#include "riemann.h"
#include "eos.h"
#include <math.h>

extern double GRAV;

static inline double i_like_grid_noise_lim(double r)
{
	double a, b;

	if (r > 0) {
		a = 1;
		b = 2;
		return fmax(fmin(a, b*r), fmin(a*r, b));
	} else {
		return 0;
	}
}

static inline double vl_lim(double r)
{
	double fabsr;

	fabsr = fabs(r);
	return (r + fabsr) / (1 + fabsr);
}

static inline double diffuse_lim(double r)
{
	return fmax(0, fmin(1, r));
}

static inline double slope_lim(double r)
{
	if (RECONSTRUCT == 3) {
		return i_like_grid_noise_lim(r);
	} else if (RECONSTRUCT == 2) {
		return vl_lim(r);
	} else {
		return 0;
	}
}

static inline double reconstruct_dq(double Lq, double Mq, double Uq)
{
	if ((Uq - Mq) == 0) {
		return 0;
	} else {
		return slope_lim((Mq - Lq) / (Uq - Mq)) * (Uq - Mq);
	}
}

static inline double lin_interp(double x0, double x1, double y0, double y1, double x)
{
	return (y1 - y0) / (x1 - x0) * (x - x0) + y0;
}

static inline double fmin3(double a, double b, double c)
{
	return fmin(fmin(a,b),c);
}

static inline double fmin4(double a, double b, double c, double d)
{
	return fmin(fmin(fmin(a,b),c),d);
}

static inline double reconstruct_flatten(double d1, double d2)
{
	double r, flatten_start_r;

	flatten_start_r = 0.92;

	if (d2 != 0) {
		r = fabs(d1 / d2);
		if (r > flatten_start_r) {
			return fmin(1, lin_interp(flatten_start_r, 1, 0, 1, r));
		} else {
			return 0;
		}
	} else {
		return 1;
	}
}

static inline void fancy_ppm(double *pql, double *pqr, double q0, double q1, double q2, double q3, double q4)
{
	double ql, qr;
	double curvl, curvr, curvc, curvf, curv;
	double C;

	C = 1.06;

	ql = (7.0*(q1 + q2) - (q0 + q3)) / 12;
	curvl = (q0 - 2*q1 + q2);
	curvr = (q1 - 2*q2 + q3);
	//curvf = 3*(q1 - 2*ql + q2);
	curvf = 4*(q1 - 2*ql + q2);
	if (SIGN(curvl) == SIGN(curvf) && SIGN(curvf) == SIGN(curvr)) {
		curv = SIGN(curvf) * fmin3(C*fabs(curvl), C*fabs(curvr), fabs(curvf));
	} else {
		curv = 0;
	}
	//ql = 0.5 * (q1 + q2) - curv / 3;
	ql = 0.5 * (q1 + q2) - curv / 6;

	qr = (7.0*(q2 + q3) - (q1 + q4)) / 12;
	curvl = (q1 - 2*q2 + q3);
	curvr = (q2 - 2*q3 + q4);
	//curvf = 3*(q2 - 2*qr + q3);
	curvf = 4*(q2 - 2*qr + q3);
	if (SIGN(curvl) == SIGN(curvf) && SIGN(curvf) == SIGN(curvr)) {
		curv = SIGN(curvf) * fmin3(C*fabs(curvl), C*fabs(curvr), fabs(curvf));
	} else {
		curv = 0;
	}
	//qr = 0.5 * (q2 + q3) - curv / 3;
	qr = 0.5 * (q2 + q3) - curv / 6;

	double test1, test2;
	test1 = 6*(qr - ql) * (q2 - 0.5*(ql + qr));
	test2 = SQR(qr - ql);

	if ((qr - q2) * (q2 - ql) <= 0 || (q3 - q2) * (q2 - q1) <= 0) {
		curvc = (q1 - 2*q2 + q3);
		curvl = (q0 - 2*q1 + q2);
		curvr = (q2 - 2*q3 + q4);
		//curvf = 6*(ql - 2*q2 + qr);
		curvf = 4*(ql - 2*q2 + qr);
		if (SIGN(curvl) == SIGN(curvc) && SIGN(curvc) == SIGN(curvr) && SIGN(curvc) == SIGN(curvf)) {
			curv = SIGN(curvf) * fmin4(C*fabs(curvl), C*fabs(curvc), C*fabs(curvr), fabs(curvf));
		} else {
			curv = 0;
		}

		if (curvf != 0) {
			ql = q2 + (ql - q2) * curv / curvf;
			qr = q2 + (qr - q2) * curv / curvf;
		} else {
			ql = q2;
			qr = q2;
		}
	} else if (test1 > test2) {
		ql = 3*q2 - 2*qr;
	} else if (-test2 > test1) {
		qr = 3*q2 - 2*ql;
	}

	*pql = ql;
	*pqr = qr;
}

void reconstruct(struct grid *g, int step, int dir)
{
	int i, j, nx, ny;
	int m, n;
	int di, dj;

	nx = g->nx;
	ny = g->ny;

	if (dir == 0) {
		di = 1;
		dj = 0;
	} else {
		di = 0;
		dj = 1;
	}

	for (n = 0; n < 4; n++) {
#if _OPENMP
#pragma omp parallel for simd private(j) num_threads(NTHREAD) schedule(THREAD_SCHEDULE)
#endif /* _OPENMP */
		for (i = 2; i < nx-2; i++) {
			for (j = 2; j < ny-2; j++) {
				double q0, q1, q2, q3, q4;
				double ql, qr;

				q0 = CEL(g->prim[n],i-2*di,j-2*dj);
				q1 = CEL(g->prim[n],i-di,j-dj);
				q2 = CEL(g->prim[n],i,j);
				q3 = CEL(g->prim[n],i+di,j+dj);
				q4 = CEL(g->prim[n],i+2*di,j+2*dj);

				if (RECONSTRUCT == 6 && (RECONSTRUCT_BOTH || step == 1)) {
					fancy_ppm(&ql, &qr, q0, q1, q2, q3, q4);
				} else if (RECONSTRUCT == 3 && (RECONSTRUCT_BOTH || step == 1)) {
					double dq1, dq2, dq3;

					dq1 = reconstruct_dq(q0, q1, q2);
					dq2 = reconstruct_dq(q1, q2, q3);
					dq3 = reconstruct_dq(q2, q3, q4);

					ql = 0.5*(q1 + q2) - (dq2 - dq1) / 6;
					qr = 0.5*(q2 + q3) - (dq3 - dq2) / 6;

					double test1, test2;
					test1 = 6*(qr - ql) * (q2 - 0.5*(ql + qr));
					test2 = SQR(qr - ql);

					if ((qr - q2) * (q2 - ql) <= 0) {
						ql = q2;
						qr = q2;
					} else if (test1 > test2) {
						ql = 3*q2 - 2*qr;
					} else if (-test2 > test1) {
						qr = 3*q2 - 2*ql;
					}

					if (RECONSTRUCT_FLATTEN) {
						double f;

						f = reconstruct_flatten(q3 - q1, q4 - q0);
						ql = (1-f)*ql + f*q2;
						qr = (1-f)*qr + f*q2;
					}
				} else if (RECONSTRUCT == 2 && (RECONSTRUCT_BOTH || step == 1)) {
					double half_step;

					half_step = 0.5 * reconstruct_dq(q1, q2, q3);

					ql = q2 - half_step;
					qr = q2 + half_step;
				} else {
					ql = q2;
					qr = q2;
				}

				FEL(g->Uprim[n],i,j) = ql;
				FEL(g->Lprim[n],i+di,j+dj) = qr;
			}
		}
	}

	for (m = 0; m < NSCALAR; m++) {
#if _OPENMP
#pragma omp parallel for simd private(j) num_threads(NTHREAD) schedule(THREAD_SCHEDULE)
#endif /* _OPENMP */
		for (i = 2; i < nx-2; i++) {
			for (j = 2; j < ny-2; j++) {
				double q0, q1, q2, q3, q4;
				double ql, qr;

				q0 = CEL(g->s[m],i-2*di,j-2*dj);
				q1 = CEL(g->s[m],i-di,j-dj);
				q2 = CEL(g->s[m],i,j);
				q3 = CEL(g->s[m],i+di,j+dj);
				q4 = CEL(g->s[m],i+2*di,j+2*dj);

				if (RECONSTRUCT == 6 && (RECONSTRUCT_BOTH || step == 1)) {
					fancy_ppm(&ql, &qr, q0, q1, q2, q3, q4);
				} else if (RECONSTRUCT == 3 && (RECONSTRUCT_BOTH || step == 1)) {
					double dq1, dq2, dq3;

					dq1 = reconstruct_dq(q0, q1, q2);
					dq2 = reconstruct_dq(q1, q2, q3);
					dq3 = reconstruct_dq(q2, q3, q4);

					ql = 0.5*(q1 + q2) - (dq2 - dq1) / 6;
					qr = 0.5*(q2 + q3) - (dq3 - dq2) / 6;

					double test1, test2;
					test1 = 6*(qr - ql) * (q2 - 0.5*(ql + qr));
					test2 = SQR(qr - ql);

					if ((qr - q2) * (q2 - ql) <= 0) {
						ql = q2;
						qr = q2;
					} else if (test1 > test2) {
						ql = 3*q2 - 2*qr;
					} else if (-test2 > test1) {
						qr = 3*q2 - 2*ql;
					}

					if (RECONSTRUCT_FLATTEN) {
						double f;

						f = reconstruct_flatten(q3 - q1, q4 - q0);
						ql = (1-f)*ql + f*q2;
						qr = (1-f)*qr + f*q2;
					}
				} else if (RECONSTRUCT == 2 && (RECONSTRUCT_BOTH || step == 1)) {
					double half_step;

					half_step = 0.5 * reconstruct_dq(q1, q2, q3);

					ql = q2 - half_step;
					qr = q2 + half_step;
				} else {
					ql = q2;
					qr = q2;
				}

				FEL(g->Us[m],i,j) = ql;
				FEL(g->Ls[m],i+di,j+dj) = qr;
			}
		}
	}

	eos_prim_floor(g->Lprim, g->nx+1, g->ny+1);
	eos_prim_to_cons(g->Lprim, g->Lcons, g->nx+1, g->ny+1);

	eos_prim_floor(g->Uprim, g->nx+1, g->ny+1);
	eos_prim_to_cons(g->Uprim, g->Ucons, g->nx+1, g->ny+1);
}

void wavespeed(struct grid *g, int step, int dir)
{
	int i, j, nx, ny;
	int di, dj; 
	double du;

	nx = g->nx;
	ny = g->ny;

	if (dir == 0) {
		du = g->dx;
		di = 1;
		dj = 0;
	} else {
		du = g->dy;
		di = 0;
		dj = 1;
	}

	for (i = 3; i < nx-2; i++) {
		for (j = 3; j < ny-2; j++) {
			double sqrt_Lrho, sqrt_Urho, Lmom, Umom, Lv, Uv;
			double Lcs, Ucs, cs;
			double avg_vel;
			double Lw, Uw;
			double trial_dt;

			if (ROE_WAVESPEED) {
				sqrt_Lrho = sqrt(FEL(g->Lcons[0],i,j));
				Lmom = FEL(g->Lcons[1+dir],i,j);
				sqrt_Urho = sqrt(FEL(g->Ucons[0],i,j));
				Umom = FEL(g->Ucons[1+dir],i,j);

				Lcs = CEL(g->cs,i-di,j-dj);
				Ucs = CEL(g->cs,i,j);
				cs = fmax(Lcs, Ucs);
				avg_vel = (Lmom/sqrt_Lrho + Umom/sqrt_Urho) / (sqrt_Lrho + sqrt_Urho);

				Lw = avg_vel - cs;
				Uw = avg_vel + cs;
				FEL(g->Lw,i,j) = Lw;
				FEL(g->Uw,i,j) = Uw;
			} else {
				Lv = FEL(g->Lprim[1+dir],i,j);
				Uv = FEL(g->Uprim[1+dir],i,j);
				Lcs = CEL(g->cs,i-di,j-dj);
				Ucs = CEL(g->cs,i,j);

				Lw = fmin(Lv - Lcs, Uv - Ucs);
				Uw = fmax(Lv + Lcs, Uv + Ucs);
				FEL(g->Lw,i,j) = Lw;
				FEL(g->Uw,i,j) = Uw;
			}

			if (step == 0) {
				trial_dt = du / fabs(Lw);
				if (g->dt > trial_dt) {
					g->dt = trial_dt;
				}
				trial_dt = du / fabs(Uw);
				if (g->dt > trial_dt) {
					g->dt = trial_dt;
				}
			}
		}
	}
}

void hlle(struct grid *g, int dir)
{
	int i, j, nx, ny;
	int m, n;
	double **J, **s_J;

	nx = g->nx;
	ny = g->ny;

	if (dir == 0) {
		J = g->Jx;
		s_J = g->s_Jx;
	} else {
		J = g->Jy;
		s_J = g->s_Jy;
	}

	for (n = 0; n < 4; n++) {
#if _OPENMP
#pragma omp parallel for simd private(j) num_threads(NTHREAD) schedule(THREAD_SCHEDULE)
#endif /* _OPENMP */
		for (i = 3; i < nx-2; i++) {
			for (j = 3; j < ny-2; j++) {
				double Lw, Uw;
				double Lv, Uv;
				double Lpress, Upress;

				Lw = FEL(g->Lw,i,j);
				Uw = FEL(g->Uw,i,j);

				Lv = FEL(g->Lprim[1+dir],i,j);
				Uv = FEL(g->Uprim[1+dir],i,j);

				Lpress = FEL(g->Lprim[3],i,j);
				Upress = FEL(g->Uprim[3],i,j);

				double Lq, Uq;
				double LJ, UJ;

				Lq = FEL(g->Lcons[n],i,j);
				Uq = FEL(g->Ucons[n],i,j);
				LJ = Lq * Lv;
				UJ = Uq * Uv;

				if (n == 1+dir) {
					LJ += Lpress;
					UJ += Upress;
				}
				if (n == 3) {
					LJ += Lpress * Lv;
					UJ += Upress * Uv;
				}

				if (Lw == 0 && Uw == 0) {
					FEL(J[n],i,j) = 0;
				} else if (Uw <= 0) {
					FEL(J[n],i,j) = UJ;
				} else if (Lw >= 0) {
					FEL(J[n],i,j) = LJ;
				} else {
					FEL(J[n],i,j) = (LJ*Uw - UJ*Lw + Uw*Lw*(Uq - Lq)) / (Uw - Lw);
				}
			}
		}
	}

	for (m = 0; m < NSCALAR; m++) {
#if _OPENMP
#pragma omp parallel for simd private(j) num_threads(NTHREAD) schedule(THREAD_SCHEDULE)
#endif /* _OPENMP */
		for (i = 3; i < nx-2; i++) {
			for (j = 3; j < ny-2; j++) {
				double Lw, Uw;
				double Lv, Uv;

				Lw = FEL(g->Lw,i,j);
				Uw = FEL(g->Uw,i,j);

				Lv = FEL(g->Lprim[1+dir],i,j);
				Uv = FEL(g->Uprim[1+dir],i,j);

				double Lq, Uq;
				double LJ, UJ;

				Lq = FEL(g->Ls[m],i,j);
				Uq = FEL(g->Us[m],i,j);
				LJ = Lq * Lv;
				UJ = Uq * Uv;

				if (Lw == 0 && Uw == 0) {
					FEL(s_J[m],i,j) = 0;
				} else if (Uw <= 0) {
					FEL(s_J[m],i,j) = UJ;
				} else if (Lw >= 0) {
					FEL(s_J[m],i,j) = LJ;
				} else {
					FEL(s_J[m],i,j) = (LJ*Uw - UJ*Lw + Uw*Lw*(Uq - Lq)) / (Uw - Lw);
				}
			}
		}
	}
}

void hllc(struct grid *g, int dir)
{
	int i, j, nx, ny;
	double **J, **s_J;

	nx = g->nx;
	ny = g->ny;

	if (dir == 0) {
		J = g->Jx;
		s_J = g->s_Jx;
	} else {
		J = g->Jy;
		s_J = g->s_Jy;
	}

#if _OPENMP
#pragma omp parallel for simd private(j) num_threads(NTHREAD) schedule(THREAD_SCHEDULE)
#endif /* _OPENMP */
	for (i = 3; i < nx-2; i++) {
		for (j = 3; j < ny-2; j++) {
			int m, n;
			double Lw, Uw, Mw;
			double Lv, Uv;
			double Lpress, Upress, Mpress;
			double Lrho, Urho, rho2, rho3;
			double Le, Ue, e2, e3;

			Lw = FEL(g->Lw,i,j);
			Uw = FEL(g->Uw,i,j);

			if (Lw == 0 && Uw == 0) {
				for (n = 0; n < 4; n++) {
					FEL(J[n],i,j) = 0;
				}
				for (m = 0; m < NSCALAR; m++) {
					FEL(s_J[m],i,j) = 0;
				}
				continue;
			}

			Lrho = FEL(g->Lprim[0],i,j);
			Urho = FEL(g->Uprim[0],i,j);
			Lv = FEL(g->Lprim[1+dir],i,j);
			Uv = FEL(g->Uprim[1+dir],i,j);
			Lpress = FEL(g->Lprim[3],i,j);
			Upress = FEL(g->Uprim[3],i,j);
			Le = FEL(g->Lcons[3],i,j);
			Ue = FEL(g->Ucons[3],i,j);

			if (Uw <= 0) {
				for (n = 0; n < 4; n++) {
					FEL(J[n],i,j) = FEL(g->Ucons[n],i,j) * Uv;
					if (n == 1+dir) {
						FEL(J[n],i,j) += Upress;
					}
					if (n == 3) {
						FEL(J[n],i,j) += Upress * Uv;
					}
				}
				for (m = 0; m < NSCALAR; m++) {
					FEL(s_J[m],i,j) = FEL(g->Us[m],i,j) * Uv;
				}
				continue;
			}
			if (Lw >= 0) {
				for (n = 0; n < 4; n++) {
					FEL(J[n],i,j) = FEL(g->Lcons[n],i,j) * Lv;
					if (n == 1+dir) {
						FEL(J[n],i,j) += Lpress;
					}
					if (n == 3) {
						FEL(J[n],i,j) += Lpress * Lv;
					}
				}
				for (m = 0; m < NSCALAR; m++) {
					FEL(s_J[m],i,j) = FEL(g->Ls[m],i,j) * Lv;
				}
				continue;
			}

			Mw = (Urho*Uv*(Uv-Uw) + Upress - Lrho*Lv*(Lv-Lw) - Lpress) / (Urho*(Uv-Uw) - Lrho*(Lv-Lw));
			rho2 = Lrho * (Lv - Lw) / (Mw - Lw);
			rho3 = Urho * (Uv - Uw) / (Mw - Uw);

			if (Mw > 0) {
				Mpress = Lrho*SQR(Lv) + Lpress - Lw*Lrho*Lv - rho2*SQR(Mw) + Lw*rho2*Mw;
			} else {
				Mpress = Urho*SQR(Uv) + Upress - Uw*Urho*Uv - rho3*SQR(Mw) + Uw*rho3*Mw;
			}

			if (Mw == 0) {
				for (n = 0; n < 4; n++) {
					FEL(J[n],i,j) = 0;
					if (n == 1+dir) {
						FEL(J[n],i,j) += Mpress;
					}
				}
				for (m = 0; m < NSCALAR; m++) {
					FEL(s_J[m],i,j) = 0;
				}
				continue;
			}

			if (Mw < 0) {
				e3 = (Uv*(Ue+Upress) - Uw*Ue - Mw*Mpress) / (Mw - Uw);

				FEL(J[0],i,j) = rho3 * Mw;
				if (dir == 0) {
					FEL(J[1],i,j) = rho3 * SQR(Mw) + Mpress;
					FEL(J[2],i,j) = rho3 * FEL(g->Uprim[2],i,j) * Mw;
				} else {
					FEL(J[1],i,j) = rho3 * FEL(g->Uprim[1],i,j) * Mw;
					FEL(J[2],i,j) = rho3 * SQR(Mw) + Mpress;
				}
				FEL(J[3],i,j) = (e3 + Mpress) * Mw;

				for (m = 0; m < NSCALAR; m++) {
					FEL(s_J[m],i,j) = (Uv - Uw) / (Mw - Uw) * FEL(g->Us[m],i,j) * Mw;
				}
			} else {
				e2 = (Lv*(Le+Lpress) - Lw*Le - Mw*Mpress) / (Mw - Lw);

				FEL(J[0],i,j) = rho2 * Mw;
				if (dir == 0) {
					FEL(J[1],i,j) = rho2 * SQR(Mw) + Mpress;
					FEL(J[2],i,j) = rho2 * FEL(g->Lprim[2],i,j) * Mw;
				} else {
					FEL(J[1],i,j) = rho2 * FEL(g->Lprim[1],i,j) * Mw;
					FEL(J[2],i,j) = rho2 * SQR(Mw) + Mpress;
				}
				FEL(J[3],i,j) = (e2 + Mpress) * Mw;

				for (m = 0; m < NSCALAR; m++) {
					FEL(s_J[m],i,j) = (Lv - Lw) / (Mw - Lw) * FEL(g->Ls[m],i,j) * Mw;
				}
			}
		}
	}
}
