#include "riemann.h"
#include "eos.h"
#include <math.h>

extern double GRAV;

static inline double grid_noise_lim(double r)
{
	double a, b;

	if (r > 0) {
		a = 1;
		b = 1.36;
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

static inline double slope_lim(double r)
{
	if (RECONSTRUCT) {
		return vl_lim(r);
	} else {
		return 0;
	}
}

void reconstruct(struct grid *g, int dir)
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
		for (i = 1; i < nx-1; i++) {
			for (j = 1; j < ny-1; j++) {
				double q0, q1, q2;
				double half_step, Lq, Uq;

				q0 = CEL(g->prim[n],i-di,j-dj);
				q1 = CEL(g->prim[n],i,j);
				q2 = CEL(g->prim[n],i+di,j+dj);

				if ((q2 - q1) == 0) {
					half_step = 0;
				} else {
					half_step = 0.5 * slope_lim((q1 - q0) / (q2 - q1)) * (q2 - q1);
				}

				Lq = q1 + half_step;
				Uq = q1 - half_step;

				FEL(g->Uprim[n],i,j) = Uq;
				FEL(g->Lprim[n],i+di,j+dj) = Lq;
			}
		}
	}

	for (m = 0; m < NSCALAR; m++) {
#if _OPENMP
#pragma omp parallel for simd private(j) num_threads(NTHREAD) schedule(THREAD_SCHEDULE)
#endif /* _OPENMP */
		for (i = 1; i < nx-1; i++) {
			for (j = 1; j < ny-1; j++) {
				double q0, q1, q2;
				double half_step, Lq, Uq;

				q0 = CEL(g->s[m],i-di,j-dj);
				q1 = CEL(g->s[m],i,j);
				q2 = CEL(g->s[m],i+di,j+dj);

				if ((q2 - q1) == 0) {
					half_step = 0;
				} else {
					half_step = 0.5 * slope_lim((q1 - q0) / (q2 - q1)) * (q2 - q1);
				}

				Lq = q1 + half_step;
				Uq = q1 - half_step;

				FEL(g->Us[m],i,j) = Uq;
				FEL(g->Ls[m],i+di,j+dj) = Lq;
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

	for (i = 2; i < nx-1; i++) {
		for (j = 2; j < ny-1; j++) {
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
		for (i = 2; i < nx-1; i++) {
			for (j = 2; j < ny-1; j++) {
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
		for (i = 2; i < nx-1; i++) {
			for (j = 2; j < ny-1; j++) {
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
	for (i = 2; i < nx-1; i++) {
		for (j = 2; j < ny-1; j++) {
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
