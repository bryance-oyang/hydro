#include "riemann.h"
#include "eos.h"
#include <math.h>

static inline double slope_lim(double r)
{
	return fmax(0, fmin(r, 1));
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
#pragma omp parallel for simd num_threads(NTHREAD) schedule(THREAD_SCHEDULE)
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
#pragma omp parallel for simd num_threads(NTHREAD) schedule(THREAD_SCHEDULE)
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
			double sqrt_Lrho, sqrt_Urho, Lmom, Umom;
			double Lcs, Ucs, cs;
			double avg_vel;
			double Lw, Uw;
			double trial_dt;

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
#pragma omp parallel for simd num_threads(NTHREAD) schedule(THREAD_SCHEDULE)
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
					CEL(J[n],i,j) = 0;
				} else if (Uw <= 0) {
					CEL(J[n],i,j) = UJ;
				} else if (Lw >= 0) {
					CEL(J[n],i,j) = LJ;
				} else {
					CEL(J[n],i,j) = (LJ*Uw - UJ*Lw + Uw*Lw*(Uq - Lq)) / (Uw - Lw);
				}
			}
		}
	}

	for (m = 0; m < NSCALAR; m++) {
#if _OPENMP
#pragma omp parallel for simd num_threads(NTHREAD) schedule(THREAD_SCHEDULE)
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
					CEL(s_J[m],i,j) = 0;
				} else if (Uw <= 0) {
					CEL(s_J[m],i,j) = UJ;
				} else if (Lw >= 0) {
					CEL(s_J[m],i,j) = LJ;
				} else {
					CEL(s_J[m],i,j) = (LJ*Uw - UJ*Lw + Uw*Lw*(Uq - Lq)) / (Uw - Lw);
				}
			}
		}
	}
}

void hllc(struct grid *g, int dir)
{
	int i, j, nx, ny;
	int m;
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

	for (i = 2; i < nx-1; i++) {
		for (j = 2; j < ny-1; j++) {
		}
	}

	for (m = 0; m < NSCALAR; m++) {
		for (i = 2; i < nx-1; i++) {
			for (j = 2; j < ny-1; j++) {
			}
		}
	}
}
