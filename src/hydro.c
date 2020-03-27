#include "hydro.h"
#include "eos.h"
#include "riemann.h"
#include "boundary.h"
#include <float.h>
#include <math.h>

extern double GRAV;

void boundary(struct grid *g)
{
	if (KH_INSTAB || SUPERSONIC || BLAST) {
		periodic_boundary_left(g);
		periodic_boundary_right(g);
		periodic_boundary_bot(g);
		periodic_boundary_top(g);
	} else if (RT_INSTAB) {
		periodic_boundary_left(g);
		periodic_boundary_right(g);
		reflecting_boundary_bot(g);
		reflecting_boundary_top(g);
	} else if (ATMOSPHERE) {
		smooth_boundary_left(g);
		smooth_boundary_right(g);
		reflecting_boundary_bot(g);
		smooth_boundary_top(g);
	} else if (BINARY) {
		smooth_boundary_left(g);
		smooth_boundary_right(g);
		smooth_boundary_bot(g);
		smooth_boundary_top(g);
	} else {
		reflecting_boundary_left(g);
		reflecting_boundary_right(g);
		reflecting_boundary_bot(g);
		reflecting_boundary_top(g);
	}
}

static double potential(double x, double y)
{
	if (BINARY) {
		double r, grav1, grav2, spin;

		r = sqrt(SQR(x) + SQR(y));
		if (r <= M1_CUTOFF) {
			grav1 = -GM1 / M1_CUTOFF;
		} else {
			grav1 = -GM1 / sqrt(SQR(x) + SQR(y));
		}
		grav2 = -GM2 / sqrt(SQR(x - BIN_SEP) + SQR(y));
		spin = -0.5 * (SQR(x - BIN_COM) + SQR(y)) * SQR(BIN_OMEGA);
		return grav1 + grav2 + spin;
	} else {
		return GRAV * y;
	}
}

static void compute_src(struct grid *g, int step)
{
	int i, j, nx, ny;
	int m;
	double dx, dy;

	nx = g->nx;
	ny = g->ny;
	dx = g->dx;
	dy = g->dy;

#if _OPENMP
#pragma omp parallel for simd private(j) num_threads(NTHREAD) schedule(THREAD_SCHEDULE)
#endif /* _OPENMP */
	for (i = 2; i < nx-2; i++) {
		for (j = 2; j < ny-2; j++) {
			double Lx, Ux, Ly, Uy, x, y;
			double div_rhov, div_rhovpot, pot_cc;

			Lx = FEL(g->x_fc,i,j);
			Ly = FEL(g->y_fc,i,j);
			Ux = FEL(g->x_fc,i+1,j);
			Uy = FEL(g->y_fc,i,j+1);
			x = CEL(g->x_cc,i,j);
			y = CEL(g->y_cc,i,j);

			div_rhov = (FEL(g->Jx[0],i+1,j) - FEL(g->Jx[0],i,j)) / dx
				+ (FEL(g->Jy[0],i,j+1) - FEL(g->Jy[0],i,j)) / dy;

			div_rhovpot = (potential(Ux, y) * FEL(g->Jx[0],i+1,j) - potential(Lx, y) * FEL(g->Jx[0],i,j)) / dx
				+ (potential(x, Uy) * FEL(g->Jy[0],i,j+1) - potential(x, Ly) * FEL(g->Jy[0],i,j)) / dy;

			pot_cc = potential(x, y);

			CEL(g->src[0],i,j) = 0;
			CEL(g->src[1],i,j) = CEL(g->prim[0],i,j) * (potential(Lx, y) - potential(Ux, y)) / dx;
			CEL(g->src[2],i,j) = CEL(g->prim[0],i,j) * (potential(x, Ly) - potential(x, Uy)) / dy;
			if (FANCY_POT_NRG) {
				CEL(g->src[3],i,j) = pot_cc * div_rhov - div_rhovpot;
			} else {
				CEL(g->src[3],i,j) = CEL(g->cons[1],i,j) * (potential(Lx, y) - potential(Ux, y)) / dx
					+ CEL(g->cons[2],i,j) * (potential(x, Ly) - potential(x, Uy)) / dy;
			}

			if (step == 0) {
				int dir;

				for (dir = 0; dir < 2; dir++) {
					double trial_dt, a, v, sqrt_part;
					double ds;

					if (dir == 0) {
						ds = dx;
					} else {
						ds = dy;
					}

					a = CEL(g->src[1+dir],i,j) / CEL(g->prim[0],i,j);
					v = CEL(g->prim[1+dir],i,j);
					sqrt_part = sqrt(2*a*ds + SQR(v));
					trial_dt = fmin(fabs((-v - sqrt_part) / a), fabs((-v + sqrt_part) / a));

					if (g->dt > trial_dt) {
#if _OPENMP
#pragma omp atomic write
#endif /* _OPENMP */
						g->dt = trial_dt;
					}
				}
			}
		}
	}

	for (m = 0; m < NSCALAR; m++) {
#if _OPENMP
#pragma omp parallel for simd private(j) num_threads(NTHREAD) schedule(THREAD_SCHEDULE)
#endif /* _OPENMP */
		for (i = 2; i < nx-2; i++) {
			for (j = 2; j < ny-2; j++) {
				CEL(g->s_src[m],i,j) = 0;
			}
		}
	}
}

static void copy_gen(struct grid *g)
{
	int i, j, nx, ny;
	int m, n;

	nx = g->nx;
	ny = g->ny;

	for (n = 0; n < 4; n++) {
#if _OPENMP
#pragma omp parallel for simd private(j) num_threads(NTHREAD) schedule(THREAD_SCHEDULE)
#endif /* _OPENMP */
		for (i = 2; i < nx-2; i++) {
			for (j = 2; j < ny-2; j++) {
				CEL(g->prim_gen[n],i,j) = CEL(g->prim[n],i,j);
				CEL(g->cons_gen[n],i,j) = CEL(g->cons[n],i,j);
			}
		}
	}

	for (m = 0; m < NSCALAR; m++) {
#if _OPENMP
#pragma omp parallel for simd private(j) num_threads(NTHREAD) schedule(THREAD_SCHEDULE)
#endif /* _OPENMP */
		for (i = 2; i < nx-2; i++) {
			for (j = 2; j < ny-2; j++) {
				CEL(g->s_gen[m],i,j) = CEL(g->s[m],i,j);
			}
		}
	}
}

static void add_flux_div_src(struct grid *g, int step)
{
	double dt;
	double xfac, yfac;
	int i, j, nx, ny;
	int m, n;

	if (step == 0) {
		dt = g->dt / 2;
	} else {
		dt = g->dt;
	}

	nx = g->nx;
	ny = g->ny;

	xfac = dt / g->dx;
	yfac = dt / g->dy;

	for (n = 0; n < 4; n++) {
#if _OPENMP
#pragma omp parallel for simd private(j) num_threads(NTHREAD) schedule(THREAD_SCHEDULE)
#endif /* _OPENMP */
		for (i = 2; i < nx-2; i++) {
			for (j = 2; j < ny-2; j++) {
				CEL(g->cons[n],i,j) = CEL(g->cons_gen[n],i,j)
					+ xfac * (FEL(g->Jx[n],i,j) - FEL(g->Jx[n],i+1,j))
					+ yfac * (FEL(g->Jy[n],i,j) - FEL(g->Jy[n],i,j+1))
					+ dt * CEL(g->src[n],i,j);
			}
		}
	}

	if (BINARY) {
#if _OPENMP
#pragma omp parallel for simd private(j) num_threads(NTHREAD) schedule(THREAD_SCHEDULE)
#endif /* _OPENMP */
		for (i = 2; i < nx-2; i++) {
			for (j = 2; j < ny-2; j++) {
				double rho, vx, vy, dtomega;

				dtomega = dt * BIN_OMEGA;
				rho = CEL(g->cons[0],i,j);
				vx = CEL(g->cons[1],i,j) / rho;
				vy = CEL(g->cons[2],i,j) / rho;
				CEL(g->cons[1],i,j) = rho
					* ((1 - SQR(dtomega))*vx + 2*dtomega*vy)
					/ (1 + SQR(dtomega));
				CEL(g->cons[2],i,j) = rho
					* (-2*dtomega*vx + (1 - SQR(dtomega))*vy)
					/ (1 + SQR(dtomega));
			}
		}
	}

	for (m = 0; m < NSCALAR; m++) {
#if _OPENMP
#pragma omp parallel for simd private(j) num_threads(NTHREAD) schedule(THREAD_SCHEDULE)
#endif /* _OPENMP */
		for (i = 2; i < nx-2; i++) {
			for (j = 2; j < ny-2; j++) {
				CEL(g->s[m],i,j) = CEL(g->s_gen[m],i,j)
					+ xfac * (FEL(g->s_Jx[m],i,j) - FEL(g->s_Jx[m],i+1,j))
					+ yfac * (FEL(g->s_Jy[m],i,j) - FEL(g->s_Jy[m],i,j+1))
					+ dt * CEL(g->s_src[m],i,j);
			}
		}
	}
}

static void compute_Jx(struct grid *g, int step)
{
	reconstruct(g, 0);
	wavespeed(g, step, 0);
	if (HLLC) {
		hllc(g, 0);
	} else {
		hlle(g, 0);
	}
}

static void compute_Jy(struct grid *g, int step)
{
	reconstruct(g, 1);
	wavespeed(g, step, 1);
	if (HLLC) {
		hllc(g, 1);
	} else {
		hlle(g, 1);
	}
}

void advance_timestep(struct grid *g)
{
	int step;

	copy_gen(g);
	g->dt = DBL_MAX;

	for (step = 0; step < 2; step++) {
		eos_sound_speed(g);
		compute_Jx(g, step);
		compute_Jy(g, step);
		compute_src(g, step);
		if (step == 0) {
			g->dt *= CFL_NUM;
		}
		add_flux_div_src(g, step);
		eos_cons_to_prim(g->cons, g->prim, g->nx, g->ny);
		boundary(g);
		eos_prim_floor(g->prim, g->nx, g->ny);
		eos_prim_to_cons(g->prim, g->cons, g->nx, g->ny);
	}

	g->time += g->dt;
}
