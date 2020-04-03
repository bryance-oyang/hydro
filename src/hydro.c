#include "hydro.h"
#include "eos.h"
#include "riemann.h"
#include "boundary.h"
#include <float.h>
#include <math.h>

#if WIND_TUNNEL == 1
#include "wind_tunnel.h"
#endif
#if BINARY == 1
#include "binary.h"
#endif

extern double GAMMA;
extern double GRAV;

void boundary(struct grid *g, int step)
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
#if BINARY == 1
		binary_boundary(g, step);
#endif
		empty_boundary_left(g);
		empty_boundary_right(g);
		empty_boundary_bot(g);
		empty_boundary_top(g);
	} else if (LINEAR_WAVE_TEST_X || LINEAR_WAVE_TEST_Y || LINEAR_WAVE_TEST_XY) {
		periodic_boundary_left(g);
		periodic_boundary_right(g);
		periodic_boundary_bot(g);
		periodic_boundary_top(g);
	} else if (WIND_TUNNEL) {
#if WIND_TUNNEL == 1
		wind_tunnel_boundary(g);
#endif
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

static double potential(double t, double x, double y)
{
#if BINARY == 1
	double r, grav1, grav2, spin, r2;

	r = sqrt(SQR(x) + SQR(y));
	if (r <= M1_CUTOFF) {
		grav1 = -GM1 / M1_CUTOFF;
	} else {
		grav1 = -GM1 / r;
	}
	if (BIN_ROT_FRAME) {
		r2 = sqrt(SQR(x - BIN_SEP*cos(BIN_ANGLE)) + SQR(y - BIN_SEP*sin(BIN_ANGLE)));
		spin = -0.5 * (SQR(x - BIN_COM*cos(BIN_ANGLE)) + SQR(y - BIN_COM*sin(BIN_ANGLE))) * SQR(BIN_OMEGA);
	} else {
		r2 = sqrt(SQR(x - BIN_SEP*cos(BIN_OMEGA*t + BIN_ANGLE)) + SQR(y - BIN_SEP*sin(BIN_OMEGA*t + BIN_ANGLE)));
		spin = 0;
	}
	if (r2 <= M2_CUTOFF) {
		grav2 = -GM2 / M2_CUTOFF;
	} else {
		grav2 = -GM2 / r2;

	}
	return grav1 + grav2 + spin;
#else
	return GRAV * y;
#endif
}

static void compute_src(struct grid *g, int step)
{
	int i, j, nx, ny;
	int m;
	double dx, dy;
	double t;

	t = g->time + step * g->dt / 2;
	nx = g->nx;
	ny = g->ny;
	dx = g->dx;
	dy = g->dy;

#if _OPENMP
#pragma omp parallel for simd private(j) num_threads(NTHREAD) schedule(THREAD_SCHEDULE)
#endif /* _OPENMP */
	for (i = 3; i < nx-3; i++) {
		for (j = 3; j < ny-3; j++) {
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

			div_rhovpot = (potential(t, Ux, y) * FEL(g->Jx[0],i+1,j) - potential(t, Lx, y) * FEL(g->Jx[0],i,j)) / dx
				+ (potential(t, x, Uy) * FEL(g->Jy[0],i,j+1) - potential(t, x, Ly) * FEL(g->Jy[0],i,j)) / dy;

			pot_cc = potential(t, x, y);

			CEL(g->src[0],i,j) = 0;
			CEL(g->src[1],i,j) = CEL(g->prim[0],i,j) * (potential(t, Lx, y) - potential(t, Ux, y)) / dx;
			CEL(g->src[2],i,j) = CEL(g->prim[0],i,j) * (potential(t, x, Ly) - potential(t, x, Uy)) / dy;
			if (FANCY_POT_NRG) {
				CEL(g->src[3],i,j) = pot_cc * div_rhov - div_rhovpot;
			} else {
				CEL(g->src[3],i,j) = CEL(g->cons[1],i,j) * (potential(t, Lx, y) - potential(t, Ux, y)) / dx
					+ CEL(g->cons[2],i,j) * (potential(t, x, Ly) - potential(t, x, Uy)) / dy;
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
					if (a == 0) {
						continue;
					}

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
		for (i = 3; i < nx-3; i++) {
			for (j = 3; j < ny-3; j++) {
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
		for (i = 3; i < nx-3; i++) {
			for (j = 3; j < ny-3; j++) {
				CEL(g->prim_gen[n],i,j) = CEL(g->prim[n],i,j);
				CEL(g->cons_gen[n],i,j) = CEL(g->cons[n],i,j);
			}
		}
	}

	for (m = 0; m < NSCALAR; m++) {
#if _OPENMP
#pragma omp parallel for simd private(j) num_threads(NTHREAD) schedule(THREAD_SCHEDULE)
#endif /* _OPENMP */
		for (i = 3; i < nx-3; i++) {
			for (j = 3; j < ny-3; j++) {
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
	double ssprk3_weight[3];

	if (step == 0) {
		ssprk3_weight[0] = 1;
		ssprk3_weight[1] = 0;
		ssprk3_weight[2] = 1;
	} else if (step == 1) {
		ssprk3_weight[0] = 3.0/4.0;
		ssprk3_weight[1] = 1.0/4.0;
		ssprk3_weight[2] = 1.0/4.0;
	} else {
		ssprk3_weight[0] = 1.0/3.0;
		ssprk3_weight[1] = 2.0/3.0;
		ssprk3_weight[2] = 2.0/3.0;
	}

	nx = g->nx;
	ny = g->ny;

	dt = g->dt;
	xfac = 1.0 / g->dx;
	yfac = 1.0 / g->dy;

	for (n = 0; n < 4; n++) {
#if _OPENMP
#pragma omp parallel for simd private(j) num_threads(NTHREAD) schedule(THREAD_SCHEDULE)
#endif /* _OPENMP */
		for (i = 3; i < nx-3; i++) {
			for (j = 3; j < ny-3; j++) {
				double deriv;

				deriv = xfac * (FEL(g->Jx[n],i,j) - FEL(g->Jx[n],i+1,j))
					+ yfac * (FEL(g->Jy[n],i,j) - FEL(g->Jy[n],i,j+1))
					+ CEL(g->src[n],i,j);

				CEL(g->cons[n],i,j) = ssprk3_weight[0] * CEL(g->cons_gen[n],i,j)
					+ ssprk3_weight[1] * CEL(g->cons[n],i,j)
					+ ssprk3_weight[2] * deriv * dt;
			}
		}
	}

	for (m = 0; m < NSCALAR; m++) {
#if _OPENMP
#pragma omp parallel for simd private(j) num_threads(NTHREAD) schedule(THREAD_SCHEDULE)
#endif /* _OPENMP */
		for (i = 3; i < nx-3; i++) {
			for (j = 3; j < ny-3; j++) {
				double deriv;

				deriv = xfac * (FEL(g->s_Jx[m],i,j) - FEL(g->s_Jx[m],i+1,j))
					+ yfac * (FEL(g->s_Jy[m],i,j) - FEL(g->s_Jy[m],i,j+1))
					+ CEL(g->s_src[m],i,j);

				CEL(g->s[m],i,j) = ssprk3_weight[0] * CEL(g->s_gen[m],i,j)
					+ ssprk3_weight[1] * CEL(g->s[m],i,j)
					+ ssprk3_weight[2] * deriv * dt;
			}
		}
	}
}

static void compute_J(struct grid *g, int step, int dir)
{
	reconstruct(g, step, dir);
	wavespeed(g, step, dir);
	if (HLLC) {
		hllc(g, dir);
	} else {
		hlle(g, dir);
	}
}

static void compute_Jx(struct grid *g, int step)
{
	compute_J(g, step, 0);
}

static void compute_Jy(struct grid *g, int step)
{
	compute_J(g, step, 1);
}

void advance_timestep(struct grid *g)
{
	int step;

	copy_gen(g);
	g->dt = DBL_MAX;

	for (step = 0; step < 3; step++) {
		eos_sound_speed(g);
		compute_Jx(g, step);
		compute_Jy(g, step);
		compute_src(g, step);
		if (step == 0) {
			g->dt *= CFL_NUM;
		}
		add_flux_div_src(g, step);
		eos_cons_to_prim(g->cons, g->prim, g->nx, g->ny);
		eos_prim_floor(g->prim, g->nx, g->ny);
		boundary(g, step);
		eos_prim_floor(g->prim, g->nx, g->ny);
		eos_prim_to_cons(g->prim, g->cons, g->nx, g->ny);
	}

	g->time += g->dt;
}
