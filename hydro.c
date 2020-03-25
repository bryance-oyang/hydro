#include "hydro.h"
#include "eos.h"
#include "riemann.h"
#include <float.h>

static void reflecting_boundary(struct grid *g)
{
	int i, j, k;
	int nx, ny;

	nx = g->nx;
	ny = g->ny;

	i = 0;
	for (k = 0; k < 2; k++) {
		for (j = 2; j < ny-2; j++) {
			CEL(g->prim[0],i+k,j) = CEL(g->prim[0],i+3-k,j);
			CEL(g->prim[1],i+k,j) = -1 * CEL(g->prim[1],i+3-k,j);
			CEL(g->prim[2],i+k,j) = CEL(g->prim[2],i+3-k,j);
			CEL(g->prim[3],i+k,j) = CEL(g->prim[3],i+3-k,j);
		}
	}

	i = nx-1;
	for (k = 0; k < 2; k++) {
		for (j = 2; j < ny-2; j++) {
			CEL(g->prim[0],i-k,j) = CEL(g->prim[0],i-3+k,j);
			CEL(g->prim[1],i-k,j) = -1 * CEL(g->prim[1],i-3+k,j);
			CEL(g->prim[2],i-k,j) = CEL(g->prim[2],i-3+k,j);
			CEL(g->prim[3],i-k,j) = CEL(g->prim[3],i-3+k,j);
		}
	}

	j = 0;
	for (i = 2; i < nx-2; i++) {
		for (k = 0; k < 2; k++) {
			CEL(g->prim[0],i,j+k) = CEL(g->prim[0],i,j+3-k);
			CEL(g->prim[1],i,j+k) = CEL(g->prim[1],i,j+3-k);
			CEL(g->prim[2],i,j+k) = -1 * CEL(g->prim[2],i,j+3-k);
			CEL(g->prim[3],i,j+k) = CEL(g->prim[3],i,j+3-k);
		}
	}

	j = ny-1;
	for (i = 2; i < nx-2; i++) {
		for (k = 0; k < 2; k++) {
			CEL(g->prim[0],i,j-k) = CEL(g->prim[0],i,j-3+k);
			CEL(g->prim[1],i,j-k) = CEL(g->prim[1],i,j-3+k);
			CEL(g->prim[2],i,j-k) = -1 * CEL(g->prim[2],i,j-3+k);
			CEL(g->prim[3],i,j-k) = CEL(g->prim[3],i,j-3+k);
		}
	}
}

static void boundary(struct grid *g)
{
	reflecting_boundary(g);
}


static void compute_src(struct grid *g)
{
	int i, j, nx, ny;
	int m, n;

	nx = g->nx;
	ny = g->ny;

#if _OPENMP
#pragma omp parallel for simd num_threads(NTHREAD) schedule(THREAD_SCHEDULE)
#endif /* _OPENMP */
	for (i = 2; i < nx-2; i++) {
		for (j = 2; j < ny-2; j++) {
			CEL(g->src[0],i,j) = 0;
			CEL(g->src[1],i,j) = 0;
			CEL(g->src[2],i,j) = 0;
			CEL(g->src[3],i,j) = 0;
		}
	}

#if _OPENMP
#pragma omp parallel for simd num_threads(NTHREAD) schedule(THREAD_SCHEDULE)
#endif /* _OPENMP */
	for (m = 0; m < NSCALAR; m++) {
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
#pragma omp parallel for simd num_threads(NTHREAD) schedule(THREAD_SCHEDULE)
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
#pragma omp parallel for simd num_threads(NTHREAD) schedule(THREAD_SCHEDULE)
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
#pragma omp parallel for simd num_threads(NTHREAD) schedule(THREAD_SCHEDULE)
#endif /* _OPENMP */
		for (i = 2; i < nx-2; i++) {
			for (j = 2; j < ny-2; j++) {
				CEL(g->cons[n],i,j) = CEL(g->cons_gen[n],i,j)
					+ xfac * (CEL(g->Jx[n],i,j) - CEL(g->Jx[n],i+1,j))
					+ yfac * (CEL(g->Jy[n],i,j) - CEL(g->Jy[n],i,j+1))
					+ dt * CEL(g->src[n],i,j);
			}
		}
	}

	for (m = 0; m < NSCALAR; m++) {
#if _OPENMP
#pragma omp parallel for simd num_threads(NTHREAD) schedule(THREAD_SCHEDULE)
#endif /* _OPENMP */
		for (i = 2; i < nx-2; i++) {
			for (j = 2; j < ny-2; j++) {
				CEL(g->s[m],i,j) = CEL(g->s_gen[m],i,j)
					+ xfac * (CEL(g->s_Jx[m],i,j) - CEL(g->s_Jx[m],i+1,j))
					+ yfac * (CEL(g->s_Jy[m],i,j) - CEL(g->s_Jy[m],i,j+1))
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
		compute_src(g);
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
