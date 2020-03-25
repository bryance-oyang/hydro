#include "eos.h"
#include "def.h"
#include <math.h>

void eos_prim_floor(double **prim, int nx, int ny)
{
	int i, j;

#if _OPENMP
#pragma omp parallel for simd num_threads(NTHREAD) schedule(THREAD_SCHEDULE)
#endif /* _OPENMP */
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			if (CEL(prim[0],i,j) < RHO_FLOOR) {
				CEL(prim[0],i,j) = RHO_FLOOR;
			}

			if (CEL(prim[3],i,j) < PRESS_FLOOR) {
				CEL(prim[3],i,j) = PRESS_FLOOR;
			}
		}
	}
}

void eos_prim_to_cons(double **prim, double **cons, int nx, int
		ny)
{
	int i, j;

#if _OPENMP
#pragma omp parallel for simd num_threads(NTHREAD) schedule(THREAD_SCHEDULE)
#endif /* _OPENMP */
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			double rho, vx, vy, press;

			rho = CEL(prim[0],i,j);
			vx = CEL(prim[1],i,j);
			vy = CEL(prim[2],i,j);
			press = CEL(prim[3],i,j);

			CEL(cons[0],i,j) = rho;
			CEL(cons[1],i,j) = rho * vx;
			CEL(cons[2],i,j) = rho * vy;
			CEL(cons[3],i,j) = 0.5*rho*(SQR(vx) + SQR(vy)) + press / (GAMMA - 1);
		}
	}
}

void eos_cons_to_prim(double **cons, double **prim, int nx, int
		ny)
{
	int i, j;

#if _OPENMP
#pragma omp parallel for simd num_threads(NTHREAD) schedule(THREAD_SCHEDULE)
#endif /* _OPENMP */
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			double rho, px, py, nrg;

			rho = CEL(cons[0],i,j);
			px = CEL(cons[1],i,j);
			py = CEL(cons[2],i,j);
			nrg = CEL(cons[3],i,j);

			CEL(prim[0],i,j) = rho;
			CEL(prim[1],i,j) = px / rho;
			CEL(prim[2],i,j) = py / rho;
			CEL(prim[3],i,j) = (nrg - 0.5*(SQR(px) + SQR(py)) / rho) * (GAMMA-1);
		}
	}
}

void eos_sound_speed(struct grid *g)
{
	int i, j, nx, ny;

	nx = g->nx;
	ny = g->ny;

	for (i = 1; i < nx-1; i++) {
		for (j = 1; j < ny-1; j++) {
			double rho, press;

			rho = CEL(g->prim[0],i,j);
			press = CEL(g->prim[3],i,j);

			CEL(g->cs,i,j) = sqrt(GAMMA * press / rho);
		}
	}
}
