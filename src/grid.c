#include "grid.h"
#include "emalloc.h"
#include "eos.h"
#include <stdlib.h>
#include <math.h>

extern double XMIN;
extern double YMIN;
extern double GRAV;

void init_grid(struct grid *g)
{
	int i, j, nx, ny;
	double x, y, dx, dy;

	nx = g->nx;
	ny = g->ny;
	dx = g->dx;
	dy = g->dy;

	for (i = 2; i < nx-1; i++) {
		for (j = 2; j < ny-1; j++) {
			CEL(g->x_cc,i,j) = XMIN + (i+0.5)*dx;
			CEL(g->y_cc,i,j) = YMIN + (j+0.5)*dy;
			FEL(g->x_fc,i,j) = XMIN + i*dx;
			FEL(g->y_fc,i,j) = YMIN + j*dy;
		}
	}

	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			x = CEL(g->x_cc,i,j);
			y = CEL(g->y_cc,i,j);

			if (KH_INSTAB) {
				if (fabs(y) > 0.25) {
					CEL(g->prim[0],i,j) = 1;
					CEL(g->prim[1],i,j) = 0.5;
					CEL(g->prim[2],i,j) = 0;
					CEL(g->prim[3],i,j) = 2.5;
				} else {
					CEL(g->prim[0],i,j) = 2;
					CEL(g->prim[1],i,j) = -0.5;
					CEL(g->prim[2],i,j) = 0;
					CEL(g->prim[3],i,j) = 2.5;
				}
				double pk2pk = 0.05;
				CEL(g->prim[1],i,j) += pk2pk * (double)rand() / RAND_MAX - pk2pk / 2;
				CEL(g->prim[2],i,j) += pk2pk * (double)rand() / RAND_MAX - pk2pk / 2;
			}

			if (RT_INSTAB) {
				double rho;
				if (y > 0) {
					rho = 2;
					CEL(g->prim[0],i,j) = rho;
					CEL(g->prim[1],i,j) = 0;
					CEL(g->prim[2],i,j) = 0.01*(1 + cos(4*PI*x))*(1 + cos(3*PI*y))/4;
					CEL(g->prim[3],i,j) = 2.5 - GRAV*rho*y;
				} else {
					rho = 1;
					CEL(g->prim[0],i,j) = rho;
					CEL(g->prim[1],i,j) = 0;
					CEL(g->prim[2],i,j) = 0.01*(1 + cos(4*PI*x))*(1 + cos(3*PI*y))/4;
					CEL(g->prim[3],i,j) = 2.5 - GRAV*rho*y;
				}
			}

			if (SOD_SHOCK) {
				if (i > nx/2) {
					CEL(g->prim[0],i,j) = 1;
					CEL(g->prim[1],i,j) = 0;
					CEL(g->prim[2],i,j) = 0;
					CEL(g->prim[3],i,j) = 1;
				} else {
					CEL(g->prim[0],i,j) = 0.125;
					CEL(g->prim[1],i,j) = 0;
					CEL(g->prim[2],i,j) = 0;
					CEL(g->prim[3],i,j) = 0.1;
				}
			}

			if (BLAST) {
				if (i == nx/2 && j == ny/2) {
					CEL(g->prim[0],i,j) = 100;
					CEL(g->prim[1],i,j) = 0;
					CEL(g->prim[2],i,j) = 0;
					CEL(g->prim[3],i,j) = 100;
				} else {
					CEL(g->prim[0],i,j) = 1;
					CEL(g->prim[1],i,j) = 0;
					CEL(g->prim[2],i,j) = 0;
					CEL(g->prim[3],i,j) = 1;
				}
			}

			if (STATIC_GRAV_TEST) {
				double rho;
				rho = 1;
				if (i == nx/2 && j == ny/2) {
					CEL(g->prim[0],i,j) = 100;
					CEL(g->prim[1],i,j) = 0;
					CEL(g->prim[2],i,j) = 0;
					CEL(g->prim[3],i,j) = 100;
				} else {
					CEL(g->prim[0],i,j) = rho;
					CEL(g->prim[1],i,j) = 0;
					CEL(g->prim[2],i,j) = 0;
					CEL(g->prim[3],i,j) = 1 - GRAV * rho * y;
				}
			}

			if (SUPERSONIC) {
				double rho;
				rho = 1;
				if (fabs(x) < 0.1 && fabs(y) < 0.1) {
					CEL(g->prim[0],i,j) = 2;
					CEL(g->prim[1],i,j) = -5;
					CEL(g->prim[2],i,j) = -3;
					CEL(g->prim[3],i,j) = 2;
				} else {
					CEL(g->prim[0],i,j) = rho;
					CEL(g->prim[1],i,j) = 0;
					CEL(g->prim[2],i,j) = 0;
					CEL(g->prim[3],i,j) = 1 - GRAV * rho * y;
				}
			}
		}
	}

	eos_prim_to_cons(g->prim, g->cons, nx, ny);
}

struct grid *alloc_grid(int nx, int ny, double dx, double dy)
{
	int m, n;
	struct grid *g;

	g = emalloc(sizeof(*g));
	g->time = 0;
	g->nx = nx;
	g->ny = ny;
	g->dx = dx;
	g->dy = dy;

#define C_ALLOC(q) g->q = emalloc(nx * ny * sizeof(*g->q))
#define F_ALLOC(q) g->q = emalloc((nx+1) * (ny+1) * sizeof(*g->q))
	for (n = 0; n < 4; n++) {
		C_ALLOC(cons[n]);
		C_ALLOC(prim[n]);
		C_ALLOC(cons_gen[n]);
		C_ALLOC(prim_gen[n]);

		F_ALLOC(Lprim[n]);
		F_ALLOC(Uprim[n]);
		F_ALLOC(Lcons[n]);
		F_ALLOC(Ucons[n]);

		F_ALLOC(Jx[n]);
		F_ALLOC(Jy[n]);

		C_ALLOC(src[n]);

	}
	C_ALLOC(x_cc);
	C_ALLOC(y_cc);
	F_ALLOC(x_fc);
	F_ALLOC(y_fc);

	C_ALLOC(cs);
	F_ALLOC(Lw);
	F_ALLOC(Uw);

	for (m = 0; m < NSCALAR; m++) {
		C_ALLOC(s[m]);
		C_ALLOC(s_gen[m]);
		F_ALLOC(Ls[m]);
		F_ALLOC(Us[m]);
		F_ALLOC(s_Jx[m]);
		F_ALLOC(s_Jy[m]);
		C_ALLOC(s_src[m]);
	}
#undef F_ALLOC
#undef C_ALLOC

	return g;
}

void free_grid(struct grid *g)
{
	int m, n;

#define FREE(q) free(g->q)
	for (n = 0; n < 4; n++) {
		FREE(cons[n]);
		FREE(prim[n]);
		FREE(cons_gen[n]);
		FREE(prim_gen[n]);

		FREE(Lprim[n]);
		FREE(Uprim[n]);
		FREE(Lcons[n]);
		FREE(Ucons[n]);

		FREE(Jx[n]);
		FREE(Jy[n]);

		FREE(src[n]);
	}
	FREE(x_cc);
	FREE(y_cc);
	FREE(x_fc);
	FREE(y_fc);

	FREE(cs);
	FREE(Lw);
	FREE(Uw);

	for (m = 0; m < NSCALAR; m++) {
		FREE(s[m]);
		FREE(s_gen[m]);
		FREE(Ls[m]);
		FREE(Us[m]);
		FREE(s_Jx[m]);
		FREE(s_Jy[m]);
		FREE(s_src[m]);
	}
#undef FREE

	free(g);
}
