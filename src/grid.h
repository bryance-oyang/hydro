#ifndef GRID_H
#define GRID_H

#include "def.h"

struct grid {
	double *cons[4];
	double *prim[4];
	double *cons_gen[4];
	double *prim_gen[4];
	double *cs;

	double *Lprim[4];
	double *Uprim[4];
	double *Lcons[4];
	double *Ucons[4];

	double *Lw;
	double *Uw;

	double *Jx[4];
	double *Jy[4];
	double *src[4];

	double *s[NSCALAR];
	double *s_gen[NSCALAR];
	double *Ls[NSCALAR];
	double *Us[NSCALAR];
	double *s_Jx[NSCALAR];
	double *s_Jy[NSCALAR];
	double *s_src[NSCALAR];

	double *x_cc;
	double *y_cc;
	double *x_fc;
	double *y_fc;

	double time;
	double dt;
	int nx;
	int ny;
	double dx;
	double dy;
};

void init_grid(struct grid *g);
struct grid *alloc_grid(int nx, int ny, double dx, double dy);
void free_grid(struct grid *g);

#endif /* GRID_H */
