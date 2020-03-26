#include "grid.h"
#include "emalloc.h"
#include "eos.h"
#include <stdlib.h>
#include <math.h>

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
