#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <math.h>

static inline void periodic_boundary_left(struct grid *g)
{
	int i, j, k, n;
	int nx, ny;

	nx = g->nx;
	ny = g->ny;

	i = 0;
	for (n = 0; n < 4; n++) {
		for (k = 0; k < 3; k++) {
			for (j = 3; j < ny-3; j++) {
				CEL(g->prim[n],i+k,j) = CEL(g->prim[n],nx-6+k,j);
			}
		}
	}
}

static inline void periodic_boundary_right(struct grid *g)
{
	int i, j, k, n;
	int nx, ny;

	nx = g->nx;
	ny = g->ny;

	i = nx-1;
	for (n = 0; n < 4; n++) {
		for (k = 0; k < 3; k++) {
			for (j = 3; j < ny-3; j++) {
				CEL(g->prim[n],i-k,j) = CEL(g->prim[n],5-k,j);
			}
		}
	}
}

static inline void periodic_boundary_bot(struct grid *g)
{
	int i, j, k, n;
	int nx, ny;

	nx = g->nx;
	ny = g->ny;

	j = 0;
	for (n = 0; n < 4; n++) {
		for (i = 3; i < nx-3; i++) {
			for (k = 0; k < 3; k++) {
				CEL(g->prim[n],i,j+k) = CEL(g->prim[n],i,ny-6+k);
			}
		}
	}
}

static inline void periodic_boundary_top(struct grid *g)
{
	int i, j, k, n;
	int nx, ny;

	nx = g->nx;
	ny = g->ny;

	j = ny-1;
	for (n = 0; n < 4; n++) {
		for (i = 3; i < nx-3; i++) {
			for (k = 0; k < 3; k++) {
				CEL(g->prim[n],i,j-k) = CEL(g->prim[n],i,5-k);
			}
		}
	}
}

static inline void smooth_boundary_left(struct grid *g)
{
	int i, j, k, n;
	int ny;

	ny = g->ny;

	i = 0;
	for (n = 0; n < 4; n++) {
		for (k = 0; k < 3; k++) {
			for (j = 3; j < ny-3; j++) {
				CEL(g->prim[n],i+k,j) = CEL(g->prim[n],i+3,j);
			}
		}
	}
}

static inline void smooth_boundary_right(struct grid *g)
{
	int i, j, k, n;
	int nx, ny;

	nx = g->nx;
	ny = g->ny;

	i = nx-1;
	for (n = 0; n < 4; n++) {
		for (k = 0; k < 3; k++) {
			for (j = 3; j < ny-3; j++) {
				CEL(g->prim[n],i-k,j) = CEL(g->prim[n],i-3,j);
			}
		}
	}
}

static inline void smooth_boundary_bot(struct grid *g)
{
	int i, j, k, n;
	int nx, ny;

	nx = g->nx;
	ny = g->ny;

	j = 0;
	for (n = 0; n < 4; n++) {
		for (i = 3; i < nx-3; i++) {
			for (k = 0; k < 3; k++) {
				CEL(g->prim[n],i,j+k) = CEL(g->prim[n],i,j+3);
			}
		}
	}
}

static inline void smooth_boundary_top(struct grid *g)
{
	int i, j, k, n;
	int nx, ny;

	nx = g->nx;
	ny = g->ny;

	j = ny-1;
	for (n = 0; n < 4; n++) {
		for (i = 3; i < nx-3; i++) {
			for (k = 0; k < 3; k++) {
				CEL(g->prim[n],i,j-k) = CEL(g->prim[n],i,j-3);
			}
		}
	}
}

static inline void reflecting_boundary_left(struct grid *g)
{
	int i, j, k;
	int ny;

	ny = g->ny;

	i = 0;
	for (k = 0; k < 3; k++) {
		for (j = 3; j < ny-3; j++) {
			CEL(g->prim[0],i+k,j) = CEL(g->prim[0],i+5-k,j);
			CEL(g->prim[1],i+k,j) = -1 * CEL(g->prim[1],i+5-k,j);
			CEL(g->prim[2],i+k,j) = CEL(g->prim[2],i+5-k,j);
			CEL(g->prim[3],i+k,j) = CEL(g->prim[3],i+5-k,j);
		}
	}
}

static inline void reflecting_boundary_right(struct grid *g)
{
	int i, j, k;
	int nx, ny;

	nx = g->nx;
	ny = g->ny;

	i = nx-1;
	for (k = 0; k < 3; k++) {
		for (j = 3; j < ny-3; j++) {
			CEL(g->prim[0],i-k,j) = CEL(g->prim[0],i-5+k,j);
			CEL(g->prim[1],i-k,j) = -1 * CEL(g->prim[1],i-5+k,j);
			CEL(g->prim[2],i-k,j) = CEL(g->prim[2],i-5+k,j);
			CEL(g->prim[3],i-k,j) = CEL(g->prim[3],i-5+k,j);
		}
	}
}

static inline void reflecting_boundary_bot(struct grid *g)
{
	int i, j, k;
	int nx, ny;

	nx = g->nx;
	ny = g->ny;

	j = 0;
	for (i = 3; i < nx-3; i++) {
		for (k = 0; k < 3; k++) {
			CEL(g->prim[0],i,j+k) = CEL(g->prim[0],i,j+5-k);
			CEL(g->prim[1],i,j+k) = CEL(g->prim[1],i,j+5-k);
			CEL(g->prim[2],i,j+k) = -1 * CEL(g->prim[2],i,j+5-k);
			CEL(g->prim[3],i,j+k) = CEL(g->prim[3],i,j+5-k);
		}
	}
}

static inline void reflecting_boundary_top(struct grid *g)
{
	int i, j, k;
	int nx, ny;

	nx = g->nx;
	ny = g->ny;

	j = ny-1;
	for (i = 3; i < nx-3; i++) {
		for (k = 0; k < 3; k++) {
			CEL(g->prim[0],i,j-k) = CEL(g->prim[0],i,j-5+k);
			CEL(g->prim[1],i,j-k) = CEL(g->prim[1],i,j-5+k);
			CEL(g->prim[2],i,j-k) = -1 * CEL(g->prim[2],i,j-5+k);
			CEL(g->prim[3],i,j-k) = CEL(g->prim[3],i,j-5+k);
		}
	}
}

static inline void empty_boundary_left(struct grid *g)
{
	int i, j, k;
	int ny;

	ny = g->ny;

	i = 0;
	for (k = 0; k < 3; k++) {
		for (j = 3; j < ny-3; j++) {
			CEL(g->prim[0],i+k,j) = RHO_FLOOR;
			CEL(g->prim[1],i+k,j) = 0;
			CEL(g->prim[2],i+k,j) = 0;
			CEL(g->prim[3],i+k,j) = PRESS_FLOOR;
		}
	}
}

static inline void empty_boundary_right(struct grid *g)
{
	int i, j, k;
	int nx, ny;

	nx = g->nx;
	ny = g->ny;

	i = nx-1;
	for (k = 0; k < 3; k++) {
		for (j = 3; j < ny-3; j++) {
			CEL(g->prim[0],i-k,j) = RHO_FLOOR;
			CEL(g->prim[1],i-k,j) = 0;
			CEL(g->prim[2],i-k,j) = 0;
			CEL(g->prim[3],i-k,j) = PRESS_FLOOR;
		}
	}
}

static inline void empty_boundary_bot(struct grid *g)
{
	int i, j, k;
	int nx, ny;

	nx = g->nx;
	ny = g->ny;

	j = 0;
	for (i = 3; i < nx-3; i++) {
		for (k = 0; k < 3; k++) {
			CEL(g->prim[0],i,j+k) = RHO_FLOOR;
			CEL(g->prim[1],i,j+k) = 0;
			CEL(g->prim[2],i,j+k) = 0;
			CEL(g->prim[3],i,j+k) = PRESS_FLOOR;
		}
	}
}

static inline void empty_boundary_top(struct grid *g)
{
	int i, j, k;
	int nx, ny;

	nx = g->nx;
	ny = g->ny;

	j = ny-1;
	for (i = 3; i < nx-3; i++) {
		for (k = 0; k < 3; k++) {
			CEL(g->prim[0],i,j-k) = RHO_FLOOR;
			CEL(g->prim[1],i,j-k) = 0;
			CEL(g->prim[2],i,j-k) = 0;
			CEL(g->prim[3],i,j-k) = PRESS_FLOOR;
		}
	}
}

#endif /* BOUNDARY_H */
