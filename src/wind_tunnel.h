#ifndef WIND_TUNNEL_H
#define WIND_TUNNEL_H

#include "def.h"
#include "grid.h"
#include "boundary.h"
#include <math.h>

#define SPEED_LIM 1e2
#define PRESS_CEIL 1e2

static inline void wind_tunnel_boundary(struct grid *g)
{
	int i, j;
	int nx, ny;

	nx = g->nx;
	ny = g->ny;

	inflow_boundary_left(g, 1, 3, 0, 1);

#if _OPENMP
#pragma omp parallel for simd private(j) num_threads(NTHREAD) schedule(THREAD_SCHEDULE)
#endif /* _OPENMP */
	for (i = 3; i < nx-3; i++) {
		for (j = 3; j < ny-3; j++) {
			double vx, vy, v;

			vx = CEL(g->prim[1],i,j);
			vy = CEL(g->prim[2],i,j);
			v = sqrt(SQR(vx) + SQR(vy));
			if (SPEED_LIM > 0) {
				if (v > SPEED_LIM) {
					CEL(g->prim[1],i,j) *= SPEED_LIM / v;
					CEL(g->prim[2],i,j) *= SPEED_LIM / v;
				}
			}

			if (CEL(g->prim[3],i,j) > PRESS_CEIL) {
				CEL(g->prim[3],i,j) = PRESS_CEIL;
			}
		}
	}
}

#endif
