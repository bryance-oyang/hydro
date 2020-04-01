#ifndef WIND_TUNNEL_H
#define WIND_TUNNEL_H

#include "def.h"
#include "grid.h"
#include "boundary.h"
#include <math.h>

#define SPEED_LIM 1e2
#define PRESS_CEIL 25
#define WIND_SPEED 3

static inline double airfoil_shape(double x)
{
	if (x >= 0 && x <= 1) {
		return 0.5 * x * (2 - x);
	} else if (x >= 1 && x <= 5) {
		return -0.125 * (x - 5);
	} else {
		return 0;
	}
}

static inline int airfoil(double x, double y)
{
	/*
	if (y >= 0 && y <= airfoil_shape(x)) {
		return 1;
	} else {
		return 0;
	}
	*/

	double rx, ry, theta;

	theta = PI/6;
	rx = x*cos(theta) + y*sin(theta);
	ry = -x*sin(theta) + y*cos(theta);
	if (SQR(rx/2) + SQR(ry/1) <= 1) {
		return 1;
	} else {
		return 0;
	}
}

static inline void wind_tunnel_boundary(struct grid *g)
{
	int i, j;
	int nx, ny;

	nx = g->nx;
	ny = g->ny;

	inflow_boundary_left(g, 1, WIND_SPEED, 0, 1);

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
