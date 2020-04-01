#ifndef INITIAL_COND_H
#define INITIAL_COND_H

#include "eos.h"
#include <stdlib.h>
#include <math.h>

#if WIND_TUNNEL == 1
#include "wind_tunnel.h"
#endif
#if BINARY == 1
#include "binary.h"
#endif

extern double XMIN;
extern double YMIN;
extern double GAMMA;
extern double GRAV;

static inline void init_prim(struct grid *g) {
	int i, j, nx, ny;
	double x, y;

	nx = g->nx;
	ny = g->ny;

	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			x = CEL(g->x_cc,i,j);
			y = CEL(g->y_cc,i,j);

			if (WIND_TUNNEL) {
#if WIND_TUNNEL == 1
				if (airfoil(x, y)) {
					CEL(g->prim[0],i,j) = 1e20;
					CEL(g->prim[1],i,j) = 0;
					CEL(g->prim[2],i,j) = 0;
					CEL(g->prim[3],i,j) = 1;
				} else {
					CEL(g->prim[0],i,j) = 1;
					CEL(g->prim[1],i,j) = 0;
					CEL(g->prim[2],i,j) = 0;
					CEL(g->prim[3],i,j) = 1;
				}
#endif
			}

			if (BINARY) {
#if BINARY == 1
				double r;

				r = sqrt(SQR(x) + SQR(y));
				CEL(g->prim[0],i,j) = bin_rho(r);
				CEL(g->prim[1],i,j) = -y/r * bin_vel(r);
				CEL(g->prim[2],i,j) = x/r * bin_vel(r);
				CEL(g->prim[3],i,j) = bin_press(r);
#endif
			}

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
				double pk2pk = 0.01;
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
				if (i < nx/2) {
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
				double r;
				r = sqrt(SQR(x) + SQR(y));
				if (r < 0.1) {
					CEL(g->prim[0],i,j) = 1;
					CEL(g->prim[1],i,j) = 0;
					CEL(g->prim[2],i,j) = 0;
					CEL(g->prim[3],i,j) = 10;
				} else {
					CEL(g->prim[0],i,j) = 1;
					CEL(g->prim[1],i,j) = 0;
					CEL(g->prim[2],i,j) = 0;
					CEL(g->prim[3],i,j) = 0.1;
				}
			}

			if (IMPLOSION) {
				if (x + y > 0.15) {
					CEL(g->prim[0],i,j) = 1;
					CEL(g->prim[1],i,j) = 0;
					CEL(g->prim[2],i,j) = 0;
					CEL(g->prim[3],i,j) = 1;
				} else {
					CEL(g->prim[0],i,j) = 0.125;
					CEL(g->prim[1],i,j) = 0;
					CEL(g->prim[2],i,j) = 0;
					CEL(g->prim[3],i,j) = 0.14;
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
					CEL(g->prim[1],i,j) = 5;
					CEL(g->prim[2],i,j) = -2;
					CEL(g->prim[3],i,j) = 2;
				} else {
					CEL(g->prim[0],i,j) = rho;
					CEL(g->prim[1],i,j) = 0;
					CEL(g->prim[2],i,j) = 0;
					CEL(g->prim[3],i,j) = 1 - GRAV * rho * y;
				}
			}

			if (ATMOSPHERE) {
				double rho0, rho, r, R0, ndens0, ndens, m0, nrg0;
				nrg0 = 2e18;
				R0 = 1e5;
				rho0 = 1.204e-3;
				ndens0 = 2.504e19;
				m0 = rho0 / ndens0;
				r = sqrt(SQR(x) + SQR(y - 4.2e5));
				ndens = ndens0 * exp(-m0 * GRAV * y / (kB * 270));
				rho = ndens * m0;
				CEL(g->prim[0],i,j) = rho;
				CEL(g->prim[1],i,j) = 0;
				CEL(g->prim[2],i,j) = 0;
				if (r < R0) {
					CEL(g->prim[3],i,j) = ndens * kB * 270 + (GAMMA - 1) * nrg0 / (PI * SQR(R0));
				} else {
					CEL(g->prim[3],i,j) = ndens * kB * 270;
				}
			}
		}
	}

	eos_prim_to_cons(g->prim, g->cons, nx, ny);
}

static inline void init_cons(struct grid *g) {
	int i, j, nx, ny;
	double x, y;

	nx = g->nx;
	ny = g->ny;

	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			x = CEL(g->x_cc,i,j);
			y = CEL(g->y_cc,i,j);

			if (LINEAR_WAVE_TEST_X || LINEAR_WAVE_TEST_Y || LINEAR_WAVE_TEST_XY) {
				double kx, ky, angle;
				double a, wave;
				double rho0, press0;

				a = 1e-6;
				if (LINEAR_WAVE_TEST_X) {
					angle = 0;
				} else if (LINEAR_WAVE_TEST_Y) {
					angle = PI/2;
				} else {
					angle = PI/3;
				}
				kx = 2*PI*cos(angle);
				ky = 2*PI*sin(angle);

				rho0 = 1;
				press0 = 1.0 / (GAMMA);
				wave = a * sin(kx*x + ky*y);

				CEL(g->cons[0],i,j) = rho0 + wave;
				CEL(g->cons[1],i,j) = 0 + cos(angle) * sqrt(GAMMA * press0 / rho0) * wave;
				CEL(g->cons[2],i,j) = 0 + sin(angle) * sqrt(GAMMA * press0 / rho0) * wave;
				CEL(g->cons[3],i,j) = press0/(GAMMA-1) + (GAMMA * press0)/((GAMMA-1) * rho0) * wave;
			}
		}
	}

	eos_cons_to_prim(g->cons, g->prim, nx, ny);
}

void init_grid(struct grid *g)
{
	int i, j, nx, ny;
	double dx, dy;

	nx = g->nx;
	ny = g->ny;
	dx = g->dx;
	dy = g->dy;

	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			CEL(g->x_cc,i,j) = XMIN + (i-2.5)*dx;
			CEL(g->y_cc,i,j) = YMIN + (j-2.5)*dy;
			FEL(g->x_fc,i,j) = XMIN + (i-3)*dx;
			FEL(g->y_fc,i,j) = YMIN + (j-3)*dy;
		}
	}

	if (LINEAR_WAVE_TEST_X || LINEAR_WAVE_TEST_Y || LINEAR_WAVE_TEST_XY) {
		init_cons(g);
	} else {
		init_prim(g);
	}
}


#endif /* INITIAL_COND_H */
