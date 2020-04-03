#ifndef BINARY_H
#define BINARY_H

#include "def.h"
#include "grid.h"
#include <math.h>

#define SATURN 0
#define AMCVN 1

#if SATURN == 1
#define GM1 (BIG_G * 5.683e29)
//#define GM2 (BIG_G * 3.7493e22)
#define GM2 (BIG_G * 3.7493e25)
#define BIN_SEP 1.85539e10
#define RGAS (kB / 5.01786569e-23)

#define BIN_ROT_FRAME 0
//#define M1_CUTOFF 5.8232e9
#define M1_CUTOFF 1e9
#define M2_CUTOFF 5e8
#define BIN_DISK_R0 7e9
#define BIN_DISK_CUTOFF_WIDTH 1e9
#define BIN_DISK_R1 1.5e10
#define BIN_RHO0 2
#define BIN_ISOTHERM 1
#define BIN_TEMP0 300
#define SPEED_LIM sqrt(2 * GM1 / M1_CUTOFF)
#endif

#if AMCVN == 1
#define GM1 1.027373356934356e5
//#define GM2 1.027373356934356e4
#define GM2 0.0
//#define BIN_SEP 36.00280088388794
#define BIN_SEP 100
#define RGAS 1

#define BIN_ROT_FRAME 0
#define M1_CUTOFF 1
#define M2_CUTOFF 0
#define BIN_DISK_R0 40
#define BIN_DISK_CUTOFF_WIDTH 10
#define BIN_DISK_R1 70
#define BIN_RHO0 2
#define BIN_ISOTHERM 1
#define BIN_TEMP0 0
#define SPEED_LIM 500
#endif

#define BIN_COM (((GM2) / ((GM1) + (GM2))) * (BIN_SEP))
#define BIN_OMEGA sqrt(((GM1) + (GM2)) / CUBE(BIN_SEP))
#define BIN_PERIOD ((2*PI) / BIN_OMEGA)
#define BIN_ANGLE 0

static inline double lin_interp(double x0, double x1, double y0, double y1, double x)
{
	return (y1 - y0) / (x1 - x0) * (x - x0) + y0;
}

static inline double alpha_temp(double r)
{
	if (BIN_ISOTHERM) {
		return BIN_TEMP0;
	} else {
		return BIN_TEMP0 * pow(r, -3.0/4.0);
	}
}

static inline double bin_temp(double r)
{
	double temp_bg;

	temp_bg = 0.1;
	if (r >= BIN_DISK_R0 - BIN_DISK_CUTOFF_WIDTH && r < BIN_DISK_R0) {
		return lin_interp(BIN_DISK_R0 - BIN_DISK_CUTOFF_WIDTH, BIN_DISK_R0, temp_bg, alpha_temp(BIN_DISK_R0), r);
	} else if (r >= BIN_DISK_R0 && r < BIN_DISK_R1) {
		return alpha_temp(r);
	} else if (r >= BIN_DISK_R1 && r < BIN_DISK_R1 + BIN_DISK_CUTOFF_WIDTH) {
		return lin_interp(BIN_DISK_R1, BIN_DISK_R1 + BIN_DISK_CUTOFF_WIDTH, alpha_temp(BIN_DISK_R1), temp_bg, r);
	} else {
		return temp_bg;
	}
}

static inline double bin_rho(double r)
{
	double rho_disk, rho_bg;

	rho_disk = BIN_RHO0;
	rho_bg = RHO_FLOOR;
	if (r >= BIN_DISK_R0 - BIN_DISK_CUTOFF_WIDTH && r < BIN_DISK_R0) {
		return lin_interp(BIN_DISK_R0 - BIN_DISK_CUTOFF_WIDTH, BIN_DISK_R0, rho_bg, rho_disk, r);
	} else if (r >= BIN_DISK_R0 && r < BIN_DISK_R1) {
		return rho_disk;
	} else if (r >= BIN_DISK_R1 && r < BIN_DISK_R1 + BIN_DISK_CUTOFF_WIDTH) {
		return lin_interp(BIN_DISK_R1, BIN_DISK_R1 + BIN_DISK_CUTOFF_WIDTH, rho_disk, rho_bg, r);
	} else {
		return rho_bg;
	}
}

static inline double bin_press(double r)
{
	return RGAS * bin_rho(r) * alpha_temp(r);
	/*
	double press_disk, press_bg;

	press_disk = 1;
	press_bg = 1e-3;
	if (r >= BIN_DISK_R0 - BIN_DISK_CUTOFF_WIDTH && r < BIN_DISK_R0) {
		return lin_interp(BIN_DISK_R0 - BIN_DISK_CUTOFF_WIDTH, BIN_DISK_R0, press_bg, press_disk, r);
	} else if (r >= BIN_DISK_R0 && r < BIN_DISK_R1) {
		return press_disk;
	} else if (r >= BIN_DISK_R1 && r < BIN_DISK_R1 + BIN_DISK_CUTOFF_WIDTH) {
		return lin_interp(BIN_DISK_R1, BIN_DISK_R1 + BIN_DISK_CUTOFF_WIDTH, press_disk, press_bg, r);
	} else {
		return press_bg;
	}
	*/
}

static inline double bin_dpdr(double r)
{
	return (bin_press(1.001*r) - bin_press(r)) / (0.001*r);
}

static inline double bin_vel(double r)
{
	double vel2;

	vel2 = fmax(0, GM1/r + r * bin_dpdr(r) / bin_rho(r));
	if (BIN_ROT_FRAME) {
		return sqrt(vel2) - r * BIN_OMEGA;
	} else {
		return sqrt(vel2);
	}
}

static inline void binary_boundary(struct grid *g, int step)
{
	int i, j;
	int nx, ny;
	double t, dt;

	nx = g->nx;
	ny = g->ny;

	if (step == 0) {
		t = g->time;
		dt = g->dt;
	} else if (step == 1) {
		t = g->time + g->dt;
		dt = g->dt / 2;
	} else {
		t = g->time + g->dt / 2;
		dt = g->dt;
	}

#if _OPENMP
#pragma omp parallel for simd private(j) num_threads(NTHREAD) schedule(THREAD_SCHEDULE)
#endif /* _OPENMP */
	for (i = 3; i < nx-3; i++) {
		for (j = 3; j < ny-3; j++) {
			double x, y, r, r2;

			x = CEL(g->x_cc,i,j);
			y = CEL(g->y_cc,i,j);
			r = sqrt(SQR(x) + SQR(y));
			if (BIN_ROT_FRAME) {
				r2 = sqrt(SQR(x - BIN_SEP*cos(BIN_ANGLE)) + SQR(y - BIN_SEP*sin(BIN_ANGLE)));
			} else {
				r2 = sqrt(SQR(x - BIN_SEP*cos(BIN_OMEGA*t + BIN_ANGLE)) + SQR(y - BIN_SEP*sin(BIN_OMEGA*t + BIN_ANGLE)));
			}

			if (r <= M1_CUTOFF || r2 <= M2_CUTOFF) {
				CEL(g->prim[0],i,j) = RHO_FLOOR;
				CEL(g->prim[1],i,j) = 0;
				CEL(g->prim[2],i,j) = 0;
				CEL(g->prim[3],i,j) = PRESS_FLOOR;
			} else {
				double rho;

				rho = CEL(g->prim[0],i,j);
				CEL(g->prim[3],i,j) = RGAS * rho * alpha_temp(r);
			}
		}
	}

	if (BIN_ROT_FRAME) {
#if _OPENMP
#pragma omp parallel for simd private(j) num_threads(NTHREAD) schedule(THREAD_SCHEDULE)
#endif /* _OPENMP */
		for (i = 3; i < nx-3; i++) {
			for (j = 3; j < ny-3; j++) {
				double vx, vy, dtomega;

				dtomega = dt * BIN_OMEGA;
				vx = CEL(g->prim[1],i,j);
				vy = CEL(g->prim[2],i,j);
				CEL(g->prim[1],i,j) = ((1 - SQR(dtomega))*vx + 2*dtomega*vy)
					/ (1 + SQR(dtomega));
				CEL(g->prim[2],i,j) = (-2*dtomega*vx + (1 - SQR(dtomega))*vy)
					/ (1 + SQR(dtomega));
			}
		}
	}

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
		}
	}
}

#endif /* BINARY_H */
