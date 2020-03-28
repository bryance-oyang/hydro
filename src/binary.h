#ifndef BINARY_H
#define BINARY_H
#if BINARY == 1

#include <math.h>

#define BIN_OMEGA 1.556166490156496
#define GM1 1.027373356934356e5
#define GM2 1.027373356934356e4
#define BIN_SEP 3.600280088388794e1
#define BIN_COM 3.272981898535267

#define M1_CUTOFF 3
#define BIN_TEMP0 4
#define BIN_DISK_R0 12
#define BIN_DISK_CUTOFF_WIDTH 3
#define BIN_DISK_R1 18
#define SPEED_LIM 500

static inline double lin_interp(double x0, double x1, double y0, double y1, double x)
{
	return (y1 - y0) / (x1 - x0) * (x - x0) + y0;
}

static inline double alpha_temp(double r)
{
	return BIN_TEMP0 * pow(r, -3.0/4.0);
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

	rho_disk = 1;
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
	return bin_rho(r) * bin_temp(r);
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
	return sqrt(vel2) - r * BIN_OMEGA;
}

static inline void binary_boundary(struct grid *g, int step)
{
	int i, j;
	int nx, ny;
	double dt;

	nx = g->nx;
	ny = g->ny;

	if (step == 0) {
		dt = g->dt / 2;
	} else if (step == 1) {
		dt = g->dt;
	} else {
		dt = 0;
	}

#if _OPENMP
#pragma omp parallel for simd private(j) num_threads(NTHREAD) schedule(THREAD_SCHEDULE)
#endif /* _OPENMP */
	for (i = 2; i < nx-2; i++) {
		for (j = 2; j < ny-2; j++) {
			double x, y, r;

			x = CEL(g->x_cc,i,j);
			y = CEL(g->y_cc,i,j);
			r = sqrt(SQR(x) + SQR(y));

			if (r <= M1_CUTOFF) {
				CEL(g->prim[0],i,j) = RHO_FLOOR;
				CEL(g->prim[1],i,j) = 0;
				CEL(g->prim[2],i,j) = 0;
				CEL(g->prim[3],i,j) = PRESS_FLOOR;
			} else {
				double rho;

				rho = CEL(g->prim[0],i,j);
				CEL(g->prim[3],i,j) = rho * alpha_temp(r);
			}
		}
	}

#if _OPENMP
#pragma omp parallel for simd private(j) num_threads(NTHREAD) schedule(THREAD_SCHEDULE)
#endif /* _OPENMP */
	for (i = 2; i < nx-2; i++) {
		for (j = 2; j < ny-2; j++) {
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

#if _OPENMP
#pragma omp parallel for simd private(j) num_threads(NTHREAD) schedule(THREAD_SCHEDULE)
#endif /* _OPENMP */
	for (i = 2; i < nx-2; i++) {
		for (j = 2; j < ny-2; j++) {
			double vx, vy, v;

			vx = CEL(g->prim[1],i,j);
			vy = CEL(g->prim[2],i,j);
			v = sqrt(SQR(vx) + SQR(vy));
			if (v > SPEED_LIM) {
				CEL(g->prim[1],i,j) *= SPEED_LIM / v;
				CEL(g->prim[2],i,j) *= SPEED_LIM / v;
			}
		}
	}
}

#endif /* BINARY == 1 */
#endif /* BINARY_H */
