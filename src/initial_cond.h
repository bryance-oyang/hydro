#ifndef INITIAL_COND_H
#define INITIAL_COND_H

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
		}
	}

	eos_prim_to_cons(g->prim, g->cons, nx, ny);
}

#endif /* INITIAL_COND_H */
