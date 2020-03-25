#ifdef CDEBUG
#define _GNU_SOURCE
#include <fenv.h>
#endif /* CDEBUG */

#include "hydro.h"
#include <stdio.h>
#include <stdlib.h>

extern int NX;
extern int NY;
extern double DX;
extern double DY;
extern double OUT_TF;

void output_array(char *filename, double *q, int nx, int ny)
{
	int i, j;
	FILE *f;

	if ((f = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "unable to open %s\n", filename);
		exit(EXIT_FAILURE);
	}

	for (i = 2; i < nx-2; i++) {
		for (j = 2; j < ny-2; j++) {
			fprintf(f, "%g ", CEL(q,i,j));
		}
		fprintf(f, "\n");
	}

	fclose(f);
}

void output(struct grid *g, int nout)
{
	char filename[BUF_LEN];
	FILE *f;

	if (nout == 0) {
		sprintf(filename, "data/param.dat");
		if ((f = fopen(filename, "w")) == NULL) {
			fprintf(stderr, "unable to open %s\n", filename);
			exit(EXIT_FAILURE);
		}

		fprintf(f, "%d\n", g->nx-4);
		fprintf(f, "%d\n", g->ny-4);
		fprintf(f, "%g\n", g->dx);
		fprintf(f, "%g\n", g->dy);

		fclose(f);
	}

	sprintf(filename, "data/rho_%05d.dat", nout);
	output_array(filename, g->prim[0], g->nx, g->ny);

	sprintf(filename, "data/press_%05d.dat", nout);
	output_array(filename, g->prim[3], g->nx, g->ny);
}

int main()
{
	long epoch, nout;
	double out_time;
	struct grid *g;

#ifdef CDEBUG
	feenableexcept(FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO);
#endif /* CDEBUG */

	global_const();

	g = alloc_grid(NX, NY, DX, DY);
	init_grid(g);
	boundary(g);

	nout = 0;
	out_time = 0;
	for (epoch = 0; epoch < MAX_EPOCH; epoch++) {
		printf("t = %g | dt = %g\n", g->time, g->dt);
		if (g->time >= out_time) {
			output(g, nout);
			out_time = g->time + OUT_DT;
			nout++;
			if (nout >= MAX_OUT) {
				break;
			}
		}

		advance_timestep(g);
	}

	free_grid(g);
	return 0;
}
