#ifdef CDEBUG
#define _GNU_SOURCE
#include <fenv.h>
#endif /* CDEBUG */

#include "hydro.h"
#include "initial_cond.h"
#include <stdio.h>
#include <stdlib.h>

extern double XMIN;
extern double XMAX;
extern double YMIN;
extern double YMAX;
extern double DX;
extern double DY;
extern double OUT_TF;
extern double OUT_DT;

void output_array(char *filename, double *q, int nx, int ny)
{
	int i, j;
	FILE *f;

	if ((f = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "unable to open %s\n", filename);
		exit(EXIT_FAILURE);
	}

	for (i = 3; i < nx-3; i++) {
		for (j = 3; j < ny-3; j++) {
			fprintf(f, "%.16e ", CEL(q,i,j));
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

		fprintf(f, "%d\n", NX);
		fprintf(f, "%d\n", NY);
		fprintf(f, "%.16e\n", g->dx);
		fprintf(f, "%.16e\n", g->dy);
		fprintf(f, "%.16e\n", XMIN);
		fprintf(f, "%.16e\n", XMAX);
		fprintf(f, "%.16e\n", YMIN);
		fprintf(f, "%.16e\n", YMAX);
		fprintf(f, "%d\n", MAX_OUT);

		fclose(f);
	}
	sprintf(filename, "data/time.dat");
	if ((f = fopen(filename, "a")) == NULL) {
		fprintf(stderr, "unable to open %s\n", filename);
		exit(EXIT_FAILURE);
	}
	fprintf(f, "%.16e\n", g->time);

	sprintf(filename, "data/rho_%05d.dat", nout);
	output_array(filename, g->prim[0], g->nx, g->ny);

	sprintf(filename, "data/vx_%05d.dat", nout);
	output_array(filename, g->prim[1], g->nx, g->ny);

	sprintf(filename, "data/vy_%05d.dat", nout);
	output_array(filename, g->prim[2], g->nx, g->ny);

	sprintf(filename, "data/press_%05d.dat", nout);
	output_array(filename, g->prim[3], g->nx, g->ny);
}

void test_lin_wave(struct grid *g)
{
	int i, j, nx, ny, n;
	double dx;
	double *orig_cons[4];
	double error[4];
	char filename[BUF_LEN];
	FILE *f;

	nx = g->nx;
	ny = g->ny;
	dx = g->dx;

	for (n = 0; n < 4; n++) {
		orig_cons[n] = emalloc(nx*ny*sizeof(*orig_cons[n]));
	}
	lin_wave(g, orig_cons);

	sprintf(filename, "data/lin_wave.dat");
	if ((f = fopen(filename, "w")) == NULL) {
		fprintf(stderr, "unable to open %s\n", filename);
		exit(EXIT_FAILURE);
	}
	fprintf(f, "%d x %d\n", NX, NY);
	fprintf(f, "t = %.16e\n", g->time);
	for (n = 0; n < 4; n++) {
		error[n] = 0;
		for (i = 0; i < nx; i++) {
			j = ny/2;
			error[n] += dx*fabs(CEL(orig_cons[n],i,j) - CEL(g->cons[n],i,j));
		}
		fprintf(f, "%.16e\n", error[n]);
	}

	fclose(f);

	for (i = 0; i < 4; i++) {
		free(orig_cons[i]);
	}
}

int main()
{
	long epoch, nout;
	double out_time;
	struct grid *g;
	int lin_wave_tested;

#ifdef CDEBUG
	feenableexcept(FE_INVALID | FE_OVERFLOW | FE_DIVBYZERO);
#endif /* CDEBUG */

	global_const();

	g = alloc_grid(NX+6, NY+6, DX, DY);
	init_grid(g);
	boundary(g, -1);

	nout = 0;
	out_time = 0;
	lin_wave_tested = 0;
	for (epoch = 0; epoch < MAX_EPOCH; epoch++) {
		printf("t = %.3e\tdt = %.3e\t%.2f%%\n", g->time, g->dt, 100*g->time/OUT_TF);
		if (LINEAR_WAVE_TEST_X == 1 && lin_wave_tested == 0 && g->time >= 1) {
			test_lin_wave(g);
			lin_wave_tested = 1;
		}
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
