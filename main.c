#include "hydro.h"
#include <stdio.h>
#include <stdlib.h>

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

	g = alloc_grid(NX, NY, DX, DY);
	init_grid(g);

	nout = 0;
	out_time = 0;
	for (epoch = 0; epoch < MAX_EPOCH; epoch++) {
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
