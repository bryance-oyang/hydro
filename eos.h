#ifndef EOS_H
#define EOS_H

#include "grid.h"

void eos_prim_floor(double **prim, int nx, int ny);
void eos_prim_to_cons(double **prim, double **cons, int nx, int
		ny);
void eos_cons_to_prim(double **cons, double **prim, int nx, int
		ny);
void eos_sound_speed(struct grid *g);

#endif /* EOS_H */
