#ifndef HYDRO_H
#define HYDRO_H

#include "grid.h"

void boundary(struct grid *g, int step);
void advance_timestep(struct grid *g);

#endif /* HYDRO_H */
