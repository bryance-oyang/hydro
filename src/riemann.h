#ifndef RIEMANN_H
#define RIEMANN_H

#include "grid.h"

void reconstruct(struct grid *g, int step, int dir);
void wavespeed(struct grid *g, int step, int dir);
void hlle(struct grid *g, int dir);
void hllc(struct grid *g, int dir);

#endif /* RIEMANN_H */
