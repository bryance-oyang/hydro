#ifndef DEF_H
#define DEF_H

#define HLLC 1
#define RECONSTRUCT 1
#define FANCY_POT_NRG 1
#define ROE_WAVESPEED 1
#define NTHREAD 4
#define THREAD_SCHEDULE static

#define NX 400
#define NY 400

#define BINARY 0
#define KH_INSTAB 0
#define RT_INSTAB 0
#define SOD_SHOCK 0
#define BLAST 0
#define IMPLOSION 1
#define STATIC_GRAV_TEST 0
#define SUPERSONIC 0
#define ATMOSPHERE 0

#define NSCALAR 0
#define RHO_FLOOR 1e-8
#define PRESS_FLOOR 1e-10

#define MAX_OUT 80
#define MAX_EPOCH 1000000000000
#define CFL_NUM 0.43

#define BUF_LEN 1024

#define PI 3.1415926535897932384626433
#define BIG_G 6.67408e-8
#define kB 1.380649e-16

#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define QUAD(x) ((x)*(x)*(x)*(x))

#define CEL(a,i,j) ((a)[(i)*ny + (j)])
#define FEL(a,i,j) ((a)[(i)*(ny+1) + (j)])

void global_const();

#endif /* DEF_H */
