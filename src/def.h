#ifndef DEF_H
#define DEF_H

#define HLLC 1
#define RECONSTRUCT 6
#define RECONSTRUCT_BOTH 1
#define RECONSTRUCT_FLATTEN 1
#define FANCY_POT_NRG 1
#define ROE_WAVESPEED 0
#define NTHREAD 4
#define THREAD_SCHEDULE static

#define NX 100
#define NY 300

#define BINARY 0
#define KH_INSTAB 0
#define RT_INSTAB 1
#define SOD_SHOCK 0
#define BLAST 0
#define IMPLOSION 0
#define STATIC_GRAV_TEST 0
#define SUPERSONIC 0
#define ATMOSPHERE 0
#define LINEAR_WAVE_TEST_X 0
#define LINEAR_WAVE_TEST_Y 0
#define LINEAR_WAVE_TEST_XY 0

#define NSCALAR 0
#define RHO_FLOOR 1e-8
#define PRESS_FLOOR 1e-10

#define MAX_OUT 40
#define MAX_EPOCH 1000000000000
#define CFL_NUM 0.43

#define BUF_LEN 1024

#define PI 3.1415926535897932384626433
#define BIG_G 6.67408e-8
#define kB 1.380649e-16

#define SIGN(x) (((x) > 0) - ((x) < 0))
#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define QUAD(x) ((x)*(x)*(x)*(x))

#define CEL(a,i,j) ((a)[(i)*ny + (j)])
#define FEL(a,i,j) ((a)[(i)*(ny+1) + (j)])

#if RECONSTRUCT != 1 && RECONSTRUCT != 2 && RECONSTRUCT != 3 && RECONSTRUCT != 6
#error "RECONSTRUCT must be in {1,2,3,6}"
#endif

void global_const();

#endif /* DEF_H */
