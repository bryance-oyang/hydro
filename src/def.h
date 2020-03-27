#ifndef DEF_H
#define DEF_H

#define HLLC 1
#define FANCY_POT_NRG 1
#define NTHREAD 4
#define THREAD_SCHEDULE static

#define NX 100
#define NY 300

#define BINARY 0
#define KH_INSTAB 0
#define RT_INSTAB 1
#define SOD_SHOCK 0
#define BLAST 0
#define STATIC_GRAV_TEST 0
#define SUPERSONIC 0
#define ATMOSPHERE 0

#define NSCALAR 1
#define RHO_FLOOR 1e-6
#define PRESS_FLOOR 1e-10
#define PRESS_RHO_CEIL 0

#define MAX_OUT 100
#define MAX_EPOCH 1000000000000
#define CFL_NUM 0.43

#define BUF_LEN 1024

#define BIN_OMEGA 1.556166490156496
#define GM1 1.027373356934356e5
#define GM2 1.027373356934356e4
#define BIN_SEP 3.600280088388794e1
#define BIN_COM 3.272981898535267
#define M1_CUTOFF 3

#define PI 3.1415926535897932384626433
#define kB 1.380649e-16
#define SQR(x) ((x)*(x))

#define CEL(a,i,j) ((a)[(i)*ny + (j)])
#define FEL(a,i,j) ((a)[(i)*(ny+1) + (j)])

void global_const();

#endif /* DEF_H */
