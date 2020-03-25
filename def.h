#ifndef DEF_H
#define DEF_H

#define HLLC 1
#define NTHREAD 4
#define THREAD_SCHEDULE static

#define KH_INSTAB 0
#define RT_INSTAB 0
#define SOD_SHOCK 0
#define BLAST 0
#define STATIC_GRAV_TEST 1

#define NSCALAR 1
#define RHO_FLOOR 1e-6
#define PRESS_FLOOR 1e-8

#define MAX_OUT 50
#define MAX_EPOCH 1000000000
#define CFL_NUM 0.43
#define OUT_DT ((double)(OUT_TF) / (MAX_OUT - 1))

#define BUF_LEN 1024

#define PI 3.1415926535897932384626433
#define SQR(x) ((x)*(x))

#define CEL(a,i,j) ((a)[(i)*ny + (j)])
#define FEL(a,i,j) ((a)[(i)*(ny+1) + (j)])

void global_const();

#endif /* DEF_H */
