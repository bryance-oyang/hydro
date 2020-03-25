#ifndef DEF_H
#define DEF_H

#define NTHREAD 4
#define THREAD_SCHEDULE static

#define HLLC 0

#define NX 200
#define NY 200
#define DX 0.1
#define DY 0.1

#define NSCALAR 1
#define RHO_FLOOR 1e-6
#define PRESS_FLOOR 1e-8

#define MAX_OUT 10
#define OUT_DT 1.0
#define MAX_EPOCH 1000000000
#define CFL_NUM 0.36

#define BUF_LEN 1024

#define GAMMA ((double)5.0/3.0)

#define SQR(x) ((x)*(x))

#define CEL(a,i,j) ((a)[(i)*ny + (j)])
#define FEL(a,i,j) ((a)[(i)*(ny+1) + (j)])

#endif /* DEF_H */
