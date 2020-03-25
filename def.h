#ifndef DEF_H
#define DEF_H

#define NTHREAD 4
#define THREAD_SCHEDULE static

#define HLLC 0

#define NX 200
#define NY 600
#define XRANGE 0.5
#define YRANGE 1.5
#define DX ((double)(XRANGE) / (NX))
#define DY ((double)(YRANGE) / (NY))

#define GRAV 0.1
//#define GRAV 0

#define NSCALAR 1
#define RHO_FLOOR 1e-6
#define PRESS_FLOOR 1e-8

#define OUT_TF 12.75
//#define OUT_TF 1
#define MAX_OUT 40
#define MAX_EPOCH 1000000000
#define CFL_NUM 0.36

#define OUT_DT ((double)(OUT_TF) / (MAX_OUT - 1))

#define BUF_LEN 1024

#define GAMMA ((double)5.0/3.0)
//#define GAMMA 1.4

#define PI 3.1415926535897932384626433

#define SQR(x) ((x)*(x))

#define CEL(a,i,j) ((a)[(i)*ny + (j)])
#define FEL(a,i,j) ((a)[(i)*(ny+1) + (j)])

#endif /* DEF_H */
