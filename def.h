#ifndef DEF_H
#define DEF_H

#define HLLC 0

#define NX 10
#define NY 10
#define DX 1
#define DY 1

#define NSCALAR 1
#define RHO_FLOOR 1e-6
#define PRESS_FLOOR 1e-8

#define MAX_OUT 3
#define OUT_DT 0.1
#define MAX_EPOCH 3
#define CFL_NUM 0.36

#define BUF_LEN 1024

#define GAMMA ((double)5.0/3.0)

#define SQR(x) ((x)*(x))

#define CEL(a,i,j) ((a)[(i)*ny + (j)])
#define FEL(a,i,j) ((a)[(i)*(ny+1) + (j)])

#endif /* DEF_H */
