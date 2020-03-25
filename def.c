#include "def.h"

int NX;
int NY;
double XRANGE;
double YRANGE;
double DX;
double DY;
double OUT_TF;
double GAMMA;
double GRAV;

void global_const()
{
#if KH_INSTAB == 1
	NX = 512;
	NY = 512;
	XRANGE = 1;
	YRANGE = 1;
	DX = XRANGE / NX;
	DY = YRANGE / NY;
	OUT_TF = 1;
	GAMMA = 1.4;
	GRAV = 0;
#endif

#if RT_INSTAB == 1
	NX = 512;
	NY = 512;
	XRANGE = 1.5;
	YRANGE = 1.5;
	DX = XRANGE / NX;
	DY = YRANGE / NY;
	OUT_TF = 12.75;
	GAMMA = 1.4;
	GRAV = 0.1;
#endif

#if SOD_SHOCK == 1
	NX = 300;
	NY = 600;
	XRANGE = 1;
	YRANGE = 1.5;
	DX = XRANGE / NX;
	DY = YRANGE / NY;
	OUT_TF = 1;
	GAMMA = 5.0/3.0;
	GRAV = 0.1;
#endif
}
