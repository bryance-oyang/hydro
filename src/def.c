#include "def.h"

double XMIN;
double XMAX;
double YMIN;
double YMAX;
double XRANGE;
double YRANGE;
double DX;
double DY;
double OUT_TF;
double GAMMA;
double GRAV;

void global_const()
{
	XMIN = -0.5;
	XMAX = 0.5;
	YMIN = -0.5;
	YMAX = 0.5;
	OUT_TF = 12;
	GAMMA = 1.4;
	GRAV = 0;

#if KH_INSTAB == 1
	XMIN = -0.5;
	XMAX = 0.5;
	YMIN = -0.5;
	YMAX = 0.5;
	OUT_TF = 1;
	GAMMA = 1.4;
	GRAV = 0;
#endif

#if RT_INSTAB == 1
	XMIN = -0.25;
	XMAX = 0.25;
	YMIN = -0.75;
	YMAX = 0.75;
	OUT_TF = 12.75;
	GAMMA = 1.4;
	GRAV = 0.1;
#endif

#if SOD_SHOCK == 1
	XMIN = -0.5;
	XMAX = 0.5;
	YMIN = -0.75;
	YMAX = 0.75;
	OUT_TF = 1;
	GAMMA = 5.0/3.0;
	GRAV = 0.1;
#endif

#if BLAST == 1
	XMIN = -0.5;
	XMAX = 0.5;
	YMIN = -0.5;
	YMAX = 0.5;
	OUT_TF = 2;
	GAMMA = 1.4;
	GRAV = 0;
#endif

#if STATIC_GRAV_TEST == 1
	XMIN = -0.5;
	XMAX = 0.5;
	YMIN = -0.5;
	YMAX = 0.5;
	OUT_TF = 2;
	GAMMA = 1.4;
	GRAV = 0.1;
#endif

#if SUPERSONIC == 1
	XMIN = -0.5;
	XMAX = 0.5;
	YMIN = -0.5;
	YMAX = 0.5;
	OUT_TF = 10;
	GAMMA = 1.4;
	GRAV = 1;
#endif

	XRANGE = XMAX - XMIN;
	YRANGE = YMAX - YMIN;
	DX = XRANGE / (NX - 4);
	DY = YRANGE / (NY - 4);
}
