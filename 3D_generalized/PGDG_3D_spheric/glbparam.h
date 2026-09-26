#ifndef GLBPARAM_H
#define GLBPARAM_H

#define SOLVER MPM
#define Particle Particle_base


#ifdef _DIM1_
#define DIM 1
#define DIM_X 1
#define DIM_Y 0
#define DIM_Z 0
#endif
#ifdef _DIM2_
#define DIM 2
#define DIM_X 1
#define DIM_Y 1
#define DIM_Z 0
#define ND_N 4 //mpm2025.10.1
#endif
#ifdef _DIM3_
#define DIM 3
#define DIM_X 1
#define DIM_Y 1
#define DIM_Z 1
#define ND_N 8 //mpm2025.10.1
#endif



#define CELL_RATIO  3.2  // 2026.7.3
#define CUT_OFF     2. //2021.11.11

#define SCALE_RATIO  2  // 

#define PI 3.1415926
#define GRAVITY  0

#define UNum 2 + DIM
#define BaseNum DIM
#define FldNum 2 //20260616

#endif
