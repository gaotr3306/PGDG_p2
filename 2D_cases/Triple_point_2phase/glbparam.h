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



#define CELL_RATIO  3.2  // 2026.7.8
#define CUT_OFF     2. //2021.11.11

#define CR_2D  1.4 // 2023.8.1
#define SCALE_RATIO  2  // the ratio between two consecutive level //mpm2025.9.11

#define PI 3.1415926
#define GRAVITY  0
#define G_P  3
#define FRIC 1.

#define Total_L 0
#define UNum 4
#define BaseNum 2
#define FldNum 2 //20260616

#endif
