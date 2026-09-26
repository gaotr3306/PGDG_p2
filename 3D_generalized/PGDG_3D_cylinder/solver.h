#ifndef SOLVER_H
#define SOLVER_H
#include "glbcls.h"
#include <cmath>

using namespace std;
using namespace boost::mpi;

class Solver{

  Real     timestep;           //current timestep size
  Real     integral_time;       //integral time in current inteval

public:  
  Real     run_time;   //current iterate time
  Solver(Initialization &Ini, communicator &world);
  void MPM_Solver(Initialization &Ini, SOLVER &mpm, int n, communicator &world);//2021.9.6
};

#endif
