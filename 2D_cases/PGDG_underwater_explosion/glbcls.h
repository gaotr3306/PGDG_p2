#ifndef GLBCLS_H
#define GLBCLS_H

//-----------------------------------------
//  Self-defined data types
//----------------------------------------
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <stdio.h>//2025.1.14
#include <algorithm>
#include <string.h>
#include <fstream>
#include <vector>
#include "glbparam.h"
#include <omp.h> //2025.1.14
#include <unordered_set> //2026.9.21

#include "tbb/tick_count.h"
// #include "tbb/atomic.h"
#include "tbb/concurrent_vector.h"
#include "tbb/concurrent_unordered_set.h"
// #include "tbb/task_scheduler_init.h"
 #include "tbb/task_arena.h" //mpm2025.9.11
#include "tbb/blocked_range.h"
#include "tbb/blocked_range2d.h"
#include "tbb/blocked_range3d.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_reduce.h"
#include "tbb/scalable_allocator.h"
#include <boost/mpi.hpp>
#include <boost/config.hpp>
#include <boost/multi_array.hpp>



#include <Eigen/Dense> //2021.11.25
#include <Eigen/Eigenvalues> //2021.11.25
#include <Eigen/Sparse> ////2025.1.13
#include <Eigen/IterativeLinearSolvers> //2025.1.13
#include <Eigen/SparseCholesky> //2025.1.13



using namespace std;
using namespace Eigen;//2021.11.25
using namespace tbb;
using namespace boost;
using namespace boost::serialization;

class   Particle_base;
class   Cell_list;
class   Level_info;
class   MPM;
class   Solver;
class   Initialization;
typedef Cell_list* p_Cell_list;
typedef Level_info* p_Level_info;
typedef char DTAG;  //define DTAG as data type to tag cells;


typedef Particle_base* p_Particle;




#ifdef _DOUBLE_
typedef double Real;
#endif

#ifdef _SINGLE_
typedef float  Real;
#endif

class my_int{
  private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & i;
    ar & j;
    ar & k;
  }
  public:
  int i; int j; int k;
  my_int(){i = j = k = 0;};
};

class my_real{
  private:
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & i;
    ar & j;
    ar & k;
  }
  public:
  Real i; Real j; Real k;
  my_real(){i = j = k = 0.;};
};



#endif
