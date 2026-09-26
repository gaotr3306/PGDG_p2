#ifndef SYSTEM_H
#define SYSTEM_H
#include "glbcls.h"

using namespace std;
using namespace boost::mpi;

// system parameters for global control
class Initialization{
  Real     start_time, end_time;
  int      n_thread;

public:
  Real     output_dt;
  int      num_output;           //num_output_file
  Initialization(communicator &world);
};

#endif
