#include <cmath>
#include "glbcls.h"
#include "mpm.h"

#include "system.h"
#include "solver.h"
#include <Eigen/Core> 

using namespace std;
using namespace tbb;
using namespace boost;
using namespace boost::mpi;
/***************************************************/
/*           Main            */
/***************************************************/

int main(int argc, char **argv)
{
  environment env(argc,argv);
  communicator world;

  Initialization Ini(world);

  task_arena init(128);//mpm2025.9.11

  Eigen::setNbThreads(16);
  initParallel();

  SOLVER mpm(Ini, world);
  mpm.Initialize_case(Ini, world);
  Solver solver(Ini, world);

  Real Time = mpm.run_time; 
  if (world.rank() == 0){
    cout<<"**********************************************************\n";
    cout<<"<<<<< Simulation starts at time: "<< Time<< "\n";
  }

  mpm.Output_plt_file(0, 0, world);
  mpm.time_total.restart();


  for(int n= 1; n<=Ini.num_output; n++){ 
    solver.run_time = Time;

    if (world.rank() == 0) cout<<"Time:  "<<Time<<"\n";

    solver.MPM_Solver(Ini, mpm, n, world);
    Time += Ini.output_dt;  
  }

  mpm.time_for_total = mpm.time_total.elapsed();

  if (world.rank() == 0){
    cout<<"**********************************************************\n";
    cout<<"Time for Simulation:       "<<mpm.time_for_simulation<<"\n";
    cout<<"Time for Mapping:          "<<mpm.time_for_mapping<<"\n";
    cout<<"Total time:                "<<mpm.time_for_total<<"\n";
    cout<<"<<<<< Simulation finished. \n";
    cout<<"**********************************************************\n";
  }
  return 0;
}
