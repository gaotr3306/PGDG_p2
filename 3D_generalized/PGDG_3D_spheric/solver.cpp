#include "solver.h"
#include "system.h"
#include "mpm.h"


/***************************************************/
/*                                                 */
/*       Functions defined in class "Solver"       */
/*                                                 */
/***************************************************/

//--------------------------------------------------
// Solver
//--------------------------------------------------
Solver::Solver(Initialization &Ini, communicator &world){
  timestep          = 0.0;
  integral_time     = 0.0;

  if (world.rank() == 0){
    cout<<"**********************************************************\n";
    cout<<"<<<<< Class Solver is initialized\n";
    cout<<"**********************************************************\n";
  }

}

//--------------------------------------------------
// MPM_Solver
//--------------------------------------------------
void Solver::MPM_Solver(Initialization &Ini, SOLVER &mpm, int n, communicator &world)
{
  integral_time = 0.0;
  while (integral_time < Ini.output_dt - 1.0e-10){

    mpm.Set_timestep(world);
    if(world.rank() == 0) cout<<"glbl_timestep="<<mpm.glbl_timestep<<endl;

    timestep = mpm.glbl_timestep;
	  int flag1=0;
    if (integral_time >= (Ini.output_dt - timestep)){
      timestep = Ini.output_dt - integral_time;
      mpm.glbl_timestep = timestep;
	    flag1=1;//
    }
    integral_time += timestep;
    mpm.run_time = run_time + integral_time;

    mpm.iterate_num ++;

    mpm.MPM_run_simulation (mpm.iterate_num, world, mpm.run_time, flag1);//


  }

  if (world.rank() == 0)
    cout<<"<<<<<<integral_time: "<<integral_time<<"\n";
  


  mpm.Output_plt_file(n, 1, world);



}
