#include <cmath>
#include "system.h"

/***************************************************/
/*                                                 */
/*   Functions defined in class "Initialization"   */
/*                                                 */
/***************************************************/
//--------------------------------------------------
// Initialization
//--------------------------------------------------
Initialization::Initialization(communicator &world){

  start_time = 0.0;
  end_time   = 1100e-6 ;
  output_dt  = 5e-6;
  num_output = int (ceil((end_time - start_time)/ output_dt));

  if (world.rank() == 0){
    cout<<"**********************************************************\n";
    cout<<"<<<<< System parameters are set\n";
    cout<<"<<<<< start time:             "<<start_time<<"\n";
    cout<<"<<<<< end time:               "<<end_time<<"\n";
    cout<<"<<<<< output dt:              "<<output_dt<<"\n";
    cout<<"<<<<< number of output file:  "<<num_output<<"\n";
    cout<<"***********************************************************\n";
  }

}