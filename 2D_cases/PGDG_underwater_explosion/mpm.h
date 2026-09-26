#ifndef MPM_H
#define MPM_H
#include "glbcls.h"
#include "Mypool.h"
#include "particle.h"

#include "cell_list.h"
#include "level_infor.h"
#include "system.h"


// #include "voro++.hh"
// using namespace voro; //2022.4.29


using namespace std;
using namespace boost;
using namespace boost::serialization;
using namespace boost::mpi;
using namespace tbb;

#include "boost/polygon/voronoi.hpp"
#include "boost/polygon/voronoi_diagram.hpp"
#include "boost/polygon/voronoi_builder.hpp"

using boost::polygon::voronoi_builder;
using boost::polygon::voronoi_diagram;
using boost::polygon::default_voronoi_builder;
using boost::polygon::x;
using boost::polygon::y;
using boost::polygon::low;
using boost::polygon::high;//2022.5.2



class Cell_dg{
public:

  concurrent_vector <Vector<Real,3>> face_center;
  concurrent_vector <Vector<Real,3>> face_normal;
  concurrent_vector <Vector<Real,3>> face_tan;
  concurrent_vector <Real> face_S;
  concurrent_vector <my_int> cell_id;
  concurrent_vector <int> bound;//2026.7.13

  Real Vol ;
  Vector<Real,3> cell_center;
};

  

//-----------------------------------------
class MPM{

public:
  int      glbl_total_num_particle;

  int      glbl_total_num_color;
  int      num_level;
  int      Lmax;
  int      Lmin;
  int      need_for_refresh;
  my_int   glbl_num_particle;
  Real     glbl_timestep;      //current timestep size
  
  Real     glbl_max_v;
  vector<int> nparticle_in_each_cpu;
  vector<int> refndnum_in_each_cpu;//2022.1.21

  int      iterate_num;        //current iterate number
  Real     local_timestep;     //current timestep size
  Real     run_time;           //current iterate time
  Real     timestep_shift;     //timestep size for redistributing particles



  Real ini_rho[2];
  Real c[2]; //2026.1.20
  Real mu[2];
  Real gamma[2];
  Real Pp0[2]; //2026.1.27
  Real Pb0 =0;//20260528


// globle variables
  Real     ini_h;       //initial smooth length
  Real     ini_scale;   //initial scale for base level
  Real     max_v;
  my_int   ini_num_cell;//number of cells at beginning
  my_real  domain;      //domain size in each direction
  my_real  box_l;       //coordinate of the left corner of the domain
  my_real  box_r;
  my_real  particle_box;
  my_real  particle_box_l;
  my_real  particle_box_r;



// local variables
  int      total_num_particle;

  p_Particle          *particle;


  Mypool<Particle>    particlepool;
  Mypool<Cell_list>   cell_listpool;
  p_Level_info        *level_info;

  vector<Real>  Gp ;
  vector<Real>  Gw ;
  vector<Real>  Gpf ;
  vector<Real>  Gwf ;//2026.7.8

  Real     time_for_simulation;
  Real     time_for_reset_neighbor;
  Real     time_for_mapping;
  Real     time_for_bulid_local_map;
  Real     time_for_map_particle_to_tree;
  Real     time_for_refresh_neighbor;
  Real     time_for_clear_cell_list;
  Real     time_for_update_every_level_info;
  Real     time_for_total;
  timer    time_simulation;
  timer    time_reset_neighbor;
  timer    time_mapping;
  timer    time_bulid_local_map;
  timer    time_map_particle_to_tree;
  timer    time_refresh_neighbor;
  timer    time_clear_cell_list;
  timer    time_update_every_level_info;
  timer    time_total;



  concurrent_vector  <p_Particle>  particle_fluid;//2025.12.25
  concurrent_vector  <p_Particle>  particle_bound;//2026.1.1


  p_Particle **Ele_cen;////2026.1.12
  Cell_dg **cell_dg;
  my_int NUM;
  my_real  p_size;      //initial_particle_size

  Real p_sz;//2026.1.12
  Real p_g;//2026.1.12

  concurrent_vector  <p_Particle> Used_Ele_cen;//2026.1.12


  int Iden=0;


// functions
  MPM(Initialization &Ini, communicator &world);

  void Initialize_case(Initialization &Ini, communicator &world);
  void MPM_run_simulation (int n, communicator &world, Real integral_time, int flag1);

  void Output_plt_file(int n, int flag, communicator &world);
  void Output_point(communicator &world, Real time1);//



  void Set_timestep (communicator &world);
  void Reset_cell_list_info(communicator &world);
  void Map_the_particle_to_tree(communicator &world);
  void Reset_particle_neighbor_info(communicator &world, int flag);//2026.7.8
  void Refresh_neighbor_info(int flag, communicator &world);

  void Used_state_update(communicator &world);//2026.1.12

  void MLS_P2G(communicator &world);//2025.12.25
  void Calculate_limiter(communicator &world);//2025.12.25
  void Grid_acceleration(communicator &world);//2025.12.25
  void Update_particle(communicator &world, int md);//mpm2025.9.11
  void Update_grid(communicator &world, int md);//2025.12.25


  Real EOS_P(Real r, int s);
  Real EOS_rhoe(Real p1, int s);//2026.1.20

  void Split_particle(communicator &world);
  void Merge_particle(communicator &world);//2026.2.7
  void Reset_refresh_merge_neighbor(communicator &world);//2026.2.
  Vector<Real, 3> MGFM(Real ul, Real ur, Real pl, Real pr, Real rhol, Real rhor, int s);//20260611

  void Split_particle_fast(communicator &world);
  void Merge_particle_fast(communicator &world);//2026.9.21
};
#endif
