#ifndef LEVEL_INFO_H
#define LEVEL_INFO_H
#include <cmath>
#include "glbcls.h"
#include "glbfunc.h"
#include "Mypool.h"
#include "boost/multi_array.hpp"


using namespace std;
using namespace boost::mpi;
using namespace boost;

typedef boost::multi_array<DTAG, 3> array_tag;
typedef boost::multi_array<p_Cell_list, 3> array_cell_list;


class Level_info{
  int                total_num_particle;
  int                total_nun_color;

public:

// variables
  int                level;
  int                Lmin;
  int                Lmax;
  int                total_num_cell;
  int                num_leaf;
  int                num_leaf_particle;
  my_int             num_cell;            //number of cells in each directions
  my_int             glbl_cell_start;
  my_int             glbl_cell_end;
  my_int             cell_start;
  my_int             cell_end;

  my_real            dcell;               //cell size in each directions
  my_real            domain;              //domain size in each direction
  my_real            box_l;               //coordinate of the left corner of the domain
  my_real            box_r;
  Real               scale;

  Real               g_p_rt;//2026.1.15


  concurrent_vector  <p_Particle> leaf_particle;
  Mypool<Cell_list>  *p_cell_listpool; 
  // every level Cell_listpool point to the main memory pool
  array_cell_list    table_cell_list;
  //point to different cell list;
  Level_info(){};  

// functions
  void Initialize (int i, MPM *mpm, communicator &world);
  void Allocate_memory();
  void Initial_cell_list_and_reset_tags();
  void Reset_cell_list_info(communicator &world);
  void Update_cell_list(Particle *current_particle, communicator &world);
  void Refresh_neighbor_info(MPM *mpm, Particle *current, my_int cell_id, int status);//mpm2029.9.11
  void Interaction(Particle *current, Particle *neighbor,  Real dist, int status, MPM *mpm);//mpm2029.9.11
  void Refresh_ele_info(MPM *mpm, Particle *current, my_int cell_id);//2026.1.12
};
#endif
