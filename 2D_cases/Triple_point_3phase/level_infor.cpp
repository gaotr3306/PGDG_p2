#include "level_infor.h"
#include "mpm.h"
#include "particle.h"

#include "cell_list.h"
#include "boost/multi_array.hpp"
typedef boost::multi_array_types::extent_range mrange;

/***************************************************/
/*                                                 */
/*     Functions defined in class "Level_info"     */
/*                                                 */
/***************************************************/

//-------------------------------------------------------
// initialze all the necessary parameters for level_info
//-------------------------------------------------------
void Level_info::Initialize
(int i, MPM *mpm, communicator &world)
{
  level              = i;
  Lmin               = mpm->Lmin;
  Lmax               = mpm->Lmax;
  total_num_particle = mpm->glbl_total_num_particle;
  total_nun_color    = mpm->glbl_total_num_color;
  p_cell_listpool    = &(mpm->cell_listpool);
  g_p_rt = mpm->p_g/mpm->p_sz;//2026.1.15

  num_leaf_particle  = 0;
  leaf_particle.clear();
  Real multi = powern(Real(SCALE_RATIO),level);
  Real box_size = mpm->domain.i;
#if DIM_Y == 1
  box_size = AMIN1(mpm->domain.j, box_size);
#endif
#if DIM_Z == 1
  box_size = AMIN1(mpm->domain.k, box_size);
#endif
  if ( mpm->ini_scale/multi > box_size){
    cout<<"ERROR: Current level cannot be initialized!!\n";
    exit(0);
  }
  my_set_data  (domain, mpm->domain);
  my_set_data  (box_l, mpm->box_l);
  my_set_data  (box_r, mpm->box_r);
  my_set_const (num_cell, 1);
  my_set_data  (num_cell, my_multiply_const(mpm->ini_num_cell, multi));
  my_self_multiply (num_cell, total_num_cell);
  dcell = my_devide_data (mpm->domain, num_cell);
  scale = AMAX1(DIM_X*dcell.i, AMAX1(DIM_Y*dcell.j, DIM_Z*dcell.k));

  cell_start.i = cell_start.j = cell_start.k = 0;
  cell_end.i = cell_end.j = cell_end.k = 1;
  my_set_data (cell_end, num_cell); 


  glbl_cell_start.i = cell_start.i;
  glbl_cell_start.j = cell_start.j;
  glbl_cell_start.k = cell_start.k;
  glbl_cell_end.i = cell_end.i;
  glbl_cell_end.j = cell_end.j;
  glbl_cell_end.k = cell_end.k;


  Allocate_memory();

  Initial_cell_list_and_reset_tags(); // reset tags to zero

    
  if (world.rank() == 0){
    cout<<"**********************************************************\n";
    cout<<"<<<<< Level: "<<level<<"I is:"<<i<<" is initialized\n";
    cout<<"<<<<< Total number of cell: "<<total_num_cell<<"\n";
    cout<<"<<<<< Scale size          : "<<scale<<"\n";
    cout<<"<<<<< cell start          : "<<cell_start.i<<" "<<cell_start.j<<" "<<cell_start.k<<"\n";
    cout<<"<<<<< cell end            : "<<cell_end.i<<" "<<cell_end.j<<" "<<cell_end.k<<"\n";//2024.11.30
    cout<<"**********************************************************\n";
  }

}

//-------------------------------------------------------
// Allocate_memory
//-------------------------------------------------------
void Level_info::Allocate_memory()
{

  table_cell_list.resize
    (boost::extents[mrange(cell_start.i,cell_end.i)][mrange(cell_start.j,cell_end.j)][mrange(cell_start.k,cell_end.k)]);


}
//-------------------------------------------------------
// reset tags to 0 for initialization
//-------------------------------------------------------
void Level_info::Initial_cell_list_and_reset_tags()
{
  static affinity_partitioner ap;
  parallel_for( blocked_range3d<int>(cell_start.i, cell_end.i, cell_start.j, cell_end.j, cell_start.k, cell_end.k),
                [&](const blocked_range3d<int>& r){
    for(int i=r.pages().begin(); i!=r.pages().end(); ++i){
      for(int j=r.rows().begin(); j!=r.rows().end(); ++j){
        for(int k=r.cols().begin(); k!=r.cols().end(); ++k){

          if (NULL == table_cell_list[i][j][k]){
            table_cell_list[i][j][k]   = p_cell_listpool->malloc();
          }
          table_cell_list[i][j][k]->Initialize(this);
        }
      }
    }
  }, ap);
}
//-------------------------------------------------------
// reset cell list info
//-------------------------------------------------------
void Level_info::Reset_cell_list_info(communicator &world)
{
  my_int start;
  my_int end;
  start.i = cell_start.i;
  start.j = cell_start.j;
  start.k = cell_start.k;
  end.i   = cell_end.i;
  end.j   = cell_end.j;
  end.k   = cell_end.k;

  static affinity_partitioner ap;
  parallel_for( blocked_range3d<int>(start.i, end.i, start.j, end.j, start.k, end.k),
           [&](const blocked_range3d<int>& r){
    for(int i = r.pages().begin(); i < r.pages().end(); i++)
      for(int j = r.rows().begin(); j < r.rows().end(); j++)
        for(int k = r.cols().begin(); k < r.cols().end(); k++)
          table_cell_list[i][j][k]->Reset_tags();
  }, ap);
}
//-------------------------------------------------------
// Update the cell_list in current level
//-------------------------------------------------------
void Level_info::Update_cell_list(Particle *current_particle, communicator &world)
{
  my_real coord_shift = my_minus_data (matrx2my_real(current_particle->coord), box_l);
  my_int  pos = get_cell_id (coord_shift, dcell, cell_start, cell_end);

  table_cell_list[pos.i][pos.j][pos.k]->Add_particle(current_particle);
  current_particle->level = level;
}


//-------------------------------------------------------
// refresh the neighbor information for each particle
//-------------------------------------------------------
void Level_info::Refresh_neighbor_info(MPM *mpm, Particle *current, my_int cell_id, int status)
{
  int radius = ceil(CELL_RATIO/g_p_rt);

  int     i_c = DIM_X==1 ? radius : 0;
  int     j_c = DIM_Y==1 ? radius : 0;
  int     k_c = DIM_Z==1 ? radius : 0;

  if(status ==3) //2026.2.7
  {
    i_c = DIM_X==1 ? 2 : 0;
    j_c = DIM_Y==1 ? 2 : 0;
    k_c = DIM_Z==1 ? 2 : 0;
  }


  for(int t=AMAX1((cell_id.i-i_c),cell_start.i); t<=AMIN1((cell_id.i+i_c),(cell_end.i-1)); t++){
    for(int s=AMAX1((cell_id.j-j_c),cell_start.j); s<=AMIN1((cell_id.j+j_c),(cell_end.j-1)); s++){
      for(int m=AMAX1((cell_id.k-k_c),cell_start.k); m<=AMIN1((cell_id.k+k_c),(cell_end.k-1)); m++){
        
          for(int jcyc=0;jcyc<int(table_cell_list[t][s][m]->particle_list.size());jcyc++){
            Particle *neighbor = table_cell_list[t][s][m]->particle_list[jcyc];

            my_real dr = matrx2my_real(current->coord -neighbor->coord);
            Real  dist = get_distance(dr);

            if(status ==2)
            {
              if(dist <= current->h*CUT_OFF) //2025.12.25
                current->neighbor.push_back(neighbor);
            }
            else if(status ==3)
            {
              if(dist <= mpm->p_sz* 2.) //2026.2.7
                current->neighbor.push_back(neighbor);
            }
  
          }
        
      }
    }
  }
}
//-------------------------------------------------------
// Refresh_ele_info
//-------------------------------------------------------
void Level_info::Refresh_ele_info(MPM *mpm, Particle *current, my_int cell_id)
{
  int     t = DIM_X==1 ? cell_id.i : 0;
  int     s = DIM_Y==1 ? cell_id.j : 0;
  int     m = DIM_Z==1 ? cell_id.k : 0;

  for(int jcyc=0;jcyc<int(table_cell_list[t][s][m]->particle_list.size());jcyc++){
    Particle *neighbor = table_cell_list[t][s][m]->particle_list[jcyc];

    if( neighbor->phase !=10 )
    {
      current->neighbor.push_back(neighbor) ;//
      neighbor->in_ele.push_back( current ) ; //2025.12.25
    }
  }

}

//-------------------------------------------------------
// interaction between two neighbor particle 
// for the neighbor information
//-------------------------------------------------------
void Level_info::Interaction
(Particle *current, Particle *neighbor, Real dist, int status, MPM *mpm) 
{
//mpm2025.9.11
  if(status ==2)
  {
    if(dist <= current->h*CUT_OFF) //2025.12.25
      current->neighbor.push_back(neighbor);
  }

}


