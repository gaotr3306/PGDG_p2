#include <cmath>
#include "glbfunc.h"
#include "mpm.h"
#include <algorithm> //2025.1.14

/***************************************************/
/*                                                 */
/*           Functions defined in class "MPM"      */
/*                                                 */
/***************************************************/


//--------------------------------------------------
// MPM incompressible class initialization
//--------------------------------------------------
MPM::MPM(Initialization &Ini, communicator &world)
{
  time_for_total                   = 0.;
  time_for_simulation              = 0.;
  time_for_reset_neighbor          = 0.;
  time_for_mapping                 = 0.;
  time_for_bulid_local_map         = 0.;
  time_for_map_particle_to_tree    = 0.;
  time_for_refresh_neighbor        = 0.;
  time_for_clear_cell_list         = 0.;
  time_for_update_every_level_info = 0.;


}
//--------------------------------------------------
// initial condition for cases
//--------------------------------------------------
void MPM::Initialize_case(Initialization &Ini, communicator &world)//2023.12.22
{


  Gp = {-1./sqrt(3.), 1./sqrt(3.)}; 
  Gw = {0.5, 0.5 }; 


  Gpf = {-1./sqrt(3.), 1./sqrt(3.)}; 
  Gwf = {0.5, 0.5  }; 
  

  Real Rt = 1 ; //
  my_real water_domain;
  water_domain.i = 4 * Rt;
  water_domain.j = 1.5* Rt;
  
  domain.i = particle_box.i = 7. *Rt ;  //5.366
  domain.j = particle_box.j = 3 *Rt;
  domain.k = particle_box.k = 1;

  num_level =  1 ; 

  p_sz = domain.j/200 * Rt ;//

  p_g = domain.j/200 * Rt ;//2026.1.12


  box_l.i  = particle_box_l.i = 0;
  box_l.j  = particle_box_l.j = 0;
  box_l.k  = particle_box_l.k = 0;

  box_r.i  = particle_box_r.i = domain.i ;
  box_r.j  = particle_box_r.j = domain.j ;
  box_r.k  = particle_box_r.k = domain.k ;

  int Nx = domain.i/p_g; Real px = domain.i/Real(Nx);
  int Ny = domain.j/p_g; Real py = domain.j/Real(Ny);
  NUM.i = Nx; NUM.j = Ny; 
  p_size.i = px;   p_size.j = py; 

  p_sz = pow(px* py, 1./DIM)/round(p_g/p_sz) ;//
  ini_h          = CELL_RATIO/CUT_OFF* p_sz;

  cout<<px<<", "<<py<<endl;

  
  ini_rho[0] = 1. ;//
  Pp0[0] = 0;
  gamma[0] = 1.4; 
  mu[0]      = 0 ; // 1e-3 
  c[0]      =  sqrt(gamma[0]* (Pp0[0] + 0)/ini_rho[0]) ;// 

  ini_rho[1] =  0.125 ;// 1.225
  Pp0[1] = 0;
  gamma[1] = 1.5;
  mu[1]      = 0  ; // 1.8e-5 
  c[1]       =  sqrt(gamma[1]* (Pp0[1] + 0)/ini_rho[1]) ; // 

  
  ini_rho[2] =  0.125 ;// 1.225
  Pp0[2] = 0.1 ;
  gamma[2] = 1.6;
  mu[2]      = 0  ; // 1.8e-5 
  c[2]       =  sqrt(gamma[2]* (Pp0[2] + 0)/ini_rho[2]) ; // 

  max_v          = 1.e-7;
  glbl_timestep  = 0.0;
  glbl_max_v     = 0.0;
  timestep_shift = 0.0;
  run_time       = 0.0;
  iterate_num    = 0;

  Lmax            = num_level-1 ;//mpm2025.9.11
  Lmin            = 0;

  ini_num_cell.i = DIM_X==1 ? int(NUM.i): 1;
  ini_num_cell.j = DIM_Y==1 ? int(NUM.j): 1;
  ini_num_cell.k = DIM_Z==1 ? int(NUM.k): 1;
  ini_scale = ini_h* CUT_OFF;//2026.7.8
//==============================================================
//==============================================================

////////////////////////////////////////////////////////////////


  static affinity_partitioner ap;

  parallel_for( blocked_range<int>(0, int(particle_fluid.size())),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i)
      {
        particle_fluid[i]->Clearup();
        particlepool.free(particle_fluid[i]);//2025.12.26
      }
  }, ap);

  parallel_for( blocked_range<int>(0, int(particle_bound.size())),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i)
      {
        particle_bound[i]->Clearup();
        particlepool.free(particle_bound[i]);//2026.1.1
      }
  }, ap);


    

//===================fluid===================================
  int nparticle= round(domain.i /p_sz);
  glbl_num_particle.i = nparticle;
  glbl_num_particle.j = round ( domain.j /p_sz );
  glbl_num_particle.k = DIM_Z==1 ? round ( domain.k /p_sz ): 1;

  int p_local_id = 0;
  my_int p_id;
  for( p_id.i= 0; p_id.i!= glbl_num_particle.i ; ++p_id.i)
    for( p_id.j= 0; p_id.j!= glbl_num_particle.j ; ++p_id.j)
      for( p_id.k= 0; p_id.k!= glbl_num_particle.k ; ++p_id.k){ 

      Vector<Real, 3> p_coord;
      p_coord(0) = (p_id.i+0.5)* p_sz;
      p_coord(1)= (p_id.j+0.5)* p_sz;
      p_coord(2) = DIM_Z==1 ? (p_id.k+0.5)* p_sz : 1;//2024.5.25


      bool f1 = (  p_coord(0) > 1 &&  p_coord(1) < 1.5 );


      p_Particle flu_particle      = particlepool.malloc();

      p_local_id ++;
      flu_particle->id          = p_local_id;

      flu_particle->color       = 0;
      flu_particle->coord = p_coord;
      flu_particle->coord_X = p_coord;
      flu_particle->coord_0 = p_coord;
      flu_particle->level      = 0; 
      flu_particle->h           = ini_h;

      flu_particle->P   = 0;//ini_rho *abs(GRAVITY) * (water_domain.j - p_coord(1)) 
      
      if ( f1 )
      {
        flu_particle->rho_0        = 1;
        flu_particle->mass        = powern(p_sz, DIM)* flu_particle->rho_0;
        flu_particle->vol         = powern(p_sz, DIM);//
        flu_particle->phase       = 1; 
        flu_particle->P       = 0.1; 

        flu_particle->rho[0]       = flu_particle->rho_0;
        flu_particle->rho0[0]      = flu_particle->rho_0;

        SetZero(flu_particle->rhov[0]);//
        SetZero(flu_particle->coef0[0]);//
        SetZero(flu_particle->coef[0]);//
        SetZero(flu_particle->dcoef[0]);//
        flu_particle->Lim_u[0].fill(1.);//
        flu_particle->mu[0]        = mu[0];
        flu_particle->E[0]       = (flu_particle->P + gamma[0]* Pp0[0])/(gamma[0]-1)  ; 
        flu_particle->c[0]         = sqrt(gamma[0]*(flu_particle->P + Pp0[0])/flu_particle->rho[0]);

      }
      else
      {

        bool f2 = (  p_coord(0) > 1);
        if(f2)
        {
          flu_particle->rho_0       = 0.125;
          flu_particle->mass        = powern(p_sz, DIM)* flu_particle->rho_0;
          flu_particle->vol         = powern(p_sz, DIM);//
          flu_particle->phase       = 2; 
          flu_particle->P       = 0.1;  //1e9
  
          flu_particle->rho[0]       = flu_particle->rho_0; //1250
          flu_particle->rho0[0]      = flu_particle->rho_0;
  
          SetZero(flu_particle->rhov[0]);//
          SetZero(flu_particle->coef0[0]);//
          SetZero(flu_particle->coef[0]);//
          SetZero(flu_particle->dcoef[0]);//
          flu_particle->Lim_u[0].fill(1.);//
          flu_particle->mu[0]        = mu[0];
          flu_particle->c[0]         = sqrt(gamma[1]*(flu_particle->P + Pp0[1])/flu_particle->rho[0]);
          flu_particle->E[0]       = (flu_particle->P + gamma[1]* Pp0[1])/(gamma[1]-1)  ; 

        }
        else{

          flu_particle->rho_0       = 1;
          flu_particle->mass        = powern(p_sz, DIM)* flu_particle->rho_0;
          flu_particle->vol         = powern(p_sz, DIM);//
          flu_particle->phase       = 3; 
          flu_particle->P       = 1; 

          flu_particle->rho[0]       = flu_particle->rho_0;
          flu_particle->rho0[0]      = flu_particle->rho_0;

          SetZero(flu_particle->rhov[0]);//
          SetZero(flu_particle->coef0[0]);//
          SetZero(flu_particle->coef[0]);//
          SetZero(flu_particle->dcoef[0]);//
          flu_particle->Lim_u[0].fill(1.);//
          flu_particle->mu[0]        = mu[2];
          flu_particle->E[0]       = (flu_particle->P + gamma[2]* Pp0[2])/(gamma[2]-1)  ; 
          flu_particle->c[0]         = sqrt(gamma[2]*(flu_particle->P + Pp0[2])/flu_particle->rho[0]);
        }

      }

      SetZero (flu_particle->a);
      SetZero (flu_particle->v);
      SetZero(flu_particle->rhov_0);//
      SetZero(flu_particle->coord_0);//


      flu_particle->neighbor.clear();
      particle_fluid.push_back(flu_particle);

    }

  total_num_particle =  p_local_id; //


////////////////////////////////////////////mpm2025.9.11


  Ele_cen = new p_Particle *[Nx ];
  for(int i=0; i< Nx ; i++)
    Ele_cen[i] = new p_Particle [Ny];

  for (int i = 0; i < Nx ; ++i) 
  for (int j = 0; j < Ny ; ++j)
  {
    Vector<Real, 3> nd_cen; SetZero(nd_cen);
    // Vector<Real,3> nd[ND_N]; 
    // for(int s=0; s < ND_N; s++)
    // {
    //   my_int tmp1 = Ele_list[i][j][s];
    //   nd[s] = node_fem[tmp1.i][tmp1.j];
    //   nd_cen += (nd[s])/ND_N;
    // }
    nd_cen(0) = px* (i +0.5);
    nd_cen(1) = py* (j +0.5);
    p_Particle ele_particle     = particlepool.malloc();

    ele_particle->coord = (nd_cen);
    ele_particle->coord_X = ele_particle->coord;
    ele_particle->id          = Nx* i + j;
    ele_particle->idI.i = i; ele_particle->idI.j = j; //2026.1.12
    ele_particle->h           = ini_h;
    ele_particle->level      = 0; 

    for (int si = 0; si < FldNum; si ++)
    {
      ele_particle->rho[si]       = ini_rho[si];
      ele_particle->rho0[si]      = ini_rho[si];

      SetZero(ele_particle->rhov[si]);//
      SetZero(ele_particle->coef0[si]);//
      SetZero(ele_particle->coef[si]);//
      SetZero(ele_particle->dcoef[si]);//
      ele_particle->Lim_u[si].fill(1.);//
      ele_particle->mu[si]        = mu[si];
      ele_particle->c[si]         = c[si];
      ele_particle->ph_state[si] = 1;
    }

    SetZero(ele_particle->rhov_0);//2025.12.25

    Ele_cen[i][j] = ele_particle ;
	}


  cell_dg = new Cell_dg *[Nx ];
  for(int i=0; i< Nx ; i++)
    cell_dg[i] = new Cell_dg [Ny];


/////////////////////////////////////
  for (int i = 0; i < Nx ; ++i) 
  for (int j = 0; j < Ny ; ++j)
  {
    cell_dg[i][j].face_center.clear();
    cell_dg[i][j].face_normal.clear();
    cell_dg[i][j].face_tan.clear();
    cell_dg[i][j].cell_id.clear();
    cell_dg[i][j].bound.clear();


    Vector<Real, 3> cen; SetZero(cen);
    Vector<Real, 3> norm; SetZero(norm);
    Vector<Real, 3> tang; SetZero(tang);
    my_int cell_j; 

    Vector<Real, 3> vol_r; SetZero(vol_r);
    vol_r = Ele_cen[i][j]->coord;
    // vol_r(0) = px* (i +0.5);
    // vol_r(1) = py* (j +0.5);

    norm(0)= 0; norm(1)= -1; 
    tang(0)= 0.5*px; tang(1)= 0; 
    cen(0) = vol_r(0); cen(1) = vol_r(1) - 0.5*py ; 
    cell_dg[i][j].face_center.push_back(cen);
    cell_dg[i][j].face_normal.push_back(norm);
    cell_dg[i][j].face_tan.push_back(tang);
    cell_dg[i][j].face_S.push_back(px);
    cell_j.i = i; cell_j.j = j-1; 
    cell_dg[i][j].cell_id.push_back(cell_j);
    if(cell_j.j >=0 ) cell_dg[i][j].bound.push_back(0);
    else cell_dg[i][j].bound.push_back(1);

    norm(0)= 1; norm(1)= 0; 
    tang(0)= 0; tang(1)= 0.5*py ; 
    cen(0) = vol_r(0) + 0.5* px; cen(1) = vol_r(1)  ; 
    cell_dg[i][j].face_center.push_back(cen);
    cell_dg[i][j].face_normal.push_back(norm);
    cell_dg[i][j].face_tan.push_back(tang);
    cell_dg[i][j].face_S.push_back(py);
    cell_j.i = i+1; cell_j.j = j; 
    cell_dg[i][j].cell_id.push_back(cell_j);
    if(cell_j.i < Nx) cell_dg[i][j].bound.push_back(0);
    else cell_dg[i][j].bound.push_back(1);

    norm(0)= 0; norm(1)= 1; 
    tang(0)= -0.5*px; tang(1)= 0; 
    cen(0) = vol_r(0) ; cen(1) = vol_r(1) + 0.5*py  ; 
    cell_dg[i][j].face_center.push_back(cen);
    cell_dg[i][j].face_normal.push_back(norm);
    cell_dg[i][j].face_tan.push_back(tang);
    cell_dg[i][j].face_S.push_back(px);
    cell_j.i = i; cell_j.j = j+1; 
    cell_dg[i][j].cell_id.push_back(cell_j);
    if(cell_j.j < Ny) cell_dg[i][j].bound.push_back(0);
    else cell_dg[i][j].bound.push_back(1);

    norm(0)= -1; norm(1)= 0; 
    tang(0)= 0; tang(1)= -0.5*py; 
    cen(0) = vol_r(0) - 0.5* px ; cen(1) = vol_r(1)  ; 
    cell_dg[i][j].face_center.push_back(cen);
    cell_dg[i][j].face_normal.push_back(norm);
    cell_dg[i][j].face_tan.push_back(tang);
    cell_dg[i][j].face_S.push_back(py);
    cell_j.i = i-1; cell_j.j = j; 
    cell_dg[i][j].cell_id.push_back(cell_j);
    if(cell_j.i >=0) cell_dg[i][j].bound.push_back(0);
    else cell_dg[i][j].bound.push_back(1);

    cell_dg[i][j].cell_center = vol_r;
    cell_dg[i][j].Vol = px* py;
  }



//////////////////////////////////////////////

  my_int pg; my_real b_p;
  pg.i = round(domain.i/p_sz); b_p.i = domain.i/pg.i;
  pg.j = round(domain.j/p_sz); b_p.j = domain.j/pg.j;


  for(int sx =0; sx <= pg.i; ++sx) 
  for(int sy=0; sy<= pg.j ; ++sy) 
  {

    Vector<Real, 3> dri; SetZero(dri);
    dri(0) = b_p.i *(sx + 0.);
    dri(1) = b_p.j *(sy + 0.);

    if(sx == 0 || sx == pg.i || sy == 0 || sy == pg.j )
    {

      p_Particle b_particle     = particlepool.malloc();

      b_particle->coord = (dri);
      b_particle->coord_X = b_particle->coord;
      b_particle->id          = pg.i* sx + sy;
      b_particle->h           = ini_h;
      b_particle->level       = 0; 
      b_particle->phase       = 10; 
      b_particle->vol         = p_sz; 

      Vector<Real, 3> normal; SetZero(normal);
      if(sx == 0) {normal(0) = -1; normal(1) = 0; }
      else if(sx == pg.i) {normal(0) = 1; normal(1) = 0; }
      else if(sy == 0) {normal(0) = 0; normal(1) = -1; }
      else if(sy == pg.j) {normal(0) = 0; normal(1) = 1; }
      b_particle->norm = normal;



      b_particle->rho_0        = ini_rho[0];
      SetZero(b_particle->rhov_0);//2025.12.25
      SetZero(b_particle->tv);//

      for (int si = 0; si < FldNum; si ++)
      {
        b_particle->mu[si]        = mu[si];
        b_particle->c[si]         = c[si];
        b_particle->rho[si]       = ini_rho[si];
        b_particle->rho0[si]      = ini_rho[si];

        SetZero(b_particle->rhov[si]);//
        SetZero(b_particle->coef0[si]);//
        SetZero(b_particle->coef[si]);//
        SetZero(b_particle->dcoef[si]);//
        b_particle->Lim_u[si].fill(1.);//
      }

      particle_bound.push_back(b_particle);
    
    }
  }


  // for(int sx =0; sx <= pg.i; ++sx) 
  // for(int sy=0; sy<= 1 ; ++sy) 
  // {
  //   Vector<Real, 3> dri; SetZero(dri);
  //   dri(0) = b_p.i *(sx + 0.5);
  //   dri(1) = domain.j *(sy + 0.);
  //   {
  //     p_Particle b_particle     = particlepool.malloc();
  //     b_particle->coord = (dri);
  //     b_particle->coord_X = b_particle->coord;
  //     b_particle->id          = pg.i* sx + sy;
  //     b_particle->h           = ini_h;
  //     b_particle->level       = 0; 
  //     b_particle->phase       = 10; 
  //     b_particle->vol         = p_sz; 
  //     Vector<Real, 3> normal; SetZero(normal);
  //     if(sy == 0) {normal(0) = 0; normal(1) = -1; }
  //     else if(sy == 1) {normal(0) = 0; normal(1) = 1; }
  //     b_particle->norm = normal;
  //     b_particle->nu          = nu;
  //     b_particle->c           = c;
  //     b_particle->rho         = ini_rho;
  //     b_particle->rho0        = ini_rho;
  //     b_particle->rho_0        = ini_rho;
  //     SetZero(b_particle->rhov);//
  //     SetZero(b_particle->coef0);//
  //     SetZero(b_particle->coef);//
  //     SetZero(b_particle->dcoef);//
  //     b_particle->Lim_u.fill(1.);//
  //     SetZero(b_particle->rhov_0);//2025.12.25
  //     SetZero(b_particle->tv);//
  //     particle_bound.push_back(b_particle);
  //   }
  // }
  // for(int sy=0; sy<= pg.j ; ++sy) 
  // for(int sx =0; sx <= 1; ++sx) 
  // {
  //   Vector<Real, 3> dri; SetZero(dri);
  //   dri(0) = domain.i *(sx + 0.);
  //   dri(1) = b_p.j *(sy + 0.5);
  //   {
  //     p_Particle b_particle     = particlepool.malloc();
  //     b_particle->coord = (dri);
  //     b_particle->coord_X = b_particle->coord;
  //     b_particle->id          = pg.i* sx + sy;
  //     b_particle->h           = ini_h;
  //     b_particle->level       = 0; 
  //     b_particle->phase       = 10; 
  //     b_particle->vol         = p_sz; 
  //     Vector<Real, 3> normal; SetZero(normal);
  //     if(sx == 0) {normal(0) = -1; normal(1) = 0; }
  //     else if(sx == 1) {normal(0) = 1; normal(1) = 0; }
  //     b_particle->norm = normal;
  //     b_particle->nu          = nu;
  //     b_particle->c           = c;
  //     b_particle->rho         = ini_rho;
  //     b_particle->rho0        = ini_rho;
  //     b_particle->rho_0        = ini_rho;
  //     SetZero(b_particle->rhov);//
  //     SetZero(b_particle->coef0);//
  //     SetZero(b_particle->coef);//
  //     SetZero(b_particle->dcoef);//
  //     b_particle->Lim_u.fill(1.);//
  //     SetZero(b_particle->rhov_0);//2025.12.25
  //     SetZero(b_particle->tv);//
  //     particle_bound.push_back(b_particle);
  //   }
  // }

////////////////////////////////////////////////////////////

  // for(int i=0; i!= int(cell_dg[12][0].cell_id.size()); ++i)
  // {
  //   cout<< cell_dg[12][0].cell_id[i].i<<","<<cell_dg[12][0].cell_id[i].j<<endl;//
  //   cout<< cell_dg[12][0].bound[i]<<endl;//
  // }

  int ed_r = 1 ;
  parallel_for( blocked_range2d<int>(0, NUM.i, 0, NUM.j), [&](const blocked_range2d<int>& r){
    for(int i= r.rows().begin(); i< r.rows().end(); i++)
      for(int j = r.cols().begin(); j < r.cols().end(); j++){ //2026.1.12
        Ele_cen[i][j]->neigh_cell.clear(); 
        for(int t=AMAX1(i- ed_r, 0); t<=AMIN1(i+ ed_r, NUM.i-1); t++)
        for(int s=AMAX1(j- ed_r, 0); s<=AMIN1(j+ ed_r, NUM.j-1); s++)
        Ele_cen[i][j]->neigh_cell.push_back(Ele_cen[t][s]) ; 

    }
  }, ap);


  level_info = new p_Level_info[num_level];
  for(int i =0; i<num_level; i++) 
  {
    level_info[i] = new Level_info;
    level_info[i]->Initialize(i, this, world);
  }


  glbl_total_num_particle = total_num_particle;
  if (world.rank() == 0){
    cout<<"**********************************************************\n";
    my_cout (glbl_num_particle, "Number of Particles:  ");
    my_cout (domain,       "Domain size:          ");
    my_cout (box_l,        "Left corner:          ");
    my_cout (box_r,        "Right corner:         ");
    cout<<"   Total number particle:  "<<glbl_total_num_particle<<"\n";
    cout<<"   Initial smooth length:  "<<ini_h<<"\n";
    cout<<"**********************************************************\n";
  }

}





//--------------------------------------------------
// mpm incompressible solver
//--------------------------------------------------
void MPM::MPM_run_simulation(int n, communicator &world, Real Time, int flag1)//2022.4.6
{
/**********************************************************/

 time_simulation.restart();




/**********************************************************/


  if(n==1  ) // || n% 20 ==0
  Used_state_update(world);//2026.1.12
  // cout<<"WRONG111"<<endl;

  Reset_cell_list_info(world);
  Map_the_particle_to_tree(world); 

  // // /**********************************************************/
  // Split_particle(world);
  // Reset_cell_list_info(world);
  // Map_the_particle_to_tree(world); 
  // Reset_refresh_merge_neighbor(world);


  Split_particle_fast(world);//2026.9.21

  for (int iteration = 0; iteration < 4; iteration++)
  {
    // Merge_particle(world);
    // Reset_cell_list_info(world);
    // Map_the_particle_to_tree(world); 
    // Reset_refresh_merge_neighbor(world);

    Merge_particle_fast(world);
  }
  /**********************************************************/

  int md =0;
  time_reset_neighbor.restart();
  Reset_particle_neighbor_info(world, 4);//
  time_for_reset_neighbor += time_reset_neighbor.elapsed();

  time_refresh_neighbor.restart();
  Refresh_neighbor_info(4, world);
  time_for_refresh_neighbor  += time_refresh_neighbor.elapsed(); //2026.7.2


  while(md< 2) 
  {
    time_reset_neighbor.restart();
    Reset_particle_neighbor_info(world, 2);//
    time_for_reset_neighbor += time_reset_neighbor.elapsed();

    time_refresh_neighbor.restart();
    Refresh_neighbor_info(2, world);
    time_for_refresh_neighbor  += time_refresh_neighbor.elapsed(); //2026.7.8
  
    ///////////////////////////////
    // if(md==0)
    // Identify_phase(world);//

    // if(md==1)
    // Free_surf_voro_particle(world);


    if(md==0)
    MLS_P2G(world);
    // Output_plt_file(100, 1, world);
    // cout<<"WRONG222"<<endl;




    Calculate_limiter(world);//
    // Output_plt_file(100, 1, world);
    // cout<<"WRONG444"<<endl;

    Grid_acceleration( world); 
    // Output_plt_file(200, 1, world);
    // cout<<"WRONG555"<<endl;

    Update_grid(world, md); 
    // Output_plt_file(250, 1, world);
    // cout<<"WRONG666"<<endl;
    // exit(0);


    Update_particle(world, md);
    // Output_plt_file(300, 1, world);
    // cout<<"WRONG777"<<endl;
    // exit(0);

    md++;
  }


  if(flag1 ==1)
  Output_point(world, Time);

    time_for_simulation += time_simulation.elapsed();
/**********************************************************/
}




//--------------------------------------------------
//  Output_point //2023.9.6
//--------------------------------------------------
void MPM::Output_point(communicator &world, Real time1)
{
	Iden += 1;

  p_Particle cp1 =  Ele_cen[NUM.i-1][0];

  Real P1 = cp1->P;


if(Iden==1){
 ofstream out("./pointa.txt",ios::trunc);  
  out<<time1<<" "<<setprecision(6)<<P1<<"\n";
  out.close();

}else{
 ofstream out("./pointa.txt",ios::out|ios::ate|ios::app);  
  out<<time1<<" "<<setprecision(6)<<P1<<"\n";
  out.close();

}

  
}




//--------------------------------------------------
// Output particle infomation
//--------------------------------------------------
void MPM::Output_plt_file(int n,  int flag, communicator &world)
{

/////////////////////////////////////////////////////////

  char    filename[256];
  sprintf(filename,"%s%d%s%d%s","./outdata/simulation_sp.",n,".",world.rank(),".plt");

  // ofstream out(filename, ios::trunc);
  // out<<"VARIABLES = \"x\",\"y\",\"phase\", \"vx\",\"vy\",\"vz\",\"rho\",\"P\" \
	//    \"surf\", \"v\",  \"n_neighbor\"\n";
  // out<<"ZONE T="<<"\"GRID 1\","<<"SOLUTIONTIME="<<run_time<<"\n";


  // for (int i = 0; i < int(particle_fluid.size()); i++)
  // {
  //   Particle *current_particle = particle_fluid[i]; //

  //   out<<current_particle->coord(0)<<' '<<current_particle->coord(1)<<' '<<current_particle->phase<<' '
  //      <<current_particle->v(0)<<' '<<current_particle->v(1)<<' '<<current_particle->vol<<' '
  //      <<current_particle->rho[0]<<' '<<current_particle->P<<' '

  //      <<(current_particle->tv).norm()<<' ' 
  //      <<(current_particle->v).norm()<<' '
  //     <<current_particle->E[0]<<endl;//
  // }
  // for (int i = 0; i < int(particle_bound.size()); i++)
  // {
  //   Particle *current_particle = particle_bound[i]; //

  //   out<<current_particle->coord(0)<<' '<<current_particle->coord(1)<<' '<<current_particle->phase<<' '
  //      <<current_particle->v(0)<<' '<<current_particle->v(1)<<' '<<current_particle->v(2)<<' '
  //      <<current_particle->rho[0]<<' '<<current_particle->P<<' '

  //      <<current_particle->interior_mark<<' ' 
  //      <<(current_particle->v).norm()<<' '
  //     <<current_particle->neighbor.size()<<endl;//
  // }
  // out.close();


  ofstream out(filename, ios::binary);
  char  versn [] =  "#!TDV112";
  // Write_file_c(out, versn, sizeof(versn));
  out.write((char*)&(versn), sizeof(versn)-1);
  Write_file(out, 1);
  Write_file(out, 0);
  Write_file(out, 0);

  vector<string>  VarName = {"x","y", "phase",  "vx", "vy", "vol", "rho", "P","tv","v","n_neighbor"}; 
  vector<int>  VarTp = {2,2, 3,  2,2,2,  2,2, 2,2,3 }; 

  Write_file(out, int(VarName.size()));

  for(string& varN: VarName)
  {
    for(char varC : varN ) 
    {
      out.write((char*)& (varC), sizeof(varC)); 
      int t=0;
      out.write((char*)& (t), sizeof(int)-sizeof(varC)); 
    }
      
    Write_file(out, 0);
  }

  Write_file(out, float(299.0));
  VarName = {"grid"}; 
  for(string& varN: VarName) 
  for(char varC : varN ) 
  {
    out.write((char*)& (varC), sizeof(varC)); 
    int t=0;
    out.write((char*)& (t), sizeof(int)-sizeof(varC)); 
  }
  Write_file(out, 0);
  
  Write_file(out, -1);
  Write_file(out, 0);//strand
  Write_file(out, double(run_time));//time
  Write_file(out, -1);//zone color
  Write_file(out, 0);//zonetype
  Write_file(out, 0);//specify var
  Write_file(out, 0);//no face
  Write_file(out, 0);//no miscell

  int ss=0;
	for (int i = 0; i < int(particle_fluid.size()); i++)
  ss ++;

  Write_file(out, int(ss) );//IMAX
  Write_file(out, 1);//jMAX
  Write_file(out, 1);//KMAX
  Write_file(out, 0);//
  Write_file(out, float(357.0));//header end
  Write_file(out, float(299.0));//data start
  for(int j=0;j< int(VarTp.size());j++) Write_file(out, VarTp[j]);//data type
  Write_file(out, 0);//no passive
  Write_file(out, 0);//no sharing
  Write_file(out, -1);// sharing info -1

  for(int j=0;j< int(VarTp.size());j++)
  {
    Write_file(out, double(-1.));// xmin
    Write_file(out, double(1.));// xmax
  }


  // int s=num_level-1; //2024.4.29
  for(int j=0;j< int(VarTp.size());j++)
	for (int i = 0; i < int(particle_fluid.size()); i++)
  {
    Particle *cp = particle_fluid[i]; //2022.2.11

    if(j==0)
      Write_file(out,  Convert_T(cp->coord(0)));
    else if(j==1)
      Write_file(out,  Convert_T(cp->coord(1) ));
    else if(j==2)
      Write_file(out,  Convert_T(cp->phase ));
    else if(j==3)
      Write_file(out,  Convert_T(cp->v(0) ));
    else if(j==4)
      Write_file(out,  Convert_T(cp->v(1) ));
    else if(j==5)
      Write_file(out,  Convert_T(cp->vol ));
    else if(j==6)
      Write_file(out,  Convert_T(cp->rho[0] ));
    else if(j==7)
      Write_file(out,  Convert_T(cp->P ));
    else if(j==8)
      Write_file(out,  Convert_T(cp->tv.head(DIM).norm() ) );
    else if(j==9)
      Write_file(out,  Convert_T(cp->v.head(DIM).norm() ));
    else if(j==10)
      Write_file(out,  Convert_T(int(cp->neighbor.size()) ));

  }
  out.close();





  char    filename1[256];
  sprintf(filename1,"%s%d%s%d%s","./outdata/ele_file.",n,".",world.rank(),".plt");

  ofstream out1(filename1, ios::trunc);
  out1<<"VARIABLES = \"x\",\"y\",\"z\", \"v1\",\"v2\", \"v\",\"rho0\",\"rho1\", \
       \"E0\",\"E1\",\"Lim0\",\"ph0\",\"ph1\", \"ph2\", \"P\"\n";
  out1<<"ZONE T="<<"\"GRID 1\","<<"SOLUTIONTIME="<<run_time<<"\n";

  for (int i = 0; i < int(NUM.i); i++)
  for (int j = 0; j < int(NUM.j); j++)
  {
    Particle *cp = Ele_cen[i][j]; //

    out1<<cp->coord(0)<<' '<<cp->coord(1)<<' '<<cp->coord(2)<<' '
    <<cp->v(0)<<' '<<cp->v(1)<<' '<<cp->v.head(DIM).norm()<<' '
    <<cp->rho[0]<<' '<<cp->rho[1]<<' '
    // <<cp->E[0]<<' '<<cp->E[1]<<' '
      <<cp->coef[1](2,1)<<' '<<cp->coef[1](2,2)<<' '
      // <<(cp->rhov[0]/cp->rho[0]).head(DIM).norm()<<' '<<(cp->rhov[1]/cp->rho[1]).head(DIM).norm()<<' '
      // <<cp->Lim_u[0](3)<<' '<<cp->Lim_u[1](3)<<' '
      <<cp->Lim_u[0](0)<<' '
      <<cp->alph[0]<<' '
      <<cp->alph[1]<<' '
      <<cp->alph[2]<<' '
      <<cp->P<<' '
      <<endl;
      
      
  }
  out1.close();



}

