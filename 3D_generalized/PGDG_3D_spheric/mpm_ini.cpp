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


  Gp = {-1./sqrt(3.), 1./sqrt(3.) }; 
  Gw = {0.5, 0.5}; 

  Gpf = {-1./sqrt(3.), 1./sqrt(3.) }; 
  Gwf = {0.5, 0.5 }; 


  Real Rt = 1 ; //
  my_real water_domain;
  water_domain.i = 4 * Rt;
  water_domain.j = 1.5* Rt;
  
  domain.i = particle_box.i = 0.350 *Rt ;  //5.366
  domain.j = particle_box.j = 0.0445 *Rt;
  domain.k = particle_box.k = 0.0445 *Rt;

  num_level =  1 ; 

  p_sz = domain.j/30 * Rt ;//

  p_g = domain.j/30  * Rt ;//2026.1.12


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

  #if DIM==3
  int Nz = domain.k/p_g; Real pz = domain.k/Real(Nz);
  NUM.k = Nz; 
  p_size.k = pz; //2026.2.18
  #endif

  #if DIM==2
  p_sz = pow(px* py, 1./DIM)/round(p_g/p_sz) ;//
  #elif DIM==3
  p_sz = pow(px* py*pz, 1./DIM)/round(p_g/p_sz) ;//
  #endif

  ini_h          = CELL_RATIO/CUT_OFF* p_sz;

  cout<<px<<", "<<py<<endl;

  
  ini_rho[0] = 1.29 ;//
  Pp0[0] = 0;
  gamma[0] = 1.4; 
  mu[0]      = 0 ; // 1e-3 
  c[0]      =  sqrt(gamma[0]* (Pp0[0] + 0)/ini_rho[0]) ;// 

  ini_rho[1] =  0.2347 ;// 1.225
  Pp0[1] = 0;
  gamma[1] = 1.648 ;
  mu[1]      = 0  ; // 1.8e-5 
  c[1]       =  sqrt(gamma[1]* (Pp0[1] + 0)/ini_rho[1]) ; // 

  

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

      
/////////////////////////////////////

  #if DIM ==2
  Ele_cen = new p_Particle *[Nx ];
  cell_dg = new Cell_dg *[Nx ];
  #elif DIM ==3
  Ele_cen = new p_Particle **[Nx ];
  cell_dg = new Cell_dg **[Nx ];
  #endif //2026.2.18
  for(int i=0;i< Nx ; i++)
  {
    #if DIM ==2
    Ele_cen[i] = new p_Particle [Ny ];
    cell_dg[i] = new Cell_dg [Ny];
    #elif DIM ==3
    Ele_cen[i] = new p_Particle *[Ny ];
    cell_dg[i] = new Cell_dg *[Ny];
    for(int j=0;j< Ny ; j++)
    {
      Ele_cen[i][j] = new p_Particle [Nz];
      cell_dg[i][j] = new Cell_dg [Nz];
    }  
    #endif //2026.2.18
  }



///////////////////////////////////////////////////////////


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


      bool f1 = (  p_coord(0) < 0.1 ); //2026.7.24


      p_Particle flu_particle      = particlepool.malloc();

      p_local_id ++;
      flu_particle->id          = p_local_id;

      flu_particle->color       = 0;
      flu_particle->coord = p_coord;
      flu_particle->coord_0 = p_coord;
      flu_particle->level      = 0; 
      flu_particle->h           = ini_h;

      flu_particle->P   = 0;//
      
      if ( f1 )
      {
        flu_particle->rho_0        =  1.697 ;
        flu_particle->mass        = powern(p_sz, DIM)* flu_particle->rho_0;
        flu_particle->vol         = powern(p_sz, DIM);//
        flu_particle->phase       = 1; 
        flu_particle->P       = 167800 ; 

        flu_particle->rho[0]       = flu_particle->rho_0;
        flu_particle->rho0[0]      = flu_particle->rho_0;


        SetZero (flu_particle->v);
        flu_particle->v(0) = 131.6 ;
        flu_particle->rhov[0] = flu_particle->v* flu_particle->rho[0]; 
        SetZero(flu_particle->coef0[0]);//
        SetZero(flu_particle->coef[0]);//
        SetZero(flu_particle->dcoef[0]);//

        flu_particle->Lim_u[0].fill(1.);//
        flu_particle->mu[0]        = mu[0];
        flu_particle->E[0]       = (flu_particle->P + gamma[0]* Pp0[0])/(gamma[0]-1)  ; 
        flu_particle->E[0]   += 0.5* flu_particle->rho[0]* flu_particle->v.head(DIM).squaredNorm() ; //202606619
        flu_particle->c[0]         = sqrt(gamma[0]*(flu_particle->P + Pp0[0])/flu_particle->rho[0]);

      }
      else
      {
#if DIM==2
        bool f2 = (  pow(p_coord(0) - 0.15, 2) + pow(p_coord(1) , 2)  > pow(0.0225, 2)); //2026.7.24
#elif DIM==3
        bool f2 = (  pow(p_coord(0) - 0.15, 2) + pow(p_coord(1) , 2) + pow(p_coord(2) , 2)  > pow(0.0225, 2)); //2026.7.24
#endif
        if(f2)
        {
          flu_particle->rho_0       = 1.179 ;
          flu_particle->mass        = powern(p_sz, DIM)* flu_particle->rho_0;//2026.7.16
          flu_particle->vol         = powern(p_sz, DIM);//
          flu_particle->phase       = 1; 
          flu_particle->P       = 101325 ;  //1e9
  
          flu_particle->rho[0]       = flu_particle->rho_0; //1250
          flu_particle->rho0[0]      = flu_particle->rho_0;
  
          SetZero (flu_particle->v);
          SetZero(flu_particle->rhov[0]);//
          SetZero(flu_particle->coef0[0]);//
          SetZero(flu_particle->coef[0]);//
          SetZero(flu_particle->dcoef[0]);//
          flu_particle->Lim_u[0].fill(1.);//
          flu_particle->mu[0]        = mu[0];
          flu_particle->c[0]         = sqrt(gamma[0]*(flu_particle->P + Pp0[0])/flu_particle->rho[0]);
          flu_particle->E[0]       = (flu_particle->P + gamma[0]* Pp0[0])/(gamma[0]-1)  ; 
          flu_particle->E[0]   += 0.5* flu_particle->rho[0]* flu_particle->v.head(DIM).squaredNorm() ; //202606619

        }
        else{

          flu_particle->rho_0       = 0.2347 ;
          flu_particle->mass        = powern(p_sz, DIM)* flu_particle->rho_0;
          flu_particle->vol         = powern(p_sz, DIM);//
          flu_particle->phase       = 2; 
          flu_particle->P       = 101325; 

          flu_particle->rho[0]       = flu_particle->rho_0;
          flu_particle->rho0[0]      = flu_particle->rho_0;

          SetZero (flu_particle->v);
          SetZero(flu_particle->rhov[0]);//
          SetZero(flu_particle->coef0[0]);//
          SetZero(flu_particle->coef[0]);//
          SetZero(flu_particle->dcoef[0]);//
          flu_particle->Lim_u[0].fill(1.);//
          flu_particle->mu[0]        = mu[1];
          flu_particle->E[0]       = (flu_particle->P + gamma[1]* Pp0[1])/(gamma[1]-1)  ; 
          flu_particle->E[0]   += 0.5* flu_particle->rho[0]* flu_particle->v.head(DIM).squaredNorm() ; //202606619
          flu_particle->c[0]         = sqrt(gamma[1]*(flu_particle->P + Pp0[1])/flu_particle->rho[0]);

        }
      }

      SetZero (flu_particle->a);
      SetZero(flu_particle->rhov_0);//
      SetZero(flu_particle->coord_0);//


      flu_particle->neighbor.clear();
      particle_fluid.push_back(flu_particle);

    }

  total_num_particle =  p_local_id; //


////////////////////////////////////////////mpm2025.9.11

  for (int i = 0; i < Nx ; ++i) 
  for (int j = 0; j < Ny ; ++j)
  #if DIM ==3
  for (int k = 0; k < Nz ; ++k)
  #endif //2026.2.18
  {
    Vector<Real, 3> nd_cen; SetZero(nd_cen);

    nd_cen(0) = px* (i +0.5);
    nd_cen(1) = py* (j +0.5);
    #if DIM ==3
    nd_cen(2) = pz* (k +0.5);
    #endif //2026.2.18
    p_Particle ele_particle     = particlepool.malloc();

    ele_particle->coord = (nd_cen);
    ele_particle->id          = Ny* i + j;
    ele_particle->idI.i = i; ele_particle->idI.j = j; 
    #if DIM ==3
    ele_particle->idI.k = k; //2026.2.18
    ele_particle->id  = Nz*Ny* i + Nz*j + k;
    #endif //2026.2.18

    ele_particle->h           = ini_h;
    ele_particle->level      = 0; 

    for (int si = 0; si < 2; si ++)
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

    #if DIM ==2
    Ele_cen[i][j] = ele_particle ;
    #elif DIM ==3
    Ele_cen[i][j][k] = ele_particle ;
    #endif //2026.2.18
	} //2026.7.24




/////////////////////////////////////
  for (int i = 0; i < Nx ; ++i) 
  for (int j = 0; j < Ny ; ++j)
  #if DIM ==3
  for (int k = 0; k < Nz ; ++k)
  #endif //2026.2.18
  {
    Cell_dg *Celli;
    #if DIM ==2
    Celli = &cell_dg[i][j];
    #elif DIM ==3
    Celli = &cell_dg[i][j][k];
    #endif //2026.2.18


    Celli->face_center.clear();
    Celli->face_normal.clear();
    Celli->face_tan.clear();
    Celli->cell_id.clear();
    Celli->bound.clear();


    Vector<Real, 3> cen; SetZero(cen);
    Vector<Real, 3> norm; SetZero(norm);
    Vector<Real, 3> tang; SetZero(tang);

    my_int cell_j; 
    int bmk=0;

    Vector<Real, 3> vol_r; SetZero(vol_r);
    // vol_r = Ele_cen[i][j]->coord;
    vol_r(0) = px* (i +0.5);
    vol_r(1) = py* (j +0.5);
    #if DIM ==3
    vol_r(2) = pz* (k +0.5);
    #endif //2026.2.18

  #if DIM ==2
    
    norm(0)= 0; norm(1)= -1; 
    tang(0)= 0.5*px; tang(1)= 0; 
    cen(0) = vol_r(0); cen(1) = vol_r(1) - 0.5*py ; 
    Celli->face_center.push_back(cen);
    Celli->face_normal.push_back(norm);
    Celli->face_tan.push_back(tang);
    Celli->face_S.push_back(px);
    cell_j.i = i; cell_j.j = j-1; 
    Celli->cell_id.push_back(cell_j);
    bmk=0;
    if(cell_j.j >=0) bmk = 0;
    else bmk = -3; //symmetric
    Celli->bound.push_back(bmk);

    norm(0)= 1; norm(1)= 0; 
    tang(0)= 0; tang(1)= 0.5*py ; 
    cen(0) = vol_r(0) + 0.5* px; cen(1) = vol_r(1)  ; 
    Celli->face_center.push_back(cen);
    Celli->face_normal.push_back(norm);
    Celli->face_tan.push_back(tang);
    Celli->face_S.push_back(py);
    cell_j.i = i+1; cell_j.j = j; 
    Celli->cell_id.push_back(cell_j);
    bmk=0;
    if(cell_j.i < Nx) bmk = 0;
    else bmk = -2; //outlet
    Celli->bound.push_back(bmk);

    norm(0)= 0; norm(1)= 1; 
    tang(0)= -0.5*px; tang(1)= 0; 
    cen(0) = vol_r(0) ; cen(1) = vol_r(1) + 0.5*py  ; 
    Celli->face_center.push_back(cen);
    Celli->face_normal.push_back(norm);
    Celli->face_tan.push_back(tang);
    Celli->face_S.push_back(px);
    cell_j.i = i; cell_j.j = j+1; 
    Celli->cell_id.push_back(cell_j);
    if(cell_j.j < Ny) Celli->bound.push_back(0);
    else Celli->bound.push_back(1);

    norm(0)= -1; norm(1)= 0; 
    tang(0)= 0; tang(1)= -0.5*py; 
    cen(0) = vol_r(0) - 0.5* px ; cen(1) = vol_r(1)  ; 
    Celli->face_center.push_back(cen);
    Celli->face_normal.push_back(norm);
    Celli->face_tan.push_back(tang);
    Celli->face_S.push_back(py);
    cell_j.i = i-1; cell_j.j = j; 
    Celli->cell_id.push_back(cell_j);
    bmk=0;
    if(cell_j.i >=0) bmk = 0;
    else bmk = -1; //inlet 20260618
    Celli->bound.push_back(bmk);

    Celli->cell_center = vol_r;
    Celli->Vol = px* py;

  #elif DIM ==3
    //yL
    norm(0)= 0; norm(1)= -1; norm(2)= 0; 
    tang(0)= 0.5*px; tang(1)= 0; tang(2)= 0; 
    cen(0) = vol_r(0); cen(1) = vol_r(1) - 0.5*py ; cen(2) = vol_r(2) ; 
    Celli->face_center.push_back(cen);
    Celli->face_normal.push_back(norm);
    Celli->face_tan.push_back(tang);
    tang(0)= 0; tang(1)= 0; tang(2)= 0.5*pz; 
    Celli->face_tan.push_back(tang);//
    Celli->face_S.push_back(px*pz);
    cell_j.i = i; cell_j.j = j-1; cell_j.k = k; 
    Celli->cell_id.push_back(cell_j);
    bmk=0;
    if(cell_j.j >=0) bmk = 0;
    else bmk = -3; //symmetric
    Celli->bound.push_back(bmk);

    //xR
    norm(0)= 1; norm(1)= 0;  norm(2)= 0; 
    tang(0)= 0; tang(1)= 0.5*py ; tang(2)= 0; 
    cen(0) = vol_r(0) + 0.5* px; cen(1) = vol_r(1) ; cen(2) = vol_r(2) ; 
    Celli->face_center.push_back(cen);
    Celli->face_normal.push_back(norm);
    Celli->face_tan.push_back(tang);
    tang(0)= 0; tang(1)= 0; tang(2)= 0.5*pz; 
    Celli->face_tan.push_back(tang);//
    Celli->face_S.push_back(py*pz);
    cell_j.i = i+1; cell_j.j = j; cell_j.k = k; 
    Celli->cell_id.push_back(cell_j);
    bmk=0;
    if(cell_j.i < Nx) bmk = 0;
    else bmk = -2; //outlet
    Celli->bound.push_back(bmk);

    //yR
    norm(0)= 0; norm(1)= 1; norm(2)= 0; 
    tang(0)= -0.5*px; tang(1)= 0; tang(2)= 0; 
    cen(0) = vol_r(0) ; cen(1) = vol_r(1) + 0.5*py  ; cen(2) = vol_r(2) ; 
    Celli->face_center.push_back(cen);
    Celli->face_normal.push_back(norm);
    Celli->face_tan.push_back(tang);
    tang(0)= 0; tang(1)= 0; tang(2)= 0.5*pz; 
    Celli->face_tan.push_back(tang);//
    Celli->face_S.push_back(px*pz);
    cell_j.i = i; cell_j.j = j+1; cell_j.k = k; 
    Celli->cell_id.push_back(cell_j);
    if(cell_j.j < Ny) Celli->bound.push_back(0);
    else Celli->bound.push_back(1);

    //xL
    norm(0)= -1; norm(1)= 0; norm(2)= 0; 
    tang(0)= 0; tang(1)= -0.5*py; tang(2)= 0; 
    cen(0) = vol_r(0) - 0.5* px ; cen(1) = vol_r(1) ; cen(2) = vol_r(2) ; 
    Celli->face_center.push_back(cen);
    Celli->face_normal.push_back(norm);
    Celli->face_tan.push_back(tang);
    tang(0)= 0; tang(1)= 0; tang(2)= 0.5*pz; 
    Celli->face_tan.push_back(tang);//
    Celli->face_S.push_back(py*pz);
    cell_j.i = i-1; cell_j.j = j; cell_j.k = k; 
    Celli->cell_id.push_back(cell_j);
    bmk=0;
    if(cell_j.i >=0) bmk = 0;
    else bmk = -1; //inlet 20260618
    Celli->bound.push_back(bmk);

    //zL
    norm(0)= 0; norm(1)= 0; norm(2)= -1; 
    tang(0)= 0.5*px; tang(1)= 0; tang(2)= 0; 
    cen(0) = vol_r(0) ; cen(1) = vol_r(1) ; cen(2) = vol_r(2) - 0.5*pz ; 
    Celli->face_center.push_back(cen);
    Celli->face_normal.push_back(norm);
    Celli->face_tan.push_back(tang);
    tang(0)= 0; tang(1)= 0.5*py; tang(2)= 0; 
    Celli->face_tan.push_back(tang);//
    Celli->face_S.push_back(py*px);
    cell_j.i = i; cell_j.j = j; cell_j.k = k-1; 
    Celli->cell_id.push_back(cell_j);
    bmk=0;
    if(cell_j.k >=0) bmk = 0;
    else bmk = -3;
    Celli->bound.push_back(bmk);

    //zR
    norm(0)= 0; norm(1)= 0; norm(2)= 1; 
    tang(0)= 0.5*px; tang(1)= 0; tang(2)= 0; 
    cen(0) = vol_r(0) ; cen(1) = vol_r(1) ; cen(2) = vol_r(2) + 0.5*pz ; 
    Celli->face_center.push_back(cen);
    Celli->face_normal.push_back(norm);
    Celli->face_tan.push_back(tang);
    tang(0)= 0; tang(1)= 0.5*py; tang(2)= 0; 
    Celli->face_tan.push_back(tang);//
    Celli->face_S.push_back(py*px);
    cell_j.i = i; cell_j.j = j; cell_j.k = k+1; 
    Celli->cell_id.push_back(cell_j);
    if(cell_j.k <Nz) Celli->bound.push_back(0);
    else Celli->bound.push_back(1);


    Celli->cell_center = vol_r;
    Celli->Vol = px* py* pz;
  #endif //2026.2.18
  }


//////////////////////////////////////////////

  for(int sx =0; sx <= NUM.i; ++sx) 
  for(int sy=0; sy<= NUM.j ; ++sy) 
  #if DIM ==3
  for(int sz=0; sz<= NUM.k ; ++sz) 
  #endif //2026.2.18
  {

    Vector<Real, 3> dri; SetZero(dri);
    dri(0) = p_size.i *(sx + 0.);
    dri(1) = p_size.j *(sy + 0.);
    #if DIM ==3
    dri(2) = p_size.k *(sz + 0.);
    #endif //2026.2.18

    bool fg = (sx == 0 || sx == NUM.i || sy == 0 || sy == NUM.j);

    #if DIM ==3
    fg = (fg || sz == 0 || sz == NUM.k);
    #endif //2026.2.18

    if( fg)
    {
      p_Particle b_particle     = particlepool.malloc();

      b_particle->coord = (dri);
      b_particle->id          = NUM.j* sx + sy;
      #if DIM ==3
      b_particle->id          = NUM.j* NUM.k* sx + sy*NUM.k + sz;
      #endif //2026.2.18
      b_particle->h           = ini_h;
      b_particle->level       = 0; 
      b_particle->phase       = 10; 
      b_particle->vol         = powern(p_sz, DIM-1);//2026.2.18

      Vector<Real, 3> normal; SetZero(normal);
      if(sx == 0) {normal(0) = -1; normal(1) = 0; normal(2) = 0;}
      else if(sx == NUM.i) {normal(0) = 1; normal(1) = 0; normal(2) = 0;}
      else if(sy == 0) {normal(0) = 0; normal(1) = -1;normal(2) = 0; }
      else if(sy == NUM.j) {normal(0) = 0; normal(1) = 1; normal(2) = 0;}
      #if DIM ==3
      else if(sz == 0) {normal(0) = 0; normal(1) = 0; normal(2) = -1;}
      else if(sz == NUM.k) {normal(0) = 0; normal(1) = 0; normal(2) = 1;} //2026.2.18
      #endif //2026.2.18
      b_particle->norm = normal;

      b_particle->rho_0        = ini_rho[0];//2026.2.18
      SetZero(b_particle->rhov_0);//2025.12.25
      SetZero(b_particle->tv);//
      for (int si = 0; si < 2; si ++)
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
  


////////////////////////////////////////////////////////////
  #if DIM ==2
  parallel_for( blocked_range2d<int>(0, NUM.i, 0, NUM.j), [&](const blocked_range2d<int>& r){
    for(int i= r.rows().begin(); i< r.rows().end(); i++)
      for(int j = r.cols().begin(); j < r.cols().end(); j++)
      { //2026.1.12
        Ele_cen[i][j]->neigh_cell.clear(); 
        for(int t=AMAX1(i- 1, 0); t<=AMIN1(i+ 1, NUM.i-1); t++)
        for(int s=AMAX1(j- 1, 0); s<=AMIN1(j+ 1, NUM.j-1); s++)
        Ele_cen[i][j]->neigh_cell.push_back(Ele_cen[t][s]) ; 
    }
  }, ap);
  #elif DIM ==3
  parallel_for( blocked_range3d<int>(0, NUM.i, 0, NUM.j, 0, NUM.k), [&](const blocked_range3d<int>& r){
    for(int i= r.pages().begin(); i< r.pages().end(); i++)
      for(int j = r.rows().begin(); j < r.rows().end(); j++)
      for(int k = r.cols().begin(); k < r.cols().end(); k++)
      { //2026.1.12
        Ele_cen[i][j][k]->neigh_cell.clear(); 
        for(int t=AMAX1(i- 1, 0); t<=AMIN1(i+ 1, NUM.i-1); t++)
        for(int s=AMAX1(j- 1, 0); s<=AMIN1(j+ 1, NUM.j-1); s++)
        for(int q=AMAX1(k- 1, 0); q<=AMIN1(k+ 1, NUM.k-1); q++)
        Ele_cen[i][j][k]->neigh_cell.push_back(Ele_cen[t][s][q]) ; 
    }
  }, ap);
  #endif //2026.7.24



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


  Reset_cell_list_info(world);
  Map_the_particle_to_tree(world); 

  // // /**********************************************************/
  // Split_particle(world);
  // Reset_cell_list_info(world);
  // Map_the_particle_to_tree(world); 
  // Reset_refresh_merge_neighbor(world);

  Split_particle_fast(world);

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
  time_for_refresh_neighbor  += time_refresh_neighbor.elapsed(); //2026.7.3

  while(md< 2) 
  {
    time_reset_neighbor.restart();
    Reset_particle_neighbor_info(world, 2);//
    time_for_reset_neighbor += time_reset_neighbor.elapsed();

    time_refresh_neighbor.restart();
    Refresh_neighbor_info(2, world);
    time_for_refresh_neighbor  += time_refresh_neighbor.elapsed(); //2026.6.30
  
    ///////////////////////////////


    if(md==0)
    MLS_P2G(world);


    Calculate_limiter(world);//


    Grid_acceleration( world);


    Update_grid(world, md);


    Update_particle(world, md);


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

  Real dw =0, jet= domain.i, ups= domain.i;
  Real mass1=0, mass2=0;
  int ss1=0, ss2=0; //2026.7.16
	for (int i = 0; i < int(particle_fluid.size()); i++)
  {
    p_Particle cp = particle_fluid[i]; 
    if(cp->phase == 2)
    {
      if(cp->coord(1) < 2.5 * p_sz && cp->coord(2) < 2.5 * p_sz ) //2025.7.24 
      {
        if(cp->coord(0) > dw ) dw = cp->coord(0);

        if(cp->coord(0) < jet ) jet = cp->coord(0);
      }
      else if(cp->coord(1) > 4* p_sz || cp->coord(2) > 4* p_sz ) // 
      {
        if(cp->coord(0) < ups ) ups = cp->coord(0);

      }
      mass1 += cp->rho[0]* cp->vol ;//2026.7.16
      ss1 ++;
    }


  }


if(Iden==1){
 ofstream out("./pointa.txt",ios::trunc);  
  out<<time1<<" "<<setprecision(6)<<dw<<" "<<jet<<" "<<ups<<"\n";
  out.close();

  ofstream out1("./mass.txt",ios::trunc);  
  out1<<time1<<" "<<setprecision(6)<<mass1<<" "<<int(particle_fluid.size())<<"\n";
  out1.close();
}else{
 ofstream out("./pointa.txt",ios::out|ios::ate|ios::app);  
  out<<time1<<" "<<setprecision(6)<<dw<<" "<<jet<<" "<<ups<<"\n";
  out.close();

  ofstream out1("./mass.txt",ios::out|ios::ate|ios::app);  
  out1<<time1<<" "<<setprecision(6)<<mass1<<" "<<int(particle_fluid.size())<<"\n";
  out1.close();
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


  ofstream out(filename, ios::binary);
  char  versn [] =  "#!TDV112";
  // Write_file_c(out, versn, sizeof(versn));
  out.write((char*)&(versn), sizeof(versn)-1);
  Write_file(out, 1);
  Write_file(out, 0);
  Write_file(out, 0);

  vector<string>  VarName = {"x","y", "z", "phase",  "vx", "vy", "vol", "rho", "P","tv","v","n_neighbor"}; 
  vector<int>  VarTp = {2,2, 2, 3,  2,2,2,  2,2, 2,2,3 }; 

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


  for(int j=0;j< int(VarTp.size());j++)
	for (int i = 0; i < int(particle_fluid.size()); i++)
  {
    Particle *cp = particle_fluid[i]; //2022.2.11

    if(j==0)
      Write_file(out,  Convert_T(cp->coord(0)));
    else if(j==1)
      Write_file(out,  Convert_T(cp->coord(1) ));
    else if(j==2)
      Write_file(out,  Convert_T(cp->coord(2) ));
    else if(j==3)
      Write_file(out,  Convert_T(cp->phase ));
    else if(j==4)
      Write_file(out,  Convert_T(cp->v(0) ));
    else if(j==5)
      Write_file(out,  Convert_T(cp->v(1) ));
    else if(j==6)
      Write_file(out,  Convert_T(cp->vol ));
    else if(j==7)
      Write_file(out,  Convert_T(cp->rho[0] ));
    else if(j==8)
      Write_file(out,  Convert_T(cp->P ));
    else if(j==9)
      Write_file(out,  Convert_T(cp->tv.head(DIM).norm() ) );
    else if(j==10)
      Write_file(out,  Convert_T(cp->v.head(DIM).norm() ));
    else if(j==11)
      Write_file(out,  Convert_T(int(cp->neighbor.size()) ));

  }
  out.close();





  char    filename1[256];
  sprintf(filename1,"%s%d%s%d%s","./outdata/ele_file.",n,".",world.rank(),".plt");

  ofstream out1(filename1, ios::trunc);
  out1<<"VARIABLES = \"x\",\"y\",\"z\", \"v1\",\"v2\", \"v\",\"rho0\",\"rho1\", \
       \"E0\",\"E1\",\"Lim0\",\"Lim1\",\"st\", \"st1\", \"P\"\n";
  out1<<"ZONE T="<<"\"GRID 1\","<<"SOLUTIONTIME="<<run_time<<"\n";

  for (int i = 0; i < int(NUM.i); i++)
  for (int j = 0; j < int(NUM.j); j++)
#if DIM==3
  for (int k = 0; k < int(NUM.k); k++)
#endif
  {
    #if DIM==2
    Particle *cp = Ele_cen[i][j]; //2020.7.6
    #elif DIM==3
    Particle *cp = Ele_cen[i][j][k]; //2020.7.6
    #endif
    
    out1<<cp->coord(0)<<' '<<cp->coord(1)<<' '<<cp->coord(2)<<' '
    <<cp->v(0)<<' '<<cp->v(1)<<' '<<cp->v.head(DIM).norm()<<' '
    <<cp->rho[0]<<' '<<cp->rho[1]<<' '
    <<cp->coef[0](1,1)<<' '<<cp->coef[0](1,2)<<' '
    <<cp->Lim_u[0](0)<<' '<<cp->Lim_u[1](0)<<' '
    <<cp->alph[0]<<' '
    <<cp->alph[1]<<' '
    <<cp->P<<' '
    <<endl;
      
  }
  out1.close();



}

