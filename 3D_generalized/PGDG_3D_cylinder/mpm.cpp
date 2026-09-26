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
// Set_timestep
//--------------------------------------------------
void MPM::Set_timestep(communicator &world)
{
  local_timestep = 1.e20;
  max_v    = -1.e20;

  static affinity_partitioner ap;

  auto& particles = (iterate_num > 1) ? Used_Ele_cen : particle_fluid;
  auto result = parallel_reduce(blocked_range<int>(0, static_cast<int>(particles.size())),
      pair<Real, Real>{local_timestep, max_v},[&](const auto& range, auto r)
      {
          for (int i = range.begin(); i != range.end(); ++i) {
              Particle* current = particles[i];

              current->Set_timestep(box_l, box_r, this);

              r.first  = AMIN1(r.first,  current->timestep);
              r.second = AMAX1(r.second, current->v.head(DIM).norm());
          }
          return r;
      },

      [](auto a, auto b)
      {
          return pair<Real, Real>{
              AMIN1(a.first,  b.first),
              AMAX1(a.second, b.second)
          };
      }
  );

  local_timestep = result.first;
  max_v          = result.second;

  all_reduce(world, local_timestep, glbl_timestep, mpi::minimum<Real>());
  all_reduce(world, max_v,          glbl_max_v,      mpi::maximum<Real>());


}




//--------------------------------------------------
// reset all the cell list
//--------------------------------------------------
void MPM::Reset_cell_list_info(communicator &world)
{
  for (int i = 0; i < num_level; i++) //2024.7.23
    level_info[i]->Reset_cell_list_info(world);

}
//--------------------------------------------------
//  map the particle 
//--------------------------------------------------
void MPM::Map_the_particle_to_tree(communicator &world)
{
  static affinity_partitioner ap;

  parallel_for( blocked_range<int>(0, int(particle_fluid.size())),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){

      int current_level = particle_fluid[i]->level;
      level_info[current_level]->Update_cell_list(particle_fluid[i], world);
    
    }
  }, ap);  

  parallel_for( blocked_range<int>(0, int(particle_bound.size())),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){

      int current_level = particle_bound[i]->level;
      level_info[current_level]->Update_cell_list(particle_bound[i], world);
    
    }
  }, ap);  

}

//--------------------------------------------------
// reset all the particle neighbor info 2026.7.3
//--------------------------------------------------
void MPM::Reset_particle_neighbor_info(communicator &world, int flag) //2026.7.1
{
  static affinity_partitioner ap;

  if(flag == 2)
  {
    parallel_for( blocked_range<int>(0, int(particle_fluid.size())),[&](const blocked_range<int>& r){
      for(int i=r.begin(); i!=r.end(); ++i){
        particle_fluid[i]->Reset_neighbor_info();
      }
    }, ap);
  
  }

  else if(flag ==4) //2026.7.1
  {
    parallel_for( blocked_range<int>(0, int(particle_fluid.size())),[&](const blocked_range<int>& r){
      for(int i=r.begin(); i!=r.end(); ++i){
        particle_fluid[i]->in_ele.clear();
      }
    }, ap);

    parallel_for( blocked_range<int>(0, int(Used_Ele_cen.size())),[&](const blocked_range<int>& r){
      for(int is=r.begin(); is!=r.end(); ++is){
        Used_Ele_cen[is]->Reset_neighbor_info();
      }
    }, ap);
  }



}

//--------------------------------------------------------
// refresh the neighbor infomation for all the particles 2026.7.3
//-------------------------------------------------------
void MPM::Refresh_neighbor_info(int flag, communicator &world)
{
  static affinity_partitioner ap;
  
  if(flag ==2)
  parallel_for( blocked_range<int>(0, int(particle_fluid.size())),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      particle_fluid[i]->Refresh_neighbor_info(this, flag);
    }
  }, ap);


  else if(flag ==4)
  {
    parallel_for( blocked_range<int>(0, int(Used_Ele_cen.size())),[&](const blocked_range<int>& r){
      for(int is=r.begin(); is!=r.end(); ++is){
          Used_Ele_cen[is]->Refresh_ele_info(this, flag);
      }
      
    }, ap);//2026.7.2

  }


}

//--------------------------------------------------------
// 2026.2.7
//-------------------------------------------------------
void MPM::Reset_refresh_merge_neighbor(communicator &world)
{
  static affinity_partitioner ap;

  parallel_for( blocked_range<int>(0, int(particle_fluid.size())),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      particle_fluid[i]->Reset_neighbor_info();
      particle_fluid[i]->Refresh_neighbor_info(this, 3);
    }
  }, ap);

}



//--------------------------------------------------
// MLS_P2G //2025.12.25
//--------------------------------------------------
void MPM::Used_state_update(communicator &world)
{
  my_int Ib; 
  Ib.i = DIM_X==1 ? 5 : 0;
  Ib.j = DIM_Y==1 ? 5 : 0;
  Ib.k = DIM_Z==1 ? 5 : 0;

  static affinity_partitioner ap;

  #if DIM ==2
  parallel_for( blocked_range2d<int>(0, NUM.i, 0, NUM.j), [&](const blocked_range2d<int>& r){
    for(int i= r.rows().begin(); i< r.rows().end(); i++)
      for(int j = r.cols().begin(); j < r.cols().end(); j++){ //2026.1.12
        Ele_cen[i][j]->used_st = 0; 
    }
  }, ap);
  parallel_for( blocked_range2d<int>(0, NUM.i, 0, NUM.j), [&](const blocked_range2d<int>& r){
    for(int i= r.rows().begin(); i< r.rows().end(); i++)
      for(int j = r.cols().begin(); j < r.cols().end(); j++){ //2026.1.12

        if( Ele_cen[i][j]->ph_state[0] == 1 || Ele_cen[i][j]->ph_state[1] == 1) //2026.2.18
        {
          for(int t=AMAX1(i- Ib.i, 0); t<=AMIN1(i+ Ib.i, NUM.i-1); t++)
          for(int s=AMAX1(j- Ib.j, 0); s<=AMIN1(j+ Ib.j, NUM.j-1); s++)
            Ele_cen[t][s]->used_st = 1; 
        }
    }
  }, ap);

  Used_Ele_cen.clear();

  parallel_for( blocked_range2d<int>(0, NUM.i, 0, NUM.j), [&](const blocked_range2d<int>& r){
    for(int i= r.rows().begin(); i< r.rows().end(); i++)
      for(int j = r.cols().begin(); j < r.cols().end(); j++){ //2026.1.12
        if(Ele_cen[i][j]->used_st == 1)
          Used_Ele_cen.push_back(Ele_cen[i][j]);
        
    }
  }, ap);

  #elif DIM ==3
  parallel_for( blocked_range3d<int>(0, NUM.i, 0, NUM.j, 0, NUM.k), [&](const blocked_range3d<int>& r){
    for(int i= r.pages().begin(); i< r.pages().end(); i++)
      for(int j = r.rows().begin(); j < r.rows().end(); j++)
      for(int k = r.cols().begin(); k < r.cols().end(); k++)
      { //2026.1.12
        Ele_cen[i][j][k]->used_st = 0; 
    }
  }, ap);

  parallel_for( blocked_range3d<int>(0, NUM.i, 0, NUM.j, 0, NUM.k), [&](const blocked_range3d<int>& r){
    for(int i= r.pages().begin(); i< r.pages().end(); i++)
    for(int j = r.rows().begin(); j < r.rows().end(); j++)
    for(int k = r.cols().begin(); k < r.cols().end(); k++)
    { //2026.1.12
      if( Ele_cen[i][j][k]->ph_state[0] == 1 || Ele_cen[i][j][k]->ph_state[1] == 1) //2026.2.18
      {
        for(int t=AMAX1(i- Ib.i, 0); t<=AMIN1(i+ Ib.i, NUM.i-1); t++)
        for(int s=AMAX1(j- Ib.j, 0); s<=AMIN1(j+ Ib.j, NUM.j-1); s++)
        for(int q=AMAX1(k- Ib.k, 0); q<=AMIN1(k+ Ib.k, NUM.k-1); q++)
          Ele_cen[t][s][q]->used_st = 1; 
      }
    }
  }, ap);

  Used_Ele_cen.clear();

  parallel_for( blocked_range3d<int>(0, NUM.i, 0, NUM.j, 0, NUM.k), [&](const blocked_range3d<int>& r){
    for(int i= r.pages().begin(); i< r.pages().end(); i++)
    for(int j = r.rows().begin(); j < r.rows().end(); j++)
    for(int k = r.cols().begin(); k < r.cols().end(); k++)
    { //2026.1.12
      if(Ele_cen[i][j][k]->used_st == 1)
        Used_Ele_cen.push_back(Ele_cen[i][j][k]);
    }
  }, ap);

  #endif //2026.2.18

}







//--------------------------------------------------
// MLS_P2G //2025.12.25
//--------------------------------------------------
void MPM::MLS_P2G(communicator &world)
{
  static affinity_partitioner ap;

  parallel_for( blocked_range<int>(0, int(Used_Ele_cen.size())),[&](const blocked_range<int>& r){
    for(int is=r.begin(); is!=r.end(); ++is){
      p_Particle cp = Used_Ele_cen[is]; 
      Vector<Real, 3> coordi = cp->coord;

      for (int t = 0; t < int(cp->neighbor.size()); t++)
      {
        cp->neighbor[t]->surf1=0;
        cp->neighbor[t]->surf=0;//2026.1.8
      }


      Real mass_t[FldNum]; 
      Vector<Real, UNum > p_avg[FldNum]; 
      Matrix<Real,UNum, BaseNum > Coef[FldNum]; 
      for(int t=0; t< FldNum; t++)
      {
        mass_t[t]=0; 
        SetZero(p_avg[t]);
        SetZero(Coef[t]) ;
      }

      for(int Id_j=0; Id_j< cp->neigh_cell.size(); Id_j++)
      {
        for (int jj = 0; jj < int(cp->neigh_cell[Id_j]->neighbor.size()); jj ++)
        {
          p_Particle neigh = cp->neigh_cell[Id_j]->neighbor[jj];
          Vector<Real, 3> xij = neigh->coord - coordi;
          Real dist = xij.head(DIM).norm();
          if( dist < 1.5 * p_sz)
          // if( dist < 2.2 * p_sz)
          {
            Real hi = p_g*1 * 3.2/CUT_OFF;  //2026.7.3
            Real ww = cp->Kernel_function(dist, hi)*neigh->vol; //2026.2.7

            Vector<Real, UNum > pj; 
            #if DIM==2
            pj<<neigh->rho[0], neigh->rhov[0](0), neigh->rhov[0](1), neigh->E[0];
            #elif DIM==3
            pj<<neigh->rho[0], neigh->rhov[0](0), neigh->rhov[0](1), neigh->rhov[0](2), neigh->E[0];
            #endif 
            

            pj += neigh->coef[0].rightCols(BaseNum)* Basis_f(-xij, p_size).tail(BaseNum);
            int ph_id = neigh->phase -1;
            if(ph_id >= 0 && ph_id < FldNum)
            {
              p_avg[ph_id] += pj*ww;
              Coef[ph_id] += neigh->coef[0].rightCols(BaseNum)* ww;
              mass_t[ph_id] += ww; 
            }


          }

        }
      }

      Real mass_total = 0;
      for(int t=0; t< FldNum; t++) mass_total += mass_t[t] ;

      for(int ss=0;ss < FldNum; ss++)
      {
        if(mass_t[ss] > 0)  cp->ph_state[ss] = 1;
        else  cp->ph_state[ss] = 0;

        cp->alph[ss] = mass_t[ss]/(mass_total + 1e-20);
  
        if( (cp->ph_state[ss] == 1 && cp->ph_state_b[ss] == 0)  ) //|| iterate_num ==1
        {
          SetZero(cp->coef[ss]);
          p_avg[ss] /= (mass_t[ss] + 1e-20);
          Coef[ss] /= (mass_t[ss] + 1e-20);//2026.1.22


          cp->coef[ss](0,0) = cp->rho[ss] = p_avg[ss](0);
          cp->coef[ss].col(0).segment<DIM>(1) = cp->rhov[ss].head(DIM) = p_avg[ss].segment<DIM>(1); //2026.2.18
          cp->coef[ss](UNum-1,0) = cp->E[ss] = p_avg[ss](UNum-1);

          cp->coef[ss].rightCols(BaseNum) = Coef[ss].rightCols(BaseNum);//2026.1.22
          cp->v = cp->rhov[ss]/cp->rho[ss];//2026.1.27

          Real e_j = cp->E[ss]  - 0.5* cp->rhov[ss].head(DIM).squaredNorm()/cp->rho[ss];
          cp->P = EOS_P(e_j, ss);//2026.1.27
          
          cp->coef0[ss] = cp->coef[ss];
        }
        cp->ph_state_b[ss] = cp->ph_state[ss] ;

        //////////////////////////////////////
        if(cp->ph_state[ss] ==0 )
        {
          SetZero(cp->coef[ss]);
          SetZero(cp->coef0[ss]);
          cp->rho[ss] = ini_rho[ss];
          SetZero(cp->rhov[ss]);
        }

      }


    }
  }, ap);


}






//--------------------------------------------------
// Calculate_limiter 2026.1.16
//--------------------------------------------------
void MPM::Calculate_limiter(communicator &world)
{
  static affinity_partitioner ap;

  parallel_for( blocked_range<int>(0, int(Used_Ele_cen.size())),[&](const blocked_range<int>& r){
    for(int is=r.begin(); is!=r.end(); ++is){
      p_Particle cp= Used_Ele_cen[is];
      int i = cp->idI.i; int j = cp->idI.j; int k = cp->idI.k; 
      Cell_dg Celli;
      #if DIM ==2
      Celli = cell_dg[i][j];
      #elif DIM ==3
      Celli = cell_dg[i][j][k];
      #endif //2026.2.18

      Vector<Real, UNum> u_max[FldNum];
      Vector<Real, UNum> u_min[FldNum];
      Vector<Real, UNum> pi[FldNum];
      for (int s = 0; s < FldNum; s ++)
      {
        u_max[s].fill(-1e20);
        u_min[s].fill(1e20);
        pi[s].setZero();//2026.1.27
      }

      Real entr_i[FldNum]; //20260521
      for (int s = 0; s < FldNum; s ++)
      {
        cp->Lim_u[s].fill(1.); //2026.6.24
        Real e_i = cp->E[s]  - 0.5* cp->rhov[s].head(DIM).squaredNorm()/cp->rho[s];
        Real P_i = EOS_P(e_i, s);//
        entr_i[s] = (P_i + Pp0[s])/pow(cp->rho[s], gamma[s]);//20260521

        #if DIM==2
        pi[s]<< cp->rho[s], cp->rhov[s](0), cp->rhov[s](1), cp->E[s]; 
        #elif DIM==3
        pi[s]<< cp->rho[s], cp->rhov[s](0), cp->rhov[s](1),cp->rhov[s](2), cp->E[s]; 
        #endif 
      }

      cp->G_v.fill(0.);
      Vector<Real, 3> vs; vs.fill(0);
      int bd_fg = 0; //2026.7.23

#if DIM ==2
      if(i== NUM.i-1)
      {
        // vs(0) = 0; vs(1) = 0; 
        vs = Ele_cen[i][j]->v;//20260618
        cp->G_v.col(0) = (vs - Ele_cen[i][j]->v).head(DIM)/(0.5* p_size.i);
        bd_fg = 1;
      }
      else if(i == 0)
      {
        // vs(0) = 0; vs(1) = 0; 
        vs = Ele_cen[i][j]->v;//20260618
        cp->G_v.col(0) = ( Ele_cen[i][j]->v - vs).head(DIM)/(0.5*p_size.i); 
        bd_fg = 1;
      }
      else
      cp->G_v.col(0) = (Ele_cen[i+1][j]->v - Ele_cen[i-1][j]->v).head(DIM)/(2*p_size.i);

      if(j == NUM.j-1)
      {
        vs(0) = 0; vs(1) = 0; 
        cp->G_v.col(1) = (vs - Ele_cen[i][j]->v).head(DIM)/(0.5*p_size.j);
        bd_fg = 1;
      }
      else if(j == 0)
      {
        vs(0) = 0; vs(1) = 0 ; 
        cp->G_v.col(1) = ( Ele_cen[i][j]->v - vs).head(DIM)/(0.5*p_size.j); 
        bd_fg = 1;
      }
      else
        cp->G_v.col(1) = (Ele_cen[i][j+1]->v - Ele_cen[i][j-1]->v).head(DIM)/(2*p_size.j);

#elif DIM ==3

      if(i== NUM.i-1)
      {
        vs = Ele_cen[i][j][k]->v;//20260618
        cp->G_v.col(0) = (vs - Ele_cen[i][j][k]->v).head(DIM)/(0.5* p_size.i);
        bd_fg = 1;
      }
      else if(i == 0)
      {
        vs = Ele_cen[i][j][k]->v;//20260618
        cp->G_v.col(0) = ( Ele_cen[i][j][k]->v - vs).head(DIM)/(0.5* p_size.i);
        bd_fg = 1;
      }
      else
      {
        Vector<Real, 3> vR = Ele_cen[i+1][j][k]->v;
        Vector<Real, 3> vL = Ele_cen[i-1][j][k]->v;

        cp->G_v.col(0) = (vR - vL).head(DIM)/(2*p_size.i);
      }


      if(j == NUM.j-1)
      {
        vs(0) = 0; vs(1) = 0; vs(2) = 0; 
        cp->G_v.col(1) = (vs - Ele_cen[i][j][k]->v).head(DIM)/(0.5*p_size.j);
        bd_fg = 1;
      }
      else if(j == 0)
      {
        // vs(0) = Ele_cen[i][j][k]->v(0); vs(1) = 0; vs(2) = Ele_cen[i][j][k]->v(2); 
        vs(0) = 0; vs(1) = 0; vs(2) = 0; 
        cp->G_v.col(1) = (Ele_cen[i][j][k]->v - vs).head(DIM)/(0.5*p_size.j);
        bd_fg = 1;
      }
      else
      {
        Vector<Real, 3> vR = Ele_cen[i][j+1][k]->v;
        Vector<Real, 3> vL = Ele_cen[i][j-1][k]->v;

        cp->G_v.col(1) = (vR - vL).head(DIM)/(2*p_size.j);
      }


      if(k == NUM.k-1)
      {
        vs(0) = 0; vs(1) = 0; vs(2) = 0; 
        cp->G_v.col(2) = (vs - Ele_cen[i][j][k]->v).head(DIM)/(0.5*p_size.k);
        bd_fg = 1;
      }
      else if(k == 0)
      {
        // vs(0) = Ele_cen[i][j][k]->v(0); vs(1) = Ele_cen[i][j][k]->v(1); vs(2) = 0; 
        vs(0) = 0; vs(1) = 0; vs(2) = 0; 
        cp->G_v.col(2) = ( Ele_cen[i][j][k]->v - vs).head(DIM)/(0.5 *p_size.k); 
        bd_fg = 1;
      }
      else
      {
        Vector<Real, 3> vR =  Ele_cen[i][j][k+1]->v;
        Vector<Real, 3> vL =  Ele_cen[i][j][k-1]->v;

        cp->G_v.col(2) = (vR - vL).head(DIM)/(2*p_size.k);
      }

      #endif //2026.2.18



      for (int jj = 0; jj < int(Celli.cell_id.size()); jj ++) //2026.2.18
      {

        my_int cell_j = Celli.cell_id[jj];
        int b_mk = Celli.bound[jj];
        Vector<Real, UNum> pij_i[FldNum]; //2026.1.27
        for (int s = 0; s < FldNum; s ++) pij_i[s].setZero();
        Vector<Real, UNum> pj; pj.setZero();

        Vector<Real, 3> fcoord_j = Celli.face_center[jj];
        Vector<Real, 3> dr =  fcoord_j - cp->coord;
        // Vector<Real, 3> drc =  cell_dg[cell_j.i][cell_j.j].cell_center - cp->coord;//2026.1.27
        Vector<Real, 3> normal = Celli.face_normal[jj];
        Real Sj = Celli.face_S[jj];


        if(b_mk == 0)
        {
          p_Particle neigh;
          #if DIM ==2
          neigh = Ele_cen[cell_j.i][cell_j.j];
          #elif DIM ==3
          neigh = Ele_cen[cell_j.i][cell_j.j][cell_j.k];
          #endif //2026.2.18
  
          for (int s = 0; s < FldNum; s ++)
          {

            pj.setZero();
            for(int s1=0; s1<FldNum; s1++)
            {
              if(neigh->alph[s1] == 0) continue;

              Real rhoj = neigh->rho[s1];
              Vector<Real, DIM> rhov_j = neigh->rhov[s1].head(DIM);
              Real Ej = neigh->E[s1];

              if(s1 != s)
              {
                Vector<Real, DIM> v_j = rhov_j/rhoj ;
                Real rhoe_j = Ej - 0.5* rhoj* v_j.squaredNorm() ;
                Real P_j = EOS_P(rhoe_j, s1);

                rhoj = cp->rho[s];
                rhov_j = rhoj* v_j;
                Ej = EOS_rhoe(P_j, s) + 0.5* rhoj* v_j.squaredNorm() ;
              }
              Vector<Real, UNum> tmp_p; 
              tmp_p(0) = rhoj;
              tmp_p.segment<DIM>(1) = rhov_j;//2026.2.18
              tmp_p(UNum-1) = Ej;

              pj += tmp_p* neigh->alph[s1];
            }
            pij_i[s] =  pj -  pi[s]; 
            for(int kk=0; kk<UNum; kk++)
            {
              u_max[s](kk) = AMAX1(u_max[s](kk), pij_i[s](kk)); 
              u_min[s](kk) = AMIN1(u_min[s](kk), pij_i[s](kk)); 
            }

          }
        }
        else
        {
          for (int s = 0; s < FldNum; s ++)
          {
            Vector<Real,BaseNum> Var = Basis_f(dr, p_size).tail(BaseNum);
            Real rhoj = cp->rho[s]  ;//
            Vector<Real, DIM> rhov_j  = cp->rhov[s].head(DIM) ;//
            Real Ej = cp->E[s]  ;//
            #if DIM==2
            pj<< rhoj, rhov_j(0), rhov_j(1), Ej; 
            #elif DIM==3
            pj<< rhoj, rhov_j(0), rhov_j(1), rhov_j(2), Ej; 
            #endif 

            if(b_mk == -1 || b_mk == -2 ) pj = pi[s]; //20260620
            if(j == NUM.j-1 || j == 0) pj.segment<DIM>(1)(1) = -rhov_j(1);  //-rhov_j(1)
            #if DIM ==3
            if(k == NUM.k-1 || k == 0) pj.segment<DIM>(1)(2) = -rhov_j(2);  //-rhov_j(2)
            #endif //2026.7.25

            pij_i[s] = pj -  pi[s]; 
            for(int kk=0; kk<UNum; kk++)
            {
              u_max[s](kk) = AMAX1(u_max[s](kk), pij_i[s](kk)); 
              u_min[s](kk) = AMIN1(u_min[s](kk), pij_i[s](kk)); 
            }
          }
        }


      }

      

      for (int s = 0; s < FldNum; s ++)
      for(int kk=0; kk<UNum; kk++)
      {
        u_max[s](kk) = AMAX1(u_max[s](kk), 0.); 
        u_min[s](kk) = AMIN1(u_min[s](kk), 0.); 
      }


      for (int jj = 0; jj < int(Celli.face_center.size()); jj ++)
      {
        Vector<Real, 3> fcoord_j = Celli.face_center[jj];
        my_int cell_j = Celli.cell_id[jj];
        Vector<Real, 3> dr =  fcoord_j - cp->coord; 
        if(Celli.bound[jj] == 0) 
        {
          #if DIM ==2
          dr = cell_dg[cell_j.i][cell_j.j].cell_center - cp->coord;
          #elif DIM ==3
          dr = cell_dg[cell_j.i][cell_j.j][cell_j.k].cell_center - cp->coord;
          #endif //2026.2.18
        }
      

        for (int s = 0; s < FldNum; s ++)
        {
          Vector<Real,BaseNum> Var = Basis_f(dr, p_size).tail(BaseNum);
          Vector<Real, UNum> dU = cp->coef[s].rightCols(BaseNum)* Var ;
          int lim_fg = bd_fg; 
          if(cp->alph[s] != 1) lim_fg = 1;//2026.7.23

          for(int kk=0; kk < UNum; kk++)
          {
            Real sd=1;


            if( dU(kk) >= 0 )
            {
              sd = u_max[s](kk)/(dU(kk)+ 1e-20);
              sd = Limiter(abs(sd), p_g, lim_fg );
            }
            else if( dU(kk) < 0 )
            {
              sd = u_min[s](kk)/(dU(kk)+ 1e-20);
              sd = Limiter(abs(sd), p_g, lim_fg );//2026.7.23
            }


            cp->Lim_u[s](kk) = AMIN1(cp->Lim_u[s](kk), sd); 
          }
        }

      }

      for (int s = 0; s < FldNum; s ++)
      {
        int fg =0;

        for (int jj = 0; jj < 4; jj ++)
        {
          Vector<Real, 3> dr; dr.setZero();
          Real theta = (jj + 0.5)*PI*0.5 ;
          dr(0) = sqrt(2.)*p_g*cos( theta); dr(1) = sqrt(2.)*p_g*sin( theta); 

          #if DIM ==3
          for(int k1=0; k1<2; k1++)
          #endif //2026.2.18
          {
          #if DIM ==3
            dr(2) = (-0.5 + k1)* p_g;
          #endif //2026.2.18

            Vector<Real,BaseNum> Var = Basis_f(dr, p_size).tail(BaseNum);
            Vector<Real, UNum> Ui = cp->coef[s].col(0);
            for(int kk=0; kk < UNum; kk++)  Ui(kk) += cp->Lim_u[s](kk)* cp->coef[s].row(kk).rightCols(BaseNum)* Var;

            Real rho_p =  Ui(0);
            Real E_p =  Ui(UNum-1);//2026.7.25
            Vector<Real, DIM> rhov_p; SetZero(rhov_p);
            rhov_p.head(DIM) = Ui.segment<DIM>(1);
            Real e_p = E_p - 0.5* rhov_p.head(DIM).squaredNorm()/rho_p;
            Real P_p = EOS_P(e_p, s);
            if(P_p < -0.05* Pp0[s]  || rho_p<0) {fg =1; goto Lim_det;}
          }


        } //2026.6.26


        Lim_det:
        {
          if(fg ==1) cp->Lim_u[s].fill(0.); //2026.6.25
        }

      }


    }
  }, ap); 



}





//--------------------------------------------------
// Grid_acceleration 2025.12.25
//--------------------------------------------------
void MPM::Grid_acceleration(communicator &world)
{
  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, int(Used_Ele_cen.size())),[&](const blocked_range<int>& r){
    for(int is=r.begin(); is!=r.end(); ++is){
      p_Particle cp = Used_Ele_cen[is]; 
      
      for (int s = 0; s < FldNum; s ++)
      {
        if(cp->ph_state[s] == 1)
        cp->Grid_acceleration(this, s);
      }

    }
  }, ap); 

}

//--------------------------------------------------
// Update_grid 2025.12.26
//--------------------------------------------------
void MPM::Update_grid(communicator &world, int md)
{
  
  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, int(Used_Ele_cen.size())),[&](const blocked_range<int>& r){
    for(int is=r.begin(); is!=r.end(); ++is){
      p_Particle cp = Used_Ele_cen[is]; 

      Vector<Real, 3> rhov_a; rhov_a.setZero();
      Real rho_a = 0;
      Real P_a = 0;

      for (int s = 0; s < FldNum ; s ++)
      {
        cp->c[s] =0;//2026.2.9
        cp->ac[s] =0;//2026.2.9

        for(int t=0;t<UNum;t++) 
        {
          cp->dcoef[s].row(t).rightCols(BaseNum) *= abs(cp->Lim_u[s](t));
          cp->coef[s].row(t).rightCols(BaseNum) *= abs(cp->Lim_u[s](t));

        }

        if(cp->ph_state[s] ==1)
        {

          if(md==0)
          {
            cp->coef0[s] = cp->coef[s];
            cp->coef[s] += cp->dcoef[s]* glbl_timestep;
          }
          else if(md==1)
          {
            cp->coef[s] = 0.5*( cp->dcoef[s]* glbl_timestep  +  cp->coef[s] + cp->coef0[s] );
          }



          cp->rho[s]     = cp->coef[s](0,0) ;
          cp->rhov[s].head(DIM) = cp->coef[s].col(0).segment<DIM>(1); //2026.2.18
          cp->E[s]       = cp->coef[s](UNum-1,0) ;
          Real e_j = cp->E[s]  - 0.5* cp->rhov[s].head(DIM).squaredNorm()/cp->rho[s];

          cp->P = EOS_P(e_j, s);//2026.1.27
          cp->c[s]  =  sqrt( gamma[s]* (Pp0[s] + cp->P)/cp->rho[s] ) ; //2026.1.27
          cp->a.head(DIM) =  (cp->dcoef[s].col(0).segment<DIM>(1) - cp->dcoef[s](0,0)* cp->v.head(DIM))/cp->rho[s] ; //2026.2.9
          cp->ac[s] = cp->a.head(DIM).norm(); 

        }

        rhov_a += cp->rhov[s]* cp->alph[s];
        rho_a += cp->rho[s]* cp->alph[s];
        P_a += cp->P * cp->alph[s];
      }
      cp->v = rhov_a/rho_a;
      cp->P = P_a ;

    }
  }, ap); 


}





//--------------------------------------------------
// G2P_interp 2025.12.26
//--------------------------------------------------
void MPM::Update_particle(communicator &world, int md)
{
  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, int(particle_fluid.size())),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      p_Particle cp= particle_fluid[i];
      p_Particle cell_i = cp->in_ele[0];

      cp->Update_state_onestep(cell_i, md, this); 
      cp->Update_PST(this); //2026.6.29
    }
  }, ap); 


}

  
  //-----------------------------------------------------
// EOS_P
//-----------------------------------------------------
Real MPM::EOS_P(Real ee, int s)
{
  Real Pp=0;
  Pp = ee* ( gamma[s] - 1) - gamma[s]* Pp0[s];

  return Pp;
}

//-----------------------------------------------------
// EOS_rho
//-----------------------------------------------------
Real MPM::EOS_rhoe(Real p1, int s)
{
  Real r1=0;
  r1 =  (p1 + gamma[s]* Pp0[s]) /(gamma[s] - 1); //2026.1.27

  return r1;
}


//-----------------------------------------------------
// MGFM
//-----------------------------------------------------
Vector<Real, 3>  MPM::MGFM(Real ul, Real ur, Real pl, Real pr, Real rhol, Real rhor, int s, int s1)
{

  pl = AMAX1(pl, 1e-10);
  pr = AMAX1(pr, 1e-10);

  rhol = AMAX1(rhol, 1e-10);
  rhor = AMAX1(rhor, 1e-10);


  Real tol = AMAX1(1e-8, 1e-3 * abs(ul - ur) ) ; //2026.6.22
  Real pstar = AMAX1(pl, pr) ; //  (pl + pr)* 0.5
  // Real pff = pstar;
  Real f1 = 1.;
  Real dx = 1e-6* AMAX1(abs(pstar), 1.);

  Real rhols=0, rhors=0, usl=0, Wl =0, Wr =0;  
  int nn=0;

  do{
     Wl = Wl_cal(rhol, pl, pstar, gamma[s], Pp0[s]);
     Wr = Wl_cal(rhor, pr, pstar, gamma[s1], Pp0[s1]);

    usl = ul - (pstar -pl)/Wl;
    rhols = RhoLs_cal(rhol, pl, pstar, gamma[s], Pp0[s]);

    f1 = ul - ur - (pstar -pl)/Wl - (pstar -pr)/Wr;

    Real Jacobi = - Jaco_cal(Wl, pl, rhol, pstar,  gamma[s], Pp0[s]) - Jaco_cal(Wr, pr, rhor, pstar, gamma[s1], Pp0[s1]) ;

    Real dp = -f1/Jacobi;
    if(abs(dp) > 0.6 *pstar)  dp = copysign(0.6 *pstar,dp);
    pstar += dp;


    pstar = AMAX1(pstar, 1e-8);

    nn++;

  }while(abs(f1) > tol && nn<50);
  Vector<Real, 3> val; val<< rhols, pstar, usl;


  return val ;
}



//-------------------------------------------------------
// Split_particle
//-------------------------------------------------------
void MPM::Split_particle(communicator &world)
{
  concurrent_vector  <p_Particle>  particle_tmp;
  particle_tmp.clear();

  Real vol_max = pow(p_sz, DIM)* 1.6 ;

  static affinity_partitioner ap;
  parallel_for( blocked_range<int>(0, int(particle_fluid.size())),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
	 	  p_Particle cp = particle_fluid[i];
		   cp->copy_p.clear();
		   cp->spl=0;
		   cp->mrg_mark = 0;
       Real dp = pow(cp->vol,1./DIM);//2026.7.24

       bool f1 = cp->vol > vol_max;
       bool f2 = cp->coord(0) - box_l.i>= dp && cp->mark_in ==1 ;
       bool f3 = cp->coord(0) > box_r.i  ;

      if(f1 || f2 || f3)
      {

        if( f1 )
        {
          cp->spl =1;
          
          for(int is=0; is< 2;is++){
            for(int js=0; js< 2 ;js++){
            #if DIM ==3
            for(int ks=0; ks< 2 ; ks++) //2026.7.24
            #endif
            {
              p_Particle partcl= particlepool.malloc();

              partcl->coord.fill(0);
              partcl->coord(0) = cp->coord(0) - dp/4 + dp/2*is; 
              partcl->coord(1) = cp->coord(1) - dp/4 + dp/2*js; 
              #if DIM ==3
              partcl->coord(2) = cp->coord(2) - dp/4 + dp/2*ks; 
              #endif
              partcl->Set_split_ptcl(cp, this);


              particle_tmp.push_back(partcl); 
            }
            }
          }
          cp->Clearup();
          particlepool.free(cp);
        }
        else if( f2) //inlet 20260618
        {
          p_Particle partcl= particlepool.malloc();
      
          partcl->coord.fill(0);
          partcl->coord = cp->coord ; 
          partcl->coord(0) = cp->coord(0) - dp; //2026.7.24
          partcl->Set_ptcl(cp, this);

          particle_tmp.push_back(cp); 
          particle_tmp.push_back(partcl); 
        }
        else if( f3 ) //outlet
        {
          cp->Clearup();
          particlepool.free(cp);
        }
      }
	 	  else
	 		  particle_tmp.push_back(cp); 
	   }
   }, ap); 



  if(particle_tmp.size()>0)
  {
	  particle_fluid.clear();

	  parallel_for( blocked_range<int>(0, particle_tmp.size()),
			   [&](const blocked_range<int>& r){
		for(int i=r.begin(); i!=r.end(); ++i){
      p_Particle cp = particle_tmp[i];
			particle_fluid.push_back(cp);

      cp->mark_in =0;
      if(cp->coord(0) - box_l.i < pow(cp->vol,1./DIM) ) cp->mark_in = 1; //2026.7.24

		}
	  }, ap);
  }



  particle_tmp.clear();
}


//-----------------------------------------------------
// Split_particle_fast
//-----------------------------------------------------
void MPM::Split_particle_fast(communicator &world)
{
  int child_N = 4;
  #if DIM == 3
  child_N = 8;
  #endif

  struct FastSplitRecord
  {
    p_Particle parent;
    my_int old_cell;
    #if DIM == 2
    p_Particle child[4-1];
    #elif DIM == 3
    p_Particle child[8-1];
    #endif
    vector<p_Particle> old_neighbors;
  };

  struct InletRecord
  {
    p_Particle particle;
  };

  struct OutletRecord
  {
    p_Particle particle;
    my_int old_cell;
    int level_index;
    vector<p_Particle> old_neighbors;
  };

  Real vol_max = powern(p_sz, DIM)*1.6;

  static affinity_partitioner ap;

  concurrent_vector<p_Particle> new_particles;
  concurrent_vector<FastSplitRecord> split_records;
  concurrent_vector<InletRecord> inlet_records;
  concurrent_vector<OutletRecord> outlet_records;
  unordered_set<p_Particle> outlet_set;

  //-----------------------------------------------------
  // 1. split / inlet / outlet
  parallel_for(blocked_range<int>(0, int(particle_fluid.size())), [&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i)
    {
      p_Particle cp = particle_fluid[i];
      cp->copy_p.clear();
      cp->spl = 0;
      cp->mrg_mark = 0;

      Real dp = pow(cp->vol, 1./DIM);

      bool f1 = cp->vol > vol_max;
      bool f2 = cp->coord(0) - box_l.i >= dp && cp->mark_in == 1;
      bool f3 = cp->coord(0) > box_r.i;

      //---------------------------------------------------
      // split
      if(f1)
      {
        cp->spl = 1;

        int level_index = cp->level - Lmin;
        p_Level_info current_level = level_info[level_index];

        my_real old_shift = my_minus_data(matrx2my_real(cp->coord), box_l);
        my_int old_cell = get_cell_id(old_shift, current_level->dcell, current_level->cell_start, current_level->cell_end);

        my_real p_ref;
        p_ref.i = cp->coord(0) - dp/4.;
        p_ref.j = cp->coord(1) - dp/4.;
        #if DIM == 3
        p_ref.k = cp->coord(2) - dp/4.;
        #endif

        #if DIM == 2
        p_Particle child[4];
        #elif DIM == 3
        p_Particle child[8];
        #endif

        int child_index = 0;

        for(int child_i=0; child_i<2; ++child_i)
        for(int child_j=0; child_j<2; ++child_j)
        #if DIM == 3
        for(int child_k=0; child_k<2; ++child_k)
        #endif
        {
          p_Particle partcl = particlepool.malloc();

          Vector<Real,3> rr; rr.setZero();
          rr(0) = p_ref.i + dp/2.*child_i;
          rr(1) = p_ref.j + dp/2.*child_j;
          #if DIM == 3
          rr(2) = p_ref.k + dp/2.*child_k;
          #endif

          partcl->coord = rr;
          partcl->Set_split_ptcl(cp, this);
          child[child_index++] = partcl;
        }

        vector<p_Particle> old_neighbors;
        old_neighbors.reserve(cp->neighbor.size());
        for(int j=0; j<int(cp->neighbor.size()); ++j) old_neighbors.push_back(cp->neighbor[j]);

        cp->Copy_split_state(child[child_N-1]);
        child[child_N-1]->Clearup();
        particlepool.free(child[child_N-1]);

        FastSplitRecord record;
        record.parent = cp;
        record.old_cell = old_cell;
        record.old_neighbors.swap(old_neighbors);

        for(int j=0; j<child_N-1; ++j)
        {
          record.child[j] = child[j];
          new_particles.push_back(child[j]);
        }

        split_records.push_back(record);
      }

      //---------------------------------------------------
      // inlet
      else if(f2)
      {
        p_Particle partcl = particlepool.malloc();

        partcl->coord.fill(0);
        partcl->coord = cp->coord;
        partcl->coord(0) = cp->coord(0) - dp;
        partcl->Set_ptcl(cp, this);

        InletRecord record;
        record.particle = partcl;
        inlet_records.push_back(record);

        new_particles.push_back(partcl);
      }

      //---------------------------------------------------
      // outlet
      else if(f3)
      {
        int level_index = cp->level - Lmin;
        p_Level_info current_level = level_info[level_index];

        my_real old_shift = my_minus_data(matrx2my_real(cp->coord), box_l);
        my_int old_cell = get_cell_id(old_shift, current_level->dcell, current_level->cell_start, current_level->cell_end);

        OutletRecord record;
        record.particle = cp;
        record.old_cell = old_cell;
        record.level_index = level_index;
        record.old_neighbors.reserve(cp->neighbor.size());
        for(int j=0; j<int(cp->neighbor.size()); ++j) record.old_neighbors.push_back(cp->neighbor[j]);

        outlet_records.push_back(record);
      }
    }
  }, ap);

  for(int i=0; i<int(outlet_records.size()); ++i) outlet_set.insert(outlet_records[i].particle);

  //-----------------------------------------------------
  // 2. new particles -> particle_fluid
  for(int i=0; i<int(new_particles.size()); ++i) particle_fluid.push_back(new_particles[i]);


  // 3. update split cell list
  for(int i=0; i<int(split_records.size()); ++i)
  {
    FastSplitRecord& record = split_records[i];

    #if DIM == 2
    p_Particle split_particles[4];
    #elif DIM == 3
    p_Particle split_particles[8];
    #endif

    split_particles[0] = record.parent;
    for(int j=0; j<child_N-1; ++j) split_particles[j+1] = record.child[j];

    int level_index = record.parent->level - Lmin;
    p_Level_info current_level = level_info[level_index];
    current_level->table_cell_list[record.old_cell.i][record.old_cell.j][record.old_cell.k]->Remove_particle(record.parent);

    for(int j=0; j<child_N; ++j)
    {
      p_Particle sp = split_particles[j];
      my_real coord_shift = my_minus_data(matrx2my_real(sp->coord), box_l);
      my_int pos = get_cell_id(coord_shift, current_level->dcell, current_level->cell_start, current_level->cell_end);

      current_level->table_cell_list[pos.i][pos.j][pos.k]->Add_particle(sp);
    }
  }


  //-----------------------------------------------------
  // 4. inlet particle -> cell list
  for(int i=0; i<int(inlet_records.size()); ++i)
  {
    p_Particle cp = inlet_records[i].particle;

    int level_index = cp->level - Lmin;
    p_Level_info current_level = level_info[level_index];

    my_real coord_shift = my_minus_data(matrx2my_real(cp->coord), box_l);
    my_int pos = get_cell_id(coord_shift, current_level->dcell, current_level->cell_start, current_level->cell_end);

    current_level->table_cell_list[pos.i][pos.j][pos.k]->Add_particle(cp);
  }


  //-----------------------------------------------------
  // 5. outlet particle -> remove from cell list
  for(int i=0; i<int(outlet_records.size()); ++i)
  {
    OutletRecord& record = outlet_records[i];
    p_Level_info current_level = level_info[record.level_index];

    current_level->table_cell_list[record.old_cell.i][record.old_cell.j][record.old_cell.k]->Remove_particle(record.particle);
  }


  //-----------------------------------------------------
  // 6. collect affected particles
  unordered_set<p_Particle> affected_set;
  unordered_set<p_Cell_list> scanned_cells;
  vector<p_Particle> affected_particles;

  auto add_affected = [&](p_Particle particle)
  {
    if(particle != NULL && !outlet_set.count(particle) && affected_set.insert(particle).second) affected_particles.push_back(particle);
  };

  int i_c = DIM_X==1 ? 2 : 0;
  int j_c = DIM_Y==1 ? 2 : 0;
  int k_c = DIM_Z==1 ? 2 : 0;


  //-----------------------------------------------------
  // 6.1 split affected particles
  for(int i=0; i<int(split_records.size()); ++i)
  {
    FastSplitRecord& record = split_records[i];

    for(int j=0; j<int(record.old_neighbors.size()); ++j) add_affected(record.old_neighbors[j]);

    #if DIM == 2
    p_Particle split_particles[4];
    #elif DIM == 3
    p_Particle split_particles[8];
    #endif

    split_particles[0] = record.parent;
    for(int j=0; j<child_N-1; ++j) split_particles[j+1] = record.child[j];

    for(int j=0; j<child_N; ++j)
    {
      p_Particle sp = split_particles[j];
      add_affected(sp);

      int level_index = sp->level - Lmin;
      p_Level_info current_level = level_info[level_index];

      my_real coord_shift = my_minus_data(matrx2my_real(sp->coord), box_l);
      my_int pos = get_cell_id(coord_shift, current_level->dcell, current_level->cell_start, current_level->cell_end);

      for(int ii=AMAX1(pos.i-i_c,current_level->cell_start.i); ii<=AMIN1(pos.i+i_c,current_level->cell_end.i-1); ++ii)
      for(int jj=AMAX1(pos.j-j_c,current_level->cell_start.j); jj<=AMIN1(pos.j+j_c,current_level->cell_end.j-1); ++jj)
      for(int kk=AMAX1(pos.k-k_c,current_level->cell_start.k); kk<=AMIN1(pos.k+k_c,current_level->cell_end.k-1); ++kk)
      {
        p_Cell_list cell = current_level->table_cell_list[ii][jj][kk];
        if(!scanned_cells.insert(cell).second) continue;
        for(int k=0; k<int(cell->particle_list.size()); ++k) add_affected(cell->particle_list[k]);
      }
    }
  }


  //-----------------------------------------------------
  // 6.2 inlet affected particles
  for(int i=0; i<int(inlet_records.size()); ++i)
  {
    p_Particle cp = inlet_records[i].particle;
    add_affected(cp);

    int level_index = cp->level - Lmin;
    p_Level_info current_level = level_info[level_index];

    my_real coord_shift = my_minus_data(matrx2my_real(cp->coord), box_l);
    my_int pos = get_cell_id(coord_shift, current_level->dcell, current_level->cell_start, current_level->cell_end);

    for(int ii=AMAX1(pos.i-i_c,current_level->cell_start.i); ii<=AMIN1(pos.i+i_c,current_level->cell_end.i-1); ++ii)
    for(int jj=AMAX1(pos.j-j_c,current_level->cell_start.j); jj<=AMIN1(pos.j+j_c,current_level->cell_end.j-1); ++jj)
    for(int kk=AMAX1(pos.k-k_c,current_level->cell_start.k); kk<=AMIN1(pos.k+k_c,current_level->cell_end.k-1); ++kk)
    {
      p_Cell_list cell = current_level->table_cell_list[ii][jj][kk];
      if(!scanned_cells.insert(cell).second) continue;
      for(int k=0; k<int(cell->particle_list.size()); ++k) add_affected(cell->particle_list[k]);
    }
  }


  //-----------------------------------------------------
  // 6.3 outlet affected particles
  for(int i=0; i<int(outlet_records.size()); ++i)
  {
    OutletRecord& record = outlet_records[i];

    for(int j=0; j<int(record.old_neighbors.size()); ++j)
    {
      p_Particle neighbor = record.old_neighbors[j];
      if(neighbor != record.particle) add_affected(neighbor);
    }
  }


  //-----------------------------------------------------
  // 7. reset + refresh neighbors
  parallel_for(blocked_range<int>(0, int(affected_particles.size())), [&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i) 
    {
      affected_particles[i]->Reset_neighbor_info();
      affected_particles[i]->Refresh_neighbor_info(this, 3);
    }
    
  }, ap);


  //-----------------------------------------------------
  // 8. remove outlet particles from particle_fluid
  concurrent_vector<p_Particle> particle_tmp;

  for(int i=0; i<int(particle_fluid.size()); ++i)
  {
    p_Particle cp = particle_fluid[i];
    if(outlet_set.count(cp)) continue;

    cp->mark_in =0;
    if(cp->coord(0) - box_l.i < pow(cp->vol,1./DIM) ) cp->mark_in = 1; //2026.7.24
    particle_tmp.push_back(cp);
  }

  particle_fluid.swap(particle_tmp);


  //-----------------------------------------------------
  // 9. outlet particles can now be freed
  for(int i=0; i<int(outlet_records.size()); ++i)
  {
    p_Particle cp = outlet_records[i].particle;
    cp->Clearup();
    particlepool.free(cp);
  }


}


//-----------------------------------------------------
// Merge_particle
//-----------------------------------------------------
void MPM::Merge_particle(communicator &world) 
{
  concurrent_vector  <p_Particle> p_temp1;
  p_temp1.clear();

  Real vol_min = powern(p_sz, DIM)* 2./3. ;
  

  static affinity_partitioner ap;
	parallel_for( blocked_range<int>(0, int(particle_fluid.size())),[&](const blocked_range<int>& r){
		for(int i=r.begin(); i!=r.end(); ++i){
		p_Particle cp = particle_fluid[i];
		cp->copy_p.clear();


		p_Particle cp1  = NULL;

		if( cp->vol < vol_min )
		{
			Real   dss = pow(vol_min, 1./DIM);//2026.7.24
			for (int j = 0; j < int(cp->neighbor.size()); j++)
			{
				Particle *neigh = cp->neighbor[j];
				if (cp != neigh)
				{
					Real  dist = (cp->coord- neigh->coord).head(DIM).norm(); 
					
					if( neigh->vol < vol_min && dist < dss && cp->phase == neigh->phase ) // to test 20251017
					{
            dss = dist;
            cp1 = neigh;
            cp->mrg =1;
					}

				}
			}
		}
		cp->copy_p.push_back(cp1);

	  }
  }, ap); // we can do iteration here

  for(int i=0; i!= int(particle_fluid.size()); ++i)
  {
    p_Particle cp = particle_fluid[i];
    if(cp->mrg == 1 )
    {
      if(cp->copy_p[0] != NULL && cp->copy_p[0]->copy_p[0]!= NULL )
      {
        if(cp->copy_p[0]->copy_p[0] == cp && cp->mrg_mark==0 )
        {
          cp->mrg_mark=1;
          cp->copy_p[0]->mrg_mark=1;
      
          p_Particle partcl= particlepool.malloc();
          partcl->Set_merge_ptcl(cp, this);
          partcl->copy_p.clear();
          p_temp1.push_back(partcl);
        }
      }
    }
  }



  parallel_for(blocked_range<int>(0, int(particle_fluid.size())),[&](const blocked_range<int>& r){
    for(int i=r.begin(); i!=r.end(); ++i){
      p_Particle cp = particle_fluid[i];

      if( cp->mrg_mark == 1){
        particle_fluid[i]->Clearup();
        particlepool.free(particle_fluid[i]);
      }
      else{
        cp->copy_p.clear();
		    p_temp1.push_back(cp);
	    }

    }
  }, ap);


  if(p_temp1.size()>0)
  {
	  particle_fluid.clear();

	  parallel_for( blocked_range<int>(0, p_temp1.size()),
			   [&](const blocked_range<int>& r){
		for(int i=r.begin(); i!=r.end(); ++i){
			p_temp1[i]->mrg_mark =0;
			p_temp1[i]->mrg = 0;			
			particle_fluid.push_back(p_temp1[i]);
      p_Particle cp = p_temp1[i];
      cp->mark_in =0;
      if(cp->coord(0) - box_l.i < pow(cp->vol, 1./DIM) ) cp->mark_in = 1; //2026.7.24
		}
	  }, ap);
  }


}



//-----------------------------------------------------
// Merge_particle_fast
//-----------------------------------------------------
void MPM::Merge_particle_fast(communicator &world)
{
  static affinity_partitioner ap;

  Real vol_min = powern(p_sz, DIM)* 2./3. ; //2026.9.21

  struct MergeRecord
  {
    p_Particle survivor;
    p_Particle removed;
    my_int survivor_old_cell;
    my_int removed_old_cell;
    int level_index;
    vector<p_Particle> old_neighbors;
  };

  vector<MergeRecord> merge_records;
  concurrent_vector  <p_Particle> p_temp1;
  p_temp1.clear();


  parallel_for(blocked_range<int>(0, int(particle_fluid.size())),[&](const blocked_range<int>& r){
    for (int i=r.begin(); i!=r.end(); ++i)
    {
      p_Particle cp = particle_fluid[i];
      cp->copy_p.clear();
      cp->mrg = 0;
      cp->mrg_mark = 0;

      p_Particle cp1  = NULL;
      if( cp->vol < vol_min )
      {
        Real   dss = pow(vol_min, 1./DIM);
        for (int j = 0; j < int(cp->neighbor.size()); j++)
        {
          Particle *neigh = cp->neighbor[j];
          if (cp != neigh)
          {
            Real  dist = (cp->coord- neigh->coord).head(DIM).norm(); 
            
            if( neigh->vol < vol_min && dist < dss && cp->phase == neigh->phase ) // to test 20251017
            {
              dss = dist;
              cp1 = neigh;
              cp->mrg =1;
            }

          }
        }
      }
      cp->copy_p.push_back(cp1);
      p_temp1.push_back(cp) ;//2026.9.21

    }
  }, ap);


  for (int i=0; i<int(particle_fluid.size()); ++i)
  {
    p_Particle cp = particle_fluid[i];

    if(cp->mrg == 1 )
    if(cp->copy_p[0] != NULL && cp->copy_p[0]->copy_p[0]!= NULL && cp->copy_p[0]->copy_p[0] == cp && cp->mrg_mark==0 )
    {
      MergeRecord record;
      record.survivor = cp;
      p_Particle partner = cp->copy_p[0] ;
      record.removed = partner;
      record.level_index = cp->level - Lmin;
      p_Level_info current_level = level_info[record.level_index];
      my_real survivor_shift = my_minus_data(matrx2my_real(cp->coord), box_l);
      my_real removed_shift = my_minus_data(matrx2my_real(partner->coord), box_l);
      record.survivor_old_cell = get_cell_id(survivor_shift, current_level->dcell, current_level->cell_start, current_level->cell_end);
      record.removed_old_cell = get_cell_id(removed_shift, current_level->dcell, current_level->cell_start, current_level->cell_end);
      for (int j=0; j<int(cp->neighbor.size()); ++j)
        record.old_neighbors.push_back(cp->neighbor[j]);
      for (int j=0; j<int(partner->neighbor.size()); ++j)
        record.old_neighbors.push_back(partner->neighbor[j]);

      cp->Set_merge_ptcl(cp, this);
      cp->mrg_mark = 0;
      partner->mrg_mark = 1;
      merge_records.push_back(record);
    }

  }
  

  unordered_set<p_Particle> affected_set;
  unordered_set<p_Cell_list> scanned_cells;
  vector<p_Particle> affected_particles;
  auto add_affected = [&](p_Particle particle)
  {
    if (particle != NULL && particle->mrg_mark != 1 && affected_set.insert(particle).second)
      affected_particles.push_back(particle);
  };

  int     i_c = DIM_X==1 ? 2 : 0;
  int     j_c = DIM_Y==1 ? 2 : 0;
  int     k_c = DIM_Z==1 ? 2 : 0;

  for (int r=0; r<int(merge_records.size()); ++r)
  {
    MergeRecord& record = merge_records[r];
    p_Level_info current_level = level_info[record.level_index];
    for (int j=0; j<int(record.old_neighbors.size()); ++j)
    {
      p_Particle neighbor = record.old_neighbors[j];
      if (neighbor != record.survivor && neighbor != record.removed)
      {
        add_affected(neighbor);
      }
    }
    current_level->table_cell_list[record.survivor_old_cell.i][record.survivor_old_cell.j][record.survivor_old_cell.k]
      ->Remove_particle(record.survivor);
    current_level->table_cell_list[record.removed_old_cell.i][record.removed_old_cell.j][record.removed_old_cell.k]
      ->Remove_particle(record.removed);

    my_real new_shift = my_minus_data(matrx2my_real(record.survivor->coord), box_l);
    my_int new_cell = get_cell_id(new_shift, current_level->dcell,  current_level->cell_start, current_level->cell_end);
    current_level->table_cell_list[new_cell.i][new_cell.j][new_cell.k]->Add_particle(record.survivor);
    add_affected(record.survivor);

    for (int ii=AMAX1(new_cell.i-i_c,current_level->cell_start.i); ii<=AMIN1(new_cell.i+i_c,current_level->cell_end.i-1); ++ii)
    for (int jj=AMAX1(new_cell.j-j_c,current_level->cell_start.j); jj<=AMIN1(new_cell.j+j_c,current_level->cell_end.j-1); ++jj)
    for (int kk=AMAX1(new_cell.k-k_c,current_level->cell_start.k); kk<=AMIN1(new_cell.k+k_c,current_level->cell_end.k-1); ++kk)
    {
      p_Cell_list cell = current_level->table_cell_list[ii][jj][kk];
      if (!scanned_cells.insert(cell).second)
        continue;
      for (int j=0; j<int(cell->particle_list.size()); ++j)
        add_affected(cell->particle_list[j]);
    }
  }

  parallel_for(blocked_range<int>(0, int(affected_particles.size())),
    [&](const blocked_range<int>& r){
      for (int i=r.begin(); i!=r.end(); ++i)
      {
        affected_particles[i]->Reset_neighbor_info();
        affected_particles[i]->Refresh_neighbor_info(this, 3);//2dx
      }

    }, ap);



  particle_fluid.clear();

  parallel_for(blocked_range<int>(0, int(p_temp1.size())),[&](const blocked_range<int>& r){
  for (int i=r.begin(); i!=r.end(); ++i)
  {
    p_Particle particle = p_temp1[i];
    if (particle->mrg_mark == 1)
    {
      particle->Clearup();
      particlepool.free(particle);
      continue;
    }
    particle->mrg = 0;
    particle->mrg_mark = 0;
    particle->copy_p.clear();
    particle->mark_in =0;
    if(particle->coord(0) - box_l.i < pow(particle->vol, 1./DIM) ) particle->mark_in = 1; //2026.9.22
    particle_fluid.push_back( particle );

  }
  }, ap);


}

