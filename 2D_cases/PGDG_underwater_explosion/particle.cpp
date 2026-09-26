#include "glbfunc.h"
#include "mpm.h"
// #include "mpm_incprs.h"
#include "level_infor.h"
#include "particle.h"

/***************************************************/
/*                                                 */
/*      Functions defined in class "Particle"      */
/*                                                 */
/***************************************************/

Particle_base::Particle_base()
{
  P           = 0.;
  mass        = 0.;
  vol         = 0.;
  h           = 0.;
  id          = 0;
  SetZero (coord);
  SetZero (v);
  SetZero (coord_X);
  SetZero (a);     
  

  level =0;
  phase =0;

  SetZero(norm);//2023.2.14 

  SetZero(rhov_0);//
  SetZero(coord_0);//
  rho_0 =0;
  SetZero(tv);//


  for(int i=0;i<2;i++)
  {
    rho[i] = 0;  
    rho0[i] = 0;  
    rhov[i].fill(0);//
    mu[i] = 0;//
    c[i] = 0;//
    coef0[i].fill(0);//
    coef[i].fill(0);//
    dcoef[i].fill(0);//
    Lim_u[i].fill(0);//
    ph_state[i] = 1;  //2026.1.16
    ph_state_b[i] = 0; //2026.1.16
    alph[i] =0;
    E[i] = 0;//2026.1.27
    aph_dt[i] = 0;
    alph0[i] =0;//2026.2.7
  }
  SetZero(G_v);
}


//-----------------------------------------------------------------
// clear particle info
//-----------------------------------------------------------------
void Particle_base::Clearup(){
  neighbor.clear();


  P           = 0.;
  mass        = 0.;
  vol         = 0.;
  h           = 0.;
  id          = 0;
  SetZero (coord);
  SetZero (v);
  SetZero (coord_X);
  SetZero (a);            //acceleration

  level =0;
  phase =0;

  SetZero(norm);//2023.2.14 

  SetZero(rhov_0);//
  SetZero(coord_0);//
  rho_0 =0;
  SetZero(tv);//


  for(int i=0;i<2;i++)
  {
    rho[i] = 0;  
    rho0[i] = 0;  
    rhov[i].fill(0);//
    mu[i] = 0;//
    c[i] = 0;//
    coef0[i].fill(0);//
    coef[i].fill(0);//
    dcoef[i].fill(0);//
    Lim_u[i].fill(0);//
    ph_state[i] = 1;  //2026.1.16
    ph_state_b[i] = 0; //2026.1.16
    alph[i] =0;
    E[i] = 0;//2026.1.27
    aph_dt[i] = 0;
    alph0[i] =0;//2026.2.7
  }
  SetZero(G_v);
  
}
//-----------------------------------------------------
// get the value of kernel function
//-----------------------------------------------------
Real Particle_base::Kernel_function(Real dist, Real h)
{
#ifdef _CUBIC_SPLINE_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

  if (DIM == 1)
    C = 0.666666667 / powern(h, 1);
  else if (DIM == 2)
    C = 0.454728408 / powern(h, 2);
  else if (DIM == 3)
    C = 0.318309886 / powern(h, 3);

  if (s < 1)
    val = C * ( 1 - 1.5*s*s + 0.75*s*s*s );
  else if (s < 2)
    val = C/4.0 * (2-s) * (2-s) * (2-s);
  else
    val = 0.0;

  return val;
#endif
#ifdef _WINDLAND_C4_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

  if (DIM == 1)
    cout<<"WendlandQuintic: Dim 1 not supported\n";
  else if (DIM == 2)
    C = 1.0 / PI / powern(h, 2) * 7. / 4.;
  else if (DIM == 3)
    C = 1.0 / PI / powern(h, 3) * 21. / 16.;

  if (s < 2)
    val = C * powern((1. - 0.5*s), 4)*(1. + 2.*s);
  else
    val = 0.0;

  return val;
#endif
#ifdef _GAUSSIAN_
  Real fac = 0.0;
  Real   s = dist / h;
  Real val = 0.0;

  fac = 0.5*2/sqrt(PI);
  if (DIM > 1)
    fac *= 0.5*2/sqrt(PI);
  else if (DIM > 2)
    fac *= 0.5*2/sqrt(PI);

  Real   C = 0.;

  if (DIM == 1)
    C = fac/h;
  else if (DIM == 2)
    C = fac / powern(h, 2);
  else if (DIM == 3)
    C = fac / powern(h, 3);

  if (s < 3)
    val = exp(-1.*s*s)*C;
  else
    val = 0.0;

  return val;
#endif
#ifdef _QUINTIC_SPLINE_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

  if (DIM == 1)
    C = 1.0 / PI * 7./478. / powern(h, 1);
  else if (DIM == 2)
    C = 1.0 / PI * 7./478. / powern(h, 2);
  else if (DIM == 3)
    C = 1.0 / PI * 7./478. / powern(h, 3);

  if (s < 1)
    val = C * (powern ((3. - s), 5) - 6.*powern ((2. - s), 5) + 15.*powern ((1. - s), 5));
  else if (s < 2)
    val = C * (powern ((3. - s), 5) - 6.*powern ((2. - s), 5));
  else if (s < 3)
    val = C * (powern ((3. - s), 5));
  else
    val = 0.0;

  return val;
#endif
}
//-----------------------------------------------------
// mpm2025.10.1
//-----------------------------------------------------
Real Particle_base::Kernel_function_lowerD(Real dist, Real h, int dim)
{
  Real fac = 0.0;
  Real   s = dist / h;
  Real val = 0.0;

  fac = 0.5*2/sqrt(PI);
  if (dim == 2)
    fac *= 0.5*2/sqrt(PI);

  Real   C = 0.;

  if (dim == 1)
    C = fac/h;
  else if (dim == 2)
    C = fac / powern(h, 2);

  
  if (s < 3)
    val = exp(-1.*s*s)*C;
  else
    val = 0.0;

  return val;
}

//-----------------------------------------------------
// mpm2025.10.1
//-----------------------------------------------------
Real Particle_base::Derivative_kernel_function_lowerD
(Real dist, my_real dr, Real h, int dim)
{
  Real val = 0.0;
  if(dim ==1)
  {
    Real fac = 0.0;
    Real   s = dist / h;

    fac = 0.5*2/sqrt(PI);
    Real   C = fac/h;
  
    if (s < 3)
      val = -2.*s*exp(-1.*s*s)*C/h/dist;
    else
      val = 0.0;
  }
  else if(dim ==2)
  {
    Real    C   = 0.0;
    Real    s   = dist / h;
    C = 1.0 / PI / powern(h, 2) * 7. / 4.;
  
    if (s < 2)
      val = -5.*C*powern((1. - 0.5*s), 3)*s/h/dist;
    else
      val = 0.0;
  }



  return val;


}

//-----------------------------------------------------
// get the derivative value of kernel function
//-----------------------------------------------------
Real Particle_base::Derivative_kernel_function
(Real dist, my_real dr, Real h)
{
#ifdef _CUBIC_SPLINE_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

  if (DIM == 1)
    C = 0.666666667 / powern(h, 2);
  else if (DIM == 2)
    C = 0.454728408 / powern(h, 3);
  else if (DIM == 3)
    C = 0.318309886 / powern(h, 4);

  if (s < 1)
    val = 3*C*(-s + 3*(s*s)/4 );
  else if (s < 2)
    val = -3/4.0*C*(2-s)*(2-s);
  else
    val = 0.0;

  val /= (dist + 1.e-20);
  return val;
#endif
#ifdef _WINDLAND_C4_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

  if (DIM == 1)
    cout<<"WendlandQuintic: Dim 1 not supported\n";
  else if (DIM == 2)
    C = 1.0 / PI / powern(h, 2) * 7. / 4.;
  else if (DIM == 3)
    C = 1.0 / PI / powern(h, 3) * 21. / 16.;

  if (s < 2)
    val = -5.*C*powern((1. - 0.5*s), 3)*s/h/dist;
  else
    val = 0.0;

  return val;
#endif
#ifdef _GAUSSIAN_
  Real fac = 0.0;
  Real   s = dist / h;
  Real val = 0.0;

  fac = 0.5*2/sqrt(PI);
  if (DIM > 1)
    fac *= 0.5*2/sqrt(PI);
  else if (DIM > 2)
    fac *= 0.5*2/sqrt(PI);

  Real   C = 0.;

  if (DIM == 1)
    C = fac/h;
  else if (DIM == 2)
    C = fac / powern(h, 2);
  else if (DIM == 3)
    C = fac / powern(h, 3);

  if (s < 3)
    val = -2.*s*exp(-1.*s*s)*C/h/dist;
  else
    val = 0.0;

  return val;
#endif
#ifdef _QUINTIC_SPLINE_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

  if (DIM == 1)
    C = 1.0 / PI * 7./478. / powern(h, 1);
  else if (DIM == 2)
    C = 1.0 / PI * 7./478. / powern(h, 2);
  else if (DIM == 3)
    C = 1.0 / PI * 7./478. / powern(h, 3);

  if (s < 1){
    val = C * (-5. * powern ((3. - s), 4) + 30.*powern ((2. - s), 4) - 75.*powern ((1. - s), 4));
    val = val/h/dist;
  }
  else if (s < 2){
    val = C * (-5. * powern ((3. - s), 4) + 30.*powern ((2. - s), 4));
    val = val/h/dist;
  }
  else if (s < 3)
    val = C * (-5. * powern ((3. - s), 4) /s/dist);
  else
    val = 0.0;

  return val;
#endif
}




//-----------------------------------------------------
// get the derivative value of kernel function regarding
// to smooth length using a quintic spline function
//-----------------------------------------------------
Real Particle_base::Derivative_h_kernel_function
(Real dist, Real h)
{
#ifdef _CUBIC_SPLINE_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

  if (DIM == 1)
    C = 0.666666667 / powern(h, 2);
  else if(DIM == 2)
    C = 0.454728408 / powern(h, 3);
  else if (DIM == 3)
    C = 0.318309886 / powern(h, 4);

  if (s < 1)
    val = C * ( -1.*DIM + 1.5*(DIM+2.)*s*s - 0.75*(DIM+3.)*s*s*s );
  else if (s < 2)
    val = C/4.0 * ( -1.*DIM*powern((2-s),3) +3.*s*powern((2.-s),2));
  else
    val = 0.0;

  return val;
#endif
#ifdef _WINDLAND_C4_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

  if (DIM == 1)
    cout<<"WendlandQuintic: Dim 1 not supported\n";
  else if (DIM == 2)
    C = 1.0 / PI / powern(h, 2) * 7. / 4.;
  else if (DIM == 3)
    C = 1.0 / PI / powern(h, 3) * 21. / 16.;

  Real w = 0.;
  Real dw = 0.;
  Real tmp = 1. - 0.5*s;
  if (s < 2){
    w   = tmp * tmp * tmp * tmp * (2.0*s + 1.0);
    dw  = -5.0 * s * tmp * tmp * tmp;
    val = -1. * C/h * ( dw*s + w*DIM );
  }
  else
    val = 0.0;

  return val;

#endif
#ifdef _GAUSSIAN_
  Real fac = 0.0;
  Real   s = dist / h;
  Real val = 0.0;

  fac = 0.5*2/sqrt(PI);
  if (DIM > 1)
    fac *= 0.5*2/sqrt(PI);
  else if (DIM > 2)
    fac *= 0.5*2/sqrt(PI);

  Real   C = 0.;

  if (DIM == 1)
    C = fac/h;
  else if (DIM == 2)
    C = fac / powern(h, 2);
  else if (DIM == 3)
    C = fac / powern(h, 3);

  Real  w = 0.;
  Real dw = 0.;
  if (s < 3){
    w   = exp(-1.*s*s);
    dw  = -2.*s*w;
    val = -1.*C/h*(dw*s + w*DIM);
  }
  else
    val = 0.0;

  return val;
#endif
#ifdef _QUINTIC_SPLINE_
  Real    C   = 0.0;
  Real    s   = dist / h;
  Real    val = 0.0;

  if (DIM == 1)
    C = 1.0 / PI * 7./478. / powern(h, 1);
  else if (DIM == 2)
    C = 1.0 / PI * 7./478. / powern(h, 2);
  else if (DIM == 3)
    C = 1.0 / PI * 7./478. / powern(h, 3);

  Real  w = 0.;
  Real dw = 0.;
  if (s < 1){
     w  = powern ((3. - s), 5);
     w -= 6. * powern ((2. - s), 5);
     w += 15. * powern ((1. - s), 5);

    dw  = -5. * powern ((3. - s), 4);
    dw += 30. * powern ((2. - s), 4);
    dw -= 75. * powern ((1. - s), 4);
  }
  else if (s < 2){
     w  = powern ((3. - s), 5);
     w -= 6. * powern ((2. - s), 5);

    dw  = -5. * powern ((3. - s), 4);
    dw += 30. * powern ((2. - s), 4);
  }
  else if (s < 3){
     w  = powern ((3. - s), 5);
    dw  = -5. * powern ((3. - s), 4);
  }
  else{
     w  = 0.;
    dw  = 0.;
  }
  val = -1. * C / h * (dw*s + w*DIM);

  return val;
#endif
}



//-----------------------------------------------------
// reset particle neighbor information
//-----------------------------------------------------
void Particle_base::Reset_neighbor_info()
{
  neighbor.clear();
}


//-----------------------------------------------------
// add neighbor for particle
//-----------------------------------------------------
void Particle_base::Add_neighbor(Particle *current_particle)
{
  neighbor.push_back(current_particle);
}

//-----------------------------------------------------
// Refresh neighbor infor
//-----------------------------------------------------
void Particle_base::Refresh_neighbor_info(MPM *mpm, int flag)
{

  p_Level_info    current_level = mpm->level_info[level - mpm->Lmin];
  my_real coord_shift = my_minus_data ( matrx2my_real(coord), mpm->box_l);
  my_int  pos = get_cell_id (coord_shift, current_level->dcell, current_level->cell_start, current_level->cell_end);
  current_level->Refresh_neighbor_info(mpm, this, pos, flag);	//2026.7.8

}

//-----------------------------------------------------
// Refresh neighbor infor //mpm2025.9.11
//-----------------------------------------------------
void Particle_base::Refresh_ele_info(MPM *mpm, int flag)
{
  
  p_Level_info    current_level = mpm->level_info[level -mpm->Lmin];
  my_real coord_shift = my_minus_data (matrx2my_real(coord), mpm->box_l);
  my_int  pos = get_cell_id (coord_shift, current_level->dcell, current_level->cell_start, current_level->cell_end);
  current_level->Refresh_ele_info(mpm, this, pos);	

}


//-----------------------------------------------------
// set timestep size
//-----------------------------------------------------
void Particle_base::Set_timestep(my_real box_l,my_real box_r, MPM *mpm)
{

  Real U = (v).head(DIM).norm();
  Real cfl= 0.3 ; //2026.7.8
  Real ds = h /CELL_RATIO*CUT_OFF;
  Real c12 = AMAX1(c[0], c[1]);
  Real a12 = AMAX1(ac[0], ac[1]);


  timestep = cfl * ds / (c12 + U );//2026.2.9
  Real tmp = cfl * sqrt(ds / a12);//

  timestep = AMIN1(timestep, tmp);//2026.2.1
}


//-----------------------------------------------------
// Set_split_ptcl
//-----------------------------------------------------
void Particle_base::Set_split_ptcl(p_Particle cp, MPM *mpm)
{
	Vector<Real, 3> dr = coord - cp->coord;
  vol = cp->vol/4.0;
  phase= cp->phase;
  mass = cp->mass/4.0;
  h = cp->h/2.0;
  level = cp->level; //2026.9.21

  coef[0] = cp->coef[0];//2026.1.22
  coef0[0] = coef[0];//2026.1.22
  rho[0] = cp->rho[0] + coef[0].row(0).rightCols(BaseNum)* Basis_f(dr, mpm->p_size).tail(BaseNum);
  rhov[0].head(DIM) = cp->rhov[0].head(DIM) + coef[0].middleRows<DIM>(1).rightCols(BaseNum)* Basis_f(dr, mpm->p_size).tail(BaseNum);//2026.7.24
  E[0] = cp->E[0]+ coef[0].row(UNum-1).rightCols(BaseNum)* Basis_f(dr, mpm->p_size).tail(BaseNum);


  a = cp->a;
  tv = v = rhov[0]/rho[0];
  Real ei = E[0]-  0.5*rho[0]* v.head(DIM).squaredNorm();
  P = mpm->EOS_P(ei, phase-1);
  c[0]  =  sqrt(mpm->gamma[phase-1]* (mpm->Pp0[phase-1] + P)/rho[0] ) ; //2026.1.27
  ee = ei;
  
  color = 0;
  mrg_mark = 0;
  mrg=0;
  spl=1;

}


//-----------------------------------------------------
// Copy_split_state //2026.9.21
//-----------------------------------------------------
void Particle_base::Copy_split_state(p_Particle cp)
{
  vol = cp->vol;
  phase= cp->phase;
  mass = cp->mass;
  h = cp->h;
  level = cp->level; //2026.9.21
  coord = cp->coord; //2026.9.21

  coef[0] = cp->coef[0];//2026.1.22
  coef0[0] = cp->coef0[0];//2026.1.22
  rho[0] = cp->rho[0] ;
  rhov[0].head(DIM) = cp->rhov[0].head(DIM) ;
  E[0] = cp->E[0] ;
  c[0] = cp->c[0] ;

  a = cp->a;
  tv = cp->tv;
  P = cp->P;
  
  color = 0;
  mrg_mark = 0;
  mrg=0;
  spl=1;
}

//-----------------------------------------------------
// Set_merge_ptcl //2026.9.21
//-----------------------------------------------------
void Particle_base::Set_merge_ptcl(p_Particle cp, MPM *mpm)
{
  p_Particle cp1 = cp->copy_p[0];

  int phase_old = cp->phase;
  int level_old = cp->level;

  Real mass0 = cp->mass;
  Real mass1 = cp1->mass;
  Real vol0 = cp->vol;
  Real vol1 = cp1->vol;

  Vector<Real,3> coord0 = cp->coord;
  Vector<Real,3> coord1 = cp1->coord;
  Matrix<Real,UNum, BaseNum+1 > coef_0 = cp->coef[0];
  Matrix<Real,UNum, BaseNum+1 > coef_1 = cp1->coef[0];
  Real rho_0 = cp->rho[0];
  Real rho_1 = cp1->rho[0];
  auto rhov_0 = cp->rhov[0];
  auto rhov_1 = cp1->rhov[0];
  Real E_0 = cp->E[0];
  Real E_1 = cp1->E[0];
  auto a_old = cp->a;


  phase = phase_old;
  level = level_old;
  mass = mass0 + mass1;
  vol = vol0 + vol1;
  h = pow(vol, 1./DIM)*CELL_RATIO/CUT_OFF;
  coord = (coord0*vol0 + coord1*vol1)/vol;
  coef[0] = (coef_0*vol0 + coef_1*vol1)/vol;
  coef0[0] = coef[0];
  rho[0] = (rho_0*vol0 + rho_1*vol1)/vol;
  rhov[0].head(DIM) = (rhov_0*vol0 + rhov_1*vol1).head(DIM)/vol;
  E[0] = (E_0*vol0 + E_1*vol1)/vol;

  a = a_old;
  tv = v = rhov[0]/rho[0];
  Real ei = E[0] - 0.5*rho[0]*v.head(DIM).squaredNorm();
  P = mpm->EOS_P(ei, phase-1);
  c[0] = sqrt(mpm->gamma[phase-1]*(mpm->Pp0[phase-1] + P)/rho[0]);

  color = 0;
  mrg_mark = 0;
  mrg = 0;
  spl = 0;
}
