#ifndef PARTICLE_H
#define PARTICLE_H
#include <cmath>
#include "glbcls.h"
#include "Mypool.h"

using namespace std;

class Particle_base{
private:

public:

// variables
  int      id;
  my_int   idI;//2026.1.12
  int      color;
  int      level;        //tree level
  int       phase; //

  Real     mass;
  Real     timestep;

  Real     P;            //pressure
  Real     vol;          //volumn
  Real     h;            //smooth length

  Vector<Real, 3>  coord;   //coordinate
  Vector<Real, 3>  coord_X;
  Vector<Real, 3>  v;      //material velocity
  Vector<Real, 3>  a;     //acceleration

  bool used_st=0;

  Real    rho[FldNum];  
  Real    rho0[FldNum];  
  Vector<Real, 3> rhov[FldNum];//
  Real    E[FldNum];//2026.1.27
  Real    mu[FldNum];//
  Real     c[FldNum];
  Matrix<Real,UNum, BaseNum+1> coef0[FldNum];//2025.12.25
  Matrix<Real,UNum, BaseNum+1 > coef[FldNum];//2025.12.25
  Matrix<Real,UNum, BaseNum+1 > dcoef[FldNum];//2025.12.26
  Vector<Real, UNum> Lim_u[FldNum];//
  Vector<Real, 3> rhov_0;//
  Vector<Real, 3> coord_0;//
  Real    rho_0;  
  Real E_0 =0;//2026.1.27
  int ph_state[FldNum] ;  //2026.1.16
  int ph_state_b[FldNum] ; //2026.1.16
  Matrix<Real,UNum, DIM > G_Ui[FldNum];//2026.1.16
  Matrix<Real,DIM, DIM > G_v;//2026.1.23
  Real     alph[FldNum];//2026.1.21
  Real vol_0 = 0;//2026.2.6


  Real ac[FldNum];//2026.6.9

  bool surf=0;
  bool surf1=0;
  int basis_l = BaseNum;
  Vector<Real, 3> tv;//2026.1.1



  int interior_mark = 0;
  concurrent_vector <p_Particle> voro_cell;
  concurrent_vector <int> voro_no;//2024.6.3
  concurrent_vector <Real> voro_edge;//2023.6.1
  Real voro_v=0;//2022.6.28
  Vector<Real, 3> del_r;//2022.7.2
  Vector<Real, 3> del_rp;//2022.7.4
  Vector<Real, 3> cen_shift;//2024.7.8





  int nd_mk = 0;//mpm2025.9.11
  Real mises = 0;//mpm2025.9.11

  concurrent_vector <p_Particle> in_ele;//2025.12.25
  concurrent_vector <std::pair<int, int>> nd_in_ele; // mpm2025.9.11

  int bound_mk=0;//mpm2025.9.20

  Real rhos_m;//
  Real rhos_b;//

  my_matrx Q0; //
  Matrix<Real, 3, 3> Q; //2025.3.5
  my_matrx n0_xi; //
  my_matrx Fns; //
  my_matrx Fvs; //
  my_matrx Fws; //


  Matrix<Real, 3, 3> sigma[G_P]; 
  Matrix<Real, 3, 3> Fn[G_P]; 

  Matrix<Real, 3, 3> N_g; //mpm2025.9.20
  Matrix<Real, 3, 3> M_g; //mpm2025.9.20
  Vector<Real, 3> f_g; //mpm2025.9.20
  Vector<Real, 3> T_g; //mpm2025.9.20
  Matrix<Real, 3, 3> xFFg;

  Real JJ=1;
  Real radi = CR_2D;//mpm2025.10.1
  Real edge= 1;//2025.11.11
  
  Vector<Real, 3> th_nn;//2023.2.14 mpi
  Vector<Real, 3> norm;//2023.2.14 
  Vector<Real, 3> th_nn0;//mpm2025.9.20

  Vector<Real, 3> w_th;//
  Vector<Real, 3> a_th;//


  Vector<Real, 3> dth_nn;//
  Vector<Real, 3> dcoord;//
  Vector<Real, 3> w_th0;//


  bool low_p=0;//mpm2025.10.4
  Real w_p=0;//mpm2025.10.4

  bool cont =0;//2025.12.11
  Real strn=0;//2025.12.11

  concurrent_vector <p_Particle> copy_p;//2026.2.7
  bool spl=0;
  bool mrg=0;
  bool mrg_mark=0;


  concurrent_vector <p_Particle> neighbor; //
  concurrent_vector <p_Particle> neigh_cell; //
  concurrent_vector <p_Particle> neighbor_mm; //
  concurrent_vector <p_Particle> neighbor_ms; //
  concurrent_vector <p_Particle> neighbor_sm; // 
  concurrent_vector <p_Particle> neighbor_ss; // 
  // concurrent_vector <p_Particle> neighbor_con;//

  Particle_base();
  ~ Particle_base(){};

  // functions
  Real    Kernel_function(Real dist, Real h);
  Real    Derivative_kernel_function(Real dist, my_real dr, Real h);
  Real    Derivative_h_kernel_function(Real dist, Real h);
  Real    Kernel_function_lowerD(Real dist, Real h, int dim);
  Real    Derivative_kernel_function_lowerD(Real dist, my_real dr, Real h, int dim);//2023.8.1


  void    Set_timestep(my_real box_l, my_real box_r, MPM *mpm);//2023.2.26
  void    Clearup();
  void    Reset_neighbor_info();
  void    Add_neighbor(Particle *current_particle);
  void    Refresh_neighbor_info(MPM *mpm, int flag);  //mpm2025.9.11
  void    Refresh_ele_info(MPM *mpm, int flag);//2025.12.25

  void    Grid_acceleration(SOLVER *mpm, int s);//2026.1.20
  void    Update_state_onestep( p_Particle cell_i, int md, SOLVER *mpm);//2025.12.25

  void    PST(SOLVER *mpm,  Vector<Real, 3>& pst, Vector<Real, 3>& ph_sharp);;//2026.1.12
  Vector<Real, 3>    Bound_pen(SOLVER *mpm);//2026.1.1
  void                Update_PST(SOLVER *mpm);//2026.1.12

  void  Set_split_ptcl(p_Particle cp, MPM *mpm);
  void  Set_merge_ptcl(p_Particle cp, MPM *mpm);//2026.2.7
  void  Copy_split_state(p_Particle cp);//2026.9.21
};

#endif
