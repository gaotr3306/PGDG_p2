#include "glbfunc.h"
#include "mpm.h"
#include "level_infor.h"
#include "particle.h"

/***************************************************/
/*                                                 */
/*      Functions defined in class "Particle"      */
/*                                                 */
/***************************************************/


//-----------------------------------------------------
// Grid_acceleration 2025.12.25
//-----------------------------------------------------
void Particle_base::Grid_acceleration(SOLVER *mpm, int ss)
{
  Vector<Real, 3> ps; SetZero(ps); 
  ps(0) = mpm->p_size.i/2; ps(1) = mpm->p_size.j/2; 

  SetZero(dcoef[ss]);
  Matrix<Real, UNum, UNum> LIM; SetZero(LIM);
  Matrix<Real, UNum, UNum> LIMj; SetZero(LIMj);
  for(int s=0;s< UNum; s++) LIM(s,s) = Lim_u[ss](s);

  Real voli = mpm->cell_dg[idI.i][idI.j].Vol; 

  Vector<Real, 6> mij; mij<< 1., pow(ps(0),2)/3., pow(ps(1),2)/3., 
  pow(ps(0)*ps(1),2)/9, pow(ps(0),4)* 4./45., pow(ps(1),4)* 4./45.; 
  mij *= voli; 


  //////cell integral
  
  for(int g1=0; g1<int(mpm->Gp.size()); g1++)
  for(int g2=0; g2<int(mpm->Gp.size()); g2++)
  #if DIM ==3
  for(int g3=0; g3<int(mpm->Gp.size()); g3++)
  #endif
  {
    Vector<Real, 3> dp; SetZero(dp);
    dp(0) = mpm->Gp[g1]* ps(0); 
    dp(1) = mpm->Gp[g2]* ps(1); 
    Real gw = mpm->Gw[g1]* mpm->Gw[g2]; 
    #if DIM ==3
    dp(2) = mpm->Gp[g3]* ps(2); 
    gw = mpm->Gw[g1]* mpm->Gw[g2]* mpm->Gw[g3]; 
    #endif


    Vector<Real, UNum> dU = LIM*(coef[ss].rightCols(BaseNum)* (Basis_f(dp, mpm->p_size).tail(BaseNum)));

    Real rho_p = rho[ss] + dU(0);
    Real E_p = E[ss] + dU(3);
    Vector<Real, DIM> rhov_p; SetZero(rhov_p);
    rhov_p.head(DIM) = rhov[ss].head(DIM) + dU.segment<DIM>(1);
    Vector<Real, DIM> v_p = rhov_p.head(DIM)/rho_p;
    Real e_p = E_p - 0.5* rhov_p.head(DIM).squaredNorm()/rho_p;
    Real P_p = mpm->EOS_P(e_p, ss);

    Vector<Real, UNum> Up; Up<<rho_p, rhov_p(0), rhov_p(1), E_p;

    Matrix<Real, DIM, DIM> tau =  ( G_v + G_v.transpose() - 2./3.* G_v.trace()*MatrixXd::Identity(DIM, DIM) )* mu[ss]   ; // 
    Matrix<Real, UNum, DIM> DF; SetZero(DF );
    DF = Up* (v_p ).transpose();
    Matrix<Real, DIM, DIM> strss = P_p* MatrixXd::Identity(DIM, DIM) - tau;

    DF.middleRows<DIM>(1) +=  strss ;
    DF.row(UNum-1) += v_p.transpose()* strss ;

    for(int bi=1; bi<= BaseNum; bi++)
    dcoef[ss].col(bi) +=  DF* Basis_fgrad(dp).row(bi-1).transpose()* gw* voli* 1; 
/////////////////////////////////////////////////////////////////////////////////////

  }
  


  ////////////////////////////face integral//////////////////////////////


  for (int j = 0; j < int(mpm->cell_dg[idI.i][idI.j].cell_id.size()); j++)
  {
    my_int cell_j = mpm->cell_dg[idI.i][idI.j].cell_id[j];
    Vector<Real, 3> Sp_j = mpm->cell_dg[idI.i][idI.j].face_center[j];
    Vector<Real, 3> normal = mpm->cell_dg[idI.i][idI.j].face_normal[j];
    int b_mk = mpm->cell_dg[idI.i][idI.j].bound[j];//2026.7.13
    Real Sj = mpm->cell_dg[idI.i][idI.j].face_S[j];
    Vector<Real, 3> dtan =  mpm->cell_dg[idI.i][idI.j].face_tan[j];
    Vector<Real, 3> dr = Sp_j - coord;


    Matrix<Real, UNum, DIM> DFij[BaseNum+1];
    for(int bi=0; bi<= BaseNum; bi++) SetZero(DFij[bi]); 
    


    for(int g1=0; g1<int(mpm->Gpf.size()); g1++)
    {
      Vector<Real, 3> dp; SetZero(dp);
      dp = mpm->Gpf[g1]* dtan; //2026.1.27
      Real gw = mpm->Gwf[g1]; 
      Vector<Real, 3> dri = dp + Sp_j - coord;
      // Real N_bi = Basis_f(dri)(bi);


      Vector<Real, UNum> dUi = LIM*(coef[ss].rightCols(BaseNum)* Basis_f(dri, mpm->p_size).tail(BaseNum));
      Real rho_i = rho[ss] + dUi(0);
      Real E_i = E[ss] + dUi(3);//2026.1.27
      Vector<Real, DIM> rhov_i; SetZero(rhov_i);
      rhov_i.head(DIM) = rhov[ss].head(DIM) + dUi.segment<DIM>(1);
      Vector<Real, DIM> v_i = rhov_i/rho_i;
      Real e_i = E_i - 0.5* rhov_i.head(DIM).squaredNorm()/rho_i;
      Real P_i = mpm->EOS_P(e_i, ss);
      Real ci = sqrt( mpm->gamma[ss]* (mpm->Pp0[ss] + P_i)/rho_i ) ; //2026.1.27

      Real entr_i = (P_i + mpm->Pp0[ss])/pow(rho_i, mpm->gamma[ss]);//20260521
      Vector<Real, UNum> Ui; Ui<<rho_i, rhov_i(0), rhov_i(1), E_i;

      Real vni = v_i.head(DIM).dot(normal.head(DIM));//2026.7.13

      Matrix<Real, UNum, DIM> DFi, DFj; SetZero(DFi), SetZero(DFj);

      Matrix<Real, DIM, DIM> tau_i =  ( G_v +  G_v.transpose() - 2./3.* G_v.trace()*MatrixXd::Identity(DIM, DIM))* mu[ss] ; // 
      Matrix<Real, DIM, DIM> strs_i = P_i* MatrixXd::Identity(DIM, DIM) - tau_i; //2026.1.27

      DFi = Ui* (v_i).head(DIM).transpose() ;//*alph[ss];
      DFi.middleRows<DIM>(1) +=  strs_i ;//*alph[ss];//
      DFi.row(UNum-1) += v_i.transpose()* strs_i ;//*alph[ss];
      
      if(b_mk == 0)
      {
        p_Particle neigh= mpm->Ele_cen[cell_j.i][cell_j.j];

        for(int s1=0;s1<2; s1++)
        {
          if(neigh->alph[s1] == 0) continue;
          
          for(int s=0;s< UNum; s++) LIMj(s,s) = neigh->Lim_u[s1](s);
          // if(ss != s1  ) LIMj.setZero(); //&& alph[ss]< 0.5
          Vector<Real, 3> drj = dp + Sp_j - neigh->coord;
          Vector<Real, UNum> dUj = LIMj* neigh->coef[s1].rightCols(BaseNum)* (Basis_f(drj, mpm->p_size).tail(BaseNum));

          Real rho_j = neigh->rho[s1]  + dUj(0);
          Real E_j = neigh->E[s1] + dUj(UNum-1);//2026.1.27
          Vector<Real, DIM> rhov_j; SetZero(rhov_j);
          rhov_j.head(DIM) = neigh->rhov[s1].head(DIM) + dUj.segment<DIM>(1);
          Vector<Real, DIM> v_j = rhov_j/rho_j;

          Real rhoe_j = E_j - 0.5* v_j.head(DIM).squaredNorm()* rho_j;
          Real P_j = mpm->EOS_P(rhoe_j, s1);
          Real cj = sqrt( mpm->gamma[s1]* (mpm->Pp0[s1] + P_j)/rho_j ) ; //2026.1.27


          Real vnj = v_j.head(DIM).dot(normal.head(DIM));
          Real vn = 0.5* (vni +  vnj ) ; 

          Vector<Real, DIM> v_ij; v_ij.setZero();
          Real P_ij= P_j; Real rho_ij = rho_j;
          Vector<Real, UNum> Uij; Uij.setZero();

          if(ss != s1)
          {
            // if(ss ==0 && s1 ==1) v_j = v_i;
            // else if(ss ==1 && s1 ==0) P_j = P_i;


            // if(abs(P_j- mpm->Pb0)/P_j < 1e-8)
            // {
            //   P_j = mpm->Pb0 ; 
            //   Vector<Real, DIM> tmpv = (v_i + v_j)* 0.5;
            //   v_i = tmpv;
            //   v_j = tmpv;
            // }

            // if((P_j) < 0)
            // {
            //   P_j = 0 ; 
            // }

            // Real Zi = rho_i* ci; Real Zj = rho_j* cj; 
            // v_ij = (Zi* v_i + Zj* v_j + (P_i - P_j)* normal.head(DIM) )/ (Zi + Zj) ;
            // P_ij = (Zi* P_i + Zj* P_j + Zi*Zj*(v_i - v_j).dot(normal.head(DIM)) )/ (Zi + Zj) ;
            // P_ij = P_j;
            // rho_ij = pow(( P_ij + mpm->Pp0[ss] )/entr_i, 1./mpm->gamma[ss]); //20260521

            // if( rho_j< 0) cout<<"WRONG111"<<endl;
            // Real cij = ( ci*pow(rho_i, 1.5) + cj*pow(rho_j, 1.5) )/(sqrt(rho_i) + sqrt(rho_j)) ;
            // // vn = 0.5* ( (vni +  vnj ) - (P_i - P_j)/cij ) ;
            // P_ij = 0.5* ( P_i + P_j + cij* (vni - vnj) ) ;
            // v_ij = 0.5* ( v_i + v_j + (P_i - P_j)/cij* normal.head(DIM) ) ;
            // if( P_ij + mpm->Pp0[ss] < 0 ) rho_ij = pow(( P_j + mpm->Pp0[ss] )/entr_i, 1./mpm->gamma[ss]); //20260521
            // else rho_ij = pow(( P_ij + mpm->Pp0[ss] )/entr_i, 1./mpm->gamma[ss]); //20260521

            // P_j = P_ij;
            // rho_j = rho_ij;
            // v_j = v_ij;
            // rhoe_j = mpm->EOS_rhoe(P_j, ss);
            // rhov_j = rho_j* v_j;
            // E_j = rhoe_j + 0.5* rho_j* v_j.head(DIM).squaredNorm();//20260521  


            // rho_j = pow(( P_j + mpm->Pp0[ss] )/entr_i, 1./mpm->gamma[ss]); //20260521
            // // rho_j = rho_i ; //20260521
            // rhoe_j = mpm->EOS_rhoe(P_j, ss);
            // rhov_j = rho_j* v_j;
            // E_j = rhoe_j + 0.5* rho_j* v_j.head(DIM).squaredNorm();//20260521  


            Vector<Real, 3> ghost_i = mpm->MGFM(vni, vnj, P_i, P_j, rho_i, rho_j, ss);
            rho_j = ghost_i(0);
            P_j = ghost_i(1);
            Real vnj_p = vnj;
            vnj = ghost_i(2);
            v_j = v_i - vni* normal.head(DIM) + vnj* normal.head(DIM); 
            rhov_j = rho_j* v_j;
            rhoe_j = mpm->EOS_rhoe(P_j, ss);
            E_j = rhoe_j + 0.5* rho_j* v_j.head(DIM).squaredNorm();//20260521  

          }


          cj = sqrt( mpm->gamma[ss]* (mpm->Pp0[ss] + P_j)/rho_j ) ; //2026.1.27

          Vector<Real, UNum> Uj; Uj<< rho_j, rhov_j(0), rhov_j(1), E_j; 
          Matrix<Real, DIM, DIM> g_vj = neigh->G_v;
          Matrix<Real, DIM, DIM> tau_j =  ( g_vj + g_vj.transpose() - g_vj.trace()*MatrixXd::Identity(DIM, DIM) )* mu[ss] ; // 
          Matrix<Real, DIM, DIM> strs_j = P_j* MatrixXd::Identity(DIM, DIM) - tau_j; //2026.1.27

          SetZero(DFj);
          DFj = Uj* (v_j).head(DIM).transpose() ;//*neigh->alph[ss];
          DFj.middleRows<DIM>(1) +=  strs_j ;//*neigh->alph[ss];//
          DFj.row(UNum-1) += v_j.transpose()* strs_j ;//*neigh->alph[ss];


          vnj = v_j.head(DIM).dot(normal.head(DIM));
          Real alph_j = neigh->alph[s1]; 
          Real alph_i = alph[s1]; 

          //////////////////////////////////////////////////////////////
          // Real SL =  AMIN1( (vni - ci), (vnj - cj) )  ;
          // Real SR =  AMAX1( (vni + ci), (vnj + cj) )  ;
          // Real S_s = ( (P_j - P_i) + rho_i * vni * (SL - vni) - rho_j * vnj * (SR - vnj) )/
          //   (  rho_i * (SL - vni) - rho_j * (SR - vnj))   ;

          // Real rho_sL = rho_i * (SL - vni)/(SL - S_s);
          // Real rho_sR = rho_j * (SR - vnj)/(SR - S_s);
          // Vector<Real, DIM> rhov_sL = rho_sL * (v_i + (S_s - vni)*normal.head(DIM)) ;
          // Vector<Real, DIM> rhov_sR = rho_sR * (v_j + (S_s - vnj)*normal.head(DIM)) ;
          // Real E_sL = rho_sL *(E_i/rho_i + (S_s - vni)*( S_s + (P_i )/(rho_i*(SL - vni)) ) );
          // Real E_sR = rho_sR *(E_j/rho_j + (S_s - vnj)*( S_s + (P_j )/(rho_j*(SR - vnj)) ) );

          // Vector<Real, UNum> U_sL; U_sL<<rho_sL, rhov_sL(0), rhov_sL(1), E_sL;
          // Vector<Real, UNum> U_sR; U_sR<<rho_sR, rhov_sR(0), rhov_sR(1), E_sR;

          // Matrix<Real, UNum, DIM> DF_hllc; SetZero(DF_hllc); 
          // // Matrix<Real, 1, DIM> alp_hllc; SetZero(alp_hllc); //2026.2.9
          // if(SL > 0)
          // {
          //   DF_hllc = DFi ;
          //   // alp_hllc = alph_i* v_i.transpose();
          // }
          // else if(SL <= 0 && S_s >=0)
          // {
          //   DF_hllc = DFi + SL* (U_sL - Ui)*  normal.head(DIM).transpose();  
          //   // alp_hllc = alph_i* S_s* normal.head(DIM).transpose();  
          //   // alp_hllc = alph_i* v_i.transpose();
          // }
          // else if(S_s <= 0 && SR >=0)
          // {
          //   DF_hllc = DFj + SR* (U_sR - Uj)*  normal.head(DIM).transpose(); 
          //   // alp_hllc = alph_j * S_s* normal.head(DIM).transpose();     
          //   // alp_hllc = alph_j * v_i.transpose();
          // }
          // else if(SR <= 0 )
          // {
          //   DF_hllc = DFj ;    
          //   // alp_hllc = alph_j* v_i.transpose();
          // }
          // //////////////////////////////////////////////////////////////
          Real SL =  AMIN1( (vni - ci), (vnj - cj) )  ;
          Real SR =  AMAX1( (vni + ci), (vnj + cj) )  ;

          Matrix<Real, UNum, DIM> DF_hllc; SetZero(DF_hllc); 
          if(SL > 0)
          {
            DF_hllc = DFi ;
          }
          else if(SL <= 0 && SR >=0)
          {
            DF_hllc = (SR* DFi - SL* DFj + SL*SR* (Uj - Ui)* normal.head(DIM).transpose() )/ (SR - SL);  
          }
          else if(SR < 0 )
          {
            DF_hllc = DFj ;    
          }
          // /////////////////////////////////////////////////////////

          Real alph_ij = alph_j ;
          Uij = 0.5*(Uj + Ui) ; 

          for(int bi=0; bi<= BaseNum; bi++) 
          {
            Real N_bi = Basis_f(dri, mpm->p_size)(bi);
            DFij[bi] += DF_hllc* N_bi* gw*  alph_ij ; 
          }
          

        }


      }
      else
      {

        Real rho_j = rho_i ;
        Vector<Real, DIM> v_j = v_i;  
        Vector<Real, DIM> rhov_j = rho_j*v_j ;
        Real P_j = P_i;
        Real E_j = E_i;//2026.1.27
        Vector<Real, UNum> Uj; Uj = Ui;

        if(b_mk == 2)
        {
          v_j(1)= -v_i(1) ;
          rhov_j = rho_j*v_j ;
        }
        // else if(b_mk == 2)
        // {
        //   v_j(1)= -v_i(1) ;
        //   rhov_j = rho_j*v_j ;
        // } //2026.7.13


        Uj<<rho_j, rhov_j(0), rhov_j(1), E_j;
        Real cj = sqrt( mpm->gamma[ss]* (mpm->Pp0[ss] + P_j)/rho_j ) ; //2026.1.27
        Real vnj = v_j.head(DIM).dot(normal.head(DIM));
        Matrix<Real, DIM, DIM> strs_j = P_j* MatrixXd::Identity(DIM, DIM) ; 
        SetZero(DFj);
        DFj = Uj* (v_j).head(DIM).transpose() ;//
        DFj.middleRows<DIM>(1) +=  strs_j ;//
        DFj.row(UNum-1) += v_j.transpose()* strs_j ;//


        Real SL =  AMIN1( (vni - ci), (vnj - cj) )  ;
        Real SR =  AMAX1( (vni + ci), (vnj + cj) )  ;
        Matrix<Real, UNum, DIM> DF_hllc; SetZero(DF_hllc); 
        if(SL > 0)
        {
          DF_hllc = DFi ;
        }
        else if(SL <= 0 && SR >=0)
        {
          DF_hllc = (SR* DFi - SL* DFj + SL*SR* (Uj - Ui)* normal.head(DIM).transpose() )/ (SR - SL);  
        }
        else if(SR < 0 )
        {
          DF_hllc = DFj ;    
        }

        for(int bi=0; bi<= BaseNum; bi++) 
        {
          Real N_bi = Basis_f(dri, mpm->p_size)(bi);
          DFij[bi] += DF_hllc* N_bi* gw ; 
        }



      }

    }

    for(int bi=0; bi<= BaseNum; bi++) 
    dcoef[ss].col(bi) += -(DFij[bi] ) *normal.head(DIM)* Sj ;
  }

  
  for(int s=0; s<= BaseNum ; s++)
    dcoef[ss].col(s) *= 1./(mij(s) + 1e-20);

  


}




//---------------------------
// 2nd TVD RK 2025.12.26
//---------------------------
void Particle_base::Update_state_onestep(p_Particle cell_i, int md, SOLVER *mpm )
{
  Real glbl_timestep = mpm->glbl_timestep;
  int ss = phase -1;

  Real divu = cell_i->G_v.trace();//2026.2.6
  if(md==0)
  {
    rho_0 = rho[0];
    coord_0 = coord;
    rhov_0 = rhov[0];
    E_0 = E[0];
    vol_0 = vol;//2026.2.6

    coord.head(DIM) += tv.head(DIM) * glbl_timestep;
    vol =  (1. + divu* glbl_timestep)* vol;//2026.2.6
  }
  else if(md ==1)
  {
    coord.head(DIM) = 0.5*( (coord_0 + coord).head(DIM) + tv.head(DIM) * glbl_timestep ) ;
    vol =  0.5*(vol_0 + vol + divu* glbl_timestep*vol );//2026.2.6
  }

  Vector<Real, 3> dri = coord - cell_i->coord;
  Vector<Real, UNum> Ui ;
  for(int t=0;t<UNum;t++)
  {
    Ui(t) = cell_i->coef[ss](t,0) + 1* cell_i->coef[ss].row(t).rightCols(BaseNum)* Basis_f(dri, mpm->p_size).tail(BaseNum);//
  }
  coef[0] = cell_i->coef[ss];//2026.1.22
  Real rho_i =  Ui(0);
  Vector<Real, 3> rhov_i; SetZero(rhov_i);
  rhov_i.head(DIM) =  Ui.segment<DIM>(1);
  Real E_i =  Ui(3);//2026.1.27
  rho[0] = rho_i;
  rhov[0] = rhov_i;
  E[0] = E_i;//2026.1.27

  v = rhov[0]/rho[0];
  Real ei = E[0]-  0.5*rho[0]* v.head(DIM).squaredNorm();
  P = mpm->EOS_P(ei, phase-1);
  c[0]  =  sqrt(mpm->gamma[phase-1]* ( mpm->Pp0[phase-1] + P)/rho[0] ) ; //2026.1.27

  if(P < -0.05 * mpm->Pp0[phase-1] || rho[0]<0 ) //2026.7.7
  {
    // cell_i->coef[phase -1].rightCols(BaseNum).fill(0);
    // coef[0] = cell_i->coef[phase -1];//2026.1.22
    coef[0].rightCols(BaseNum).fill(0); //2026.6.29
    rho[0] = cell_i->rho[phase-1];
    rhov[0] = cell_i->rhov[phase-1];
    E[0] = cell_i->E[phase-1];//2026.1.27
    c[0]  =  cell_i->c[phase-1];//2026.7.1

  }


  mass = vol* rho[0]; //2026.2.6
  h = sqrt(vol)* CELL_RATIO/CUT_OFF;



}


//---------------------------
// 2nd TVD RK 2025.12.26
//---------------------------
void Particle_base::Update_PST(SOLVER *mpm)
{
  Real glbl_timestep = mpm->glbl_timestep;

  Vector<Real, 3>  dr; SetZero(dr);

  
  SetZero(tv);
  Vector<Real, 3> pst, ph_sharp; SetZero(pst); SetZero(ph_sharp);
  PST(mpm, pst, ph_sharp);
  Vector<Real, 3> dv; SetZero(dv);
  Real px = mpm->p_sz;
  // SetZero(ph_sharp);

  Vector<Real, 3> dx = pst* glbl_timestep* glbl_timestep/2.;
  Real dis = dx.head(DIM).norm();
  Real lim = AMIN1(px/10, dis)/(dis+ 1e-20) ; 
  dx *=  lim;

  Vector<Real, 3> dx1 = ph_sharp* glbl_timestep* glbl_timestep/2.;
  Real dis1 = dx1.head(DIM).norm();
  Real lim1 = AMIN1(px/10, dis1)/(dis1+ 1e-20) ; 
  dx1 *=  lim1;

  Vector<Real, 3> dx2 = Bound_pen(mpm) * glbl_timestep* glbl_timestep/2.;
  Real dis2 = dx2.head(DIM).norm();
  Real lim2 = AMIN1(px/2, dis2)/(dis2+ 1e-20) ; //2026.7.8
  dx2 *=  lim2;

  dv =  (dx + dx1 + dx2)/glbl_timestep ;

  
  tv += dv;

}

//-----------------------------------------------------
// PST
//-----------------------------------------------------
void Particle_base::PST(SOLVER *mpm,  Vector<Real, 3>& pst, Vector<Real, 3>& ph_sharp)
{
  Real ca = mpm->glbl_max_v* 5. ;
  Real P_b = rho[0]* ca* ca ; // 
  SetZero(pst);
  SetZero(ph_sharp);
  Real v_ww=0;

  Real dx = h/CELL_RATIO* CUT_OFF;//2026.7.3

  for (int i = 0; i < int(neighbor.size()); i++){
    Particle *neigh = neighbor[i];
    Vector<Real, 3>   dr = coord - neigh->coord;
    my_real dr1 = matrx2my_real(dr);
    Real    dist = dr.head(DIM).norm();


    if ( dist != 0) //
    {
      Real    vol_j = neigh->vol; //
      Real    aw   = Derivative_kernel_function(dist, dr1, h);
      my_real awax = my_multiply_const (dr1, aw);
      Real  Wij = Kernel_function(dist, h); 
      Real  W_dp = Kernel_function( dx, h);//2026.7.8
      Real R_Wij_W_dp = 0.2* pow( Wij/W_dp, 4);
      Real    tmp1  = -1.*P_b*(1./rho[0] ) *(1 + R_Wij_W_dp);

      Real tmp2 = 1./(rho[0] + neigh->rho[0]);

      if( neigh->phase== 1 || neigh->phase== 2) 
      {
        pst +=  my_real2matrx(awax)* tmp1* vol_j;

        if(neigh->phase!= phase)
          // ph_sharp += -0.1 *tmp2 *(abs(P) + abs(neigh->P))* my_real2matrx(awax)* vol_j;
          ph_sharp += -0.1 *tmp2 *(c[0] * c[0]* rho[0] + neigh->c[0]* neigh->c[0]* neigh->rho[0])* my_real2matrx(awax)* vol_j;

      }
      else if( neigh->phase== 10 )
      { 
        Vector<Real, 3>  normal = neigh->norm;
        pst +=  normal* Wij* tmp1* neigh->vol;
      }

    }
    if(neigh->surf ==1) surf1=1;
  }

  tv = v;

}


//-----------------------------------------------------
// Bound_pen
//-----------------------------------------------------
Vector<Real, 3> Particle_base::Bound_pen(SOLVER *mpm)
{
  Vector<Real, 3> pst; SetZero(pst);

  p_Particle cell_i = in_ele[0];
  my_int idc = cell_i->idI;
  for (int j = 0; j < int(mpm->cell_dg[idc.i][idc.j].cell_id.size()); j++)
  {
    Vector<Real, 3> Sp_j = mpm->cell_dg[idc.i][idc.j].face_center[j];
    Vector<Real, 3> normal = mpm->cell_dg[idc.i][idc.j].face_normal[j];
    bool b_mk = mpm->cell_dg[idc.i][idc.j].bound[j];

    if(b_mk == 1)
    {
      Vector<Real, 3> dr =  coord - Sp_j ;
      Real ds = dr.head(DIM).dot(normal.head(DIM));
      if(ds < 0) normal = -1* normal;
  
      Real dis = abs(ds);
      if(dis< 0.4* mpm->p_sz )
      {
        pst += 0.02/dis * pow(c[0],2)* normal;
      }
    }
  }
  return pst;

}

