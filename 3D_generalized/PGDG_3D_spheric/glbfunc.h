#ifndef GLBFUNC_H
#define GLBFUNC_H

#include "glbparam.h"
#include "glbcls.h"
#include<cstdlib>
#include <cmath>
#include <boost/serialization/string.hpp>
#include <boost/serialization/utility.hpp>


using namespace std;
using namespace tbb;
using namespace boost;



template<class T>  void SetZero(T &A){
  int row= A.rows();
  int col= A.cols();
  for(int j=0;j<row;j++)
  for(int k=0;k<col;k++)
  A(j,k)=0;

}; //2024.6.10

//  the nth power
template<class T> T powern(T a, int n)
{
  T res = 1;
  if (n >= 0)
    for (int i=0; i<n; i++)
      res *= a;
  else if (n < 0){
    for (int i=0; i<(-1*n); i++)
      res *= a;
    res = 1/res;
  }
  return res;
};
//  Get the maximum
template<class T> T AMAX1(T a, T b)
{
  return (a >= b ? a : b);
};
//  Get the minimum
template<class T> T AMIN1(T a, T b)
{
  return (a <= b ? a : b);
};
//  multiply constant value to my_int or my_real
template<class T1, class T2> T1 my_multiply_const(T1 a, T2 b)
{
  T1 c;
  c.i = DIM_X==1 ? a.i*b : a.i;
  c.j = DIM_Y==1 ? a.j*b : a.j;
  c.k = DIM_Z==1 ? a.k*b : a.k;
  return c;
};
//  multiply self-defined data type to my_int or my_real
template<class T1, class T2> T1 my_multiply_data(T1 a, T2 b)
{
  T1 c;
  c.i = DIM_X==1 ? a.i*b.i : a.i;
  c.j = DIM_Y==1 ? a.j*b.j : a.j;
  c.k = DIM_Z==1 ? a.k*b.k : a.k;
  return c;
};
//  devide self-defined data type
template<class T1, class T2> T1 my_devide_data(T1 a, T2 b)
{
  T1 c;
  c.i = DIM_X==1 ? a.i/b.i : a.i;
  c.j = DIM_Y==1 ? a.j/b.j : a.j;
  c.k = DIM_Z==1 ? a.k/b.k : a.k;
  return c;
};
//  add constant value to my_int or my_real
template<class T1, class T2> T1 my_add_const(T1 a, T2 b)
{
  T1 c;
  c.i = DIM_X==1 ? a.i+b : a.i;
  c.j = DIM_Y==1 ? a.j+b : a.j;
  c.k = DIM_Z==1 ? a.k+b : a.k;
  return c;
};
//  add self-defined data type to my_int or my_real
template<class T1, class T2> T1 my_add_data(T1 a, T2 b)
{
  T1 c;
  c.i = DIM_X==1 ? a.i+b.i : a.i;
  c.j = DIM_Y==1 ? a.j+b.j : a.j;
  c.k = DIM_Z==1 ? a.k+b.k : a.k;
  return c;
};
//  minus self-defined data type
template<class T1, class T2> T1 my_minus_data(T1 a, T2 b)
{
  T1 c;
  c.i = DIM_X==1 ? a.i-b.i : a.i;
  c.j = DIM_Y==1 ? a.j-b.j : a.j;
  c.k = DIM_Z==1 ? a.k-b.k : a.k;
  return c;
};
//  set the value of my_int or my_real according to the same data type
template<class T> void my_set_data(T &a, T b)
{
  a.i = DIM_X==1 ? b.i : a.i;
  a.j = DIM_Y==1 ? b.j : a.j;
  a.k = DIM_Z==1 ? b.k : a.k;
};
//  set the value of my_int or my_real according to another constant
template<class T1, class T2> void my_set_const(T1 &a, T2 b)
{
  a.i = DIM_X==1 ? b : 0;
  a.j = DIM_Y==1 ? b : 0;
  a.k = DIM_Z==1 ? b : 0;//2025.3.15
};
//  multiply all the members
template<class T1, class T2> void my_self_multiply(T1 a, T2 &b)
{
  b = 1;
  b *= (DIM_X==1 ? a.i : 1);
  b *= (DIM_Y==1 ? a.j : 1);
  b *= (DIM_Z==1 ? a.k : 1);
};
//  add all the members
template<class T1, class T2> void my_self_add(T1 a, T2 &b)
{
  b = 0;
  b += (DIM_X==1 ? a.i : 0);
  b += (DIM_Y==1 ? a.j : 0);
  b += (DIM_Z==1 ? a.k : 0);
};
//  output
template<class T> void my_cout(T a, char *title)
{
  cout<<"   "<<title<<"  "<<a.i<<"  "<<a.j<<"  "<<a.k<<"  "<<"\n";
};
//allocate & delete 3d martrix
template<class T> void allocate_3d_matrix(T*** &matrix, my_int res)
{
  matrix = new T**[res.i];
  for(int i=0; i<res.i; i++){
    matrix[i] = new T*[res.j];
    for(int j=0; j<res.j; j++){
      matrix[i][j] = new T[res.k];
    }
  }
};
template <typename T> //switched order for deduction
boost::multi_array<T, 3> allocate_3d_multi_array(int xmin, int xmax, int ymin, int ymax, int zmin, int zmax, T value)
{
  boost::multi_array<T, 3> arr(boost::extents[xmax-xmin][ymax-ymin][zmax-zmin]);
  std::fill(arr.data(), arr.data() + arr.num_elements(), value);
  return arr;
};
//d. Sign of the first value is determined by the secomax_species_number value's sign


inline my_int get_cell_id(my_real coord, my_real dcell, my_int cell_start, my_int cell_end){
  my_int pos;

  pos.i = DIM_X==1 ? int(floor(coord.i/dcell.i)) : 0;
  pos.j = DIM_Y==1 ? int(floor(coord.j/dcell.j)) : 0;
  pos.k = DIM_Z==1 ? int(floor(coord.k/dcell.k)) : 0;

  pos.i = AMAX1(cell_start.i,AMIN1(pos.i,cell_end.i-1));
  pos.j = AMAX1(cell_start.j,AMIN1(pos.j,cell_end.j-1));
  pos.k = AMAX1(cell_start.k,AMIN1(pos.k,cell_end.k-1));

  return pos;
};

inline Real get_distance(my_real dr){
  Real dist = 0.0;
  dist += DIM_X==1 ? dr.i*dr.i : 0.0;
  dist += DIM_Y==1 ? dr.j*dr.j : 0.0;
  dist += DIM_Z==1 ? dr.k*dr.k : 0.0;

  dist = sqrt(dist);
  return dist;
};



inline my_real matrx2my_real(Vector<Real, 3> A){
  my_real val;

  val.i = A(0);
  val.j = A(1);
  val.k = A(2);

  return val;
}; //2023.2.16

inline Vector<Real, 3> my_real2matrx(my_real A){
  Vector<Real, 3> val;

  val(0) = A.i;
  val(1) = A.j ;
  val(2) = A.k;

  return val;
}; //2023.2.17





inline Vector<Real, BaseNum+1> Basis_f(Vector<Real,3> rr, my_real ps) //2026.6.20
{
  Vector<Real,BaseNum+1> V;
  Real x=rr(0);
  Real y=rr(1);
  Real z=rr(2);

#if DIM==2
  if(BaseNum==2)
  V<< 1., x, y;
  else if(BaseNum==3)
  V<< 1., x, y, x*y;
  else if(BaseNum==5)
  V<< 1., x, y, x*y, x*x, y*y;
#elif DIM ==3
  if(BaseNum==3)
  V<< 1., x, y, z; //2026.2.18
#endif


  return V;
} 

inline Matrix<Real, BaseNum, DIM> Basis_fgrad(Vector<Real,3> rr)
{
  Matrix<Real, BaseNum, DIM> V;
  Real x=rr(0);
  Real y=rr(1);
  Real z=rr(2);


#if DIM==2
  if(BaseNum==2)
  V<< 1, 0,
    0, 1;
  else if(BaseNum==3)
  V<< 1, 0,
    0, 1,
    y, x;
  else if(BaseNum==5)
  V<< 1, 0,
    0, 1,
    y, x,
    2*x, 0,
    0, 2*y;
#elif DIM ==3
  if(BaseNum==3)
  V<< 1, 0, 0,
    0, 1, 0,
    0, 0, 1.; //2026.2.18
#endif

  return V;
} 


inline Vector<Real, UNum> Cell_avg(Matrix<Real, UNum, BaseNum+1> Coef, my_real ps) //20260620
{
  Vector<Real,UNum> Var = Coef.col(0);
  // Real dx= ps.i;
  // Real dy= ps.j;

  // if(BaseNum >=4)
  //   Var += pow(dx, 2)/12 * Coef.col(3);

  // if(BaseNum >=5)
  //   Var += pow(dy, 2)/12 * Coef.col(4);
  
  return Var;
} 


inline Real Limiter(Real ds, Real h, int bd_fg ) //2026.7.23
{

  // Real val = (abs(ds) + ds )/(abs(ds) + 1 ) ;


  Real  val = 1;
  if(ds>0 && ds<0.5) val = 2* ds;
  else if(ds>=0.5 ) val = 1; //2026.7.23
  if(bd_fg == 1)
  {
    val = (abs(ds) + ds )/(abs(ds) + 1 ) ;
  }

  
  return val;
} 


template<class T> void Write_file(ofstream &out, T ss)
{
  out.write((char*)&(ss), sizeof(ss));
};//2026.2.23

template<class T> T Convert_T(T num)
{
  return T(num);
};



inline Real RhoLs_cal(Real rhol, Real pl, Real pstar, Real gam, Real Pp0)
{

  Real valN = (gam + 1)* pstar + (gam - 1)* pl + 2* gam * Pp0; 
  Real valD = (gam - 1)* pstar + (gam + 1)* pl + 2* gam * Pp0; 

  return rhol* valN/ valD; 
} //20260611

inline Real Wl_cal(Real rhol, Real pl, Real pstar, Real gam, Real Pp0)
{
  Real Wl = sqrt(  0.5*rhol*  ( (gam + 1.0)*pstar  + (gam - 1.0)*pl + 2.0*gam * Pp0  )  );

  return Wl;
} //20260611


inline Real Jaco_cal(Real W, Real pl, Real rhol, Real pstar, Real gam, Real Pp0)
{

  Real dK1dp = 1./W - (pstar - pl)/pow(W,2)* rhol* (gam + 1)/4./W;

  return dK1dp;
} //20260611


#endif