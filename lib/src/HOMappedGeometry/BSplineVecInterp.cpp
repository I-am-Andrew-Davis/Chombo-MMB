#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "BSplineVecInterp.H"
#include "CHMatrixOps.H"
#include "math.h"

#include "CH_Timer.H"

#include "NamespaceHeader.H"
using namespace BSpline;
/******************************************************************************/
/**
 * \file BSplineVecInterp.cpp
 *
 * \brief Class Non-inline definitions for classes in bsplineVecInterp.H
 *
 *//*+*************************************************************************/

/*--------------------------------------------------------------------*/
//  Default constructor
/*--------------------------------------------------------------------*/
template <int Order>
BSplineVecInterp<Order>::BSplineVecInterp()
{
}

/*--------------------------------------------------------------------*/
//  destructor
/*--------------------------------------------------------------------*/
template <int Order>
BSplineVecInterp<Order>::~BSplineVecInterp()
{
}

/*--------------------------------------------------------------------*/
//  Constructor with a baseFab
/** \param[in]  a_interp  array of point to interpolate between
 *  \param[in]  a_size    number of element in array
 *  \param[in]  a_xstart  initial x location
 *  \param[in]  a_dx      x spacing
 *  \param[in]  a_stride  spacing ( def-ult 1 for contiguous mem)
 *//*-----------------------------------------------------------------*/
template <int Order>
BSplineVecInterp<Order>::BSplineVecInterp(const Real* a_interp,
                                          const int   a_size,
                                          const Real  a_xstart,
                                          const Real  a_dx,
                                          const int   a_stride)
{
  define(a_interp, a_size, a_xstart, a_dx, a_stride);
}

/*--------------------------------------------------------------------*/
//  Weak constructor with a baseFab
/** \param[in]  a_interp  array of point to interpolate between
 *  \param[in]  a_size    number of element in array
 *  \param[in]  a_xstart  initial x location
 *  \param[in]  a_dx      x spacing
 *  \param[in]  a_stride  spacing ( default 1 for contiguous mem)
 *//*-----------------------------------------------------------------*/
template <int Order>
void
BSplineVecInterp<Order>::define(const Real* a_interp,
                                const int   a_size,
                                const Real  a_xstart,
                                const Real  a_dx,
                                const int   a_stride)
{
  CH_assert(a_size >= 7); // stencil size for ends, should fix this
  static_assert(Order == 1 || Order == 3 || Order == 5,
                "Requires Order = 1,3,5");
  // set spacing
  m_dx = a_dx;
  m_x0 = a_xstart;
  // alias to interpolation data
  m_interp = a_interp;
  m_size = a_size;
  m_str = a_stride;
  // memory for basis funcs
  stc::RVec<Order+1> bf;
  // size the control point vector
  m_ctrlPts.resize(m_size + Order - 1);
  // create knot vector
  knots(m_T, m_dT, m_size);
  // solve the control points
  solveCP<Order>(bf.dataPtr(), m_interp, m_ctrlPts.data(), m_size,
                 Order, m_dx, m_T,
                 m_str, 1);
}

/*--------------------------------------------------------------------*/
//  Interpolate point in domain
/** \param[in]  a_xi     location to interp
 *                       must be within the size of the baseFab box
 *  \return              interpolated value
 *//*-----------------------------------------------------------------*/
template <int Order>
Real
BSplineVecInterp<Order>::interp(const Real a_xi) const
{
  CH_TIME("BSplineVecInterp::interp");
  // check interp point is in bounds
  CH_assert(a_xi >= m_x0 - 0.00000001);
  CH_assert(a_xi <= m_x0 + (m_size-1)*m_dx + 0.00000001);
  // normalize point
  Real t = (a_xi-m_x0)/m_dx;
  // get knot index
  int Tidx = getIndx(t, m_T, m_dT);
  // last point is technically out of bounds, usability fix
  // also allows for extrapolation, but no guarantee what you get
  if((a_xi-m_x0) >= (m_size-1)*m_dx)
    {
      // Tidx = m_T[m_T.size() - Order - 1];
      Tidx = m_T.size() - Order - 2;
    }
  stc::RVec<Order+1> bf;
  basisFuncs<Order>(bf.dataPtr(), Tidx, t, m_T);
  return evalPoint<Order>(&(m_ctrlPts.data()[Tidx-Order]),
                          1, bf.dataPtr());
}

/*--------------------------------------------------------------------*/
//  derivative of interpolated point in domain
/** \param[in]  a_xi     location to interp
 *                       must be within the size of the baseFab box
 *  \param[in]  a_deriv  which derivative to find
 *                       only valid for derivatives < spline order
 *  \return              interpolated value
 *//*-----------------------------------------------------------------*/
template <int Order>
Real
BSplineVecInterp<Order>::interpD(const Real a_xi,
                                 const int  a_deriv) const
{
  CH_TIME("BSplineVecInterp::interpD");
  // check interp point is in bounds
  CH_assert(a_xi >= m_x0 - 0.00000001);
  CH_assert(a_xi <= m_x0 + (m_size-1)*m_dx + 0.00000001);
  CH_assert(a_deriv > 0);
  // normalize point
  Real t = (a_xi-m_x0)/m_dx;
  // get knot index
  int Tidx = getIndx(t, m_T, m_dT);
  // last point is technically out of bounds, usability fix
  if((a_xi-m_x0) >= (m_size-1)*m_dx)
    {
      Tidx = m_T.size() - Order - 2;
    }
  stc::RVec<Order+1> bf;
  stc::RVec<Order+1> dCP;
  switch (a_deriv)
    {
    case 0:
      BSpline::BasisFuncs<Order, 0>::eval(bf.dataPtr(), Tidx, t, m_T);
      return BSpline::PointDeriv<Order, 0>::eval(
        Tidx,
        m_T,
        m_dx,
        &(m_ctrlPts.data()[Tidx-Order]),
        1,
        bf.dataPtr());
    case 1:
      BSpline::BasisFuncs<Order, 1>::eval(bf.dataPtr(), Tidx, t, m_T);
      return BSpline::PointDeriv<Order, 1>::eval(
        Tidx,
        m_T,
        m_dx,
        &(m_ctrlPts.data()[Tidx-Order]),
        1,
        bf.dataPtr());
    default:
      BSpline::BasisFuncs<Order>::eval(bf.dataPtr(), Tidx, t, m_T, a_deriv);
      return BSpline::PointDeriv<Order>::eval(
        Tidx,
        m_T,
        m_dx,
        &(m_ctrlPts.data()[Tidx-Order]),
        1,
        bf.dataPtr(),
        a_deriv);
    }
}

/*--------------------------------------------------------------------*/
//  Get the index of knot vector for a given value t bounded by T
/** \param[in]  t  point to interp index, in range of T  
 *  \param[in]  T  knot vector
 *  \param[in]  dt knot spacing
 *  \return        knot index
 *//*-----------------------------------------------------------------*/
template <int Order>
int
BSplineVecInterp<Order>::getIndx(const Real t,
                          const std::vector<Real>& T,
                          const Real dt) const
{
  return int((t - fmod(t,dt))/dt) + Order;
}

/*--------------------------------------------------------------------*/
//  Create the knot vector, uniformly spaced of unit stride
//  There are repeated knots at the end points as needed
/** \param[out] T     The knot vector to size and fill
 *  \param[out] dT    Knot spacing to set
 *  \param[in]  size  Number of point to interpolate
 *//*-----------------------------------------------------------------*/
template <int Order>
void
BSplineVecInterp<Order>::knots(std::vector<Real>& T, Real& dT, const int size)
{
  // set number of knot vects, then set each vector along FAB directions
  T.resize(size + 2*Order);
  dT = 1.0;
  // fill vector with uniform points
  for(int j=0; j!=Order; j++)
    {
      T[j] = 0.0;
    }
  for(int j=Order; j!=size+Order; j++)
    {
      T[j] = dT*(j-Order);
    }
  for(int j=size+Order; j!=T.size(); j++)
    {
      T[j] = dT*(size - 1);
    }
}



/*******************************************************************************
 */
/// Set of functions needed for general basis spline functionality
/// Part of BSpline namespace
/**
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Interpolate a point at a given location
/** \param[in]  t          Independent variable value to interp at
 *  \param[in]  a_order    polynomial order
 *  \param[in]  a_Tidx     index in knot vect for evaluation
 *  \param[in]  a_T        knot vector
 *  \param[in]  a_ctrlPts  array of control points to use
 *  \param[in]  cp_str     stride of control points
 *  \return                Interpolated y value
 *//*-----------------------------------------------------------------*/
// Real BSpline::evalPoint(const Real               t,
//                         const int                a_order,
//                         const int                a_Tidx,
//                         const std::vector<Real>& a_T,
//                         const Real*              a_ctrlPts,
//                         const int                cp_str,
// //**Ideally make this an stc::Vector
//                         const Real *const        a_bf)
// {
//   CH_TIME("BSpline::evalPoint");
//   Real y = 0;
// //**a_order is constexpr, can unroll loop
//   for(int i=0; i!=a_order+1; i++)
//     {
//       //y += a_ctrlPts[(i+a_Tidx-a_order)*cp_str]*a_bf[i];
//       y += a_ctrlPts[i*cp_str]*a_bf[i];
//     }
//   return y;
// }

/*--------------------------------------------------------------------*/
//  Interpolate a derivative at a given location
/** \param[in]  t          Independent variable value to interp at
 *  \param[in]  a_order    polynomial order
 *  \param[in]  a_Tidx     index in knot vect for evaluation
 *  \param[in]  a_T        knot vector
 *  \param[in]  a_dx       physical spacing
 *  \param[in]  a_ctrlPts  array of control points to use
 *  \param[in]  cp_str     stride of control points
 *  \return                Interpolated y value
 *//*-----------------------------------------------------------------*/
// Real BSpline::evalPointDeriv(const Real               t,
//                              const int                a_deriv,
//                              const int                a_order,
//                              const int                a_Tidx,
//                              const std::vector<Real>& a_T,
//                              const Real               a_dx,
//                              const Real*              a_ctrlPts,
//                              const int                cp_str,
// //**Ideally make this an stc::Vector
//                              const Real *const        a_bf,
//                              Real*                    a_dCP)
// {
//   CH_TIME("BSpline::evalPointDeriv");
//   CH_assert(a_deriv <= a_order);
//   Real y = 0;
//   int cp_size = a_order + 1 - a_deriv;
// //**Precompute just like bf
//   derCtrlPts(a_dCP, a_ctrlPts, cp_str, a_order, a_deriv,
//              a_T, a_Tidx-cp_size, a_dx);
//   for(int i=0; i!=cp_size; i++)
//     {
//       y += a_dCP[i]*a_bf[i];
//     }
//   return y;
// }

/*--------------------------------------------------------------------*/
//  Solve the control points for a 1D slice of m_interp
//  FIXME -- ineffecient copy of data to matrix
//  FIXME -- Uniform knots only
//  FIXME -- Try using not-a-knot end conditions
/** \param[in]  Q         Array containing interpolation points
 *  \param[out] P         Array to store control points
 *  \param[in]  q_size    Number of point in Q
 *  \param[in]  a_order   Polynomial order
 *  \param[in]  a_dx      Uniform x spacing
 *  \param[in]  a_T       Knot vector
 *  \param[in]  q_str     Stride of Q data (defult 1)
 *  \param[in]  p_str     Stride of P data (default 1)
 *//*-----------------------------------------------------------------*/
// template <int Order>
// void BSpline::solveCP(Real*                    a_bf,
//                       const Real* const        Q,
//                       Real* const              P,
//                       const int                q_size,
//                       const int                a_order,
//                       const Real               a_dx,
//                       const std::vector<Real>& a_T,
//                       const int                q_str,
//                       const int                p_str)
// {
//   CH_TIME("BSpline::solveCP");
//   // matrices for solving CP
//   int n = q_size - 1; // number of interpolation elements
//   int pn = n + a_order - 1; // number of elements of P
//   int in = q_size - 2; // number of interior cp to solve 
//   CHVector R(in); // stores modified interpolation points
//   int kul = (a_order-1)/2;
//   int L = 3*kul+1;
//   int dI = L-kul-1;
//   CHMatrix B(L, in); // banded matrix of basis funcions at interp points
//   //std::vector<Real> bf(a_order);
//   for(int i=0; i!=in; i++)
//     {
//       // set R vector equal to Q
//       R(i) = Q[(i+1)*q_str];
//       // basis function eval at knot
//       basisFuncs<Order>(a_bf, a_order+1+i, a_T[a_order+1+i], a_T);
//       // fill in values if valid
//       for(int j=0; j!=a_order; j++)
//         {
//           int sft = j-(a_order-1)/2;
//           if((i+sft >= 0) && (i+sft < in))
//             {
//               B(dI-sft, i+sft) = a_bf[j];
//             }
//         }
//     }
  
//   // solve P end point given boundary conditions
//   // end point interpolation, linear
//   P[0] = Q[0];
//   P[pn*p_str] = Q[n*q_str];
//   Real PdBeg, PdEnd;
//   Real D_hi, D_lo; 
//   // end derivatives, cubic
//   if(a_order >= 3)
//     {
//       derivOneSided(D_hi, D_lo, Q, q_size, q_str, a_dx);
//       PdBeg = P[0] + (a_dx/a_order)*D_lo;
//       PdEnd = P[pn*p_str] - (a_dx/a_order)*D_hi;
//     }
//   // end second derivatives, quintic
//   if(a_order >= 5)
//     {
//       secDerivOneSided(D_hi, D_lo, Q, q_size, q_str, a_dx);
//       P[2*p_str] = D_lo*a_dx*a_dx/10.0 + 3.0*PdBeg - 2.0*P[0];
//       P[(pn-2)*p_str] = D_hi*a_dx*a_dx/10.0
//         + 3.0*PdEnd - 2.0*P[pn*p_str];
//     }
//   // update P values after the calcuation, allows for inplace solve
//   if(a_order >= 3)
//     {
//       P[1*p_str] = PdBeg;
//       P[(pn-1)*p_str] = PdEnd;
//     }
  
//   // adjust R but solved end derivatives
//   // lo side
//   if(a_order >= 3) // cubic
//     {
//       basisFuncs<Order>(a_bf, a_order+1, a_T[a_order+1], a_T);
//       R(0) -= a_bf[0]*P[1*p_str];
//     }
//   if(a_order >= 5) // quartic
//     {
//       R(0) -= a_bf[1]*P[2*p_str];
//       basisFuncs<Order>(a_bf, a_order+2, a_T[a_order+2], a_T);
//       R(1) -= a_bf[0]*P[2*p_str];
//     }

//   // hi side
//   if(a_order >= 3) // cubic
//     {
//       basisFuncs<Order>(a_bf, pn, a_T[pn], a_T);
//       R(in-1) -= a_bf[a_order-1]*P[(pn-1)*p_str];
//     }
//   if(a_order >= 5) // quartic
//     {
//       R(in-1) -= a_bf[a_order-2]*P[(pn-2)*p_str];
//       basisFuncs<Order>(a_bf, pn-1, a_T[pn-1], a_T);
//       R(in-2) -= a_bf[a_order-1]*P[(pn-2)*p_str];
//     }
    
//   // solve digonal system of eqs
//   CHgbsv(B, kul, kul, R);
//   // copy appropriate portion into P
//   for(int i=0; i!=in; i++) P[(i+(a_order-1)/2+1)*p_str] = R(i); 
// }

/*--------------------------------------------------------------------*/
//  Solve for the basis functions at a given point
/** \param[in]  i  knot index  
 *  \param[in]  t  evaluation point, within range of T
 *  \param[in]  T  knot vector
 *  \param[in]  d  order of basis function
 *  \return        vector of basis functions, of size d+1
 *//*-----------------------------------------------------------------*/
// void
// BSpline::basisFuncs(Real*                    bf,
//                     const int                i,
//                     const Real               t,
//                     const std::vector<Real>& T,
//                     const int                d)
// {
//   CH_TIME("BSpline::basisFuncs");
//   // some temp values for basis func evals
//   Real lt[d+1],  rt[d+1];
//   Real saved = 0.0;
//   Real temp = 0.0;
//   // eval basis funcs
//   bf[0] = 1.0;
// //**d is constexpr, can unroll loop
//   for(int j=1; j!=d+1; j++)
//     {
//       lt[j] = t - T[i+1-j];
//       rt[j] = T[i+j] - t;
//       saved = 0.0;
// //**j is constexpr, can unroll loop
//       for(int r=0; r!=j; r++)
//         {
//           temp = bf[r]/(rt[r+1]+lt[j-r]);
//           bf[r] = saved + rt[r+1]*temp;
//           saved = lt[j-r]*temp;
//         }
//       bf[j] = saved;
//     }
// }

/*--------------------------------------------------------------------*/
//  Solve derivative control points, need for spline derivatives
/** \param[out] dCP        control point derivative specified
 *  \param[in]  a_ctrlPts  control points 
 *  \param[in]  cp_str     control point stride
 *  \param[in]  a_order    order of basis func
 *  \param[in]  a_deriv    what deriv to go up to
 *  \param[in]  a_T        knot vector
 *  \param[in]  a_Tidx     knot vector eval point
 *  \param[in]  a_dx       interp point spacing
 *//*-----------------------------------------------------------------*/
void
BSpline::derCtrlPts(Real*                    a_dCP,
                    const Real*              a_ctrlPts,
                    const int                cp_str,
                    const int                a_order,
                    const int                a_deriv,
                    const std::vector<Real>& a_T,
                    const int                a_Tidx,
                    const Real               a_dx)
{
  CH_TIME("BSpline::derCtrlPts");
  CH_assert(a_deriv >= 0);
  CH_assert(a_deriv <= a_order);
  int r = a_order+1;
  if(a_deriv > 0)
    {
      Real dCtrlPts [a_deriv+1][r];
      // zeroth deriv
      for(int i=0; i!=r; i++)
        {
          dCtrlPts[0][i] = a_ctrlPts[i*cp_str];
        }
      // higher derivs
      for(int k=1; k<a_deriv; k++)
        {
          Real factor = r-k;
          for(int i=0; i<r-k; i++)
            {
              dCtrlPts[k][i] = factor*(dCtrlPts[k-1][i+1] - dCtrlPts[k-1][i])
                /(a_dx*(a_T[a_Tidx+i+a_order+1] - a_T[a_Tidx+i+k]));
            }
        }
      // the desired (highest) deriv
      int k = a_deriv;
      Real factor = r-k;
      for(int i=0; i<r-k; i++)
        {
          a_dCP[i] = factor*(dCtrlPts[k-1][i+1] - dCtrlPts[k-1][i])
            /(a_dx*(a_T[a_Tidx+i+a_order+1] - a_T[a_Tidx+i+k]));
        }
    }
  else
    {
      for(int i=0; i!=r; i++)
        {
          a_dCP[i] = a_ctrlPts[i*cp_str];
        }
    }
}

/*--------------------------------------------------------------------*/
//  Solve one sided derivatives at ends of array
/** \param[in]  U      array of point to solve end derivs
 *  \param[in]  size   number of elements in array
 *  \param[in]  stride spacing of element ( 1 if contiguous mem)
 *  \param[in]  dx     spacing between points
 *  \return            two element vector with the derivatives
 *                     Solves using O(x^6) stencils
 *//*-----------------------------------------------------------------*/
void
BSpline::derivOneSided(Real&       hi_deriv,
                       Real&       lo_deriv,
                       const Real* U,
                       const int   size,
                       const int   stride,
                       const Real  dx)
{
  CH_assert(size >= 7); // stencil size for ends
  // hi side
  int s = size-1;
  hi_deriv = (10.0*U[(s-6)*stride] -72.0*U[(s-5)*stride]
              +225.0*U[(s-4)*stride] -400.0*U[(s-3)*stride]
              +450.0*U[(s-2)*stride] -360.0*U[(s-1)*stride]
              +147.0*U[(s-0)*stride])/(60.0*dx);
  // lo side
  s = 0;
  lo_deriv = (-147.0*U[(s+0)*stride] +360.0*U[(s+1)*stride]
              -450.0*U[(s+2)*stride] +400.0*U[(s+3)*stride]
              -225.0*U[(s+4)*stride] +72.0*U[(s+5)*stride]
              -10.0*U[(s+6)*stride])/(60.0*dx);
}

/*--------------------------------------------------------------------*/
//  Solve one sided second derivatives at ends of array
/** \param[in]  U      array of point to solve end derivs
 *  \param[in]  size   number of elements in array
 *  \param[in]  stride spacing of element ( 1 if contiguous mem)
 *  \param[in]  dx     spacing between points
 *  \return            two element vector with the derivatives
 *                     Solves using O(x^6) stencils
 *//*-----------------------------------------------------------------*/
void
BSpline::secDerivOneSided(Real&       hi_deriv,
                          Real&       lo_deriv,
                          const Real* U,
                          const int   size,
                          const int   stride,
                          const Real  dx)
{
  CH_assert(size >= 8); // stencil size for ends
  // hi side
  int s = size-1;
  hi_deriv = (- 126.0*U[(s-7)*stride] +1019.0*U[(s-6)*stride]
              -3618.0*U[(s-5)*stride] +7380.0*U[(s-4)*stride]
              -9490.0*U[(s-3)*stride] +7911.0*U[(s-2)*stride]
              -4014.0*U[(s-1)*stride] + 938.0*U[(s-0)*stride])/(180.0*dx*dx);
  // lo side
  s = 0;
  lo_deriv = (  938.0*U[(s+0)*stride] -4014.0*U[(s+1)*stride]
                +7911.0*U[(s+2)*stride] -9490.0*U[(s+3)*stride]
                +7380.0*U[(s+4)*stride] -3618.0*U[(s+5)*stride]
                +1019.0*U[(s+6)*stride] - 126.0*U[(s+7)*stride])/(180.0*dx*dx);
}


/*******************************************************************************
 *
 * Explicit instantiations of class BSplineVecInterp
 *
 ******************************************************************************/

template class BSplineVecInterp<1>;
template class BSplineVecInterp<3>;
template class BSplineVecInterp<5>;

#include "NamespaceFooter.H"
