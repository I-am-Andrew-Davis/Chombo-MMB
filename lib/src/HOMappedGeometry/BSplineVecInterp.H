#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

/******************************************************************************/
/**
 * \file BSplineVecInterp.cpp
 *
 * \brief Class for interpolation between points in array
 *
 *//*+*************************************************************************/

#ifndef _BSPLINEVECINTERP_H_
#define _BSPLINEVECINTERP_H_

#include "REAL.H"
#include "CH_Timer.H"
#include "SPACE.H"
#include "StcVector.H"
#include "CHMatrixOps.H"
#include <vector>

#include "NamespaceHeader.H"

/*******************************************************************************
 */
/// Set of functions needed for general basis spline creation and evaluation
/**
 ******************************************************************************/
namespace BSpline
{
  /// solve for control point derivitives
  void derCtrlPts(Real*                    a_dCP,
                  const Real*              a_ctrlPts,
                  const int                c_str,
                  const int                a_order,
                  const int                a_deriv,
                  const std::vector<Real>& a_T,
                  const int                beg,
                  const Real               a_dx);

  /// solve interpolated value from a 1D set of control points
  template <int Order>
  inline Real
  evalPoint(const Real*              a_ctrlPts,
            const int                cp_str,
            const Real *const        a_bf)
  {
    CH_TIME("BSpline::evalPoint");
    return stc::sumEachElement<Order+1>(
      [&]
      (const stc::array_size_type i)
        {
          return a_ctrlPts[i*cp_str]*a_bf[i];
        });
  }

  /// solve interpolated derivitive from a 1D set of control points
  template <int Order>
  Real evalPointDeriv(const Real               t,
                      const int                a_deriv,
                      const int                a_Tidx,
                      const std::vector<Real>& a_T,
                      const Real               a_dx,
                      const Real*              a_ctrlPts,
                      const int                cp_str,
                      const Real *const        a_bf,
                      Real*                    a_dCP)
  {
    CH_TIME("BSpline::evalPointDeriv");
    CH_assert(a_deriv <= Order);
    Real y = 0;
    int cp_size = Order + 1 - a_deriv;
    derCtrlPts(a_dCP, a_ctrlPts, cp_str, Order, a_deriv,
               a_T, a_Tidx-cp_size, a_dx);
    for(int i=0; i!=cp_size; i++)
      {
        y += a_dCP[i]*a_bf[i];
      }
    return y;
  }

  template <int Order, int Deriv = -1>
  struct PointDeriv
  {
    static Real eval(const int                a_Tidx,
                     const std::vector<Real>& a_T,
                     const Real               a_dx,
                     const Real*              a_ctrlPts,
                     const int                cp_str,
                     const Real *const        a_bf,
                     const int                a_deriv)
      {
        CH_assert(a_deriv < Order);
        stc::RVec<Order+1> dCP;
        const int cp_size = Order + 1 - a_deriv;
        derCtrlPts(dCP.dataPtr(), a_ctrlPts, cp_str, Order, a_deriv, a_T,
                   a_Tidx - cp_size, a_dx);
        Real y = 0;
        for(int i = 0; i != cp_size; ++i)
          {
            y += dCP[i]*a_bf[i];
          }
        return y;
      }
  };

  template <int Order>
  struct PointDeriv<Order, 0>
  {
    static Real eval(const int                a_Tidx,
                     const std::vector<Real>& a_T,
                     const Real               a_dx,
                     const Real*              a_ctrlPts,
                     const int                cp_str,
                     const Real *const        a_bf)
      {
        return evalPoint<Order>(a_ctrlPts, cp_str, a_bf);
      }
  };

  template <int Order>
  struct PointDeriv<Order, 1>
  {
    static Real eval(const int                a_Tidx,
                     const std::vector<Real>& a_T,
                     const Real               a_dx,
                     const Real*              a_ctrlPts,
                     const int                cp_str,
                     const Real *const        a_bf)
      {
        constexpr int r = Order + 1;
        Real dCtrlPts[r];
        // Zeroth deriv
        for(int i = 0; i != r; ++i)
          {
            dCtrlPts[i] = a_ctrlPts[i*cp_str];
          }
        // The desired (highest) deriv
        return stc::sumEachElement<Order>(
          [&]
          (const int i)
            {
              return a_bf[i]*Order*(dCtrlPts[i+1] - dCtrlPts[i])/
                (a_dx*(a_T[a_Tidx + 1 + i] - a_T[a_Tidx + 1 + i - Order]));
            });
      }
  };

  /// solve basis functions, any order
  template <int Order>
  inline void
  basisFuncs(Real*                    bf,
             const int                i,
             const Real               t,
             const std::vector<Real>& T)
  {
    CH_TIME("BSpline::basisFuncs");
    // some temp values for basis func evals
    Real lt[Order+1], rt[Order+1];
    Real saved = 0.0;
    Real temp = 0.0;
    // eval basis funcs
    bf[0] = 1.0;
    for(int j=1; j!=Order+1; j++)
      {
        lt[j] = t - T[i+1-j];
        rt[j] = T[i+j] - t;
        saved = 0.0;
        for(int r=0; r!=j; r++)
          {
            temp = bf[r]/(rt[r+1]+lt[j-r]);
            bf[r] = saved + rt[r+1]*temp;
            saved = lt[j-r]*temp;
          }
        bf[j] = saved;
      }
  }

  template <int Order, int Deriv = -1>
  struct BasisFuncs
  {
  private:
    // Dispatch if Deriv is not compile-time
    static void eval(std::true_type,
                     Real*                    bf,
                     const int                i,
                     const Real               t,
                     const std::vector<Real>& T,
                     const int                a_deriv)
      {
        CH_TIME("BSpline::BasisFuncs::eval(run-time)");
        // some temp values for basis func evals
        Real lt[Order+1], rt[Order+1];
        const int d = Order - a_deriv;
        Real saved = 0.0;
        Real temp = 0.0;
        // eval basis funcs
        bf[0] = 1.0;
        for(int j = 1; j != d+1; ++j)
          {
            lt[j] = t - T[i+1-j];
            rt[j] = T[i+j] - t;
            saved = 0.0;
            for(int r = 0; r != j; ++r)
              {
                temp = bf[r]/(rt[r+1]+lt[j-r]);
                bf[r] = saved + rt[r+1]*temp;
                saved = lt[j-r]*temp;
              }
            bf[j] = saved;
          }
      }

    // Dispatch if Deriv is compile-time 0 or 1
    static void eval(std::false_type,
                     Real*                    bf,
                     const int                i,
                     const Real               t,
                     const std::vector<Real>& T,
                     const int                a_deriv)
      {
        CH_TIME("BSpline::BasisFuncs::eval(compile-time)");
        // some temp values for basis func evals
        constexpr int d = Order - Deriv;
        Real lt[d+1], rt[d+1];
        Real saved = 0.0;
        Real temp = 0.0;
        // eval basis funcs
        bf[0] = 1.0;
        for(int j = 1; j != d+1; ++j)
          {
            lt[j] = t - T[i+1-j];
            rt[j] = T[i+j] - t;
            saved = 0.0;
            for(int r = 0; r != j; ++r)
              {
                temp = bf[r]/(rt[r+1]+lt[j-r]);
                bf[r] = saved + rt[r+1]*temp;
                saved = lt[j-r]*temp;
              }
            bf[j] = saved;
          }
      }

  public:
    static void eval(Real*                    bf,
                     const int                i,
                     const Real               t,
                     const std::vector<Real>& T,
                     const int a_deriv = -1)
      {
        return eval(std::integral_constant<bool, Deriv == -1>{},
                    bf, i, t, T, a_deriv);
      }
  };

  /// solve derivitive at ends
  void derivOneSided(Real&       hi_deriv,
                     Real&       lo_deriv,
                     const Real* U,
                     const int   size,
                     const int   stride,
                     const Real  dx);

  /// solve second derivitive at ends
  void secDerivOneSided(Real&       hi_deriv,
                        Real&       lo_deriv,
                        const Real* U,
                        const int   size,
                        const int   stride,
                        const Real  dx);

  /// solve for control points in a 1D problem
  template <int Order>
  void solveCP(Real*                    a_bf,
               const Real*              Q,
               Real*                    P,
               const int                q_size,
               const int                a_order,
               const Real               a_dx,
               const std::vector<Real>& a_T,
               const int                q_str=1,
               const int                p_str=1)
  {
    CH_TIME("BSpline::solveCP");
    // matrices for solving CP
    int n = q_size - 1; // number of interpolation elements
    int pn = n + a_order - 1; // number of elements of P
    int in = q_size - 2; // number of interior cp to solve 
    CHVector R(in); // stores modified interpolation points
    int kul = (a_order-1)/2;
    int L = 3*kul+1;
    int dI = L-kul-1;
    CHMatrix B(L, in); // banded matrix of basis funcions at interp points
    //std::vector<Real> bf(a_order);
    for(int i=0; i!=in; i++)
      {
        // set R vector equal to Q
        R(i) = Q[(i+1)*q_str];
        // basis function eval at knot
        basisFuncs<Order>(a_bf, a_order+1+i, a_T[a_order+1+i], a_T);
        // fill in values if valid
        for(int j=0; j!=a_order; j++)
          {
            int sft = j-(a_order-1)/2;
            if((i+sft >= 0) && (i+sft < in))
              {
                B(dI-sft, i+sft) = a_bf[j];
              }
          }
      }
  
    // solve P end point given boundary conditions
    // end point interpolation, linear
    P[0] = Q[0];
    P[pn*p_str] = Q[n*q_str];
    Real PdBeg, PdEnd;
    Real D_hi, D_lo; 
    // end derivatives, cubic
    if(a_order >= 3)
      {
        derivOneSided(D_hi, D_lo, Q, q_size, q_str, a_dx);
        PdBeg = P[0] + (a_dx/a_order)*D_lo;
        PdEnd = P[pn*p_str] - (a_dx/a_order)*D_hi;
      }
    // end second derivatives, quintic
    if(a_order >= 5)
      {
        secDerivOneSided(D_hi, D_lo, Q, q_size, q_str, a_dx);
        P[2*p_str] = D_lo*a_dx*a_dx/10.0 + 3.0*PdBeg - 2.0*P[0];
        P[(pn-2)*p_str] = D_hi*a_dx*a_dx/10.0
          + 3.0*PdEnd - 2.0*P[pn*p_str];
      }
    // update P values after the calcuation, allows for inplace solve
    if(a_order >= 3)
      {
        P[1*p_str] = PdBeg;
        P[(pn-1)*p_str] = PdEnd;
      }
  
    // adjust R but solved end derivatives
    // lo side
    if(a_order >= 3) // cubic
      {
        basisFuncs<Order>(a_bf, a_order+1, a_T[a_order+1], a_T);
        R(0) -= a_bf[0]*P[1*p_str];
      }
    if(a_order >= 5) // quartic
      {
        R(0) -= a_bf[1]*P[2*p_str];
        basisFuncs<Order>(a_bf, a_order+2, a_T[a_order+2], a_T);
        R(1) -= a_bf[0]*P[2*p_str];
      }

    // hi side
    if(a_order >= 3) // cubic
      {
        basisFuncs<Order>(a_bf, pn, a_T[pn], a_T);
        R(in-1) -= a_bf[a_order-1]*P[(pn-1)*p_str];
      }
    if(a_order >= 5) // quartic
      {
        R(in-1) -= a_bf[a_order-2]*P[(pn-2)*p_str];
        basisFuncs<Order>(a_bf, pn-1, a_T[pn-1], a_T);
        R(in-2) -= a_bf[a_order-1]*P[(pn-2)*p_str];
      }
    
    // solve digonal system of eqs
    CHgbsv(B, kul, kul, R);
    // copy appropriate portion into P
    for(int i=0; i!=in; i++) P[(i+(a_order-1)/2+1)*p_str] = R(i); 
  }
}


/*******************************************************************************
 */
///  One dimensional b-spline interpolation between data of an array
///  Data must be uniformly spaced in current implementation
///  Supports use of linear(1), cubic(3), and quintic(5) interpolation
///  Solves for continuous fit through first m_order-1 derivatives
///  providing for O( x^(m_order+1) )
///  Requires minimum of 6 points for cubic, 7 for quintic
/**
 *   // Example usage with set of x,y points
 *     int size = 6;
 *     Real x[size] = {1, 2, 3, 4, 5, 6};       // don't actually need this data
 *     Real y[size] = {1, 8, 27, 64, 125, 216}; 
 *     Real dx = 1.0;       // x spacing
 *     int stride = 1;      // y array indexing stride
 *     Real xstart = x[0];, // first 
 *     int polyOrder = 3;   // create a cubic polynomial
 * 
 *     BSplineVecInterp Spline(y, size, xstart, dx, polyOrder, stride);
 *   // solve y point at x = 2.4
 *     Real pt1 = Spline.interp(2.4);
 *
 *   // Same example using std::vector
 *     std::vector<Real> y_vec = {1, 8, 27, 64, 125, 216};
 *     BSplineVecInterp SplineV(y_vec.data(), 
 *                              y_vec.size(), 
 *                              x_start, 
 *                              dx, 
 *                              polyOrder,
 *                              1);     // stride is 1 for std::vectors 
 *   // solve y point at x = 3.7
 *     Real pt2 = SplineV.interp(3.7);
 *
 *  \tparam[Order]      Polynomial order (must be odd, default 5)
 */
/**
 ******************************************************************************/
template <int Order = 5>
class BSplineVecInterp
{
/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:
  /// default constructor
  BSplineVecInterp();
  
  /// constructor
  BSplineVecInterp(const Real* a_interp,
                   const int   a_size,
                   const Real  a_xstart,
                   const Real  a_dx,
                   const int   a_stride = 1);
  
  /// weak constructor
  void define(const Real* a_interp,
              const int   a_size,
              const Real  a_xstart,
              const Real  a_dx,
              const int   a_stride = 1);
  
  /// copy constructor disabled
  BSplineVecInterp(const BSplineVecInterp&) = delete;

  /// move constructor disabled
  BSplineVecInterp(BSplineVecInterp&&) = delete;

  /// assignment constructor disabled
  BSplineVecInterp& operator=(const BSplineVecInterp&) = delete;

  /// move assignment constructor disabled
  BSplineVecInterp& operator=(BSplineVecInterp&&) = delete;

  /// destructor
  ~BSplineVecInterp();


/*==============================================================================
 * Public member functions
 *============================================================================*/

  /// interpolate at given point
  Real interp(const Real a_xi) const;

  /// derivitive at given point
  Real interpD(const Real a_xi, const int a_deriv=1) const;


/*==============================================================================
 * Protected member functions
 *============================================================================*/
protected: 

  /// fill knot vector
  void knots(std::vector<Real>& T, Real& dT, const int size);

  /// get index in knot vector for any t
  int getIndx(const Real t, const std::vector<Real>& T, const Real dT) const;
  

/*==============================================================================
 * Data Members
 *============================================================================*/
  
protected:
  const Real*       m_interp;   ///< pointer to interpolation data
  int               m_size;     ///< number of interp points
  int               m_str;      ///< data stride ( 1 if contiguous )
  std::vector<Real> m_ctrlPts;  ///< store control points used for interpolation
  std::vector<Real> m_T;        ///< knot vector
  Real              m_dT;       ///< uniform knot spacing
  Real              m_dx;       ///< physical spacing
  Real              m_x0;       ///< physical start location
  Real*             m_bf;       ///< storage for basis function evaluation
  Real*             m_dCP;      ///< storage for control point derivatives
};


#include "NamespaceFooter.H"
#endif
