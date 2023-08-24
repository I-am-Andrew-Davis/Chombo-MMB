#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _ROOTSOLVER_H_
#define _ROOTSOLVER_H_

/******************************************************************************/
/**
 * \file
 *
 * \brief Solvers for numerical integration
 *
 *//*+*************************************************************************/

#include <cmath>
#include <algorithm>

#include "CH_assert.H"

#include "BaseNamespaceHeader.H"

namespace NumericIntegral
{
  template <typename T> struct NITr
  {
  };
  // Default epsilon and tolerance for floats
  template <> struct NITr<float>
  {
    enum
    {
      maxIter = 100
    };

    static float eps()
      {
        return 3.0E-7;
      }

    static float tolerance()
      {
        return 1.0e-6;
      }
  };

// Default epsilon and tolerance for doubles
  template <> struct NITr<double>
  {
    enum
    {
      maxIter = 100
    };

    static double eps()
      {
        return 3.0E-15;
      }

    static double tolerance()
      {
        return 1.0e-12;
      }
  };
  
/** Adaptive Simpson's Rule, Recursive function */
  template <typename T, typename Func>
  T adaptiveSimpsonsAux(const Func& f,
                        T           a,
                        T           b,
                        T           tol,
                        T           whole,
                        T           fa,
                        T           fb,
                        T           fm,
                        int         rec)
  {
    T m   = (a + b)/2,  h   = (b - a)/2;
    T lm  = (a + m)/2,  rm  = (m + b)/2;
    // serious numerical trouble: it won't converge
    if ((tol/2 == tol) || (a == lm))
      {
        std::cerr  << "NumericIntegral::adaptiveSimpsons: convergence failed"
                   << std::endl;
        return whole;
      }
    T flm = f(lm),      frm = f(rm);
    T left  = (h/6) * (fa + 4*flm + fm);
    T right = (h/6) * (fm + 4*frm + fb);
    T delta = left + right - whole;

    // if (rec <= 0)
    //   {
    //     std::cerr  << "NumericIntegral::adaptiveSimpsons: maximum number of iterations reached"
    //                << std::endl;
    //   }
    // Lyness 1969 + Richardson extrapolation; see article
    if (rec <= 0 || std::abs(delta) <= 15.*tol)
      {
        return left + right + (delta)/15.;
      }
    return (adaptiveSimpsonsAux(f, a, m, tol/2., left,  fa, fm, flm, rec-1) +
            adaptiveSimpsonsAux(f, m, b, tol/2., right, fm, fb, frm, rec-1));
  }

/** Adaptive Simpson's Rule Wrapper
 *  (fills in cached function evaluations) */
  template <typename T, typename Func>
  T adaptiveSimpsons(const Func& f,     // function ptr to integrate
                     T           a,
                     T           b,               // interval [a,b]
                     T           tol = NITr<T>::tolerance(),  // error tolerance
                     int         maxIter = 1000)
  { // recursion cap
    CH_assert(b > a);
    T h = b - a;
    if (h == 0)
      {
        return 0;
      }
    T fa = f(a), fb = f(b), fm = f((a + b)/2);
    T S = (h/6)*(fa + 4*fm + fb);
    // now that its set up, call the recursive solver
    return adaptiveSimpsonsAux(f, a, b, tol, S, fa, fb, fm, maxIter);
  }
  
}  // End of namespace NumericIntegration

#include "BaseNamespaceFooter.H"

#endif
