#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef __FACTORIAL_H_
#define __FACTORIAL_H_

#include "MayDay.H"
#include "RealVect.H"
#include "IndexTM.H"
#include "NamespaceHeader.H"

/**
 *  These helper functions calculate factorials, binomial coeff, etc.
 *  for multinomials.
 */

#ifdef USE_STCVECTOR
using Sz = stc::array_size_type;
#else
using Sz = int;
#endif

/// Calculates factorial for an integer
inline Real factorial(const int n)
{
    CH_assert(n >= 0);
    Real nfact = 1.0;
    for (int i = 2; i <= n; nfact*= i, i++);
    return nfact;
}
/// computes x^p
inline Real POW(const Real& a_x,  const int& a_p)
{
    Real retval = 1;
    for(int iexp = 0; iexp < Abs(a_p); iexp++)
      {
        if(a_p >= 0)
          {
            retval *= a_x;
          }
        else
          {
            retval /= a_x;
          }
      }
    return retval;
}

/// Calculates the binomial coefficient, "n choose k" 
inline Real nCk(const int n, const int k)
{
    CH_assert((n >= k)&&(k >= 0)&&(n >= 0));
    Real nck = 1.0;
    for (int i = n-k+1; i <= n; nck*= i, i++);
    return nck / factorial(k);
}

/// Calculates factorials for a multinomial
template <Sz Dim>
inline Real pfactorial(const IndexTM<int, Dim>& p)
{
    Real pfact = 1;
    for(int idir = 0; idir < Dim; idir++)
      {
        CH_assert(p[idir] >= 0);
        pfact *=  factorial(p[idir]); 
      }
    return pfact;
}

/// Calculates the multinomial coefficient, "p choose k" 
template <Sz Dim>
inline Real pCk(const IndexTM<int, Dim>& p, const IndexTM<int, Dim>& k)
{
  //CH_assert((p >= IntVect::Zero)&&(k >= IntVect::Zero)&&(p >= k));

    Real pfact = 1;
    for(int idir = 0; idir < Dim; idir++)
      {
        pfact *= nCk(p[idir],k[idir]); 
      }
    return pfact;
}


/// calculate x^p
template <Sz Dim>
inline Real power(const IndexTM<Real, Dim>& a_x, const IndexTM<int, Dim>& a_p)
{
  Real retval = 1;
  for(int idir = 0; idir < Dim; idir++)
    {
      if((a_p[idir] >= 0))
        {
          for(int ipow = 0; ipow < a_p[idir]; ipow++)
            {
              retval *= a_x[idir];
            }
        }
      else
        {
          for(int ipow = 0; ipow < a_p[idir]; ipow++)
            {
              retval /= a_x[idir];
            }
        }
    }
    return retval;
}


#include "NamespaceFooter.H"
#endif
