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
 * \file StretchingCosSmlCtr.H
 *
 * \brief Provides grid stretching based on the cosine function
 *
 *//*+*************************************************************************/

#ifndef _STRETCHINGCOSSMLCTR_H_
#define _STRETCHINGCOSSMLCTR_H_

//----- Standard Library -----//

#include <cmath>

//----- Chombo -----//

#include "REAL.H"
#include "CONSTANTS.H"

#include "NamespaceHeader.H"


/*******************************************************************************
 */
/// StretchingCosSmlCtr
/**
 *  1-D mapping class that implements stretching based on a cosine function
 *
 *  Typically used to concentrate grid points at the center of a domain.  The
 *  advantage of the cos function is the ability to have constant cell sizes
 *  at the domain extents and center.
 *
 *  Given a variable xi on a domain [0,1], this class provides x also on a
 *  domain [0,1].  Whereas xi is expected to be on a Cartesian grid (uniform
 *  distribution), x is mapped based on a cosine function.
 *
 *  \note
 *  <ul>
 *    <li> Use the accompanying wxmaxima file to explore the effect of
 *         parameters on your function
 *    <li> This function can be closely approximated with a polynomial which is
 *         probably faster to evaluate.  But having about 11 terms is unwieldy
 *         to write and alpha must still be compile-time constant.
 *  </ul>
 *
 *//*+*************************************************************************/

template <int alpha>
class StretchingCosSmlCtr
{
public:

  /// Constructor
  /** \param[in]  a_delta Ratio of smallest cell to largest cell
   */
  StretchingCosSmlCtr(const Real a_delta)
    :
    m_delta(a_delta)
    { }

  /// Get physical coordinate from computational coordinate
  /** \param[in]  a_xi    Computational coordinate
   */
  Real x(const Real a_xi) const noexcept;

  /// Get derivative dx/dxi at computational coordinate
  /** \param[in]  a_xi    Computational coordinate
   */
  Real dxdxi(const Real a_xi) const noexcept;

protected:

  const Real m_delta;
};

//--For alpha = 1

template <>
inline Real
StretchingCosSmlCtr<1>::x(const Real a_xi) const noexcept
{
  return (1. - m_delta)*std::sin(2*Pi*a_xi)/(2*Pi*(m_delta + 1.)) + a_xi;
}
template <>
inline Real
StretchingCosSmlCtr<1>::dxdxi(const Real a_xi) const noexcept
{
  return (1. - m_delta)*std::cos(2*Pi*a_xi)/(m_delta + 1.) + 1.;
}

//--For alpha = 2

template <>
inline Real
StretchingCosSmlCtr<2>::x(const Real a_xi) const noexcept
{
  return (1. - m_delta)*(std::sin(4*Pi*a_xi) + 8*std::sin(2*Pi*a_xi))/
    (4*Pi*(5*m_delta + 3.)) + a_xi;
}
template <>
inline Real
StretchingCosSmlCtr<2>::dxdxi(const Real a_xi) const noexcept
{
  return (1. - m_delta)*(std::cos(4*Pi*a_xi) + 4*std::cos(2*Pi*a_xi))/
    (5*m_delta + 3.) + 1.;
}

//--For alpha = 3

template <>
inline Real
StretchingCosSmlCtr<3>::x(const Real a_xi) const noexcept
{
  const Real c1 = std::sin(2*Pi*a_xi);
  return (1. - m_delta)*(9*std::sin(4*Pi*a_xi) + 48*c1 - 4*c1*c1*c1)/
    (12*Pi*(11*m_delta + 5.)) + a_xi;
}
template <>
inline Real
StretchingCosSmlCtr<3>::dxdxi(const Real a_xi) const noexcept
{
  const Real c1 = std::cos(2*Pi*a_xi);
  return (1. - m_delta)*(3*std::cos(4*Pi*a_xi) + 6*c1 + 2*c1*c1*c1)/
    (11*m_delta + 5.) + 1.;
}

#include "NamespaceFooter.H"

#endif  /* ! defined _STRETCHINGCOSSMLCTR_H_ */
