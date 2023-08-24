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
 * \file StretchingTanh.H
 *
 * \brief Provides grid stretching based on the tanh function
 *
 *//*+*************************************************************************/

#ifndef _STRETCHINGTANH_H_
#define _STRETCHINGTANH_H_

//----- Standard Library -----//

#include <cmath>

//----- Chombo -----//

#include "REAL.H"
#include "CONSTANTS.H"

#include "NamespaceHeader.H"


/*******************************************************************************
 */
/// StretchingTanh
/**
 *  1-D mapping class that implements stretching based on tanh
 *
 *  Typically used to concentrate grid points at one extent of a domain.  The
 *  advantage of the tanh function is the ability to have constant cell sizes
 *  at the domain extents and transitions in the center.
 *
 *  Given a variable xi on a domain [0,1], this class provides x also on a
 *  domain [0,1].  Whereas xi is expected to be on a Cartesian grid (uniform
 *  distribution), x is mapped based on a tanh function.
 *
 *  \note
 *  <ul>
 *    <li> Use the accompanying wxmaxima file to explore the effect of
 *         parameters on your function
 *  </ul>
 *
 *//*+*************************************************************************/

class StretchingTanh
{
public:

  /// Constructor
  /** The default function transitions from small cell to large cell over the
   *  range 0.2 < xi < 0.8.  Setting a_alpha to a larger value compresses the
   *  range of the transition.  The left and right bounds can be adjusted using
   *  a_blo and a_bhi to offset the transition to one side.
   *  \param[in]  a_delta Ratio of smallest cell to largest cell
   *  \param[in]  a_alpha Compression of range of transition.  Set to negative
   *                      to reverse small and large cells (default 1)
   *  \param[in]  a_blo   Lower boundary (default 0)
   *  \param[in]  a_bhi   Upper boundary (default 1)
   */
  StretchingTanh(const Real a_delta,
                 const Real a_alpha = 1.,
                 const Real a_blo   = 0.,
                 const Real a_bhi   = 1.)
    :
    m_delta(a_delta),
    m_alpha(a_alpha),
    m_blo(a_blo),
    m_bhi(a_bhi),
    m_cxi1(2*Pi*a_alpha*(a_blo - a_bhi)),
    m_cxi2(2*Pi*a_alpha*(a_delta + 1.)*(a_blo - a_bhi)),
    m_cnum1((1 - a_delta)*std::log(std::cosh(Pi*a_alpha*(2*a_blo - 1.)))),
    m_den((a_delta - 1)*
          std::log(std::cosh(Pi*a_alpha*(2*a_blo - 1.))/
                   std::cosh(Pi*a_alpha*(2*a_bhi - 1.))) - m_cxi2)
    { }

  /// Get physical coordinate from computational coordinate
  /** \param[in]  a_xi    Computational coordinate
   */
  Real x(const Real a_xi) const noexcept
    {
      const Real cxi1 = std::cosh((m_cxi1)*a_xi - Pi*m_alpha*(2*m_blo - 1.));
      return ((1. - m_delta)*std::log(cxi1) - m_cxi2*a_xi - m_cnum1)/m_den;
    }

  /// Get derivative dx/dxi at computational coordinate
  /** \param[in]  a_xi    Computational coordinate
   */
  Real dxdxi(const Real a_xi) const noexcept
    {
      const Real cxi2 = m_cxi1*a_xi - Pi*m_alpha*(2*m_blo - 1.);
      const Real cxi3 = std::cosh(cxi2);
      return ((1. - m_delta)*m_cxi1*std::sinh(cxi2) - m_cxi2*cxi3)/(m_den*cxi3);
    }

protected:

  const Real m_delta;
  const Real m_alpha;
  const Real m_blo;
  const Real m_bhi;
  // Cached coefficients
  const Real m_cxi1;
  const Real m_cxi2;
  const Real m_cnum1;
  const Real m_den;
};

#include "NamespaceFooter.H"

#endif  /* ! defined _STRETCHINGTANH_H_ */
