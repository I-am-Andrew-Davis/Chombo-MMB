#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LOGSCHWARZCHRISTOFFELRAMPCS_H_
#define _LOGSCHWARZCHRISTOFFELRAMPCS_H_

#include "CONSTANTS.H"
#include "NewFourthOrderCoordSys.H"
#include "NamespaceHeader.H"


/*******************************************************************************
 */
/// Conformal Schwarz-Christoffel mapping of a ramp with log stretching
/**                                                                    \verbatim
 *  Physical-space (z)                Computational Space
 *
 *      ^                               ^
 *    i |      ,-+                    i |
 *         ,--'  |
 *     ,--'       \                     +-----------+
 *  +-'          ,+                     |           |
 *   \        ,-'                       |           |
 *    |    ,-'  alpha                   |           |
 *    +--*'- - - - - - >                +--*--------+ - - - >
 *       0                                 \hat(xi)_0
 *                                                                  \endverbatim
 *  <p>
 *  \f$ z = k(xi - \hat(xi))^c e^(i\alpha)\f$</br/>
 *  where</br>
 *  \f$alpha\f$ = angle of ramp</br>
 *  \f$\hat(xi)_0\f$ = location of corner</br>
 *  \f$c = 1 - \alpha/\pi\f$</br>
 *  \f$k = len/(\hat(xi)_0^c + (1 - \hat(xi)_0)^c)\f$ where len is the length
 *  of the bottom boundary (both the horizontal and following the ramp) in
 *  physical space
 *
 *  \note
 *  <ul>
 *    <li> In computational space, the domain extends from 0 to 1 in the x
 *         direction and the corner is anchored at m_xiCorner.  The cells
 *         are square and the domain extends from 0 in the y direction (height
 *         is not specified)
 *    <li> In the physical space, the corner is anchored at (0, 0) and the
 *         length of the bottom boundary is specified.
 *    <li> If a third dimension is hanging around, it is scaled by k
 *    <li> Tolerances tested up to 1024 cells in X-Dir.
 *  </ul>
 *
 ******************************************************************************/

class LogSchwarzChristoffelRampCS : public NewFourthOrderCoordSys
{
public:
  /// Constructor
  LogSchwarzChristoffelRampCS(const RealVect& a_dx,
                              const int       a_nCellBaseDir0,  // See .cpp for
                              const Real      a_alpha = Pi/6.,  // details
                              const Real      a_XLeadLength = 0.5,
                              const Real      a_XRampLength = 4.5,
                              const RealVect& a_shift = RealVect::Unit,
                              const RealVect& a_domL  = RealVect::Unit);

  /// Destructor
  virtual ~LogSchwarzChristoffelRampCS();

  /// Given coordinate in mapped space, return its location in real space
  virtual RealVect realCoord(const RealVect& a_hatXi) const;

  /// Given coordinate in real space, return its location in the mapped space
  virtual RealVect mappedCoord(const RealVect& a_x) const;

  /// Calculate the derivative of each coordinate vector
  virtual Real dXdXi(const RealVect& a_Xi, int a_dirX, int a_dirXi) const;

  /// Return the angle of the ramp in radians
  Real rampAngle() const;

  /// Return the location of the corner point in mapped space
  const RealVect& mappedCorner() const;

  /// Return the x-location of lower left domain corner in physical space
  Real physicalDomainX0Min() const;

protected:

  /// Get a complex in polar coordinates
  RealVect getPolar(const RealVect& a_Xi) const;

  const Real m_alpha;               ///< Ramp angle (0 <= m_alpha < Pi)
  const Real m_c;                   ///< Exponent coefficient
  Real m_k;                         ///< Scaling between mapped and physical
                                    ///< space
  RealVect m_XiCorner;              ///< Corner in mapped space
  RealVect m_shift;                 ///< Vector of the stretching
  RealVect m_eShift;                ///< exp(m_shift) - 1
  RealVect m_diffShift;             ///< m_eShift/m_shift
  RealVect m_domL;                  ///< computational domain length
};


/*******************************************************************************
 */
/// Factory for LogSchwarzChristoffelRampCS
/*
 ******************************************************************************/

class LogSchwarzChristoffelRampCSFactory : public NewCoordSysFactory
{
public:
  /// Constructor
  LogSchwarzChristoffelRampCSFactory(
    const int       a_nCellBaseDir0, // See .cpp for details
    const Real      a_alpha = Pi/6.,
    const Real      a_XLeadLength = 0.5,
    const Real      a_XRampLength = 4.0,
    const RealVect& a_shift = RealVect::Unit);

  ///
  virtual NewCoordSys* getCoordSys(const ProblemDomain& a_domain,
                                   const RealVect& a_dx) const;

protected:

  const int  m_nCellBaseDir0;
  const Real m_alpha;
  const Real m_XLeadLength;
  const Real m_XRampLength;
  const RealVect m_shift;           ///< Vector of the stretching
};


/*******************************************************************************
 *
 * Class LogSchwarzChristoffelRampCS: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Return the angle of the ramp in radians
/*--------------------------------------------------------------------*/

inline Real
LogSchwarzChristoffelRampCS::rampAngle() const
{
  return m_alpha;
}

/*--------------------------------------------------------------------*/
//  Return the location of the corner point in mapped space
/*--------------------------------------------------------------------*/

inline const RealVect&
LogSchwarzChristoffelRampCS::mappedCorner() const
{
  return m_XiCorner;
}

/*--------------------------------------------------------------------*/
//  Return the x-location of lower left domain corner in physical
//  space
/*--------------------------------------------------------------------*/

inline Real
LogSchwarzChristoffelRampCS::physicalDomainX0Min() const
{
  return -m_k*std::pow(std::fabs(m_XiCorner[0]), m_c);
}

#include "NamespaceFooter.H"

#endif
