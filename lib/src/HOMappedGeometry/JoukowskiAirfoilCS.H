#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _JOUKOWSKIAIRFOILCS_H_
#define _JOUKOWSKIAIRFOILCS_H_

//----- Standard Library -----//

#include <cmath>
#include <complex>

//----- Chombo -----//

#include "CONSTANTS.H"
#include "NewFourthOrderCoordSys.H"
#include "StretchingTanh.H"
#include "StretchingCosSmlCtr.H"

#include "NamespaceHeader.H"


/*******************************************************************************
 */
/// Joukowski mapping of an airfoil
/**                                                                    \verbatim
 *  \note
 *  <ul>
 *    <li> In computational space, the domain extends from 0 to 1 in all
 *         directions
 *    <li> If a third dimension is hanging around, it is scaled by the airfoil
 *         span
 *  </ul>
 *
 ******************************************************************************/

class JoukowskiAirfoilCS : public NewFourthOrderCoordSys
{
public:

  /// Constructor
  JoukowskiAirfoilCS(
    const RealVect& a_dx,
    const Real      a_chord            = 1.,
    const Real      a_span             = 1.,
    const Real      a_domainRatio      = 4.,
    const Real      a_camberRatio      = 0.05,
    const Real      a_thicknessRatio   = 0.15,
    const Real      a_alpha            = 0.,
    const Real      a_cellRatioRadial  = 0.1,
    const Real      a_cellRatioAzimuth = 0.2);

  /// Destructor
  virtual ~JoukowskiAirfoilCS();

  /// Given coordinate in mapped space, return its location in real space
  virtual RealVect realCoord(const RealVect& a_hatXi) const;

  /// Given coordinate in real space, return its location in the mapped space
  virtual RealVect mappedCoord(const RealVect& a_x) const;

  /// Calculate the derivative of each coordinate vector in a Jacobian
  virtual void dXdXiJacobian(const RealVect& a_Xi,
                             Real*           a_J,
                             int             a_ldJ) const;

  /// Calculate the derivative of each coordinate vector
  virtual Real dXdXi(const RealVect& a_Xi, int a_dirX, int a_dirXi) const;

protected:

  const Real m_chord;                 ///< Airfoil chord
  const Real m_span;                  ///< Airfoil span
  const Real m_camberRatio;           ///< Camber ratio (h/chord)
  const Real m_thicknessRatio;        ///< Thickness ratio (t/chord)
  const Real m_alpha;                 ///< Angle of attack
  const Real m_c;                     ///< Critical point (for potential flow
                                      ///< around a Joukowski airfoil)
  const Real m_kappa;                 ///< Shift due to camber
  const Real m_epsilon;               ///< Shift due to thickness
  const Real m_r0;                    ///< Radius of circle mapped to airfoil
  const Real m_r1;                    ///< Radius of domain before mapping to
                                      ///< airfoil
  const std::complex<Real> m_mu;      ///< Origin of circle mapped to airfoil
  const StretchingCosSmlCtr<2> m_stretchCosX;
                                      ///< Stretching in azimuthal direction
  const StretchingTanh m_stretchTanhY;
                                      ///< Stretching in radial direction
  RealVect m_XTE;                     ///< Location of trailing edge
  std::vector<RealVect> m_mappedGuessXi;
  std::vector<RealVect> m_mappedGuessX;

public:

  mutable RealVect m_X;               ///< For evaluating residuals
};


/*******************************************************************************
 */
/// Factory for JoukowskiAirfoilCS
/*
 ******************************************************************************/

class JoukowskiAirfoilCSFactory : public NewCoordSysFactory
{
public:
  /// Constructor
  JoukowskiAirfoilCSFactory(
    const Real a_chord            = 1.,
    const Real a_span             = 1.,
    const Real a_domainRatio      = 4.,
    const Real a_camberRatio      = 0.05,
    const Real a_thicknessRatio   = 0.15,
    const Real a_alpha            = 0.,
    const Real a_cellRatioRadial  = 0.1,
    const Real a_cellRatioAzimuth = 0.2);

  ///
  virtual NewCoordSys* getCoordSys(const ProblemDomain& a_domain,
                                   const RealVect&      a_dx) const;

protected:

  const Real m_chord;
  const Real m_span;
  const Real m_domainRatio;
  const Real m_camberRatio;
  const Real m_thicknessRatio;
  const Real m_alpha;
  const Real m_cellRatioRadial;
  const Real m_cellRatioAzimuth;
};


// /*******************************************************************************
//  *
//  * Class JoukowskiAirfoilCS: inline member definitions
//  *
//  ******************************************************************************/

// /*--------------------------------------------------------------------*/
// //  Return the angle of the ramp in radians
// /*--------------------------------------------------------------------*/

// inline Real
// JoukowskiAirfoilCS::rampAngle() const
// {
//   return m_alpha;
// }

// /*--------------------------------------------------------------------*/
// //  Return the location of the corner point in mapped space
// /*--------------------------------------------------------------------*/

// inline const RealVect&
// JoukowskiAirfoilCS::mappedCorner() const
// {
//   return m_XiCorner;
// }

// /*--------------------------------------------------------------------*/
// //  Return the x-location of lower left domain corner in physical
// //  space
// /*--------------------------------------------------------------------*/

// inline Real
// JoukowskiAirfoilCS::physicalDomainX0Min() const
// {
//   return -m_k*std::pow(std::fabs(m_XiCorner[0]), m_c);
// }

#include "NamespaceFooter.H"

#endif
