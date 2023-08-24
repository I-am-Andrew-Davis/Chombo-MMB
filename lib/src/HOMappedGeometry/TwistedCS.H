#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _TWISTEDCS_H_
#define _TWISTEDCS_H_

#include "NewFourthOrderCoordSys.H"

#include "NamespaceHeader.H"

/*******************************************************************************
 */
/// Calhoun-Leveque "twisted" Coordinate mapping (constant Jacobian)
/**
 *//*+*************************************************************************/

class TwistedCS : public NewFourthOrderCoordSys
{
public:
  /// Constructor
  TwistedCS(const RealVect& a_dX,
            const Real& a_R,
            const Real& a_twist,
            const RealVect& a_domainLength);

  /// Destructor
  virtual ~TwistedCS();

  /// Given coordinate in mapped space, return its location in real space
  RealVect realCoord(const RealVect& a_Xi) const;

  /// Given coordinate in real space, return its location in mapped space
  RealVect mappedCoord(const RealVect& a_x) const;

  // Calculate the derivative of each coordinate vector
  virtual Real dXdXi(const RealVect& a_Xi, int a_dirX, int a_dirXi) const;


protected:
  RealVect m_domainLength;            ///< physical size of the domain
  Real m_R;                           ///< radius for n-sphere of twisted grid
  Real m_theta;                       ///< Amount of twist in the grid

  // convenience
  Real m_Pi;
  Real m_scale;

};

/*******************************************************************************
 */
///  Factory object to create a twisted coordinate system
/**
 *//*+*************************************************************************/
class
TwistedCSFactory : public NewCoordSysFactory
{

public:
  /// Constructor
  TwistedCSFactory(const Real& a_R,
                   const Real& a_twist);

  /// Return a pointer to the TwistedCS object
  virtual NewCoordSys* getCoordSys(const ProblemDomain& a_domain,
                                   const RealVect& a_dx) const;

protected:
  Real m_R;                           ///< radius for n-sphere of twisted grid
  Real m_twist;                       ///< Amount of twist in the grid
};

#include "NamespaceFooter.H"

#endif
