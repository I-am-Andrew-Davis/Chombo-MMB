#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _WARPEDCS_H_
#define _WARPEDCS_H_

#include "NewFourthOrderCoordSys.H"
#include "NamespaceHeader.H"

#ifdef CH_USE_FLOAT
#  define RTOL 1.e-6
#  define ATOL 1.e-7
#else
#  define RTOL 1.e-12
#  define ATOL 1.e-13
#endif

/*******************************************************************************
 */
/// A "warped" coordinate mapping. This has a non-constant Jacobian and grid
/// metrics, but does preserve a rectangular domain.
/// This mapping is primarily useful for testing
/**
 *//*+*************************************************************************/
class WarpedCS : public NewFourthOrderCoordSys
{
public:
  /// Constructor
  WarpedCS(const RealVect& a_dX,
           const RealVect& a_scale,
           const RealVect a_domainLength,
           const Real a_rtol=RTOL,
           const Real a_atol=ATOL,
           const int a_imax=100);

  /// Destructor
  virtual ~WarpedCS();

  /// Given coordinate in mapped space, return its location in real space
  RealVect realCoord(const RealVect& a_Xi) const;

  /// Given coordinate in real space, return its location in mapped space
  RealVect mappedCoord(const RealVect& a_x) const;

  // Calculate the derivative of each coordinate vector
  virtual Real dXdXi(const RealVect& a_Xi, int a_dirX, int a_dirXi) const;

  /// Sets the time for moving grids and adjusts m_scale
  virtual void setTime(const Real& a_time) const;

  /// Initialize scale
  void scale(const RealVect& a_scale);

protected:
  // scale for displacement
  mutable RealVect m_scale;           ///< The 'c'
  const RealVect m_staticScale;       ///< Preserves initial scale argument
                                      ///< given to constructor for moving grids

  // parameters for fixed-point iteration for inverse mapping
  Real m_rtol;
  Real m_atol;
  int m_imax;

  // convenience
  const Real m_twoPi;
  mutable int m_rootDir;
  RealVect m_domainLength;

  // RHS of fixed-point iteration for inverse mapping
  Real g( const Real xi, const RealVect& x ) const;
};

/*******************************************************************************
 */
///  Factory object to create an warped coordinate system
/**
 *//*+*************************************************************************/
class
WarpedCSFactory : public NewCoordSysFactory
{

public:
  /// Constructor
  WarpedCSFactory(const RealVect& a_scale,
                  const Real a_rtol=RTOL,
                  const Real a_atol=ATOL,
                  const int a_imax=100);

  /// Return a pointer to the WarpedCS object
  virtual NewCoordSys* getCoordSys(const ProblemDomain& a_domain,
                                   const RealVect& a_dx) const;

protected:
  RealVect m_scale;                   ///< The scale for displacement
  Real m_rtol;                        ///< Tolerance 
  Real m_atol;
  int m_imax;
};


#include "NamespaceFooter.H"

#endif
