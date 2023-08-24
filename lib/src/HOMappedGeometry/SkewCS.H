#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SKEWCS_H_
#define _SKEWCS_H_

#include "NewFourthOrderCoordSys.H"
#include "NamespaceHeader.H"

/******************************************************************************/
/**
 * A parallelogram coordinate system. This has a constant Jacobian and grid
 *  metrics. The grid metrics are a triangular matrix.
 *
 *          -----------------   -
 *         /               /    |
 *        /  angle        /     | height
 *       / )             /      |
 *       ----------------       -
 *        length
 *//***************************************************************************/

class SkewCS : public NewFourthOrderCoordSys
{
public:
  /// Constructor
  SkewCS(const RealVect& a_dX,
         const Real&     a_angle,
         const RealVect& a_physLength,
         const RealVect& a_domainLength);
  
  /// Destructor
  virtual ~SkewCS();

  /// Given coordinate in mapped space, return its location in real space
  RealVect realCoord(const RealVect& a_Xi) const;

  /// Given coordinate in real space, return its location in mapped space
  RealVect mappedCoord(const RealVect& a_x) const;

  // Calculate the derivative of each coordinate vector
  virtual Real dXdXi(const RealVect& a_Xi, int a_dirX, int a_dirXi) const;

protected:
  Real m_angle;                       ///< angle of skewness
  RealVect m_physLeng;                ///< height and length to transform to
  RealVect m_domainLeng;              ///< size of computational domain
};

/*******************************************************************************
 */
///  Factory object to create an skew coordinate system
/**
 ******************************************************************************/
class
SkewCSFactory : public NewCoordSysFactory
{

public:
  /// Constructor
  SkewCSFactory(const Real&     a_angle,
                const RealVect& a_physicalLength);

  /// Return a pointer to the SkewCS object
  virtual NewCoordSys* getCoordSys(const ProblemDomain& a_domain,
                                   const RealVect&      a_dx) const;

protected:
  Real m_angle;                       ///< angle of skewness
  RealVect m_physLeng;                ///< height and length to transform to
};


#include "NamespaceFooter.H"

#endif
