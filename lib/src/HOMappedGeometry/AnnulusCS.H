#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _ANNULUSCS_H_
#define _ANNULUSCS_H_

#include "NewFourthOrderCoordSys.H"
#include "NamespaceHeader.H"

/*******************************************************************************
 */
///  A mapping for an annulus geometry. The grid is mapped it fill the space
///   between concentric cylinders.
/**
 * 
 * The annulus transformation is illustrated here.
 *
 *         (th)
 *      -----------                   ______          
 *      |         |              rh /        \       
 *      |         |                /    _     \   
 * (rl) |         | (rh)          |    / \  tl |
 *      |         |               |   |   |----|
 *      |         |               |  rl\_/  th |
 *      |         |                \          /      
 *      |         |                 \ ______ /           
 *      -----------                  
 *         (tl)
 *   computational space         physical space
 *
 * tl - theta at the low side
 * th - theta at high side
 * rl - inner radius
 * rh - outer radius
 ******************************************************************************/

class AnnulusCS : public NewFourthOrderCoordSys
{
public:
  /// Constructor
  AnnulusCS(const RealVect& a_dX,
            const Real&     a_inRadius,
            const Real&     a_outRadius,
            const RealVect& a_domainLen,
            const Real&     a_radialStretch = 1.0);

  /// Destructor
  virtual ~AnnulusCS();

  /// Given coordinate in mapped space, return its location in real space
  RealVect realCoord(const RealVect& a_Xi) const;

  /// Given coordinate in real space, return its location in mapped space
  RealVect mappedCoord(const RealVect& a_x) const;

  // Calculate the derivative of each coordinate vector
  virtual Real dXdXi(const RealVect& a_Xi, int a_dirX, int a_dirXi) const;

protected:
  const Real m_inRadius;              ///< radius of the inner cylinder
  const Real m_outRadius;             ///< radius of the outer cylinder
  const RealVect m_domainLen;         ///< length of the computational domain
  const Real m_Cs;                    ///< radial stretching factor

  // convenience
  const Real m_twoPi;
};

/*******************************************************************************
 */
///  Factory object to create an annulus coordinate system
/**
 ******************************************************************************/
class
AnnulusCSFactory : public NewCoordSysFactory
{

public:
  /// Constructor
  AnnulusCSFactory(const Real& a_inRadius,
                   const Real& a_outRadius,
                   const Real& m_radialStretch);

  /// Return a pointer to the AnnulusCS object
  virtual NewCoordSys* getCoordSys(const ProblemDomain& a_domain,
                                   const RealVect& a_dx) const;

protected:
  const Real m_inRadius;              ///< radius of the inner cylinder
  const Real m_outRadius;             ///< radius of the outer cylinder
  const Real m_radialStretch;         ///< stretch in the radial direction
};


#include "NamespaceFooter.H"

#endif
