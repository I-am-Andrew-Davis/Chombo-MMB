#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LOGCS_H_
#define _LOGCS_H_

#include "NewFourthOrderCoordSys.H"
#include "NamespaceHeader.H"

/*******************************************************************************
 */
///  A log stretched grid mapping, that maintains the boundary box
///  Stretching of 0 in any direction results in a Cartesian grid.
/**
 ******************************************************************************/

class LogCS : public NewFourthOrderCoordSys
{
public:
  /// Constructor
  LogCS(const RealVect& a_dX,
        const RealVect& a_shift,
        const RealVect a_domL);

  /// Destructor
  virtual ~LogCS();

  /// Given coordinate in mapped space, return its location in real space
  RealVect realCoord(const RealVect& a_Xi) const;

  /// Given coordinate in real space, return its location in mapped space
  RealVect mappedCoord(const RealVect& a_x) const;

  // Calculate the derivative of each coordinate vector
  virtual Real dXdXi(const RealVect& a_Xi, int a_dirX, int a_dirXi) const;

protected:
  // scale for displacement
  RealVect m_shift;                   ///< Vector of the stretching
  RealVect m_eShift;                  ///< exp(m_shift) - 1
  RealVect m_diffShift;               ///< m_eShift/m_shift
  RealVect m_domL;                    ///< computational domain length
};

/*******************************************************************************
 */
///  Factory object to create an log coordinate system
/**
 ******************************************************************************/
class
LogCSFactory : public NewCoordSysFactory
{

public:
  /// Constructor
  LogCSFactory(const RealVect& a_shift);

  /// Return a pointer to the LogCS object
  virtual NewCoordSys* getCoordSys(const ProblemDomain& a_domain,
                                   const RealVect&      a_dx) const;

protected:
  RealVect m_shift;
};


#include "NamespaceFooter.H"
#endif
