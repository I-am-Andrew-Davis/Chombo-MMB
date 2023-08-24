#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CUBEDSPHERESHELL3DCS_H_
#define _CUBEDSPHERESHELL3DCS_H_

#include <utility>
using std::pair;

#include "REAL.H"
#include "DisjointBoxLayout.H"
#include "RealVect.H"
#include "IntVect.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "Vector.H"
#include "Spline1DMapping.H"
#include "MultiBlockCoordSys.H"

#include "NamespaceHeader.H"
// #include "UsingNamespace.H"

/// Multiblock coordinate system for cubed sphere shell
/// with variable-thickness spherical layers 
/**
   The CubedSphereShell3DCS coordinate system class is a multiblock example
   that takes place on a set of cubes wrapping around a sphere, and
   creates a shell with r-direction layers of variable thickness.
*/
class CubedSphereShell3DCS : public MultiBlockCoordSys
{
public:
  /// default constructor
  CubedSphereShell3DCS();

  /**
     Destructor.
   */
  virtual ~CubedSphereShell3DCS();

  virtual void define(const ProblemDomain& a_levelDomain,
                      const RealVect&      a_dx,
                      const Real           a_radiusInner,
                      const Real           a_radiusOuter,
                      const Real           a_alpha,
                      const Real           a_beta);

  /// block mapping conversion function
  /** given a_xiSrc in block a_nSrc, returns which block (a_n_valid)
      contains that location in its valid region and its
      location (a_xi_valid) in that block's computational coordinate space.
  */
  virtual void blockRemapping(RealVect&            a_xi_valid,
                              int&                 a_n_valid,
                              bool&                a_validExists,
                              RigidTransformation& a_extraDispl,
                              const RealVect&      a_xiSrc,
                              const IntVect&       a_iSrc,
                              int                  a_nSrc) const;

protected:

  int  m_numBlocks;                   ///< Number of blocks -- should be 6

  // define m_boundaries, called by define()
  void defineBoundaries();
};


class CubedSphereShell3DCSFactory : public MultiBlockCoordSysFactory
{
public:

  CubedSphereShell3DCSFactory(const Real a_radiusInner,
                              const Real a_radiusOuter,
                              const Real a_alpha,
                              const Real a_beta)
    :
    m_radiusInner(a_radiusInner),
    m_radiusOuter(a_radiusOuter),
    m_alpha(a_alpha),
    m_beta(a_beta)
    { }

  virtual ~CubedSphereShell3DCSFactory()
  {
  }

  virtual MultiBlockCoordSys* getCoordSys(const ProblemDomain& levelDomain,
                                          const RealVect& a_dx) const;

protected:

  Real m_radiusInner;
  Real m_radiusOuter;
  Real m_alpha;
  Real m_beta;
};

#include "NamespaceFooter.H"

#endif
