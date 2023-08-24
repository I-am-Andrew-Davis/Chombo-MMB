#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SMOOTHRAMPCS_H_
#define _SMOOTHRAMPCS_H_

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

/// Multiblock coordinate system for smooth-ramp case
/// with variable, linear mesh stretching
/**
   The SmoothRamp coordinate system class is a multiblock example
   that takes place on 3 blocks: an inlet block, a ramp block, and an outlet
   block. The case is uniform in the z-direction. The inlet and outlet blocks
   are uniform in the x-direction. User-defined linear-stretching is allowed
   in the y-direction (wall-normal direction). This case assumes that the
   top-boundary is a slip-wall. Ghost cells above this wall may have some
   strange behavior, but they should still be well defined in general. The
   ghost cells below the ramp should be a little more well behaved as long as
   the mesh stretching isn't too much.
*/
class SmoothRampCS : public MultiBlockCoordSys
{
public:
  /// default constructor
  SmoothRampCS();

  /**
     Destructor.
   */
  virtual ~SmoothRampCS();

  virtual void define(const ProblemDomain& a_levelDomain,
                      const RealVect&      a_dx,
                      const Real           a_stretch,
                      const Real           a_rampHeightRelax,
                      const Real           a_zWidth,
                      const int            a_onlyInletBlock,
                      const int            a_numXCellsInletBlock,
                      const int            a_numXCellsRampBlockPreRamp,
                      const int            a_numXCellsRampBlockRamp,
                      const int            a_numXCellsRampBlockPostRamp,
                      const int            a_numXCellsOutletBlock);

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

  Real m_zWidth;
  int  m_numBlocks;                   ///< Number of blocks -- could be 1 to 3

  // define m_boundaries, called by define()
  void defineBoundaries();
};


class SmoothRampCSFactory : public MultiBlockCoordSysFactory
{
public:

  SmoothRampCSFactory(const Real a_stretch,
                      const Real a_rampHeightRelax,
                      const Real a_zWidth,
                      const int  a_onlyInletBlock,
                      const int  a_numXCellsInletBlock,
                      const int  a_numXCellsRampBlockPreRamp,
                      const int  a_numXCellsRampBlockRamp,
                      const int  a_numXCellsRampBlockPostRamp,
                      const int  a_numXCellsOutletBlock)
    :
    m_stretch(a_stretch),
    m_rampHeightRelax(a_rampHeightRelax),
    m_zWidth(a_zWidth),
    m_onlyInletBlock(a_onlyInletBlock),
    m_numXCellsInletBlock(a_numXCellsInletBlock),
    m_numXCellsRampBlockPreRamp(a_numXCellsRampBlockPreRamp),
    m_numXCellsRampBlockRamp(a_numXCellsRampBlockRamp),
    m_numXCellsRampBlockPostRamp(a_numXCellsRampBlockPostRamp),
    m_numXCellsOutletBlock(a_numXCellsOutletBlock)
    { }

  virtual ~SmoothRampCSFactory()
  {
  }

  virtual MultiBlockCoordSys* getCoordSys(const ProblemDomain& levelDomain,
                                          const RealVect& a_dx) const;

protected:

  Real m_stretch;
  Real m_rampHeightRelax;
  Real m_zWidth;
  int  m_onlyInletBlock;
  int  m_numXCellsInletBlock;
  int  m_numXCellsRampBlockPreRamp;
  int  m_numXCellsRampBlockRamp;
  int  m_numXCellsRampBlockPostRamp;
  int  m_numXCellsOutletBlock;
};

#include "NamespaceFooter.H"

#endif
