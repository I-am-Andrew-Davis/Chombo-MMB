#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SMOOTHRAMPBLOCKSCS_H_
#define _SMOOTHRAMPBLOCKSCS_H_

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
#include "NewFourthOrderCoordSys.H"
#include "Spline1DMapping.H"

#include "NamespaceHeader.H"
// #include "UsingNamespace.H"

/// Multi-block coordinate-system for smooth-ramp case 
/**
   The SmoothRampBlocks class implements the NewFourthOrderCoordSys class
   for the smooth-ramp case.
*/
class SmoothRampBlocksCS : public NewFourthOrderCoordSys
{
public:
  SmoothRampBlocksCS(const RealVect&      a_dX,
                     const ProblemDomain& a_domain,
                     const int            a_blkIdx,
                     const IntVect&       a_evalGhost,
                     IntVect&             a_ix,
                     const Real           a_stretch,
                     const Real           a_rampHeightRelax,
                     const Real           a_zWidth,
                     const Real           a_xLengthInletBlock,
                     const Real           a_xLengthRampBlockPreRamp,
                     const Real           a_xLengthRampBlock,
                     const Real           a_xLengthOutletBlock,
                     const int            a_numRampXCells);

  /**
     Destructor.
   */
  virtual ~SmoothRampBlocksCS();

  /// given coordinate in mapped space, return its location in real space
  RealVect realCoord(const RealVect& a_Xi) const;

  /// given coordinate in real space, return its location in the mapped space
  RealVect mappedCoord(const RealVect& a_x) const;

  // Calculate the derivative of each Cartesian coordinate vector
  virtual Real dXdXi(const RealVect& a_Xi, int a_dirX, int a_dirXi) const;

protected:

  ProblemDomain m_domain;             ///< Problem domain for this block
  IntVect       m_evalGhost;          ///< Num ghost cells to extrapolate
  int           m_blkIdx;             ///< Block index (panel ID)
  Real          m_stretch;            ///< Stretching factor for y-direction
  Real          m_rampHeightRelax;    ///< Relaxation factor for y-direction
                                      ///< for ramp-height
  Real          m_zWidth;             ///< Spanwise physical-space width
  Real          m_xLengthInletBlock;  ///< X-length of inlet block
  Real          m_xLengthRampBlockPreRamp;
                                      ///< X-length of ramp block before ramp
  Real          m_xLengthRampBlock;   ///< X-length of total ramp block
  Real          m_xLengthOutletBlock; ///< X-length of outlet block
  int           m_numRampXCells;      ///< Number of cells in x-direction
                                      ///< along ramp
  IntVect       m_ix;                 ///< Index of the bottom-left edge of
                                      ///< this panel in index space
  RealVect      m_xiOrigin;           ///< Start location of xi
};

#include "NamespaceFooter.H"

#endif
