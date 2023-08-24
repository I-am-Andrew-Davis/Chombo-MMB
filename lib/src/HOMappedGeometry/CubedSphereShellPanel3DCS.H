#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CUBEDSPHERESHELLPANEL3DCS_H_
#define _CUBEDSPHERESHELLPANEL3DCS_H_

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

/// Basic single-panel cubed sphere Coordinates
/// Multiblock coordinate system for single-panel cubed sphere shell
/// with variable-thickness spherical layers 
/**
   The CubedSphereShellPanelCS class implements the NewFourthOrderCoordSys class
   for cubed sphere panels in equiangular coordinates.
*/
class CubedSphereShellPanel3DCS : public NewFourthOrderCoordSys
{
public:
  /// constructor
  /**
     nPanel
       (0-3 Equatorial, 4 North Pole, 5 South Pole) and defines the
       panel on the cubed sphere described by this object.
     dDX
       The cell spacing on this panel.
     ix
       Index of the bottom-left edge of this panel in index space.
  */
  CubedSphereShellPanel3DCS(const RealVect&      a_dX,
                            const ProblemDomain& a_domain,
                            const int            a_blkIdx,
                            const IntVect&       a_evalGhost,
                            IntVect&             a_ix,
                            const Real           a_radiusInner,
                            const Real           a_radiusOuter,
                            const Real           a_alpha,
                            const Real           a_beta);

  /**
     Destructor.
   */
  virtual ~CubedSphereShellPanel3DCS();

  //void panel(int nPanel) {m_nPanel = nPanel;}

  /// given coordinate in mapped space, return its location in real space
  RealVect realCoord(const RealVect& a_Xi) const;

  /// given coordinate in real space, return its location in the mapped space
  RealVect mappedCoord(const RealVect& a_x) const;

  // Calculate the derivative of each Cartesian coordinate vector
  virtual Real dXdXi(const RealVect& a_Xi, int a_dirX, int a_dirXi) const;

/*--------------------------------------------------------------------*/
//  Non-linear function for radial mapping for cubed sphere
/**
 *  \param[in]  a_r     Physical-space radial-distance from center
 *  \param[in]  a_radiusInner
 *                      Inner radius of sphere
 *  \param[in]  a_radiusOuter
 *                      Outer radius of sphere
 *  \param[in]  a_alpha Stretching constant
 *  \param[in]  a_beta  Stretching constant
 *//*-----------------------------------------------------------------*/

  struct RadiusFunc
  {
    RadiusFunc(const Real a_r,
               const Real a_radiusInner,
               const Real a_radiusOuter,
               const Real a_alpha,
               const Real a_beta
      )
      :
      m_r(a_r),
      m_radiusInner(a_radiusInner),
      m_radiusOuter(a_radiusOuter),
      m_alpha(a_alpha),
      m_beta(a_beta)
      {}
    Real operator()(const Real& a_etaHat) const
      {
        Real alpha_eta = (m_radiusOuter - m_radiusInner)
          - m_alpha*std::tanh(m_beta);
        return m_radiusInner + alpha_eta*a_etaHat
          + m_alpha*std::tanh(m_beta*a_etaHat) - m_r;
      }
    const Real m_r;
    const Real m_radiusInner;
    const Real m_radiusOuter;
    const Real m_alpha;
    const Real m_beta;
  };

protected:

  ProblemDomain m_domain;             ///< Problem domain for this block
  IntVect       m_evalGhost;          ///< Num ghost cells to extrapolate
  int           m_blkIdx;             ///< Block index (panel ID)
  Real          m_radiusInner;        ///< Inner radius of the physical-space
                                      ///< spherical shell mapping
  Real          m_radiusOuter;        ///< Outer radius of the physical-space
                                      ///< spherical shell mapping
  Real          m_alpha;              ///< Radial stretching coefficient
  Real          m_beta;               ///< Radial stretching coefficient
  IntVect       m_ix;                 ///< Index of the bottom-left edge of
                                      ///< this panel in index space
  RealVect      m_xiOrigin;           ///< Start location of xi
};

#include "NamespaceFooter.H"

#endif
