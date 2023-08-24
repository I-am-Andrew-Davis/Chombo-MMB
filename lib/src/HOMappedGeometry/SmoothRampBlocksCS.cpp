#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "SmoothRampBlocksCS.H"
#include "BoxIterator.H"
#include <cmath>
#include "MBMiscUtil.H"
#include "csignal"
#include "RootSolver.H"

#include "FourthOrderUtil.H"
#include "FourthOrderUtilF_F.H"
#include "FCDivergenceF_F.H"
#include "AdvectOpF_F.H"

#include "NamespaceHeader.H"
// #include "UsingNamespace.H"

SmoothRampBlocksCS::SmoothRampBlocksCS(
  const RealVect&      a_dX,
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
  const int            a_numRampXCells)
  :
  m_domain(a_domain),
  m_evalGhost(a_evalGhost),
  m_blkIdx(a_blkIdx),
  m_stretch(a_stretch),
  m_rampHeightRelax(a_rampHeightRelax),
  m_zWidth(a_zWidth),
  m_xLengthInletBlock(a_xLengthInletBlock),
  m_xLengthRampBlockPreRamp(a_xLengthRampBlockPreRamp),
  m_xLengthRampBlock(a_xLengthRampBlock),
  m_xLengthOutletBlock(a_xLengthOutletBlock),
  m_numRampXCells(a_numRampXCells),
  m_ix(a_ix)
{
  CH_assert((a_blkIdx >= 0) && (a_blkIdx < 3));

  m_dx = a_dX;
  // This is for AMR
  m_xiOrigin = RealVect{m_dx*m_domain.domainBox().smallEnd()};
}

SmoothRampBlocksCS::~SmoothRampBlocksCS()
{
}

RealVect
SmoothRampBlocksCS::realCoord(const RealVect& a_Xi) const
{
  CH_TIME("SmoothRampBlocksCS::realCoord");

#if CH_SPACEDIM == 3
  // Getting the z-direction coordinate is easy
  const Real zMin = 0.;
  const Real zMax = m_zWidth;
  const Real zetaMin = 0.;
  const Real zetaMax = (m_dx[2])*(m_domain.domainBox().size()[2]);
  const Real zVal = a_Xi[2]*((zMax - zMin)/(zetaMax - zetaMin))
    + ((zMin*zetaMax - zMax*zetaMin)/(zetaMax - zetaMin));
#endif

  // Getting x-direction coordinate is relatively easy --- the ramp length
  // is always set to unit length and physDx is calculated based on that.
  //**NOTE: The beginning of the ramp is set to the zero location
  const Real refLoc = 0.; // Located at start of ramp
  Real xMin = 0.;
  Real xMax = 1.;
  if (m_blkIdx == 0) // Inlet block
    {
      xMin = refLoc - m_xLengthRampBlockPreRamp - m_xLengthInletBlock;
      xMax = refLoc - m_xLengthRampBlockPreRamp;
    }
  else if (m_blkIdx == 1) // Ramp block
    {
      xMin = refLoc - m_xLengthRampBlockPreRamp;
      xMax = refLoc - m_xLengthRampBlockPreRamp + m_xLengthRampBlock;
    }
  else if (m_blkIdx == 2) // Outlet block
    {
      xMin = refLoc - m_xLengthRampBlockPreRamp + m_xLengthRampBlock;
      xMax = refLoc - m_xLengthRampBlockPreRamp + m_xLengthRampBlock
        + m_xLengthOutletBlock;
    }
  const Real xiMin = m_dx[0]*(m_domain.domainBox().smallEnd()[0]);
  const Real xiMax = m_dx[0]*(m_domain.domainBox().bigEnd()[0] + 1);
  const Real xVal = a_Xi[0]*((xMax - xMin)/(xiMax - xiMin))
    + ((xMin*xiMax - xMax*xiMin)/(xiMax - xiMin));

  // Getting y-direction coordinate is a little harder --- start with the
  // linear version just for the sake of debugging
  const Real Hval = 0.22; // Fixed value
  const Real yMin = 0.;   // Fixed value
  const Real yMax = 0.74; // Fixed value
  const Real etaMin = 0.;
  const Real etaMax = (m_dx[1])*(m_domain.domainBox().size()[1]);
  const Real hValScaleMin = 1.; // At the wall, Hval needs to be 0.22
  const Real hValScaleMax = 0.; // At the top boundary, Hval needs to be 0
  // All these terms are just being used here to collapse variables together
  // so that the resulting equations are simplified. There's really no rhyme
  // or reason to the naming.
  const Real hValStretch = m_rampHeightRelax;

  // These terms are for extruding the 2D plane in the wall-normal direction
  const Real A = m_stretch;
  const Real B = A*((etaMin*etaMin - etaMax*etaMax)/(etaMax - etaMin))
    + (yMax - yMin)/(etaMax - etaMin);
  const Real C = A*etaMin*etaMax
    + ((yMin*etaMax - yMax*etaMin)/(etaMax - etaMin));

  // This is the ramp-profile
  Real D = 1. - 10.*xVal*xVal*xVal + 15.*xVal*xVal*xVal*xVal
    - 6.*xVal*xVal*xVal*xVal*xVal;
  // Before and after the ramp, the wall is flat (i.e. D is constant)
  if (xVal < 0.) { D = 1.; }
  else if (xVal > 1.) { D = 0.; }

  // These terms are for relaxing the ramp-profile to a straight line
  const Real a = hValStretch;
  const Real b = a*((etaMin*etaMin - etaMax*etaMax)/(etaMax - etaMin))
    + (hValScaleMax - hValScaleMin)/(etaMax - etaMin);
  const Real c = a*etaMin*etaMax
    + ((hValScaleMin*etaMax - hValScaleMax*etaMin)/(etaMax - etaMin));

  const Real E = Hval*D*a + A;
  const Real F = Hval*D*b + B;
  const Real G = Hval*D*c + C;
  const Real yVal = E*a_Xi[1]*a_Xi[1] + F*a_Xi[1] + G;

  RealVect xVect(D_DECL(xVal, yVal, zVal));

  return xVect;
}

// given coordinate in real space, return its location in the mapped space
RealVect
SmoothRampBlocksCS::mappedCoord(const RealVect& a_x) const
{
  CH_TIME("SmoothRampBlocksCS::mappedCoord");

#if CH_SPACEDIM == 3
  // Getting the zeta-direction coordinate should be fairly easy
  const Real zMin = 0.;
  const Real zMax = m_zWidth;
  const Real zetaMin = 0.;
  const Real zetaMax = (m_dx[2])*(m_domain.domainBox().size()[2]);
  const Real zetaVal =
    (a_x[2] - ((zMin*zetaMax - zMax*zetaMin)/(zetaMax - zetaMin)))*(
      (zetaMax - zetaMin)/(zMax - zMin));
#endif

  // Getting x-direction coordinate shouldn't be too hard --- the ramp length
  // is always set to unit length and physDx is calculated based on that.
  //**NOTE: The beginning of the ramp is set to the zero location
  const Real refLoc = 0.; // Located at start of ramp
  Real xMin = 0.;
  Real xMax = 1.;
  const Real xVal = a_x[0];
  if (m_blkIdx == 0) // Inlet block
    {
      xMin = refLoc - m_xLengthRampBlockPreRamp - m_xLengthInletBlock;
      xMax = refLoc - m_xLengthRampBlockPreRamp;
    }
  else if (m_blkIdx == 1) // Ramp block
    {
      xMin = refLoc - m_xLengthRampBlockPreRamp;
      xMax = refLoc - m_xLengthRampBlockPreRamp + m_xLengthRampBlock;
    }
  else if (m_blkIdx == 2) // Outlet block
    {
      xMin = refLoc - m_xLengthRampBlockPreRamp + m_xLengthRampBlock;
      xMax = refLoc - m_xLengthRampBlockPreRamp + m_xLengthRampBlock
        + m_xLengthOutletBlock;
    }
  const Real xiMin = m_dx[0]*(m_domain.domainBox().smallEnd()[0]);
  const Real xiMax = m_dx[0]*(m_domain.domainBox().bigEnd()[0] + 1);
  const Real xiVal = (a_x[0] - ((xMin*xiMax - xMax*xiMin)/(xiMax - xiMin)))*(
    (xiMax - xiMin)/(xMax - xMin));

  // Getting y-direction coordinate is a little harder --- start with the
  // linear version just for the sake of debugging
  const Real Hval = 0.22; // Fixed value
  const Real yMin = 0.;   // Fixed value
  const Real yMax = 0.74; // Fixed value
  const Real etaMin = 0.;
  const Real etaMax = (m_dx[1])*(m_domain.domainBox().size()[1]);
  const Real hValScaleMin = 1.; // At the wall, Hval needs to be 0.22
  const Real hValScaleMax = 0.; // At the top boundary, Hval needs to be 0
  // All these terms are just being used here to collapse variables together
  // so that the resulting equations are simplified. There's really no rhyme
  // or reason to the naming.
  const Real hValStretch = m_rampHeightRelax;

  // These terms are for extruding the 2D plane in the wall-normal direction
  const Real A = m_stretch;
  const Real B = A*((etaMax*etaMax-etaMin*etaMin)/(etaMin - etaMax))
    + (yMin-yMax)/(etaMin-etaMax);
  const Real C = A*etaMin*etaMax
    + ((yMax*etaMin - yMin*etaMax)/(etaMin-etaMax));

  // This is the ramp-profile
  Real D = 1. - 10.*xVal*xVal*xVal + 15.*xVal*xVal*xVal*xVal
    - 6.*xVal*xVal*xVal*xVal*xVal;
  // Before and after the ramp, the wall is flat (i.e. D is constant)
  if (xVal < 0.) { D = 1.; }
  else if (xVal > 1.) { D = 0.; }

  // These terms are for relaxing the ramp-profile to a straight line
  const Real a = hValStretch;
  const Real b = a*((etaMin*etaMin - etaMax*etaMax)/(etaMax - etaMin))
    + (hValScaleMax - hValScaleMin)/(etaMax - etaMin);
  const Real c = a*etaMin*etaMax
    + ((hValScaleMin*etaMax - hValScaleMax*etaMin)/(etaMax - etaMin));

  const Real E = Hval*D*a + A;
  const Real F = Hval*D*b + B;
  const Real G = Hval*D*c + C;
  const Real I = G - a_x[1];

  Real etaVal = (a_x[1] - G)/F;
  if ((m_stretch || m_rampHeightRelax) > 0.)
    {
      etaVal = (-F + std::sqrt(F*F - 4.*E*I))/(2.*E);
    }

  RealVect xiVect(D_DECL(xiVal, etaVal, zetaVal));

  return xiVect;
}

Real
SmoothRampBlocksCS::dXdXi(const RealVect& a_Xi,
                        int             a_dirX,
                        int             a_dirXi) const
{
  CH_TIME("SmoothRampBlocksCS::dXdXi");
  CH_assert((m_blkIdx >= 0) && (m_blkIdx <= 2));

  // There are some relatively easy derivatives
  // dx/dxi is a constant --- (xMax - xMin)/(xiMax - xiMin)
  // dx/deta = 0
  // dx/dzeta = 0
  // dy/dzeta = 0
  // dz/dxi = 0
  // dz/deta = 0
  // dz/dzeta = constant --- (zMax - zMin)/(zetaMax - zetaMin)

  // There are two terms that are a little more work
  // dy/deta = G
  // dy/dxi = (dy/dx)*(dx/dxi)
  // dy/dx = -30*M*x*x + 60*M*x*x*x - 30*M*x*x*x*x
  // where M = (Hval*(etaMax - eta)/A)
  // This is just in the ramp region, everywhere else, it's zero

#if CH_SPACEDIM == 3
  // Getting the z-direction coordinate is easy
  const Real zMin = 0.;
  const Real zMax = m_zWidth;
  const Real zetaMin = 0.;
  const Real zetaMax = (m_dx[2])*(m_domain.domainBox().size()[2]);
  const Real zSlope = ((zMax - zMin)/(zetaMax - zetaMin));
#endif

  // Getting x-direction coordinate is relatively easy --- the ramp length
  // is always set to unit length and physDx is calculated based on that.
  //**NOTE: The beginning of the ramp is set to the zero location
  const Real refLoc = 0.; // Located at start of ramp
  Real xMin = 0.;
  Real xMax = 1.;
  if (m_blkIdx == 0) // Inlet block
    {
      xMin = refLoc - m_xLengthRampBlockPreRamp - m_xLengthInletBlock;
      xMax = refLoc - m_xLengthRampBlockPreRamp;
    }
  else if (m_blkIdx == 1) // Ramp block
    {
      xMin = refLoc - m_xLengthRampBlockPreRamp;
      xMax = refLoc - m_xLengthRampBlockPreRamp + m_xLengthRampBlock;
    }
  else if (m_blkIdx == 2) // Outlet block
    {
      xMin = refLoc - m_xLengthRampBlockPreRamp + m_xLengthRampBlock;
      xMax = refLoc - m_xLengthRampBlockPreRamp + m_xLengthRampBlock
        + m_xLengthOutletBlock;
    }
  const Real xiMin = m_dx[0]*(m_domain.domainBox().smallEnd()[0]);
  const Real xiMax = m_dx[0]*(m_domain.domainBox().bigEnd()[0] + 1);
  const Real xVal = a_Xi[0]*((xMax - xMin)/(xiMax - xiMin))
    + ((xMin*xiMax - xMax*xiMin)/(xiMax - xiMin));
  const Real xSlope = ((xMax - xMin)/(xiMax - xiMin));

  // Getting y-direction coordinate is a little harder --- start with the
  // linear version just for the sake of debugging
  const Real Hval = 0.22; // Fixed value
  const Real yMin = 0.;   // Fixed value
  const Real yMax = 0.74; // Fixed value
  const Real etaMin = 0.;
  const Real etaMax = (m_dx[1])*(m_domain.domainBox().size()[1]);
  const Real hValScaleMin = 1.; // At the wall, Hval needs to be 0.22
  const Real hValScaleMax = 0.; // At the top boundary, Hval needs to be 0
  // All these terms are just being used here to collapse variables together
  // so that the resulting equations are simplified. There's really no rhyme
  // or reason to the naming.
  const Real hValStretch = m_rampHeightRelax;

  // These terms are for extruding the 2D plane in the wall-normal direction
  const Real A = m_stretch;
  const Real B = A*((etaMax*etaMax-etaMin*etaMin)/(etaMin - etaMax))
    + (yMin-yMax)/(etaMin-etaMax);

  // This is the ramp-profile
  Real D = 1. - 10.*xVal*xVal*xVal + 15.*xVal*xVal*xVal*xVal
    - 6.*xVal*xVal*xVal*xVal*xVal;
  // Before and after the ramp, the wall is flat (i.e. D is constant)
  if (xVal < 0.) { D = 1.; }
  else if (xVal > 1.) { D = 0.; }

  // These terms are for relaxing the ramp-profile to a straight line
  const Real a = hValStretch;
  const Real b = a*((etaMin*etaMin - etaMax*etaMax)/(etaMax - etaMin))
    + (hValScaleMax - hValScaleMin)/(etaMax - etaMin);
  const Real c = a*etaMin*etaMax
    + ((hValScaleMin*etaMax - hValScaleMax*etaMin)/(etaMax - etaMin));

  const Real E = Hval*D*a + A;
  const Real F = Hval*D*b + B;

  const Real M = Hval*(a*a_Xi[1]*a_Xi[1] + b*a_Xi[1] + c);
  Real derivD = -30.*xVal*xVal + 60.*xVal*xVal*xVal - 30.*xVal*xVal*xVal*xVal;
  if ((xVal < 0.) || (xVal > 1.)) { derivD = 0.; }

  Real value = 0.0;

  if (a_dirX == 0)
    {
      if (a_dirXi == 0) { value = xSlope; }
      else if (a_dirXi == 1) { value = 0.; }
      else if (a_dirXi == 2) { value = 0.; }
    }
  else if (a_dirX == 1)
    {
      if (a_dirXi == 0) { value = M*derivD*xSlope; }
      else if (a_dirXi == 1) { value = F + 2.*E*a_Xi[1]; }
      else if (a_dirXi == 2) { value = 0.; }
    }
#if CH_SPACEDIM == 3
  else if (a_dirX == 2)
    {
      if (a_dirXi == 0) { value = 0.; }
      else if (a_dirXi == 1) { value = 0.; }
      else if (a_dirXi == 2) { value = zSlope; }
    }
#endif

  return value;
}

#include "NamespaceFooter.H"
