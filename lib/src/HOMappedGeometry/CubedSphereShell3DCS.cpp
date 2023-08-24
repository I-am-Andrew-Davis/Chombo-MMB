#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "CubedSphereShell3DCS.H"
#include "CubedSphereShellPanel3DCS.H"
#include "CubedSphereShellF_F.H"
#include "MBMiscUtil.H"
#include "CONSTANTS.H"

#include <cfloat>

#include "NamespaceHeader.H"
// #include "UsingNamespace.H"

///////////////////////////////////////////////////////////////////////////////

CubedSphereShell3DCS::CubedSphereShell3DCS()
{
}

CubedSphereShell3DCS::~CubedSphereShell3DCS()
{
  if (m_gotCoordSysVect)
    {
      for (int i=0; i != m_numBlocks; ++i)
        {
          delete m_coordSysVect[i];
        }
    }
}

void
CubedSphereShell3DCS::define(const ProblemDomain& a_levelDomain,
                             const RealVect&      a_dx,
                             const Real           a_radiusInner,
                             const Real           a_radiusOuter,
                             const Real           a_alpha,
                             const Real           a_beta)
{
  // Must always be SpaceDim*2 blocks in this mapping
  m_numBlocks = SpaceDim*2;

  // set coordSys for each block before heading off to the big define
  m_coordSysVect.resize(m_numBlocks, NULL);
  m_mappingBlocks.resize(m_numBlocks);

  // Computational grid spacing
  m_dxVect = a_dx;

  // !!!FIXME!!! this should be an argument
  IntVect numGhost = 8*IntVect::Unit;
  // Origin of block -- will be moved around throughout initialization
  IntVect blockOrigin(IntVect::Zero);

  // Check that alpha_eta > |a_alpha*a_beta|
  Real alpha_eta = (a_radiusOuter - a_radiusInner) - a_alpha*std::tanh(a_beta);
  if (alpha_eta <= std::abs(a_alpha*a_beta))
    {
      pout() << "alpha_eta = " << alpha_eta << "; a_alpha = " << a_alpha
             << "; a_beta = " << a_beta << "; a_alpha*a_beta = "
             << std::abs(a_alpha*a_beta) << endl;
      MayDay::Error(
        "CubedSphereShell3DCS::define: alpha_eta must be > "
        "abs(a_alpha*a_beta)");
    }

  for (int blkIdx = 0; blkIdx != m_numBlocks; ++blkIdx)
    {
      // Construct the box defining the block
      // First, we get the size of the current block on this level
      IntVect domainSize = a_levelDomain.domainBox().size();
      // Then, we create an equivalent box with the low corner at zero
      Box blockBox(IntVect_zero, domainSize - 1);
      // Finally, we can set the box for this block by offsetting
      // blockBox as necessary
      m_mappingBlocks[blkIdx] = blockBox + blockOrigin;
      pout() << "  adding block(" << blkIdx << "): "
             << m_mappingBlocks[blkIdx] << std::endl;

      IntVect small = m_mappingBlocks[blkIdx].smallEnd();

      //**FIXME: give inner/outer radius to CubedSphereShellPanel3DCS
      CubedSphereShellPanel3DCS* thisBlockPtr = 
        new CubedSphereShellPanel3DCS(m_dxVect,
                                      m_mappingBlocks[blkIdx],
                                      blkIdx,
                                      numGhost,
                                      small,
                                      a_radiusInner,
                                      a_radiusOuter,
                                      a_alpha,
                                      a_beta);
      m_coordSysVect[blkIdx] = reinterpret_cast<NewCoordSys*>(thisBlockPtr);

      // Increment next computational start point
      // Note that the smallest block we can have is size-8, and we need
      // at least 8 ghost cells in index space for each block, so we just
      // set the distance to the next block to be 3 times the block size
      blockOrigin[0] += blockBox.size()[0] + 3*blockBox.size()[0];
    }
  m_gotCoordSysVect = true;
  m_gotMappingBlocks = true;

  // Define boundaries
  defineBoundaries();

  // Initialize block transformations
  initializeBlockTransformations();

  // Set cache for MultiBlockCoordSys
  setCache();
}

void
CubedSphereShell3DCS::defineBoundaries()
{
  // Start by setting all block boundaries to be of type BOUNDARY.
  setAllBoundaries(BlockBoundary::BOUNDARY);

  // length in dimension 0 of central block (index 0)
  IntVect baseLength = m_mappingBlocks[0].size();

  if (SpaceDim == 3)
    {
      constexpr int xLo = 0;
      constexpr int xHi = SpaceDim;
      constexpr int yLo = 1;
      constexpr int yHi = 1 + SpaceDim;
      constexpr int zLo = 2;
      constexpr int zHi = 2 + SpaceDim;
      const IntVect xShift = baseLength[0]*BASISV(0);
      const IntVect xGapShift = 3*baseLength[0]*BASISV(0);
      const IntVect yShift = baseLength[1]*BASISV(1);
      const IntVect zShift = baseLength[2]*BASISV(2);
      constexpr IntVect permSame = IntVect{0, 1, 2, 0, 0, 0};
      constexpr IntVect permDiff = IntVect{2, 1, 0, 0, 0, 0};
      constexpr IntVect clockwise = IntVect{-1, 1, 1, 0, 0, 0};
      constexpr IntVect anticlockwise = IntVect{1, 1, -1, 0, 0, 0};
      constexpr IntVect flip = IntVect{-1, 1, -1, 0, 0, 0};
      IndicesTransformation it;

      // Equatorial panel 0
      it.defineFromTranslation(xGapShift);
      m_boundaries[0][xHi].define(it, 1);
      it.defineFromTranslation(4*xShift + 3*xGapShift);
      m_boundaries[0][xLo].define(it, 3);
      it.defineFromPivot(zShift,
                         6*xShift + 5*xGapShift + zShift,
                         permSame,
                         flip);
      m_boundaries[0][zHi].define(it, 5);
      it.defineFromPivot(IntVect_zero,
                         5*xShift + 4*xGapShift + zShift,
                         permDiff,
                         clockwise);
      m_boundaries[0][zLo].define(it, 4);

      // Equatorial panel 1
      it.defineFromTranslation(xGapShift);
      m_boundaries[1][xHi].define(it, 2);
      it.defineFromTranslation(-xGapShift);
      m_boundaries[1][xLo].define(it, 0);
      it.defineFromPivot(xShift + xGapShift + zShift,
                         5*xShift + 5*xGapShift + zShift,
                         permDiff,
                         clockwise);
      m_boundaries[1][zHi].define(it, 5);
      it.defineFromPivot(xShift + xGapShift,
                         5*xShift + 4*xGapShift,
                         permSame,
                         flip);
      m_boundaries[1][zLo].define(it, 4);

      // Equatorial panel 2
      it.defineFromTranslation(xGapShift);
      m_boundaries[2][xHi].define(it, 3);
      it.defineFromTranslation(-xGapShift);
      m_boundaries[2][xLo].define(it, 1);
      it.defineFromPivot(2*xShift + 2*xGapShift + zShift,
                         5*xShift + 5*xGapShift,
                         permSame,
                         IntVect_unit);
      m_boundaries[2][zHi].define(it, 5);
      it.defineFromPivot(2*xShift + 2*xGapShift,
                         4*xShift + 4*xGapShift,
                         permDiff,
                         anticlockwise);
      m_boundaries[2][zLo].define(it, 4);

      // Equatorial panel 3
      it.defineFromTranslation(-4*xShift - 3*xGapShift);
      m_boundaries[3][xHi].define(it, 0);
      it.defineFromTranslation(-xGapShift);
      m_boundaries[3][xLo].define(it, 2);
      it.defineFromPivot(3*xShift + 3*xGapShift + zShift,
                         6*xShift + 5*xGapShift,
                         permDiff,
                         anticlockwise);
      m_boundaries[3][zHi].define(it, 5);
      it.defineFromPivot(3*xShift + 3*xGapShift,
                         4*xShift + 4*xGapShift + zShift,
                         permSame,
                         IntVect_unit);
      m_boundaries[3][zLo].define(it, 4);

      // North polar panel (4)
      it.defineFromPivot(5*xShift + 4*xGapShift + zShift,
                         IntVect_zero,
                         permDiff,
                         anticlockwise);
      m_boundaries[4][xHi].define(it, 0);
      it.defineFromPivot(4*xShift + 4*xGapShift,
                         2*xShift + 2*xGapShift,
                         permDiff,
                         clockwise);
      m_boundaries[4][xLo].define(it, 2);
      it.defineFromPivot(4*xShift + 4*xGapShift + zShift,
                         3*xShift + 3*xGapShift,
                         permSame,
                         IntVect_unit);
      m_boundaries[4][zHi].define(it, 3);
      it.defineFromPivot(5*xShift + 4*xGapShift,
                         xShift + xGapShift,
                         permSame,
                         flip);
      m_boundaries[4][zLo].define(it, 1);

      // South polar panel (5)
      it.defineFromPivot(6*xShift + 5*xGapShift,
                         3*xShift + 3*xGapShift + zShift,
                         permDiff,
                         clockwise);
      m_boundaries[5][xHi].define(it, 3);
      it.defineFromPivot(5*xShift + 5*xGapShift + zShift,
                         xShift + xGapShift + zShift,
                         permDiff,
                         anticlockwise);
      m_boundaries[5][xLo].define(it, 1);
      it.defineFromPivot(6*xShift + 5*xGapShift + zShift,
                         zShift,
                         permSame,
                         flip);
      m_boundaries[5][zHi].define(it, 0);
      it.defineFromPivot(5*xShift + 5*xGapShift,
                         2*xShift + 2*xGapShift + zShift,
                         permSame,
                         IntVect_unit);
      m_boundaries[5][zLo].define(it, 2);
    }
  else if (SpaceDim == 2)
    {
      int xLo = 0;
      int xHi = SpaceDim;
      int yLo = 1;
      int yHi = 1 + SpaceDim;
      IntVect xShift = baseLength[0]*BASISV(0);
      IntVect xGapShift = 3*baseLength[0]*BASISV(0);
      IntVect yShift = baseLength[1]*BASISV(1);
      IntVect permSame = IntVect(D_DECL6(0, 1, 2, 0, 0, 0));
      IntVect permDiff = IntVect(D_DECL6(1, 0, 2, 0, 0, 0));
      IntVect clockwise = IntVect(D_DECL6(-1, 1, 1, 0, 0, 0));
      IntVect anticlockwise = IntVect(D_DECL6(1, -1, 1, 0, 0, 0));
      IntVect flip = IntVect(D_DECL6(-1, -1, 1, 0, 0, 0));
      IndicesTransformation it;

      // Equatorial panel 0
      it.defineFromTranslation(xGapShift);
      m_boundaries[0][xHi].define(it, 1);
      it.defineFromTranslation(4*xShift + 3*xGapShift);
      m_boundaries[0][xLo].define(it, 3);

      // Equatorial panel 1
      it.defineFromTranslation(xGapShift);
      m_boundaries[1][xHi].define(it, 2);
      it.defineFromTranslation(-xGapShift);
      m_boundaries[1][xLo].define(it, 0);

      // Equatorial panel 2
      it.defineFromTranslation(xGapShift);
      m_boundaries[2][xHi].define(it, 3);
      it.defineFromTranslation(-xGapShift);
      m_boundaries[2][xLo].define(it, 1);

      // Equatorial panel 3
      it.defineFromTranslation(-4*xShift - 3*xGapShift);
      m_boundaries[3][xHi].define(it, 0);
      it.defineFromTranslation(-xGapShift);
      m_boundaries[3][xLo].define(it, 2);
    }

  m_gotBoundaries = true;
}

void CubedSphereShell3DCS::blockRemapping(RealVect&            a_xi_valid,
                                          int&                 a_n_valid,
                                          bool&                a_validExists,
                                          RigidTransformation& a_extraDispl,
                                          const RealVect&      a_xiSrc,
                                          const IntVect&       a_iSrc,
                                          int                  a_nSrc) const
{
  // We can use which block we're in to get the Cartesian coordinates
  RealVect xSrc = m_coordSysVect[a_nSrc]->realCoord(a_xiSrc);

#if CH_SPACEDIM == 3
  // Now we can tell which block we should be going to
  if ((std::fabs(xSrc[0]) > std::fabs(xSrc[1])) &&
      (std::fabs(xSrc[0]) > std::fabs(xSrc[2])))
    {
      if (xSrc[0] > 0.)
        {
          a_n_valid = 0;
        }
      else
        {
          a_n_valid = 2;
        }
    }
  else if (std::fabs(xSrc[1]) > std::fabs(xSrc[2]))
    {
      if (xSrc[1] > 0.)
        {
          a_n_valid = 4;
        }
      else
        {
          a_n_valid = 5;
        }
    }
  else
    {
      if (xSrc[2] > 0.)
        {
          a_n_valid = 3;
        }
      else
        {
          a_n_valid = 1;
        }
    }
#elif CH_SPACEDIM == 2
  if (std::fabs(xSrc[0]) > std::fabs(xSrc[1]))
    {
      if (xSrc[0] > 0.)
        {
          a_n_valid = 0;
        }
      else
        {
          a_n_valid = 2;
        }
    }
  else
    {
      if (xSrc[1] > 0.)
        {
          a_n_valid = 3;
        }
      else
        {
          a_n_valid = 1;
        }
    }
#endif

  // Now we can tell where we should end up in mapped space
  a_xi_valid = m_coordSysVect[a_n_valid]->mappedCoord(xSrc);

  a_validExists = true;
}

// -- begin factory implementations ---------------------------

MultiBlockCoordSys*
CubedSphereShell3DCSFactory::getCoordSys(const ProblemDomain& a_levelDomain,
                                         const RealVect& a_dx) const
{
  CubedSphereShell3DCS* coordSysPtr = new CubedSphereShell3DCS();
  coordSysPtr->define(a_levelDomain,
                      a_dx,
                      m_radiusInner,
                      m_radiusOuter,
                      m_alpha,
                      m_beta);
  // setAllPhysical(coordSysPtr);
  return (static_cast <MultiBlockCoordSys*> (coordSysPtr));
}

#include "NamespaceFooter.H"
