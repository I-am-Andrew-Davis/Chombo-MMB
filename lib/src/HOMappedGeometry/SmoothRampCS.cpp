#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "SmoothRampCS.H"
#include "SmoothRampBlocksCS.H"
#include "MBMiscUtil.H"
#include "CONSTANTS.H"

#include <cfloat>

#include "NamespaceHeader.H"

//Offset between blocks in computational space # of ghost cells + some cushion
static constexpr IntVect c_blockOffset = IntVect(D_DECL(8*3,0,0));

SmoothRampCS::SmoothRampCS()
{
}

SmoothRampCS::~SmoothRampCS()
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
SmoothRampCS::define(const ProblemDomain& a_levelDomain,
                     const RealVect&      a_dx,
                     const Real           a_stretch,
                     const Real           a_rampHeightRelax,
                     const Real           a_zWidth,
                     const int            a_onlyInletBlock,
                     const int            a_numXCellsInletBlock,
                     const int            a_numXCellsRampBlockPreRamp,
                     const int            a_numXCellsRampBlockRamp,
                     const int            a_numXCellsRampBlockPostRamp,
                     const int            a_numXCellsOutletBlock)
{
  m_numBlocks = 3; // Inlet + ramp + outlet
  if (a_onlyInletBlock)
    {
      m_numBlocks = 1; // Inlet only
    }

  // set coordSys for each block before heading off to the big define
  m_coordSysVect.resize(m_numBlocks, NULL);
  m_mappingBlocks.resize(m_numBlocks);

  // Computational grid spacing
  m_dxVect = a_dx;

  // Set the domain width
  m_zWidth = a_zWidth;

  // !!!FIXME!!! this should be an argument
  IntVect numGhost = 8*IntVect::Unit;
  // Origin of block -- will be moved around throughout initialization
  IntVect blockOrigin(IntVect::Zero);

  // The base level must always be 1, as set in the constructor
  int refRatio = 1./a_dx[0];
  pout() << "REFRATIO: " << refRatio << std::endl;
  pout() << "DX: " << a_dx << std::endl;

  // We can take the number of cells in the y-direction and the z-direction
  // and use these for the entire domain. The x-direction number of cells
  // is specific to each block and must be specified by the user.

  // The domain height is fixed and the domain width is uniform. The block
  // lengths are individually variable, except that we need this to be as
  // easy for the user to specify and keep track of as possible. There are a
  // few possibilities here:
  // (1) the lengths of the individual blocks with the constraint that the
  //     ramp itself must be unit length
  // (2) the dx values of the individual blocks except that we require all
  //     blocks to have the same dx for convenience and quality here
  //     (further changes could allow for a stretching, but this will not
  //      be included in the current setup)
  // (3) the number of cells in each block
  // With user-convenience in mind, it would seem that the best approach would
  // be to allow the user to specify the number of cells covering the length
  // of the ramp. This will then define the dx value in all of the blocks.
  // After that, it would be natural for the user to define the number of cells
  // in the small start region of the ramp block and the small end region of the
  // ramp block. Finally, they can set the number of cells in the inlet and
  // outlet blocks.

  for (int blkIdx = 0; blkIdx != m_numBlocks; ++blkIdx)
    {
      // Construct the box defining the block
      // First, we get the size of the current block on this level
      IntVect domainSize = a_levelDomain.domainBox().size();
      //**NOTE: This is currently set up assuming that the inlet block is
      //        always block-0. If there's no ramp, there's no outlet.
      //        If there's a ramp, there's an outlet. So there's only two cases
      //        and we don't have to check which one this is
      if (blkIdx == 0)
        {
          domainSize[0] = refRatio*a_numXCellsInletBlock;
        }
      else if (blkIdx == 1)
        {
          domainSize[0] = refRatio*a_numXCellsRampBlockRamp;
        }
      else if (blkIdx == 2)
        {
          domainSize[0] = refRatio*a_numXCellsOutletBlock;
        }
      // Then, we create an equivalent box with the low corner at zero
      Box blockBox(IntVect_zero, domainSize - 1);
      // Finally, we can set the box for this block by offsetting
      // blockBox as necessary
      m_mappingBlocks[blkIdx] = blockBox + blockOrigin;
      pout() << "  adding block(" << blkIdx << "): "
             << m_mappingBlocks[blkIdx] << std::endl;

      IntVect small = m_mappingBlocks[blkIdx].smallEnd();

      const Real physDx = 1./(refRatio*a_numXCellsRampBlockRamp
                              - refRatio*a_numXCellsRampBlockPreRamp
                              - refRatio*a_numXCellsRampBlockPostRamp);
      const Real xLengthInletBlock = physDx*refRatio*a_numXCellsInletBlock;
      const Real xLengthRampBlockPreRamp =
        physDx*refRatio*a_numXCellsRampBlockPreRamp;
      const Real xLengthRampBlock = refRatio*physDx*(a_numXCellsRampBlockRamp);
      const Real xLengthOutletBlock = refRatio*physDx*a_numXCellsOutletBlock;

      SmoothRampBlocksCS* thisBlockPtr = 
        new SmoothRampBlocksCS(m_dxVect,
                               m_mappingBlocks[blkIdx],
                               blkIdx,
                               numGhost,
                               small,
                               a_stretch,
                               a_rampHeightRelax,
                               a_zWidth,
                               xLengthInletBlock,
                               xLengthRampBlockPreRamp,
                               xLengthRampBlock,
                               xLengthOutletBlock,
                               (refRatio*a_numXCellsRampBlockRamp));
      m_coordSysVect[blkIdx] = reinterpret_cast<NewCoordSys*>(thisBlockPtr);

      // Increment next computational start point
      blockOrigin[0] += m_mappingBlocks[blkIdx].size()[0]
        + refRatio*c_blockOffset[0];
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
SmoothRampCS::defineBoundaries()
{
  // Start by setting all block boundaries to be of type BOUNDARY.
  setAllBoundaries(BlockBoundary::BOUNDARY);

  IndicesTransformation it;
  RigidTransformation rigidTransform;
  const int refRatio = 1./m_dxVect[0];
  if (m_numBlocks == 1) // Only define the inlet block
    {
      if (SpaceDim == 3) // Only the spanwise is periodic
        {
          IntVect baseLength = m_mappingBlocks[0].size();
          // Low z-side
          BlockBoundary& zLo = boundary(0, 2, Side::Lo);
          it.defineFromTranslation(baseLength[2]*BASISV(2));
          zLo.define(it, 0);
          zLo.defineConformal(0, true);
          RealVect physShift = RealVect_zero;
          physShift[2] = m_zWidth;
          rigidTransform.define(RealVect_zero, physShift);
          zLo.definePhysicalTransformation(rigidTransform);

          // High z-side
          BlockBoundary& zHi = boundary(0, 2, Side::Hi);
          it.defineFromTranslation(-baseLength[2]*BASISV(2));
          zHi.define(it, 0);
          zHi.defineConformal(0, true);
          physShift = RealVect_zero;
          physShift[2] = -m_zWidth;
          rigidTransform.define(RealVect_zero, physShift);
          zHi.definePhysicalTransformation(rigidTransform);
        }
    }
  else // Define the inlet, ramp, and outlet blocks
    {
      // Inlet block (block 0)
      IntVect baseLength = m_mappingBlocks[0].size();
      {
        const int block = 0;
        // High x-side
        it.defineFromTranslation(refRatio*c_blockOffset*BASISV(0));
        BlockBoundary& xHi = boundary(block, 0, Side::Hi);
        xHi.define(it, 1);
        xHi.defineConformal(1);

#if CH_SPACEDIM == 3
        // Low z-side
        BlockBoundary& zLo = boundary(block, 2, Side::Lo);
        it.defineFromTranslation(baseLength[2]*BASISV(2));
        zLo.define(it, block);
        zLo.defineConformal(block, true);
        RealVect physShift = RealVect_zero;
        physShift[2] = m_zWidth;
        rigidTransform.define(RealVect_zero, physShift);
        zLo.definePhysicalTransformation(rigidTransform);

        // High z-side
        BlockBoundary& zHi = boundary(block, 2, Side::Hi);
        it.defineFromTranslation(-baseLength[2]*BASISV(2));
        zHi.define(it, block);
        zHi.defineConformal(block, true);
        physShift = RealVect_zero;
        physShift[2] = -m_zWidth;
        rigidTransform.define(RealVect_zero, physShift);
        zHi.definePhysicalTransformation(rigidTransform);
#endif
      }

      // Ramp block (block 1)
      baseLength = m_mappingBlocks[1].size();
      {
        const int block = 1;
        // Low x-side
        it.defineFromTranslation(-refRatio*c_blockOffset*BASISV(0));
        BlockBoundary& xLo = boundary(block, 0, Side::Lo);
        xLo.define(it, 0);
        xLo.defineConformal(0);

        // High x-side
        it.defineFromTranslation(refRatio*c_blockOffset*BASISV(0));
        BlockBoundary& xHi = boundary(block, 0, Side::Hi);
        xHi.define(it, 2);
        xHi.defineConformal(2);

#if CH_SPACEDIM == 3
        // Low z-side
        BlockBoundary& zLo = boundary(block, 2, Side::Lo);
        it.defineFromTranslation(baseLength[2]*BASISV(2));
        zLo.define(it, block);
        zLo.defineConformal(block, true);
        RealVect physShift = RealVect_zero;
        physShift[2] = m_zWidth;
        rigidTransform.define(RealVect_zero, physShift);
        zLo.definePhysicalTransformation(rigidTransform);

        // High z-side
        BlockBoundary& zHi = boundary(block, 2, Side::Hi);
        it.defineFromTranslation(-baseLength[2]*BASISV(2));
        zHi.define(it, block);
        zHi.defineConformal(block, true);
        physShift = RealVect_zero;
        physShift[2] = -m_zWidth;
        rigidTransform.define(RealVect_zero, physShift);
        zHi.definePhysicalTransformation(rigidTransform);
#endif
      }

      // Outlet block (block 2)
      baseLength = m_mappingBlocks[2].size();
      {
        const int block = 2;
        // Low x-side
        it.defineFromTranslation(-refRatio*c_blockOffset*BASISV(0));
        BlockBoundary& xLo = boundary(block, 0, Side::Lo);
        xLo.define(it, 1);
        xLo.defineConformal(1);

#if CH_SPACEDIM == 3
        // Low z-side
        BlockBoundary& zLo = boundary(block, 2, Side::Lo);
        it.defineFromTranslation(baseLength[2]*BASISV(2));
        zLo.define(it, block);
        zLo.defineConformal(block, true);
        RealVect physShift = RealVect_zero;
        physShift[2] = m_zWidth;
        rigidTransform.define(RealVect_zero, physShift);
        zLo.definePhysicalTransformation(rigidTransform);

        // High z-side
        BlockBoundary& zHi = boundary(block, 2, Side::Hi);
        it.defineFromTranslation(-baseLength[2]*BASISV(2));
        zHi.define(it, block);
        zHi.defineConformal(block, true);
        physShift = RealVect_zero;
        physShift[2] = -m_zWidth;
        rigidTransform.define(RealVect_zero, physShift);
        zHi.definePhysicalTransformation(rigidTransform);
#endif
      }
    }

  m_gotBoundaries = true;
}

void SmoothRampCS::blockRemapping(RealVect&            a_xi_valid,
                                  int&                 a_n_valid,
                                  bool&                a_validExists,
                                  RigidTransformation& a_extraDispl,
                                  const RealVect&      a_xiSrc,
                                  const IntVect&       a_iSrc,
                                  int                  a_nSrc) const
{
  const int refRatio = 1./m_dxVect[0];
  // Check if we're in the current block
  if (m_mappingBlocks[a_nSrc].contains(a_iSrc))
    {
      a_xi_valid = a_xiSrc;
      a_n_valid = a_nSrc;
      pout() << "How often does this happen?" << std::endl;
      return;
    }
  // Check if the x-index is in the current block
  if ((a_iSrc[0] >= (m_mappingBlocks[a_nSrc].smallEnd()[0])) &&
      (a_iSrc[0] <= (m_mappingBlocks[a_nSrc].bigEnd()[0])))
    {
      IntVect destIntVect = a_iSrc;
#if CH_SPACEDIM == 3 // Only address the periodic boundaries
      Box blockBox = m_mappingBlocks[a_nSrc];
      destIntVect[2] = (a_iSrc[2] + blockBox.size()[2]) % blockBox.size()[2];
#endif
      // Now we need to transform this into a xi value
      a_xi_valid = destIntVect*m_dxVect;
      a_n_valid = a_nSrc;
      return;
    }
  // Check if the x-index is lower than the current block
  if (a_iSrc[0] < (m_mappingBlocks[a_nSrc].smallEnd()[0]))
    {
      // If this is the inlet block, don't shift it in x, but check if it is
      // across a periodic boundary
      if (a_nSrc == 0)
        {
          IntVect destIntVect = a_iSrc;
#if CH_SPACEDIM == 3
          Box blockBox = m_mappingBlocks[a_nSrc];
          destIntVect[2] =
            (a_iSrc[2] + blockBox.size()[2]) % blockBox.size()[2];
#endif
          // Now we need to transform this into a xi value
          a_xi_valid = destIntVect*m_dxVect;
          a_n_valid = a_nSrc;
          return;
        }
      else // If this is the ramp or outlet block, shift it one block lower
        {
          IntVect destIntVect = a_iSrc;
          destIntVect[0] = a_iSrc[0] - refRatio*c_blockOffset[0];
#if CH_SPACEDIM == 3 // Also check if it's across a periodic boundary
          Box blockBox = m_mappingBlocks[a_nSrc];
          destIntVect[2] =
            (a_iSrc[2] + blockBox.size()[2]) % blockBox.size()[2];
#endif
          // Now we need to transform this into a xi value
          a_xi_valid = destIntVect*m_dxVect;
          a_n_valid = a_nSrc - 1;
          return;
        }
    }
  // Check if the x-index is higher than the current block
  if (a_iSrc[0] > (m_mappingBlocks[a_nSrc].bigEnd()[0]))
    {
      // If this is the outlet block, don't shift it in x, but check if it is
      // across a periodic boundary
      if (a_nSrc == 2)
        {
          IntVect destIntVect = a_iSrc;
#if CH_SPACEDIM == 3
          Box blockBox = m_mappingBlocks[a_nSrc];
          destIntVect[2] =
            (a_iSrc[2] + blockBox.size()[2]) % blockBox.size()[2];
#endif
          // Now we need to transform this into a xi value
          a_xi_valid = destIntVect*m_dxVect;
          a_n_valid = a_nSrc;
          return;
        }
      else // If this is the ramp or outlet block, shift it one block lower
        {
          IntVect destIntVect = a_iSrc;
          destIntVect[0] = a_iSrc[0] + refRatio*c_blockOffset[0];
#if CH_SPACEDIM == 3 // Also check if it's across a periodic boundary
          Box blockBox = m_mappingBlocks[a_nSrc];
          destIntVect[2] =
            (a_iSrc[2] + blockBox.size()[2]) % blockBox.size()[2];
#endif
          // Now we need to transform this into a xi value
          a_xi_valid = destIntVect*m_dxVect;
          a_n_valid = a_nSrc + 1;
          return;
        }
    }

  a_validExists = true;
}

// -- begin factory implementations ---------------------------

MultiBlockCoordSys*
SmoothRampCSFactory::getCoordSys(const ProblemDomain& a_levelDomain,
                                 const RealVect& a_dx) const
{
  SmoothRampCS* coordSysPtr = new SmoothRampCS();
  coordSysPtr->define(a_levelDomain,
                      a_dx,
                      m_stretch,
                      m_rampHeightRelax,
                      m_zWidth,
                      m_onlyInletBlock,
                      m_numXCellsInletBlock,
                      m_numXCellsRampBlockPreRamp,
                      m_numXCellsRampBlockRamp,
                      m_numXCellsRampBlockPostRamp,
                      m_numXCellsOutletBlock);
  return (static_cast <MultiBlockCoordSys*> (coordSysPtr));
}

#include "NamespaceFooter.H"
