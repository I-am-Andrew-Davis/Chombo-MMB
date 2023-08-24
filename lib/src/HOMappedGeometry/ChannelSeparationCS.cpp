#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "ChannelSeparationCS.H"
#include "CartesianBlockCS.H"
#include "CartesianCS.H"
#include "BoxIterator.H"
#include "ReadCGNS.H"
#include "ExternalCS.H"
#include "DebugOut.H"
#include <queue>

#include "NamespaceHeader.H"

//Offset between blocks in computational space # of ghost cells + some cushin
static constexpr IntVect c_blockOffset = IntVect(D_DECL(8*3,0,0));

ChannelSeparationCS::ChannelSeparationCS()
{
 // m_coordSysVect.resize(NUMBLOCKS, NULL);

 // m_mappingBlocks.resize(NUMBLOCKS);
}

ChannelSeparationCS::~ChannelSeparationCS()
{
}

void
ChannelSeparationCS::define(const ProblemDomain& a_levelDomain,
                            const RealVect& a_dx,
                            const Real a_z_dx,
                            const Real a_l_dx,
                            const Real a_r_dx,
                            const int a_z_length,
                            const int a_l_length,
                            const int a_r_length,
                            const int a_blocksAreConformal,
                            const int a_numRampBlocks)
{
  CH_TIME("ChannelSeparationCS::define");
  if (m_verbosity > 1)
    {
      pout() << "ChannelSeparationCS begin defining coordinate system"
             << std::endl;
    }
  CH_assert(m_gotExternalInfo);

#ifdef CH_USE_CGNS
  ReadCGNS cgns;
  if (cgns.openFile(m_fileName))
    {
      pout() << "Unable to open CGNS file " << m_fileName << "!  Aborting."
          << std::endl;
      MayDay::Error("Require CGNS file");
    }
  if (cgns.numBase() < 1)
    {
      MayDay::Error("Require at least 1 base node in CGNS file");
    }
  // We always select base 1.  Here the map from zone names to block indices is
  // also constructed
  m_zoneNameMap.clear();
  cgns.selectBase(1, m_zoneNameMap);
  m_numBlocks = a_numRampBlocks;
  NUMBLOCKS = a_numRampBlocks;
  m_coordSysVect.resize(m_numBlocks, nullptr);
  m_mappingBlocks.resize(m_numBlocks);
  m_dxVect = a_dx;
  //pout()<<"a_dx: " << a_dx << std::endl;
  m_blocksAreConformal = a_blocksAreConformal;
  m_numZCells = a_z_length;
  m_zDx = a_z_dx;

  // the base level must always be 1, as set in the constructor
  int refRatio = 1./a_dx[0];
  pout() << "REFRATIO: " << refRatio << std::endl;
  pout() << "DX: " << a_dx << std::endl;

  IntVect numGhost = 8*IntVect::Unit;

  // Select procs to distribute ExternalCS (single block CS) construction
  int procPerBlock = numProc()/m_numBlocks;
  if (m_numBlocks == 1)
    {
      procPerBlock = std::max((procPerBlock - 1), 1);
    }
  int startProc = 0;
  //NOTE: the y cell num and dx are always 1 because the values used do not come
  //from these int/realvects they come from the cgns file
  //left extDx
  const RealVect extDx_l  =RealVect(D_DECL(a_l_dx, 1, a_z_dx));
  //middle extDx
  const RealVect extDx_m  =RealVect(D_DECL(1, 1, a_z_dx)); //note: the dx spacing for the first component does not matter
  //Right extDx
  const RealVect extDx_r  =RealVect(D_DECL(a_r_dx, 1, a_z_dx));
  //left ext distance
  const IntVect extCell_l =IntVect (D_DECL(a_l_length + 1, 1, a_z_length + 1));
  //middle ext distance
  const IntVect extCell_m =IntVect (D_DECL(1, 1, a_z_length + 1)); //number of points ex (9,0) == 8 cells
  //right ext distance
  const IntVect extCell_r =IntVect (D_DECL(a_r_length + 1, 1, a_z_length + 1));
  const int constIdx = 0;

  for (int idxBlk = 0; idxBlk != m_numBlocks; ++idxBlk)
    {
      pout() << "idxBlk: " << idxBlk<< std::endl;
      // Construct the box defining the block
      std::string zoneName;
      IntVect zoneNumCell;
      cgns.readZoneInfo(1, zoneName, zoneNumCell);
      if (idxBlk == 0)
        {
         zoneName = "blk-0";
         IntVect blockOrigin(IntVect::Zero);
         Box blockBox(IntVect::Zero,
                      IntVect(D_DECL(extCell_l[0]-2,
                                     zoneNumCell[1]-1,
                                     extCell_l[2]-2)));
         blockBox.refine(refRatio);
         m_mappingBlocks[idxBlk] = blockBox + blockOrigin;

         pout() << "  adding block(" << zoneName << ", " << idxBlk << "): "
                << m_mappingBlocks[idxBlk] << std::endl;
         m_coordSysVect[idxBlk] = new ExtrudedExternalCS(cgns,
                                                         a_dx,
                                                         m_scale,
                                                         m_mappingBlocks[idxBlk],
                                                         m_fileName,
                                                         zoneName,
                                                         constIdx,
                                                         numGhost,
                                                         startProc,
                                                         procPerBlock,
                                                         idxBlk,
                                                         extDx_l,
                                                         extCell_l);
        }

      if (idxBlk == 1)
      {
       zoneName = "blk-1";
       // Construct the box defining the block
       IntVect start = IntVect::Zero;
       start[0] = m_mappingBlocks[0].bigEnd(0)+1;
       Box blockBox(IntVect::Zero,
                    IntVect(D_DECL(zoneNumCell[0] -1,zoneNumCell[1]-1, extCell_m[2]-2 )));
       IntVect blockOrigin(start + refRatio*c_blockOffset);
       blockBox.refine(refRatio);
       m_mappingBlocks[idxBlk] = blockBox + blockOrigin;

       pout() << "  adding block(" << zoneName << ", " << idxBlk << "): "
              << m_mappingBlocks[idxBlk] << std::endl;
       m_coordSysVect[idxBlk] = new ExtrudedExternalCS(cgns,
                                                       a_dx,
                                                       m_scale,
                                                       m_mappingBlocks[idxBlk],
                                                       m_fileName,
                                                       zoneName,
                                                       constIdx,
                                                       numGhost,
                                                       startProc,
                                                       procPerBlock,
                                                       idxBlk,
                                                       extDx_m,
                                                       extCell_m);
      }
      if (idxBlk == 2)
      {
       zoneName = "blk-2";
       // Construct the box defining the block
       IntVect start = IntVect::Zero;
       start[0] = m_mappingBlocks[1].bigEnd(0)+1;
       Box blockBox(IntVect::Zero,
                    IntVect(D_DECL(extCell_r[0]-2, zoneNumCell[1] -1, extCell_r[2]-2)));
       // pout() << "box size: " << blockBox.size() << std::endl;
       IntVect blockOrigin(start + refRatio*c_blockOffset);
       blockBox.refine(refRatio);
       m_mappingBlocks[idxBlk] = blockBox + blockOrigin;

       pout() << "  adding block(" << zoneName << ", " << idxBlk << "): "
               << m_mappingBlocks[idxBlk] << std::endl;
       m_coordSysVect[idxBlk] = new ExtrudedExternalCS(cgns,
                                                       a_dx,
                                                       m_scale,
                                                       m_mappingBlocks[idxBlk],
                                                       m_fileName,
                                                       zoneName,
                                                       constIdx,
                                                       numGhost,
                                                       startProc,
                                                       procPerBlock,
                                                       idxBlk,
                                                       extDx_r,
                                                       extCell_r);
      }
    }
  m_gotMappingBlocks = true;

  defineBoundaries();
  initializeBlockTransformations();

  m_gotCoordSysVect = true;
  setCache();

  if (m_verbosity > 1)
  {
      pout() << "ExternalMultiBlockCS finished defining coordinate system" << std::endl;
  }
  cgns.closeFile();
#else
  MayDay::Error("ExternalMultiBlockCS requires CGNS support");
#endif // #ifdef USE_CGNS
}

void
ChannelSeparationCS::defineBoundaries()
{
 CH_assert(gotMappingBlocks());
 m_boundaries.resize(NUMBLOCKS);
 IndicesTransformation it;
 // The base level must always be 1, as set in the constructor
 int refRatio = 1./m_dxVect[0];
  for (int iblock = 0; iblock < NUMBLOCKS; ++iblock)
  {
   //First Block
   if (iblock == 0)
   {
    //The x-low direction
    BlockBoundary& bbxl = boundary(iblock, 0, Side::Lo);
    bbxl.define(iblock);
    //The x-hi direction
    if (NUMBLOCKS == 3)
      {
        it.defineFromTranslation(refRatio*(c_blockOffset)*BASISV(0));
        BlockBoundary& bb = boundary(iblock, 0, Side::Hi);
        bb.define(it, 1);
        if (m_blocksAreConformal)
          {
            bb.defineConformal(1);
          }
      }
    else
      {
        BlockBoundary& bbxh = boundary(iblock, 0, Side::Hi);
        bbxh.define(iblock);
      }
    //The y-low direction
    BlockBoundary& bbyl = boundary(iblock, 1, Side::Lo);
    bbyl.define(iblock);
    //The y-hi direction
          BlockBoundary& bbyh = boundary(iblock, 1, Side::Hi);
          bbyh.define(iblock);
#if CH_SPACEDIM  == 3

          //The z-low direction
          BlockBoundary& bbzl = boundary(iblock, 2, Side::Lo);
          it.defineFromTranslation(refRatio*m_numZCells*BASISV(2));
          bbzl.define(it, iblock);
          bbzl.defineConformal(iblock, true);
          RigidTransformation rigidTransform;
          RealVect physShift = RealVect_zero;
          physShift[2] = m_numZCells*m_zDx;
          rigidTransform.define(RealVect_zero, physShift);
          bbzl.definePhysicalTransformation(rigidTransform);
          //The z-hi direction
          BlockBoundary& bbzh = boundary(iblock, 2, Side::Hi);
          it.defineFromTranslation(-refRatio*m_numZCells*BASISV(2));
          bbzh.define(it, iblock);
          bbzh.defineConformal(iblock, true);
          physShift[2] = -m_numZCells*m_zDx;
          rigidTransform.define(RealVect_zero, physShift);
          bbzh.definePhysicalTransformation(rigidTransform);
#endif
   }

   //Second Block
   if (iblock == 1)
   {
    //The x-lo direction
          it.defineFromTranslation((-1 * refRatio*c_blockOffset)*BASISV(0));
          BlockBoundary& bbl = boundary(iblock, 0, Side::Lo);
          bbl.define(it, 0);
          bbl.definePhysicalTransformation
          (RigidTransformation(RealVect_zero, RealVect_zero));
          if (m_blocksAreConformal)
            {
              bbl.defineConformal(0);
            }
          it.defineFromTranslation((1 * refRatio*c_blockOffset)*BASISV(0));
          //The x-hi direction
          BlockBoundary& bbh = boundary(iblock, 0, Side::Hi);
          bbh.define(it, 2);
          if (m_blocksAreConformal)
            {
              bbh.defineConformal(2);
            }
          //The y-low direction
          BlockBoundary& bbyl = boundary(iblock, 1, Side::Lo);
          bbyl.define(iblock);
          //The y-hi direction
          BlockBoundary& bbyh = boundary(iblock, 1, Side::Hi);
          bbyh.define(iblock);
#if CH_SPACEDIM  == 3
          //The z-low direction
          BlockBoundary& bbzl = boundary(iblock, 2, Side::Lo);
          it.defineFromTranslation(refRatio*m_numZCells*BASISV(2));
          bbzl.define(it, iblock);
          bbzl.defineConformal(iblock, true);
          RigidTransformation rigidTransform;
          RealVect physShift = RealVect_zero;
          physShift[2] = m_numZCells*m_zDx;
          rigidTransform.define(RealVect_zero, physShift);
          bbzl.definePhysicalTransformation(rigidTransform);
          //The z-hi direction
          BlockBoundary& bbzh = boundary(iblock, 2, Side::Hi);
          it.defineFromTranslation(-refRatio*m_numZCells*BASISV(2));
          bbzh.define(it, iblock);
          bbzh.defineConformal(iblock, true);
          physShift[2] = -m_numZCells*m_zDx;
          rigidTransform.define(RealVect_zero, physShift);
          bbzh.definePhysicalTransformation(rigidTransform);
#endif
   }

   //Third and final block
   if (iblock == 2)
   {
    //The x-lo direction
    it.defineFromTranslation((-1 * refRatio*c_blockOffset)*BASISV(0));
    BlockBoundary& bb = boundary(iblock, 0, Side::Lo);
    bb.define(it, 1);
    bb.definePhysicalTransformation
    (RigidTransformation(RealVect_zero, RealVect_zero));
    if (m_blocksAreConformal)
      {
        bb.defineConformal(1);
      }
    //The x-hi direction
    BlockBoundary& bbxh = boundary(iblock, 0, Side::Hi);
    bbxh.define(iblock);
    //The y-low direction
    BlockBoundary& bbyl = boundary(iblock, 1, Side::Lo);
    bbyl.define(iblock);
    //The y-hi direction
    BlockBoundary& bbyh = boundary(iblock, 1, Side::Hi);
    bbyh.define(iblock);
#if CH_SPACEDIM  == 3
    //The z-low direction
    BlockBoundary& bbzl = boundary(iblock, 2, Side::Lo);
    it.defineFromTranslation(refRatio*m_numZCells*BASISV(2));
    bbzl.define(it, iblock);
    bbzl.defineConformal(iblock, true);
    RigidTransformation rigidTransform;
    RealVect physShift = RealVect_zero;
    physShift[2] = m_numZCells*m_zDx;
    rigidTransform.define(RealVect_zero, physShift);
    bbzl.definePhysicalTransformation(rigidTransform);
    //The z-hi direction
    BlockBoundary& bbzh = boundary(iblock, 2, Side::Hi);
    it.defineFromTranslation(-refRatio*m_numZCells*BASISV(2));
    bbzh.define(it, iblock);
    bbzh.defineConformal(iblock, true);
    physShift[2] = -m_numZCells*m_zDx;
    rigidTransform.define(RealVect_zero, physShift);
    bbzh.definePhysicalTransformation(rigidTransform);
#endif
   }
  }

  m_gotBoundaries = true;
}
/*--------------------------------------------------------------------*/
///  Return the ChannelCS object
/**
 *   \param[in]  a_domain
 *                      Problem domain the coordinate system is in
 *   \return            Pointer to the coordinate system
 *//*-----------------------------------------------------------------*/

ExternalMultiBlockCS*
ChannelSeparationCSFactory::getCoordSys(const ProblemDomain& a_levelDomain,
                                        const RealVect& a_dx) const
{
 ChannelSeparationCS* coordSysPtr = new ChannelSeparationCS();
 coordSysPtr->setupForExternalGrid(m_fileName, m_scale);
 coordSysPtr->setVerbosity(m_verbosity);
 coordSysPtr->define( a_levelDomain, a_dx, m_z_dx, m_l_dx, m_r_dx,
                      m_z_length, m_l_length, m_r_length,
                      m_blocksAreConformal, m_numRampBlocks);
 // coordSysPtr->initExt(a_extDx, a_extCell, a_readDir, a_readSide)

 return (static_cast <ExternalMultiBlockCS*> (coordSysPtr));
}

#include "NamespaceFooter.H"
