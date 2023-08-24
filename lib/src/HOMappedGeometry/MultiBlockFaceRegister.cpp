#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <map>

#include "MultiBlockFaceRegister.H"
#include "NewFourthOrderCoordSys.H"
#include "BoxIterator.H"
#include "CH_assert.H"

#include "PointwiseDotProdF_F.H"

#include "NamespaceHeader.H"


/*--------------------------------------------------------------------*/
//  Constructor
/**
*//*------------------------------------------------------------------*/

MultiBlockFaceRegister::MultiBlockFaceRegister(
  const MultiBlockCoordSys *const a_CrCoordSys,
  const MultiBlockRegions&        a_mbRegions)
  :
  m_coordSys(a_CrCoordSys),
  m_mbRegions(a_mbRegions),
  m_isDefined(FaceRegUndefined),
  m_scaleFineFaces(true)
{
  CH_assert(a_CrCoordSys != nullptr);
}

/*--------------------------------------------------------------------*/
//  Define everything
/**
*//*------------------------------------------------------------------*/

void
MultiBlockFaceRegister::define(const DisjointBoxLayout& a_FnGrid,
                               int                      a_nRefine,
                               int                      a_nComp,
                               bool                     a_scaleFineFaces)
{
  define(a_FnGrid,
         a_nRefine*IntVect_unit,
         a_nComp,
         a_scaleFineFaces);
}


void
MultiBlockFaceRegister::define(const DisjointBoxLayout& a_FnGrid,
                               IntVect                  a_nRefine,
                               int                      a_nComp,
                               bool                     a_scaleFineFaces)
{
  CH_assert(a_nRefine > IntVect_zero);

  m_nRefine = a_nRefine;
  m_scaleFineFaces = a_scaleFineFaces;

  DisjointBoxLayout crFnGrid;
  coarsen(crFnGrid, a_FnGrid, a_nRefine);

  m_fineFlux.define(crFnGrid, a_nComp, IntVect::Unit);
  m_isDefined = FaceRegDefined;
}

/*--------------------------------------------------------------------*/
//  Define the coarse register and copiers
/**
*//*------------------------------------------------------------------*/

//**Not required if using MultiBlockRegions
#if 0
void
MultiBlockFaceRegister::defineCoarse(const DisjointBoxLayout& a_CrGrid)
{
  CH_TIME("MultiBlockFaceRegister::defineCoarse");
  const Vector<Box>& blocks = m_coordSys->mappingBlocks();
  CH_assert(blocks.size() > 1);
  CH_assert(m_isDefined & FluxRegFineDefined);
  m_isDefined |= ( FluxRegCoarseDefined | FluxRegDefined );

//--Setup of base class (since we don't actually call it's define function)

  const DisjointBoxLayout& crFnGrid = m_fineFlux.getBoxes();
  // Note that m_coarFlux and the reverse copier are *not* defined
  for (const int dir : EachDir)
    {
      m_coarseLocations[dir].define(a_CrGrid);
      m_coarseLocations[dir+SpaceDim].define(a_CrGrid);
    }

//--Iterate over all boxes and gather information on those whose boundaries
//--coincide with block boundaries.

  // We construct data structures for accumulating fine faces from remote
  // blocks in each direction.
  Vector<Box> loSrcBoxes, loDstBoxes, hiSrcBoxes, hiDstBoxes;
  Vector<int> loSrcProcs, loDstProcs, hiSrcProcs, hiDstProcs;
  std::map<Box, Box> loSrcToDstMap, hiSrcToDstMap;
  for (const int dir : EachDir)
    {
      loSrcBoxes.clear();
      loDstBoxes.clear();
      hiSrcBoxes.clear();
      hiDstBoxes.clear();
      loSrcProcs.clear();
      loDstProcs.clear();
      hiSrcProcs.clear();
      hiDstProcs.clear();
      loSrcToDstMap.clear();
      hiSrcToDstMap.clear();
      for (LayoutIterator litCr(a_CrGrid.layoutIterator()); litCr.ok(); ++litCr)
        {
          // Find the block containing this box.
          Box crBox = a_CrGrid[litCr];
          int dstProc = a_CrGrid.procID(litCr());
          int blockID = m_coordSys->whichBlock(crBox);
          CH_assert(blockID != -1);
          const Box& blockBox = blocks[blockID];

          // Now find out whether this face abuts the block boundary and set up
          // any block-block face transfers that are needed.
          Box crFaceBox = crBox;
          crFaceBox.grow(dir, 1);

          // Find correction locations at both the low and high sides of the
          // box simultaneously to avoid multiple loops
          bool haveCorrection = false;
          
          // This would be a correction to coarse cells/faces at the low side of
          // crBox.  Note that this equals the high side of the fine box.
          // int nbrBlockIDLo = -1;
          IndicesTransformation transformationLo;
          Box crFaceBoxTrLo;
          Vector<Box>* coarseLocationsLo = nullptr;
          if (crFaceBox.smallEnd(dir) < blockBox.smallEnd(dir))
            {
              haveCorrection = true;
              const BlockBoundary& boundary =
                m_coordSys->boundary(blockID, dir, Side::Lo);
              // nbrBlockIDLo = boundary.neighbor();
              transformationLo = boundary.getTransformation();
              
              // Find any fine boxes on the far side of the boundary.
              // Transform the test box to the remote space
              crFaceBoxTrLo = transformationLo.transformFwd(crFaceBox);

              // This is with respect to the CrFn boxes (following conventions
              // in LevelFluxRegister)
              if (dstProc == procID())
                {
                  DataIndex didxCr(litCr());
                  coarseLocationsLo =
                    &(getCoarseLocations(dir, Side::Hi)[didxCr]);
                }
            }
          
          // This would be a correction to coarse cells/faces at the high side
          // of crBox.  Note that this equals the low side of the fine box.
          // int nbrBlockIDHi = -1;
          IndicesTransformation transformationHi;
          Box crFaceBoxTrHi;
          Vector<Box>* coarseLocationsHi = nullptr;
          if (crFaceBox.bigEnd(dir) > blockBox.bigEnd(dir))
            {
              haveCorrection = true;
              const BlockBoundary& boundary =
                m_coordSys->boundary(blockID, dir, Side::Hi);
              // nbrBlockIDHi = boundary.neighbor();
              transformationHi = boundary.getTransformation();
              
              // Find any fine boxes on the far side of the boundary.
              // Transform the test box to the remote space
              crFaceBoxTrHi = transformationHi.transformFwd(crFaceBox);

              // This is with respect to the CrFn boxes (following conventions
              // in LevelFluxRegister)
              if (dstProc == procID())
                {
                  DataIndex didxCr(litCr());
                  coarseLocationsHi = 
                    &(getCoarseLocations(dir, Side::Lo)[didxCr]);
                }
            }

          // The painful N^2 loop (but we bail and truncate where possible)
          if (haveCorrection)
            {
              int small0 = std::numeric_limits<int>::max();
              int big0   = std::numeric_limits<int>::min();
              bool haveLo = false;
              bool haveHi = false;
              if (!crFaceBoxTrLo.isEmpty())
                {
                  haveLo = true;
                  small0 = crFaceBoxTrLo.smallEnd(0);
                  big0   = crFaceBoxTrLo.bigEnd(0);
                }
              if (!crFaceBoxTrHi.isEmpty())
                {
                  haveHi = true;
                  small0 = std::min(small0, crFaceBoxTrHi.smallEnd(0));
                  big0   = std::max(big0,   crFaceBoxTrLo.bigEnd(0));
                }
              for (LayoutIterator litCrFn = crFnGrid.layoutIterator();
                   litCrFn.ok(); ++litCrFn)
                {
                  Box crFnBox = crFnGrid[litCrFn];
                  int srcProc = crFnGrid.procID(litCrFn());

                  // No intersection
                  if (crFnBox.bigEnd(0)   < small0) continue;
                  // Due to sorting no other boxes intersect
                  if (crFnBox.smallEnd(0) > big0) break;

                  // Intersections on the low side
                  if (haveLo && crFnBox.intersectsNotEmpty(crFaceBoxTrLo))
                    {
                      // Now in the space of crBox
                      const Box crFnBoxLo =
                        transformationLo.transformBack(crFnBox);
                      // Intersection should be same in both local and remote
                      // index spaces
                      CH_assert(crFnBoxLo.intersectsNotEmpty(crFaceBox));

                      Box intersect = crFnBoxLo & crFaceBox;
                      // Find the cells just inside the original coarse box
                      intersect = adjCellHi(intersect, dir);

                      // Local settings
                      if (dstProc == procID())
                        {
                          CH_assert(coarseLocationsLo);
                          coarseLocationsLo->push_back(intersect);
                        }

                      // Global settings
                      loDstBoxes.push_back(intersect);
                      Box otherSide = transformationLo.transformFwd(intersect);
                      loSrcBoxes.push_back(otherSide);
                      // Map from src (otherSide) to dst (intersect)
                      auto ins = loSrcToDstMap.insert({ otherSide, intersect });
                      CH_assert(ins.second);  //  Must be unique insertion

                      loSrcProcs.push_back(srcProc);
                      loDstProcs.push_back(dstProc);
                    }

                  // Intersections on the high side
                  if (haveHi && crFnBox.intersectsNotEmpty(crFaceBoxTrHi))
                    {
                      // Now in the space of crBox
                      const Box crFnBoxHi =
                        transformationHi.transformBack(crFnBox);
                      // Intersection should be same in both local and remote
                      // index spaces
                      CH_assert(crFnBoxHi.intersectsNotEmpty(crFaceBox));

                      Box intersect = crFnBoxHi & crFaceBox;
                      // Find the cells just inside the original coarse box
                      intersect = adjCellLo(intersect, dir);

                      // Local settings
                      if (dstProc == procID())
                        {
                          CH_assert(coarseLocationsHi);
                          coarseLocationsHi->push_back(intersect);
                        }

                      // Global settings
                      hiDstBoxes.push_back(intersect);
                      Box otherSide = transformationHi.transformFwd(intersect);
                      // Map from src (otherSide) to dst (intersect)
                      auto ins = hiSrcToDstMap.insert({ otherSide, intersect });
                      CH_assert(ins.second);  //  Must be unique insertion
                      hiSrcBoxes.push_back(otherSide);

                      hiSrcProcs.push_back(srcProc);
                      hiDstProcs.push_back(dstProc);
                    }
                }  // Loop over coarsened-fine layout
            }  // If have correction
        }  // Loop over coarse layout

      // Now that we have source/destination boxes for coarsened fine faces
      // in this direction, we create a copier that will traffic the data
      // from one to the other.

      // Low side
      {
        DisjointBoxLayout srcGrid(loSrcBoxes, loSrcProcs);
        DisjointBoxLayout dstGrid(loDstBoxes, loDstProcs);
        m_remoteCopiers[dir][Side::Lo].define(crFnGrid, // a_from
                                              a_CrGrid, // a_to
                                              srcGrid, // src
                                              dstGrid, // dest
                                              loSrcToDstMap,
                                              IntVect::Zero);
      }
      // High side
      {
        DisjointBoxLayout srcGrid(hiSrcBoxes, hiSrcProcs);
        DisjointBoxLayout dstGrid(hiDstBoxes, hiDstProcs);
        m_remoteCopiers[dir][Side::Hi].define(crFnGrid, // a_from
                                              a_CrGrid, // a_to
                                              srcGrid, // src
                                              dstGrid, // dest
                                              hiSrcToDstMap,
                                              IntVect::Zero);
      }
    }  // Loop over directions
}
#endif

/*--------------------------------------------------------------------*/
//  Increment fine
/**
*//*------------------------------------------------------------------*/

void
MultiBlockFaceRegister::incrementFine(const LevelData<FluxBox>& a_fnFlux,
                                      const Real                a_scale,
                                      const Interval&           a_srcIntv,
                                      const Interval&           a_dstIntv)
{
  CH_assert(isDefined());

  RealVect denom(1.0);
  if (m_scaleFineFaces)
    {
      for (const int dir : EachDir)
        {
          IntVect nRef(m_nRefine);
          nRef[dir] = 1;
          denom[dir] = stc::product(nRef);
        }
    }
  // Since we are not implementing a correction, we do not alter the sign based
  // on the side as is done in LevelFluxRegister
  const RealVect scale = a_scale/denom;

  for (DataIterator dit = a_fnFlux.dataIterator(); dit.ok(); ++dit)
    {
      for (const int dir : EachDir)
        {
          for (const auto side : EachSide)
            {
              const Vector<Box>& amrLocations =
                m_mbRegions.getFnLocations(dit, dir, side);
              if (amrLocations.size() > 0)
                {
                  const FluxBox& fnFlux = a_fnFlux[dit];
                  FArrayBox& fnFluxDir = const_cast<FArrayBox&>(fnFlux[dir]);
                  // To outer cells
                  fnFluxDir.shiftHalf(dir, sign(side));

                  // Zero registers
                  FArrayBox& crFnFlux = m_fineFlux[dit];
                  const Box& disjointFnBox = a_fnFlux.disjointBoxLayout()[dit];
                  const Box disjointCrFnBox = coarsen(disjointFnBox, m_nRefine);
                  const Box zeroBox = adjCellBox(disjointCrFnBox, dir, side, 1);
                  crFnFlux.setVal(
                    0., zeroBox, a_dstIntv.begin(), a_dstIntv.size());

                  // Fill registers
                  for (int cDst = a_dstIntv.begin(),
                         cDst_end = a_dstIntv.end() + 1,
                         cSrc = 0; cDst != cDst_end; ++cDst, ++cSrc)
                    {
                      for (const Box& box : amrLocations)
                        {
                          CH_assert(fnFluxDir.contains(box));
                          MD_BOXLOOP(box, i)
                            {
                              const IntVect ivCr =
                                stc::coarsen(MD_GETIV(i), m_nRefine);
                              crFnFlux[MD_IV(ivCr, cDst)] +=
                                scale[dir]*fnFluxDir[MD_IX(i, cSrc)];
                            }
                        }
                    }

                  // Shift back
                  fnFluxDir.shiftHalf(dir, -sign(side));
                }
            }
        }
    }

  // Increment fine was called
  m_isDefined |= FaceRegIncFine;
}

/*--------------------------------------------------------------------*/
//  Operator applied during copyTo
/**
*//*------------------------------------------------------------------*/

class MBAddFaceOp: public LDOperator<FArrayBox>
{
public:

//--Constructors/destructors

#if 0  /* old implementation */
  // Construct
  MBAddFaceOp(const int a_dir,
              const Side::LoHiSide a_side,
              const MultiBlockCoordSys* const a_coordSysPtr,
              const Real a_scale = 1.,
              const bool a_directional = true)
    :
    m_dir(a_dir),
    m_side(a_side),
    m_coordSysPtr(a_coordSysPtr),
    m_directional(a_directional),
    m_scale(a_scale)
  {  }
#endif

  /// Constructor
  MBAddFaceOp(LevelData<FluxBox>&             a_dstLvlFxb,
              const MultiBlockCoordSys* const a_coordSysPtr,
              const Real                      a_scale = 1.)
    :
    m_dstLvlFxb(a_dstLvlFxb),
    m_coordSysPtr(a_coordSysPtr),
    m_scale(a_scale)
    { }

  // Use synthesized destructor, copy, and assignment

//--Member functions

  /// Use API 1 (which gives the dataIndex for the destination)
  virtual int API() const override
    { return 1; }

/*
  There's a lot of weird stuff in the old implementation.  E.g., I'm not 
  convinced both vector and tensor transforms are required.  Right now, the
  class is optimized for the volume flux which is a scalar flux (see linearIn1
  and op1).  If scriptN is not used, this class might be used to replace N on
  coarser meshes with that from finer meshes.  If you want to use that, you have
  to reimplement the functionality for copying N (sorry).  The old
  implementation is here in linearIn, op, and backTransform.
*/

#if 0  /* old implementation */
  /// Operation for distributed memory (src is the buffer)
  virtual void linearIn(FArrayBox& arg,
                        void* buf,
                        const Box& regionTo,
                        const Interval& comps) const
    {
      Real* buffer = (Real*)buf;

      IndicesTransformation itBack = backTransform(arg);
      Box regionFrom = itBack.transformFwd(regionTo);

      FArrayBox tmp(regionFrom, comps.size(), buffer);
      op(arg, regionFrom, comps, regionTo, tmp, comps);
    }


  /// Operation for local memory
  void op(FArrayBox&       a_dest,
          const Box&       a_regionFrom,
          const Interval&  a_Cdest,
          const Box&       a_regionTo,
          const FArrayBox& a_src,
          const Interval&  a_Csrc) const
    {
      CH_assert(a_regionFrom.numPts() == a_regionTo.numPts());
      // Shift to align the face with cells.
      // If low side, then shift each m_dir-face up to the cell above.
      // If high side, then shift each m_dir-face down to the cell below.
      // We'll shift back when we're done.
      a_dest.shiftHalf(m_dir, -sign(m_side));

      CH_assert(a_dest.box().contains(a_regionTo));
      CH_assert(a_src.box().contains(a_regionFrom));

      IndicesTransformation itBack = backTransform(a_dest);

      // Get an understanding of sign transformations in the 'to space'.  Signs
      // are not considered if 'm_directional' is set to false
      const IntVect signTrfm = (m_directional) ?
        itBack.transformVectorBack(IntVect::Unit) : IntVect::Unit;

      CH_assert(a_Csrc.size() == a_Cdest.size());
      int ncomp = a_Csrc.size();
      if (ncomp % SpaceDim == 0)
        {

//--Not only do we have to transform the to/from locations, but the components
//--themselves are assumed to represent vector data and must be transformed.

          const int compPerDir = ncomp/SpaceDim;
          IntVect compIVDest0;
          for (const int dir : EachDir)
            {
              compIVDest0[dir] = dir*compPerDir;
            }
          IntVect compIVSrc0 = itBack.transformVectorFwd(compIVDest0);

          if (compPerDir == SpaceDim)
            {

//--Tensor, meaning we also have to transform each range of components

              IntVect compIVDest1{ 0, 1, 2, 3, 4, 5 };
              IntVect compIVSrc1 = itBack.transformVectorFwd(compIVDest1);

              // For each direction of components
              for (const int dir : EachDir)
                {
                  // Over all components in a given direction
                  for (int comp = 0; comp != compPerDir; ++comp)
                    {
                      const int compCdest =
                        a_Cdest.begin() + compIVDest0[dir] + compIVDest1[comp];
                      // Sign changes taken from the inner vector transform
                      const Real scale = signTrfm[comp]*m_scale;
                      const int compCsrc = a_Csrc.begin() +
                        abs(compIVSrc0[dir]) + abs(compIVSrc1[comp]);
                      MD_BOXLOOP(a_regionTo, iTo)
                        {
                          IntVect ivFrom = itBack.transformFwd(MD_GETIV(iTo));
                          a_dest[MD_IX(iTo, compCdest)] +=
                            scale * a_src[MD_IV(ivFrom, compCsrc)];
                        }
                    }
                }
            }
          else
            {

//--Vector only (no additional vector meaning within compPerDir)

              // I'm not sure what the use case is nor if the layout of
              // components is being interpreted correctly.  This code is okay
              // only if you have compPerDir components and we are rearranging
              // those blocks.
              CH_assert(false);

              // For each direction of components
              for (const int dir : EachDir)
                {
                  // Over all components in a given direction
                  int compCdest = a_Cdest.begin() + compIVDest0[dir];
                  // Sign changes
                  const Real scale = signTrfm[dir]*m_scale;
                  int compCsrc = a_Csrc.begin() + abs(compIVSrc0[dir]);
                  for (int comp = 0; comp != compPerDir; ++comp)
                    {
                      MD_BOXLOOP(a_regionTo, iTo)
                      for (BoxIterator bit(a_regionTo); bit.ok(); ++bit)
                        {
                          IntVect ivFrom = itBack.transformFwd(MD_GETIV(iTo));
                          a_dest[MD_IX(iTo, compCdest)] +=
                            scale * a_src[MD_IV(ivFrom, compCsrc)];
                        }
                      ++compCsrc;
                      ++compCdest;
                    }
                }
            }
        }
      else
        {

//--The components are assumed to be scalar fluxes associated with 'm_dir'.
//--We can directly copy from 'from' direction to 'to direction' and only have
//--to worry about the sign change associated with 'm_dir'.
//--Locations are transformed.

          int compCsrc = a_Csrc.begin();
          int compCdest = a_Cdest.begin();
          // Sign changes
          const Real scale = signTrfm[m_dir]*m_scale;
          for (int comp = 0; comp < ncomp; comp++)
            {
              MD_BOXLOOP(a_regionTo, iTo)
                {
                  IntVect ivFrom = itBack.transformFwd(MD_GETIV(iTo));
                  a_dest[MD_IX(iTo, compCdest)] +=
                    scale * a_src[MD_IV(ivFrom, compCsrc)];
                }
              ++compCsrc;
              ++compCdest;
            }
        }

      // Shift back to original
      a_dest.shiftHalf(m_dir, sign(m_side));
    }
#endif

  /// Operation for distributed memory
  virtual void linearIn1(FArrayBox&       a_dstFab,
                         void*            a_srcBuf,
                         const Box&       a_dstRegion,
                         const Interval&  a_comps,
                         const DataIndex& a_dstDidx) const override
    {
      Real* buffer = static_cast<Real*>(a_srcBuf);

      MultiBlockRegions::FaceTag faceTag(-1, a_dstDidx);
      IndicesTransformation trfm = getTransform(a_dstRegion, faceTag);
      const int dir = faceTag.dir();
      const Side::LoHiSide side = faceTag.side();

      // The copy is from outer cells to outer cells.  To find outer on the
      // src side, we transform inner on the destination side.
      const int shiftIn = -sign(side);
      Box inner(a_dstRegion);
      inner.shift(dir, shiftIn);
      const Box srcRegion = trfm.transformFwd(inner);

      FArrayBox tmp(srcRegion, a_comps.size(), buffer);
      trfmOp(a_dstDidx, a_comps, a_dstRegion, tmp, a_comps, dir, side, trfm);
    }

  /// Operation for shared memory
  virtual void op1(FArrayBox&       a_dstFab,  // Unused because of redirect
                   const Box&       a_srcRegion,
                   const Interval&  a_dstComps,
                   const Box&       a_dstRegion,
                   const FArrayBox& a_srcFab,
                   const Interval&  a_srcComps,
                   const DataIndex& a_dstDidx) const override
    {
      CH_assert(a_srcRegion.numPts() == a_dstRegion.numPts());
      CH_assert(a_srcComps.size()    == a_dstComps.size());

      MultiBlockRegions::FaceTag faceTag(-1, a_dstDidx);
      IndicesTransformation trfm = getTransform(a_dstRegion, faceTag);
#ifndef NDEBUG
      Box testRegion(a_dstRegion);
      testRegion.shift(faceTag.dir(), -sign(faceTag.side()));
      CH_assert(trfm.transformFwd(testRegion) == a_srcRegion);
#endif

      trfmOp(a_dstDidx, a_dstComps, a_dstRegion, a_srcFab, a_srcComps,
             faceTag.dir(), faceTag.side(), trfm);
    }

private:

#if 0  /* old implementation */
  IndicesTransformation backTransform(const FArrayBox&  a_destFab) const
  {
    // Call whichBlockOverlap() on the FABs' boxes,
    // instead of whichBlock() on a_regionFrom and a_regionTo,
    // because we're going to have cells OUTSIDE valid blocks.
    // int blockFrom = m_coordSysPtr->whichBlockOverlap(a_src.box());
    int blockTo = m_coordSysPtr->whichBlockOverlap(a_destFab.box());
    // petermc, 14 Nov 2011: Do not obtain the IndicesTransformation
    // this way, because there may be more than one way to get from
    // source block to dest block!
    // const IndicesTransformation& it =
    // m_coordSysPtr->blockTransformation(blockFrom, blockTo);
    const BlockBoundary& bb = m_coordSysPtr->boundary(blockTo, m_dir, m_side);
    IndicesTransformation itBack = bb.getTransformation();
    return itBack;
  }
#endif

  /// The op for copying from a src to dst
  /** Destinations are redirected to m_dstLvlFxb
   */
  void trfmOp(const DataIndex&             a_dstDidx,
              const Interval&              a_dstComps,
              const Box&                   a_dstRegion,
              const FArrayBox&             a_srcFab,
              const Interval&              a_srcComps,
              const int                    a_dir,
              const Side::LoHiSide         a_side,
              const IndicesTransformation& a_trfm) const
    {
      // Grab the fab and shift.  Fine cell values are in layers of cells
      // outside the coarse box.  Store to faces.
      FArrayBox& dstFab = m_dstLvlFxb[a_dstDidx][a_dir];
      dstFab.shiftHalf(a_dir, sign(a_side));
      CH_assert(dstFab.contains(a_dstRegion));

      // Directions are in the new index space so we need the inverse transform
      // to query a_dir in destination space.
      const IndicesTransformation invTr = a_trfm.inverse();
      const int fluxSign = invTr.getSign()[a_dir];

      const Real scale = fluxSign*m_scale;

      // Box srcRegion = trmf.transformFwd(
      // dumpFAB2DSlicePretty(&a_srcFab, 0, a_trmf.transformFws

      const int shiftIn = -sign(a_side);
      for (int cDst = a_dstComps.begin(), cDst_end = a_dstComps.end() + 1,
             cSrc = a_srcComps.begin(); cDst != cDst_end; ++cDst, ++cSrc)
        {
          MD_BOXLOOP(a_dstRegion, i)
            {
              // a_dstRegion is the outer cells, the inner cells, when
              // transformed, match the source region.
              const IntVect ivDstIn = stc::shift(MD_GETIV(i), a_dir, shiftIn);
              const IntVect ivSrc = a_trfm.transformFwd(ivDstIn);
              dstFab[MD_IX(i, cDst)] += scale*a_srcFab[MD_IV(ivSrc, cSrc)];
            }
        }

      // Shift back to original
      dstFab.shiftHalf(a_dir, -sign(a_side));
    }

  /// Get the transformation for this copy
  /** \param[in]  a_dstDidx
   *                      The dataIndex for data we are copying into
   *  \param[out] a_faceTag
   *                      The face we are working on
   */
  IndicesTransformation getTransform(
    const Box&                  a_dstRegion,
    MultiBlockRegions::FaceTag& a_faceTag) const
    {
      int dstIdxBlk =
        m_dstLvlFxb.disjointBoxLayout().blockIndex(a_faceTag.m_srcDidx);
      CH_assert(dstIdxBlk >= 0);
      Box block = m_coordSysPtr->mappingBlocks()[dstIdxBlk];
      // Query the face.  The destination region will have dimension size 1 in
      // this direction and is outside the mapping block.
      a_faceTag.m_idxFace = -1;
      for (const int dir : EachDir)
        {
          const int idxDir = a_dstRegion.smallEnd(dir);
          if (idxDir == a_dstRegion.bigEnd(dir))
            {
              for (const auto side : EachSide)
                {
                  if (idxDir == (block.sideEnd(side)[dir] + sign(side)))
                    {
                      a_faceTag.m_idxFace =
                        MultiBlockRegions::FaceTag::indexFace(dir, side);
                      goto haveFace;
                    }
                }
            }
        }
      haveFace: ;
      CH_assert(a_faceTag.m_idxFace != -1);
      const BlockBoundary& bb =
        m_coordSysPtr->boundary(dstIdxBlk, a_faceTag.dir(), a_faceTag.side());
      return bb.getTransformation();
    }

//--Data members

#if 0  /* old implementation */
  int m_dir;                          ///< The dir this operation affects
  Side::LoHiSide m_side;              ///< The side this operation affects
  const MultiBlockCoordSys* m_coordSysPtr;
  bool m_directional;                 ///< T - Vector data on faces has a
                                      ///<     direction (normal & default)
                                      ///< F - The vector data does not have
                                      ///<     a direction.  The best example of
                                      ///<     this is <N> which has vector
                                      ///<     components but really just scales
                                      ///<     areas.

public:
  Real m_scale;                       ///< Extra scaling beyond what is normally
                                      ///< needed to integrated the fine fluxes
#endif

  LevelData<FluxBox>& m_dstLvlFxb;    ///< CopyTo is hijacked and writes
                                      ///< redirected to here.  This permits
                                      ///< copying from a LevelData<FArrayBox>
                                      ///< to a LevelData<FluxBox>.  The
                                      ///< original destination provided as an
                                      ///< argument to copyTo is unused.
  const MultiBlockCoordSys* m_coordSysPtr;
                                      ///< Coordinate system
  Real m_scale;                       ///< Extra scaling beyond what is normally
                                      ///< needed to integrated the fine fluxes
  
};

/*--------------------------------------------------------------------*/
//  Overwrite solution on the coarse faces
/** \param[out] a_CrF   Coarse information updated with face
 *                      corrections from the finer level
 *  \param[in]  a_CrFIntv
 *                      Interval in a_CrF to update
 *  \param[in]  a_faceIntv
 *                      Interval in interior structures
 *  \param[in]  a_scale Additional scaling
 *  \param[in]  a_directional
 *                      T - Vector data on faces has a direction
 *                          (normal & default)
 *                      F - The vector data does not have a direction.
 *                          The best example of this is <N> which has
 *                          vector components but really just scales
 *                          areas.
 *  \param[in]  a_N    Metrics on the faces.  Define if a_CrF is a
 *                     volume flux and we need to jump across periodic
 *                     boundaries.
 *  \param[in]  a_dx   Mesh spacing for the level.  Only needed if a_N
 *                     is defined.
*//*------------------------------------------------------------------*/

void
MultiBlockFaceRegister::reface(LevelData<FluxBox>&       a_CrF,
                               const Interval&           a_CrFIntv,
                               const Interval&           a_faceIntv,
                               const Real                a_scale,
                               const bool                a_directional,
                               const LevelData<FluxBox>* a_N,
                               const RealVect&           a_dx)
{
  // Note: Some of these variables could be named better.  CrF is not
  //       coarsened-fine but coarse flux.
  CH_TIME("MultiBlockFaceRegister::reface");
  CH_assert(m_isDefined & FaceRegIncFine);
  // Have to call incrementFine before next access to reface
  m_isDefined &= ~FaceRegIncFine;

  // Initialize coarse locations to zero, because the flux update following is
  // an additive operation that needs to replace the current values
  for (DataIterator dit = a_CrF.dataIterator(); dit.ok(); ++dit)
    {
      const Box& disjointBox = a_CrF.disjointBoxLayout()[dit];
      FluxBox& CrFFlbx = a_CrF[dit];
      const int blockID = m_coordSys->whichBlock(disjointBox);
      const NewFourthOrderCoordSys* blockCoordSys =
        static_cast<const NewFourthOrderCoordSys*>(
          m_coordSys->getCoordSys(blockID));
      // Side of the coarse box
      for (const auto side : EachSide)
        {
          const int sideSign = sign(side);
          for (const int dir : EachDir)
            {
              const Vector<Box>& intersect =
                m_mbRegions.getCrLocations(dit, dir, side);
              if (intersect.size() > 0)
                {
                  FArrayBox& CrFFabDir = CrFFlbx[dir];
                  if (a_N)
                    {
                      CH_assert(a_CrFIntv.size() == 1);
                      // Cells just inside the boundary
                      Box deltaBox =
                        adjCellBox(disjointBox, dir, side, -1);
                      // Shift onto faces
                      deltaBox.shiftHalf(dir, sideSign);
                      FArrayBox delta(deltaBox, SpaceDim);
                      const FArrayBox& Ndir = a_N->operator[](dit)[dir];
                      // Offset to get a Cart. Coord. from an IntVect
                      const RealVect offset(
                        0.5*(RealVect_unit - RealVect_basis(dir)));
                      RigidTransformation displ =
                        m_coordSys->boundary(blockID, dir, side)
                        .getPhysTransformation();
                  
                      for (int b = 0, bEnd = intersect.size(); b != bEnd; ++b)
                        {
                          Box faceIntersect(intersect[b]);
                          faceIntersect.shiftHalf(dir, sideSign);
                          MD_BOXLOOP(faceIntersect, i)
                            {
                              // We assume the diplacement is not a constant so
                              // we can also handle axisymmetric periodicity
                              IntVect iv = MD_GETIV(i);
                              const RealVect cartLoc = a_dx*(iv + offset);
                              const RealVect physLoc =
                                blockCoordSys->realCoord(cartLoc);
                              // Need displacement from CrFn to Cr.  Displace X
                              // from Cr to CrFn and then negate.
                              const RealVect deltaPnt =
                                (physLoc - displ.transformFwd(physLoc))/
                                SpaceDim;
                              for (const int c : EachDir)
                                {
                                  delta[MD_IX(i, c)] = deltaPnt[c];
                                }
                            }
                          // Dot product of row of N^T and delta
                          int cBeg = 0;
                          int cNum = 1;
                          int cBegN = blockCoordSys->getNcomponent(0, dir);
                          CH_assert(CrFFabDir.box().contains(faceIntersect));
                          CH_assert(delta.box().contains(faceIntersect));
                          CH_assert(Ndir.box().contains(faceIntersect));
                          FORT_POINTFDOTN(CHF_FRA(CrFFabDir),
                                          CHF_CONST_INT(cBeg),
                                          CHF_CONST_INT(cNum),
                                          CHF_CONST_FRA(delta),
                                          CHF_CONST_INT(cBeg),
                                          CHF_CONST_INT(cNum),
                                          CHF_CONST_INT(cNum),
                                          CHF_CONST_FRA(Ndir),
                                          CHF_CONST_INT(cBegN),
                                          CHF_CONST_INT(SpaceDim),
                                          CHF_BOX(faceIntersect));
                        }
                    }
                  else
                    {
                      for (int b = 0, bEnd = intersect.size(); b != bEnd; ++b)
                        {
                          Box faceIntersect(intersect[b]);
                          faceIntersect.shiftHalf(dir, sideSign);
                          CrFFabDir.setVal(0.,
                                           faceIntersect,
                                           a_CrFIntv.begin(),
                                           a_CrFIntv.size());
                        }
                    }
                }
            }
        }
    }

  /*
    Our source information is in FArrayBoxes.  The destination is ultimately a
    FluxBox.  You cannot copyTo between a FAB and FxB, so instead we fake the
    destination as an alias.  As long as a_CrF has at least one layer of ghost
    cells, we are good to go in terms of valid regions since we are copying into
    the outer cells.  The actual copy is redirected in MBAddFaceOp.  This may
    work even if there are zero ghost cells as only an assert not aware of the
    redirect may cause problems.
  */

  {
    LevelData<FArrayBox> CrFAlias(
      a_CrF.getBoxes(),
      a_CrF.nComp(),
      IntVect_unit,
      FABAliasFlBxDataFactory(&a_CrF, a_CrF.interval(), 0));
    MBAddFaceOp op(a_CrF, m_coordSys, a_scale);
    m_fineFlux.copyTo(a_faceIntv,
                      CrFAlias,
                      a_CrFIntv,
                      m_mbRegions.amrFaceRemoteCopier(),
                      op);
  }

#if 0
  // Copy all the fine faces from their "native" locations to those abutting
  // the coarse faces.
  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      // Alias a direction of the flux box with the original interval
      LevelData<FArrayBox> CrFDirAlias(
        a_CrF.getBoxes(),
        a_CrF.nComp(),
        IntVect_unit,
        FABAliasFlBxDataFactory(&a_CrF, a_CrF.interval(), dir));

      for (const auto side : EachSide)
        {
          if (m_remoteCopiers[dir][side].isDefined())
            {
              // Copy the coarsened fine faces to the corresponding coarse faces
              // on their neighboring blocks.
              MBAddFaceOp op(dir, side, m_coordSys, -a_scale, a_directional);
              m_fineFlux.copyTo(a_faceIntv,
                                CrFDirAlias,
                                a_CrFIntv,
                                m_remoteCopiers[dir][side],
                                op);
            }
        } // end iteration over side
    } // end iteration over dir
#endif

}

#include "NamespaceFooter.H"
