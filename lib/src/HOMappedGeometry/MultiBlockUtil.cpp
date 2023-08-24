#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <stack>

#include "LoHiSide.H"
#include "MultiBlockUtil.H"
#include "BoxIterator.H"
#include "MBStencilIterator.H"
#include "MBVectorStencilIterator.H"
#include "BoxCollapser.H"
#include "BoxFixedOff.H"
#include "MBMiscUtil.H"
#include "FourthOrderUtil.H"
#include "LAPACKMatrix.H"
#include "MultiBlockLevelGeom.H"

#include "NamespaceHeader.H"

/// constructor
MultiBlockUtil::MultiBlockUtil()
  :
  m_isDefined(false)
{}

/// destructor
MultiBlockUtil::~MultiBlockUtil()
{
  undefine();
}

/// full constructor
MultiBlockUtil::MultiBlockUtil(const MultiBlockCoordSys*  a_coordSysPtr,
                               Interval                   a_fixedDims,
                               Vector<int>                a_fixedPt)
  :
  m_isDefined(false)
{
  define(a_coordSysPtr, a_fixedDims, a_fixedPt);
}

void 
MultiBlockUtil::undefine()
{
  if (m_isDefined)
    {
      m_interpDimsVect.clear();
      m_fixedDimsVect.clear();
      for (BoxIterator bitOff(m_offsetAllBox); bitOff.ok(); ++bitOff)
        {
          IntVect ivOff = bitOff();
          for (int srcBlock = 0; srcBlock < m_nblocks; srcBlock++)
            {
              delete m_transformsAll(ivOff, srcBlock);
            }
        }
    }
  m_isDefined = false;
}

void 
MultiBlockUtil::define(const MultiBlockCoordSys*  a_coordSysPtr,
                       Interval                   a_fixedDims,
                       Vector<int>                a_fixedPt)
{ CH_TIME("MultiBlockUtil::define");
  undefine();
  m_coordSysPtr = (MultiBlockCoordSys*) a_coordSysPtr;

  m_mappingBlocks = m_coordSysPtr->mappingBlocks();
  m_nblocks = m_mappingBlocks.size();

#if 0
  m_boundaries = m_coordSysPtr->boundaries();
#else
  // workaround until we figure out why the above breaks COGENT
  const Vector < Tuple< BlockBoundary, 2*SpaceDim > >& src_boundaries =
    m_coordSysPtr->boundaries();
  m_boundaries.resize(m_nblocks);
  for (int i=0; i<m_nblocks; ++i)
    {
      m_boundaries[i] = src_boundaries[i];
    }
#endif

  m_fixedDims = a_fixedDims;
  m_interpDimsVect.clear();
  m_fixedDimsVect.clear();
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if ( a_fixedDims.contains(idir) )
        { // fix dimension idir
          m_interpUnit[idir] = 0;
          m_fixedDimsVect.push_back(idir);
        }
      else
        { // interpolate in dimension idir
          m_interpUnit[idir] = 1;
          m_interpDimsVect.push_back(idir);
        }
    }
  m_fixedPt = Vector<int>(a_fixedPt);

  // Use these in displacementPowers(), in cell-averaged case.
  int nInterpDims = m_interpDimsVect.size();
  m_convSize = 1+2*nInterpDims;
  m_convBaseCells.resize(m_convSize);
  m_convBaseCells.assign(IntVect::Zero);
  // Shifts are [0,0], [-1,0], [0,-1], [+1,0], [0,+1].
  // WAS for (int idir = 0; idir < SpaceDim; idir++)
  for (int ind = 0; ind < nInterpDims; ind++)
    {
      int idir = m_interpDimsVect[ind];
      m_convBaseCells[1 + ind][idir] -= 1;
      m_convBaseCells[1 + nInterpDims + ind][idir] += 1;
    }
  Real convScaling = 1./24.;
  m_convWeight.resize(m_convSize);
  m_convWeight.assign(convScaling);
  m_convWeight[0] = 1. - (2*nInterpDims)*convScaling;

  // Box m_offsetAllBox;
  m_offsetAllBox = Box(-IntVect::Unit, IntVect::Unit);

  //  m_transformsAll.resize(m_nblocks);
  m_transformsAll.define(m_offsetAllBox, m_nblocks);

  for (BoxIterator bitOff(m_offsetAllBox); bitOff.ok(); ++bitOff)
    {
      IntVect ivOff = bitOff();
      for (int srcBlock = 0; srcBlock < m_nblocks; srcBlock++)
        {
          m_transformsAll(ivOff, srcBlock) =
            new Vector<IndicesTransformation>(m_nblocks);
        }
    }

  for (int dstBlock = 0; dstBlock < m_nblocks; dstBlock++)
    {
      /*
        For each offset IntVect ivOff,
        validBlocks(ivOff, 0:validNum(ivOff, 0)-1)
        will hold the block numbers of valid blocks containing
        ghost cells in direction ivOff from block dstBlock,
        and validTransformations(ivOff, 0:validNum(ivOff, 0)-1)
        will hold the transformations from dstBlock to valid blocks
        in the same order.
      */
      const BaseFab<int>& validBlocks =
        m_coordSysPtr->validBlocks(dstBlock);
      const BaseFab<IndicesTransformation>& validTransformations =
        m_coordSysPtr->validTransformations(dstBlock);
      // validNum gives number of components of
      // validBlocks and validTransformations.
      const BaseFab<int>& validNum =
        m_coordSysPtr->validNum(dstBlock);
      for (BoxIterator bitOff(m_offsetAllBox); bitOff.ok(); ++bitOff)
        {
          IntVect ivOff = bitOff();
          int validNumOff = validNum(ivOff, 0);
          for (int i = 0; i < validNumOff; i++)
            {
              int srcBlock = validBlocks(ivOff, i);

              // Transformation from dstBlock to srcBlock.
              IndicesTransformation tfm = validTransformations(ivOff, i);
              // (*m_transformsAll[srcBlock])(ivOff, dstBlock) = tfm;
              (*m_transformsAll(ivOff, srcBlock))[dstBlock] = tfm;
            }
        }
    }

  m_isDefined = true;
}


IntVect
MultiBlockUtil::whichOffset(int             a_block,
                            const IntVect&  a_cell) const
{
  CH_assert(m_isDefined);
  const Box& blockBox = m_mappingBlocks[a_block];
  IntVect offset = IntVect::Zero;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (a_cell[idir] < blockBox.smallEnd(idir))
        offset[idir] = -1;
      else if (a_cell[idir] > blockBox.bigEnd(idir))
        offset[idir] = +1;
    }
  return offset;
}


const IndicesTransformation&
MultiBlockUtil::blockTransform(int a_srcBlock,
                               int a_dstBlock,
                               const IntVect& a_offset) const
{
  CH_assert(m_isDefined);
  const Vector<IndicesTransformation>& tfmSrcOff =
    transformsSrcOff(a_srcBlock, a_offset);
  return tfmSrcOff[a_dstBlock];
}


const Vector<IndicesTransformation>&
MultiBlockUtil::transformsSrcOff(int a_srcBlock,
                                 const IntVect& a_offset) const
{
  CH_assert(m_isDefined);
  const Vector<IndicesTransformation>& tfmSrcOff =
    *m_transformsAll(a_offset, a_srcBlock);
  return tfmSrcOff;
}


void
MultiBlockUtil::growIVS(IntVectSet& a_ivs,
                        int a_bufferSize) const
{
  CH_assert(m_isDefined);
  CH_TIME("MultiBlockUtil::growIVS");
  IntVectSet ivsGrown(a_ivs);
  ivsGrown.grow(a_bufferSize);
  for (IVSIterator ivsit(ivsGrown); ivsit.ok(); ++ivsit)
    {
      IntVect iv = ivsit();
      if ( !a_ivs.contains(iv) )
        {
          int blockValid = m_coordSysPtr->whichBlock(iv);
          if (blockValid != -1)
            { // iv is a valid cell of blockValid
              a_ivs |= iv;
            }
          else
            { // iv is not a valid cell of any block
              int dstBlock = m_coordSysPtr->whichBlockBuffered(iv, a_bufferSize);
              IntVect offset = whichOffset(dstBlock, iv);
              CH_assert(offset != IntVect::Zero);
              for (int srcBlock = 0; srcBlock < m_nblocks; srcBlock++)
                {
                  // tfm transforms from dstBlock to srcBlock,
                  // taking offset from dstBlock.
                  const IndicesTransformation& tfm = blockTransform(srcBlock,
                                                                    dstBlock,
                                                                    offset);
                  if ( tfm.isDefined() )
                    {
                      IntVect ivValid = tfm.transform(iv);
                      a_ivs |= ivValid;
                    }
                }
            }
        }
    }
}

void
MultiBlockUtil::commonVertexCells(Vector<IntVect>& a_neighborCells,
                                  Vector<int>& a_neighborBlocks,
                                  const IntVect& a_baseCell,
                                  int a_baseBlock) const
{
  CH_assert(m_isDefined);
  CH_TIME("MultiBlockUtil::commonVertexCells");
  Box baseBlockBox = m_mappingBlocks[a_baseBlock];

  // Check for the easy case: a_baseCell is away from all block boundaries.
  Box innerBlockBox = grow(baseBlockBox, -1);
  if (innerBlockBox.contains(a_baseCell))
    {
      a_neighborCells.clear();
      a_neighborBlocks.clear();
      Box neighborsBox(a_baseCell - IntVect::Unit,
                       a_baseCell + IntVect::Unit);
      for (BoxIterator bitNbr(neighborsBox); bitNbr.ok(); ++bitNbr)
        {
          IntVect nbr = bitNbr();
          a_neighborCells.push_back(nbr);
          a_neighborBlocks.push_back(a_baseBlock);
        }
    }
  else
    {
      Box offsetBox(IntVect::Zero, IntVect::Unit);
      Vector<IntVect> allVertices;
      Vector<int> allVertexBlocks;
      for (BoxIterator bitOff(offsetBox); bitOff.ok(); ++bitOff)
        {
          IntVect offset = bitOff();
          IntVect vertex = a_baseCell + offset;
          // Note that as a NODE, vertex is still valid within a_baseBlock.
          // First, find all aliases of vertex.
          Vector<IntVect> aliasVertices;
          Vector<int> aliasVertexBlocks;
          vertexAliases(aliasVertices, aliasVertexBlocks,
                        vertex, a_baseBlock);
          allVertices.append(aliasVertices);
          allVertexBlocks.append(aliasVertexBlocks);
        }
      validCellsFromVertices(a_neighborCells, a_neighborBlocks,
                             allVertices, allVertexBlocks);
    }
  return;
}

void
MultiBlockUtil::commonVertexStencilElements(
  StencilMap&                a_neighbors,
  const MBStencilElement&    a_base,
  const RigidTransformation& a_baseExtraDispl) const
{
  CH_assert(m_isDefined);
  CH_TIME("MultiBlockUtil::commonVertexCells");
  int baseBlock = a_base.block();
  const IntVect& baseCell = a_base.cell();
  Box baseBlockBox = m_mappingBlocks[baseBlock];

  Box innerBlockBox = grow(baseBlockBox, -1);
  if (innerBlockBox.contains(baseCell))
    { // easy case: baseCell is away from all block boundaries.
      Box neighborsBox(baseCell - IntVect_unit,
                       baseCell + IntVect_unit);
      MD_BOXLOOP(neighborsBox, nbr)
        {
          MBStencilElement element(MD_GETIV(nbr), baseBlock);
          a_neighbors.insert({ element, a_baseExtraDispl });
        }
    }
  else
    { // harder case: baseCell is adjacent to a block boundary.
      Box shiftCellToVertexBox(IntVect_zero, IntVect_unit);
      // Fill allVertexElements with all vertex stencil elements that are
      // vertices of cell stencil element a_base.
      StencilMap allVertexElements;
      MD_BOXLOOP(shiftCellToVertexBox, i)
        {
          const IntVect vertex = baseCell + MD_GETIV(i);
          MBStencilElement elt(vertex, baseBlock);
          // Note that as a NODE, vertex is still valid within baseBlock.
          // First, find all aliases of vertex.
          vertexStencilElementAliases(allVertexElements, elt, a_baseExtraDispl);
        }
      validStencilElementsFromVertices(a_neighbors, allVertexElements);
    }
  return;
}

void
MultiBlockUtil::vertexAliases(Vector<IntVect>& a_aliasVertices,
                              Vector<int>& a_aliasBlocks,
                              const IntVect& a_baseVertex,
                              int a_baseBlock) const
{
  // It might be better to have a struct consisting of IntVect and int.
  List<IntVect> aliasVerticesList;
  List<int> aliasBlocksList;
  List<IntVect> toVisitVertices;
  List<int> toVisitBlocks;
  toVisitVertices.add(a_baseVertex);
  toVisitBlocks.add(a_baseBlock);
  while (toVisitVertices.isNotEmpty())
    {
      // pop this vertex from toVisit list
      IntVect thisVertex = toVisitVertices.firstElement();
      int thisBlock = toVisitBlocks.firstElement();
      toVisitVertices.removeFirst();
      toVisitBlocks.removeFirst();

      // add this vertex to aliases
      aliasVerticesList.add(thisVertex);
      aliasBlocksList.add(thisBlock);

      Box thisBlockNodesBox = surroundingNodes(m_mappingBlocks[thisBlock]);
      int faceID = 0;
      for (SideIterator sit; sit.ok(); sit.next())
        {
          Side::LoHiSide side = sit();
          IntVect sideBlockEnd = thisBlockNodesBox.sideEnd(side);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              if (thisVertex[idir] == sideBlockEnd[idir])
                {
                  const BlockBoundary& bb = m_boundaries[thisBlock][faceID];
                  if (bb.isConformal() || bb.isMapped())
                    {
                      IndicesTransformation tfm = bb.getTransformation();
                      int nbrBlock = bb.neighbor();
                      IntVect nbrVertex = tfm.transformNode(thisVertex);
                      if ( !aliasVerticesList.includes(nbrVertex) )
                        if ( !toVisitVertices.includes(nbrVertex) )
                          {
                            toVisitVertices.add(nbrVertex);
                            toVisitBlocks.add(nbrBlock);
                          }
                    }
                }
              faceID++;
            }
        }
    }

  a_aliasVertices.clear();
  a_aliasBlocks.clear();
  while (aliasVerticesList.isNotEmpty())
    {
      a_aliasVertices.push_back(aliasVerticesList.firstElement());
      a_aliasBlocks.push_back(aliasBlocksList.firstElement());
      aliasVerticesList.removeFirst();
      aliasBlocksList.removeFirst();
    }
  return;
}

void
MultiBlockUtil::vertexStencilElementAliases(
  StencilMap&                a_aliases,
  const MBStencilElement&    a_base,
  const RigidTransformation& a_baseExtraDispl) const
{
  std::stack<StencilMap::iterator> toVisit;

  // Insert the base into aliasList and add iterator to toVisit
  toVisit.push(a_aliases.insert({a_base, a_baseExtraDispl}).first);
  while (!toVisit.empty())
    {
      // Pop this stencil element from toVisit list
      const MBStencilElement& thisElement = toVisit.top()->first;
      const RigidTransformation& thisExtraDispl = toVisit.top()->second;
      toVisit.pop();

      const int thisBlock      = thisElement.block();
      const IntVect thisVertex = thisElement.cell();

      Box thisBlockNodesBox = surroundingNodes(m_mappingBlocks[thisBlock]);
      for (const auto side : EachSide)
        {
          IntVect sideBlockEnd = thisBlockNodesBox.sideEnd(side);
          for (const int dir : EachDir)
            {
              if (thisVertex[dir] == sideBlockEnd[dir])
                {
                  const BlockBoundary& bb = m_coordSysPtr->boundary(thisBlock,
                                                                    dir,
                                                                    side);
                  if (bb.isConformal() || bb.isMapped())
                    {
                      IndicesTransformation tfm = bb.getTransformation();
                      int nbrBlock = bb.neighbor();
                      IntVect nbrVertex = tfm.transformNode(thisVertex);
                      MBStencilElement nbr(nbrVertex, nbrBlock);
                      auto insert = a_aliases.insert(
                        { nbr, thisExtraDispl + bb.getPhysTransformation() });
                      if (insert.second)  // Means insertion succeeded because
                        {                 // key was unique
                          // Save the iterator to the insertion for future visit
                          toVisit.push(insert.first);
                        }
                    }
                }
            }  // Loop over directions
        }  // Loop over sides
    }
}

void
MultiBlockUtil::validCellsFromVertices(Vector<IntVect>& a_cellIndices,
                                       Vector<int>& a_cellBlocks,
                                       const Vector<IntVect>& a_vertexIndices,
                                       const Vector<int>& a_vertexBlocks) const
{
  Box nodeToCellBox(-IntVect::Unit, IntVect::Zero);
  List<IntVect> cellIndicesList;
  List<int> cellBlocksList;
  for (int i = 0; i < a_vertexIndices.size(); i++)
    {
      IntVect thisVertex = a_vertexIndices[i];
      int thisBlock = a_vertexBlocks[i];
      for (BoxIterator bitOff(nodeToCellBox); bitOff.ok(); ++bitOff)
        {
          IntVect nodeToCell = bitOff();
          IntVect cell = thisVertex + nodeToCell;
          if (m_mappingBlocks[thisBlock].contains(cell))
            if ( !cellIndicesList.includes(cell) )
              {
                cellIndicesList.add(cell);
                cellBlocksList.add(thisBlock);
              }
        }
    }

  a_cellIndices.clear();
  a_cellBlocks.clear();
  while (cellIndicesList.isNotEmpty())
    {
      a_cellIndices.push_back(cellIndicesList.firstElement());
      a_cellBlocks.push_back(cellBlocksList.firstElement());
      cellIndicesList.removeFirst();
      cellBlocksList.removeFirst();
    }
  return;
}

void
MultiBlockUtil::validStencilElementsFromVertices(
  StencilMap&       a_stencilElements,
  const StencilMap& a_vertices) const
{
  Box nodeToCellBox(-IntVect_unit, IntVect_zero);
  for (const auto& v : a_vertices)
    {
      const int thisBlock = v.first.block();
      const IntVect& thisVertex = v.first.cell();
      const Box& thisBlockCellBox = m_mappingBlocks[thisBlock];
      MD_BOXLOOP(nodeToCellBox, nodeToCell)
        {
          IntVect cell = thisVertex + MD_GETIV(nodeToCell);
          if (thisBlockCellBox.contains(cell))
            {
              MBStencilElement cellElement(cell, thisBlock);
              // Only inserts if unique
              a_stencilElements.insert({ cellElement, v.second });
            }
        }
    }
}

bool
MultiBlockUtil::validCellShift(IntVect&        a_shiftedCell,
                               int&            a_shiftedBlock,
                               const IntVect&  a_origCell,
                               int             a_origBlock,
                               int             a_dir,
                               int             a_shift) const
{
  CH_assert(m_isDefined);
  a_shiftedCell = a_origCell;
  Box origBox = m_mappingBlocks[a_origBlock];
  CH_assert(origBox.contains(a_origCell));
  a_shiftedCell.shift(a_dir, a_shift);
  if (origBox.contains(a_shiftedCell))
    { // still in the same block
      a_shiftedBlock = a_origBlock;
      return true;
    }
  else
    {
      int faceID = a_dir;
      if (a_shift > 0) faceID += SpaceDim;

      const BlockBoundary& bb = m_boundaries[a_origBlock][faceID];
      if (bb.isDomainBoundary())
        {
          return false;
        }
      else
        {
          a_shiftedBlock = bb.neighbor();
          IndicesTransformation tfm = bb.getTransformation();
          a_shiftedCell = tfm.transform(a_shiftedCell);
          return true;
        }
    }
}

bool
MultiBlockUtil::validCellShift(MBStencilElement&       a_shiftedElt,
                               RigidTransformation&    a_extraDispl,
                               const MBStencilElement& a_origElt,
                               int                     a_dir,
                               int                     a_shift) const
{
  CH_assert(m_isDefined);
  const IntVect& origCell = a_origElt.cell();
  int origBlock = a_origElt.block();
  Box origBox = m_mappingBlocks[origBlock];
  CH_assert(origBox.contains(origCell));

  IntVect shiftedCell = origCell;
  shiftedCell.shift(a_dir, a_shift);
  int shiftedBlock = origBlock;
  bool retval = true;
  if (!origBox.contains(shiftedCell))
    {
      const Side::LoHiSide side = static_cast<Side::LoHiSide>(a_shift > 0);
      const BlockBoundary& bb = m_coordSysPtr->boundary(origBlock, a_dir, side);
      if (bb.isDomainBoundary())
        {
          retval = false;
        }
      else
        {
          shiftedBlock = bb.neighbor();
          IndicesTransformation tfm = bb.getTransformation();
          shiftedCell = tfm.transform(shiftedCell);
          a_extraDispl += bb.getPhysTransformation();
          retval = true;
        }
    }
  a_shiftedElt = MBStencilElement(shiftedCell, shiftedBlock);
  return retval;
}

//  Find a valid cell index for each ghost cell
/** The main technology used is blockRemapping
 */
void
MultiBlockUtil::getValid(IVSFAB<MBLGeomValidInfo>& a_validInfo,
                         const IntVectSet&         a_ghostCellsIVS,
                         const int                 a_baseBlockNum) const
{
  for (IVSIterator ivsit(a_ghostCellsIVS); ivsit.ok(); ++ivsit)
    {
      IntVect thisGhostCell = ivsit();

      MBLGeomValidInfo& thisValidInfo = a_validInfo(thisGhostCell, 0);
      getValid(thisValidInfo.m_idxCell,
               thisValidInfo.m_idxBlk,
               thisValidInfo.m_validExists,
               thisValidInfo.m_extraDispl,
               thisGhostCell,
               a_baseBlockNum);
    }
}

void
MultiBlockUtil::getValid(BaseFab<IntVect>& a_validIndices,
                         BaseFab<int>& a_validBlock,
                         BaseFab<bool>& a_validExists,
                         BaseFab<RigidTransformation>& a_extraDispl,
                         const Box& a_bx,
                         const int a_baseBlockNum) const
{
  CH_assert(m_isDefined);
  for (BoxIterator bit(a_bx); bit.ok(); ++bit)
    {
      IntVect thisGhostCell = bit();
      getValid(a_validIndices(thisGhostCell, 0),
               a_validBlock(thisGhostCell, 0),
               a_validExists(thisGhostCell, 0),
               a_extraDispl(thisGhostCell, 0),
               thisGhostCell,
               a_baseBlockNum);
    }
}

void
MultiBlockUtil::getValid(IntVect&             a_validIndices,
                         int&                 a_validBlock,
                         bool&                a_validExists,
                         RigidTransformation& a_extraDispl,
                         const IntVect&       a_ghostCell,
                         const int            a_baseBlock) const
{
  CH_assert(m_isDefined);
  const NewCoordSys* baseCoordSysPtr =
    m_coordSysPtr->getCoordSys(a_baseBlock);

  // baseCenter:  center of a_ghostCell in
  // mapped coordinates of block a_baseBlock
  RealVect baseCenter =
    baseCoordSysPtr->centerMappedCoordinates(a_ghostCell);

  // Point baseCenter in a_baseBlock's mapped coordinates becomes
  // validMappedCenter in a_validBlock's mapped coordinates, where
  // a_validBlock is the index of the block where the point is valid.
  RealVect validMappedCenter;
  // a_extraDipls is set to zero here in case particular CS do not consider it
  a_extraDispl.reset();
  m_coordSysPtr->blockRemapping(validMappedCenter,
                                a_validBlock,
                                a_validExists,
                                a_extraDispl,
                                baseCenter,
                                a_ghostCell,
                                a_baseBlock);
  CH_assert((a_validBlock >= 0) && (a_validBlock < m_nblocks));
  // a_validIndices is in index space of a_validBlock, transformed
  // from a_ghostCell in index space of baseBlock.
  const NewCoordSys* validCoordSysPtr =
    m_coordSysPtr->getCoordSys(a_validBlock);
  const RealVect& validDxi = validCoordSysPtr->dx();
  // WAS for (int idir = 0; idir < SpaceDim; idir++)
  for (int ind = 0; ind < m_interpDimsVect.size(); ind++)
    {
      int idir = m_interpDimsVect[ind];
      a_validIndices[idir] = floor(validMappedCenter[idir] / validDxi[idir]);
    }
  for (int ind = 0; ind < m_fixedDimsVect.size(); ind++)
    {
      int idir = m_fixedDimsVect[ind];
      a_validIndices[idir] = a_ghostCell[idir];
    }
  // You're in trouble if you go outside the domain.
  // So find the closest valid cell that is in the domain.
  const Box& validBox = m_mappingBlocks[a_validBlock];
  if (!validBox.contains(a_validIndices))
    {
      // This should happen only when directly on a block boundary.
      // WAS for (int idir = 0; idir < SpaceDim; idir++)
      for (int ind = 0; ind < m_interpDimsVect.size(); ind++)
        {
          int idir = m_interpDimsVect[ind];
          if (a_validIndices[idir] < validBox.smallEnd(idir))
            {
              a_validIndices[idir] = validBox.smallEnd(idir);
            }
          else if (a_validIndices[idir] > validBox.bigEnd(idir))
            {
              a_validIndices[idir] = validBox.bigEnd(idir);
            }
        }
    }
}

/** \param[out] a_stencils
 *                      Stencil to fill the ghost cell.  For each ghost cell,
 *                      this is a vector of cells.
 *  \param[out] a_stencilExtraDispl
 *                      Extra displacement used to get to a particular stencil
 *                      cell from the ghost cell.
 *  \param[out] a_minBox
 *                      The minimum box containing stencil elements for each
 *                      block in the domain.
 *  \param[in]  a_ghostCellsIVS
 *                      The ghost cells that must be filled during an exchange
 *  \param[in]  a_validIndices
 *                      The integer location of the valid cell that best
 *                      overlaps the ghost cell.  If the ghost cell is outside
 *                      the domain, this might be the nearest valid cell.
 *  \param[in]  a_validBlock
 *                      The block containing the valid cell
 *  \param[in]  a_order Order of interpolation (choose 1, 2, or 4)
 *  \param[in]  a_radius
 *                      Requested separation from domain boundaries for the
 *                      center cell of the inner stencil (1 for 4th-order)
 */
void
MultiBlockUtil::getStencilCells(
  IVSFAB<MBStencil>&                        a_stencils,
  IVSFAB<std::vector<RigidTransformation>>& a_stencilExtraDispl,
  Vector<Box>&                              a_minBox,
  const IntVectSet&                         a_ghostCellsIVS,
  const IVSFAB<MBLGeomValidInfo>&           a_validInfo,
  const int                                 a_order,
  const int                                 a_radius) const
{
  CH_assert(m_isDefined);
  // Start with empty a_minBox.
  a_minBox.assign(m_nblocks, Box{});
  for (IVSIterator ivsit(a_ghostCellsIVS); ivsit.ok(); ++ivsit)
    {
      IntVect thisGhostCell = ivsit();

      MBStencil stencil;

      // Information about the valid cell
      const MBLGeomValidInfo& thisValidInfo = a_validInfo(thisGhostCell, 0);

      getStencilCells(stencil,
                      a_stencilExtraDispl(thisGhostCell, 0),
                      a_minBox,
                      thisGhostCell,
                      thisValidInfo.m_idxCell,
                      thisValidInfo.m_idxBlk,
                      thisValidInfo.m_validExists,
                      thisValidInfo.m_extraDispl,
                      a_order,
                      a_radius);
      a_stencils(thisGhostCell, 0) = stencil;
    }
}

/** \param[out] a_stencil
 *                      Vector of elements in the stencil including integer
 *                      location and block.  The vector is resized herein.
 *  \param[out] a_extraDispl
 *                      Extra displacements for each element of the stencil
 *                      primarily due to periodic boundaries.  The vector is
 *                      resized herein.  The displacements reference the ghost
 *                      cell.
 *  \param[out] a_minBox
 *                      The minimum box containing stencil elements for each
 *                      block in the domain.
 *  \param[in]  a_ghostIndices
 *                      The integer location of the ghost cell.
 *  \param[in]  a_validIndices
 *                      The integer location of the valid cell that best
 *                      overlaps the ghost cell.  If the ghost cell is outside
 *                      the domain, this might be the nearest valid cell.
 *  \param[in]  a_validBlockNum
 *                      The block containing the valid cell
 *  \param[in]  a_order Order of interpolation (choose 1, 2, or 4)
 *  \param[in]  a_radius
 *                      Requested separation from domain boundaries for the
 *                      center cell of the inner stencil (1 for 4th-order)
 */
void
MultiBlockUtil::getStencilCells(
  MBStencil&                        a_stencil,
  std::vector<RigidTransformation>& a_stencilExtraDispl,
  Vector<Box>&                      a_minBox,
  const IntVect&                    a_ghostIndices,
  const IntVect&                    a_validIndices,
  const int                         a_validBlockNum,
  const bool                        a_validExists,
  const RigidTransformation&        a_validExtraDispl,
  const int                         a_order,
  const int                         a_radius) const
{
  MBStencilElement validStencilElement(a_validIndices, a_validBlockNum);

  StencilMap stencilMap;

  if (a_order <= 1 || !a_validExists)
    {
      // // if (a_ghostIndices == IntVect{279, 32})
      // //   {
      // //     validStencilElement.define(IntVect{255, 31}, 2);
      // //   }
      // // if (a_ghostIndices == IntVect{278, 31})
      // //   {
      // //     validStencilElement.define(IntVect{254, 31}, 2);
      // //   }
      // if (a_ghostIndices == IntVect{256, 31})
      //   {
      //     validStencilElement.define(IntVect{280, 31}, 3);
      //   }
      // if (a_ghostIndices == IntVect{257, 31})
      //   {
      //     validStencilElement.define(IntVect{281, 31}, 3);
      //   }
      stencilMap.insert({ validStencilElement, a_validExtraDispl });
      // // Box queryBoxProb(IntVect{280, 31}, IntVect{280, 31});
      // // queryBoxProb.grow(2);
      // // Box queryBoxProbA1(IntVect{278, 31}, IntVect{280, 31});
      // // Box queryBoxProbA2(IntVect{280, 31}, IntVect{280, 33});
      // Box queryBoxProbA1(IntVect{1471, 0}, IntVect{1473, 0});
      // Box queryBoxProbA2(IntVect{1471, -2}, IntVect{1471, 0});
      // Box queryBoxProbB1(IntVect{255, 31}, IntVect{257, 31});
      // Box queryBoxProbB2(IntVect{255, 31}, IntVect{255, 33});
      // // Box queryBoxCtrl(IntVect{904, 15}, IntVect{904, 15});
      // // queryBoxCtrl.grow(2);
      // // Box queryBoxCtrl1(IntVect{902, 15}, IntVect{904, 15});
      // // Box queryBoxCtrl2(IntVect{904, 15}, IntVect{904, 17});
      // Box queryBoxCtrl1(IntVect{687, 15}, IntVect{689, 15});
      // Box queryBoxCtrl2(IntVect{687, 15}, IntVect{687, 17});
      // if (queryBoxProbA1.contains(a_ghostIndices) ||
      //     queryBoxProbA2.contains(a_ghostIndices) ||
      //     queryBoxProbB1.contains(a_ghostIndices) ||
      //     queryBoxProbB2.contains(a_ghostIndices) ||
      //     queryBoxCtrl1.contains(a_ghostIndices) ||
      //     queryBoxCtrl2.contains(a_ghostIndices))
      //   {
      //     std::cout << "Stencil for " << a_ghostIndices
      //               << ": " << validStencilElement.cell()
      //               << " blk " << validStencilElement.block()
      //               << " exists " << a_validExists << std::endl;
      //   }
    }
  else
    {
      // Set centralCell to be the central cell of the neighborhood,
      // in a_validBlockNum's index space.  It is the cell closest to
      // thisValidCell that is separated from the external boundary
      // by at least a_radius cells.
      // (So centralCell == a_validIndices except near external boundary.)
      // We use only 4th or 5th order, and for those, a_radius == 1.
      IntVect centralCell =
        m_coordSysPtr->cellAvoidingBoundaries(a_validIndices,
                                              a_validBlockNum,
                                              a_radius);
      if (m_fixedDimsVect.size() > 0)
        {
          for (int ind = 0; ind < m_fixedDimsVect.size(); ind++)
            {
              int idir = m_fixedDimsVect[ind];
              centralCell[idir] = a_ghostIndices[idir];
            }
        }

      MBStencilElement baseStencilElement(centralCell, a_validBlockNum);
      // Stencil's inner set:
      // All valid cells having a vertex in common with centralCell.
      commonVertexStencilElements(stencilMap,
                                  baseStencilElement,
                                  a_validExtraDispl);
      if (m_fixedDimsVect.size() > 0)
        { // Remove all stencil cells that do not have 
          // the fixed dimensions set to what they should be.
          for (StencilMap::iterator cell = stencilMap.begin(),
                 cell_end = stencilMap.end(); cell != cell_end;)
            {
              const MBStencilElement& elt = cell->first;
              IntVect thisCellIndices = elt.cell();
              bool keepCell = true;
              for (int idxFxD = 0, idxFxD_end = m_fixedDimsVect.size();
                   idxFxD != idxFxD_end; ++idxFxD)
                {
                  int idxDir = m_fixedDimsVect[idxFxD];
                  if (thisCellIndices[idxDir] != a_ghostIndices[idxDir])
                    {
                      keepCell = false;
                    }
                }
              if (keepCell)
                {
                  ++cell;
                }
              else
                {
                  cell = stencilMap.erase(cell);
                }
            }
        }

      if (a_order > 2)
        {
          // To stencilVec, append
          // a_validIndices +/- 2*e, whenever within physical boundary.
          // Then we'll have either 4 or 5 in a row in every dimension.
          // Note this is from a_validIndices, not centralCell.
          IntVect offCenter = centralCell - a_validIndices;
          { CH_TIME("add to stencil");
            for (int ind = 0; ind < m_interpDimsVect.size(); ind++)
              {
                int idir = m_interpDimsVect[ind];
                int offCenterDir = offCenter[idir];
                for (const auto side : EachSide)
                  {
                    int sgn = sign(side);
                    // When offCenterDir == 0, can shift from a_validIndices
                    // by either -2 or +2.
                    
                    // When offCenterDir == -1, we have shifted left
                    // from a_validIndices to centralCell,
                    // so no more room on the right; do iff sgn == -1,
                    // and shift from a_validIndices by -3.
                    
                    // When offCenterDir == +1, we have shifted right
                    // from a_validIndices to centralCell,
                    // so no more room on the left; do iff sgn == +1,
                    // and shift from a_validIndices by +3.
                    if (sgn * offCenterDir >= 0)
                      {
                        int shiftAmount = sgn*2 + offCenterDir;
                        MBStencilElement shiftedElt;
                        RigidTransformation extraDispl(a_validExtraDispl);
                        if (validCellShift(shiftedElt,
                                           extraDispl,
                                           validStencilElement,
                                           idir,
                                           shiftAmount))
                          {
                            stencilMap.insert({ shiftedElt, extraDispl });
                          }
                      }
                  }
              }
          }
        }
    }

  // So much for a_order == 4.
  // If a_order == 5, need more stencil elements in stencilVec.

//--Now update minbox and separate the output into calling API

  a_stencil = RefCountedPtr< Vector<MBStencilElement> >
    (new Vector<MBStencilElement>(stencilMap.size()));
  a_stencilExtraDispl.clear();
  a_stencilExtraDispl.reserve(stencilMap.size());
  {
    int i = 0;
    for (const auto& elem : stencilMap)
      {
        // Fill output
        a_stencil[i++] = elem.first;
        a_stencilExtraDispl.emplace_back(elem.second);

        // Set minbox
        IntVect thisStencilCell = elem.first.cell();
        int thisStencilBlock = elem.first.block();
        Box& minBoxSrc = a_minBox[thisStencilBlock];
        if ( !minBoxSrc.contains(thisStencilCell) )
          {
            Box thisStencilCellBox(thisStencilCell, thisStencilCell);
            minBoxSrc.minBox(thisStencilCellBox);
          }
      }
  }
}


void
MultiBlockUtil::getVectorStencilCells(IVSFAB<MBVectorStencil>& a_vectorstencils,
                                      const IVSFAB<MBStencil>& a_stencils,
                                      const IntVectSet& a_ghostCellsIVS) const
{
  for (IVSIterator ivsit(a_ghostCellsIVS); ivsit.ok(); ++ivsit)
    {
      IntVect thisGhostCell = ivsit();

      const MBStencil& stencil = a_stencils(thisGhostCell, 0);
      int npts = stencil.size();

      Vector<MBVectorStencilElement> vectorstencilVec(npts);
      for (int ipt = 0; ipt < npts; ipt++)
        {
          const MBStencilElement& stencilElement = stencil[ipt];
          const IntVect& cell = stencilElement.cell();
          int block = stencilElement.block();
          vectorstencilVec[ipt] = MBVectorStencilElement(cell, block);
        }

      MBVectorStencil vectorstencil =
        RefCountedPtr< Vector<MBVectorStencilElement> >
        (new Vector<MBVectorStencilElement>(vectorstencilVec));

      a_vectorstencils(thisGhostCell, 0) = vectorstencil;
    }
}

void
MultiBlockUtil::displacementPowers(
  Vector<Real>&              a_powers,
  const IntVect&             a_cell,
  const int                  a_cellBlockNum,
  const RigidTransformation& a_cellExtraDispl,
  const RealVect&            a_commonXBasePoint,
  const Real                 a_avgDistance,
  const Vector<IntVect>&     a_exponents,
  XCoordCache&               a_xCache) const
{ CH_TIME("displacementPowers");
  // This function will return
  // a_powers[pvec] =
  // <((Xvec() - Xvec(a_commonMappedBasePoint))/a_avgDistance)^pvec>_(a_cell),
  // where a_cell is in index space of block a_thisBlockNum.
  const NewCoordSys* coordSysCellBlockPtr =
    m_coordSysPtr->getCoordSys(a_cellBlockNum);

  // Initialize powers
  a_powers.assign(0.);

  // The loop here covers all cells in the convolution, including the center
  for (int idxConv = 0; idxConv != m_convSize; ++idxConv)
    { CH_TIME("powers of displacement");
      const IntVect convCell = a_cell + m_convBaseCells[idxConv];
      const RealVect convXi =
        coordSysCellBlockPtr->centerMappedCoordinates(convCell);
      // This is the physical location of the cell as expressed in the
      // neighborhood of the ghost cell.
      const RealVect convX = a_cellExtraDispl.transformBack(
        m_coordSysPtr->cachedXCoord(
          convCell, a_cellBlockNum, convXi, a_xCache));
      // The displacement is from the center of the common base to the center
      // of the cell in the convolution.
      const RealVect convDisplacementNormalized =
        (convX - a_commonXBasePoint)/a_avgDistance;

      // Increment a_powers += m_convWeight[iconv] *
      // ((Xvec() - Xvec(ghostcenter))/avgDistance)^pvec_(convCell),
      // where Xvec() has real coordinates,
      // and convCell is in index space of thisStencilCellBlock.
      addPowersPoint(a_powers, convDisplacementNormalized,
                     a_exponents, m_convWeight[idxConv]);      
    }
}
      
void
MultiBlockUtil::displacementPowersTransformed(Vector<Real>&  a_powers,
                                              const BaseFab<VectorTransformation>& a_vectorStencilTransformations,
                                              const IntVect&     a_cell,
                                              const RealVect&    a_commonMappedBasePoint,
                                              Real               a_avgDistance,
                                              int                a_thisBlockNum,
                                              int                a_commonBlockNum,
                                              int                a_ghostBlockNum,
                                              const Vector<IntVect>&   a_exponents) const
{ CH_TIME("displacementPowersTransformed");
  // This function will return
  // a_powers[pvec] =
  // <T(Xvec) * ((Xvec() - Xvec(a_commonMappedBasePoint))/a_avgDistance)^pvec>_(a_cell),
  // where a_cell is in index space of block a_thisBlockNum.
  // We do this for all SpaceDim*SpaceDim components of T(Xvec).
  const NewCoordSys* coordSysThisBlockPtr =
    m_coordSysPtr->getCoordSys(a_thisBlockNum);

  Vector<RealVect> convCellMappedCenters(m_convSize); // 1+2*nInterpDims
  for (int iconv = 0; iconv < m_convSize; iconv++)
    {
      IntVect convCell = a_cell + m_convBaseCells[iconv];
      convCellMappedCenters[iconv] =
        coordSysThisBlockPtr->centerMappedCoordinates(convCell);
    }

  Vector<int> convSrcBlocks(m_convSize, a_thisBlockNum);

  // For each cell i in the convolution stencil,
  // RealVect convDisplacements[i] is displacement
  // from point a_commonMappedBasePoint in block a_commonBlockNum
  // to point convCellMappedCenters[i] in block convSrcBlocks[i].
  Vector<RealVect> convDisplacements =
    m_coordSysPtr->displacements(convCellMappedCenters,
                                 convSrcBlocks,
                                 a_commonMappedBasePoint,
                                 a_commonBlockNum);

  // int numTaylorCoeffs = a_powers.size() / (SpaceDim*SpaceDim);
  a_powers.assign(0.); // length SpaceDim*SpaceDim * numTaylorCoeffs

  for (int iconv = 0; iconv < m_convSize; iconv++)
    { CH_TIME("powers of displacement");
      RealVect convDisplacementsNormalized =
        convDisplacements[iconv] / a_avgDistance;

      // Set basePowers = m_convWeight[iconv] *
      // ((Xvec() - Xvec(ghostcenter))/avgDistance)^pvec_(convCell),
      // where Xvec() has real coordinates,
      // and convCell is in index space of thisStencilCellBlock.
      Real weight = m_convWeight[iconv];
      // Multiply weight by T(Xvec), where we get T(XVec) from
      // m_coordSysPtr function that uses
      // convCellMappedCenters, convSrcBlocks, a_commonBlockNum.

      // For each cell i in the convolution stencil,
      // RealVect convDisplacements[i] is displacement
      // from point a_commonMappedBasePoint in block a_commonBlockNum
      // to point convCellMappedCenters[i] in block convSrcBlocks[i].
      
      // VectorTransformation vectorBlockTransformation(int a_nDst,
      // const RealVect& a_xiSrc, int a_nSrc) const;
      // VectorTransformation vtToGhost transforms components of a vector
      // in basis of block thisVectorStencilCellBlock
      // to components in a basis in block a_ghostBlockNum.
      // The vector is situated at point thisVectorStencilCellCenter
      // of block thisVectorStencilCellBlock.
      IntVect convCell = a_cell + m_convBaseCells[iconv];
      VectorTransformation vt = a_vectorStencilTransformations(convCell, 0);

      // Now set a_powers += basePowers * T(Xvec).
      addPowersVector(a_powers, convDisplacementsNormalized,
                      a_exponents, weight, vt.dataPtr());
    }
}
      
void
MultiBlockUtil::getVectorTransformations(Vector< BaseFab<VectorTransformation>* >& a_vectorTransformations,
                                         const IVSFAB<MBVectorStencil>& a_vectorstencilsFab,
                                         const IndexType& a_type,
                                         const IntVectSet& a_ghostCellsIVS,
                                         int a_ghostBlockNum) const
{ CH_TIME("MultiBlockUtil::getVectorTransformations");
  Vector< BaseFab<bool>* > gotVectorTransformations(m_nblocks);
  for (int srcBlock = 0; srcBlock < m_nblocks; srcBlock++)
    {
      const Box& bxSrc = a_vectorTransformations[srcBlock]->box();
      gotVectorTransformations[srcBlock] = new BaseFab<bool>(bxSrc, 1);
      gotVectorTransformations[srcBlock]->setVal(false);
    }
  for (IVSIterator ivsit(a_ghostCellsIVS); ivsit.ok(); ++ivsit)
    {
      IntVect thisGhostCell = ivsit();
      const MBVectorStencil& thisGhostStencil = a_vectorstencilsFab(thisGhostCell, 0);
      int vectorstencilSize = thisGhostStencil.size();
      for (int icell = 0; icell < vectorstencilSize; icell++)
        {
          const MBVectorStencilElement& vectorstencilElement =
            thisGhostStencil[icell];
          const IntVect& thisVectorStencilCellIndices =
            vectorstencilElement.cell();
          int thisVectorStencilCellBlock =
            vectorstencilElement.block();
          BaseFab<VectorTransformation>& vectorTransformationsSrc =
            *(a_vectorTransformations[thisVectorStencilCellBlock]);
          BaseFab<bool>& gotVectorTransformationsSrc =
            *(gotVectorTransformations[thisVectorStencilCellBlock]);
          // lambda for setting up transformations
          auto setupVectorTransformation = [&] (const IntVect& a_cellIndex)
          {
            if (!gotVectorTransformationsSrc(a_cellIndex, 0))
              {
                VectorTransformation vtToGhost =
                  m_coordSysPtr->vectorBlockTransformationCenter(
                    a_ghostBlockNum,
                    a_cellIndex,
                    thisVectorStencilCellBlock);
                vectorTransformationsSrc(a_cellIndex, 0) = vtToGhost.inverse();
                gotVectorTransformationsSrc(a_cellIndex, 0) = true;
              }
          };
          // set transform depending on index type
          if (a_type == IndexType::TheNodeType())
            {
              setupVectorTransformation(thisVectorStencilCellIndices);
            }
          else if (a_type == IndexType::TheCellType())
            {
              for (int iconv = 0; iconv < m_convSize; iconv++)
                {
                  IntVect convCellIndices = 
                    thisVectorStencilCellIndices + m_convBaseCells[iconv];
                  setupVectorTransformation(convCellIndices);
                }
            }
        }
    }
  for (int srcBlock = 0; srcBlock < m_nblocks; srcBlock++)
    {
      delete gotVectorTransformations[srcBlock];
    }
}

void
MultiBlockUtil::getWeights(
  IVSFAB<MBStencil>&                              a_stencilsFab,
  const IVSFAB<std::vector<RigidTransformation>>& a_stencilExtraDispl,
  const IndexType&                                a_type,
  const IntVectSet&                               a_ghostCellsIVS,
  const int                                       a_ghostBlockNum,
  const IVSFAB<MBLGeomValidInfo>&                 a_validInfo,
  const Vector<IntVect>&                          a_exponents,
  XCoordCache&                                    a_xCache) const
{ CH_TIME("MultiBlockUtil::getWeights");
  CH_assert(m_isDefined);
  for (IVSIterator ivsit(a_ghostCellsIVS); ivsit.ok(); ++ivsit)
    {
      IntVect thisGhostCell = ivsit();
      MBStencil& thisGhostStencil = a_stencilsFab(thisGhostCell, 0);

      // Information about the valid cell
      const MBLGeomValidInfo& thisValidInfo = a_validInfo(thisGhostCell, 0);

      getWeights(thisGhostStencil,         // Stencil
                 a_stencilExtraDispl(thisGhostCell, 0),
                 a_type,
                 thisGhostCell,            // Ghost cell
                 a_ghostBlockNum,
                 thisValidInfo.m_idxCell,  // Valid cell
                 thisValidInfo.m_idxBlk,
                 thisValidInfo.m_validExists,
                 thisValidInfo.m_extraDispl,
                 a_exponents,
                 a_xCache);
    }
}

/** \param[in]  a_stencil
 *                      Vector of elements in the stencil including integer
 *                      location and block.
 *  \param[out] a_stencil
 *                      Weights added
 *  \param[in]  a_stencilExtraDispl
 *                      Extra displacements for each element of the stencil
 *                      w.r.t. the ghost cell, primarily due to periodic
 *                      boundaries.
 *  \param[in]  a_ghostIndices
 *                      The integer location of the ghost cell
 *  \param[in]  a_ghostBlockNum
 *                      Block number of the ghost cell
 *  \param[in]  a_validIndices
 *                      The integer location of the valid cell that best
 *                      overlaps the ghost cell.  If the ghost cell is outside
 *                      the domain, this might be the nearest valid cell.
 *  \param[in]  a_validBlockNum
 *                      The block containing the valid cell
 *  \param[in]  a_validExists
 *                      T - There is a valid cell overlapping the ghost cell
 *                      F - The ghost cell is outside the domain
 *  \param[in]  a_validExtraDispl
 *                      Extra displacements for the valid cell w.r.t. the ghost
 *                      cell, primarily due to periodic boundaries.
 *  \param[in]  a_exponents
 *                      Exponents from the desired interpolation
 *  \param[in]  a_xCache
 *                      Includes methods for evaluating xi->x and caching the
 *                      results for multiple lookup
 */
void
MultiBlockUtil::getWeights(
  MBStencil&                              a_stencil,
  const std::vector<RigidTransformation>& a_stencilExtraDispl,
  const IndexType&                        a_type,
  const IntVect&                          a_ghostIndices,
  const int                               a_ghostBlockNum,
  const IntVect&                          a_validIndices,
  const int                               a_validBlockNum,
  const bool                              a_validExists, 
  const RigidTransformation&              a_validExtraDispl,
  const Vector<IntVect>&                  a_exponents,
  XCoordCache&                            a_xCache) const
{ CH_TIME("MultiBlockUtil::getWeights on one");
  const int stencilSize = a_stencil.size();

  // There is no stencil for exterior ghost cells, and instead it is just
  // a direct copy
  if (!a_validExists)
    {
      CH_assert(stencilSize == 1);
      a_stencil[0].setWeight(1.0);
      return; // skip out of the function
    }

  // Get coordinate of valid cell in computational space
  const RealVect validXiCenter = m_coordSysPtr->getCoordSys(a_validBlockNum)
    ->centerMappedCoordinates(a_validIndices);
  // This is the physical location of the valid cell as expressed in the
  // neighborhood of the ghost cell.
  const RealVect validXCenter = a_validExtraDispl.transformBack(
    m_coordSysPtr->cachedXCoord(a_validIndices,
                                a_validBlockNum,
                                validXiCenter,
                                a_xCache));
  
  // Set to total of the distances, in real space,
  // from center of ghost cell to center of each stencil cell.
  Real totalDistance = 0.;
  // Mapped-space coordinates of center of each valid cell in the stencil.
  Vector<RealVect> stencilXiCenter(stencilSize);
  // Block of each valid cell in the stencil.
  Vector<int> stencilBlocks(stencilSize);
  { CH_TIME("total distance");
    for (int icell = 0; icell < stencilSize; icell++)
      {
        const MBStencilElement& stencilElement = a_stencil[icell];
        const IntVect& thisStencilCellIndices = stencilElement.cell();
        int thisStencilCellBlock = stencilElement.block();
        stencilBlocks  [icell] = thisStencilCellBlock;
        stencilXiCenter[icell] =
          m_coordSysPtr->getCoordSys(thisStencilCellBlock)
          ->centerMappedCoordinates(thisStencilCellIndices);
        // Center of ghost cell in mapped coordinates is
        // a_validMappedCenter in block a_validBlockNum.
//Well, this is what we have to fix
//Should only query valid locations and add in (or substract) extra displacements
        Real distance = m_coordSysPtr->distance(
          validXCenter,
          thisStencilCellIndices,
          thisStencilCellBlock,
          stencilXiCenter[icell],
          a_stencilExtraDispl[icell],
          a_xCache);
        totalDistance += distance;
      }
  }
  Real avgDistance = totalDistance / Real(stencilSize);

  // a_exponents contains the exponents, SpaceDim at a time.
  // With 2D and degree 3 and no fixed dimensions, these are:
  // 00 10 20 30 01 11 21 02 12 03.
  // These are Vector<int> instead of Vector<IntVect> because
  // we need to send them to Chombo Fortran.
  int numTaylorCoeffs = a_exponents.size();

  // To hold displacement powers (centered or averaged) of stencil cells,
  // numTaylorCoeffs at a time.
  Vector<Real> coeffs;
  if (a_type == IndexType::TheNodeType())
    { CH_TIME("centered displacement powers"); // cell-centered only
      for (int icell = 0; icell < stencilSize; icell++)
        {
          const MBStencilElement& stencilElement = a_stencil[icell];
          // This is the physical location of the stencil cell as expressed in
          // the neighborhood of the ghost cell.
          const RealVect stencilXCenter =
            a_stencilExtraDispl[icell].transformBack(
              m_coordSysPtr->cachedXCoord(stencilElement.cell(),
                                          stencilElement.block(),
                                          stencilXiCenter[icell],
                                          a_xCache));
          const RealVect dispSrcCell =
            (stencilXCenter - validXCenter)/avgDistance;
          // Set coeffsCell[pvec] =
          // ((Xvec() - Xvec(ghostcenter))/avgDistance)^pvec_(thisStencilCellIndices),
          // where Xvec() has real coordinates,
          // a_exponents holds the powers pvec,
          // and thisStencilCellIndices in index space of thisStencilCellBlock.
          Vector<Real> coeffsCell(numTaylorCoeffs, 0.);
          addPowersPoint(coeffsCell, dispSrcCell, a_exponents, 1.);
          coeffs.append(coeffsCell);
        } // end loop over cells in stencil
    }
  else if (a_type == IndexType::TheCellType())
    { CH_TIME("averaged displacement powers"); // cell-averaged only
      Vector<Real> dispAvgPowers(numTaylorCoeffs, 0.);
      for (int icell = 0; icell < stencilSize; icell++)
        {
          const MBStencilElement& stencilElement = a_stencil[icell];
          displacementPowers(dispAvgPowers,
                             stencilElement.cell(),
                             stencilElement.block(),
                             a_stencilExtraDispl[icell],
                             validXCenter,
                             avgDistance,
                             a_exponents,
                             a_xCache);
          coeffs.append(dispAvgPowers);
        } // end loop over cells in stencil
    }
  else
    {
      MayDay::Error("Bad index type");
    }

  // Now using these powers, find weights.

  // stencilSize equations or function evaluations.
  // numTaylorCoeffs variables as displacement powers or coefficients.

  // Vector<Real> coeffs has length stencilSize * numTaylorCoeffs.
  // - On entry, coeffs holds powers of displacement,
  // numTaylorCoeffs at a time for each stencil cell.
  // - On exit, coeffs is replaced by coefficients,
  // numTaylorCoeffs at a time for each stencil cell.
  LAPACKMatrix A(stencilSize, numTaylorCoeffs);
  {
    int icoeff = 0;
    for (int icell = 0; icell < stencilSize; icell++)
      {
        for (int ico = 0; ico < numTaylorCoeffs; ico++)
          {
            A(icell, ico) = coeffs[icoeff];
            icoeff++;
          }
      }
  }
  A.pseudoInvertUsingQR();
  {
    int icoeff = 0;
    for (int icell = 0; icell < stencilSize; icell++)
      {
        for (int ico = 0; ico < numTaylorCoeffs; ico++)
          {
            coeffs[icoeff] = A(ico, icell);
            icoeff++;
          }
      }
  }

  // From coeffs(numTaylorCoeffs*icell+(0:numTaylorCoeffs-1))
  // need to find Real wt,
  // and then do a_stencil[icell].setWeight(wt);
  // This is the weight in the interpolation stencil.
  if (a_type == IndexType::TheNodeType())
    { // cell-centered only:
      // At point Xvec in the ghost cell,
      // weight of icell in thisStencilCellIndices is
      // coeffs[icell, pvec] *
      // ((Xvec() - Xvec(ghostcenter))/avgDistance)^pvec_icell.
      // At cell center, Xvec() = Xvec(ghostcenter),
      // so we can ignore the terms with nonzero pvec, and
      // this is simply coeffs[icell, 0].
      int stencilBaseIndex = 0;
      for (int icell = 0; icell < stencilSize; icell++)
        {
          const Real wt = coeffs[stencilBaseIndex];
          a_stencil[icell].setWeight(wt);
          stencilBaseIndex += numTaylorCoeffs;
        }
    }
  else if (a_type == IndexType::TheCellType())
    { // cell-averaged only:
      /*
        thisDispPowers[pvec] =
        <((Xvec() - Xvec(ghostcenter))/avgDistance)^pvec>_g
        where Xvec() has real coordinates.

        If we're interpolating cell centers instead of cell averages,
        then we don't need thisDispPowers, because
        the only Xvec() will be Xvec(ghostcenter), and
        ((Xvec(ghostcenter) - Xvec(ghostcenter))/avgDistance)^pvec_g
        == delta(pvec, 0).
        In other words, when interpolating cell centers, we use
        only the first coefficient of the Taylor series, and
        ignore the other coefficients.
      */
      // CORRECTION: This is the location of the ghost cell relative to the
      //             valid cell!
      Vector<Real> thisDispPowers(numTaylorCoeffs);
      displacementPowers(thisDispPowers,
                         a_ghostIndices,
                         a_ghostBlockNum,
                         // The displacements reference the ghost cell itself
                         RigidTransformation{},
                         validXCenter,
                         avgDistance,
                         a_exponents,
                         a_xCache);

      // At point Xvec in the ghost cell,
      // weight of icell in thisStencilCellIndices is
      // coeffs[icell, pvec] *
      // <((Xvec() - Xvec(ghostcenter))/avgDistance)^pvec>_icell
      // where the average is over the ghost cell.
      int stencilBaseIndex = 0;
      for (int icell = 0; icell < stencilSize; icell++)
        {
          const Real wt = dotSubvectors(coeffs,
                                        stencilBaseIndex,
                                        thisDispPowers,
                                        0,
                                        numTaylorCoeffs);
          a_stencil[icell].setWeight(wt);
          stencilBaseIndex += numTaylorCoeffs;
        }
    }
}

void
MultiBlockUtil::getVectorWeights(IVSFAB<MBVectorStencil>& a_vectorstencilsFab,
                                 const Vector< BaseFab<VectorTransformation>* >& a_vectorStencilTransformations,
                                 const IndexType& a_type,
                                 const IntVectSet& a_ghostCellsIVS,
                                 const IVSFAB<IntVect>& a_validIndices,
                                 const IVSFAB<int>& a_validBlock,
                                 const IVSFAB<RealVect>& a_validMappedCenter,
                                 const Vector<IntVect>& a_exponents,
                                 int a_ghostBlockNum) const
{ CH_TIME("MultiBlockUtil::getVectorWeights");
  CH_assert(m_isDefined);
  for (IVSIterator ivsit(a_ghostCellsIVS); ivsit.ok(); ++ivsit)
    {
      IntVect thisGhostCell = ivsit();
      MBVectorStencil& thisGhostVectorStencil =
        a_vectorstencilsFab(thisGhostCell, 0);
      // cell in block validBlockNum where validMappedCenterHere is
      const IntVect& validIndices = a_validIndices(thisGhostCell, 0);
      // which valid block contains thisGhostCell
      int validBlockNum = a_validBlock(thisGhostCell, 0);
      // center of thisGhostCell in valid block's mapped coordinates
      const RealVect& validMappedCenterHere = a_validMappedCenter(thisGhostCell, 0);
      getVectorWeights(thisGhostVectorStencil, a_vectorStencilTransformations,
                       a_type, thisGhostCell,
                       validIndices, validBlockNum, validMappedCenterHere,
                       a_exponents, a_ghostBlockNum);
    }
}


void
MultiBlockUtil::getVectorWeights(MBVectorStencil& a_vectorstencil,
                                 const Vector< BaseFab<VectorTransformation>* >& a_vectorStencilTransformations,
                                 const IndexType& a_type,
                                 const IntVect& a_ghostIndices,
                                 const IntVect& a_validIndices,
                                 int a_validBlockNum,
                                 const RealVect& a_validMappedCenter,
                                 const Vector<IntVect>& a_exponents,
                                 int a_ghostBlockNum) const
{ CH_TIME("MultiBlockUtil::getVectorWeights on one");
  /*
    FIXME: this is a copy of a section of getWeights().
   */
  MayDay::Error("Is this used?");
  int vectorstencilSize = a_vectorstencil.size();

  // Set to total of the distances, in real space,
  // from center of ghost cell to center of each stencil cell.
  Real totalDistance = 0.;
  // Indices of each valid cell in the stencil.
  Vector<IntVect> vectorstencilCellIndices(vectorstencilSize);
  // Mapped-space coordinates of center of each valid cell in the stencil.
  Vector<RealVect> vectorstencilCellCenters(vectorstencilSize);
  // Block of each valid cell in the stencil.
  Vector<int> vectorstencilBlocks(vectorstencilSize);
  { CH_TIME("total distance");
    for (int icell = 0; icell < vectorstencilSize; icell++)
      {
        const MBVectorStencilElement& vectorstencilElement = a_vectorstencil[icell];
        const IntVect& thisVectorStencilCellIndices = vectorstencilElement.cell();
        vectorstencilCellIndices[icell] = thisVectorStencilCellIndices;
        int thisVectorStencilCellBlock = vectorstencilElement.block();
        vectorstencilBlocks[icell] = thisVectorStencilCellBlock;
        vectorstencilCellCenters[icell] =
          m_coordSysPtr->getCoordSys(thisVectorStencilCellBlock)
          ->centerMappedCoordinates(thisVectorStencilCellIndices);
        totalDistance +=
          m_coordSysPtr->distance(vectorstencilCellCenters[icell],
                                  thisVectorStencilCellBlock,
                                  a_validMappedCenter,
                                  a_validBlockNum);
      }
  }
  Real avgDistance = totalDistance / Real(vectorstencilSize);

  // a_exponents contains the exponents, SpaceDim at a time.
  // With 2D and degree 3 and no fixed dimensions, these are:
  // 00 10 20 30 01 11 21 02 12 03.
  // These are Vector<int> instead of Vector<IntVect> because
  // we need to send them to Chombo Fortran.
  int numTaylorCoeffs = a_exponents.size();
  int numVectorAllCoeffs = SpaceDim*SpaceDim * numTaylorCoeffs;

  /*
    This is used only if cell-centered.
    FIXME: this is a copy of a section of getWeights().
  */
  Vector<RealVect> disp;
  if (a_type == IndexType::TheNodeType())
    { // cell-centered only
      // Return Vector<RealVect> disp, where for each stencil cell i,
      // RealVect disp[i] is displacement
      // from the point a_validMappedCenter in block a_validBlockNum 
      // to the point vectorstencilCellCenters[i] in block vectorstencilBlocks[i].
      // These are displacements in real space
      // from the center of the ghost cell to the center of each stencil cell.
      // Center of ghost cell in mapped coordinates is
      // a_validMappedCenter in block a_validBlockNum.
      disp = m_coordSysPtr->displacements(vectorstencilCellCenters,
                                          vectorstencilBlocks,
                                          a_validMappedCenter,
                                          a_validBlockNum);
      // For each i, we will get monomials from disp[i],
      // and multiply the monomials by all the components of the
      // basis transform matrix.
    }

  /*
    If cell-averaged:

    thisDispPowersTransformed[dst, src, pvec] =
    T(Xvec())[dst,src] * <((Xvec() - Xvec(ghostcenter))/avgDistance)^pvec>_g
    where Xvec() has real coordinates.

    If we're interpolating cell centers instead of cell averages,
    then we don't need thisDispPowersTransforemd, because
    the only Xvec() will be Xvec(ghostcenter), and
    ((Xvec(ghostcenter) - Xvec(ghostcenter))/avgDistance)^pvec_g
    == delta(pvec, 0).
    In other words, when interpolating cell centers, we use
    only the first coefficient of the Taylor series, and
    ignore the other coefficients.
  */
  Vector<Real> thisDispPowers;
  if (a_type == IndexType::TheCellType())
    { CH_TIME("center displacement powers"); // cell-averaged only
      thisDispPowers.resize(numTaylorCoeffs);
      // displacementPowers(thisDispPowers,
      //                    a_ghostIndices,
      //                    a_validMappedCenter,
      //                    avgDistance,
      //                    a_ghostBlockNum,
      //                    a_validBlockNum,
      //                    a_exponents);
    }

  // To hold transformed displacement powers (centered or averaged) of stencil cells,
  // numVectorAllCoeffs = SpaceDim*SpaceDim * numTaylorCoeffs at a time.
  // For each stencil cell, have SpaceDim*SpaceDim groups, each group
  // having numTaylorCoeffs at a time corresponding to a matrix element.
  Vector<Real> coeffs;
  if (a_type == IndexType::TheNodeType())
    { CH_TIME("centered displacement powers"); // cell-centered only
      for (int icell = 0; icell < vectorstencilSize; icell++)
        {
          // Recall: disp[icell] is displacement in real space
          // from the center of the valid cell containing the ghost cell
          // to the center of stencil cell icell.
          RealVect dispSrcCell = disp[icell] / avgDistance;
          // Set baseCoeffsCell[pvec] =
          // ((Xvec() - Xvec(ghostcenter))/avgDistance)^pvec_(thisVectorStencilCellIndices),
          // where Xvec() has real coordinates,
          // a_exponents holds the powers pvec,
          // and thisVectorStencilCellIndices in index space of thisVectorStencilCellBlock.
          Vector<Real> baseCoeffsCell(numTaylorCoeffs, 0.);
          addPowersPoint(baseCoeffsCell, dispSrcCell,
                         a_exponents, 1.);

          // Now set coeffs = baseCoeffsCell * T(Xvec).
          // Need vectorstencilCellCenters[icell] and vectorstencilBlocks[icell].

          // Store in this order for each stencil cell (say 2D, degree 3):
          // [1, x, x^2, x^3, y, x*y, x*2*y, y^2, x*y^2, y^3]*T00,
          // [1, x, x^2, x^3, y, x*y, x*2*y, y^2, x*y^2, y^3]*T01,
          // [1, x, x^2, x^3, y, x*y, x*2*y, y^2, x*y^2, y^3]*T10,
          // [1, x, x^2, x^3, y, x*y, x*2*y, y^2, x*y^2, y^3]*T11.
          // The first two rows above make up one row in the matrix:
          // coefficients of 1st function component evaluated at stencil cell.
          // The second two rows above make up one row in the matrix:
          // coefficients of 2nd function component evaluated at stencil cell.
          // Here (x,y) is the normalized displacement in real space
          // from the ghost cell center to the stencil cell center.
          // ( T00 T01 )
          // ( T10 T11 ) is the transformation matrix at stencil cell center
          // from ghost cell block to stencil cell block.
          // The overdetermined system is (P < N):
          // ( x1^p*y1^q*T00(x1,y1) x1^p*y1^q*T01(x1,y1) )( a1 )   ( f(x1,y1) )
          // ( x1^p*y1^q*T10(x1,y1) x1^p*y1^q*T11(x1,y1) )( :  )   ( g(x1,y1) )
          // (          :                    :           )( aP ) = (  :       )
          // (          :                    :           )( b1 ) = (  :       )
          // ( xN^p*yN^q*T00(xN,yN) xN^p*yN^q*T01(xN,yN) )( :  )   ( f(xN,yN) )
          // ( xN^p*yN^q*T10(xN,yN) xN^p*yN^q*T11(xN,yN) )( bP )   ( g(xN,yN) )
          // because for each stencil cell i,
          // f(xi,yi) = sum_{p,q}
          //   a_{p,q}*xi^p*yi^q*T00(xi,yi) + b_{p,q}*xi^p*yi^q*T01(xi,yi) ;
          // g(xi,yi) = sum_{p,q}
          //   a_{p,q}*xi^p*yi^q*T10(xi,yi) + b_{p,q}*xi^p*yi^q*T11(xi,yi) .
          RealVect thisVectorStencilCellCenter = vectorstencilCellCenters[icell];

          int thisVectorStencilCellBlock = vectorstencilBlocks[icell];
          // VectorTransformation vectorBlockTransformation(int a_nDst,
          // const RealVect& a_xiSrc, int a_nSrc) const;
          // VectorTransformation vtToGhost transforms components of a vector
          // in basis of block thisVectorStencilCellBlock
          // to components in a basis in block a_ghostBlockNum.
          // The vector is situated at point thisVectorStencilCellCenter
          // of block thisVectorStencilCellBlock.
          IntVect thisVectorStencilCellIndices = vectorstencilCellIndices[icell];
          VectorTransformation vt =
            (*a_vectorStencilTransformations[thisVectorStencilCellBlock])(thisVectorStencilCellIndices, 0);
          for (int dst = 0; dst < SpaceDim; dst++)
            {
              for (int src = 0; src < SpaceDim; src++)
                {
                  Real dstTransComp = vt.component(dst, src);
                  Vector<Real> transformedCoeffsCell(numTaylorCoeffs);
                  for (int ico = 0; ico < numTaylorCoeffs; ico++)
                    {
                      transformedCoeffsCell[ico] = dstTransComp *
                        baseCoeffsCell[ico];
                    }
                  coeffs.append(transformedCoeffsCell);
                }
            }
        } // end loop over cells in stencil
    }
  else if (a_type == IndexType::TheCellType())
    { CH_TIME("averaged displacement powers"); // cell-averaged only
      for (int icell = 0; icell < vectorstencilSize; icell++)
        {
          const MBVectorStencilElement& vectorstencilElement = a_vectorstencil[icell];
          const IntVect& thisVectorStencilCellIndices = vectorstencilElement.cell();
          int thisVectorStencilCellBlock = vectorstencilElement.block();
          
          Vector<Real> dispAvgPowers(numVectorAllCoeffs, 0.);
          // There are SpaceDim*SpaceDim components of T(x,y).
          displacementPowersTransformed(dispAvgPowers,
                                        *(a_vectorStencilTransformations[thisVectorStencilCellBlock]),
                                        thisVectorStencilCellIndices,
                                        a_validMappedCenter,
                                        avgDistance,
                                        thisVectorStencilCellBlock,
                                        a_validBlockNum,
                                        a_ghostBlockNum,
                                        a_exponents);
          coeffs.append(dispAvgPowers);
        } // end loop over cells in stencil
    }
  else
    {
      MayDay::Error("Bad index type");
    }

  // Now using these powers, find weights.

  // numFunctionValues = SpaceDim * vectorstencilSize equations or
  // function evaluations.
  // numVectorCoeffs = SpaceDim * numTaylorCoeffs variables as
  // displacement powers or coefficients.

  // Vector<Real> coeffs has length vectorstencilSize * numVectorAllCoeffs.
  // - On entry, coeffs holds a row of transformed powers of displacement,
  // numVectorCoeffs = SpaceDim * numTaylorCoeffs at a time
  // for each function evaluation, and within each function evaluation,
  // grouped by matrix column, numTaylorCoeffs at a time.
  // - On exit, coeffs is replaced by coefficients,
  // numVectorCoeffs = SpaceDim * numTaylorCoeffs at a time
  // for each function evaluation, and within each function evaluation,
  // grouped by matrix column, numTaylorCoeffs at a time.
  int numFunctionValues = SpaceDim * vectorstencilSize;
  int numVectorCoeffs = SpaceDim * numTaylorCoeffs;
  LAPACKMatrix A(numFunctionValues, numVectorCoeffs);
  {
    int icoeff = 0;
    for (int icell = 0; icell < numFunctionValues; icell++)
      {
        for (int ico = 0; ico < numVectorCoeffs; ico++)
          {
            A(icell, ico) = coeffs[icoeff];
            icoeff++;
          }
      }
  }
  A.pseudoInvertUsingQR();
  {
    int icoeff = 0;
    for (int icell = 0; icell < numFunctionValues; icell++)
      {
        for (int ico = 0; ico < numVectorCoeffs; ico++)
          {
            coeffs[icoeff] = A(ico, icell);
            icoeff++;
          }
      }
  }

  int vectorstencilBaseIndex = 0;
  for (int icell = 0; icell < vectorstencilSize; icell++)
    {
      Tuple<Real, SpaceDim*SpaceDim> wt;
      // From coeffs(numVectorAllCoeffs*icell+(0:numVectorAllCoeffs-1))
      // need to find Tuple<Real, SpaceDim*SpaceDim> wt,
      // and then do a_vectorstencil[icell].setWeight(wt);
      // This is the weight in the interpolation stencil.
      int tupleIndex = 0;
      for (int dst = 0; dst < SpaceDim; dst++)
        {
          for (int src = 0; src < SpaceDim; src++)
            {
              // tupleIndex = dst*SpaceDim + src
              if (a_type == IndexType::TheNodeType())
                { // cell-centered only:
                  // At point Xvec in the ghost cell,
                  // weight of icell in thisVectorStencilCellIndices is
                  // coeffs[icell, pvec] *
                  // ((Xvec() - Xvec(ghostcenter))/avgDistance)^pvec_icell.
                  // At cell center, Xvec() = Xvec(ghostcenter),
                  // so we can ignore the terms with nonzero pvec, and
                  // this is simply coeffs[icell, 0].
                  // In 2D:
                  // wt[0] is weight of src=0 for dst=0
                  // wt[1] is weight of src=1 for dst=0
                  // wt[2] is weight of src=0 for dst=1
                  // wt[3] is weight of src=1 for dst=1
                  wt[tupleIndex] = coeffs[vectorstencilBaseIndex];
                }
              else if (a_type == IndexType::TheCellType())
                { // cell-averaged only:
                  // At point Xvec in the ghost cell,
                  // weight of icell in thisVectorStencilCellIndices is
                  // coeffs[icell, pvec] *
                  // <((Xvec() - Xvec(ghostcenter))/avgDistance)^pvec>_icell
                  // where the average is over the ghost cell.
                  wt[tupleIndex] =
                    dotSubvectors(coeffs, vectorstencilBaseIndex,
                                  thisDispPowers, 0,
                                  numTaylorCoeffs);
                  // If I set wt[tupleIndex] = coeffs[vectorstencilBaseIndex];
                  // then I get 2nd-order convergence.
                }
              tupleIndex++;
              vectorstencilBaseIndex += numTaylorCoeffs;
            } // end loop over source components
        } // end loop over dest components
      a_vectorstencil[icell].setWeight(wt);
      // numVectorAllCoeffs = SpaceDim*SpaceDim * numTaylorCoeffs
      // vectorstencilBaseIndex += numVectorAllCoeffs;
    } // end loop over all cells in the stencil

  // dummy statement in order to get around gdb bug
  int dummy_unused = 0; dummy_unused++;
}


void
MultiBlockUtil::copyStencilsFromBlock(IVSFAB<MBStencil>& a_stencilsFab,
                                      const IntVectSet& a_ghostCellsIVS,
                                      const IVSFAB<MBStencil>& a_blockStencilsFab) const
{
  CH_assert(m_isDefined);
  for (IVSIterator ivsit(a_ghostCellsIVS); ivsit.ok(); ++ivsit)
    {
      // on the collapsed layout
      IntVect thisGhostCell = ivsit();
      IntVect fixedCell(thisGhostCell);
      for (int ind = 0; ind < m_fixedDimsVect.size(); ind++)
        {
          int idir = m_fixedDimsVect[ind];
          fixedCell[idir] = m_fixedPt[ind];
        }
      copyStencilsFromBlock(a_stencilsFab(thisGhostCell, 0),
                            thisGhostCell,
                            a_blockStencilsFab(fixedCell, 0));
    }
}      


void
MultiBlockUtil::copyStencilsFromBlock(MBStencil& a_stencil,
                                      const IntVect& a_iv,
                                      const MBStencil& a_blockStencil) const
{
  // a_stencil will be the same as a_blockStencil except that
  // cell[m_fixedDimsVect] for each stencil element will be set to
  // a_iv[m_fixedDimsVect] instead of m_fixedPt.
  Vector<MBStencilElement> stencilVec;
  MBStencilIterator stencilit(a_blockStencil);
  for (stencilit.begin(); stencilit.ok(); ++stencilit)
    {
      const MBStencilElement& blockStencilElement = stencilit();
      const IntVect& blockStencilCell = blockStencilElement.cell();
      int block = blockStencilElement.block(); // keep this
      Real weight = blockStencilElement.weight(); // keep this
      IntVect thisStencilCell = blockStencilCell; // modify this in m_fixedDims
      for (int idir = m_fixedDims.begin(); idir <= m_fixedDims.end(); idir++)
        {
          thisStencilCell[idir] = a_iv[idir];
        }
      MBStencilElement elt(thisStencilCell, block, weight);
      stencilVec.push_back(elt);
    }
  a_stencil = RefCountedPtr< Vector<MBStencilElement> >
    (new Vector<MBStencilElement>(stencilVec));
}


void
MultiBlockUtil::copyVectorStencilsFromBlock(IVSFAB<MBVectorStencil>& a_vectorstencilsFab,
                                            const IntVectSet& a_ghostCellsIVS,
                                            const IVSFAB<MBVectorStencil>& a_blockVectorStencilsFab) const
{
  CH_assert(m_isDefined);
  for (IVSIterator ivsit(a_ghostCellsIVS); ivsit.ok(); ++ivsit)
    {
      // on the collapsed layout
      IntVect thisGhostCell = ivsit();
      IntVect fixedCell(thisGhostCell);
      for (int ind = 0; ind < m_fixedDimsVect.size(); ind++)
        {
          int idir = m_fixedDimsVect[ind];
          fixedCell[idir] = m_fixedPt[ind];
        }
      copyVectorStencilsFromBlock(a_vectorstencilsFab(thisGhostCell, 0),
                                  thisGhostCell,
                                  a_blockVectorStencilsFab(fixedCell, 0));
    }
}      


void
MultiBlockUtil::copyVectorStencilsFromBlock(MBVectorStencil& a_vectorstencil,
                                            const IntVect& a_iv,
                                            const MBVectorStencil& a_blockVectorStencil) const
{
  // a_vectorstencil will be the same as a_blockVectorStencil except that
  // cell[m_fixedDimsVect] for each vector stencil element will be set to
  // a_iv[m_fixedDimsVect] instead of m_fixedPt.
  Vector<MBVectorStencilElement> vectorstencilVec;
  MBVectorStencilIterator vectorstencilit(a_blockVectorStencil);
  for (vectorstencilit.begin(); vectorstencilit.ok(); ++vectorstencilit)
    {
      const MBVectorStencilElement& blockVectorStencilElement = vectorstencilit();
      const IntVect& blockVectorStencilCell = blockVectorStencilElement.cell();
      int block = blockVectorStencilElement.block(); // keep this
      Tuple<Real, SpaceDim*SpaceDim> weight = blockVectorStencilElement.weight(); // keep this
      IntVect thisVectorStencilCell = blockVectorStencilCell; // modify this in m_fixedDims
      for (int idir = m_fixedDims.begin(); idir <= m_fixedDims.end(); idir++)
        {
          thisVectorStencilCell[idir] = a_iv[idir];
        }
      MBVectorStencilElement elt(thisVectorStencilCell, block, weight);
      vectorstencilVec.push_back(elt);
    }
  a_vectorstencil = RefCountedPtr< Vector<MBVectorStencilElement> >
    (new Vector<MBVectorStencilElement>(vectorstencilVec));
}


Box
MultiBlockUtil::boxFixed(const Box& a_bx)
{
  Box returnBox(a_bx);
  for (int ind = 0; ind < m_fixedDimsVect.size(); ind++)
    {
      int idir = m_fixedDimsVect[ind];
      returnBox.setRange(idir, m_fixedPt[ind]);
    }
  return returnBox;
}


IntVectSet
MultiBlockUtil::extraBlockGhosts(const Box& a_baseBox,
                                 int a_ghostLayer,
                                 int a_baseBlockNum) const
{
  CH_assert(m_isDefined);
  const Vector<Box>& mappingBlocksAll = m_coordSysPtr->mappingBlocks();
  const Box& baseBlockBox = mappingBlocksAll[a_baseBlockNum];

  // Without m_fixedDims: Box grownBox = grow(a_baseBox, a_ghostLayer);
  // Grow box in the interpolating dimensions only.
  IntVect ghostLayerVect = IntVect::Zero;
  for (int ind = 0; ind < m_interpDimsVect.size(); ind++)
    {
      int idir = m_interpDimsVect[ind];
      ghostLayerVect[idir] = a_ghostLayer;
    }
  Box grownBox = grow(a_baseBox, ghostLayerVect);
  
  IntVectSet ivsReturn = IntVectSet(); // empty
  if ( !baseBlockBox.contains(grownBox) ) // otherwise, empty
    {
      // We'll set ivs to extra-block ghost cells of baseBox.
      // Start with ivs being the full grownBox.
      DenseIntVectSet ivs = DenseIntVectSet(grownBox, true);
      // Note baseBlockBox includes all of baseBox,
      // so after this, ivs contains ghost cells only.
      ivs -= baseBlockBox;
      if ( !ivs.isEmpty() )
        {
          // Now remove ghost cells that are outside domain.
          ivs.recalcMinBox();
          // WAS for (int idir = 0; idir < SpaceDim; idir++)
          for (int ind = 0; ind < m_interpDimsVect.size(); ind++)
            {
              int idir = m_interpDimsVect[ind];
              // Check idir Lo face.
              int endBlockLo = baseBlockBox.smallEnd(idir);
              if (grownBox.smallEnd(idir) < endBlockLo)
                { // Remove ghost cells from idir Lo face.
                  if (m_boundaries[a_baseBlockNum][idir].isDomainBoundary())
                    { // Remove ghost cells beyond this face.
                      Box grownBoxFace(grownBox);
                      grownBoxFace.setBig(idir, endBlockLo-1);
                      ivs -= grownBoxFace;
                    }
                }
              // Check idir Hi face.
              int endBlockHi = baseBlockBox.bigEnd(idir);
              if (grownBox.bigEnd(idir) > endBlockHi)
                { // Remove ghost cells from idir Hi face.
                  if (m_boundaries[a_baseBlockNum][idir + SpaceDim].isDomainBoundary())
                    { // Remove ghost cells beyond this face.
                      Box grownBoxFace(grownBox);
                      grownBoxFace.setSmall(idir, endBlockHi+1);
                      ivs -= grownBoxFace;
                    }
                }
            }
          ivs.recalcMinBox();
          // Now ivs is what we want it to be.
          // Find the valid block and valid cell of each ghost
          // cell in ivs.
          // I have to convert from DenseIntVectSet.
          ivsReturn = IntVectSet(ivs);
        }
    }
  return ivsReturn;
}

bool
MultiBlockUtil::allGridsHaveFixedPt(const BoxLayout& a_layout)
{
  if (m_fixedDimsVect.size() == 0)
    {
      return true;
    }
  else
    {
      for (LayoutIterator lit = a_layout.layoutIterator(); lit.ok(); ++lit)
        {
          const Box& bx = a_layout[lit];
          for (int ind = 0; ind < m_fixedDimsVect.size(); ind++)
            {
              int idir = m_fixedDimsVect[ind];
              int val = m_fixedPt[ind];
              if ( ! ( ( bx.smallEnd(idir) <= val ) &&
                       ( val <= bx.bigEnd(idir) ) ) )
                { // val is not in range of bx[idir]
                  return false;
                }
            }
        }
      // every box passed the test in every fixed dimension
      return true;
    }
}


void
MultiBlockUtil::getCollapsedLayout(BoxLayout& a_layoutCollapsed,
                                   const BoxLayout& a_layoutFull)
{
  if (m_fixedDimsVect.size() == 0)
    {
      a_layoutCollapsed = a_layoutFull;
    }
  else
    {
      BoxCollapser collapser(m_fixedDims);
      a_layoutCollapsed.deepCopy(a_layoutFull);
      a_layoutCollapsed.transform(collapser);
      a_layoutCollapsed.closeNoSort();
    }
}

void
MultiBlockUtil::getCollapsedLayout(DisjointBoxLayout& a_layoutCollapsed,
                                   const DisjointBoxLayout& a_layoutFull)
{
  if (m_fixedDimsVect.size() == 0)
    {
      a_layoutCollapsed = a_layoutFull;
    }
  else
    {
      BoxCollapser collapser(m_fixedDims);
      a_layoutCollapsed.deepCopy(a_layoutFull);
      a_layoutCollapsed.transform(collapser);
      a_layoutCollapsed.closeNoSort();
    }
}

void
MultiBlockUtil::getFixedOffLayout(BoxLayout& a_layoutFixedOff,
                                  const BoxLayout& a_layoutFull)
{
  BoxFixedOff fixer(m_fixedDims);
  a_layoutFixedOff.deepCopy(a_layoutFull);
  a_layoutFixedOff.transform(fixer);
  a_layoutFixedOff.closeNoSort();
}

void
MultiBlockUtil::order2grad(LevelData<FArrayBox>&        a_gradData,
                           const LevelData<FArrayBox>&  a_data)
{
  CH_assert(m_isDefined);
  // This function fills in a_gradData on valid cells ONLY.
  const DisjointBoxLayout& layout = a_data.disjointBoxLayout();
  int ncomp = a_data.nComp();
  LevelData<FArrayBox> dataGhosted(layout, ncomp, m_interpUnit);
  DataIterator dit = a_data.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& baseBox = layout[dit];

      const FArrayBox& dataFab = a_data[dit];
      FArrayBox& dataGhostedFab = dataGhosted[dit];

      dataGhostedFab.copy(dataFab);
      // In each interpolating dimension, fill a layer of ghost cells
      // with extrapolant.
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          if ( !m_fixedDims.contains(idir) )
            {
              secondOrderCellExtrap(dataGhostedFab, baseBox, idir);
            }
        }
    }
  // Where it exists, valid data from neighboring boxes overwrites
  // the extrapolated data in dataGhosted.
  // If we have data covering the whole domain, then all extrapolated
  // data will be overwritten.  But the usefulness here comes when
  // we do not have data covering the whole domain.  It works out
  // that we apply the second-order one-sided finite-difference
  // derivative formula on the edge of where we have data.
  dataGhosted.exchange();

  // Now find gradient from dataGhosted.
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& baseBox = layout[dit];
      ProblemDomain blockDomain = m_coordSysPtr->blockDomainOfBox(baseBox);

      const FArrayBox& dataGhostedFab = dataGhosted[dit];
      FArrayBox& gradFab = a_gradData[dit];

      // Don't use any of dataGhostedFab from outside blockDomain.
      Box dataBox(dataGhostedFab.box());
      dataBox &= blockDomain;
      
      // This function is in MBMiscUtil.
      order2gradient(gradFab, baseBox,
                     dataGhostedFab, dataBox,
                     m_interpDimsVect);
    }

  // dummy statement in order to get around gdb bug
  int dummy_unused = 0; dummy_unused++;
}

#include "NamespaceFooter.H"
