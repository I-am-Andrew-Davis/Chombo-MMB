#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MULTIBLOCKCOORDSYS_H_
#define _MULTIBLOCKCOORDSYS_H_

#include "NewCoordSys.H"
#include "Vector.H"
#include "BlockBoundary.H"
#include "Tuple.H"
#include "IndicesTransformation.H"
#include "VectorTransformation.H"
#include "VectorAvgTransformation.H"
#include "LevelData.H"
#include "FluxBox.H"
#include "ComputeCache.H"
#include "NamespaceHeader.H"

/// Additional information cached per block
struct BlockInfo
{
  ProblemDomain m_domain;
  std::string   m_name;
};

/// Methods for evaluating xi->x and caching the results for multiple lookup
/** Use this object as the cache for member function cachedXCoord
 */
using XCoordCache = ComputeCacheUMap<
  stc::IVec<SpaceDim+1>,                      // Key
  RealVect,                                   // Cached value
  std::function<RealVect(const NewCoordSys&,  // CoordSys object
                         const RealVect&)>>;  // Coord in computational space

/// virtual interface class encapsulating multi-block mapping API
/**
   The MultiBlockCoordSys is designed to represent multiblock mappings.
   Each mapping is represented by a separate NewCoordSys which has an
   associated valid region defined by a Box in index space.

   Implementations of specific multiblock coordinate systems will be
   provided by derived classes.
*/
class MultiBlockCoordSys
{

public:

  /// default constructor
  MultiBlockCoordSys();

  /// destructor (should this delete the members of coordSysVect?
  virtual ~MultiBlockCoordSys();

  /// is this multiblock?  Default is yes.  SingleBlockCSAdaptor overrides
  virtual bool isMultiBlock() const
    { return true; }

  /// returns block number in which a_box is contained in the global index space
  /** returns -1 if a_box isn't fully contained in the valid region of
      any block
  */
  virtual int whichBlock(const Box& a_box) const;

  /// returns block number in which a_iv is contained in the global index space
  /** returns -1 if a_iv isn't fully contained in the valid region of
      any block
  */
  virtual int whichBlock(const IntVect& a_iv) const;

  /// returns block number that overlaps a_box in the global index space
  /** returns -1 if either
      a_box doesn't overlap the valid region of any block, or
      a_box overlaps the valid region of more than one block.
  */
  virtual int whichBlockOverlap(const Box& a_box) const;

  /// returns block number in which a_iv is contained within buffer of a_bufferSize
  /** returns -1 if a_iv isn't fully contained in the buffered region of
      any block
  */
  virtual int whichBlockBuffered(const IntVect& a_iv,
                                 int a_bufferSize) const;

  /// returns a block number given a block name
  virtual int whichBlock(const std::string& a_name) const;

  /// block mapping conversion function
  /** given a_xiSrc in block a_nSrc, returns which block (a_n_valid)
      contains that location in its valid region and its
      location a_xi_valid in that block's computational coordinate space.
  */
  virtual void blockRemapping(RealVect&            a_xi_valid,
                              int&                 a_n_valid,
                              bool&                a_validExists,
                              RigidTransformation& a_extraDispl,
                              const RealVect&      a_xiSrc,
                              const IntVect&       a_iSrc,
                              int                  a_nSrc) const
    {
      MayDay::Error("Override this function");
    }
  // This is unused now and will not be called.  But I'm too lazy to change
  // all the CS classes for now.
  virtual void blockRemapping(RealVect&            a_xi_valid,
                              int&                 a_n_valid,
                              bool&                a_validExists,
                              const RealVect&      a_xiSrc,
                              int                  a_nSrc) const
    { }

  // Purpose of this function is that, relative to calling
  // blockRemapping one point at a time, things will run
  // faster if you start by assuming each point is in the same
  // block as the previous point.
  //  virtual void blockRemapping(FArrayBox& a_xi_valid,
  //                              BaseFab<int>& a_n_valid,
  //                              const FArrayBox& a_xiSrc,
  //                              int a_nSrc) const;

  /// general block mapping conversion function
  /** Given a_xiSrc in mapped space of block a_nSrc,
      and block index a_nDst,
      returns coordinates in mapped space of block a_nDst.

      Default implementation goes through real space.
  */
  virtual RealVect blockRemappingGeneral(int a_nDst,
                                         const RealVect& a_xiSrc,
                                         int a_nSrc) const;

  /// return a cell that is far enough away from external boundaries
  /** Given a_cell in block a_blockNum, return another cell that
      is a_cell shifted as necessary away from the external boundaries
      so as to leave distance of at least a_radius to all external boundaries.
  */
  virtual IntVect cellAvoidingBoundaries(const IntVect&   a_cell,
                                         int              a_blockNum,
                                         int              a_radius) const;

  /// modify a_bx so that it lies entirely within the physical domain
  /** Modify a_bx by removing parts of it that lie off edges of
      block a_blockNum that correspond to physical boundaries.
      Cells lying off the block in other directions are retained.
   */
  virtual void keepInDomain(Box& a_bx,
                            int a_blockNum) const;

  /// return ProblemDomain within block, with box and ghost cells
  /** Return non-periodic ProblemDomain with base box
      grow(a_baseBox, a_numGhost), cropped at the physical boundaries.
      Cells lying off the block in other directions are retained.
   */
  //**This is obsolete and replaced by problemDomain(Box) or problemDomain(int)
  virtual ProblemDomain problemDomain(const Box&  a_baseBox,
                                      int         a_numGhost) const;

  /// Returns a problem domain for a block containing a_box
  /** Unlike the above the domain box covers all core cells in the block and
   *  no ghost cells.  Connectivity, which works the same as periodic, is
   *  appropriately marked in the problem domain.
   *  Warning: uses whichBlock which may be expensive.
   */
  const ProblemDomain& problemDomain(const Box& a_box) const;

  /// Returns a problem domain for a block number a_idxBlk
  /** The domain box covers all core cells in the block and no ghost cells.
   *  Connectivity, which works the same as periodic, is appropriately marked
   *  in the problem domain.
   *  Use whichBlock to get a_idxBlk for a box
   */
  const ProblemDomain& problemDomain(const int a_idxBlk) const;

  /// return ProblemDomain of block containing a box
  /** Return non-periodic ProblemDomain consisting of the whole block box
      that contains a_baseBox.
  */
  virtual ProblemDomain blockDomainOfBox(const Box& a_baseBox) const;

  /// A non-periodic problem domain covering the index space of all blocks.
  /** This can be used to construct DisjointBoxLayouts for this level.
   *  Single block adapter should override this to return actual problem domain
   */
  virtual ProblemDomain levelDomain() const;

  /// return vector transformation from one block to another at a point.
  /** Return object that transforms the components of a vector
      situated at a_xiSrc from basis of block a_nSrc
      to components in basis in block a_nDst.

      Default implementation goes through real space.
   */
  virtual VectorTransformation vectorBlockTransformation(int a_nDst,
                                                         const RealVect& a_xiSrc,
                                                         int a_nSrc) const;

  /// return vector transformation from one block to another at center of a cell.
  /** Return object that transforms the components of a vector
      situated at center of cell a_ivSrc from basis of block a_nSrc
      to components in basis in block a_nDst.

      Calls vectorBlockTransformation.
   */
  virtual VectorTransformation vectorBlockTransformationCenter(int a_nDst,
                                                               const IntVect& a_ivSrc,
                                                               int a_nSrc) const;

  virtual void vectorBlockTransformationCenter(BaseFab<VectorTransformation>& a_vtFab,
                                               int a_nDst,
                                               const Box& a_srcBox,
                                               int a_nSrc) const;

  /// return average vector transformation from one block to another at a cell
  /** Return average of object that transforms the components of a vector
      on the cell a_ivSrc from basis of block a_nSrc
      to components in basis in block a_nDst.

      Calls vectorBlockTransformationCenter.
   */
  virtual VectorTransformation vectorBlockTransformationAvg(int a_nDst,
                                                            const IntVect& a_ivSrc,
                                                            int a_nSrc) const;

  /// return average vector transformation from one block to another at a cell, with averaging over specified dimensions only
  /** Return average of object that transforms the components of a vector
      on the cells of a_srcBox from basis of block a_nSrc
      to components in basis in block a_nDst,
      where averaging is done in all dimensions that are NOT in a_fixedDims.

      This is returned in a_vtFab(a_srcBox, 0).

      Calls vectorBlockTransformationCenter.
   */
  virtual void vectorBlockTransformationAvg(BaseFab<VectorTransformation>& a_vtFab,
                                            int a_nDst,
                                            const Box& a_srcBox,
                                            int a_nSrc,
                                            Interval a_fixedDims = Interval() ) const;

  /// return gradient of vector transformation from one block to another at a cell
  /** Return gradient (or more precisely, difference, or gradient multiplied
      by mesh spacing) of object that transforms the components of a vector
      at center of cell a_ivSrc from basis of block a_nSrc
      to components in basis in block a_nDst.

      Calls vectorBlockTransformationCenter.
   */
  virtual Tuple<VectorTransformation, SpaceDim> vectorBlockTransformationGrad(int a_nDst,
                                                                              const IntVect& a_ivSrc,
                                                                              int a_nSrc) const;

  /// return gradient of vector transformation from one block to another at a cell, over specified dimensions only
  /** Return gradient (or more precisely, difference, or gradient multiplied
      by mesh spacing) of object that transforms the components of a vector
      at center of cells in a_srcBox from basis of block a_nSrc
      to components in basis in block a_nDst,
      on all dimensions that are NOT in a_fixedDims.

      This is returned in a_vtFab(a_srcBox, 0:a_interpDimsVect.size()-1).

      Calls vectorBlockTransformationCenter.
   */
  virtual void vectorBlockTransformationGrad(BaseFab<VectorTransformation>& a_vtGradFab,
                                             int a_nDst,
                                             const Box& a_srcBox,
                                             int a_nSrc,
                                             Interval a_fixedDims = Interval() ) const;

  virtual void vectorBlockAvgTransformation(BaseFab<VectorAvgTransformation>& a_vtAvgFab,
                                            int a_nDst,
                                            const Box& a_srcBox,
                                            int a_nSrc,
                                            Interval a_fixedDims = Interval() ) const;

  /// Return physical coordinates given computational coordinates
  /** This routine maintains a cache of coordinates
   *  \param[in]  a_idxCell
   *                      A unique index on the lattice related to the
   *                      computational coordinate.  It is used as part of the
   *                      key for lookup in the cache.
   *  \param[in]  a_idxBlk
   *                      Index of the block.  Used to retrieve the block
   *                      CS and together with 'a_idxCell', this makes up the
   *                      key for lookup in the cache.
   *  \param[in]  a_compCoord
   *                      The coordinate in computational space
   *  \param[in]  a_cache The cache of coordinates in physical space
   *  \param[out] a_cache The cache possibly updated with the new key and
   *                      value.
   *
   *  \note
   *  <ul>
   *    <li> The cache must have been constructed by binding to
   *         NewCoordSys::realCoord(const RealVect&) in order to get
   *         physical-space coordinates from computation-space coordinates
   *  </ul>
   */
  RealVect cachedXCoord(const IntVect&  a_idxCell,
                        const int       a_idxBlk,
                        const RealVect& a_xiCoord,
                        XCoordCache&    a_cache) const
    {
      // Make a key { a_idxCell, a_idxBlk }
      const stc::IVec<SpaceDim+1> key = stc::IVec<SpaceDim+1>(a_idxCell) +
        a_idxBlk*stc::make_IVec<SpaceDim+1>::basis(SpaceDim);
      return a_cache(key, *getCoordSys(a_idxBlk), a_xiCoord);
    }

  /// return displacements
  /**
     Return a Vector of same length as a_dstCoords and a_dstBlocks,
     with component i of this Vector being the displacement
     from point a_srcCoords in mapped coordinates of block a_srcBlock
     to point a_dstCoords[i] in mapped coordinates of block a_dstBlocks[i].

     If this method is not overwritten, then component i of the returned
     Vector is simply the displacement in physical coordinates from
     a_srcCoords to a_dstCoords[i].
   */
  virtual Vector<RealVect> displacements(
    const Vector<RealVect>& a_dstCoords,
    const Vector<int>&      a_dstBlocks,
    const RealVect&         a_srcCoords,
    int                     a_srcBlock) const;

  /// return distance between two points
  virtual Real distance(const RealVect& a_srcCoords,
                        int             a_srcBlock,
                        const RealVect& a_dstCoords,
                        int             a_dstBlock) const;

  /// return distance between two points
  Real distance(const RealVect&            a_srcX,
                const IntVect&             a_dstCell,
                const int                  a_dstBlk,
                const RealVect&            a_dstXi,
                /// Extra displacement going from src to dst
                const RigidTransformation& a_extraDisplX,
                XCoordCache&               a_xCache) const;

  /// Return the mapping which corresponds to the block with the given index.
  virtual const NewCoordSys* getCoordSys(int a_n) const
  {
    return m_coordSysVect[a_n];
  }

  /// Return the mapping which corresponds to the given box.
  /// \param a_box The Box for which the coordinate system is retrieved. This
  ///              Box must be contained completely within a block or this
  ///              method will fail.
  virtual const NewCoordSys* getCoordSys(const Box& a_box) const
  { return m_coordSysVect[whichBlock(a_box)]; }

  /// const access to the mapping blocks
  const Vector<Box>& mappingBlocks() const
  {return m_mappingBlocks;}

  /// number of blocks
  int numBlocks() const
  {return m_mappingBlocks.size();}

  /// const access to block information
  const BlockInfo& blockInfo(const int a_idxBlk) const
    { return m_blockInfo[a_idxBlk]; }

  // /// access to block information
  // BlockInfo& blockInfo(const int a_idxBlk)
  //   { return m_blockInfo[a_idxBlk]; }

  /// const access to the block boundaries
  const Vector< Tuple<BlockBoundary, 2*SpaceDim> >& boundaries() const
  {return m_boundaries; }

  /// const access to a block boundary
  const BlockBoundary& boundary(const int            a_idxBlk,
                                const int            a_dir,
                                const Side::LoHiSide a_side) const
    { return m_boundaries[a_idxBlk][a_side*SpaceDim + a_dir]; }

  /// access to a block boundary
  BlockBoundary& boundary(const int            a_idxBlk,
                          const int            a_dir,
                          const Side::LoHiSide a_side)
    { return m_boundaries[a_idxBlk][a_side*SpaceDim + a_dir]; }

  /// const access to m_validNum[a_block]
  const BaseFab<int>& validNum(int a_block) const
  {return *(m_validNum[a_block]); }

  /// const access to m_validBlocks[a_block]
  const BaseFab<int>& validBlocks(int a_block) const
  {return *(m_validBlocks[a_block]); }

  /// const access to m_validTransformations[a_block]
  const BaseFab<IndicesTransformation>& validTransformations(int a_block) const
  {return *(m_validTransformations[a_block]); }

  /// const access to m_allBlockTransformations[blockSrc][blockDst]
  const IndicesTransformation& blockTransformation(int a_blockSrc,
                                                   int a_blockDst) const
  {return m_allBlockTransformations[a_blockSrc][a_blockDst]; }

  /// Need to have m_boundaries for this.
  void initializeBlockTransformations();

  /// Return list of valid blocks and cells corresponding to a_cell in a_block.
  void validBlocksCells(Vector<int>& a_validBlocks,
                        Vector<IntVect>& a_validCells,
                        int a_origBlock,
                        const IntVect& a_origCell);

  /// Have separate volume flux on each block.  Normally does nothing, but cubed-sphere coordinate system needs it.
  virtual void separateVolFlux(LevelData<FluxBox>& a_flux) const;

  int realDim() const
  {
    return m_realDim;
  }

  RealVect dxVect() const
  {
    return m_dxVect;
  }

  bool gotBoundaries() const
  {
    return m_gotBoundaries;
  }

  bool gotBlockTransformations() const
  {
    return m_gotBlockTransformations;
  }

  bool gotMappingBlocks() const
  {
    return m_gotMappingBlocks;
  }

  bool gotCoordSysVect() const
  {
    return m_gotCoordSysVect;
  }

  bool useScriptN() const
  {
    return m_useScriptN;
  }

  bool isDefined() const;

  void setCache();

  void setVerbosity(int a_verbosity)
    {
      m_verbosity = a_verbosity;
    }

  int verbosity() const
    {
      return m_verbosity;
    }

protected:

  void permute(Vector<int>& a_vec,
               int a_start,
               int a_length,
               Vector< Vector<int> >& a_allPathsDir,
               const Vector<int>& a_changedDirs);

  void setAllBoundaries(BlockBoundary::btype a_btype);

  void setBoundaryFromFaces(int a_srcBlock,
                            int a_srcFaceID,
                            int a_dstBlock,
                            int a_dstFaceID,
                            IntVect a_sign = IntVect::Unit);

  /// Do we have m_coordSysVect?
  bool m_gotCoordSysVect;

  /// Vector of NewCoordSys -- one for each mapping block
  Vector<NewCoordSys*> m_coordSysVect;

  /// Do we have m_mappingBlocks?
  bool m_gotMappingBlocks;

  /// Associated regions in index space over which each mapping is valid
  Vector<Box> m_mappingBlocks;

  /// Have we filled in m_boundaries?
  /// covers m_boundaries
  bool m_gotBoundaries;

  /// The grid spacing, must be the same across all mapping blocks
  RealVect m_dxVect;

  /// Number of dimensions in real space; default SpaceDim
  int m_realDim;

  /// Description of each boundary face of each block; set these in derived class and then set m_gotBoundaries true
  /** These have (C) ordering [blkIdx][side][dir]
   */
  Vector< Tuple<BlockBoundary, 2*SpaceDim> > m_boundaries;

  /// have we filled in transformations in initializeBlockTransformations()?
  /// covers m_allBlockTransformations, m_validNum, m_validBlocks, m_validTransformations
  bool m_gotBlockTransformations;

  /// [blockSrc][blockDst] gives transformation from blockSrc to blockDst
  Vector< Vector<IndicesTransformation> > m_allBlockTransformations;

  /// number of valid blocks containing ghost cells of each block in each direction:
  /// this is the number of components in m_validBlocks and m_validTransformations.
  /// m_validNum lives on [-1:1]^SpaceDim and has 1 component.
  Vector< BaseFab<int>* > m_validNum;

  /// for each block, block numbers of valid blocks containing its ghost cells in each direction.
  /// m_validBlocks[iblock] lives on [-1:1]^SpaceDim and has SpaceDim! components, of which we use m_validNum[iblock](offset, 0) for each offset in [-1:1]^SpaceDim.
  Vector< BaseFab<int>* > m_validBlocks;

  /// for each block, transformation from each block to valid blocks containing its ghost cells in each direction.
  /// m_validTransformation[iblock] lives on [-1:1]^SpaceDim and has SpaceDim! components, of which we use m_validNum[iblock](offset, 0) for each offset in [-1:1]^SpaceDim.
  Vector< BaseFab<IndicesTransformation>* > m_validTransformations;

  /// Whether scriptN is used in calculating N
  //  By default, this is set to true.  If it should be set to false, do so
  //  in the constructors of derived classes
  bool m_useScriptN;

  /// Maps each disjoint box to a block for quick lookup
  // //**FIXME replace with IV lookup
  // std::unordered_map<Box, int, CH_Hash::google_CityHash<Box>> m_boxMap;

  /// Information about a block
  Vector<BlockInfo> m_blockInfo;

  /// How much output do you want?
  int m_verbosity = 0;
};


/// pure-virtual base class defining factory interface for MultiBlockCoordSys
class
MultiBlockCoordSysFactory
{
public:

  MultiBlockCoordSysFactory()
  {
  }

  virtual ~MultiBlockCoordSysFactory()
  {
  }

  virtual MultiBlockCoordSys* getCoordSys(const ProblemDomain& levelDomain,
                                          const RealVect& a_dx) const = 0;
};

#include "NamespaceFooter.H"
#endif
