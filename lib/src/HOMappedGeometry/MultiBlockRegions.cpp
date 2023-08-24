#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


/******************************************************************************/
/**
 * \file MultiBlockRegions.cpp
 *
 * \brief Member functions for MultiBlockRegions
 *
 *//*+*************************************************************************/

//----- Standard Library -----//

//----- Chombo Library -----//

#include "LevelFluxRegister.H"
#include "MultiBlockCoordSys.H"
#include "IndicesTransformation.H"
#include "BlockBoundary.H"
#include "MultiBlockRegions.H"

#define WINLOCKMETHOD 2


/*******************************************************************************
 *
 * Class MultiBlockRegions: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Weak construction
/** Two regions are defined.  The first is for a single level across
 *  MMB connections.  These are denoted by 'lvl' and are used to
 *  average fluxes.  The second is where an AMR interfaced overlaps a
 *  MMB boundary.  These are denoted by 'amr' and are used for face
 *  and flux registers.
 *  \param[in]  a_coordSys
 *                      Multiblock coordinate system
 *  \param[in]  a_CrGrid
 *                      Disjoint box layout for the coarse grid
 *  \param[in]  a_FnGrid
 *                      Disjoint box layout for the fine grid if it
 *                      exists (otherwise set to nullptr)
 *  \param[in]  a_nRef  Refinement ratio between the coarse and fine
 *                      levels
 *
 *  In most cases, codimension 1 of the blocks is a small subset of
 *  the domain.  For same level operations (which may be on the finest
 *  level), compact layouts are created.  For multilevel operations,
 *  the original layouts are retained on the coarse-grid resolution
 *  for storing data, one for the coarse grid and one for a
 *  coarsenening of the fine grid.  For the latter, see
 *  LevelFluxRegister to understand usage.
*//*------------------------------------------------------------------*/

void
MultiBlockRegions::define(const MultiBlockCoordSys&      a_coordSys,
                          const DisjointBoxLayout&       a_CrGrid,
                          const DisjointBoxLayout *const a_FnGrid,
                          const int                      a_nRef)
{
  CH_TIME("MultiBlockRegions::define");
  if (a_coordSys.verbosity() >= 3)
    {
      pout() << "MultiBlockRegions::define " << a_FnGrid << std::endl;
    }
  const int lclProc = procID();  (void)lclProc;

  // lvl -> block adjacency with no resolution change (flux average)
  // amr -> block adjacency with AMR interface (face and flux registers)

  DisjointBoxLayout crFnGrid;
  if (a_FnGrid != nullptr)
    {
      CH_assert(a_nRef > 0);
      coarsen(crFnGrid, *a_FnGrid, a_nRef);
    }

  // Define the locations
  // lvl is used for face averaging
  m_lvlLocations.define(a_CrGrid);
  // amr is used for face and flux registers
  m_crLocations.define(a_CrGrid);
  if (a_FnGrid != nullptr)
    {
      m_fnLocations.define(crFnGrid);
    }
  else
    {
      m_fnLocations.define(BoxLayout(Vector<Box>{}, Vector<int>{}));
    }

  // No grid on the coarse level: clear remaining structures and bail
  if (a_CrGrid.size() == 0)
    {
      CH_assert(a_FnGrid == nullptr);
      clear();
      return;
    }

  // Open the level layout for redefinition
  m_lvlStoLayout = BoxLayout{};

  // Clear AMR structures as required.
  if (a_FnGrid == nullptr)
    {
      m_amrFaceRemoteCopier.clear();
      m_amrFluxRemoteCopier.clear();
    }
  // m_amrCrRcvToOrigin.clear();

  // Motions items for creating the remote copiers
  Vector<RemoteMotion> lvlMotions;
  lvlMotions.reserve(64);
  Vector<RemoteMotion> amrFaceMotions;
  amrFaceMotions.reserve(64);
  Vector<RemoteMotion> amrFluxMotions;
  amrFluxMotions.reserve(64);

/*
  Where AMR interfaces align with block boundaries:
  -------------------------------------------------
  We loop through the coarse boxes and for every box where an AMR interface
  aligns with a block boundary, we identify the coarse cells that will be
  updated because of a flux correction.  We also prepare the copier to post a
  receive to obtain the fine grid information.  But we also need the remote
  processor to know that it has to post the matchine send (the remote processor
  has the fine grid).  There are two ways to find that information, a) run
  another nested loop on the crFnGrid and crGrid and find, from the fine
  perspective, the AMR interfaces that align with block boundaries or b) send a
  message to the remote processor letting it know it needs to post a send.
  Copier performs the former when copying between layouts.  Here we attempt the
  second.  I'm not sure which is best but because we use RMA, the messages are
  intersperesed within a fairly heavy loop and hopefully hidden.  Here, we
  define the counters and buffers necessary to figure out the fine locations and
  the matching sends that need to be posted.

  Where the grid is at the same level on block boundaries:
  --------------------------------------------------------
  Note that this problem does not exist for finding exchanges on the same level
  because the same layout is used for both loops.  If we visit one side of a
  block boundary (and post a send and receive), we know the remote process is
  visiting the same boundary and posting matching send and receives.

  Layouts:
  --------
  Finally, we create compact layouts for averages on the same level.  This is
  because the layout for the fine level is very costly.  For averaging between
  layouts (AMR aligns with BB), data is stored on the coarse level and this is
  tolerable.  The original layouts of the coarse and coarsened-fine are used for
  data (see allocations in LevelFluxRegister)
*/

#ifdef CH_MPI
  int mpierr;
  volatile int* volatile idxBufData;    // For counters into the buffer
  MPI_Win idxBufWin;
  volatile char* volatile bufferData;   // The buffer for receiving data from remotes
  MPI_Win bufferWin;
  // Size of a buffer objects
  const int sizeBoxBytes = linearSize(Box{});
  //                         4 integers            2 boxes
  const int sizeDisplBytes = 4*linearSize(int{}) + 2*sizeBoxBytes;
  // Above is 48 B in 2D and 64 B in 3D
  Vector<char> lclBuffer;   // A buffer for sending
  if (a_FnGrid != nullptr)  // Only need if we have a fine level
    {
      // mpierr = MPI_Alloc_mem(2*sizeof(int), MPI_INFO_NULL, (void*)&idxBufData);
      // if (mpierr)
      //   {
      //     MayDay::Error("MultiBlockRegions::define: MPI_Alloc_mem failed");
      //   }
      mpierr = MPI_Win_allocate(2*sizeof(int), sizeof(int), MPI_INFO_NULL,
                                Chombo_MPI::comm, (void*)&idxBufData,
                                &idxBufWin);
      if (mpierr)
        {
          MayDay::Error("MultiBlockRegions::define: MPI_Win_allocate failed on "
                        "idxBufWin");
        }
      // We get an exclusive lock while initializing the data in idxBufData
      MPI_Win_lock(MPI_LOCK_EXCLUSIVE, lclProc, 0, idxBufWin);
      idxBufData[0] = 0;  // Current index in buffer
//    idxBufData[1] =     // Maximum index in buffer (found next)

/*
  A lot of this is trying to figure out a sufficient buffer size for the RMA
  window.  I think going with a static approach is preferrable and these
  estimates should be quite conservative.  Dynamics windows are supported in
  MPI-3 but it is not particularly straightforward to tell a remote process
  that it needs to allocate more memory
*/

      lclBuffer.resize(sizeDisplBytes);
      // Safety factors for estimating sufficient buffers size
      constexpr int c_safety1 = 2*SpaceDim;  // max face
      constexpr int c_safety2 = 2*(SpaceDim-1);
                                             // max neighbors per face.
                                             // Refinement probably leads to
                                             // less neighbors per face (from
                                             // the fine side)

      // Below assumes each face of each box is on a box boundary and adjacent
      // to c_safety2 neighbors.  /2 assumes adjacency on only half the faces.
      // This permits ~42 boxes in 3-D (512 messages) before going to a full
      // search.
      long long bufferSize = 0;
      constexpr long long bufferDefaultSize = 32768ll;
      if (static_cast<long long>(crFnGrid.dataIterator().size())*c_safety1*
          c_safety2*sizeDisplBytes/2 < bufferDefaultSize)  // Avoid counting
        {
          bufferSize = bufferDefaultSize;
        }
      else
        {
          // Count the faces adjacent to block boundaries
          for (DataIterator ditCrFn(crFnGrid.dataIterator()); ditCrFn.ok();
               ++ditCrFn)
            {
              // Find the block containing this box.
              const Box crFnBox = crFnGrid[ditCrFn];
              const int idxBlk = crFnGrid.blockIndex(ditCrFn);
              CH_assert(idxBlk != -1);
              const Box& blockBox = a_coordSys.mappingBlocks()[idxBlk];
              for (const auto side : EachSide)
                {
                  for (const int dir : EachDir)
                    {
                      if (crFnBox.sideEnd(side)[dir] ==
                          blockBox.sideEnd(side)[dir])
                        {
                          const BlockBoundary& boundary =
                            a_coordSys.boundary(idxBlk, dir, side);
                          if (boundary.isInterface())
                            {
                              ++bufferSize;
                            }
                        }
                    }
                }
            }
          // if (bufferSize != 0)
            {
              bufferSize *= c_safety2*sizeDisplBytes;
              // Additional safety in case of a few pathological faces
              bufferSize = std::max(bufferSize, bufferDefaultSize);
            }
        }
      // This is very important, the total number of locations in the buffer.
      // All processes should compare their given index with this max value and
      // abort if it is exceeded.
      idxBufData[1] = bufferSize/sizeDisplBytes;  // Total items in buffer
      // mpierr = MPI_Win_create((void*)idxBufData, 2*sizeof(int), sizeof(int),
      //                         MPI_INFO_NULL, Chombo_MPI::comm, &idxBufWin);
      // if (mpierr)
      //   {
      //     MayDay::Error("MultiBlockRegions::define: MPI_Win_create failed on "
      //                   "idxBufWin");
      //   }
      MPI_Win_unlock(lclProc, idxBufWin);
      // Ensure idxBufData is initialized for all processors
      mpierr = MPI_Barrier(Chombo_MPI::comm);
      CH_assert(mpierr == MPI_SUCCESS);
      mpierr = MPI_Win_allocate(idxBufData[1]*sizeDisplBytes, sizeDisplBytes,
                                MPI_INFO_NULL, Chombo_MPI::comm, (void*)&bufferData,
                                &bufferWin);
      if (mpierr)
        {
          MayDay::Error("MultiBlockRegions::define: MPI_Win_allocate failed on "
                        "bufferWin");
        }
      // Note, destruction of window will destory bufferData as well in the
      // above
    }  // Have finer level
#endif  /* CH_MPI */

  // Boxes on which we define more compact storage.  FaceTag identifies the
  // origin (a DataIndex in a_CrGrid and a face).
  std::unordered_map<Box, FaceTag, CH_Hash::google_CityHash<Box>> stoBoxes;

  int idxBoFace[2*SpaceDim];            // List of face indices that align with
                                        // block boundaries
  IndicesTransformation transformation[2*SpaceDim];
                                        // Transformation across block
                                        // boundaries
  Box crBox1Tr[2*SpaceDim];             // A layer of cells outside the box,
                                        // transformed into the remote space

/*----------------------------------------------------------------------------*
 * Start the outer coarse loop
 *----------------------------------------------------------------------------*/

#ifdef CH_MPI
  // Windows exposed to RMA
  if (a_FnGrid != nullptr)  // Only need if we have a fine level
    {
#if WINLOCKMETHOD==0
      mpierr = MPI_Win_lock_all(0, idxBufWin);
      CH_assert(mpierr == MPI_SUCCESS);
#elif WINLOCKMETHOD<=1
      mpierr = MPI_Win_lock_all(0, bufferWin);
      CH_assert(mpierr == MPI_SUCCESS);
//#else WINLOCKMETHOD==2 meaning use individual locks everywhere
#endif
    }
#endif  /* CH_MPI */

  for (DataIterator ditCr = a_CrGrid.dataIterator(); ditCr.ok(); ++ditCr)
    {
      // Find the block containing this box.
      const Box crBox = a_CrGrid[ditCr];
      const int idxBlk = a_CrGrid.blockIndex(ditCr);
      CH_assert(idxBlk != -1);
      const Box& blockBox = a_coordSys.mappingBlocks()[idxBlk];
      stc::Vector<Vector<Box>, 2*SpaceDim>& lvlLocations =
        m_lvlLocations[ditCr];
      stc::Vector<Vector<Box>, 2*SpaceDim>& crLocations =
        m_crLocations[ditCr];
      int numBoFace = 0;  // Number of faces on the boundary of the block
      for (const auto side : EachSide)
        {
          for (const int dir : EachDir)
            {
              if (crBox.sideEnd(side)[dir] == blockBox.sideEnd(side)[dir])
                {
                  const BlockBoundary& boundary =
                    a_coordSys.boundary(idxBlk, dir, side);
                  if (boundary.isInterface())
                    {
                      const int idxFace = FaceTag::indexFace(dir, side);
                      idxBoFace[numBoFace++] = idxFace;
                      // The intersection is a single layer of cells external
                      // and adjacent to crBox
                      const Box crBox1 = adjCellBox(crBox, dir, side, 1);
                      // Get the intersection box in the remote space.
                      transformation[idxFace] = boundary.getTransformation();
                      crBox1Tr[idxFace] =
                        transformation[idxFace].transformFwd(crBox1);
                    }
                }
            }
        }

      if (numBoFace == 0) continue;

      // Define a window, outside of which we no longer have to search boxes
      // because of sorting.  If crBox is adjacent to multiple blocks, its
      // transformation could span a large window which may mitigate this
      // optimization.  When crBox is adjacent to only one block, this
      // optimization should be quite effective.
      int small0 = std::numeric_limits<int>::max();
      int big0   = std::numeric_limits<int>::min();
      for (int iBo = 0; iBo != numBoFace; ++iBo)
        {
          const int idxFace = idxBoFace[iBo];
          small0 = std::min(small0, crBox1Tr[idxFace].smallEnd(0));
          big0   = std::max(big0,   crBox1Tr[idxFace].bigEnd(0));
        }

/*----------------------------------------------------------------------------*
 * Start the inner coarse loop
 *----------------------------------------------------------------------------*/

//--Another loop over the coarse layout to determine grid connections across
//--block interfaces

      for (LayoutIterator litCr = a_CrGrid.layoutIterator(); litCr.ok();
           ++litCr)
        {
          const Box crBoxB = a_CrGrid[litCr];
          const int rmtProc = a_CrGrid.procID(litCr());

          // No intersection
          if (crBoxB.bigEnd(0)   < small0) continue;
          // Due to sorting no other boxes intersect
          if (crBoxB.smallEnd(0) > big0) break;

          // Check intersections
          for (int iBo = 0; iBo != numBoFace; ++iBo)
            {
              const int idxFace = idxBoFace[iBo];
              if (crBoxB.intersectsNotEmpty(crBox1Tr[idxFace]))
                {

/*
  The local and remote "storage boxes" need to be canonical, because we need to
  be able to find their definition based on the original DataIndex and the
  location of the face.  One would typically expect that the full face of a box
  to be on the block boundary so this should be reasonable.  These canonical
  boxes are then used to find the DataIndex in the compact layout.
 */

                  // Compact local storage box
                  const int lclDir = FaceTag::dir(idxFace);
                  const Side::LoHiSide lclSide = FaceTag::side(idxFace);
                  // Outer cells on local side
                  const Box crBox1 = adjCellBox(crBox, lclDir, lclSide, 1);
                  // Add inner cells.  This is the local storage box
                  Box lclSto(crBox1);
                  lclSto.growDir(lclDir, Side::flip(lclSide), 1);

                  // Compact remote storage box
                  const IndicesTransformation invTr =
                    transformation[idxFace].inverse();
                  const int rmtDir = invTr.getPermutation()[lclDir];
                  const Side::LoHiSide rmtSide = Side::flip(
                    static_cast<Side::LoHiSide>(
                      ((Side::sign(lclSide)*invTr.getSign()[lclDir]) + 1)/2));
                  // Outer cells on remote side
                  Box rmtSto = adjCellBox(crBoxB, rmtDir, rmtSide, 1);
                  // Add inner cells.  This is the remote storage box
                  rmtSto.growDir(rmtDir, Side::flip(rmtSide), 1);

                  // From hereon, work in local space
                  const Box crBoxBTr =
                    transformation[idxFace].transformBack(crBoxB);

                  // The intersect should be a layer of cells outside crBox
                  const Box outer = crBoxBTr & crBox1;
                  CH_assert(!outer.isEmpty());
                  // Find the layer of cells just inside crBox
                  const Box inner =
                    adjCellBox(outer, lclDir, Side::flip(lclSide), 1);
                  // Note: Local flux is stored to inner and communications fill
                  // outer from the remote side.  Locations refer to inner (the
                  // cells that need to be updated with an averaged flux)
                  lvlLocations[idxFace].push_back(inner);
                  auto ins = stoBoxes.insert(
                    { lclSto, FaceTag(idxFace, ditCr()) });
                  if (!ins.second)  // Storage already exists.  Make sure
                    {               // FaceTag, and DataIndex match
                      CH_assert(ins.first->second.m_idxFace == idxFace);
                      CH_assert(ins.first->second.m_srcDidx == ditCr());
                    }
                  // In the notation here, rmtOuter is the transform of local
                  // outer, not the outer cells on the remote side.
                  Box rmtOuter = transformation[idxFace].transformFwd(outer);
                  Box rmtInner = transformation[idxFace].transformFwd(inner);
                  // Lcl inner ->send-to---> Rmt inner
                  // Lcl outer ->recv-from-> Rmt outer
                  // Note that lclSto and rmtSto are not used for communication,
                  // but saving them here allows for finding the DataIndex on
                  // the compact layout which is used for storage.
                  //                     Storage      Send   Receive
                  lvlMotions.emplace_back(lclSto,    inner,    outer,  // Local
                                          rmtSto, rmtOuter, rmtInner,  // Remote
                                          rmtProc);
                }
            }
        }  // Loop over all boxes on coarse layout

/*----------------------------------------------------------------------------*
 * Start the inner coarsened-fine loop
 *----------------------------------------------------------------------------*/

//--A loop over the fine layout to determine grid connections across
//--block interfaces.

      if (a_FnGrid != nullptr)
        {
          for (LayoutIterator litCrFn = crFnGrid.layoutIterator(); litCrFn.ok();
               ++litCrFn)
            {
              const Box crFnBox = crFnGrid[litCrFn];
              const int rmtProc = crFnGrid.procID(litCrFn());

              // No intersection
              if (crFnBox.bigEnd(0)   < small0) continue;
              // Due to sorting no other boxes intersect
              if (crFnBox.smallEnd(0) > big0) break;

              // Check intersections
              for (int iBo = 0; iBo != numBoFace; ++iBo)
                {
                  const int idxFace = idxBoFace[iBo];
                  if (crFnBox.intersectsNotEmpty(crBox1Tr[idxFace]))
                    {
                      const int lclDir = FaceTag::dir(idxFace);
                      const Side::LoHiSide lclSide = FaceTag::side(idxFace);
                      // Outer cells on local side
                      const Box crBox1 = adjCellBox(crBox, lclDir, lclSide, 1);

                      // Get the remote face direction and side
                      const IndicesTransformation invTr =
                        transformation[idxFace].inverse();
                      const int rmtDir = invTr.getPermutation()[lclDir];
                      const Side::LoHiSide rmtSide = Side::flip(
                        Side::LoHiSide(
                          (Side::sign(lclSide)*invTr.getSign()[lclDir] + 1)/2));
                      const int rmtIdxFace =
                        FaceTag::indexFace(rmtDir, rmtSide);

                      // From hereon, work in local space
                      const Box crFnBoxTr =
                        transformation[idxFace].transformBack(crFnBox);

                      // The intersect should be a layer of cells outside crBox
                      const Box outer = crFnBoxTr & crBox1;
                      CH_assert(!outer.isEmpty());
                      // Find the layer of cells just inside crBox
                      const Box inner =
                        adjCellBox(outer, lclDir, Side::flip(lclSide), 1);

//**FIXME This doesn't work for FluxRegister.  There, we receive into the inner
//**      cells for a direct update.  Using inner for a flux, versus a flux
//**      correction, is not unique for all faces and directions.  Probably need
//**      both amrOuter and amrInner copiers
                      // Local flux correction is stored to inner and
                      // communications fill outer from the remote side.
                      // Locations refer to inner (the cells that need to be
                      // updated with a corrected flux).
                      // WARNING! In LevelFluxRegister, the sides of crLocations
                      // are flipped.  E.g. the high side in a crLocations means
                      // the high side of the fine box, which is the low side of
                      // the coarse box.  This is taken into account when using
                      // addCrLocationsToFluxReg
                      crLocations[idxFace].push_back(inner);
                      // In the notation here, rmtOuter is the transform of
                      // local outer, not the outer cells on the remote side.
                      // Box rmtOuter =
                      //   transformation[idxFace].transformFwd(outer);
                      Box rmtInner =
                        transformation[idxFace].transformFwd(inner);
                      // Face motions are from outer cells to outer cells.
                      // I.e., the cells do not overlap in physical space.
                      // Lcl outer ->recv-from-> Rmt inner
                      //                        Storage      Send Receive
                      amrFaceMotions.emplace_back(RemoteMotion::receive_only{},
                                                  Box{},           outer, // Lcl
                                                  Box{}, rmtInner,        // Rmt
                                                  rmtProc);
                      RemoteMotion& amrFaceMotion = amrFaceMotions.back();
                      amrFaceMotion.m_lclStoDidx = ditCr();
                      amrFaceMotion.m_rmtStoDidx = DataIndex(litCrFn());
                      // Flux motions are from outer cells to inner cells.
                      // I.e., the cells overlap in physical space.
                      // Lcl outer ->recv-from-> Rmt outer
                      //                        Storage      Send Receive
                      amrFluxMotions.emplace_back(RemoteMotion::receive_only{},
                                                  Box{},           inner, // Lcl
                                                  Box{}, rmtInner,        // Rmt
                                                  rmtProc);
                      RemoteMotion& amrFluxMotion = amrFluxMotions.back();
                      amrFluxMotion.m_lclStoDidx = ditCr();
                      amrFluxMotion.m_rmtStoDidx = DataIndex(litCrFn());
                      // // Outer is the receive or destination location from a
                      // // copyTo.  For multiblock, it is also unique.  Relate
                      // // it to the DataIndex, face, and side of the coarse
                      // // box
                      // auto ins = m_amrCrRcvToOrigin.insert(
                      //   { outer, FaceTag(idxFace, ditCr()) });
                      // // Insertion must have succeeded or this is not unique
                      // CH_assert(ins.second);
                      // If local transfer, save the fine locations.
#ifdef CH_MPI
                      if (rmtProc == lclProc)
#endif
                        {
                          Box fnBox = a_FnGrid->operator[](litCrFn);
                          // Find the layer of cells just outside fnBox
                          fnBox = adjCellBox(fnBox, rmtDir, rmtSide, 1);
                          // Remember rmtInner is outer cells on remote side
                          Box fnOuter = refine(rmtInner, a_nRef);
                          fnOuter &= fnBox;
                          CH_assert(!fnOuter.isEmpty());
                          m_fnLocations[DataIndex(litCrFn())][rmtIdxFace]
                            .push_back(fnOuter);
                        }

//--Receives also cover local copy.  However, if using MPI, we need to inform
//--the remote process that it must send a message.

#ifdef CH_MPI
                      if (rmtProc != lclProc)
                        {
                          int idxBuf;
                          int maxBuf;
                          int one = 1;
                          // Get a location to store data into the remote.
                          // This is atomic and thread safe
#if WINLOCKMETHOD>=1
                          mpierr = MPI_Win_lock(MPI_LOCK_SHARED, rmtProc, 0,
                                                idxBufWin);
                          CH_assert(mpierr == MPI_SUCCESS);
#endif
                          mpierr = MPI_Get(&maxBuf, 1, MPI_INT, rmtProc, 1, 1,
                                           MPI_INT, idxBufWin);
                          CH_assert(mpierr == MPI_SUCCESS);
                          mpierr = MPI_Fetch_and_op(&one, &idxBuf, MPI_INT,
                                                    rmtProc, 0, MPI_SUM,
                                                    idxBufWin);
                          CH_assert(mpierr == MPI_SUCCESS);
#if WINLOCKMETHOD==0
                          mpierr = MPI_Win_flush(rmtProc, idxBufWin);
                          if (mpierr)
                            {
                              MayDay::Error("MultiBlockRegions::define: "
                                            "MPI_Win_flush_local failed on "
                                            "idxBufWin");
                            }
#elif WINLOCKMETHOD>=1
                          mpierr = MPI_Win_unlock(rmtProc, idxBufWin);
                          if (mpierr)
                            {
                              MayDay::Error("MultiBlockRegions::define: "
                                            "MPI_Win_unlock failed on "
                                            "idxBufWin");
                            }
#endif
                          if (idxBuf >= maxBuf)  // Insufficient mem on remote
                            {
                              // This isn't a run out of memory error, rather
                              // we did not budget enough memory for faces.
                              pout() << "MBR: ERROR: Save from " << lclProc
                                     << " to " << rmtProc << " with index: "
                                     << idxBuf << " and max: " << maxBuf
                                     << std::endl;
                              MayDay::Error("MultiBlockRegions::define: "
                                            "insufficient memory on remote to "
                                            "save send information");
                            }
                          // Make sure the local buffer is free
#if WINLOCKMETHOD<=1
                          mpierr = MPI_Win_flush_local_all(bufferWin);
                          CH_assert(mpierr == MPI_SUCCESS);
#else  /* WINLOCKMETHOD=2 */
                          mpierr = MPI_Win_lock(MPI_LOCK_EXCLUSIVE, rmtProc, 0,
                                                bufferWin);
                          CH_assert(mpierr == MPI_SUCCESS);
#endif
                          // Pack the information into a local buffer
                          char *b = lclBuffer.data();
                          linearOut(b, litCrFn().intCode()); // Rmt layout index
                          b += sizeof(int);
                          linearOut(b, ditCr().intCode());   // Lcl layout index
                          b += sizeof(int);
                          linearOut(b, rmtIdxFace);          // Rmt face index
                          b += sizeof(int);
                          linearOut(b, idxFace);             // Lcl face index
                          b += sizeof(int);
                          linearOut(b, rmtInner);            // Rmt inner box
                          b += sizeBoxBytes;
                          linearOut(b, inner);               // Lcl inner box
                          b += sizeBoxBytes;
                          mpierr = MPI_Put(lclBuffer.data(), sizeDisplBytes,
                                           MPI_BYTE, rmtProc, idxBuf,
                                           sizeDisplBytes, MPI_BYTE, bufferWin);
                          CH_assert(mpierr == MPI_SUCCESS);
#if WINLOCKMETHOD==2
                          mpierr = MPI_Win_unlock(rmtProc, bufferWin);
                          CH_assert(mpierr == MPI_SUCCESS);
#endif
                        }
#endif  /* CH_MPI */
                    }  // Have intersection with crFnBox
                }  // Loop over faces adjacent to a block boundary
            }  // Loop over all boxes onn coarsened-fine layout
        }  // Have a finer level
    }  // Loop over local boxes on coarse level

/*----------------------------------------------------------------------------*
 * Finish construction
 *----------------------------------------------------------------------------*/

//--Build the compact layout for same-level operation

  m_lvlStoLayout.define<FaceTag>(stoBoxes,
                                 &m_lvlSto2Origin,
                                 &m_lvlMapOrigin2Sto);
  // We need the local blocks on the layout
  m_lvlStoLayout.defineLocalBlocks(a_coordSys.mappingBlocks(), false);

  // Populate the DataIndex for the motion items on the same level
  {
    // Relate boxes to DataIndex.  We need this because we are going to do a
    // random lookup by box.
    std::unordered_map<Box, DataIndex, CH_Hash::google_CityHash<Box>> stoDidx;
    stoDidx.max_load_factor(0.7);
    stoDidx.reserve(m_lvlStoLayout.size());
    for (LayoutIterator lit = m_lvlStoLayout.layoutIterator(); lit.ok(); ++lit)
      {
        auto ins = stoDidx.insert({ m_lvlStoLayout[lit], DataIndex(lit()) });
        CH_assert(ins.second);  // The boxes must all be unique
      }
    // Find DataIndex for the boxes in the remote motions
    for (RemoteMotion& lvlMotion : lvlMotions)
      {
        const auto lclIter = stoDidx.find(lvlMotion.m_lclSto);
        CH_assert(lclIter != stoDidx.end());
        lvlMotion.m_lclStoDidx = lclIter->second;
        const auto rmtIter = stoDidx.find(lvlMotion.m_rmtSto);
        CH_assert(rmtIter != stoDidx.end());
        lvlMotion.m_rmtStoDidx = rmtIter->second;
      }
  }

//--Build the remote copier for the level

  m_lvlRemoteCopier.define(lvlMotions);

//--Create motion items for sends from local fine level to remote coarse level
//--These have been placed in our buffer by remote processes

#ifdef CH_MPI
  if (a_FnGrid != nullptr)  // Only need if we have a fine level
    {
      // Complete RMA in all windows (and synchronize in time)
#if WINLOCKMETHOD==0
      mpierr = MPI_Win_unlock_all(idxBufWin);
      if (mpierr)
        {
          MayDay::Error("MultiBlockRegions::define: MPI_Win_unlock_all failed "
                        "on idxBufWin");
        }
#elif WINLOCKMETHOD<=1
      mpierr = MPI_Win_unlock_all(bufferWin);
      if (mpierr)
        {
          MayDay::Error("MultiBlockRegions::define: MPI_Win_unlock_all failed "
                        "on bufferWin");
        }
//#else WINLOCKMETHOD==2 meaning use individual locks everywhere
#endif
      // Make sure unlocks from remote procs are complete.  I.e., we don't want
      // this proc to lock before MPI_Win_lock_all (which is not collective!)
      mpierr = MPI_Barrier(Chombo_MPI::comm);
      CH_assert(mpierr == MPI_SUCCESS);
      // Start epoch for local access
      mpierr = MPI_Win_lock(MPI_LOCK_EXCLUSIVE, lclProc, 0, bufferWin);
      CH_assert(mpierr == MPI_SUCCESS);

      const int numSend = idxBufData[0];
      if (numSend > idxBufData[1])  // Insufficient local memory
        {
          // This isn't a run out of memory error, rather we did not budget
          // enough memory for faces.
          MayDay::Error("MultiBlockRegions::define: insufficient memory on "
                        "local for saved send information");
        }
      volatile char *b = bufferData;
      for (int idxSend = 0; idxSend != numSend; ++idxSend)
        {
          int lclIntCode, rmtIntCode, idxFace, rmtIdxFace;
          Box outer, rmtOuter;
          linearIn(lclIntCode, (char*)b);  // Lcl (coarsened-fine) layout int code
          b += sizeof(int);
          linearIn(rmtIntCode, (char*)b);  // Rmt (coarse) layout int code
          b += sizeof(int);
          linearIn(idxFace, (char*)b);     // Lcl face index
          b += sizeof(int);
          linearIn(rmtIdxFace, (char*)b);  // Rmt face index
          b += sizeof(int);
          linearIn(outer, (char*)b);       // Lcl outer box
          b += sizeBoxBytes;
          linearIn(rmtOuter, (char*)b);    // Rmt outer box
          b += sizeBoxBytes;
          DataIndex didxCrFn = DataIndex(crFnGrid.layoutIterator()[lclIntCode]);
          DataIndex didxCr   = DataIndex(a_CrGrid.layoutIterator()[rmtIntCode]);
          const int rmtProc = a_CrGrid.procID(didxCr);
          CH_assert(rmtProc != lclProc);
          // In the notation here, rmtInner is the transform of local inner, not
          // the inner cells on the remote side.  So rmtInner is the outer cells
          // on the remote side.
          const Box rmtInner = adjCellBox(rmtOuter,
                                          FaceTag::dir(rmtIdxFace),
                                          FaceTag::side(rmtIdxFace),
                                          1);
          // Face motions are from outer cells to outer cells.
          // I.e., the cells do not overlap in physical space.
          // Lcl outer ->recv-from-> Rmt inner
          //                        Storage   Send   Receive
          amrFaceMotions.emplace_back(RemoteMotion::send_only{},
                                      Box{}, outer,           // Local
                                      Box{},        rmtInner, // Remote
                                      rmtProc);
          RemoteMotion& amrFaceMotion = amrFaceMotions.back();
          amrFaceMotion.m_lclStoDidx = didxCrFn;
          amrFaceMotion.m_rmtStoDidx = didxCr;
          // Flux motions are from outer cells to inner cells.
          // I.e., the cells overlap in physical space.
          // Lcl outer ->recv-from-> Rmt outer
          //                        Storage   Send   Receive
          amrFluxMotions.emplace_back(RemoteMotion::send_only{},
                                      Box{}, outer,           // Local
                                      Box{},        rmtOuter, // Remote
                                      rmtProc);
          RemoteMotion& amrFluxMotion = amrFluxMotions.back();
          amrFluxMotion.m_lclStoDidx = didxCrFn;
          amrFluxMotion.m_rmtStoDidx = didxCr;
          Box fnBox = a_FnGrid->operator[](didxCrFn);
          // Find the layer of cells just outside fnBox
          const int lclDir = FaceTag::dir(idxFace);
          const Side::LoHiSide lclSide = FaceTag::side(idxFace);
          fnBox = adjCellBox(fnBox, lclDir, lclSide, 1);
          outer.refine(a_nRef);
          outer &= fnBox;
          CH_assert(!outer.isEmpty());
          m_fnLocations[didxCrFn][idxFace].push_back(outer);
        }
      // End epoch for local access
      mpierr = MPI_Win_unlock(lclProc, bufferWin);
      CH_assert(mpierr == MPI_SUCCESS);

      // Release buffers and window
      mpierr = MPI_Win_free(&idxBufWin);
      CH_assert(mpierr == MPI_SUCCESS);
      // mpierr = MPI_Free_mem((void*)idxBufData);
      // CH_assert(mpierr == MPI_SUCCESS);
      mpierr = MPI_Win_free(&bufferWin);
      CH_assert(mpierr == MPI_SUCCESS);
      // Do not free bufferData since bufferWin was created with
      // MPI_Win_allocate
    }  // Have finer level
#endif  /* CH_MPI */

//--Build the remote copier for AMR

  if (a_FnGrid != nullptr)  // Only needed if we have a fine level
    {
      // pout() << "Fn Loc: \n";
      // for (DataIterator dit(*a_FnGrid); dit.ok(); ++dit)
      //   {
      //     const Box& disjointBox = a_FnGrid->operator[](dit);
      //     for (const int dir : EachDir)
      //       {
      //         for (const auto side : EachSide)
      //           {
      //             pout() << "  " << disjointBox << " dir" << dir << ",side" << side << ":\n";
      //             const Vector<Box>& locs = getFnLocations(dit, dir, side);
      //             for (const auto& box : locs)
      //               {
      //                 pout() << "    " << box << std::endl;
      //               }
      //           }
      //       }
      //   }
      // pout() << "FaceMot:\n";
      // for (const auto& i : amrFaceMotions)
      //   {
      //     pout() << i.m_lclSto << " " << i.m_lclSnd << " " << i.m_lclRcv << " "
      //            << i.m_rmtSto << " " << i.m_rmtSnd << " " << imm_rmtRcv
      //            << std::endl;
      //   }
      m_amrFaceRemoteCopier.define(amrFaceMotions);
      m_amrFluxRemoteCopier.define(amrFluxMotions);
    }
}

/*--------------------------------------------------------------------*/
//  Clear the data structures except locations
/** 
 *//*-----------------------------------------------------------------*/

void
MultiBlockRegions::clear()
{
  m_lvlStoLayout = BoxLayout(Vector<Box>{}, Vector<int>{});
  m_lvlSto2Origin.define(m_lvlStoLayout);
  m_lvlMapOrigin2Sto.clear();
  m_lvlRemoteCopier.clear();
  m_amrFaceRemoteCopier.clear();
  m_amrFluxRemoteCopier.clear();
  // m_amrCrRcvToOrigin.clear();
}

/*--------------------------------------------------------------------*/
//  Add coarse locations from block boundaries to the
//  LevelFluxRegister
/** \param[in]  a_fluxReg
 *                      The level flux register already fully defined
 *                      for intra-block corrections.
 *  \param[out] a_fluxReg
 *                      Coarse locations updated for inter-block
 *                      corrections
 *//*-----------------------------------------------------------------*/

void
MultiBlockRegions::addCrLocationsToFluxReg(LevelFluxRegister& a_fluxReg) const
{
  CH_assert(a_fluxReg.isAllDefined());
  for (DataIterator dit(m_lvlLocations.dataIterator()); dit.ok(); ++dit)
    {
      const stc::Vector<Vector<Box>, 2*SpaceDim>& crLocations =
        m_crLocations[dit];
      for (int idxFace = 0; idxFace != 2*SpaceDim; ++idxFace)
        {
          if (crLocations[idxFace].size() > 0)
            {
              // Find direction and side.  We need to reverse the side for
              // the flux register because conventions there is that the
              // interface is on the side of the fine mesh.  Don't worry about
              // that actual orientation of the remote face, that is taken care
              // of during the copy.
              const int dir = FaceTag::dir(idxFace);
              const Side::LoHiSide side = FaceTag::side(idxFace);
              Vector<Box>& flxLvlCrLocations =
                a_fluxReg.getCoarseLocations(dir, Side::flip(side))[dit];
              for (const Box& box : crLocations[idxFace])
                {
                  flxLvlCrLocations.push_back(box);
                }
            }
        }
    }
}
