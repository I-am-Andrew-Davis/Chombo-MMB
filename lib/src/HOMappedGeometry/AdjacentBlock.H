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
 * \file AdjacentBlock.H
 *
 * \brief Technique to visit adjacent block for some task
 *
 *//*+*************************************************************************/

#ifndef _ADJACENTBLOCK_H_
#define _ADJACENTBLOCK_H_

//----- Chombo Library -----//

#include "Box.H"
#include "MultiBlockCoordSys.H"
#include "IndicesTransformation.H"

#include "NamespaceHeader.H"

namespace AdjacentBlock
{

  constexpr int c_vecSize = 16;       ///< Size of reservation for new vectors


/*******************************************************************************
 */
///  AdjacentBlockApply
/**
 *//*+*************************************************************************/

template <typename F>
struct AdjacentBlockApply
{
  static void op(const MultiBlockCoordSys&    a_coordSys,
                 const IndicesTransformation& a_hereTrfm,
                 const Box&                   a_box,
                 const Box&                   a_nodeBox,
                 const int                    a_idxBlk,
                 F&&                          a_op,
                 const IntVect&               a_dirs,
                 Vector<int>&                 a_history,
                 int                          a_hops)
    {
      if (++a_hops > 64)
        {
          MayDay::Error("Improbable number of hops in AdjacentBlockApply::op.  "
                        "Recursive chain terminated.");
        }
//**  You may want to return if this block connects to the original block across
//**  a face and a_hops > 1.  But my guess is that we do want access to fully
//**  wrap around edges and corners, even for a Cartesian arrangement of blocks.
      a_history.push_back(a_idxBlk);
      // Save the current history and size
      const Vector<int> hereHistory = a_history;
      const int hereHistorySize = a_history.size();
      // Anchor for a new history chain.  This is only used for new directions.
      Vector<int> anchorHistoryChain;

      const Box blockBox = a_coordSys.mappingBlocks()[a_idxBlk];
      const Box opBox = a_box & blockBox;
      // Execute the operator
      a_op(opBox, a_idxBlk, a_hereTrfm);
      if (opBox == a_box && hereHistorySize == 2)
        {
          // We hopped across a face ond there is nothing left to transfer
          return;
        }
      // else, if history.size() > 2 we shift opBox into the next domain and
      // try to visit all blocks somehow connected to the nodeBox

      // Transfer to adjacent blocks
      for (const int dir : EachDir)
        {
          for (const auto side : EachSide)
            {
              int idxNbrBlk;   // Indexof neighbor block
              Box blockNodes;  // Nodes in the block
              Box box;      // This will be the box transferred to the neighbor
              Box nodeBox;  // This is the nodes transferred to the neighbor
              IndicesTransformation nextTrfm;       // From this to neighbor
              IndicesTransformation compositeTrfm;  // All
              IntVect dirs;              // Visited dirs passed to neighbor
              Vector<int>* nextHistory;  // History passed to neighbor
              
              const BlockBoundary& boundary =
                a_coordSys.boundary(a_idxBlk, dir, side);
              if (!boundary.isInterface()) goto checkNextBlock;
              idxNbrBlk = boundary.neighbor();
              blockNodes = a_coordSys.mappingBlocks()[idxNbrBlk];
              // If this is a direction that has been visited, do not revisit
              // neighbors.  You can revisit neighbors in new directions
              // because these might be connections to the same block - think
              // of a singe block that is fully periodic.
              if (a_dirs[dir] != 0)
                {
                  for (const int jNbr : a_history)
                    {
                      if (jNbr == idxNbrBlk) goto checkNextBlock;
                    }
                }
              // To continue, the nodes of the neighbor block must intersect
              // the nodeBox
              nextTrfm = boundary.getTransformation();
              nodeBox = nextTrfm.transformFwd(a_nodeBox);
              blockNodes.surroundingNodes();
              if (!blockNodes.intersectsNotEmpty(nodeBox)) goto checkNextBlock;

//If a new direction, use pure history
// Create and use new copy of hereHistory

//If an old direction, use updated history
// Use a_history

              if (a_dirs[dir] == 0)
                // If this is a newly visited direction, transfer leftovers of
                // (a_box - opBox).  Anchor a new history chain.
                {
                  box = a_box;
                  // Crop the part of a_box that is inside the current domain
                  switch (side)
                    {
                    case Side::Lo:
                      if (blockBox.smallEnd(dir) - 1 < box.smallEnd(dir))
                        {
                          // Box would be empty
                          goto checkNextBlock;
                        }
                      box.setBig(dir, blockBox.smallEnd(dir) - 1);
                      break;
                    default:
                      if (blockBox.bigEnd(dir) + 1 > box.bigEnd(dir))
                        {
                          // Box would be empty
                          goto checkNextBlock;
                        }
                      box.setSmall(dir, blockBox.bigEnd(dir) + 1);
                      break;
                    }
                  anchorHistoryChain = hereHistory;
                  anchorHistoryChain.reserve(c_vecSize);
                  nextHistory = &anchorHistoryChain;
                }
              else
                // Shift the remainder of opBox into the neighbor block.  This
                // must be accomplished by a transformation because we always
                // need to be able to related back to the origin block.  History
                // is accumulated in current chain.
                {
                  box = opBox;
                  IntVect translate = IntVect_zero;
                  translate[dir] = Side::sign(side)*opBox.size(dir);
                  IndicesTransformation transTrfm;
                  transTrfm.defineFromTranslation(translate);
                  nextTrfm = transTrfm.compose(boundary.getTransformation());
                  nextHistory = &a_history;
                }

              box = nextTrfm.transformFwd(box);
              compositeTrfm = a_hereTrfm.compose(nextTrfm);
              dirs = a_dirs;
              dirs[dir] = 1;
              dirs = nextTrfm.transformVectorFwd(dirs);
              AdjacentBlockApply::op(a_coordSys,
                                     compositeTrfm,
                                     box,
                                     nodeBox,
                                     idxNbrBlk,
                                     std::forward<F>(a_op),
                                     dirs,
                                     *nextHistory,
                                     a_hops);
              checkNextBlock : ;
            }
        }
    }
};

/*--------------------------------------------------------------------*/
///
/**
 *//*-----------------------------------------------------------------*/

template <typename F>
void
hereAndForEachApply(const MultiBlockCoordSys& a_coordSys,
                    const Box&                a_box,
                    const int                 a_idxBlk,
                    F&&                       a_op)
{
  // Apply op to the current block and then go to adjacent blocks
  const Box blockBox = a_coordSys.mappingBlocks()[a_idxBlk];
  const Box opBox = a_box & blockBox;
  // Execute the operator
  a_op(opBox, a_idxBlk, IndicesTransformation::Identity);
  if (opBox == a_box)
    {
      // No cells left to transfer
      return;
    }
  forEachApply(a_coordSys,
               a_box,
               a_idxBlk,
               std::forward<F>(a_op));
}

/*--------------------------------------------------------------------*/
///
/**
 *//*-----------------------------------------------------------------*/

template <typename F>
void
forEachApply(const MultiBlockCoordSys& a_coordSys,
             const Box&                a_box,
             const int                 a_idxBlk,
             F&&                       a_op)
{
  CH_TIME("AdjacentBlock:forEachApply");
  const Box blockBox = a_coordSys.mappingBlocks()[a_idxBlk];
  if (blockBox.contains(a_box)) return;
  Vector<int> history(c_vecSize);
  for (const int dir : EachDir)
    {
      for (const auto side : EachSide)
        {
          if (Side::sign(side)*(blockBox.sideEnd(side)[dir] -
                                a_box.sideEnd(side)[dir]) < 0)
            {
              Box box = a_box;
              // Crop the part of a_box that is inside the current domain
              switch (side)
                {
                case Side::Lo:
                  box.setBig(dir, blockBox.smallEnd(dir) - 1);
                  break;
                default:
                  box.setSmall(dir, blockBox.bigEnd(dir) + 1);
                  break;
                }
              // The node box.  Only adjacent blocks that intersect this are of
              // interest.
              Box nodeBox = bdryBox(box, dir, Side::flip(side), 1);
              nodeBox.surroundingNodes();
              // Get the boundary and transform both boxes
              const BlockBoundary& boundary =
                a_coordSys.boundary(a_idxBlk, dir, side);
              IndicesTransformation nextTrfm = boundary.getTransformation();
              box = nextTrfm.transformFwd(box);
              nodeBox = nextTrfm.transformFwd(nodeBox);
              const int idxNbrBlk = boundary.neighbor();
              // Note that each face is treated as a new chain
              history.clear();
              history.push_back(a_idxBlk);
              IntVect dirs = IntVect_basis(dir);
              dirs = nextTrfm.transformVectorFwd(dirs);
              AdjacentBlockApply<F>::op(a_coordSys,
                                        nextTrfm,
                                        box,
                                        nodeBox,
                                        idxNbrBlk,
                                        std::forward<F>(a_op),
                                        dirs,
                                        history,
                                        0);
            }
        }
    }
}

/*--------------------------------------------------------------------*/
///
/**
 *//*-----------------------------------------------------------------*/

template <typename F>
void
forEachApply(const MultiBlockCoordSys& a_coordSys,
             const Box&                a_box,
             const int                 a_idxBlk,
             const int                 a_dir,
             const Side::LoHiSide      a_side,
             F&&                       a_op)
{
  const Box blockBox = a_coordSys.mappingBlocks()[a_idxBlk];
  Box box = a_box;
  // Crop the part of a_box that is inside the current domain
  switch (a_side)
    {
    case Side::Lo:
      box.setBig(a_dir, blockBox.smallEnd(a_dir) - 1);
      break;
    default:
      box.setSmall(a_dir, blockBox.bigEnd(a_dir) + 1);
      break;
    }
  if (box.isEmpty())
    {
      // a_box did not extend across the argument face
      return;
    }
  // The node box.  Only adjacent blocks that intersect this are of interest.
  Box nodeBox = bdryBox(box, a_dir, Side::flip(a_side), 1);
  nodeBox.surroundingNodes();
  // Get the boundary and transform both boxes
  const BlockBoundary& boundary = a_coordSys.boundary(a_idxBlk, a_dir, a_side);
  IndicesTransformation nextTrfm = boundary.getTransformation();
  box = nextTrfm.transformFwd(box);
  nodeBox = nextTrfm.transformFwd(nodeBox);
  const int idxNbrBlk = boundary.neighbor();
  Vector<int> history(c_vecSize);
  history.clear();
  history.push_back(a_idxBlk);
  IntVect dirs = IntVect_basis(a_dir);
  dirs = nextTrfm.transformVectorFwd(dirs);
  int hops = 0;
  AdjacentBlockApply<F>::op(a_coordSys,
                            nextTrfm,
                            box,
                            nodeBox,
                            idxNbrBlk,
                            std::forward<F>(a_op),
                            dirs,
                            history,
                            0);
}

}  // namespace AdjacentBlock

#include "NamespaceFooter.H"

#endif  /* ! defined _ADJACENTBLOCK_H_ */
