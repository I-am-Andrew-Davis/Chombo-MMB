#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "MultiBlockLevelGeom.H"
#include "BoxFixedOff.H"

#include "NamespaceHeader.H"

/// constructor
MultiBlockLevelGeom::MultiBlockLevelGeom(const MultiBlockCoordSys*  a_coordSysPtr,
                                         const DisjointBoxLayout&   a_grids,
                                         int                        a_ghosts,
                                         Interval                   a_fixedDims,
                                         Vector<int>                a_fixedPt)
  :
  m_isDefined(false)
{
  define(a_coordSysPtr, a_grids, a_ghosts,
         a_fixedDims, a_fixedPt);
}

/// destructor
MultiBlockLevelGeom::~MultiBlockLevelGeom()
{
  undefine();
}

void MultiBlockLevelGeom::undefine()
{
  if (m_isDefined)
    {
      delete m_mbUtil;
      m_interpDimsVect.clear();
      m_fixedDimsVect.clear();
      DataIterator dit = m_grids.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          delete m_validInfo[dit];
        }
    }
  m_isDefined = false;
}

//  Define knowing only ghosts required for the solution update procedure. DEPRECATED.
/** \param[in]  a_ghosts
 *                      The number of ghosts that are required in the solution
 *                      update procedure can be interpolated.
 *  \param[in]  a_spaceOrder
 *                      Spatial order of the scheme
 */
void
MultiBlockLevelGeom::define(const MultiBlockCoordSys* a_coordSysPtr,
                            const DisjointBoxLayout&  a_grids,
                            const int                 a_solUpdateGhosts,
                            const int                 a_spaceOrder,
                            Interval                  a_fixedDims,
                            Vector<int>               a_fixedPt)
{
  // CH_assert(a_spaceOrder >= 2); // Why the assert for something unused?
  // petermc, 20 Nov 2012: earlier had a_solUpdateGhosts + a_spaceOrder/2
  // instead of a_solUpdateGhosts, but that was with a previous version
  // of MultiBlockLevelExchange.
  define(a_coordSysPtr, a_grids, a_solUpdateGhosts,
         a_fixedDims, a_fixedPt);
}

//  Define with explicit specification of number of multiblock ghosts
/** This function should be avoided in user code unless you want to explicitly
 *  maintain the required number of ghosts.
 *  \param[in]  a_ghosts
 *                      The number of ghosts that must exchanged across block
 *                      boundaries so that the number of ghosts required
 *                      in the solution advance procedure can be interpolated.
 */
void
MultiBlockLevelGeom::define(const MultiBlockCoordSys*  a_coordSysPtr,
                            const DisjointBoxLayout&   a_grids,
                            int                        a_ghosts,
                            Interval                   a_fixedDims,
                            Vector<int>                a_fixedPt)
{
  // pout() << "Begin MBLG define" << std::endl;
  CH_TIMERS("MultiBlockLevelGeom::define");
  CH_TIMER("getValid()", t1_getValid);
  CH_TIMER("wait for other process", t2_waitOther);
  undefine();
  CH_assert(a_coordSysPtr->isDefined());
  m_coordSysPtr = (MultiBlockCoordSys*) a_coordSysPtr;
  m_ghosts = a_ghosts;

  // MultiBlockUtil has a new in it.
  m_mbUtil = new MultiBlockUtil(m_coordSysPtr, a_fixedDims, a_fixedPt);

  m_fixedDims = a_fixedDims;
  m_interpDimsVect.clear();
  m_fixedDimsVect.clear();
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if ( a_fixedDims.contains(idir) )
        {
          m_fixedDimsVect.push_back(idir);
        }
      else
        {
          m_interpDimsVect.push_back(idir);
        }
    }
  m_fixedPt = Vector<int>(a_fixedPt);

  m_gridsFull = a_grids;

  // Set BoxLayout m_gridsFixedOff to DisjointBoxLayout m_gridsFull
  // modified so that every Box
  // has range 0 in m_interpDimsVect, and
  // shifted down by m_gridsFull.smallEnd(m_fixedDimsVect) in m_fixedDimsVect.
  // Purpose: all cells in m_ghostCells[dit] + m_gridsFixedOff[dit]
  // will have the same interpolation stencil.
  m_allGridsHaveFixedPt = m_mbUtil->allGridsHaveFixedPt(m_gridsFull);
  m_mbUtil->getCollapsedLayout(m_grids, m_gridsFull);
  m_mbUtil->getFixedOffLayout(m_gridsFixedOff, m_gridsFull);

  CH_assert(m_coordSysPtr->gotBoundaries());
  // access boundaries
  // m_coordSysPtr->boundaries();

  DataIterator dit = m_grids.dataIterator();

  // Get m_block.
  m_block.define(m_grids);
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& baseBox = m_grids[dit];
      // LayoutData<int> m_block;
      int thisBlock = m_coordSysPtr->whichBlock(baseBox);
      m_block[dit] = thisBlock;
    }

  // Get m_ghostCells.
  extraBlockGhosts(m_ghostCells, m_ghosts);
  m_validInfo.define(m_grids);

  CH_START(t1_getValid);
  for (dit.begin(); dit.ok(); ++dit)
    {
      const IntVectSet& ghostCellsIVS = m_ghostCells[dit];
      m_validInfo[dit] = new IVSFAB<MBLGeomValidInfo>(ghostCellsIVS, 1);
      m_mbUtil->getValid(*m_validInfo[dit], ghostCellsIVS, m_block[dit]);
    }
  CH_STOP(t1_getValid);
  CH_START(t2_waitOther);
#ifdef CH_MPI
  // This barrier is just for timings.  If any proc takes "forever" to get
  // valid, every processor timing will report it in
  // MultiBlockLevelGeom::define.
  MPI_Barrier(Chombo_MPI::comm);
#endif
  CH_STOP(t2_waitOther);
  m_isDefined = true;
  // pout() << "Finished MBLG define" << std::endl;
}


void
MultiBlockLevelGeom::extraBlockGhosts(LayoutData< IntVectSet >& a_ghostCells,
                                      int a_ghostLayer) const
{
  CH_TIME("MultiBlockLevelGeom::extraBlockGhosts");
  a_ghostCells.define(m_grids);
  DataIterator dit = m_grids.dataIterator();
  if (a_ghostLayer == 0)
    { // empty
      for (dit.begin(); dit.ok(); ++dit)
        {
          a_ghostCells[dit] = IntVectSet();
        }
    }
  else
    {
      for (dit.begin(); dit.ok(); ++dit)
        {
          const Box& baseBox = m_grids[dit];
          int baseBlockNum = m_block[dit];
          a_ghostCells[dit] =
            m_mbUtil->extraBlockGhosts(baseBox, a_ghostLayer, baseBlockNum);
        }
    }
}

#include "NamespaceFooter.H"
