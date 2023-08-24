#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "ExtrudedExternalCS.H"

#include "NamespaceHeader.H"

#ifdef USE_SPLINE

/// Pull a mapping from an external source
/*--------------------------------------------------------------------*/
///  Constructor
/**
 *   \param[in]  a_dX   The computational cell size
 *   \param[in]  a_scale
 *                      Scaling to apply to the mesh; should always be 1
 *   \param[in]  a_domain
 *                      The problem domain
 *   \param[in]  a_file_name
 *                      The name of the external file to read
 *   \param[in]  a_zoneName
 *                      The name of the 1 unique zone to be read
 *   \param[in]  a_constIdx
 *                      Idx that insures only the 1 correct zone will be read
 *   \param[in]  a_blkIdx
 *                      The idx of the block that is being created
 *   \param[in]  a_extDx
 *                      The physical spacing of the block
 *   \param[in]  a_extCell
 *                      The number of cells that will be created in each dir
 *//*-----------------------------------------------------------------*/
ExtrudedExternalCS::ExtrudedExternalCS(const ReadCGNS& a_cgns,
                                       const RealVect& a_dX,
                                       const RealVect& a_scale,
                                       const ProblemDomain& a_domain,
                                       const std::string a_file_name,
                                       const std::string& a_zoneName,
                                       const int a_constIdx, //this is so the correct zone in the .cgns file is read everytime
                                       const IntVect& a_evalGhost,
                                       const int a_startProc,
                                       const int a_numProc,
                                       const int a_blkIdx,
                                       const RealVect a_extDx,
                                       const IntVect a_extCell
)

:
  ExternalCS(a_cgns,
             a_dX,
             a_scale,
             a_domain,
             a_file_name,
             a_zoneName,
             a_constIdx, //this is so the correct zone in the .cgns file is read everytime
             a_evalGhost,
             a_startProc,
             a_numProc),

  m_blkIdx(a_blkIdx),
  m_extDx(a_extDx),
  m_extCell(a_extCell)
{
 CH_assert(m_dx <= m_baseDx);
 m_dx = a_dX;
 m_refRatio = m_baseDx[0]/m_dx[0]; // this should work out exactly
 CH_assert(m_refRatio >= 1);
 CH_assert(trunc(m_baseDx[0]/m_dx[0]) == m_refRatio);
 m_xiOrigin = RealVect{m_dx*m_domain.domainBox().smallEnd()};
 // parallel construction
 m_startProc = a_startProc;
 m_numProc = a_numProc;
 CH_assert((m_startProc >= 0) && (m_startProc < numProc()));
 CH_assert((m_numProc >= 0) && (m_numProc <= numProc()));


 // setup caching for the inverse
 {
  std::function<RealVect(const IntVect&)> coordFunc =
  [&]
  (const IntVect& a_key) -> RealVect
  {
   RealVect key(a_key);
   key *= m_dx/2;
   return realCoordRaw(key);
  };
  m_coordCache = std::make_unique<coordCache> (std::move(coordFunc));
 }
 disableCache(); // default empty

 // Setting up the mapping, this has some expense

 readGridFile(a_cgns);

 // Create the KD-tree adapter, using the read in point data
 {
  CH_TIME("KD-tree construction");
  m_KDadaptor = std::make_unique<CHArrayKDadaptor> (m_extrudedGridNodes);
  // Create the KD-tree, using leaf size of 50
  m_KDtree = std::make_unique<KDtree> (SpaceDim,
                                       *m_KDadaptor,
                                       nanoflann::KDTreeSingleIndexAdaptorParams(50));
  m_KDtree->buildIndex();
 }

 /// Solve domain offset
 m_domainOffset.resize(SpaceDim);
 // get face box
 Box bbox = m_domain.domainBox();
 bbox.convert(IntVect::Unit);
 for (int d = 0; d != SpaceDim; ++d)
   {
     // get a Xi domain length vector

     RealVect shift(RealVect::Zero);
     shift[d] = 1.0;
     RealVect xi_lo = shift*m_scaleXi*m_dx*bbox.smallEnd();
     RealVect xi_hi = shift*m_scaleXi*m_dx*bbox.bigEnd();
     // map to X space
     RealVect dom_lo = realCoordNonAdjst(xi_lo);
     RealVect dom_hi = realCoordNonAdjst(xi_hi);
     m_domainOffset[d] = dom_hi - dom_lo;
   }
}

ExtrudedExternalCS::~ExtrudedExternalCS()
{
}

//Read grid file function
/*--------------------------------------------------------------------*/
///  Reads in a cgns file and then creates a new block based on the block
///  read in
/**
 *   \param[in]  a_cgns
 *                      The CGNS file to be read in
 *   \param[out]  m_extrudedGridNodes
 *                      The CHArray of grid nodes making up the new block
 *//*-----------------------------------------------------------------*/
void ExtrudedExternalCS::readGridFile(const ReadCGNS& a_cgns)
{
 CH_TIME("ExternalCS::readGridFile");
#if !(defined(CH_USE_CGNS) || defined(CH_USE_OGEN))
 MayDay::Error("Do not have ability to read in a grid!\nTry including an additional libary, such as CGNS or Ogen");
#endif
 // Read in a CGNS file
#ifdef CH_USE_CGNS
 a_cgns.readZoneCoords(m_idxZone, m_gridNodes, m_gridSize);
 m_scaleXi = RealVect::Unit*m_dx;
#endif

 CH_assert(m_blkIdx >= 0 && m_blkIdx < 3);
//this is the "left" or first block
 if (m_blkIdx == 0){
  Box sizeBox(IntVect(D_DECL(0,0,0)),
              IntVect(D_DECL(m_extCell[0]-1,m_gridSize[1]-1, m_extCell[2]-1)));
  Box faceBox(IntVect(D_DECL(0,0,0)), IntVect(D_DECL(0,m_gridSize[1]-1,0)));
  m_extrudedGridNodes.define(sizeBox.size(),SpaceDim);
#if CH_SPACEDIM == 3
  for (int k=0; k<m_extCell[2]; ++k)
#endif
  {
   for (int j=0; j<m_extCell[0]; ++j)
   {
    MD_BOXLOOP(faceBox, i)
    {
     m_extrudedGridNodes(MD_GETIV(i)+IntVect(D_DECL(j,0,k)), 1) =
     m_gridNodes(MD_GETIV(i), 1);

     m_extrudedGridNodes(MD_GETIV(i)+IntVect(D_DECL(j,0,k)), 0) =
     (m_gridNodes(MD_GETIV(i), 0) - m_extDx[0]*(m_extCell[0] -1 -j));

#if CH_SPACEDIM == 3
     m_extrudedGridNodes(MD_GETIV(i)+IntVect(D_DECL(j,0,k)), 2) =
     (m_gridNodes(MD_GETIV(i), 2) + m_extDx[2]*(k));
#endif
    }
   }
  }
  // Solve for metrics
  m_gridInterp.define(m_extrudedGridNodes.begin(), // interp data
                      sizeBox.size(),      // data size
                      SpaceDim,        // data dimension
                      ColStorage,      // data storage organization
                      m_dx,            // step spacing
                      m_xiOrigin,       // origin
                      m_startProc,
                      m_numProc);
 }
//this is the "middle" or second block
 if(m_blkIdx == 1){

  Box sizeBox(IntVect(D_DECL(0,0,0)),
              IntVect(D_DECL(m_gridSize[0]-1,m_gridSize[1]-1, m_extCell[2]-1)));

  Box faceBox(IntVect(D_DECL(0,0,0)),
              IntVect(D_DECL(m_gridSize[0]-1,m_gridSize[1]-1,0)));
  m_extrudedGridNodes.define(sizeBox.size(),SpaceDim);

#if CH_SPACEDIM == 3
  for (int k=0; k<m_extCell[2]; ++k)
#endif
  {
   MD_BOXLOOP(faceBox, i)
   {
    IntVect gridPoint = IntVect(MD_GETIV(i));
    IntVect iter = MD_GETIV(i);

#if CH_SPACEDIM == 3
    iter[2] = k;
#endif
    m_extrudedGridNodes(iter, 0)=m_gridNodes(MD_GETIV(i), 0);
    m_extrudedGridNodes(iter, 1)=m_gridNodes(MD_GETIV(i), 1);

#if CH_SPACEDIM == 3
    m_extrudedGridNodes(iter, 2)=m_gridNodes(MD_GETIV(i), 2)+m_extDx[2]*(k);
#endif
   }
  }
  // Solve for metrics
  m_gridInterp.define(m_extrudedGridNodes.begin(), // interp data
                      sizeBox.size(),      // data size
                      SpaceDim,        // data dimension
                      ColStorage,      // data storage organization
                      m_dx,            // step spacing
                      m_xiOrigin,       // origin
                      m_startProc,
                      m_numProc);
 }

//this is the "right" or third block
 if (m_blkIdx == 2)
 {
  Box sizeBox(IntVect(D_DECL(0,0,0)),
              IntVect(D_DECL(m_extCell[0]-1,m_gridSize[1]-1, m_extCell[2]-1)));

  Box faceBox(IntVect(D_DECL(m_gridSize[0]-1,0,0)),
              IntVect(D_DECL(m_gridSize[0]-1,m_gridSize[1]-1,0)));
  m_extrudedGridNodes.define(sizeBox.size(),SpaceDim);

#if CH_SPACEDIM == 3
  for (int k=0; k<m_extCell[2]; ++k)
#endif
  {
   for (int j=0; j<m_extCell[0]; ++j)
   {
    MD_BOXLOOP(faceBox, i)
    {
     m_extrudedGridNodes(MD_GETIV(i)-IntVect(D_DECL(m_gridSize[0]-1,0,0))
                         +IntVect(D_DECL(j,0,k)), 0)=
     m_extDx[0]*(j) + m_gridNodes(MD_GETIV(i), 0);

     m_extrudedGridNodes(MD_GETIV(i)-IntVect(D_DECL(m_gridSize[0]-1,0,0))
                         +IntVect(D_DECL(j,0,k)), 1)=
     m_gridNodes(MD_GETIV(i), 1);

#if CH_SPACEDIM == 3
     m_extrudedGridNodes(MD_GETIV(i)-IntVect(D_DECL(m_gridSize[0]-1,0,0))
                         +IntVect(D_DECL(j,0,k)), 2)=
     m_gridNodes(MD_GETIV(i), 2) + m_extDx[2]*(k);
#endif
    }
   }
  }

  // Solve for metrics
  m_gridInterp.define(m_extrudedGridNodes.begin(), // interp data
                      sizeBox.size(),      // data size
                      SpaceDim,        // data dimension
                      ColStorage,      // data storage organization
                      m_dx,            // step spacing
                      m_xiOrigin,       // origin
                      m_startProc,
                      m_numProc);
 }

}
#endif // end #ifdef USE_SPLINE

#include "NamespaceFooter.H"
