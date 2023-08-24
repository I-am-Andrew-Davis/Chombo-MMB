#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "ExternalCS.H"

#include <algorithm>
#include <vector>

#include "BoxIterator.H"
#include "LoHiCenter.H"
#include "LoHiSide.H"
#include "RootSolver.H"
#include "ReadCGNS.H"

#include "NamespaceHeader.H"

#ifdef USE_SPLINE
/// Pull a mapping from an external source
/*--------------------------------------------------------------------*/
///  Constructor
/**
 *   \param[in]  a_dX   The computational cell size
 *   \param[in]  a_scale
 *                      Scaling to apply to the mesh
 *   \param[in]  a_domain
 *                      The problem domain
 *   \param[in]  a_file_name
 *                      The name of the external file to read
 *//*-----------------------------------------------------------------*/
ExternalCS::ExternalCS(const ReadCGNS& a_cgns,
                       const RealVect& a_dX,
                       const RealVect& a_scale,
                       const ProblemDomain& a_domain,
                       const std::string a_file_name,
                       const std::string& a_zoneName,
                       const int a_blkIdx,
                       const IntVect& a_evalGhost,
                       const int a_startProc,
                       const int a_numProc)
  :
  m_grid_file(a_file_name),
  m_zoneName(a_zoneName),
  m_idxZone(a_blkIdx + 1),
  m_domain(a_domain),
  m_scale(a_scale),
  m_baseDx(RealVect::Unit), // The coarse grid must always have unit size!!
  m_evalGhost(a_evalGhost)
{
  CH_assert(m_dx < m_baseDx);
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
    m_KDadaptor = std::make_unique<CHArrayKDadaptor> (m_gridNodes);
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

/*--------------------------------------------------------------------*/
/// Destructor
/**
 *//*-----------------------------------------------------------------*/
ExternalCS::~ExternalCS()
{
}

/*--------------------------------------------------------------------*/
/// Allow caching of repeatedly used realCoord values
/**
 *//*-----------------------------------------------------------------*/
void
ExternalCS::enableCache()
{
  m_useCache = true;
}

/*--------------------------------------------------------------------*/
/// Disallow caching of realCoord values, and clear existing cache
/**
 *//*-----------------------------------------------------------------*/
void
ExternalCS::disableCache()
{
  m_useCache = false;
  m_coordCache->clear();
}

/*--------------------------------------------------------------------*/
/// Given a coordinate in computational space, adjust it for periodic 
/// and scale as needed to fit into mapping
/**
 *   \param[in]  a_Xi   Location possibly outside of domain
 *   \return            Location shifted inside domain, if periodic
 *//*-----------------------------------------------------------------*/
IntVect
ExternalCS::coordAdj(RealVect& a_Xi) const
{
  IntVect shift(IntVect::Zero); // vector of domain shifts
  if(m_domain.isPeriodic())
    {
      a_Xi /= m_dx; // integer index
      Box bbox = m_domain.domainBox();
      bbox.convert(IntVect::Unit);
      // adjust for periodic, cycle point around
      for (int d=0; d!=SpaceDim; d++)
        {
          if(m_domain.isPeriodic(d))
            {
              if(a_Xi[d] > (Real)bbox.bigEnd(d))
                {
                  a_Xi[d] -= bbox.size(d)-1;
                  shift[d] = 1;
                }
              if(a_Xi[d] < (Real)bbox.smallEnd(d))
                {
                  a_Xi[d] += bbox.size(d)-1;
                  shift[d] = -1;
                }
            }
        }
      a_Xi *= m_dx; // switch back to normal index
    }
  //a_Xi *= m_scaleXi; // scale as needed
  a_Xi = m_xiOrigin + (a_Xi-m_xiOrigin)*m_scaleXi;
  return shift;
}

/*--------------------------------------------------------------------*/
/// Given coordinate in mapped space, return its location in real space
/// This does gives the raw transformation without accounting for
/// periodic boundaries.
/**
 *   \param[in]  a_Xi   Location in computational space
 *   \return            Location in physical space
 *//*-----------------------------------------------------------------*/
RealVect
ExternalCS::realCoordNonAdjst(const RealVect& a_Xi) const
{
  RealVect realLoc;
  // evaluate the metrics
  for (int d=0; d!=SpaceDim; d++)
    {
      realLoc[d] = m_gridInterp.interp(a_Xi, d);
    }
  return realLoc*m_scale;
}

/*--------------------------------------------------------------------*/
/// Given coordinate in mapped space, return its location in real space
/**
 *   \param[in]  a_Xi   Location in computational space
 *   \return            Location in physical space
 *//*-----------------------------------------------------------------*/
RealVect
ExternalCS::realCoordRaw(const RealVect& a_Xi) const
{
  CH_TIME("ExernalCS::realCoordRaw");
  // shift the coordinate
  RealVect Xi(a_Xi);
  IntVect shift = coordAdj(Xi);
  
  // solve for periodic offset
  RealVect offset(RealVect::Zero);
  if(m_domain.isPeriodic() && (shift != IntVect::Zero))
    {
      offset = D_TERM(  (shift[0]*m_domainOffset[0]),
                      + (shift[1]*m_domainOffset[1]),
                      + (shift[2]*m_domainOffset[2]));
    }
  return realCoordNonAdjst(Xi) + offset;
}

/*--------------------------------------------------------------------*/
/// Given coordinate in mapped space, return its location in real space
/**
 *   \param[in]  a_Xi   Location in computational space
 *   \return            Location in physical space
 *//*-----------------------------------------------------------------*/
RealVect
ExternalCS::realCoord(const RealVect& a_Xi) const
{
  CH_TIME("ExernalCS::realCoord");
  if (m_useCache)
    {
      // attempt conversion to the near cell/face point
      RealVect key(2*a_Xi/m_dx);
      if (std::all_of(key.dataPtr(), key.dataPtr() + (SpaceDim-1),
                      [&](Real x){
                        Real intPart;
                        Real relEps = 2*std::numeric_limits<Real>::epsilon()*std::abs(x);
                        return std::abs(std::modf(x, &intPart)) <= relEps;}))
        {
          return (*m_coordCache)(key);
        }
    }
  // else
  return realCoordRaw(a_Xi);
}

/*--------------------------------------------------------------------*/
/// A global brute force search for the coarse inverse in a region
/// This is guaranteed to find the closes grid point, but is slow
/**
 *   \return            Approximate location in computational space
 *   \param[in]  a_x    Location in physical space
 *   \param[in]  a_dxi  grid spacing
 *   \param[in]  a_search
 *                      Grid box to search in
 *//*-----------------------------------------------------------------*/
RealVect
ExternalCS::globalInverseSearch(const RealVect& a_x,
                                const RealVect& a_dxi,
                                const Box& a_search) const
{
  CH_TIME("ExernalCS::globalInverseSearch");
  // mapping function, relative to solution
  auto func = [&](RealVect xi){return realCoord(xi) - a_x;};
  // initial guess, middle of the domain?
  RealVect xi(m_domain.size()/2);
  RealVect xi_new(xi);
  auto funcErr = [&](RealVect xi){auto fxi = func(xi);
                                  return dot(fxi, fxi);};
  Real err = funcErr(xi);
  Real err_new(err);
  MD_BOXLOOP(a_search, i)
    {
      // xi_new = a_dxi*(MD_GETIV(i) + 0.5*RealVect::Unit);
      xi_new = centerMappedCoordinates(MD_GETIV(i));
      err_new = funcErr(xi_new);
      if (err_new < err)
        {
          err = err_new;
          xi = xi_new;
        }
    }
  return xi;
}

/*--------------------------------------------------------------------*/
/// A global search for the coarse inverse given an initial guess
/// An efficient method for finding the closes grid point
/**
 *   \return            Approximate location in computational space
 *   \param[in]  a_x    Location in physical space
 *   \param[in]  a_dxi  grid spacing
 *   \param[in]  a_search
 *                      Grid box to search in.
                        One cell past this is evaluated
 *   \param[in]  a_xiEst
 *                      Initial guess
 *//*-----------------------------------------------------------------*/
RealVect
ExternalCS::inverseSearch(const RealVect& a_x,
                          const RealVect& a_dxi,
                          const Box& a_search,
                          const RealVect& a_xiEst) const
{
  CH_TIME("ExernalCS::inverseSearch");
  // mapping function at cell center, relative to solution
  // euclidean distance squared
  auto funcErr = [&](IntVect cellIdx)
    {
      //auto xi = centerMappedCoordinates(cellIdx);
      RealVect xi = cellIdx*a_dxi;
      RealVect fxi = realCoord(xi) - a_x;
      return dot(fxi, fxi);
    };
  // force a point to be inside a box
  auto projectIntoBox = [] (const Box& a_box, const IntVect& a_iv)
    {
      auto newIv = a_iv;
      newIv = stc::max(newIv, a_box.smallEnd());
      newIv = stc::min(newIv, a_box.bigEnd());
      CH_assert(a_box.contains(newIv));
      return newIv;
    };
  // initial cell index for xi
  IntVect xiEst = a_xiEst/a_dxi;
  xiEst = projectIntoBox(a_search, xiEst);
  // information for regions to search
  Box boundRadius(-IntVect::Unit, IntVect::Unit);
  Box interestBox(xiEst, xiEst);
  Box evalBox(grow(interestBox,1));
  Box discardBox(interestBox);
  // Store the min on each side of the search region, with access function
  std::array<Real, SpaceDim*2> faceMinVal;
  std::array<IntVect, SpaceDim*2> faceMinIdx;
  std::array<bool, SpaceDim*2> faceInBound;
  auto indexFace = [&](const int a_dir, const Side::LoHiSide a_side)
    {
      int shift = (a_side == Side::Lo) ? 0 : 1;
      return a_dir + SpaceDim*shift;
    };
  auto getIndexFaceMin = [&](int& a_dir,
                             Side::LoHiSide& a_side,
                             bool a_checkBound = true)
    {
      int minI = -1;
      for (int i = 0; i!=faceMinVal.size(); ++i)
        {
          if (!a_checkBound || faceInBound[i])
            {
              if (minI < 0)
                {
                  minI = i;
                }
              else
                {
                  minI = (faceMinVal[i] < faceMinVal[minI]) ? i : minI;
                }
            }
        }
      CH_assert(faceInBound[minI] || !a_checkBound);
      CH_assert(minI >= 0);
      CH_assert(minI < faceMinVal.size());
      a_dir = minI%SpaceDim;
      a_side = (minI < SpaceDim) ? Side::Lo : Side::Hi;
      return indexFace(a_dir, a_side);
    };
  // setup the initial region
  {
    FArrayBox evalRegion(evalBox, 1);
    MD_ARRAY_RESTRICT(evalArr, evalRegion);
    MD_BOXLOOP(evalBox, i)
      {
        evalArr[MD_IX(i, 0)] = funcErr(MD_GETIV(i));
      }
    auto boundBox = boundRadius + xiEst;
    Real ivDist = evalRegion(xiEst, 0);
    bool ivContained = true;
    MD_BOXLOOP(boundBox, j)
      {
        ivContained &= (evalArr[MD_IX(j, 0)] >= ivDist);
      }
    if (ivContained)
      {
        // solution found, search done
        goto inverseFound;
      }
    
    // set the face min values
    for (auto dir : EachDir)
      {
        for (auto side : EachSide)
          {
            auto sideIdx = indexFace(dir, side);
            IntVect sideIv = xiEst + sign(side)*BASISV(dir);
            faceMinIdx[sideIdx] = sideIv;
            faceMinVal[sideIdx] = evalRegion(sideIv, 0);
            faceInBound[sideIdx] = a_search.contains(sideIv);
          }
      }
  }
  
  // Search until either the inverse is found, or the region is exhausted
  // the check is extra, but helps save infinite loop on error
  while ((a_search != discardBox) &&
         (a_search.contains(discardBox)))
    {
      // find the new search direction
      Side::LoHiSide searchSide;
      int searchDir;
      getIndexFaceMin(searchDir, searchSide);
      // setup the region of interest
      interestBox = adjCellBox(discardBox, searchDir, searchSide, 1);
      CH_assert(a_search.contains(interestBox));
      // eval the region
      evalBox = grow(interestBox, 1);
      CH_assert(evalBox.contains(evalBox));
      FArrayBox evalRegion(evalBox, 1);
      MD_ARRAY_RESTRICT(evalArr, evalRegion);
      MD_BOXLOOP(evalBox, i)
        {
          evalArr[MD_IX(i, 0)] = funcErr(MD_GETIV(i));
        }
      // check if any cell in the region of interest is bounded
      // and track any closer surrounding cells
      CH_assert(evalRegion.contains(interestBox));
      MD_BOXLOOP(interestBox, i)
        {
          IntVect iv(MD_GETIV(i));
          auto boundBox = boundRadius + iv;
          CH_assert(evalRegion.contains(boundBox));
          Real ivDist = evalRegion(iv, 0);
          bool ivContained = true;
          MD_BOXLOOP(boundBox, j)
            {
              ivContained &= (evalArr[MD_IX(j, 0)] >= ivDist);
            }
          // a bounded region found, exit
          if (ivContained)
            {
              xiEst = iv;
              goto inverseFound;
            }
        }
      // update the tracked minimum from the new eval region
      for (auto dir : EachDir)
        {
          for (auto side : EachSide)
            {
              if (!(dir == searchDir && side == flip(searchSide)))
                {
                  bool requireUpdate(dir == searchDir && side == searchSide);
                  auto sideIdx = indexFace(dir, side);
                  auto oldMin = faceMinVal[sideIdx];
                  auto sideBox = adjCellBox(interestBox, dir, side, 1);
                  CH_assert(!discardBox.contains(sideBox));
                  CH_assert(evalRegion.contains(sideBox));
                  MD_BOXLOOP(sideBox, k)
                    {
                      IntVect ivk(MD_GETIV(k));
                      auto newMin = evalArr[MD_IX(k, 0)];
                      // force update in search directions
                      if ((newMin < oldMin) || requireUpdate)
                        {
                          faceMinVal[sideIdx] = newMin;
                          faceMinIdx[sideIdx] = ivk;
                          faceInBound[sideIdx] = a_search.contains(ivk);
                          requireUpdate = false;
                        }
                    }
                  CH_assert(!discardBox.contains(faceMinIdx[sideIdx]));
                }
            }
        }
      // expand the discarded region
      discardBox.growDir(searchDir, searchSide, 1);
    }
  // no inverse found in the search region, so return the closest point
  CH_assert(discardBox == a_search);
  Side::LoHiSide minSide;
  int minDir;
  xiEst = faceMinIdx[getIndexFaceMin(minDir, minSide, false)];
  xiEst = projectIntoBox(a_search, xiEst);
  
  inverseFound:
  CH_assert(a_search.contains(xiEst));
  return (RealVect(xiEst)+0.5)*a_dxi;
}

/*--------------------------------------------------------------------*/
/// A global search for the coarse inverse
/// An efficient method for finding the closes grid point
/// This uses a KD-tree to find the closest node value in O(log n) time
/**
 *   \return            Approximate location in computational space
 *   \param[in]  a_x    Location in physical space
 *//*-----------------------------------------------------------------*/
RealVect
ExternalCS::inverseSearch(const RealVect& a_x) const
{
  CH_TIME("KD-tree lookup");
  size_t ret_index;
  Real out_dist_sqr;
  int num_results = 1;
  nanoflann::KNNResultSet<Real> resultSet(num_results);
  resultSet.init(&ret_index, &out_dist_sqr);
  m_KDtree->findNeighbors(resultSet, a_x.dataPtr(), nanoflann::SearchParams());
  // Get the index of the nearest node
  Box nodeBox(IntVect::Zero, m_gridSize-IntVect::Unit);
  BoxIterator bit(nodeBox);
  IntVect KDidx(bit.at(ret_index));
  // convert to index on coarse domain
  Box dom = m_domain.domainBox();
  dom.coarsen(m_refRatio);
  KDidx += dom.smallEnd();
  // convert to a xi value
  return RealVect(KDidx)*m_baseDx;
}


/*--------------------------------------------------------------------*/
/// Check if an inverse is inside the domain
/**
 *   \return             If a_x is inside the domain
 *   \param[in]  a_x     Location in physical space
 *   \param[out] a_xiEst A point near the inverse, and always internal
 *                       This should always be an adequate initial guess
 *                       for iterative solvers if an inverse exists
 *//*-----------------------------------------------------------------*/
bool
ExternalCS::inBlockDomain(const RealVect &a_x,
                          RealVect & a_xiEst) const
{
  CH_TIME("ExernalCS::inBlockDomain");
  // like Box.contains(iv) except for with RealVects and node bases boxes
  auto containsPoint = [&](const RealVect& a_xi,
                           const Box& a_box)
  {
    RealVect point = (a_xi/m_baseDx);
    return ((RealVect(a_box.bigEnd()) >= point) &&
            (RealVect(a_box.smallEnd()) <= point));
  };

  // Get the nearest node point on the coarse grid from the KD-tree
  a_xiEst = inverseSearch(a_x);
  // Define the coarse domain
  Box searchBox = surroundingNodes(m_domain.domainBox());
  searchBox.coarsen(m_refRatio);
  // Some points are clearly bounded in the domain
  Box inside = grow(searchBox, -1);
  if (containsPoint(a_xiEst, inside))
    {
      return true;
    }
  // Points near the boundary are hard to detect
  else
    {
      // Check if the point is feasibly in the domain
      // If the coarse inverse is not interior or adjacent to the domain
      // boundary, it can not be possibly be contained in the domain.

      // find the directions of the boundary from nearIdx
      std::vector<RealVect> boundaryDirs;
      for (const auto dir : EachDir)
        {
          for (const auto side : EachSide)
            {
              const RealVect testDir = RealVect_basis(dir)*sign(side)*m_baseDx;
              if (!containsPoint(a_xiEst + testDir, searchBox))
                {
                  boundaryDirs.push_back(testDir);
                }
            }
        }
      // nearIdx must be on at least one block boundary
      CH_assert(boundaryDirs.size() > 0);
      // find the bounded face of nearIdx which may contain an inverse
      // track the boundary adjacent point nearest to the inverse
      RealVect adj = a_xiEst;
      Real minAdjDist = stc::mag(realCoord(adj) - a_x);
      Real minOut2Dist = std::numeric_limits<Real>::max();
      for (const auto bdDir : boundaryDirs)
        {
          // update at 2 cells out
          RealVect testAdj2 = adj + 2*bdDir;
          minOut2Dist = std::min(minOut2Dist,
                                 stc::mag(realCoord(testAdj2) - a_x));
          // update at 1 cell out
          RealVect testAdj = adj + bdDir;
          Real testAdjDist = stc::mag(realCoord(testAdj) - a_x);
          if (testAdjDist < minAdjDist)
            {
              adj = testAdj;
              minAdjDist = testAdjDist;
            }
        }

      // The nearest adjacent point must be closer to a_x than any other external
      // point for the point to feasibly be inside
      if (minAdjDist < minOut2Dist)
        {
          // initial guess is confirmed to be within 2 coarse cells of the inverse
          // which allows an iterative solve to get the exact inverse
          RealVect xiExact = ExternalCS::mappedCoord(a_x, a_xiEst);
          // see if the solution is insdie the block
          bool bounded = containsPoint(xiExact, searchBox);
          // Save the exact inverse if it's the solution we want and exit
          if (bounded)
            {
              a_xiEst = xiExact;
              return true;
            }
        }
      return false;
    }
}

/*--------------------------------------------------------------------*/
/// Get Xi of the cell center point nearest to X
/**
 *   \param[in]  a_x    Location in physical space
 *   \return            Cell center in mapped space
 *//*-----------------------------------------------------------------*/
RealVect
ExternalCS::nearestDomainCell(const RealVect& a_x) const
{
  // The inverse falls outside the domain if this point is reached, so
  // update a_xiEst to the nearest cell center.
  CH_TIME("ExernalCS::nearestDomainCell");
  // A function for searching the nearest cells
  auto nearestCell = [&](RealVect& a_xi,
                         const RealVect& a_x,
                         const Box& a_box,
                         const RealVect& a_dx)
                       {
                         Real nearestDist = std::numeric_limits<Real>::max();
                         IntVect nearCellIdx;
                         MD_BOXLOOP(a_box, i)
                           {
                             IntVect cellIdx = MD_GETIV(i);
                             RealVect xiCell = (RealVect(cellIdx)+0.5)*a_dx;
                             Real dist = stc::mag(realCoord(xiCell) - a_x);
                             if (dist < nearestDist)
                               {
                                 nearCellIdx = cellIdx;
                                 nearestDist = dist;
                                 a_xi = xiCell;
                               }
                           }
                         return nearCellIdx;
                       };

  // Get the nearest node point on the coarse grid from the KD-tree
  RealVect xiCell = inverseSearch(a_x);
  // Define the coarse domain
  Box searchBox = m_domain.domainBox();
  searchBox.coarsen(m_refRatio);
  // And the index of the nearest node
  IntVect nodeIdx = xiCell;
  stc::forEachElement<SpaceDim>(
      [&] (const stc::array_size_type a_idx){
        nodeIdx[a_idx] = std::round(xiCell[a_idx]/m_baseDx[a_idx]);
      });

  // Grow a box internal to the domain around the nearest node
  Box nearBox(nodeIdx-IntVect::Unit, nodeIdx);
  nearBox &= searchBox;
  // Check all (4 at most in 3D) cell centers for the min distance
  // updates xiCell
  IntVect cCellIdx = nearestCell(xiCell, a_x, nearBox, m_baseDx);
  // Refine the nearest cell, and check all fine cells that are along
  // the boundaries
  if (m_refRatio > 1)
    {
      // Refine the cell
      Box fnearBox(cCellIdx, cCellIdx);
      fnearBox.refine(m_refRatio);
      // Check all cell centers, update xiCell
      nearestCell(xiCell, a_x, fnearBox, m_baseDx/m_refRatio);
    }
  return xiCell;
}

/*--------------------------------------------------------------------*/
/// Given coordinate in real space, return its location in mapped space
/**
 *   \param[in]  a_x    Location in physical space
 *   \return            Location in mapped space
 *//*-----------------------------------------------------------------*/
RealVect
ExternalCS::mappedCoord(const RealVect& a_x) const
{
  CH_TIME("ExernalCS::mappedCoord");
  RealVect mapLoc;
  // Get a initial guess by doing a brute force search, ghost cells included
  //  this should be within a dxi of the inverse
  auto searchBox = m_domain.domainBox();
  searchBox.grow(m_evalGhost);
  auto Xi_est = globalInverseSearch(a_x, m_dx, searchBox);
  // Use the newton solver to get an exact inverse
  mapLoc = ExternalCS::mappedCoord(a_x, Xi_est);
  return mapLoc;
}

/*--------------------------------------------------------------------*/
/// Given coordinate in real space, return its location in mapped space
/**
 *   \param[in]  a_x    Location in physical space
 *   \param[in]  a_xiEst Initial guess in mapped space
 *   \return            Location in mapped space
 *//*-----------------------------------------------------------------*/
RealVect
ExternalCS::mappedCoord(const RealVect& a_x,
                        const RealVect& a_xiEst) const
{
  RealVect mapLoc;
  int numIter;
  // mapping function and derivative
  // auto func = [&](RealVect xi, int comp){
  //   return m_gridInterp.interp(xi, comp) - a_x[comp];};
  // auto funcD = [&](RealVect xi, int comp, RealVect d){
  //   return m_gridInterp.interpD(xi, comp, d);};
  auto func = [&](RealVect xi, int comp){
    auto x = realCoordRaw(xi);
    return x[comp] - a_x[comp];};
  auto funcD = [&](RealVect xi, int comp, RealVect d){
    int deriv;
    for (auto dir : EachDir)
      {
        if (d[dir] == 1)
          {
            deriv = dir;
          }
      }
    return dXdXi(xi, comp, deriv);};
  // Use the newton solver to get an exact inverse
  Real invTolerance = 1.0e-10;
  mapLoc = RootSolver::Newton(numIter, func, funcD, a_xiEst, invTolerance);
  return mapLoc;
}

/*--------------------------------------------------------------------*/
/// Calculate the derivative of each coordinate vector
/**
 *   \param[in]  a_Xi   Location in computational space
 *   \param[in]  a_dirX The component in physical space (x, y)
 *   \param[in]  a_dirXi
 *                      The component in computational space (xi, eta)  
 *   \return            derivative of dX/dXi
 *//*-----------------------------------------------------------------*/
Real
ExternalCS::dXdXi(const RealVect& a_Xi, int a_dirX, int a_dirXi) const
{
  CH_TIME("ExernalCS::dXdXi");
  // set something invalid, which should get overwritten - returns NaN
  Real value = 1.0/0.0;;
  // make the periodic adjustment
  RealVect Xi(a_Xi);
  coordAdj(Xi);
  
  IntVect deriv(IntVect::Zero);
  deriv[a_dirXi] = 1;
  value = m_gridInterp.interpD(Xi, a_dirX, deriv)
    *m_scale[a_dirXi]*m_scaleXi[a_dirXi];
  return value;
}

/*--------------------------------------------------------------------*/
/// Read in and store grid values. Compute all metrics needed
/**
 *//*-----------------------------------------------------------------*/

void ExternalCS::readGridFile(const ReadCGNS& a_cgns)
{
  CH_TIME("ExternalCS::readGridFile");
#if !(defined(CH_USE_CGNS) || defined(CH_USE_OGEN))
  MayDay::Error("Do not have ability to read in a grid!\nTry including an additional libary, such as CGNS or Ogen"); 
#endif
  // Read in a CGNS file
#ifdef CH_USE_CGNS 
  a_cgns.readZoneCoords(m_idxZone, m_gridNodes, m_gridSize);
  m_scaleXi = RealVect(m_gridSize - 1)/m_domain.size();
#endif

  // Solve for metrics
  m_gridInterp.define(m_gridNodes.begin(), // interp data
                      m_gridSize,      // data size
                      SpaceDim,        // data dimension
                      ColStorage,      // data storage organization
                      m_dx,            // step spacing
                      m_xiOrigin,       // origin
                      m_startProc,
                      m_numProc);
  // check interpolation over the correct domain
  // pout() << "Origin" << m_xiOrigin <<std::endl;
  // pout() << "Comparing " << *m_gridA.begin();
  // pout() << " to "<< realCoord(m_dx*m_domain.domainBox().smallEnd())[0] << std::endl;
  // pout() << "Comparing " << *m_gridA.end();
  // pout() << " to "<< realCoord(m_dx*m_domain.domainBox().bigEnd())[0] << std::endl;
  // setCaching(true); // cache spline values
}
#endif // end #ifdef USE_SPLINE


// -- begin factory implementations ---------------------------

/*--------------------------------------------------------------------*/
///  Constructor
/**
 *   \param[in]  a_file_name
 *                      The name of the external file to read
 *   \param[in]  a_scale
 *                      Scaling to apply to the mesh
 *//*-----------------------------------------------------------------*/
ExternalCSFactory::ExternalCSFactory(const std::string& a_grid_file,
                                     const RealVect& a_scale):
  m_grid_file(a_grid_file),
  m_scale(a_scale)
{
}

/*--------------------------------------------------------------------*/
///  Return the ExternalCS object
/**
 *   \param[in]  a_domain
 *                      Problem domain the coordinate system is in
 *   \return            Pointer to the coordinate system
 *//*-----------------------------------------------------------------*/
NewCoordSys*
ExternalCSFactory::getCoordSys(const ProblemDomain& a_domain,
                               const RealVect& a_dx) const
{
  CH_assert(false);
  IntVect domSize = a_domain.size();
  RealVect domL;
  for (auto dir : EachDir)
    {
      domL[dir] = domSize[dir]*a_dx[dir];
    }
  // !! FIXME !!
  // ExternalCS* newCSPtr = new ExternalCS(a_dx, m_scale, a_domain, m_grid_file, 0);
  // return static_cast< NewCoordSys* >(newCSPtr);
  return nullptr;
}

#include "NamespaceFooter.H"
