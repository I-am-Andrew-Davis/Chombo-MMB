#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "OgenCS.H"

#include <algorithm>

#include "BoxIterator.H"
#include "LoHiCenter.H"
#include "LoHiSide.H"
#include "RootSolver.H"

#include "NamespaceHeader.H"

#ifdef CH_USE_OGEN
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
OgenCS::OgenCS(const RealVect& a_dX,
               const RealVect& a_scale,
               const ProblemDomain& a_domain,
               const std::string a_file_name,
               const int a_blkIdx,
               const IntVect& a_evalGhost)
  :
  m_grid_file(a_file_name),
  m_idxZone(a_blkIdx + 1),
  m_domain(a_domain),
  m_scale(a_scale),
  m_evalGhost(a_evalGhost)
{
  m_dx = a_dX;
  m_xiOrigin = RealVect{m_dx*m_domain.domainBox().smallEnd()};
  
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
  disableCache();
  
  /// Setting up the mapping, this has some expense
  readGridFile();
  
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
OgenCS::~OgenCS()
{
}

/*--------------------------------------------------------------------*/
/// Allow caching of repeatedly used realCoord values
/**
 *//*-----------------------------------------------------------------*/
void
OgenCS::enableCache()
{
  m_useCache = true;
}

/*--------------------------------------------------------------------*/
/// Disallow caching of realCoord values, and clear existing cache
/**
 *//*-----------------------------------------------------------------*/
void
OgenCS::disableCache()
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
OgenCS::coordAdj(RealVect& a_Xi) const
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
OgenCS::realCoordNonAdjst(const RealVect& a_Xi) const
{
  RealVect realLoc;
  // evaluate the metrics
  realLoc = ogen.realCoord(a_Xi);
  return realLoc*m_scale;
}

/*--------------------------------------------------------------------*/
/// Given coordinate in mapped space, return its location in real space
/**
 *   \param[in]  a_Xi   Location in computational space
 *   \return            Location in physical space
 *//*-----------------------------------------------------------------*/
RealVect
OgenCS::realCoordRaw(const RealVect& a_Xi) const
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
OgenCS::realCoord(const RealVect& a_Xi) const
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
OgenCS::globalInverseSearch(const RealVect& a_x,
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
OgenCS::inverseSearch(const RealVect& a_x,
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
/// Check if inverse is in domain
/**
 *   \param[in]  a_x     Location in physical space
 *   \param[in]  a_xiEst Initial guess for the inverse if feasible
 *   \param[out] a_xiEst A point near the inverse
 *//*-----------------------------------------------------------------*/
bool
OgenCS::inBlockDomain(const RealVect &a_x,
                          RealVect & a_xiEst) const
{
  CH_TIME("ExernalCS::inBlockDomain");
  RealVect baseDx(RealVect::Unit); // the factory always specifies unit base levels
  int refRatio = baseDx[0]/m_dx[0]; // this should work out exactly
  // convert a real valued point within a cell to an integer cell index
  auto convertToCellIdx = [&](const RealVect& a_xi)
  {
    //RealVect cellVal = (a_xi/m_dx);
    RealVect cellVal = (a_xi/baseDx);
    IntVect cellIdx;
    stc::forEachElement<SpaceDim>(
      [&] (const stc::array_size_type a_idx){
        cellIdx[a_idx] = std::floor(cellVal[a_idx]);
      });
    return cellIdx;
  };

  // do a coarse global search for the inverse, using an initial guess when able
  Box searchBox = m_domain.domainBox();
  searchBox.coarsen(refRatio);
  if (!searchBox.contains(convertToCellIdx(a_xiEst)))
    {
      a_xiEst = m_dx*0.5*(m_domain.domainBox().smallEnd() + m_domain.domainBox().bigEnd());
    }
  a_xiEst = inverseSearch(a_x, baseDx, searchBox, a_xiEst);
  IntVect nearIdx = convertToCellIdx(a_xiEst);
  // Some points are clearly bounded in the domain
  Box inside = grow(searchBox, -1);
  if (inside.contains(nearIdx))
    {
      return true;
    }
  // Points near the boundary are hard to detect
  else
    {
      // check if the point is feasibly in the domain
      // if the coarse inverse is not interior or adjacent to the boundary
      // region then it can not be contained.
      Box searchBoundary(nearIdx, nearIdx);
      searchBoundary.grow(2);
      auto xiBnd = inverseSearch(a_x, baseDx, searchBoundary, a_xiEst);
      searchBoundary.grow(-1);
      // Coarse point is adjacent to the boundary must be fully resolved
      if (searchBoundary.contains(convertToCellIdx(xiBnd)))
        {
          // initial guess is confirmed to be within 2 cells of the inverse
          a_xiEst = OgenCS::mappedCoord(a_x, xiBnd);
          // make sure the solution still lies within the bounded region
          // if this fails something has gone wrong with the initial point
          // for the newton solver, or the newton solve itself
          CH_assert(stc::mag(a_xiEst - xiBnd) < 2*stc::mag(baseDx));
          return searchBox.contains(convertToCellIdx(a_xiEst));
        }
      else
        {
          return false;
        }
    }
}

/*--------------------------------------------------------------------*/
/// Given coordinate in real space, return its location in mapped space
/**
 *   \param[in]  a_x    Location in physical space
 *   \return            Location in mapped space
 *//*-----------------------------------------------------------------*/
RealVect
OgenCS::mappedCoord(const RealVect& a_x) const
{
  CH_TIME("ExernalCS::mappedCoord");
  RealVect mapLoc;
  mapLoc = ogen.mappedCoord(a_x);
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
OgenCS::mappedCoord(const RealVect& a_x,
                        const RealVect& a_xiEst) const
{
  MayDay::Error("Real -> Mapped not defined!");
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
OgenCS::dXdXi(const RealVect& a_Xi, int a_dirX, int a_dirXi) const
{
  CH_TIME("ExernalCS::dXdXi");
  // set something invalid, which should get overwritten - returns NaN
  Real value = 1.0/0.0;;
  // make the periodic adjustment
  RealVect Xi(a_Xi);
  coordAdj(Xi);
  
  value = ogen.dXdXi(Xi, a_dirX, a_dirXi)*m_scale[a_dirXi]*m_scaleXi[a_dirXi];
  return value;
}

/*--------------------------------------------------------------------*/
/// Read in and store grid values. Compute all metrics needed
/**
 *//*-----------------------------------------------------------------*/

void OgenCS::readGridFile()
{
  CH_TIME("OgenCS::readGridFile");
  // point to the data base
  ogen.openFile(m_grid_file);
  // scale from Chombo domain size to Overture unit domain
  m_scaleXi = 1.0/(m_dx*m_domain.size());
}


// -- begin factory implementations ---------------------------

/*--------------------------------------------------------------------*/
///  Constructor
/**
 *   \param[in]  a_file_name
 *                      The name of the external file to read
 *   \param[in]  a_scale
 *                      Scaling to apply to the mesh
 *//*-----------------------------------------------------------------*/
OgenCSFactory::OgenCSFactory(const std::string& a_grid_file,
                                     const RealVect& a_scale):
  m_grid_file(a_grid_file),
  m_scale(a_scale)
{
}

#endif // #ifdef CH_USE_OGEN

/*--------------------------------------------------------------------*/
///  Return the OgenCS object
/**
 *   \param[in]  a_domain
 *                      Problem domain the coordinate system is in
 *   \return            Pointer to the coordinate system
 *//*-----------------------------------------------------------------*/
NewCoordSys*
OgenCSFactory::getCoordSys(const ProblemDomain& a_domain,
                               const RealVect& a_dx) const
{
  CH_assert(false);
  IntVect domSize = a_domain.size();
  RealVect domL;
  for (auto dir : EachDir)
    {
      domL[dir] = domSize[dir]*a_dx[dir];
    }
  OgenCS* newCSPtr = new OgenCS(a_dx, m_scale, a_domain, m_grid_file, 0);
  return static_cast< NewCoordSys* >(newCSPtr);
}

#include "NamespaceFooter.H"
