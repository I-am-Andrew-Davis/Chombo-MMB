#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "ExternalMultiBlockCS.H"
#include "ExternalCS.H"
#include "DebugOut.H"
#include "ReadCGNS.H"
#include <queue>

#include "NamespaceHeader.H"

ExternalMultiBlockCS::ExternalMultiBlockCS()
{
  m_gotExternalInfo = false;
}

ExternalMultiBlockCS::~ExternalMultiBlockCS()
{
  if (m_gotCoordSysVect)
    {
      for (int iblock = 0; iblock < m_numBlocks; iblock++)
        {
          delete m_coordSysVect[iblock];
        }
    }
}

void
ExternalMultiBlockCS::define(const ProblemDomain& a_levelDomain,
                             const RealVect& a_dx)
{
  CH_TIME("ExernalMultiBlockCS::define");
  if (m_verbosity > 1)
    {
      pout() << "ExternalMultiBlockCS begin defining coordinate system" << std::endl;
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
  m_numBlocks = cgns.numZone();
  m_coordSysVect.resize(m_numBlocks, nullptr);
  m_mappingBlocks.resize(m_numBlocks);
  m_dxVect = a_dx;

  // the base level must always be 1, as set in the constructor
  int refRatio = 1./a_dx[0]; 
  pout() << "REFRATIO: " << refRatio << std::endl;
  pout() << "DX: " << a_dx << std::endl;

  // Set periodicity from problem domain
  for (auto dir : EachDir)
    {
      m_periodicDir[dir] = a_levelDomain.isPeriodic(dir);
    }

  // !!!FIXME!!! this should be an argument
  IntVect numGhost = 8*IntVect::Unit;
  
  IntVect blockOrigin(IntVect::Zero);

  // Select procs to distribute ExternalCS (single block CS) construction
  const int procPerBlock = numProc()/m_numBlocks;
  int startProc = 0;

  for (int idxBlk = 0; idxBlk != m_numBlocks; ++idxBlk)
    {
      // Construct the box defining the block
      std::string zoneName;
      IntVect zoneNumCell;
      cgns.readZoneInfo(idxBlk + 1, zoneName, zoneNumCell);
      Box blockBox(IntVect_zero, zoneNumCell - 1);
      blockBox.refine(refRatio);
      m_mappingBlocks[idxBlk] = blockBox + blockOrigin;

      pout() << "  adding block(" << zoneName << ", " << idxBlk << "): "
             << m_mappingBlocks[idxBlk] << std::endl;
      // Get the mapping for each block (note: idxBlk = idxZone - 1)
      m_coordSysVect[idxBlk] = new ExternalCS(cgns,
                                              a_dx,
                                              m_scale,
                                              m_mappingBlocks[idxBlk],
                                              m_fileName,
                                              zoneName,
                                              idxBlk,
                                              numGhost,
                                              startProc,
                                              procPerBlock);
      // allow solution caching for the coarse grid only
      // provides speed up, but memory overhead is likely too high for finer grids
      // TODO Enabling for all levels requires the ability to clear cache during
      // regrid, but mapping in uninvolved in regrid
      if (refRatio == 1)
        {
          auto map = dynamic_cast<ExternalCS*>(m_coordSysVect[idxBlk]);
          map->enableCache();
        }
      
      // increment next computational start point
      blockOrigin[0] += blockBox.size()[0] + refRatio*3*numGhost[0];

      // increment proc
      if (procPerBlock > 0)
        {
          startProc += procPerBlock;
        }
      else
        {
          startProc++;
          startProc = startProc % numProc();
        }

    }
  
  m_gotCoordSysVect = true;
  m_gotMappingBlocks = true;
  
  defineBoundaries(cgns);
  initializeBlockTransformations();
  
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
ExternalMultiBlockCS::defineBoundaries(const ReadCGNS& a_cgns)
{
  CH_TIME("ExernalMultiBlockCS::defineBoundaries");
  CH_assert(gotMappingBlocks());
  m_boundaries.resize(m_numBlocks);

#ifdef CH_USE_CGNS
  IndicesTransformation idxTfm;
  RigidTransformation physTfm;
  // loop over each block, and corresponding zone 
  for (int idxBlk = 0; idxBlk != m_numBlocks; ++idxBlk)
    {
      // index for cgns
      int idxZone = idxBlk + 1;
      // start with setting all boundaries to type BOUNDARY
      for (const int dir : EachDir)
        {
          for (auto side : EachSide)
            {
              boundary(idxBlk, dir, side).define(0); // a physical boundary
            }
        }
      // check each face for connections and update accordingly
      std::vector<int> connectedIdxBlk;
      std::vector<Box> facesFrom;
      std::vector<Box> facesTo;
      std::vector<IntVect> trfms;
      a_cgns.readZoneConnections(
        idxZone, m_zoneNameMap, connectedIdxBlk, facesFrom, facesTo, trfms);

      // count the connected boundaries - there should be a least one if this is a multiblock grid
      if (connectedIdxBlk.size() == 0)
        {
          std::string zoneName;
          IntVect zoneNumCell;
          a_cgns.readZoneInfo(idxZone, zoneName, zoneNumCell);
          std::string errMsg = "Block(" + zoneName + ", " + std::to_string(idxBlk) + ") has no connected boundaries - check the CGNS grid";
          MayDay::Error(errMsg.c_str());
        }
      else if (m_verbosity > 2)
        {
          std::string zoneName;
          IntVect zoneNumCell;
          a_cgns.readZoneInfo(idxZone, zoneName, zoneNumCell);
          pout() << "Block(" << zoneName << ", " << idxBlk << ") has " << connectedIdxBlk.size() << " connections." << std::endl;
        }
      // Set up the connections
      for (int con = 0; con!=connectedIdxBlk.size(); ++con)
        {
          const int connIdx = connectedIdxBlk[con];
          // set up the index transformation
          //   separate out the permutation and sign
          IntVect permute; // The index that connects the other block (0:SpaceDim)
          IntVect sign; // The direction of each index (-1 or 1)
          const IntVect& T = trfms[con];
          stc::forEachElement<SpaceDim>(
            [&](const stc::array_size_type a_idx)
            {
              // permute[a_idx] = abs(T[a_idx]) - 1;
              int p = abs(T[a_idx]) - 1;
              permute[p] = a_idx;
              sign[p] = (T[a_idx] > 0) ? 1 : ((T[a_idx] < 0) ? -1 : 0);
            });
          
          // get the side and direction from the cgns face index boxes
          //  cgns index starts at 1. Connections must be over the entire face
          auto getDirAndSide = [&]
            (const Box& a_box, int& a_dir, Side::LoHiSide& a_side)
            {
              IntVect diff = a_box.bigEnd() - a_box.smallEnd();
              for (a_dir = 0; a_dir < SpaceDim; ++a_dir)
                {
                  // find the a_direction of the face
                  if (diff[a_dir] == 0)
                    {
                      // get the side
                      a_side = ((a_box.smallEnd(a_dir) == 0)
                                ? Side::Lo : Side::Hi);
                      return; //  and exit loop
                    }
                }
              MayDay::Error("getDirandSide failed, check that a flat box was given");
            };

          int thisDir, otherDir;
          Side::LoHiSide thisSide, otherSide;
          getDirAndSide(facesFrom[con], thisDir, thisSide);
          getDirAndSide(facesTo[con], otherDir, otherSide);
          // fix the permutation sign in the normal direction, apparently
          // grid generation sources can not be trusted for this
          if (thisSide != otherSide)
            {
              sign[otherDir] = 1;
            }
          else
            {
              sign[otherDir] = -1;
            }

          if (m_verbosity > 2)
            {
              pout() << "  connection dir " <<  thisDir
                     << " side " << thisSide
                     << " to block " << connIdx
                     << " dir " << otherDir
                     << " side " << otherSide
                     << " :" << std::endl;
            }

          // The connected faces
          Box thisFace = (m_mappingBlocks[idxBlk] - m_mappingBlocks[idxBlk].smallEnd()).surroundingNodes();
          thisFace.setRange(thisDir,
                            (thisSide==Side::Hi) ? thisFace.size(thisDir)-1 : 0);
          Box otherFace = (m_mappingBlocks[connIdx] - m_mappingBlocks[connIdx].smallEnd()).surroundingNodes();
          otherFace.setRange(otherDir,
                             (otherSide==Side::Hi) ? otherFace.size(otherDir)-1 : 0);
          // start of each block in global space
          IntVect thisBlkOrg = m_mappingBlocks[idxBlk].smallEnd();
          IntVect otherBlkOrg = m_mappingBlocks[connIdx].smallEnd();
          // translation portion from block origin to connected face
          IntVect thisBlkRelIdx = thisFace.smallEnd();
          IntVect otherBlkRelIdx = otherFace.smallEnd();
          for (const auto dir : EachDir)
            {
              if (sign[dir] == -1)
                otherBlkRelIdx[dir] = otherFace.bigEnd(dir);
              //otherBlkRelIdx[dir] = facesTo[con].bigEnd(dir);
            }
          // the global shift to access each connected face
          IntVect thisOffset = thisBlkOrg + thisBlkRelIdx;
          IntVect otherOffset = otherBlkOrg + otherBlkRelIdx;
          // define translation for only pivot - this will get fully defined later
          idxTfm.define(permute, sign, RealVect::Zero);
          // solve the translate two contributions
          //   "undo" the offset of for the face of this block
          //   shift to the face of the other block
          IntVect translate = otherOffset - idxTfm.transformNode(thisOffset);

          // define the MB transformation
          auto & thisboundary = boundary(idxBlk, thisDir, thisSide);
          idxTfm.define(permute, sign, translate);
          thisboundary.define(idxTfm, connIdx);
          if (m_verbosity > 2)
            pout() << "  Index permutation " << permute
                   << " sign " << sign
                   << " translation " << translate << std::endl;

          // set prescribed periodic boundaries
          if (m_periodicDir[thisDir])
            {
              // This means the the repeated domain follows the same mapping.
              // 'true' means it is also periodic (the domain repeats)
              thisboundary.defineConformal(connIdx, true);
            }
          
          // Get Xi on this face and its counterpart, and measure the distance
          // This should be zero for non-periodic connections - but MMB doesn't
          // distinguish periodic boundaries from normal connected ones
          // Fist get the connected face
          Box thisBox = thisFace + thisBlkOrg;
          RealVect thisXi = thisBox.smallEnd(); // just need any point in thisBox
          RealVect otherXi = idxTfm.transformNode(thisXi);
          //   get X on this face and its counterpart
          RealVect thisX = m_coordSysVect[idxBlk]->realCoord(m_dxVect*thisXi);
          RealVect otherX = m_coordSysVect[connIdx]->realCoord(m_dxVect*otherXi);
          //    compute the total physical offset
          RealVect physShift = otherX - thisX;
          if (m_verbosity > 2)
            pout() << "  Physical shift " << physShift << std::endl;

          // the connections are setup correctly, they should match exactly on the face
          std::vector<IntVect> checkDirs(SpaceDim+1);
          for (auto testD : EachDir)
            {
              checkDirs[testD] = BASISV(testD);
              // checkDirs[testD] = sign(thisSide)*BASISV(testD);
            }
          checkDirs[SpaceDim] = IntVect::Zero;
          for (auto shiftIV : checkDirs)
          {
            RealVect thisTestXi = thisXi + shiftIV;
            RealVect otherTestXi = idxTfm.transformNode(thisTestXi);
            //   get X on this face and its counterpart
            RealVect thisTestX = m_coordSysVect[idxBlk]->realCoord(m_dxVect*thisTestXi);
            RealVect otherTestX = m_coordSysVect[connIdx]->realCoord(m_dxVect*otherTestXi);
            // check the points match within a tolerance

            //pout() << "  Bnd xi " << thisTestXi << " -> " << otherTestXi << std::endl;
            //pout() << "  Bnd x  " << thisTestX << " == " << otherTestX << std::endl;
            const RealVect shiftErr = abs(thisTestX - (otherTestX - physShift));
            const RealVect shiftTol = max(abs(thisTestX), abs(otherTestX))*1e-14*RealVect::Unit;
            const bool shiftInterior = (BASISV(thisDir) == shiftIV) || (m_periodicDir[thisDir]);
            if ((shiftErr > shiftTol) && (!shiftInterior))
              {
                pout() << " tol " << shiftErr
                       << " > " << shiftTol
                       << std::endl;
                pout() << "  Bnd xi " << thisTestXi << " -> " << otherTestXi << std::endl;
                pout() << "  Bnd x  " << thisTestX << " == " << otherTestX << std::endl;
                MayDay::Error("Connected block boundary points do not align at physical locations");
                // If you get here something with the mapping has broken check:
                // 1) Does your grid actually connect properly
                // 2) the tolerances on this may be too tight
                // 3) the Chombo internal mapping is borked, good luck debugging
                //    - first check the permutations/sign
                //    - hope the mapping didn't break
              }
            // matching on the interior implies they conform
            // good starting point, but do something more robust than check 1 cell
            // else if ((shiftErr < shiftTol) && (shiftInterior))
            //   {
            //     MayDay::Warning("Automatically detected a conformal boundary");
            //     thisboundary.defineConformal(connIdx, true);
            //   }
          }

          // FIXME!! This is only valid for connected or simple rectangular periodic domains
          //         need to compute the physical rotation for more complicated geometry
          RealVect physRotate(RealVect::Zero);
          physTfm.define(physRotate, physShift);
          thisboundary.definePhysicalTransformation(physTfm);
        }
    }

  m_gotBoundaries = true;
#else
  MayDay::Error("ExternalMultiBlockCS requires CGNS support");
#endif // #ifdef USE_CGNS
}

/**
 *  \param[out] a_xi_valid
 *                      Mapped coordinate of valid cell
 *  \param[out] a_n_valid
 *                      Block containing valid cell
 *  \param[out] a_validExists
 *                      T : The ghost cells overlaps a valid cell
 *                      F : The ghost cells is either outside the
 *                          domain or the boundary is conformal (the
 *                          latter is a hack).  For the former, a
 *                          nearest valid cell is provided.
 *  \param[out] a_extraDispl
 *                      Physical displacements to get to the valid
 *                      cell
 *  \param[in]  a_xiSrc Mapped coordinate of the ghost cell
 *  \param[in]  a_iSrc  IntVect of the ghost cell
 *  \param[in]  a_nSrc  Block containing the ghost cell
 */
void
ExternalMultiBlockCS::blockRemapping(RealVect&            a_xi_valid,
                                     int&                 a_n_valid,
                                     bool&                a_validExists,
                                     RigidTransformation& a_extraDispl,
                                     const RealVect&      a_xiSrc,
                                     const IntVect&       a_iSrc,
                                     int                  a_nSrc) const
{
  CH_TIME("ExternalMultiBlockCS::blockRemapping");
  if (m_verbosity > 4)
    {
      pout() << "  ExternalMultiBlockCS blockRemapping: from block ="
             << a_nSrc << ", xi = " << a_xiSrc << std::endl;
    }

#if 0 // begin grid checking 
    // check if cell is valid
    Real ptJ = m_coordSysVect[a_nSrc]->pointwiseJ(a_xiSrc);
    if (ptJ <= 0)
      {
        pout() << "Bad J at ghost " << a_xiSrc << " of block " << a_nSrc << std::endl;
        MayDay::Error("Grid has invalid jacobian");
      }
#endif // end grid checking 
    
  bool followGridLines = false;
  a_n_valid = -1;
  // get source location in physical space
  RealVect X_orig = m_coordSysVect[a_nSrc]->realCoord(a_xiSrc);
  // adjust for "periodic" boundaries. I.e. those with connections
  // but different physical locations
  //   Determine the side shifted over, and apply the corresponding
  //   physical transformation. Null transformation if simply connected
  Box cellDomain = problemDomain(a_nSrc).domainBox();
  std::vector<BlockBoundary::btype> bdrysCrossed;
  for (auto side : EachSide)
    {
      // RealVect boxSideXi = nodeDomain.sideEnd(side)*m_dxVect;
      for (auto dir : EachDir)
        {
          // check if a connected boundary is crossed
          auto & srcBoundary = boundary(a_nSrc, dir, side);
          if (sign(side)*(a_iSrc[dir] - cellDomain.sideEnd(side)[dir]) > 0 &&
              srcBoundary.isInterface())
            {
              bdrysCrossed.push_back(srcBoundary.type());
              // shift the physical vector over the boundary
              a_extraDispl += srcBoundary.getPhysTransformation();
              X_orig = srcBoundary.getPhysTransformation().transformFwd(X_orig);
              // give an initial guess for the remapping
//**FIXME check with Nate
              const NewCoordSys* baseCoordSysPtr = getCoordSys(a_nSrc);
              a_xi_valid = baseCoordSysPtr->centerMappedCoordinates(
                srcBoundary.getTransformation().transformFwd(a_iSrc));
              // a_xi_valid = srcBoundary.getTransformation().transformFwd(a_xiSrc/m_dxVect);
              // a_xi_valid += a_xiSrc - m_dxVect*IntVect(a_xiSrc/m_dxVect);
              a_n_valid = srcBoundary.neighbor();
              if (m_followGridLinesOffBC)
                {
                  // Boundary in the original block othogonal to 'dir' and
                  // normal to a physical boundary (we prefix with 'tan')
                  for (const int tanDir : EachDir)
                    {
                      if (tanDir != dir)
                        {
                          for (const auto tanSide : EachSide)
                            {
                              const BlockBoundary& bo =
                                boundary(a_nSrc, tanDir, tanSide);
                              // (1) Must be domain boundary (should actually be
                              // a wall but we don't know the BC type in the
                              // current implementation)
                              if (bo.isDomainBoundary())
                                {
                                  // (2) Must be within specified normal
                                  // distance of physical boundary
                                  const int iNrmDist = Side::sign(tanSide)*(
                                    cellDomain.sideEnd(tanSide)[tanDir] -
                                    a_iSrc[tanDir]);
                                  if (iNrmDist >= 0 &&
                                      iNrmDist < m_followGridLinesNumNormal)
                                    {
                                      // (3) Ghost must be within specified
                                      // distance from source block
                                      const int iGhostDist = Side::sign(side)*(
                                        a_iSrc[dir] -
                                        cellDomain.sideEnd(side)[dir]);
                                      if (iGhostDist > 0 &&
                                          iGhostDist <=
                                          m_followGridLinesNumGhost)
                                        {
                                          // (4) The same boundary in the
                                          // neighbor must be an interface
                                          // (i.e., the BC went away)
                                          const int nbrTanDir =
                                            srcBoundary.dirOther(tanDir);
                                          const Side::LoHiSide nbrTanSide =
                                            static_cast<Side::LoHiSide>(
                                              ((Side::sign(tanSide)*
                                                srcBoundary
                                                .reorientFace(tanDir)) + 1)/2);
                                          const BlockBoundary& nbrBo =
                                            boundary(a_n_valid,
                                                     nbrTanDir,
                                                     nbrTanSide);
                                          if (nbrBo.isInterface())
                                            {
                                              followGridLines = true;
                                            }
                                        }
                                    }
                                }
                            }  // Loop over tangential sides
                        }
                    }  // Loop over tangential  directions
                }

            }
        }
    }
  // This is the point X we need an inverse for
  const RealVect X = X_orig;
  
  bool mapConformal = false;
  if (bdrysCrossed.size() == 1)
    {
      mapConformal = (bdrysCrossed[0] == BlockBoundary::CONFORMAL);
    }
  if (m_verbosity > 4)
    {
      if (mapConformal)
        {
          pout() << "     -> is a conformal boundary" << std::endl;
        }
      else
        {
          pout() << "     -> is a mapped boundary" << std::endl;
        }
    }
  // conformal boundaries do not need to be searched
  if (mapConformal)
    {
      // a_n_valid and a_xi_valid should already be correct
      // pout() << "Conformal init " << a_n_valid << " " << a_xi_valid << std::endl;
      // Hack to get conformal boundaries using the high order machinery by treating
      // them as first order
      a_validExists = false;
      return;
    }

  if (bdrysCrossed.size() == 1 && followGridLines)
    // Requiring bdrysCrossed.size() == 1 means that ghost cells across block
    // edges will not follow grids lines.  The impact of this is not fully
    // understood but most problems seem to arise from face interpolation in
    // normal directions.
    {
      // To follow grids lines, the valid cell should _not_ be in the direct
      // neighbor block.  The problem is that it is across an edge or corner in
      // a block where we do not want to source ghost cell information.
      // Note: a_n_valid and a_xi_valid were set when followGridLines was set
      //       true
      auto map = static_cast<ExternalCS*>(m_coordSysVect[a_n_valid]);
      if (!map->inBlockDomain(X, a_xi_valid))
        {
          a_xi_valid = map->nearestDomainCell(X);
          a_validExists = true;
          // if (m_verbosity >= 4)
            {
              pout() << "Following grid lines for MB ghost " << a_iSrc
                     << " to valid " << a_xi_valid << " in block "
                     << a_n_valid << std::endl;
            }
          return;  
        }
    }
      
  // Create a queue of neighbor blocks to search, with no duplicate searches
  std::queue<int> searchBlockQueue;
  std::vector<bool> searchedBlocks(m_numBlocks, false);
  // Lambda to pop the top of queue, mark as searched, and add new neighbors
  auto incrementQueue = [&]()
    {
      auto nextBlock = searchBlockQueue.front();
      searchBlockQueue.pop();
      searchedBlocks[nextBlock] = true;
      for (auto side: EachSide)
        {
          for (auto dir : EachDir)
            {
              auto & nextBoundary = boundary(nextBlock, dir, side);
              auto neighbor = nextBoundary.neighbor();
              if ((!nextBoundary.isDomainBoundary()) && (!searchedBlocks[neighbor]))
                {
                  searchBlockQueue.push(neighbor);
                }
            }
        }
      if (searchBlockQueue.empty())
      {
      for (auto i=0; i < searchedBlocks.size(); ++i)
        {
          if (!searchedBlocks[i])
            {
             searchBlockQueue.push(i);
            }
        }
      }
      return nextBlock;
    };

  // initialize the queue, with either the
  if (a_n_valid > -1)
    {
      // add the connected neighbor, this is usually where remapping is
      searchBlockQueue.push(a_n_valid);
    }
  else
    {
      // start with the current block if nothing better is known
      searchBlockQueue.push(a_nSrc);
    }
  // searchBlockQueue.push(0);
  // searchBlockQueue.push(1);
  // searchBlockQueue.push(2);

  // store the nearest neighbor for the case when the remapping is invalid
  RealVect nearestValid_xi;
  int nearestBlock = -1;
  Real nearestDist = std::numeric_limits<Real>::max();
  
  // Solve which neighboring block contains the point in physical space
  //   Check each neighboring block, if the point is in block bounds
  //   we found the block and are done. Otherwise try the next neighbor,
  //   and add second neighbors to the search space
  a_n_valid = -1;
  while (!searchBlockQueue.empty())
    {
      auto blockIdx = incrementQueue();
      // casting is bad, but in this case it is known this is the mapping type
      auto map = static_cast<ExternalCS*>(m_coordSysVect[blockIdx]);
      // For mapped boundaries, provide a initial guess X for the valid cell
      if (map->inBlockDomain(X, a_xi_valid))
        {
          a_validExists = true;
          a_n_valid = blockIdx;
          // map the point back to computational space in the neighboring block
          // use the guess from checking the domain to speed this up
          a_xi_valid = map->mappedCoord(X, a_xi_valid);
          // the point was found, so we are done here
          return;
        }
      else 
        {
          auto this_X = map->realCoord(a_xi_valid);
          Real this_dist(stc::mag(this_X - X));
          if (this_dist < nearestDist)
            {
              nearestDist = this_dist;
              nearestValid_xi = a_xi_valid;
              nearestBlock = blockIdx;
            }
        }
    }
  
  // If a remapping was not found, then a_xiSrc must be a ghost cell
  // that extends outside of the physical domain.
  // Set a_xi_valid and a_n_valid to the nearest location that does lie
  // inside the domain.
  if (a_n_valid == -1)
    {
      // make sure the search queue did not get messed up
      for (auto e : searchedBlocks)
        {
          if (!e)
            {
              pout() << e << endl;
              MayDay::Error("Block remapping not found, with some blocks unsearched");
            }
        }
      // indicate the point is invalid
      // set to nearest physical neighbor in the domain
      a_validExists = false;
      a_n_valid = nearestBlock;
      auto map = static_cast<ExternalCS*>(m_coordSysVect[a_n_valid]);
      a_xi_valid = map->nearestDomainCell(X);
      if (m_verbosity > 3)
        {
          pout() << "No block remapping for "
                 << a_xiSrc << " blk " << a_nSrc
                 << ", using "
                 << a_xi_valid << " blk " << a_n_valid
                 << std::endl;
        }
      // // the simplest result is to use closest point in a_nSrc
      // a_n_valid = a_nSrc;
      // // and project the point into the box
      // a_xi_valid = a_xiSrc;
      // a_xi_valid = stc::min(a_xi_valid, domain.bigEnd()*m_dxVect);
      // a_xi_valid = stc::max(a_xi_valid, domain.smallEnd()*m_dxVect);
    }
}


Vector<RealVect>
ExternalMultiBlockCS::displacements(const Vector<RealVect>&   a_dstCoords,
                                    const Vector<int>&        a_dstBlocks,
                                    const RealVect&           a_srcCoords,
                                    int                       a_srcBlock) const
{
  auto len = a_dstCoords.size();
  CH_assert(a_dstBlocks.size() == len);

  Vector<RealVect> returnVec(len);
  // the source location
  auto srcX = m_coordSysVect[a_srcBlock]->realCoord(a_srcCoords);
  // copy and move the source over all adjusted boundaries, allows periodic
  Vector<RealVect> srcShifts{};
  for (auto side : EachSide)
    {
      for (auto dir : EachDir)
        {
          auto & srcBoundary = boundary(a_srcBlock, dir, side);
          auto X = srcBoundary.getPhysTransformation().transformFwd(srcX);
          if (X != srcX)
            {
              srcShifts.push_back(X);
            }
        }
    }
  // compute each of the the displacements
  for (int i = 0; i < len; i++)
    {
      auto dstX = m_coordSysVect[a_dstBlocks[i]]->realCoord(a_dstCoords[i]);
      // compute distance
      RealVect distVect = dstX - srcX;
      // check if one of the periodic shifts is smaller
      for (auto & X : srcShifts)
        {
          RealVect distVectShift = dstX - X;
          if (stc::mag(distVectShift) < stc::mag(distVect))
            {
              distVect = distVectShift;
            }
        }
      returnVec[i] = distVect;
    }
  return returnVec;
}

// void
// ExternalMultiBlockCS::separateVolFlux(LevelData<FluxBox>& a_flux) const
// {
//   // I don't think this needs to do anything
//   // const DisjointBoxLayout& layout = a_flux.disjointBoxLayout();
//   // DataIterator dit = layout.dataIterator();
//   // for (dit.begin(); dit.ok(); ++dit)
//   //   {
//   //     Box bx = layout[dit];
//   //     int blockNum = whichBlock(bx);

//   //     FluxBox& flub = a_flux[dit];
//   //     for (int idir = 0; idir < SpaceDim; idir++)
//   //       {
//   //         //if (blockOrigin[idir] > 0)
//   //           { // not the first
//   //             flub[idir].negate();
//   //           }
//   //       }
//   //   }
// }

void
ExternalMultiBlockCS::setupForExternalGrid(const std::string a_fileName,
                                           const RealVect a_scale)
{
  m_fileName = a_fileName;
  m_scale = a_scale;
  m_gotExternalInfo = true;
}

// -- begin factory implementations ---------------------------

/*--------------------------------------------------------------------*/
/// Create an ExternalMultiBlockCS multiblock coordinate system
/**
 *  \param[in]  a_levelDomain
 *                      Only used to get periodicity
 *  \param[in]  a_dx    Computational-grid mesh spacing on this level.
 *                      For this particular class, it is fixed that
 *                      a_dx = 1 on the base grid.
 *//*-----------------------------------------------------------------*/

MultiBlockCoordSys*
ExternalMultiBlockCSFactory::getCoordSys(const ProblemDomain& a_levelDomain,
                                         const RealVect&      a_dx) const
{
  ExternalMultiBlockCS* coordSysPtr = new ExternalMultiBlockCS();
  coordSysPtr->setupForExternalGrid(m_fileName, m_scale);
  coordSysPtr->setVerbosity(m_verbosity);
  coordSysPtr->define(a_levelDomain, a_dx);
    
  return (static_cast <MultiBlockCoordSys*> (coordSysPtr));
}

#include "NamespaceFooter.H"
