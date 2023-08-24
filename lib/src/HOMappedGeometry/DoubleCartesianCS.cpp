#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "DoubleCartesianCS.H"
#include "CartesianBlockCS.H"
#include "BoxIterator.H"

#include "NamespaceHeader.H"

DoubleCartesianCS::DoubleCartesianCS()
{
  m_coordSysVect.resize(NUMBLOCKS, NULL);

  m_mappingBlocks.resize(NUMBLOCKS);
}

DoubleCartesianCS::~DoubleCartesianCS()
{
  if (m_gotCoordSysVect)
    {
      for (int iblock = 0; iblock < NUMBLOCKS; iblock++)
        {
          delete m_coordSysVect[iblock];
        }
    }
}

void
DoubleCartesianCS::define(const ProblemDomain& a_levelDomain,
                          const RealVect& a_dx)
{
  const Box& levelBox = a_levelDomain.domainBox();
  int baseLength = levelBox.size(0) / 3;
  CH_assert(levelBox.size() == 3 * baseLength * IntVect::Unit);
  CH_assert(levelBox.smallEnd() == IntVect::Zero);

  IntVect baseLo = IntVect::Zero;
  IntVect baseHi = (baseLength - 1) * IntVect::Unit;
  Box baseBox(baseLo, baseHi);

  int blockIndex = 0;
  m_blockIndicesBox.define(IntVect::Zero, IntVect::Unit);
  m_blockIndices.define(m_blockIndicesBox, 1);
  BoxIterator bit(m_blockIndicesBox);
  for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      m_origin[blockIndex] = (2*baseLength) * iv;
      m_blockIndices(iv, 0) = blockIndex;
      blockIndex++;
    }

  for (int iblock = 0; iblock < NUMBLOCKS; iblock++)
    {
      m_mappingBlocks[iblock] = baseBox + m_origin[iblock];
    }
  m_gotMappingBlocks = true;

  /*
    Define block boundaries.
   */
  defineBoundaries();

  initializeBlockTransformations();

  for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      int iblock = m_blockIndices(iv, 0);
      const Box& bx = m_mappingBlocks[iblock];
      m_coordSysVect[iblock] = new CartesianBlockCS(iblock, iv, a_dx, bx);
    }

  m_gotCoordSysVect = true;
  setCache();
}

void
DoubleCartesianCS::defineBoundaries()
{
  CH_assert(gotMappingBlocks());
  // length in dimension 0 of central block (index 0)
  int baseLength = m_mappingBlocks[0].size(0);
  m_boundaries.resize(NUMBLOCKS);

  IndicesTransformation it;
  BoxIterator bit(m_blockIndicesBox);
  for (bit.begin(); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      int iblock = m_blockIndices(iv, 0);

      for (const auto dir : EachDir)
        {
          for (const auto side : EachSide)
            {
              BlockBoundary& bb = boundary(iblock, dir, side);
              const int sideSign = Side::sign(side);
              const int blockSideSign = 2*iv[dir] - 1;
              const int blockOther =
                m_blockIndices(iv - blockSideSign*BASISV(dir), 0);
              if (sideSign == blockSideSign)
                {
                  // Big jump across domain
                  it.defineFromTranslation((sideSign * -3 * baseLength)*
                                           BASISV(dir));
                  bb.define(it, blockOther);
                  // Periodic displacement (unit jump)
                  bb.definePhysicalTransformation(
                    RigidTransformation(
                      RealVect_zero,
                      RealVect((-sideSign)*RealVect_basis(dir))));
                }
              else
                {
                  // Next to each other but cross spacing
                  it.defineFromTranslation((sideSign * baseLength)*
                                           BASISV(dir));
                  bb.define(it, blockOther);
                  bb.definePhysicalTransformation(
                    RigidTransformation(RealVect_zero, RealVect_zero));
                }
            }
        }
    }

  m_gotBoundaries = true;
}

/*--------------------------------------------------------------------*/
//  Block mapping conversion function
/** \param[out] a_xi_valid
 *                      Same coordinate as a_xiSrc from the valid cell
 *                      in computational space
 *  \param[out] a_n_valid
 *                      Index of block containing valid cell
 *  \param[out] a_validExists
 *                      T - the ghost cell is within the domain and
 *                          has a valid cell
 *                      F - the ghost cell is outside the domain and
 *                          has a "nearest" valid cell
 *  \param[out] a_extraDispl
 *                      Extra displacement going from ghost cell to
 *                      valid cell, primarily due to periodic
 *                      boundaries
 *  \param[in]  a_xiSrc Coordinates of ghost cell in computational
 *                      space
 *  \param[in]  a_iSrc  IntVect of the ghost cell
 *  \param[in]  a_nSrc  Index of block containing the ghost cell
 *//*-----------------------------------------------------------------*/

void
DoubleCartesianCS::blockRemapping(RealVect&            a_xi_valid,
                                  int&                 a_n_valid,
                                  bool&                a_validExists,
                                  RigidTransformation& a_extraDispl,
                                  const RealVect&      a_xiSrc,
                                  const IntVect&       a_iSrc,
                                  int                  a_nSrc) const
{
  // This will be set to the X from the perspective of the valid cell
  RealVect X_valid = m_coordSysVect[a_nSrc]->realCoord(a_xiSrc);

  RealVect translation = RealVect_zero;
  // X_valid in [0:1]^D.
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (X_valid[idir] < 0.)
        {
          X_valid[idir] += 1.;
          translation[idir] = 1.;
        }
      else if (X_valid[idir] > 1.)
        {
          X_valid[idir] -= 1.;
          translation[idir] = -1.;
        }
    }
  a_extraDispl.define(RealVect_zero, translation);

  a_xi_valid = X_valid;
  IntVect iv_valid;  // Gives location of the block
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (X_valid[idir] > 0.5)
        {
          iv_valid[idir] = 1;
          a_xi_valid[idir] = X_valid[idir] + 0.5;
        }
      else
        {
          iv_valid[idir] = 0;
          a_xi_valid[idir] = X_valid[idir];
        }
    }
  a_n_valid = m_blockIndices(iv_valid, 0);
  a_validExists = true;
}

Vector<RealVect>
DoubleCartesianCS::displacements(const Vector<RealVect>&   a_dstCoords,
                                 const Vector<int>&        a_dstBlocks,
                                 const RealVect&           a_srcCoords,
                                 int                       a_srcBlock) const
{
  MayDay::Error("Is this used?");
  Vector<RealVect> disps =
    MultiBlockCoordSys::displacements(a_dstCoords, a_dstBlocks,
                                      a_srcCoords, a_srcBlock);
  int len = disps.size();
  for (int i = 0; i < len; i++)
    {
      for (int idir = 0; idir < SpaceDim; idir++)
        { // allow for wraparound in dimension idir
          if (disps[i][idir] > 1.)
            disps[i][idir] -= 1.;
          else if (disps[i][idir] < -1.)
            disps[i][idir] += 1.;
        }
    }
  return disps;
}


void
DoubleCartesianCS::separateVolFlux(LevelData<FluxBox>& a_flux) const
{
  CH_assert(false);
  const DisjointBoxLayout& layout = a_flux.disjointBoxLayout();
  DataIterator dit = layout.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      Box bx = layout[dit];
      int blockNum = whichBlock(bx);
      IntVect blockOrigin = m_origin[blockNum];

      FluxBox& flub = a_flux[dit];
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          if (blockOrigin[idir] > 0)
            { // not the first
              flub[idir].negate();
            }
        }
    }
}


MultiBlockCoordSys*
DoubleCartesianCSFactory::getCoordSys(const ProblemDomain& a_levelDomain,
                                      const RealVect& a_dx) const
{
  DoubleCartesianCS* coordSysPtr = new DoubleCartesianCS();
  if (m_doSaneDefine)
    {
      // Scale the problem domain to 3/2 of current to account for spacing
      // between blocks
      int baseLength = a_levelDomain.size(0)/2;
      CH_assert(a_levelDomain.domainBox().size() == 2*baseLength*IntVect_unit);
      ProblemDomain scaledDomain(Box(IntVect_zero,
                                     (3*baseLength - 1)*IntVect_unit));
      coordSysPtr->define(scaledDomain, a_dx);
    }
  else
    {
      coordSysPtr->define(a_levelDomain, a_dx);
    }
  return (static_cast <MultiBlockCoordSys*> (coordSysPtr));
}

#include "NamespaceFooter.H"
