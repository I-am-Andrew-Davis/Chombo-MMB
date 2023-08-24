#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "SingleBlockCSAdaptor.H"
#include "CartesianBlockCS.H"

#include "NamespaceHeader.H"

//-----------------------------------------------------------------------
SingleBlockCSAdaptor::
SingleBlockCSAdaptor(const NewCoordSys* const a_coordSys,
                     const ProblemDomain&     a_problemDomain,
                     const RealVect&          a_domainLength,
                     const Topology&          a_topology)
  :
  MultiBlockCoordSys(),
  m_topology(a_topology)
{
  CH_assert(a_coordSys != NULL);

  m_coordSysVect.resize(1);
  m_coordSysVect[0] = const_cast<NewCoordSys*>(a_coordSys);
  m_gotCoordSysVect = true;

  m_mappingBlocks.resize(1);
  m_mappingBlocks[0] = a_problemDomain.domainBox();
  m_gotMappingBlocks = true;

  m_boundaries.resize(1);
  Tuple<BlockBoundary, 2*SpaceDim>& blockBoundaries = m_boundaries[0];
  switch (a_topology)
    {
    case TopologyPeriodic:
    {
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          if (a_problemDomain.isPeriodic(idir))
            {
              IndicesTransformation idxTfm;
              BlockBoundary& boLo = blockBoundaries[idir];
              idxTfm.defineFromTranslation(
                IntVect_basis(idir)*a_problemDomain.size(idir));
              boLo.define(idxTfm, 0);
              boLo.defineConformal(0, true);
              const RealVect shiftLo =
                RealVect_basis(idir)*a_domainLength[idir];
              const RigidTransformation physTfmLo(RealVect::Zero, shiftLo);
              boLo.definePhysicalTransformation(physTfmLo);
              BlockBoundary& boHi = blockBoundaries[idir + SpaceDim];
              idxTfm.defineFromTranslation(
                IntVect_basis(idir)*(-a_problemDomain.size(idir)));
              boHi.define(idxTfm, 0);
              boHi.defineConformal(0, true);
              const RealVect shiftHi =
                RealVect_basis(idir)*(-a_domainLength[idir]);
              const RigidTransformation physTfmHi(RealVect::Zero, shiftHi);
              boHi.definePhysicalTransformation(physTfmHi);
            }
          else
            {
              blockBoundaries[idir].define(0);
              blockBoundaries[idir + SpaceDim].define(0);
            }
        }

      m_gotBoundaries = true;

      // Initialize block transformations
      initializeBlockTransformations();

      // Set the cache of block information.  We do this explicitly instead of
      // with the setCache function because we want a single-block problem
      // domain
      m_blockInfo.resize(1);
      BlockInfo& blockInfo = m_blockInfo[0];
      blockInfo = { a_problemDomain };
      break;
    }
    case TopologyDonutX:
    {
      // x-dir
      {
        IndicesTransformation idxTfm;
        BlockBoundary& boLo = blockBoundaries[0];
        idxTfm.defineFromTranslation(IntVect_basis(0)*a_problemDomain.size(0));
        boLo.define(idxTfm, 0);
        boLo.defineConformal(0, true);
        const RigidTransformation physTfm(RealVect_zero, RealVect_zero);
        boLo.definePhysicalTransformation(physTfm);
        BlockBoundary& boHi = blockBoundaries[0 + SpaceDim];
        idxTfm.defineFromTranslation(
          IntVect_basis(0)*(-a_problemDomain.size(0)));
        boHi.define(idxTfm, 0);
        boHi.defineConformal(0, true);
        boHi.definePhysicalTransformation(physTfm);
      }
      // y-dir
      blockBoundaries[1].define(0);
      blockBoundaries[1 + SpaceDim].define(0);
      // z-dir (periodic)
      if (SpaceDim > 2)
        {
          IndicesTransformation idxTfm;
          BlockBoundary& boLo = blockBoundaries[2];
          idxTfm.defineFromTranslation(
            IntVect_basis(2)*a_problemDomain.size(2));
          boLo.define(idxTfm, 0);
          boLo.defineConformal(0, true);
          const RealVect shiftLo = RealVect_basis(2)*a_domainLength[2];
          const RigidTransformation physTfmLo(RealVect::Zero, shiftLo);
          boLo.definePhysicalTransformation(physTfmLo);
          BlockBoundary& boHi = blockBoundaries[2 + SpaceDim];
          idxTfm.defineFromTranslation(
            IntVect_basis(2)*(-a_problemDomain.size(2)));
          boHi.define(idxTfm, 0);
          boHi.defineConformal(0, true);
          const RealVect shiftHi = RealVect_basis(2)*(-a_domainLength[2]);
          const RigidTransformation physTfmHi(RealVect::Zero, shiftHi);
          boHi.definePhysicalTransformation(physTfmHi);
        }

      m_gotBoundaries = true;

      // Initialize block transformations
      initializeBlockTransformations();

      setCache();
      break;
    }
    }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
SingleBlockCSAdaptor::
~SingleBlockCSAdaptor()
{
  delete m_coordSysVect[0];
}
//-----------------------------------------------------------------------

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
SingleBlockCSAdaptor::blockRemapping(RealVect&            a_xi_valid,
                                     int&                 a_n_valid,
                                     bool&                a_validExists,
                                     RigidTransformation& a_extraDispl,
                                     const RealVect&      a_xiSrc,
                                     const IntVect&       a_iSrc,
                                     int                  a_nSrc) const
{
  switch (m_topology)
    {
    case TopologyPeriodic:
      MayDay::Error("SingleBlockCSAdaptor with periodic topology is not "
                    "treated as multiblock.  So why was blockRemapping "
                    "called?");
      break;
    case TopologyDonutX:
      a_xi_valid = a_xiSrc;
      if (a_xi_valid[0] > 1.)
        {
          CH_assert(a_xi_valid[0] < 2.);
          a_xi_valid[0] -= 1.;
        }
      else if (a_xi_valid[0] < 0.)
        {
          CH_assert(a_xi_valid[0] > -1.);
          a_xi_valid[0] += 1.;
        }
      a_validExists = true;
      if (a_xi_valid[1] < 0.)
        {
          a_validExists = false;
          a_xi_valid[1] = 0.;
        }
      else if (a_xi_valid[1] > 1.)
        {
          a_validExists = false;
          a_xi_valid[1] = 1.;
        }
      a_n_valid = 0;
      a_extraDispl.defineAsIdentity();
      break;
    }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
SingleBlockCSAdaptorFactory::
SingleBlockCSAdaptorFactory(
  const NewCoordSysFactory* const      a_coordSysFactory,
  const RealVect                       a_domainLength,
  const SingleBlockCSAdaptor::Topology a_topology)
  :
  m_coordSysFactory(a_coordSysFactory),
  m_domainLength(a_domainLength),
  m_topology(a_topology)
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
SingleBlockCSAdaptorFactory::
~SingleBlockCSAdaptorFactory()
{
  delete m_coordSysFactory;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
MultiBlockCoordSys*
SingleBlockCSAdaptorFactory::
getCoordSys(const ProblemDomain& a_levelDomain,
            const RealVect& a_dx) const
{
  SingleBlockCSAdaptor* coordSysPtr = new SingleBlockCSAdaptor(
    m_coordSysFactory->getCoordSys(a_levelDomain, a_dx),
    a_levelDomain,
    m_domainLength,
    m_topology);
  return (static_cast <MultiBlockCoordSys*> (coordSysPtr));
}
//-----------------------------------------------------------------------

#include "NamespaceFooter.H"
