#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#if defined(CH_USE_FFTW) && defined(CH_MPI) && CH_SPACEDIM==3

// #define DEBUG_SPECTRALFILTER


/******************************************************************************/
/**
 * \file SpectralUtil.cpp
 *
 * \brief Member functions for SpectralUtil
 *
 *//*+*************************************************************************/

//----- Chombo -----//

#include "SpectralUtil.H"
#include "ShapeArray.H"


/*==============================================================================
 * Local helper functions
 *============================================================================*/


/*******************************************************************************
 *
 * Class LDOperatorAoSIn: member definitions
 *
 ******************************************************************************/

LDOperatorAoSIn<FArrayBox>::~LDOperatorAoSIn()
  { }

/// Linear in switches from SoA to AoS
void
LDOperatorAoSIn<FArrayBox>::linearIn(FArrayBox&      a_arg,
                                     void*           a_buf,
                                     const Box&      a_R,
                                     const Interval& a_comps) const
{
  CH_assert(Misc::isAligned<Real>(a_buf));
  const Real* buffer = (Real*)a_buf;
  const int n0 = a_arg.box().size(0);
  const int n1 = a_arg.box().size(1);
  const int n2 = a_arg.box().size(2);
  const int nC = a_arg.nComp();
  CH_assert(a_comps.begin() >= 0 && a_comps.end() < nC);
  // Access pattern
  // ptr + (i2-s2)*n1*n0*nC + (i1-s1)*n0*nC + (i0-s0)*nC + iC;
  // (-s2)*n1*n0*nC + (-s1)*n0*nC + (-s0)*nC
  auto arr = shape::make_array(a_arg.dataPtr(), n2, n1, n0, nC);
  arr.setRelativeLB(a_arg.smallEnd()[2],
                    a_arg.smallEnd()[1],
                    a_arg.smallEnd()[0],
                    0);
  // Real *const _arrdataPtr = a_arg.dataPtr() -
  //   ((a_arg.smallEnd()[2] *n1 +
  //     a_arg.smallEnd()[1])*n0 +
  //     a_arg.smallEnd()[0])*nC;
  // auto arr = (Real (*)[n1][n0][nC])_arrdataPtr;
  // // Assert this memory is used
  // assert(&(a_arg.operator()(a_arg.smallEnd(), 0)) ==
  //        &(arr[a_arg.smallEnd()[2]]
  //             [a_arg.smallEnd()[1]]
  //             [a_arg.smallEnd()[0]]
  //             [0]));
  // volatile Real arr_dbgval0 = arr[a_arg.smallEnd()[2]]
  //                                [a_arg.smallEnd()[1]]
  //                                [a_arg.smallEnd()[0]]
  //                                [0];
  // (void)arr_dbgval0;
  // volatile Real arr_dbgval1 = arr[a_arg.bigEnd()[2]]
  //                                [a_arg.bigEnd()[1]]
  //                                [a_arg.bigEnd()[0]]
  //                                [nC-1];
  // (void)arr_dbgval1;
  // End assert
  for (int iC = a_comps.begin(), iC_end = a_comps.end() + 1;
       iC != iC_end; ++iC)
    {
      MD_BOXLOOP(a_R, i)
        {
          arr[i2][i1][i0][iC] = *buffer++;
        }
    }
}

/// Local copy switches from SoA to AoS
void
LDOperatorAoSIn<FArrayBox>::op(FArrayBox&       a_dst,
                               const Box&       a_regionFrom,
                               const Interval&  a_Cdst,
                               const Box&       a_regionTo,
                               const FArrayBox& a_src,
                               const Interval&  a_Csrc) const
{
  // MD_ARRAY_RESTRICT(arrSrc, a_src);
  const int n0 = a_dst.box().size(0);
  const int n1 = a_dst.box().size(1);
  const int n2 = a_dst.box().size(2);
  const int nC = a_dst.nComp();
  CH_assert(a_Cdst.size() == a_Csrc.size() && a_regionFrom == a_regionTo);
  // Access pattern
  // ptr + (i2-s2)*n1*n0*nC + (i1-s1)*n0*nC + (i0-s0)*nC + iC;
  // (-s2)*n1*n0*nC + (-s1)*n0*nC + (-s0)*nC
  auto arrDst = shape::make_array(a_dst.dataPtr(), n2, n1, n0, nC);
  arrDst.setRelativeLB(a_dst.smallEnd()[2],
                       a_dst.smallEnd()[1],
                       a_dst.smallEnd()[0],
                       0);
  // Real *const _arrDstdataPtr = a_dst.dataPtr() -
  //   ((a_dst.smallEnd()[2] *n1 +
  //     a_dst.smallEnd()[1])*n0 +
  //     a_dst.smallEnd()[0])*nC;
  // auto arrDst = (Real (*)[n1][n0][nC])_arrDstdataPtr;
  // // Assert this memory is used
  // assert(&(a_dst.operator()(a_dst.smallEnd(), 0)) ==
  //        &(arrDst[a_dst.smallEnd()[2]]
  //                [a_dst.smallEnd()[1]]
  //                [a_dst.smallEnd()[0]]
  //                [0]));
  // volatile Real arrDst_dbgval0 = arrDst[a_dst.smallEnd()[2]]
  //                                      [a_dst.smallEnd()[1]]
  //                                      [a_dst.smallEnd()[0]]
  //                                      [0];
  // (void)arrDst_dbgval0;
  // volatile Real arrDst_dbgval1 = arrDst[a_dst.bigEnd()[2]]
  //                                      [a_dst.bigEnd()[1]]
  //                                      [a_dst.bigEnd()[0]]
  //                                      [nC-1];
  // (void)arrDst_dbgval1;
  // End assert
  for (int iCdst = a_Cdst.begin(), iCdst_end = a_Cdst.end() + 1,
         iCsrc = a_Csrc.begin(); iCdst != iCdst_end; ++iCdst, ++iCsrc)
    {
      MD_BOXLOOP(a_regionTo, i)
        {
          arrDst[i2][i1][i0][iCdst] = a_src[MD_IX(i, iCsrc)];
        }
    }
}


/*******************************************************************************
 *
 * Class LDOperatorAoSIn: member definitions
 *
 ******************************************************************************/

/// Destructor
LDOperatorAoSOut<FArrayBox>::~LDOperatorAoSOut()
  { }

/// Linear out switches from AoS to SoA
void
LDOperatorAoSOut<FArrayBox>::linearOut(const FArrayBox& a_arg,
                                       void*            a_buf,
                                       const Box&       a_R,
                                       const Interval&  a_comps) const
{
  CH_assert(Misc::isAligned<Real>(a_buf));
  Real* buffer = (Real*)a_buf;
  const int n0 = a_arg.box().size(0);
  const int n1 = a_arg.box().size(1);
  const int n2 = a_arg.box().size(2);
  const int nC = a_arg.nComp();
  CH_assert(a_comps.begin() >= 0 && a_comps.end() < nC);
  // Access pattern
  // ptr + (i2-s2)*n1*n0*nC + (i1-s1)*n0*nC + (i0-s0)*nC + iC;
  // (-s2)*n1*n0*nC + (-s1)*n0*nC + (-s0)*nC
  auto arr = shape::make_array(a_arg.dataPtr(), n2, n1, n0, nC);
  arr.setRelativeLB(a_arg.smallEnd()[2],
                    a_arg.smallEnd()[1],
                    a_arg.smallEnd()[0],
                    0);
  // const Real *const _arrdataPtr = a_arg.dataPtr() -
  //   ((a_arg.smallEnd()[2] *n1 +
  //     a_arg.smallEnd()[1])*n0 +
  //     a_arg.smallEnd()[0])*nC;
  // auto arr = (const Real (*)[n1][n0][nC])_arrdataPtr;
  // // Assert this memory is used
  // assert(&(a_arg.operator()(a_arg.smallEnd(), 0)) ==
  //        &(arr[a_arg.smallEnd()[2]]
  //             [a_arg.smallEnd()[1]]
  //             [a_arg.smallEnd()[0]]
  //             [0]));
  // volatile Real arr_dbgval0 = arr[a_arg.smallEnd()[2]]
  //                                [a_arg.smallEnd()[1]]
  //                                [a_arg.smallEnd()[0]]
  //                                [0];
  // (void)arr_dbgval0;
  // volatile Real arr_dbgval1 = arr[a_arg.bigEnd()[2]]
  //                                [a_arg.bigEnd()[1]]
  //                                [a_arg.bigEnd()[0]]
  //                                [nC-1];
  // (void)arr_dbgval1;
  // End assert
  for (int iC = a_comps.begin(), iC_end = a_comps.end() + 1;
       iC != iC_end; ++iC)
    {
      MD_BOXLOOP(a_R, i)
        {
          *buffer++ = arr[i2][i1][i0][iC];
        }
    }
}

/// Local copy switches from AoS to SoA
void
LDOperatorAoSOut<FArrayBox>::op(FArrayBox&       a_dst,
                                const Box&       a_regionFrom,
                                const Interval&  a_Cdst,
                                const Box&       a_regionTo,
                                const FArrayBox& a_src,
                                const Interval&  a_Csrc) const
{
  // MD_ARRAY_RESTRICT(arrDst, a_dst);
  const int n0 = a_src.box().size(0);
  const int n1 = a_src.box().size(1);
  const int n2 = a_src.box().size(2);
  const int nC = a_src.nComp();
  CH_assert(a_Cdst.size() == a_Csrc.size() && a_regionFrom == a_regionTo);
  // Access pattern
  // ptr + (i2-s2)*n1*n0*nC + (i1-s1)*n0*nC + (i0-s0)*nC + iC;
  // (-s2)*n1*n0*nC + (-s1)*n0*nC + (-s0)*nC
  auto arrSrc = shape::make_array(a_src.dataPtr(), n2, n1, n0, nC);
  arrSrc.setRelativeLB(a_src.smallEnd()[2],
                       a_src.smallEnd()[1],
                       a_src.smallEnd()[0],
                       0);
  // const Real *const _arrSrcdataPtr = a_src.dataPtr() -
  //   ((a_src.smallEnd()[2] *n1 +
  //     a_src.smallEnd()[1])*n0 +
  //     a_src.smallEnd()[0])*nC;
  // auto arrSrc = (const Real (*)[n1][n0][nC])_arrSrcdataPtr;
  // // Assert this memory is used
  // assert(&(a_src.operator()(a_src.smallEnd(), 0)) ==
  //        &(arrSrc[a_src.smallEnd()[2]]
  //                [a_src.smallEnd()[1]]
  //                [a_src.smallEnd()[0]]
  //                [0]));
  // volatile Real arrSrc_dbgval0 = arrSrc[a_src.smallEnd()[2]]
  //                                      [a_src.smallEnd()[1]]
  //                                      [a_src.smallEnd()[0]]
  //                                      [0];
  // (void)arrSrc_dbgval0;
  // volatile Real arrSrc_dbgval1 = arrSrc[a_src.bigEnd()[2]]
  //                                      [a_src.bigEnd()[1]]
  //                                      [a_src.bigEnd()[0]]
  //                                      [nC-1];
  // (void)arrSrc_dbgval1;
  // End assert
  for (int iCsrc = a_Csrc.begin(), iCsrc_end = a_Csrc.end() + 1,
         iCdst = a_Cdst.begin(); iCsrc != iCsrc_end; ++iCsrc, ++iCdst)
    {
      MD_BOXLOOP(a_regionFrom, i)
        {
          a_dst[MD_IX(i, iCdst)] = arrSrc[i2][i1][i0][iCsrc];
        }
    }
}

/*******************************************************************************
 *
 * Class FFTWDataFactory: member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
/// Destructor
/**
 *//*-----------------------------------------------------------------*/

FFTWDataFactory<FArrayBox>::~FFTWDataFactory()
{
  if (m_numCreate >= 0)  // Means define was called
    {
      fftw_free(m_data);
    }
}

/*--------------------------------------------------------------------*/
/// Weak construction
/** \param[in]  a_l_alloc
 *                      Number of complex variable points to allocate
 *                      with size fftw_complex
 *//*-----------------------------------------------------------------*/

void FFTWDataFactory<FArrayBox>::define(const int a_l_alloc)
{
  m_data = static_cast<Real*>(fftw_malloc(sizeof(fftw_complex)*a_l_alloc));
  m_numCreate = 0;
}

/*--------------------------------------------------------------------*/
/// Create an FArrayBox for the LevelData
/**
 *//*-----------------------------------------------------------------*/

FArrayBox* FFTWDataFactory<FArrayBox>::create(const Box&       a_box,
                                              int              a_ncomps,
                                              const DataIndex& a_datInd) const
{
  // Make sure the data is not aliased twice.  You could in theory clear the
  // LevelData and reallocate but we do not neet this.  So in practice data
  // for the LevelData can only be defined once.
  CH_assert(m_numCreate == 0);
  return new FArrayBox(a_box, a_ncomps, m_data);
  ++m_numCreate;
}

/*--------------------------------------------------------------------*/
/// Grab the date pointer (used when no box on this processor)
/** You cannot call both dataPtr and create
 *//*-----------------------------------------------------------------*/

Real* FFTWDataFactory<FArrayBox>::dataPtr() const
{
  CH_assert(m_numCreate == 0);
  return m_data;
  ++m_numCreate;
}


/*******************************************************************************
 *
 * Class SpectralFilter: member definitions
 *
 ******************************************************************************/

/*==============================================================================
 * Member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Define the class (weak construction)
/** \param[in]  a_problemDomain
 *                      Problem domain must be single-level periodic
 *                      rectangular domain
 *  \param[in]  a_boxes Disjoint box layout of the state to filter
 *  \param[in]  a_numComp
 *                      The FFTW plans are built for this many
 *                      components.  Generally, this should equal the
 *                      number of components in the LevelData passed
 *                      to the apply function.  Note that in apply
 *                      you can restrict filtering to a subset of
 *                      the components with the a_intv argument.
 *//*-----------------------------------------------------------------*/

void
SpectralFilter::define(const ProblemDomain&        a_problemDomain,
                       const DisjointBoxLayout&    a_boxes,
                       const int                   a_numComp,
                       const SpectralFilterProfile a_filterProfile)
{
  CH_TIMERS("SpectralFilter::define");
  // Spectral filtering might work with 32-bit floats but has not been tested
  assert(sizeof(Real) == 8);

  // Requirements:
  // 1) Single level with rectangular domain (we do not check that a_boxes
  //    cover the domain although this is a requirement).
  // 2) Periodic domain in all directions.
  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      CH_assert(a_problemDomain.isPeriodic(dir));
    }

  m_numComp = a_numComp;
  m_filterProfile = a_filterProfile;

  // FFTW can only transform AoS which does not match our standard data
  // layout that does not have components with unit stride.  Re-using a plan
  // (per component) requires matching alignments and we must accomodate
  // extra space needed by FFTW (as returned from *local_size*).

  // Approaches:
  // 1) MPI copy components 1 by 1 and transform (probably slowest
  //    alogorithms but only requires domain x 1 component of extra storage)
  // 2) Same as 1 but create a pipeline of copies and transforms.  Requires
  //    3x memory of 1): inbound, transform, outbound.  Gather messages may
  //    induce contention with FFTW interal messages.  Or maybe it is a
  //    good way to hide communication behind FFTW transform.
  // 3) Accommodate mpi_local_size in expanded k-dir and make sure either
  //    i or j have size 8 (for alignment).  MPI copy all components and
  //    then transform each component.  Requires extra 1x current state and
  //    no apparent advantage over 2).
  // 4) MPI copy all components and then do local copy to new array
  //    switching from SoA to AoS (probably fastest algorithms but requires
  //    extra 2x current state)
  // 5) MPI copy all components and switch from SoA to AoS in linear-in
  //    and vice-versa in linear-out (probably fastest algorithms, requires
  //    extra 1x current state, but quite confusing and all Chombo debugging
  //    tools go out the window for looking at the new layout)
  // The current algorithm implements 5) above because we are a bit masochistic.

  // fftw use C (row-ordered) storage.  We simply reverse the indices of our
  // array.

  // dims gives dimensions of the complex array.  Note modification to
  // dims[SpaceDim-1] to account for complex vars.  For the real array cast in
  // the same memory, r_n0 = 2*dims[SpaceDim-1] and others remain the same.

  // Prefix c means applied to complex variables only and r for real only.
  // If c and r are both absent, it means the value it is the same for both.
  m_c_numI = a_problemDomain.size(0)/2 + 1;
  m_numJ   = a_problemDomain.size(1);
  m_numK   = a_problemDomain.size(2);
  ptrdiff_t c_dims[SpaceDim];
  c_dims[SpaceDim - 1] = m_c_numI;
  for (int dir = 1; dir != SpaceDim; ++dir)
    {
      c_dims[SpaceDim - dir - 1] = a_problemDomain.size(dir);
    }

  // Except for dims, in the following, we refer to our normal ordering of
  // dimensions.  E.g., 0 = x-direction.  The last direction, SpaceDim - 1,
  // is denoted by K and the second last by J.
  // Remember that l_alloc is the number of complex variables that must be
  // locally allocated
  ptrdiff_t lc_alloc = fftw_mpi_local_size_many(
    SpaceDim, c_dims, a_numComp,
    FFTW_MPI_DEFAULT_BLOCK, MPI_COMM_WORLD, &m_l_numK, &m_l_idxKbeg);
  // Swap 0 and 1 indices to get transpose shape
  {
    ptrdiff_t tmp = c_dims[0];
    c_dims[0] = c_dims[1];
    c_dims[1] = tmp;
    ptrdiff_t dummy   = fftw_mpi_local_size_many(
      SpaceDim, c_dims, a_numComp,
      FFTW_MPI_DEFAULT_BLOCK, MPI_COMM_WORLD, &m_lt_numJ, &m_lt_idxJbeg);
    (void)dummy;
  }
  // Construct fftw boxes
  std::vector<int> all_numK(numProc());
  MPI_Allgather(&m_l_numK,       1, MPI_INT,
                all_numK.data(), 1, MPI_INT, Chombo_MPI::comm);
  // Count the number of non-empty boxes
  int fftwNumBox = 0;
  for (int iProc = 0, iProc_end = numProc(); iProc != iProc_end; ++iProc)
    {
      if (all_numK[iProc] > 0) ++fftwNumBox;
    }
  Vector<int> procIDs(fftwNumBox);
  Vector<Box> fftwBoxes(fftwNumBox);
  const int r_loI = a_problemDomain.domainBox().smallEnd(0);
  const int   loK = a_problemDomain.domainBox().smallEnd(SpaceDim - 1);
  int r_maxI = r_loI;
  int idxK = loK;
  for (int iProc = 0; iProc != fftwNumBox; ++iProc)
    {
      CH_assert(all_numK[iProc] > 0);
      procIDs[iProc] = iProc;
      Box box = a_problemDomain.domainBox();
      box.setBig(0, r_loI + 2*m_c_numI - 1);
      box.setSmall(SpaceDim - 1, idxK);
      idxK += all_numK[iProc];
      box.setBig  (SpaceDim - 1, idxK - 1);
      fftwBoxes[iProc] = box;
      if (procID() == iProc)
        {
          CH_assert(box.numPts() < 2*lc_alloc);
        }
      r_maxI = std::max(r_maxI, box.bigEnd(0));
    }
  CH_assert(idxK == m_numK);

  // We have to make the problem domain match the fftwBoxes.  No need for
  // periodic
  Box fftwDomainBox = a_problemDomain.domainBox();
  fftwDomainBox.setBig(0, r_maxI);
  ProblemDomain fftwProblemDomain(fftwDomainBox);
  DisjointBoxLayout dbl(fftwBoxes, procIDs, fftwProblemDomain);
  dbl.close();
  // l_alloc is used for the fabs in the LevelData and memory is allocated
  // with fftw_malloc (although our own alignment support would have been
  // sufficient)
  m_fftwFilterDataFactory.define(lc_alloc);
  m_fftwFilterData.define(dbl,
                          a_numComp,
                          IntVect::Zero,
                          m_fftwFilterDataFactory);
  Real* l_data;
  if (m_l_numK > 0)
    {
      // Grab the data from the fab
      l_data = m_fftwFilterData[m_fftwFilterData.dataIterator()].dataPtr();
    }
  else
    {
      // Directly get the data from the factory (there is no box on this
      // processor but fftw may want us to allocate some memory anyways).
      l_data = m_fftwFilterDataFactory.dataPtr();
    }
  fftw_complex *lc_data = reinterpret_cast<fftw_complex*>(l_data);

  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      // This is the size of the transform, both for real and complex.  Just
      // reusing c_dims.
      c_dims[SpaceDim - dir - 1] = a_problemDomain.size(dir);
    }
  CH_TIMER("SpectralFilter::fftwPlan", tPlan);
  CH_START(tPlan);
  m_fftwFilterForwardPlan  = fftw_mpi_plan_many_dft_r2c(
    SpaceDim, c_dims, a_numComp,
    FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, l_data, lc_data,
    MPI_COMM_WORLD, FFTW_MEASURE | FFTW_MPI_TRANSPOSED_OUT);
  m_fftwFilterBackwardPlan = fftw_mpi_plan_many_dft_c2r(
    SpaceDim, c_dims, a_numComp,
    FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, lc_data, l_data,
    MPI_COMM_WORLD, FFTW_MEASURE | FFTW_MPI_TRANSPOSED_IN);
  CH_STOP(tPlan);
  // In copier
  ProblemDomain notPeriodicDomain(a_problemDomain.domainBox());
  m_fftwFilterInCopier.define(a_boxes, dbl, notPeriodicDomain);
  // Out copier
  m_fftwFilterOutCopier = m_fftwFilterInCopier;
  m_fftwFilterOutCopier.reverse();
}

/*--------------------------------------------------------------------*/
//  Filter the state
/** \param[in]  a_gridResolution
 *                      Filterered grid resolution
 *  \param[in]  a_U     State vector to filter
 *  \param[in]  a_intv  Interval in state to filter
 *
 *  Note that the code transforms the entire state vector; using a
 *  smaller interval will not increase performance.
 *//*-----------------------------------------------------------------*/

void
SpectralFilter::apply(const IntVect&        a_gridResolution,
                      LevelData<FArrayBox>& a_U,
                      const Interval&       a_intv)
{
  CH_TIME("SpectralFilter::apply");

  // If not a sharp cutoff filter, filter parameter must be defined
  CH_assert(m_filterProfile == SpectralFilterSharp ||
            m_filterProfile == SpectralFilterSharpIsotropic ||
            m_filterParam != std::numeric_limits<Real>::infinity());

  CH_assert(a_intv.begin() >= 0 && a_intv.end() < a_U.nComp());
  const Box domainBox = a_U.disjointBoxLayout().physDomain().domainBox();
  const int numTotal = domainBox.numPts();
  const IntVect nD = domainBox.size();
  IntVect fw(D_DECL6(nD[0] - a_gridResolution[0] + 1,
                     nD[1] - a_gridResolution[1] + 1,
                     nD[2] - a_gridResolution[2] + 1,
                     nD[3] - a_gridResolution[3] + 1,
                     nD[4] - a_gridResolution[4] + 1,
                     nD[5] - a_gridResolution[5] + 1));
  CH_assert(fw <= domainBox.size());
  // There should only be 1 (or zero) fftw box per process
  CH_assert(m_fftwFilterData.dataIterator().size() <= 1);

#ifdef DEBUG_SPECTRALFILTER
  {
    Vector<int> oneProcIDs(1, 0);
    Vector<Box> oneBoxes(1, domainBox);
    ProblemDomain oneDomain(domainBox);
    DisjointBoxLayout oneDbl(oneBoxes, oneProcIDs, oneDomain);
    oneDbl.close();
    LevelData<FArrayBox> oneData(oneDbl, a_U.nComp(), IntVect::Zero);
    a_U.copyTo(oneData);
    if (procID() == 0)
      {
        FArrayBox& Ufab = oneData[oneData.dataIterator()];
        IntVect loBox(Ufab.smallEnd());
        IntVect hiBox(Ufab.bigEnd());
        IntVect lo(loBox);
        IntVect hi(hiBox);
        lo[2] = 0;
        hi[2] = 0;
        printf("original u, c=0:\n");
        dumpFAB2DSlicePretty(&Ufab, 0, lo, hi, 2, std::cout);
        printf("original u, c=1:\n");
        dumpFAB2DSlicePretty(&Ufab, 1, lo, hi, 2, std::cout);
        lo = loBox;
        hi = hiBox;
        lo[0] = 0;
        hi[0] = 0;
        printf("original u, c=2:\n");
        dumpFAB2DSlicePretty(&Ufab, 2, lo, hi, 2, std::cout);
        lo = loBox;
        hi = hiBox;
        lo[1] = 0;
        hi[1] = 0;
        printf("original u, c=3:\n");
        dumpFAB2DSlicePretty(&Ufab, 3, lo, hi, 2, std::cout);
        lo = loBox;
        hi = hiBox;
        lo[2] = 0;
        hi[2] = 0;
        printf("original u, c=4:\n");
        dumpFAB2DSlicePretty(&Ufab, 4, lo, hi, 2, std::cout);
      }
  }
#endif

  // Zero data in fftwFilterData if not filtering all components
  if (a_intv.size() < m_numComp && m_l_numK > 0)
    {
      m_fftwFilterData[m_fftwFilterData.dataIterator()].setVal(0.);
    }

  // Copy SoA layout into AoS layout for fftw
  a_U.copyTo(a_intv,
             m_fftwFilterData,
             a_intv,
             m_fftwFilterInCopier,
             LDOperatorAoSIn<FArrayBox>{});
  // Forward
  fftw_execute(m_fftwFilterForwardPlan);

//--Filter

  // You cannot use the BaseFab directly since the components have unit stride
  // This is a transposed array.  In 3-D, the c-storage-order complex-variable
  // dimensions are [m_lt_numJ][numK][c_numI][nC]
  const int c_numI = m_c_numI;  // Must be on stack for VLA-reshape
  const int numK = m_numK;
  const int nC = m_fftwFilterData.nComp();

  // Only processes with data participate in the filtering
  if (m_l_numK > 0)
    {
      // Reshape data into transposed form
      FArrayBox& fftwFilterFab =
        m_fftwFilterData[m_fftwFilterData.dataIterator()];
      auto ltc_data3d = shape::make_array(fftwFilterFab.dataPtr(),
                                          (int)m_lt_numJ, numK, c_numI, nC, 2);
      // auto ltc_data3d = (Real (*)[numK][c_numI][nC][2])fftwFilterFab.dataPtr();
      // // Assert this memory is used
      // assert(&(fftwFilterFab(fftwFilterFab.smallEnd(), 0)) ==
      //        &(ltc_data3d[0][0][0][0][0]));
      // volatile Real ltc_data3d_dbgval0 = (ltc_data3d[0][0][0][0][0]);
      // (void)ltc_data3d_dbgval0;
      // volatile Real ltc_data3d_dbgval1 =
      //   (ltc_data3d[m_lt_numJ-1][numK-1][c_numI-1][nC-1][1]);
      // (void)ltc_data3d_dbgval1;

#ifdef DEBUG_SPECTRALFILTER
      if (procID() == 0)
        {
          for (int i2 = 0; i2 != 1; ++i2)
            {
              printf("Proc 0, F(u) k=%d, c=%d:\n", i2, 1);
              // We print out the transposed array as is
              for (int i1 = 0; i1 != m_lt_numJ; ++i1)
                {
                  for (int i0 = 0; i0 != c_numI; ++i0)
                    {
                      printf(" (%6.2f,%6.2f)",
                             ltc_data3d[i1][i2][i0][1][0]/numTotal,
                             ltc_data3d[i1][i2][i0][1][1]/numTotal);
                    }
                  printf("\n");
                }
              printf("Proc 0, F(u) k=%d, c=%d:\n", i2, 4);
              // We print out the transposed array as is
              for (int i1 = 0; i1 != m_lt_numJ; ++i1)
                {
                  for (int i0 = 0; i0 != c_numI; ++i0)
                    {
                      printf(" (%6.2f,%6.2f)",
                             ltc_data3d[i1][i2][i0][4][0]/numTotal,
                             ltc_data3d[i1][i2][i0][4][1]/numTotal);
                    }
                  printf("\n");
                }
            }
        }
#endif

      if (m_filterProfile == SpectralFilterTanh)
        {
          if (((fw[0]/2) > 0) && ((fw[1]/2) > 0) && ((fw[2]/2) > 0))
            {
              const Real alpha = m_filterParam;
              const int beta = a_gridResolution[0]/2;
              for (int i1 = 0; i1 != m_lt_numJ; ++i1)
                for (int i2 = 0; i2 != numK; ++i2)
                  for (int i0 = 0; i0 != c_numI; ++i0)
                    for (int c = 0; c != nC; ++c)
                      {
                        int i1Actual = i1 + m_lt_idxJbeg;
                        int i1Shifted = std::abs(i1Actual - m_numJ);
                        int i2Shifted = std::abs(i2 - m_numK);
                        int waveNum0 = i0;
                        int waveNum1 = std::min(i1Actual, i1Shifted);
                        int waveNum2 = std::min(i2, i2Shifted);
                        const Real waveNumMag =
                          std::sqrt(waveNum0*waveNum0
                                    + waveNum1*waveNum1
                                    + waveNum2*waveNum2);
                        Real scalingFact =
                          0.5*std::tanh(-(waveNumMag - beta)/alpha) + 0.5;
                        if (waveNumMag <= 0.)
                          {
                            scalingFact = 1.;
                          }
                        ltc_data3d[i1][i2][i0][c][0] *= scalingFact;
                        ltc_data3d[i1][i2][i0][c][1] *= scalingFact;
                      }
            }
        }
      else if (m_filterProfile == SpectralFilterGaussian)
        {
          if (((fw[0]/2) > 0) && ((fw[1]/2) > 0) && ((fw[2]/2) > 0))
            {
              // beta and rolloff work together to shift half-Gaussian
              const int beta = a_gridResolution[0]/2;
              const Real rolloff = m_filterParam;
              const Real pi = std::atan(1.)*4.;
              const Real delta = pi/rolloff;
              const Real gamma = 3.56; // 0.5 spectral response at beta
              for (int i1 = 0; i1 != m_lt_numJ; ++i1)
                for (int i2 = 0; i2 != numK; ++i2)
                  for (int i0 = 0; i0 != c_numI; ++i0)
                    for (int c = 0; c != nC; ++c)
                      {
                        int i1Actual = i1 + m_lt_idxJbeg;
                        int i1Shifted = std::abs(i1Actual - m_numJ);
                        int i2Shifted = std::abs(i2 - m_numK);
                        int waveNum0 = i0;
                        int waveNum1 = std::min(i1Actual, i1Shifted);
                        int waveNum2 = std::min(i2, i2Shifted);
                        const Real waveNumMag =
                          std::sqrt(waveNum0*waveNum0
                                    + waveNum1*waveNum1
                                    + waveNum2*waveNum2);
                        Real scalingFact = std::exp(
                          -delta*delta*(waveNumMag - beta + rolloff)*
                          (waveNumMag - beta + rolloff)/(4.*gamma));
                        if (waveNumMag <= (beta - rolloff))
                          {
                            scalingFact = 1.;
                          }
                        ltc_data3d[i1][i2][i0][c][0] *= scalingFact;
                        ltc_data3d[i1][i2][i0][c][1] *= scalingFact;
                      }
            }
        }
      else if (m_filterProfile == SpectralFilterSharpIsotropic)
        {
          if (((fw[0]/2) > 0) && ((fw[1]/2) > 0) && ((fw[2]/2) > 0))
            {
              int isoCutoff = a_gridResolution[0]/2;
              for (int i1 = 0; i1 != m_lt_numJ; ++i1)
                for (int i2 = 0; i2 != numK; ++i2)
                  for (int i0 = 0; i0 != c_numI; ++i0)
                    for (int c = 0; c != nC; ++c)
                      {
                        int i1Actual = i1 + m_lt_idxJbeg;
                        int i1Shifted = std::abs(i1Actual - m_numJ);
                        int i2Shifted = std::abs(i2 - m_numK);
                        int waveNum0 = i0;
                        int waveNum1 = std::min(i1Actual, i1Shifted);
                        int waveNum2 = std::min(i2, i2Shifted);
                        const Real waveNumMag =
                          std::sqrt(waveNum0*waveNum0
                                    + waveNum1*waveNum1
                                    + waveNum2*waveNum2);
                        Real scalingFact = 1.;
                        if (waveNumMag > isoCutoff)
                          {
                            scalingFact = 0.;
                          }
                        ltc_data3d[i1][i2][i0][c][0] *= scalingFact;
                        ltc_data3d[i1][i2][i0][c][1] *= scalingFact;
                      }
            }
        }
      else if (m_filterProfile == SpectralFilterSharp)
        {
          // Direction 0 (half-symmetry)
          if ((fw[0]/2) > 0)
            {
              for (int i1 = 0; i1 != m_lt_numJ; ++i1)
                for (int i2 = 0; i2 != numK; ++i2)
                  for (int i0 = c_numI - (fw[0]/2); i0 != c_numI; ++i0)
                    for (int c = 0; c != nC; ++c)
                      {
                        ltc_data3d[i1][i2][i0][c][0] = 0.;
                        ltc_data3d[i1][i2][i0][c][1] = 0.;
                      }
            }
          // Direction 1 (full-symmetry, is partitioned in transposed form)
          if ((fw[1]/2) > 0)
            {
              const int isEven = 1 - m_numJ%2;
              const int i1_beg = std::max(
                (ptrdiff_t)0, (m_numJ/2) - (fw[1]/2-1) - m_lt_idxJbeg);
              const int i1_lst = std::min(
                m_lt_numJ  - 1, (m_numJ/2) + (fw[1]/2-isEven) - m_lt_idxJbeg);
              for (int i1 = i1_beg; i1 <= i1_lst; ++i1)
                for (int i2 = 0; i2 != numK; ++i2)
                  for (int i0 = 0; i0 != c_numI; ++i0)
                    for (int c = 0; c != nC; ++c)
                      {
                        ltc_data3d[i1][i2][i0][c][0] = 0.;
                        ltc_data3d[i1][i2][i0][c][1] = 0.;
                      }
            }
          // Direction 2 (full-symmetry)
          if ((fw[2]/2) > 0)
            {
              const int isEven = 1 - m_numK%2;
              const int i2_beg = (m_numK/2) - (fw[2]/2-1);
              const int i2_lst = (m_numK/2) + (fw[2]/2-isEven);
              for (int i1 = 0; i1 != m_lt_numJ; ++i1)
                for (int i2 = i2_beg; i2 <= i2_lst; ++i2)
                  for (int i0 = 0; i0 != c_numI; ++i0)
                    for (int c = 0; c != nC; ++c)
                      {
                        ltc_data3d[i1][i2][i0][c][0] = 0.;
                        ltc_data3d[i1][i2][i0][c][1] = 0.;
                      }
            }
        }

#ifdef DEBUG_SPECTRALFILTER
      if (procID() == 0)
        {
          for (int i2 = 0; i2 != 1; ++i2)
            {
              printf("Proc 0, filtered F(u) k=%d, c=%d:\n", i2, 1);
              // We print out the transposed array as is
              for (int i1 = 0; i1 != m_lt_numJ; ++i1)
                {
                  for (int i0 = 0; i0 != c_numI; ++i0)
                    {
                      printf(" (%6.2f,%6.2f)",
                             ltc_data3d[i1][i2][i0][1][0]/numTotal,
                             ltc_data3d[i1][i2][i0][1][1]/numTotal);
                    }
                  printf("\n");
                }
              printf("Proc 0, filtered F(u) k=%d, c=%d:\n", i2, 4);
              // We print out the transposed array as is
              for (int i1 = 0; i1 != m_lt_numJ; ++i1)
                {
                  for (int i0 = 0; i0 != c_numI; ++i0)
                    {
                      printf(" (%6.2f,%6.2f)",
                             ltc_data3d[i1][i2][i0][4][0]/numTotal,
                             ltc_data3d[i1][i2][i0][4][1]/numTotal);
                    }
                  printf("\n");
                }
            }
        }
#endif
    }  // Back to using all processors

  // Reverse
  fftw_execute(m_fftwFilterBackwardPlan);

  // Normalize
  if (m_l_numK > 0)
    {
      const int r_numI = 2*m_c_numI;  // Padded 0 dimension
      const int numJ = m_numJ;
      FArrayBox& fftwFilterFab =
        m_fftwFilterData[m_fftwFilterData.dataIterator()];
      auto l_data3d = shape::make_array(fftwFilterFab.dataPtr(),
                                        (int)m_l_numK, numJ, r_numI, nC);
      // auto l_data3d = (Real (*)[numJ][r_numI][nC])fftwFilterFab.dataPtr();
      // // Assert this memory is used
      // assert(&(fftwFilterFab(fftwFilterFab.smallEnd(), 0)) ==
      //        &(l_data3d[0][0][0][0]));
      // volatile Real l_data3d_dbgval0 = (l_data3d[0][0][0][0]);
      // (void)l_data3d_dbgval0;
      // volatile Real l_data3d_dbgval1 =
      //   (l_data3d[m_l_numK-1][numJ-1][r_numI-1][nC-1]);
      // (void)l_data3d_dbgval1;
      for (int i2 = 0; i2 != m_l_numK; ++i2)
        for (int i1 = 0; i1 != numJ; ++i1)
          // The following loop is only over used 0 dimension (not padded)
          for (int i0 = 0, i0_end = domainBox.size(0); i0 != i0_end; ++i0)
            for (int c = 0; c != nC; ++c)
              {
                l_data3d[i2][i1][i0][c] /= numTotal;
              }
    }

  // Copy AoS layout into standard SoA layout
  m_fftwFilterData.copyTo(a_intv,
                          a_U,
                          a_intv,
                          m_fftwFilterOutCopier,
                          LDOperatorAoSOut<FArrayBox>{});

#ifdef DEBUG_SPECTRALFILTER
  {
    Vector<int> oneProcIDs(1, 0);
    Vector<Box> oneBoxes(1, domainBox);
    ProblemDomain oneDomain(domainBox);
    DisjointBoxLayout oneDbl(oneBoxes, oneProcIDs, oneDomain);
    oneDbl.close();
    LevelData<FArrayBox> oneData(oneDbl, a_U.nComp(), IntVect::Zero);
    a_U.copyTo(oneData);
    if (procID() == 0)
      {
        FArrayBox& Ufab = oneData[oneData.dataIterator()];
        IntVect loBox(Ufab.smallEnd());
        IntVect hiBox(Ufab.bigEnd());
        IntVect lo(loBox);
        IntVect hi(hiBox);
        lo[2] = 0;
        hi[2] = 0;
        printf("filtered u, c=0:\n");
        dumpFAB2DSlicePretty(&Ufab, 0, lo, hi, 2, std::cout);
        printf("filtered u, c=1:\n");
        dumpFAB2DSlicePretty(&Ufab, 1, lo, hi, 2, std::cout);
        lo = loBox;
        hi = hiBox;
        lo[0] = 0;
        hi[0] = 0;
        printf("filtered u, c=2:\n");
        dumpFAB2DSlicePretty(&Ufab, 2, lo, hi, 2, std::cout);
        lo = loBox;
        hi = hiBox;
        lo[1] = 0;
        hi[1] = 0;
        printf("filtered u, c=3:\n");
        dumpFAB2DSlicePretty(&Ufab, 3, lo, hi, 2, std::cout);
        lo = loBox;
        hi = hiBox;
        lo[2] = 0;
        hi[2] = 0;
        printf("filtered u, c=4:\n");
        dumpFAB2DSlicePretty(&Ufab, 4, lo, hi, 2, std::cout);
      }
  }
#endif
}

#endif  /* ! (defined(CH_USE_FFTW) && defined(CH_MPI) && CH_SPACEDIM==3) */
