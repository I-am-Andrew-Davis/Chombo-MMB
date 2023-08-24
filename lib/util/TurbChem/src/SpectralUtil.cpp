#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#if defined(CH_USE_FFTW) && CH_SPACEDIM==3

// #define DEBUG_ESD


/******************************************************************************/
/**
 * \file SpectralUtil.cpp
 *
 * \brief Member functions for SpectralUtil
 *
 *//*+*************************************************************************/

//----- Standard Library -----//

#include <chrono>
#include <fstream>

//----- Chombo -----//

#include "RealVect.H"
#include "ShapeArray.H"

//----- Internal -----//

#include "SpectralUtil.H"


/*==============================================================================
 * Local helper functions
 *============================================================================*/


/*******************************************************************************
 */
/// Operator that switches from Struct-of-Arrays (SoA) to Array-of-Structs (AoS)
/**
 *  Chombo normally features components at largest stride in an array.  This
 *  switches to components having unit stride.  The resulting memory is used
 *  by fftw.
 *
 *//*+*************************************************************************/

template <typename T>
class LDOperatorAoSIn;
template <>
class LDOperatorAoSIn<FArrayBox> : public LDOperator<FArrayBox>
{
public:
  /// Destructor
  virtual ~LDOperatorAoSIn()
    { }
  /// Linear in switches from SoA to AoS
  virtual void linearIn(FArrayBox&      a_arg,
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
  virtual void op(FArrayBox&       a_dst,
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
};


/*******************************************************************************
 */
/// Operator that switches from Array-of-Structs (AoS) to Struct-of-Arrays (SoA)
/**
 *  Chombo normally features components at largest stride in an array.  This
 *  switches from components having unit stride to a normal Chombo layout.  The
 *  AoS layout was used by fftw.
 *
 *//*+*************************************************************************/

template <typename T>
class LDOperatorAoSOut;
template <>
class LDOperatorAoSOut<FArrayBox> : public LDOperator<FArrayBox>
{
public:
  /// Destructor
  virtual ~LDOperatorAoSOut()
    { }
  /// Linear out switches from AoS to SoA
  virtual void linearOut(const FArrayBox& a_arg,
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
  virtual void op(FArrayBox&       a_dst,
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
};


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
  if (m_numCreate >= 0)  // Clear any previous definition
    {
      fftw_free(m_data);
    }
  std::cout << "Local allocation size: "
            << static_cast<float>(sizeof(fftw_complex)*a_l_alloc)/1024
            << " kiB" << std::endl;;
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
 * Class EnergySpectralDensity: member definitions
 *
 ******************************************************************************/

/*==============================================================================
 * Constructors and Destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

EnergySpectralDensity::~EnergySpectralDensity()
{
  if (m_r_numI != -1)
    {
      fftw_destroy_plan(m_fftwForwardPlan);
    }
}

/*==============================================================================
 * Member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Define the class (weak construction)
/** \param[in]  a_problemDomain
 *                      Problem domain must be single-level periodic
 *                      rectangular domain
 *  \param[in]  a_boxes Disjoint box layout of the state to filter
 *//*-----------------------------------------------------------------*/

void
EnergySpectralDensity::defineSerial(const ProblemDomain& a_problemDomain)
{
  if (!m_use) return;
  // Spectral filtering might work with 32-bit floats but has not been tested
  assert(sizeof(Real) == 8);

  // Do nothing if problem domain size has not changed
  if (a_problemDomain.domainBox() == m_problemDomain.domainBox()) return;
  m_problemDomain = a_problemDomain;

  // Clear any previous definitions
  if (m_r_numI != -1)
    {
      fftw_destroy_plan(m_fftwForwardPlan);
    }

  // Requirements (but neither are checked for):
  // 1) Single level with rectangular domain.
  // 2) Periodic domain in all directions.

  // FFTW can only transform AoS which does not match our standard data
  // layout that does not have components with unit stride.  Re-using a plan
  // (per component) requires matching alignments.

  // FFTW use C (row-ordered) storage.  We simply reverse the indices of our
  // array.

  // dims gives dimensions of the complex array.  Note modification to
  // dims[SpaceDim-1] to account for complex vars.  For the real array cast in
  // the same memory, r_n0 = 2*dims[SpaceDim-1] and others remain the same.

  // Prefix c means applied to complex variables only and r for real only.
  // If c and r are both absent, it means the value it is the same for both.
  m_r_numI = a_problemDomain.size(0);
  m_c_numI = a_problemDomain.size(0)/2 + 1;
  m_numJ   = a_problemDomain.size(1);
  m_numK   = a_problemDomain.size(2);

  // Construct fftw boxes
  const int fftwNumBox = 1;
  Vector<int> procIDs(fftwNumBox);
  Vector<Box> fftwBoxes(fftwNumBox);
  const int r_loI = a_problemDomain.domainBox().smallEnd(0);
  int r_maxI = r_loI;
  {
    procIDs[0] = 0;
    Box box = a_problemDomain.domainBox();
    box.setBig(0, r_loI + 2*m_c_numI - 1);
    fftwBoxes[0] = box;
    r_maxI = std::max(r_maxI, box.bigEnd(0));
  }
  // We have to make the problem domain match the fftwBoxes.  No need for
  // periodic
  Box fftwDomainBox = a_problemDomain.domainBox();
  fftwDomainBox.setBig(0, r_maxI);
  ProblemDomain fftwProblemDomain(fftwDomainBox);
  DisjointBoxLayout dbl(fftwBoxes, procIDs, fftwProblemDomain);
  dbl.close();
  m_fftwDataFactory.define(m_c_numI*m_numJ*m_numK*SpaceDim);
  m_fftwData.define(dbl, SpaceDim, IntVect::Zero, m_fftwDataFactory);

  // Grab the data from the fab
  Real* l_data = m_fftwData[m_fftwData.dataIterator()].dataPtr();
  fftw_complex *lc_data = reinterpret_cast<fftw_complex*>(l_data);

  int dims[SpaceDim];
  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      // This is the size of the transform, both for real and complex.
      dims[SpaceDim - dir - 1] = a_problemDomain.size(dir);
    }

  using Clock = std::chrono::steady_clock;
  using TimePoint = typename Clock::time_point;
  TimePoint startTime = Clock::now();
  std::cout << "Developing FFTW plans..." << std::endl;
  /* The arrays are stored with the component index has unit stride (AoS
   * format).  For a single plan, this is the only option for MPI.  We use the
   * same approach here even though a serial plan could be created with a
   * contiguous array for each component (SoA format), consistent with the
   * standard Chombo layout
   */
  m_fftwForwardPlan  = fftw_plan_many_dft_r2c(
    SpaceDim, dims, SpaceDim,
    //             +-- between elements of same array
    //             |         +-- to the next array
    //             v         v
    l_data,  NULL, SpaceDim, 1,
    lc_data, NULL, SpaceDim, 1, FFTW_MEASURE);
  std::cout << "Time to develop plans\n"
            << std::chrono::duration_cast<std::chrono::duration<double>>(
              Clock::now() - startTime).count() << " s" << std::endl;
}

/*--------------------------------------------------------------------*/
//  Define the class (weak construction)
/** \param[in]  a_problemDomain
 *                      Problem domain must be single-level periodic
 *                      rectangular domain
 *  \param[in]  a_boxes Disjoint box layout of the state to filter
 *//*-----------------------------------------------------------------*/

void
EnergySpectralDensity::load(const Box&       a_box,
                            const FArrayBox& a_dataFab,
                            const Interval&  a_intv)
{
  if (!m_use) return;
  LDOperatorAoSIn<FArrayBox> op{};
  // Copy switches from SoA to AoS (i.e., components having unit stride)
  op.op(m_fftwData[m_fftwData.dataIterator()],
        a_box,
        Interval(0, SpaceDim-1),
        a_box,
        a_dataFab,
        a_intv);
}

/*--------------------------------------------------------------------*/
//  Filter the state
/** \param[in]  a_physDomain
 *                      Physical domain dimensions
 *  \param[in]  a_numSubSample
 *                      Sub-sampling in frequency domain.  Each
 *                      frequency bin is divided into this many cells
 *                      (in each direction) before accumulating into
 *                      energy shells.  Must be even number > 1
 *                      (default 4);
 *//*-----------------------------------------------------------------*/

void
EnergySpectralDensity::compute(const RealVect&      a_physDomain,
                               const int            a_numSubSample)
{
  if (!m_use) return;
  // There should only be 1 (or zero) fftw box per process
  CH_assert(m_fftwData.dataIterator().size() <= 1);
  // Sub-sample resolution must be > 1 and even
  CH_assert(a_numSubSample > 1 && a_numSubSample%2 == 0);

#ifdef DEBUG_ESD
  {
    const int r_numI = 2*m_c_numI;  // Must be on stack for VLA-reshape
    const int numJ = m_numJ;
    const int nC = SpaceDim;

    // Reshape data
    FArrayBox& fftwFab = m_fftwData[m_fftwData.dataIterator()];
    auto l_data3d = shape::make_array(fftwFab.dataPtr(),
                                      m_numK, numJ, r_numI, nC);
    // auto l_data3d = (Real (*)[numJ][r_numI][nC])fftwFab.dataPtr();
    // // Assert this memory is used
    // assert(&(fftwFab(fftwFab.smallEnd(), 0)) == &(l_data3d[0][0][0][0]));
    // volatile Real l_data3d_dbgval0 = (l_data3d[0][0][0][0]);
    // (void)l_data3d_dbgval0;
    // volatile Real l_data3d_dbgval1 =
    //   (l_data3d[m_numK-1][numJ-1][r_numI-1][nC-1]);
    // (void)l_data3d_dbgval1;

    for (int i2 = 0; i2 != 1; ++i2)
      {
        printf("u k=%d:\n", i2);
        for (int i1 = 0; i1 != numJ; ++i1)
          {
            for (int i0 = 0; i0 != m_r_numI; ++i0)
              {
                printf(" %6.2f",
                       l_data3d[i2][i1][i0][0]);
              }
            printf("\n");
          }
      }
    for (int i2 = 0; i2 != 1; ++i2)
      {
        printf("v k=%d:\n", i2);
        for (int i1 = 0; i1 != numJ; ++i1)
          {
            for (int i0 = 0; i0 != m_r_numI; ++i0)
              {
                printf(" %6.2f",
                       l_data3d[i2][i1][i0][1]);
              }
            printf("\n");
          }
      }
  }
#endif

  // Forward FFT
  fftw_execute(m_fftwForwardPlan);

  // In 3-D, the c-storage-order complex-variable dimensions are
  // [m_numK][m_numJ][m_c_numI][SpaceDim]
  const int c_numI = m_c_numI;  // Must be on stack for VLA-reshape
  const int numJ = m_numJ;
  const int nC = SpaceDim;

  // Reshape data
  FArrayBox& fftwFab = m_fftwData[m_fftwData.dataIterator()];
  auto lc_data3d = shape::make_array(fftwFab.dataPtr(),
                                     m_numK, numJ, c_numI, nC, 2);
  // auto lc_data3d = (Real (*)[numJ][c_numI][nC][2])fftwFab.dataPtr();
  // // Assert this memory is used
  // assert(&(fftwFab(fftwFab.smallEnd(), 0)) == &(lc_data3d[0][0][0][0][0]));
  // volatile Real lc_data3d_dbgval0 = (lc_data3d[0][0][0][0][0]);
  // (void)lc_data3d_dbgval0;
  // volatile Real lc_data3d_dbgval1 =
  //   (lc_data3d[m_numK-1][numJ-1][c_numI-1][nC-1][1]);
  // (void)lc_data3d_dbgval1;

#ifdef DEBUG_ESD
  {
    const int numTotal = m_r_numI*m_numJ*m_numK;
    for (int i2 = 0; i2 != 1; ++i2)
      {
        printf("F(u) k=%d:\n", i2);
        for (int i1 = 0; i1 != numJ; ++i1)
          {
            for (int i0 = 0; i0 != c_numI; ++i0)
              {
                printf(" (%6.2f,%6.2f)",
                       lc_data3d[i2][i1][i0][0][0]/numTotal,
                       lc_data3d[i2][i1][i0][0][1]/numTotal);
              }
            printf("\n");
          }
      }
    for (int i2 = 0; i2 != 1; ++i2)
      {
        printf("F(v) k=%d:\n", i2);
        for (int i1 = 0; i1 != numJ; ++i1)
          {
            for (int i0 = 0; i0 != c_numI; ++i0)
              {
                printf(" (%6.2f,%6.2f)",
                       lc_data3d[i2][i1][i0][1][0]/numTotal,
                       lc_data3d[i2][i1][i0][1][1]/numTotal);
              }
            printf("\n");
          }
      }
  }
#endif

//--Cartesian summation of energy and normalization

  const int c_numJ = m_numJ/2 + 1;  // Number of non-mirrored complex components
  const int c_numK = m_numK/2 + 1;

  const Real cnorm = a_physDomain.product()/
    std::pow(((Real)m_r_numI)*((Real)m_numJ)*((Real)m_numK), 2);

  for (int i2 = 0; i2 != c_numK; ++i2)
    for (int i1 = 0; i1 != c_numJ; ++i1)
      for (int i0 = 0; i0 != m_c_numI; ++i0)
        {
          Real umag2 = (std::pow(lc_data3d[i2][i1][i0][0][0], 2) +
                        std::pow(lc_data3d[i2][i1][i0][0][1], 2));
          Real vmag2 = (std::pow(lc_data3d[i2][i1][i0][1][0], 2) +
                        std::pow(lc_data3d[i2][i1][i0][1][1], 2));
          Real wmag2 = (std::pow(lc_data3d[i2][i1][i0][2][0], 2) +
                        std::pow(lc_data3d[i2][i1][i0][2][1], 2));
          lc_data3d[i2][i1][i0][0][0] =  // Now holds KE
            0.5*(umag2 + vmag2 + wmag2)*cnorm;
        }

//--Spherical summation into of energy (and spectrum).  Use sub-sample
//--resolution of a_numSubSample per frequency (think of sub-cells per cell).

  const int c_numIJKmin = std::min(std::min(m_c_numI, c_numJ), c_numK);

  // Number of sub-samples (sub-cells) in each direction
  const int halfNumSubSample = a_numSubSample/2;
  const int n0s = (m_c_numI-1)*a_numSubSample + (m_r_numI%2)*halfNumSubSample;
  const int n1s = (  c_numJ-1)*a_numSubSample + (  m_numJ%2)*halfNumSubSample;
  const int n2s = (  c_numK-1)*a_numSubSample + (  m_numK%2)*halfNumSubSample;

  // Sub-samples based on a sample (or wave) separation of 1
  const Real subSampleLength = 1./a_numSubSample;
  const Real subSampleVolume8 = std::pow(subSampleLength, SpaceDim)*8;

  m_ESD.resize(c_numIJKmin + 1);
  for (auto& esd : m_ESD) esd = 0.;

  for (int i2s = 0; i2s != n2s; ++i2s)  // Sub-sample index
    {
      const Real w2 = (i2s + 0.5)*subSampleLength;  // Wave position 2
#ifdef DEBUG_ESD
      const Real w2sq = 0.;
#else
      const Real w2sq = w2*w2;
#endif
      const int i2 = (i2s + halfNumSubSample)/a_numSubSample;  // Wave index
      for (int i1s = 0; i1s != n1s; ++i1s)
        {
          const Real w1 = (i1s + 0.5)*subSampleLength;
          const Real w1sq = w1*w1;
          const int i1 = (i1s + halfNumSubSample)/a_numSubSample;
          for (int i0s = 0; i0s != n0s; ++i0s)
            {
              const Real w0 = (i0s + 0.5)*subSampleLength;
              // radius to sub-cell center
              const Real r = std::sqrt(w0*w0 + w1sq + w2sq);
              // (r < 0.5 : idxBin = 0, 0.5 <= r < 1.5 : idxBin = 1, ...)
              const int idxBin = std::min((int)(r + 0.5), c_numIJKmin);
              // Wave index from sub-cell index
              const int i0 = (i0s + halfNumSubSample)/a_numSubSample;
              m_ESD[idxBin] += lc_data3d[i2][i1][i0][0][0]*subSampleVolume8;
            }
        }
    }
#ifdef DEBUG_ESD
  {
    Real tke = 0.;
    for (auto& esd : m_ESD) tke += esd;
    std::cout << "tke3: " << tke << std::endl;
  }
#endif
}

/*--------------------------------------------------------------------*/
//  Write the ESD to a CSV file
/**
 *//*-----------------------------------------------------------------*/

void
EnergySpectralDensity::write() const
{
  if (!m_use) return;
  std::ofstream fout(m_fileName);
  fout.setf(std::ios_base::scientific, std::ios_base::floatfield);
  fout.precision(std::numeric_limits<Real>::digits10);
  const int numSample = m_ESD.size() - 1;
  fout << "# Extra higher frequencies: " << m_ESD[numSample] << std::endl;
  fout << "# Samples: " << numSample << std::endl;
  for (int idx = 0, idx_end = m_ESD.size() - 1; idx != idx_end; ++idx)
    {
      fout << idx << ", " << m_ESD[idx] << std::endl;
    }
}

#endif  /* ! (defined(CH_USE_FFTW) && CH_SPACEDIM==3) */
