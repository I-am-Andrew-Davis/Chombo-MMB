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
 * \file SpectralUtil.cpp
 *
 * \brief Spectral operations for use by Chombo applications
 *
 * Makes use of fftw for FFTs.
 *
 *//*+*************************************************************************/

#ifndef _SPECTRALUTIL_H_
#define _SPECTRALUTIL_H_

//----- Chombo -----//

#include "FArrayBox.H"
#include "LevelData.H"
#include "BaseFabMacros.H"

//----- External -----//

#if defined(CH_USE_FFTW) && defined(CH_MPI) && CH_SPACEDIM==3
#include "fftw3-mpi.h"

/*****************************************************************************/
/**
 * \brief Operator that switches from Struct-of-Arrays (SoA) to
 *        Array-of-Structs (AoS)
 *
 *  Chombo normally features components at largest stride in an array.  This
 *  switches to components having unit stride.  The resulting memory is used
 *  by FFTW.
 *
 *//*+************************************************************************/
template <typename T>
class LDOperatorAoSIn;
template <>
class LDOperatorAoSIn<FArrayBox> : public LDOperator<FArrayBox>
{
public:

  /// Destructor
  virtual ~LDOperatorAoSIn();

  /// Linear in switches from SoA to AoS
  virtual void linearIn(FArrayBox&      a_arg,
                        void*           a_buf,
                        const Box&      a_R,
                        const Interval& a_comps) const;

  /// Local copy switches from SoA to AoS
  virtual void op(FArrayBox&       a_dst,
                  const Box&       a_regionFrom,
                  const Interval&  a_Cdst,
                  const Box&       a_regionTo,
                  const FArrayBox& a_src,
                  const Interval&  a_Csrc) const;
};


/*****************************************************************************/
/**
 * \brief Operator that switches from Array-of-Structs (AoS) to
 *        Struct-of-Arrays (SoA)
 *
 *  Chombo normally features components at largest stride in an array.
 *  This switches from components having unit stride to a normal Chombo layout.
 *  The AoS memory layout is used by FFTW.
 *
 *//*+************************************************************************/
template <typename T>
class LDOperatorAoSOut;
template <>
class LDOperatorAoSOut<FArrayBox> : public LDOperator<FArrayBox>
{
public:

  /// Destructor
  virtual ~LDOperatorAoSOut();

  /// Linear out switches from AoS to SoA
  virtual void linearOut(const FArrayBox& a_arg,
                         void*            a_buf,
                         const Box&       a_R,
                         const Interval&  a_comps) const;

  /// Local copy switches from AoS to SoA
  virtual void op(FArrayBox&       a_dst,
                  const Box&       a_regionFrom,
                  const Interval&  a_Cdst,
                  const Box&       a_regionTo,
                  const FArrayBox& a_src,
                  const Interval&  a_Csrc) const;
};


/*******************************************************************************
 */
/// A factory for LevelData objects that have a layout optimized for fftw
/**
 *//*+*************************************************************************/
template <typename T>
class FFTWDataFactory;
template <>
class FFTWDataFactory<FArrayBox> : public DataFactory<FArrayBox>
{
public:

  /// Constructor
  FFTWDataFactory()
    :
    m_numCreate(-1)
    { }

  /// Destructor
  virtual ~FFTWDataFactory();

  /// Weak construction
  void define(const int a_l_alloc);

  /// Create an FArrayBox for the LevelData
  virtual FArrayBox* create(const Box&       a_box,
                            int              a_ncomps,
                            const DataIndex& a_datInd) const;

  /// Grab the date pointer (used when no box on this processor)
  Real* dataPtr() const;

//--Data

protected:

  Real* m_data;
  mutable int m_numCreate;
};

#endif


/*******************************************************************************
 */
/// Class that implements spectral filtering using fftw
/**
 *//*+*************************************************************************/

class SpectralFilter
{

/*--------------------------------------------------------------------*/
/// Type of spectral filter profiles
/*--------------------------------------------------------------------*/

public:

enum SpectralFilterProfile
{
  SpectralFilterSharp = 0,            ///< Sharp spectral cutoff (rectangular)
  // SharpIsotropic is the default
  SpectralFilterSharpIsotropic,       ///< Sharp spectral cutoff (spherical)
  SpectralFilterTanh,                 ///< Hyperbolic-tangent filter
  SpectralFilterGaussian              ///< Gaussian filter
};

/*==============================================================================
 * Constructors and destructors
 *============================================================================*/

  public:

  /// Construct with undefined state
  SpectralFilter()
    :
    m_filterParam(std::numeric_limits<Real>::infinity())
    { }

#if defined(CH_USE_FFTW) && defined(CH_MPI) && CH_SPACEDIM==3

  /// Full construction
  SpectralFilter(
    const ProblemDomain&        a_problemDomain,
    const DisjointBoxLayout&    a_boxes,
    const int                   a_numComp,
    const SpectralFilterProfile a_filterProfile = SpectralFilterSharpIsotropic)
    :
    m_filterParam(std::numeric_limits<Real>::infinity())
    { define(a_problemDomain, a_boxes, a_numComp, a_filterProfile); }

  /// Destructor
  ~SpectralFilter() = default;

  /// Copy constructor and copy assignment are deleted
  SpectralFilter(const SpectralFilter&) = delete;
  SpectralFilter& operator=(const SpectralFilter&) = delete;
  /// Move constructor and move assignment are default
  //  Beware that underlying objects may not yet have proper move semantics.
  //  Test extensively if you use these.
  SpectralFilter(SpectralFilter&&) = default;
  SpectralFilter& operator=(SpectralFilter&&) = default;

/*==============================================================================
 * Member functions
 *============================================================================*/

  /// Weak construction
  void define(
    const ProblemDomain&        a_problemDomain,
    const DisjointBoxLayout&    a_boxes,
    const int                   a_numComp,
    const SpectralFilterProfile a_filterProfile = SpectralFilterSharpIsotropic);

  /// Parameter if tanhFilter or Gaussian filter selected
  /** a_filterParam       Not required for either sharp cutoff filter
   *                      tanh     - smoothing width
   *                      Gaussian - rolloff width
   */
  void setFilterParam(const Real a_filterParam)
    { m_filterParam = a_filterParam; }

  /// Apply the filter
  void apply(const IntVect&        a_gridResolution,
             LevelData<FArrayBox>& a_U,
             const Interval&       a_intv);
#else  /* Null operation if FFTW not available */
  void apply(const IntVect&        a_gridResolution,
             LevelData<FArrayBox>& a_U,
             const Interval&       a_intv)
    { }
#endif

/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  Real m_filterParam;                 ///< Not required for either sharp cutoff
                                      ///< filter
                                      ///< tanh     - smoothing width
                                      ///< Gaussian - rolloff width

#if defined(CH_USE_FFTW) && defined(CH_MPI) && CH_SPACEDIM==3
public:

  static constexpr int c_verbosity = 1;

protected:

  int m_numComp;                      ///< Number of components expected in the
                                      ///< state to filter.  The FFTW plans
                                      ///< are built using this many components.
  SpectralFilterProfile m_filterProfile;
                                      ///< Shape of the filter

  // Prefix c means applied to complex variables only and r for real only.
  // If c and r are both absent, it means the value it is the same for both.
  // Prefix l means local to this processor
  // Prefix t means this value applies to transposed form
  ptrdiff_t m_c_numI;                 ///< Number of complex variables in unit
                                      ///< stride direction
  ptrdiff_t m_numJ;                   ///< Number of real/complex variables in
                                      ///< next to least-contiguous direction
  ptrdiff_t m_numK;                   ///< Number of real/complex variables in
                                      ///< least-contiguous direction
  ptrdiff_t m_l_numK;                 ///< Number of local real/complex
                                      ///< variables in least-contiguous
                                      ///< direction on this processor
  ptrdiff_t m_l_idxKbeg;              ///< Begin 'K' index on this processor for
                                      ///< local partition of array
  ptrdiff_t m_lt_numJ;                ///< Number of local real/complex
                                      ///< variables in next to least-contiguous
                                      ///< direction on this processor when in
                                      ///< transposed form
  ptrdiff_t m_lt_idxJbeg;             ///< Begin 'J' index on this processor
                                      ///< when in transposed form for local
                                      ///< partition of array

  FFTWDataFactory<FArrayBox> m_fftwFilterDataFactory;
                                      ///< Factory using fftw_malloc
  LevelData<FArrayBox> m_fftwFilterData;
                                      ///< Layouts that match fftw's
                                      ///< requirements
  fftw_plan m_fftwFilterForwardPlan;  ///< Plans do to the FFTs (forward)
  fftw_plan m_fftwFilterBackwardPlan; ///< Plans do to the FFTs (reverse)
  Copier m_fftwFilterInCopier;        ///< Copy from m_boxes to fftw layout
  Copier m_fftwFilterOutCopier;       ///< Copy from fftw layout to m_boxes

#endif  /* ! (defined(CH_USE_FFTW) && defined(CH_MPI) && CH_SPACEDIM==3) */
};

#endif  /* ! defined _SPECTRALUTIL_H_ */
