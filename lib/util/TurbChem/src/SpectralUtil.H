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

#if defined(CH_USE_FFTW) && CH_SPACEDIM==3
#include "fftw3.h"
#endif


/*******************************************************************************
 */
/// A factory for LevelData objects that have a layout optimized for fftw
/**
 *//*+*************************************************************************/

#if defined(CH_USE_FFTW) && CH_SPACEDIM==3
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
/// Class that computes energy spectral density using fftw
/**
 *//*+*************************************************************************/

class EnergySpectralDensity
{

/*==============================================================================
 * Constructors and destructors
 *============================================================================*/

  public:

  /// Construct with undefined state
  EnergySpectralDensity()
    :
    m_use(false),
    m_problemDomain(),
    m_r_numI(-1)
    { }

#if defined(CH_USE_FFTW) && CH_SPACEDIM==3

  /// Destructor
  ~EnergySpectralDensity();

  /// Full construction
  EnergySpectralDensity(const ProblemDomain& a_problemDomain)
    { defineSerial(a_problemDomain); }

  /// Copy constructor and copy assignment are deleted
  EnergySpectralDensity(const EnergySpectralDensity&) = delete;
  EnergySpectralDensity& operator=(const EnergySpectralDensity&) = delete;
  /// Move constructor and move assignment are default
  //  Beware that underlying objects may not yet have proper move semantics.
  //  Test extensively if you use these.
  EnergySpectralDensity(EnergySpectralDensity&&) = default;
  EnergySpectralDensity& operator=(EnergySpectralDensity&&) = default;
#endif

/*==============================================================================
 * Member functions
 *============================================================================*/

  /// Set to compute ESD
  void use()
    { m_use = true; }
  void unuse()
    { m_use = false; }
  bool isUsed() const
    { return m_use; }

  /// Set output file name
  void outputFileName(const std::string& a_fileName)
    { m_fileName = a_fileName; }

#if defined(CH_USE_FFTW) && CH_SPACEDIM==3

  /// Weak construction
  void defineSerial(const ProblemDomain& a_problemDomain);

  /// Copy data into local storage
  void load(const Box&       a_box,
            const FArrayBox& a_dataFab,
            const Interval&  a_intv);

  /// Apply the transform and compute the ESD
  void compute(const RealVect&    a_physDomain,
               const int          a_numSubSample = 4);

  /// Write the ESD to a CSV file
  void write() const;

#else  /* Null operation if FFTW not available */
  void defineSerial(const ProblemDomain& a_problemDomain) 
    { }
  void load(const Box&       a_box,
            const FArrayBox& a_dataFab,
            const Interval&  a_intv)
    { }
  void compute(const RealVect&    a_physDomain,
               const int          a_numSubSample = 4)
    { }
  void write() const
    { }
#endif

  const std::vector<Real>& getESDvec() const
    { return m_ESD; }

/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  bool m_use;                         ///< T - Compute ESD
  std::string m_fileName;             ///< Output csv file name
  
  ProblemDomain m_problemDomain;      ///< Last problem domain that was defined
                                      ///< If new define has same dimensions,
                                      ///< can avoid re-defining

  // Prefix c means applied to complex variables only and r for real only.
  // If c and r are both absent, it means the value it is the same for both.
  int m_r_numI;                       ///< Number of real variables in unit
                                      ///< stride direction
  int m_c_numI;                       ///< Number of complex variables in unit
                                      ///< stride direction
  int m_numJ;                         ///< Number of real/complex variables in
                                      ///< next to least-contiguous direction
  int m_numK;                         ///< Number of real/complex variables in
                                      ///< least-contiguous direction

  std::vector<Real> m_ESD;            ///< The computed energy spectral density

#if defined(CH_USE_FFTW) && CH_SPACEDIM==3

  FFTWDataFactory<FArrayBox> m_fftwDataFactory;
                                      ///< Factory using fftw_malloc
  LevelData<FArrayBox> m_fftwData;    ///< Layouts that match fftw's
                                      ///< requirements
  fftw_plan m_fftwForwardPlan;        ///< Plan do to the FFTs (forward)

#endif  /* ! (defined(CH_USE_FFTW) && CH_SPACEDIM==3) */
};

#endif  /* ! defined _SPECTRALUTIL_H_ */
