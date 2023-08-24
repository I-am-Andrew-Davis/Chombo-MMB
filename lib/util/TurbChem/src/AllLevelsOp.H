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
 * \file AllLevelsOp.H
 *
 * \brief Structure for defining and operating on all AMR levels
 *
 *//*+*************************************************************************/

#ifndef _ALLLEVELSOP_H_
#define _ALLLEVELSOP_H_

//----- Standard Library -----//

//----- Chombo Library -----//

//----- Internal -----//

#include "LevelTCOp.H"

//--Forward declarations

class EnergySpectralDensity;


/*******************************************************************************
 */
///  Structure for defining and operating on all AMR levels
/**
 *//*+*************************************************************************/

class AllLevelsOp
{

/*====================================================================*
 * Types
 *====================================================================*/

public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Constructor
  AllLevelsOp(
    const std::vector<std::string> a_outTimePntVars,
    const std::vector<std::string> a_outTimeSumVars,
    const std::vector<std::string> a_outGlobalVars,
    const std::vector<std::string> a_BCs,
    const std::vector<Real>        a_varScale,
    const Real                     a_timeMult,
    const int                      a_computeMethod,
    const int                      a_numOutComps,
    const int                      a_deconvOrder,
    const int                      a_compOrder,
    const int                      a_numReqGhost,
    const int                      a_numOutGhost,
    const int                      a_verbosity,
    const bool                     a_compress,
    const bool                     a_thermPerf,
    const bool                     a_mapped);

  /// Copy constructor not permitted
  AllLevelsOp(const AllLevelsOp&) = delete;

  /// Assignment constructor not permitted
  AllLevelsOp& operator=(const AllLevelsOp&) = delete;

  /// Destructor -- default destructor

/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Definition of the level hierarchy structure
  void define(const std::string&      a_fileName,
              const std::vector<int>& a_periodicity,
              Real&                   a_timeBaseLevel);

  /// Run through all levels present in the hierarchy
  void runLevels(std::vector<Real>&     a_domSums,
                 EnergySpectralDensity& a_ESD);

  /// Write the resulting data
  void writeLevels(const std::string a_fileName);
  

/*==============================================================================
 * Protected member functions
 *============================================================================*/

protected:

  /// Any protected functions

/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  std::vector<std::unique_ptr<LevelTCOp>> m_levelObj;
                                      ///< Vector of Levels
  std::vector<std::string> m_outTimePntVars;
                                      ///< List of variables to compute
                                      ///< and print at 1 space pnt over time
  std::vector<std::string> m_outTimeSumVars;
                                      ///< List of variables to compute, sum
                                      ///< over domain, and print over time
  std::vector<std::string> m_outGlobalVars;
                                      ///< List of variables to compute and
                                      ///< plot to HDF5 file (entire domain)
  std::vector<Real> m_varScale;       ///< Variable scaling (multiplication)
  Real              m_timeMult;       ///< Time scaling factor (multiplication)
  int               m_numLevels;      ///< Number of levels (starts at 1)
  int               m_maxLevel;       ///< Maximum level (starts at 0)
  int               m_compMethod;     ///< Method of computing the variables
                                      ///< (1) is standard pnt based method
  int               m_numComps;       ///< Number of components in input file
  int               m_numOutComps;    ///< Number of output components
  int               m_deconvOrder;    ///< Order of accuracy of deconvolution
  int               m_compOrder;      ///< Order of accuracy of computation
  int               m_numReqGhost;    ///< Number of ghost cells required for
                                      ///< computations in domain
  int               m_numOutGhost;    ///< Number of output ghost cells
  int               m_iteration;      ///< Iteration number for current file
  int               m_verbosity;      ///< Level of output during processing
  bool              m_compress;       ///< Flow is compressible
  bool              m_thermPerf;      ///< Fluid is thermally perfect
  bool              m_mapped;         ///< Simulation uses mapping
};

#endif  /* ! defined _ALLLEVELSOP_H_ */
