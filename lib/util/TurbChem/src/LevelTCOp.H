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
 * \file LevelTCOp.H
 *
 * \brief Level operator for post processing turbulent chemistry simulations
 *
 *//*+*************************************************************************/

#ifndef _LEVELTCOP_H_
#define _LEVELTCOP_H_

//----- Standard Library -----//

//----- Chombo Library -----//

//----- Internal -----//

#include "PatchOps.H"

//--Forward declarations

class EnergySpectralDensity;


/*******************************************************************************
 */
///  Level operator for post processing turbulent chemistry simulations
/**
 *//*+*************************************************************************/

class LevelTCOp
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
  LevelTCOp(const int a_level);

  /// Copy constructor not permitted
  LevelTCOp(const LevelTCOp&) = delete;

  /// Assignment constructor not permitted
  LevelTCOp& operator=(const LevelTCOp&) = delete;

  /// Destructor -- default destructor

/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Define the level operator (weak construction)
  void define(HDF5Handle&            a_handle,
              const std::string      a_fileName,
              const int              a_numGhostCells,
              const int              a_numOutGhostCells,
              const int              a_numOutComps,
              const std::vector<int> a_periodicity,
              Real&                  a_timeBaseLevel);

  /// Compute the quantities of interest
  void postProcessCompute(const std::vector<std::string> a_outTimePntVars,
                          const std::vector<std::string> a_outTimeSumVars,
                          const std::vector<std::string> a_outGlobalVars,
                          std::vector<Real>&             a_domSums,
                          const int                      a_deconvOrder,
                          const int                      a_compOrder,
                          const int                      a_numReqGhost,
                          const int                      a_numComps,
                          const int                      a_numOutComps,
                          EnergySpectralDensity&         a_ESD);

  /// Write global data to HDF5 file
  int writeDataHDF5(HDF5Handle&                    a_handle,
                    const std::vector<std::string> a_outGlobalVars,
                    const int                      a_numOutComps,
                    const int                      a_numOutGhost);

/*==============================================================================
 * Protected member functions
 *============================================================================*/

protected:

  /// Any protected functions

/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  LevelData<FArrayBox> m_dataOld;     ///< Existing data read from HDF5 file
  LevelData<FArrayBox> m_dataNew;     ///< Post-processed data
  PatchOps             m_patchOp;     ///< Operators for patches of level
  ProblemDomain        m_probDomain;  ///< Problem domain
  DisjointBoxLayout    m_dbl;         ///< Layout of boxes for level
  Box                  m_domainBox;   ///< Box covering entire domain
  RealVect             m_dx;          ///< Cell size on current level
  IntVect              m_refRatio;    ///< AMR refinement ratio on this level
  Real                 m_time;        ///< Simulation time for the level
  Real                 m_dt;          ///< Time step size on current level
  int                  m_level;       ///< Level number (starting at 0)
};

#endif  /* ! defined _LEVELTCOP_H_ */
