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
 * \file AllLevelsOp.cpp
 *
 * \brief Member functions for AllLevelsOp
 *
 *//*+*************************************************************************/

//----- Standard Library -----//

//----- Chombo Library -----//

#include "ParmParse.H"
#include "AMR.H"
#include "AMRLevel.H"
#include "CONSTANTS.H"
#include "AMRIO.H"
#include "CH_HDF5.H"

//----- Internal -----//

#include "AllLevelsOp.H"
#include "SpectralUtil.H"


/*******************************************************************************
 *
 * Class AllLevelsOp: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** \param[in]  a_something 
 *                      Something
 *//*-----------------------------------------------------------------*/

AllLevelsOp::AllLevelsOp(
  //**FIXME Why are the vectors not references?
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
  const bool                     a_mapped)
:
  m_numLevels(0),
  m_maxLevel(1),
  m_numComps(0),
  m_iteration(0)
{
  //**FIXME Why are these assignments not in the initializer list?
  m_outTimePntVars = a_outTimePntVars;
  m_outTimeSumVars = a_outTimeSumVars;
  m_outGlobalVars  = a_outGlobalVars;
  m_varScale       = a_varScale;
  m_timeMult       = a_timeMult;
  m_compMethod     = a_computeMethod;
  m_numOutComps    = a_numOutComps;
  m_deconvOrder    = a_deconvOrder;
  m_compOrder      = a_compOrder;
  m_numReqGhost    = a_numReqGhost;
  m_numOutGhost    = a_numOutGhost;
  m_verbosity      = a_verbosity;
  m_compress       = a_compress;
  m_thermPerf      = a_thermPerf;
  m_mapped         = a_mapped;
}


/*==============================================================================
 * Public member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Define the level operator (weak construction)
/** \param[in]  a_level The level index
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *//*-----------------------------------------------------------------*/
void
AllLevelsOp::define(const std::string&      a_fileName,
                    const std::vector<int>& a_periodicity,
                    Real&                   a_timeBaseLevel)
{
  // Read in the file
  HDF5Handle handle(a_fileName, HDF5Handle::OPEN_RDONLY);
  HDF5HeaderData header;
  header.readFromFile(handle);
  // Read the main header from the file
  m_numComps  = header.m_int["num_components"];
  m_numLevels = header.m_int["num_levels"];
  m_maxLevel  = header.m_int["max_level"];
  m_iteration = header.m_int["iteration"];
  // Read the level headers
  for (int lvl = 0; lvl != m_numLevels; ++lvl)
    {
      std::unique_ptr<LevelTCOp> lvlPtr(new LevelTCOp(lvl));
      lvlPtr->define(handle,
                     a_fileName,
                     m_numReqGhost,
                     m_numOutGhost,
                     m_numOutComps,
                     a_periodicity,
                     a_timeBaseLevel);
      m_levelObj.push_back(std::move(lvlPtr));
    }
  handle.close();
}

/*--------------------------------------------------------------------*/
//  Define the level operator (weak construction)
/** \param[in]  a_level The level index
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *//*-----------------------------------------------------------------*/
void
AllLevelsOp::runLevels(std::vector<Real>&     a_domSums,
                       EnergySpectralDensity& a_ESD)
{
  if (a_ESD.isUsed() && m_numLevels != 1)
    {
      std::cout << "Computing energy spectral density requires 1 AMR level\n";
      std::cout << "ABORTING\n";
      exit(1);
    }
  for (int lvl = 0; lvl != m_numLevels; ++lvl)
    {
      m_levelObj[lvl]->postProcessCompute(m_outTimePntVars,
                                          m_outTimeSumVars,
                                          m_outGlobalVars,
                                          a_domSums,
                                          m_deconvOrder,
                                          m_compOrder,
                                          m_numReqGhost,
                                          m_numComps,
                                          m_numOutComps,
                                          a_ESD);
    }
}

/*--------------------------------------------------------------------*/
//  Define the level operator (weak construction)
/** \param[in]  a_level The level index
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *//*-----------------------------------------------------------------*/
void
AllLevelsOp::writeLevels(const std::string a_fileName)
{
  // FIXME: handle should be obtained as a function argument a_handle
  //        this will allow for writing of multiple levels of data
  HDF5Handle handle(a_fileName,HDF5Handle::CREATE);
  HDF5HeaderData header;
  // FIXME: option should allow plotting all data (old and new)
  //        together in single file
  header.m_int["max_level"]  = m_maxLevel;
  header.m_int["num_levels"] = m_numLevels;
  header.m_int["iteration"]  = m_iteration;
  header.writeToFile(handle);

  for (int lvl = 0; lvl != m_numLevels; ++lvl)
    {
      m_levelObj[lvl]->writeDataHDF5(handle,
                                     m_outGlobalVars,
                                     m_numOutComps,
                                     m_numOutGhost);
    }
  handle.close();
}
