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
 *  \file ThreadTopology.cpp
 *
 *  \brief Non-inline definitions for classes in ThreadTopology.H
 *
 *//*+*************************************************************************/


//----- Standard Library -----//

//----- Internal -----//

#include "ThreadTopology.H"
#include "ThreadParameters.H"

#include "NamespaceHeader.H"

using namespace ThreadTools;

//-- define static members
bool ThreadTopology::s_hasRTPrio = false;

/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Construct the thread topology
/** \param[in]  a_numWorkers
 *                      The number of worker threads
 *  \param[in]  a_numProcComms
 *                      The number of process communicator threads
 *  \param[in]  a_numGpuComms
 *                      The number of GPU communicator threads
 *//*-----------------------------------------------------------------*/

ThreadTopology::ThreadTopology(const unsigned a_numWorkers,
                               const unsigned a_numProcComms,
                               const unsigned a_numGpuComms)
: m_numWorkers {a_numWorkers},
  m_numProcComms {a_numProcComms},
  m_numGpuComms {a_numGpuComms}
{
  ThreadSupport::setRealtimePriorities();
}

/*==============================================================================
 * Public member functions
 *============================================================================*/

#include "NamespaceFooter.H"
