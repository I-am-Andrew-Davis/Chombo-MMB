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
 *  \file ThreadTeamArchitect.cpp
 *
 *  \brief Non-inline definitions for classes in ThreadTeamArchitect.H
 *
 *//*+*************************************************************************/

//----- Internal -----//

#include "ThreadTeamArchitect.H"
#include "ThreadTopology.H"
#include "ThreadTopologyBuilder.H"
#include "ThreadParameters.H"

//----- System -----//

#include <sys/resource.h>
#include <sys/types.h>
#include <sys/syscall.h>

//----- Standard Library -----//

#include <sstream>

#include "NamespaceHeader.H"

using namespace ThreadTools;

/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** \param[in]  a_hwTopology
 *//*-----------------------------------------------------------------*/

ThreadTeamArchitect::ThreadTeamArchitect(
  std::shared_ptr<const HardwareTopology> a_hwTopology)
:
ThreadArchitect(a_hwTopology)
{
}

/*==============================================================================
 * Public member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Setup the thread topology
//  The number of worker threads can not be 0. If less than zero, the
//  number of worker threads will equal the number of cores this
//  process is bound to
/** \param[in]  a_numWorkers
 *                      The number of worker threads
 *  \param[in]  a_threadsPerTeam
 *                      The number of threads per team
 *  \param[in]  a_numProcComms
 *                      The number of process communication threads to use
 *  \param[in]  a_numGpuComms
 *                      The number of GPU communication threads to use
 *  \param[in]  a_extraStackSize_MiB
 *                      Extra MiB on top of the default stack size
 *  \param[in]  a_useHyperThreads
 *                      Whether to use hyperthreads
 *  \param[in]  a_bindThread
 *                      Whether to bind threads to hardware
 *//*-----------------------------------------------------------------*/

void ThreadTeamArchitect::setupTopology(const int      a_numWorkers,
                                        const int      a_threadsPerTeam,
                                        const int      a_numProcComms,
                                        const int      a_numGpuComms,
                                        const size_t   a_extraStackSize_MiB,
                                        const bool     a_useHyperThreads,
                                        const bool     a_bindThread)
{
  // Set up the thread topology using the base class
  ThreadArchitect::setupTopology(a_numWorkers,
                                 a_numProcComms,
                                 a_numGpuComms,
                                 a_extraStackSize_MiB,
                                 a_useHyperThreads,
                                 a_bindThread);

  const auto numWorkers = m_threadTopology->numWorkers();

  TH_assert(a_threadsPerTeam != 0
            && "Thread team size can not be 0. If less than zero,\n"
               " the thread team size will be equal to the number\n"
               " of threads per core.");

  setThreadsPerTeam(a_threadsPerTeam, a_useHyperThreads);

  const int numTeams = numWorkers / m_threadsPerTeam;
  TH_assert(m_threadsPerTeam * numTeams == static_cast<int>(numWorkers));

  TH_assert(numWorkers % m_threadsPerTeam == 0
            && "The number of workers must be divisible by the"
               " threads per team");

  TH_assert(m_threadsPerTeam != 0
            && "Thread team size can not be 0. If less than zero,\n"
               " the thread team size will be equal to the number\n"
               " of threads per core.");

  if (s_verbosity >= Verb::v1)
    {
      std::stringstream ss;
      ss << "Num Teams                             : " << numTeams << "\n"
         << "Threads per team                      : " << m_threadsPerTeam << "\n\n";
      s_thrOut << ss.rdbuf();
    }

  m_teams.reserve(numTeams + 1);
  m_teams.emplace_back(1); // master
  for (int team = 0; team < numTeams; ++team)
    {
      m_teams.emplace_back(m_threadsPerTeam);
    }
  // TODO teams should have variable sizes
}

/*--------------------------------------------------------------------*/
//  Set the number of threads per team
/** \param[in]  a_threadsPerTeam
 *                      The number of threads per team. If < 0, the
 *                      number of threads per team equals 1 if
 *                      hyperthreads not present and 2 if hyper threads
 *                      are present.
 *  
 *//*-----------------------------------------------------------------*/

void ThreadTeamArchitect::setThreadsPerTeam(const int  a_threadsPerTeam,
                                            const bool a_useHyperThreads)
{
  if (a_useHyperThreads && m_hwTopology->hasHyperThreads())
    {
      setThreadsPerTeamWithHyperThreads(a_threadsPerTeam);
    }
  else
    {
      setThreadsPerTeamNoHyperThreads(a_threadsPerTeam);
    }
}

/*--------------------------------------------------------------------*/
//  Set the number of threads per team when using hyperthreads
/** \param[in]  a_threadsPerTeam
 *                      The number of threads per team. If < 0, the
 *                      number of threads per team equals 1 if
 *                      hyperthreads not present and 2 if hyper threads
 *                      are present.
 *//*-----------------------------------------------------------------*/

void ThreadTeamArchitect::setThreadsPerTeamWithHyperThreads(
  const int  a_threadsPerTeam)
{
  m_useHyperThreads = true;

  if (a_threadsPerTeam < 0)
    {
      m_threadsPerTeam = m_hwTopology->hasHyperThreads() ? 2 : 1;
    }
  else
    {
      m_threadsPerTeam = a_threadsPerTeam;
    }
}

/*--------------------------------------------------------------------*/
//  Set the number of threads per team when not using hyperthreads
/** \param[in]  a_threadsPerTeam
 *                      The number of threads per team. If < 0, the
 *                      number of threads per team equals 1.
 *//*-----------------------------------------------------------------*/

void ThreadTeamArchitect::setThreadsPerTeamNoHyperThreads(
  const int a_threadsPerTeam)
{
  m_useHyperThreads = false;

  if (a_threadsPerTeam > 0)
    {
      m_threadsPerTeam = a_threadsPerTeam;
    }
  else
    {
      m_threadsPerTeam = 1;
    }
}

#include "NamespaceFooter.H"
