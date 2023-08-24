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
 *  \file ThreadArchitect.cpp
 *
 *  \brief Non-inline definitions for classes in ThreadArchitect.H
 *
 *//*+*************************************************************************/


//----- Internal -----//

#include "ThreadArchitect.H"
#include "ThreadParameters.H"
#include "ThreadTopology.H"
#include "ThreadTopologyBuilder.H"

//----- System -----//

#include <sys/resource.h>
#include <sys/syscall.h>
#include <sys/types.h>

//----- Standard Library -----//

#include <sstream>

#include "NamespaceHeader.H"

using namespace ThreadTools;

std::chrono::time_point<std::chrono::steady_clock> ThreadArchitect::s_startTime;

/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** \param[in]  a_hwTopology
 *                      Hardware topology
 *//*-----------------------------------------------------------------*/

ThreadArchitect::ThreadArchitect(
  std::shared_ptr<const HardwareTopology> a_hwTopology)
:
m_hwTopology{a_hwTopology ?
               a_hwTopology : std::make_shared<const HardwareTopology>()},
m_threadBuilder{m_hwTopology}
{
  s_startTime = std::chrono::steady_clock::now();
  setupMaster();
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

ThreadArchitect::~ThreadArchitect()
{
  waitToFinish();
}

/*==============================================================================
 * Public member functions
 *============================================================================*/


/*--------------------------------------------------------------------*/
//  Setup the master thread
/**
 *//*-----------------------------------------------------------------*/

void ThreadArchitect::setupMaster()
{
  if (s_verbosity >= Verb::v1)
    {
      s_thrOut << "Master\n";
      s_thrOut << "hasRT: " << ThreadSupport::hasRTPrio() << std::endl;
    }
  int a_policy;
  int a_prio;
  getThreadPrio(Role::master, a_policy, a_prio);
  ThreadSupport::setThreadPrio(a_policy, a_prio);
  ThreadSupport::printThreadCPUBinding(m_hwTopology->topology());
  ThreadSupport::printPThreadAttr(pthread_self());
  pid_t tid = syscall(SYS_gettid);
  errno = 0;
  std::stringstream ss;
  int niceness = getpriority(PRIO_PROCESS, tid);
  if (errno)
    {
      s_thrErr << "getpriority(PRIO_PROCESS, tid)" << std::endl;
      niceness = 0;
    }
  if (s_verbosity >= Verb::v1)
    {
      ss << "  Niceness                            : " << niceness << std::endl;
    }
  s_thrOut << ss.rdbuf();
  // Do not bind the master thread. It is the responsibility of the
  // user to bind processes and the master thread will keep the same cpuset
  // as the process. The communication and gpu thread will also use the
  // cpuset from the process.
}

/*--------------------------------------------------------------------*/
//  Setup the thread topology
//  The number of worker threads can not be 0. If less than zero, the
//  number of worker threads will equal the number of cores this
//  process is bound to
/** \param[in]  a_numWorkers
 *                      The number of worker threads
 *  \param[in]  a_numProcComms
 *                      The number of process communication threads to use
 *  \param[in]  a_numGpuComms
 *                      The number of GPU communication threads to use
 *  \param[in]  a_extraStackSize_MiB
 *                      Extra MiB on top of the default stack size
 *  \param[in]  a_affinity
 *
 *  \param[in]  a_useHyperThreads
 *                      Whether to use hyperthreads
 *  \param[in]  m_bindThreads
 *                      Should the thread be bound to a hardware obj
 *//*-----------------------------------------------------------------*/

void ThreadArchitect::setupTopology(const int      a_numWorkers,
                                    const int      a_numProcComms,
                                    const int      a_numGpuComms,
                                    const size_t   a_extraStackSize_MiB,
                                    const bool     a_useHyperThreads,
                                    const bool     a_bindThread)
{
  getThreadStackSize(a_extraStackSize_MiB);

  m_bindThread = a_bindThread;

  ThreadTopologyBuilder builder(m_hwTopology);
  m_threadTopology =
    std::make_shared<ThreadTopology>(builder.build(a_numWorkers,
                                                   a_useHyperThreads,
                                                   a_numProcComms,
                                                   a_numGpuComms));

  const auto numWorkers = m_threadTopology->numWorkers();

  assert(numWorkers != 0 &&
         "The number of worker threads can not be 0. If less\n"
         " than zero, the number of worker threads will equal\n"
         " the number of cores this process is bound to.");

  if (s_verbosity >= Verb::v1)
    {
      std::stringstream ss;
      ss << std::boolalpha
         << "Number of workers                     : " << numWorkers << "\n"
         << "Using hyper threads                   : " << m_threadTopology->usingHyperThreads() << "\n"
         << "Number of comm threads                : " << m_threadTopology->numProcComms() << "\n"
         << "Number of GPU threads                 : " << m_threadTopology->numGpuComms() << "\n"
         << "\n";
      s_thrOut << ss.rdbuf();
    }
}

/*--------------------------------------------------------------------*/
//  Get the thread stack size
/** \param[in] a_extraStackSize_MiB
 *                      Extra stack size to add to the default size
 *//*-----------------------------------------------------------------*/

void ThreadArchitect::getThreadStackSize(const size_t a_extraStackSize_MiB)
{
  struct rlimit limits;
  if (getrlimit(RLIMIT_STACK, &limits))
    {
      s_thrErr << "getrlimit(RLIMIT_STACK)" << std::endl;
    }
  if (limits.rlim_cur == RLIM_INFINITY)
    {
      if (s_verbosity >= Verb::v1)
        {
          s_thrOut << "Default stack size                    : "
                   << "infinite (using 8192 kiB for thread)" << std::endl;
        }
      limits.rlim_cur = 8192 * 1024;
    }
  else
    {
      if (s_verbosity >= Verb::v1)
        {
          s_thrOut << "Default stack size                    : "
                   << limits.rlim_cur / 1024 << " kiB" << std::endl;
        }
    }
  m_stackSize_B = m_hwTopology->alignToPage(limits.rlim_cur +
                                            a_extraStackSize_MiB * 1024 * 1024);
  m_guardSize_B = 1 * m_hwTopology->getPageSize(); // For guard
  assert(m_stackSize_B >= static_cast<size_t>(m_hwTopology->getPageSize()));
  if (s_verbosity >= Verb::v1)
    {
      s_thrOut << "Per worker thread stack size          : "
               << m_stackSize_B / 1024 << " kiB" << std::endl;
    }
}

/*--------------------------------------------------------------------*/
//  Set thread priority based on role
/** \param[in]  a_role  Role of the thread
 *  \return             0 : success
 *                      1 : error
 *                      2 : insufficient privileges
 *  This routine simply calls the more general version using
 *  hard-coded priorities for different thread roles.  Note two
 *  different approaches:
 *    1. If we can use realtime scheduling, the worker threads are
 *       not modified but special threads are set to a realtime
 *       scheduler (with typically low priority).  However, realtime
 *       threads have priority over all other threads.
 *    2. If we do not have privileges for realtime threads, all worker
 *       threads are set to a high nice level (giving them lower
 *       priority).
 *  It is preferred to striclty bind the worker threads but let the
 *  special threads float around the hardware to find an
 *  underutilized core.
 *//*-----------------------------------------------------------------*/
int ThreadArchitect::setThreadPrio(const Role a_role) const noexcept
{
  if (ThreadSupport::hasRTPrio())
    {
      switch (a_role)
        {
        case Role::master:
          return ThreadSupport::setThreadPrio(SCHED_FIFO, 2);
        case Role::gpu:
          return ThreadSupport::setThreadPrio(SCHED_FIFO, 1);
        case Role::comm:
          return ThreadSupport::setThreadPrio(SCHED_FIFO, 1);
        case Role::worker:
        default:
          return 0;
        }
    }
  else
    {
      switch (a_role)
        {
        case Role::worker:
          return ThreadSupport::setThreadPrio(SCHED_OTHER, 15);
        default:
          return 0;
        }
    }
  return 0;
}

#include "NamespaceFooter.H"
