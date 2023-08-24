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
 *  \file ThreadArchitect.H
 *
 *  \brief Definition of class ThreadArchitect
 *
 *//*+*************************************************************************/

#ifndef _THREADARCHITECT_H_
#define _THREADARCHITECT_H_

//----- Standard Library -----//

#include <thread>

//----- Internal -----//

#include "ThreadSupport.H"
#include "ThreadTopology.H"
#include "HardwareTopology.H"
#include "ThreadBuilder.H"

#include "NamespaceHeader.H"

namespace ThreadTools
{

/// Roles of threads
  enum class Role: int8_t
  {
    master,
    gpu,
    comm,
    worker
  };

/*******************************************************************************
 */
/// Builds and places threads
/**
 *  Builds worker, gpu comm, and proc comm threads and binds them to hardware
 *
 ******************************************************************************/

  class ThreadArchitect
  {


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

    /// Constructor
    ThreadArchitect(std::shared_ptr<const HardwareTopology> hwTopology = nullptr);

    /// Destructor
    virtual ~ThreadArchitect();

//--Use synthesized copy, copy assignment, move, and move assignment


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

    /// Setup the thread topology
    void setupTopology(int numThreads,
                       int numCommThreads,
                       int numGpuThreads,
                       size_t extraStackSize_MiB,
                       bool useHyperThreading,
                       bool bindThread = true);

    /// Build the next worker thread
    int buildNextWorker(void* (*func)(void*), void* args);

    /// Build the next GPU thread
    int buildNextGpuCommThread(void* (*func)(void*), void* args);

    /// Build the next proc communication thread
    int buildNextProcCommThread(void* (*func)(void*), void* args);

    /// Get the thread topology
    std::shared_ptr<const ThreadTopology> getThreadTopology() const;

    /// Get the hardware topology
    std::shared_ptr<const HardwareTopology> getHardwareTopology() const;

    /// Wait until all threads are finished
    virtual void waitToFinish();

    /// Are all the threads built?
    bool isFull() const;

  protected:

    /// Configures the calling thread as master
    void setupMaster();

    /// Get the thread stack size
    void getThreadStackSize(size_t a_extraStackSize_MiB);

    /// Set the next worker thread
    void setNextWorkerThread(std::thread&& thread);

    /// Set the next GPU communicator thread
    void setNextGpuCommThread(std::thread&& thread);

    /// Set the next process communicator thread
    void setNextProcCommThread(std::thread&& thread);

    /// Set thread priority based on role
    int setThreadPrio(Role a_role) const noexcept;

    /// Get the thread policy and priority for a Role
    void getThreadPrio(Role a_role, int& a_policy, int& a_prio) const noexcept;

/*==============================================================================
 * Data members
 *============================================================================*/

  public:

  static std::chrono::time_point<std::chrono::steady_clock> s_startTime;

  protected:
  
    // TODO use unique for ThreadTopology. probably not needed
    std::shared_ptr<ThreadTopology> m_threadTopology;
                                      ///<
    std::shared_ptr<const HardwareTopology> m_hwTopology;
                                      ///<

    ThreadBuilder m_threadBuilder;    ///<
    std::vector<std::thread> m_workers;
                                      ///<
    std::vector<std::thread> m_procComms;
                                      ///<
    std::vector<std::thread> m_gpuComms;
                                      ///<

    // std::vector of teams
    // map of internal tid to thread
    // map of hardware tid to thread

    bool m_useHyperThreads;           ///<

    static thread_local int tls_thrID;
                                      ///<

    size_t m_stackSize_B;             ///<
    size_t m_guardSize_B;             ///<
    bool m_bindThread;
  };

/*******************************************************************************
 *
 * Class ThreadArchitect: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Get the thread topology
/** \return             The thread topology
 *//*-----------------------------------------------------------------*/

  inline
  std::shared_ptr<const ThreadTopology>
  ThreadArchitect::getThreadTopology() const
  {
    return m_threadTopology;
  }

/*--------------------------------------------------------------------*/
//  Get the hardware topology
/** \return             The hardware topology
 *//*-----------------------------------------------------------------*/

  inline
  std::shared_ptr<const HardwareTopology>
  ThreadArchitect::getHardwareTopology() const
  {
    return m_hwTopology;
  }

/*--------------------------------------------------------------------*/
//  Build the next worker thread
/** \param[in]  a_func  The function the thread will execute
 *  \param[in]  a_args  The args for the function
 *  \return             The status of the system calls
 *//*-----------------------------------------------------------------*/

  inline
  int ThreadArchitect::buildNextWorker(void* (*func)(void*), void* args)
  {
    int status = 0;
    int policy;
    int prio;
    const int thrId = m_workers.size() + 1;
    getThreadPrio(Role::worker, policy, prio);

    std::thread thr;

    if (m_bindThread)
      {
        hwloc_obj_t pu = m_threadTopology->getObjOfWorker(m_workers.size());
        thr = m_threadBuilder.createThread(pu,
                                           m_stackSize_B,
                                           m_guardSize_B,
                                           policy,
                                           prio,
                                           thrId,
                                           func,
                                           args,
                                           status);
      }
    else
      {
        thr = m_threadBuilder.createThreadNoBind(policy,
                                                 prio,
                                                 thrId,
                                                 func,
                                                 args,
                                                 status);
      }

    setNextWorkerThread(std::move(thr));

    return status;
  }

/*--------------------------------------------------------------------*/
//  Build the next GPU communicator thread
/** \param[in]  a_func  The function the thread will execute
 *  \param[in]  a_args  The args for the function
 *  \return             The status of the system calls
 *//*-----------------------------------------------------------------*/

  inline
  int ThreadArchitect::buildNextGpuCommThread(void* (*a_func)(void*),
                                              void* a_args)
  {
    assert(false);
    int status = 0;
    const int thrId = 1 + m_threadTopology->numWorkers() + m_gpuComms.size();
    (void)thrId;
    return status;
  }

/*--------------------------------------------------------------------*/
//  Build the next process communicator thread
/** \param[in]  a_func  The function the thread will execute
 *  \param[in]  a_args  The args for the function
 *  \return             The status of the system calls
 *//*-----------------------------------------------------------------*/

  inline
  int ThreadArchitect::buildNextProcCommThread(void* (*a_func)(void*),
                                               void* a_args)
  {
    assert(false);
    int status = 0;
    const int thrId = 1 + m_threadTopology->numWorkers()
                        + m_threadTopology->numGpuComms()
                        + m_procComms.size();
    (void)thrId;
    return status;
  }

/*--------------------------------------------------------------------*/
//  Set the next worker thread
/** \param[in]  a_thread The thread
 *//*-----------------------------------------------------------------*/

  inline
  void ThreadArchitect::setNextWorkerThread(std::thread&& a_thread)
  {
    assert(m_threadTopology != nullptr);
    assert(m_threadTopology->numWorkers() > m_workers.size());
    m_workers.emplace_back(std::forward<std::thread>(a_thread));
  }

/*--------------------------------------------------------------------*/
//  Set the next GPU communicator thread
/** \param[in]  a_thread The thread
 *//*-----------------------------------------------------------------*/

  inline
  void ThreadArchitect::setNextGpuCommThread(std::thread&& a_thread)
  {
    assert(m_threadTopology != nullptr);
    assert(m_threadTopology->numGpuComms() > m_gpuComms.size());
    m_gpuComms.emplace_back(std::forward<std::thread>(a_thread));
  }

/*--------------------------------------------------------------------*/
//  Set the next process communicator thread
/** \param[in]  a_thread The thread
 *//*-----------------------------------------------------------------*/

  inline
  void ThreadArchitect::setNextProcCommThread(std::thread&& a_thread)
  {
    assert(m_threadTopology != nullptr);
    assert(m_threadTopology->numProcComms() > m_procComms.size());
    m_procComms.emplace_back(std::forward<std::thread>(a_thread));
  }


// TODO put this someplace else
/*--------------------------------------------------------------------*/
//  Get the thread policy and priority for a Role
/** \param[in]  a_role  The thread role
 *  \param[in]  a_policy
 *                      The scheduling polocy of the Role
 *  \param[in]  a_prio  The priority of the role
 *//*-----------------------------------------------------------------*/

  inline
  void
  ThreadArchitect::getThreadPrio(const Role a_role,
                                 int&       a_policy,
                                 int&       a_prio) const noexcept
  {
    a_policy = 0;
    a_prio = -1;
    if (ThreadSupport::hasRTPrio())
      {
        switch (a_role)
          {
          case Role::master:
            a_policy = SCHED_FIFO;
            a_prio = 2;
            break;
          case Role::gpu:
            a_policy = SCHED_FIFO;
            a_prio = 1;
            break;
          case Role::comm:
            a_policy = SCHED_FIFO;
            a_prio = 1;
            break;
          case Role::worker:
            a_policy = SCHED_FIFO;
            a_prio = 3;
          default:
            break;
          }
      }
    else
      {
        switch (a_role)
          {
          case Role::worker:
            a_policy = SCHED_OTHER;
            a_prio = 15;
            break;
          default:
            break;
          }
      }
    return;
  }

/*--------------------------------------------------------------------*/
//  Wait for all the threads to finish
/** 
 *//*-----------------------------------------------------------------*/

  inline
  void ThreadArchitect::waitToFinish()
  {
    for (auto& thr : m_workers)
      {
        if (thr.joinable())
          {
            thr.join();
          }
      }
    for (auto& thr : m_procComms)
      {
        if (thr.joinable())
          {
            thr.join();
          }
      }
    for (auto& thr : m_gpuComms)
      {
        if (thr.joinable())
          {
            thr.join();
          }
      }
  }

/*--------------------------------------------------------------------*/
//  Return are all the threads set
/** \return             Wether all the threads have been set
 *//*-----------------------------------------------------------------*/

  inline
  bool ThreadArchitect::isFull() const
    {
      return m_threadTopology ?
        m_workers.size() == m_threadTopology->numWorkers()
        && m_procComms.size() == m_threadTopology->numProcComms()
        && m_gpuComms.size() == m_threadTopology->numGpuComms()
        : false;
    }
}

#include "NamespaceFooter.H"

#endif  /* ! defined _THREADARCHITECT_H_ */
