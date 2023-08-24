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
 *  \file ThreadTeam.H
 *
 *  \brief Definition of class ThreadTeam
 *
 *//*+*************************************************************************/

#ifndef _THREADTEAM_H_
#define _THREADTEAM_H_

//----- Standard Library -----//

#include <mutex>
#include <condition_variable>

//----- Internal -----//

#include "ThreadBarrier.H"
#include "ThreadHold.H"

#include "NamespaceHeader.H"

namespace ThreadTools
{

/*******************************************************************************
 */
/// A thread team that uses the same function for each member
/**
 *
 ******************************************************************************/

  class ThreadTeam
  {


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  public:

    /// Constructor
    ThreadTeam(unsigned a_numMembers)
      :
      m_numMembers(a_numMembers),
      m_barrier(std::make_unique<ThreadBarrier>(a_numMembers)),
      m_hold(std::make_unique<ThreadHold>())
      {
      }
    
    /// Explicitly delete copy constructor
    ThreadTeam(const ThreadTeam&) = delete;

    /// Move constructor
    ThreadTeam(ThreadTeam&& other)
      :
      m_numMembers(other.m_numMembers),
      m_barrier(std::move(other.m_barrier)),
      m_task{other.m_task.load()}
      {
      }

    /// Destructor
    ~ThreadTeam() = default;


/*==============================================================================
 * Public member functions
 *============================================================================*/

  public:

    /// Wait for all threads
    void waitForAllThreads()
      {
        m_barrier->wait();
      }

    /// Thread(s) will wait until released
    void hold()
      {
        m_hold->wait();
      }

    /// Release the threads waiting at hold
    void release()
      {
        m_hold->release();
      }

    /// Reset the hold to be used again
    void reset()
      {
        m_hold->reset();
      }

    /// Get the task
    void* task()
      {
        return m_task;
      }
    
    /// Set the task
    void setTask(void* a_task)
      {
        m_task = a_task;
      }

    /// Get the number of members in the team
    unsigned size()
      {
        return m_numMembers;
      }

/*==============================================================================
 * Data members
 *============================================================================*/

  protected:

    unsigned m_numMembers;            ///< The number of threads in the team
    std::unique_ptr<ThreadBarrier> m_barrier;
                                      ///< Barrier for threads wait until all
                                      ///< threads reach it.  Use a pointer
                                      ///< to allow this class to be moved.
                                      ///< TODO barrier or spinbarrier
    std::unique_ptr<ThreadHold> m_hold;
                                      ///< Hold for threads to wait until
                                      ///< explicitly released.  Use a pointer
                                      ///< to allow this class to be moved
    std::atomic<void*> m_task;        ///< Task that each thread will execute
  };
}

#include "NamespaceFooter.H"

#endif  /* ! defined _THREADTEAM_H_ */
