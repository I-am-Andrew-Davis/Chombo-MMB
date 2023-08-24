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
 *  \file ThreadSpinBarrier.H
 *
 *  \brief Spin barrier for threads
 *
 *//*+*************************************************************************/


#ifndef _THREADSPINBARRIER_H_
#define _THREADSPINBARRIER_H_

//----- Standard Library -----//

#include <atomic>
#include <cassert>
#include <condition_variable>
#include <mutex>
#include <thread>

#include "NamespaceHeader.H"

namespace ThreadTools
{

/*******************************************************************************
 */
/// Thread spin barrier
/**
 *  This barrier may be used in sequence without fear of conflicts.  Also,
 *  different threads can be safely used at subsequent barriers.
 *
 ******************************************************************************/

class ThreadSpinBarrier
{


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:
  /// Default constructor
  ThreadSpinBarrier();

  /// Construct with parameters
  ThreadSpinBarrier(const int a_numThread);

  /// Destructor
  ~ThreadSpinBarrier() = default;

//--Use synthesized copy, copy assignment, move, and move assignment


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:
  /// Set the number of threads.  Weak construction.
  void define(const int a_numThread);

  /// Wait at the barrier
  void spinWait();

/*==============================================================================
 * Data members
 *============================================================================*/

private:

  std::atomic<int> m_count;           ///< Downwards count of threads at barrier
  std::atomic<int> m_clear;           ///< Downwards count as threads clear
                                      ///< barrier
  int m_numThread;                    ///< Number of threads in the barrier
  std::atomic_flag m_lock = ATOMIC_FLAG_INIT;
                                      ///< lock
};


/*******************************************************************************
 *
 * Class ThreadSpinBarrier: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Constructor
/** \param[in]  a_numThread
 *                      Number of threads using the barrier
 *//*-----------------------------------------------------------------*/

inline ThreadSpinBarrier::ThreadSpinBarrier()
  :
  m_count(0),
  m_clear(0),
  m_numThread(0)
{
}

/*--------------------------------------------------------------------*/
//  Constructor
/** \param[in]  a_numThread
 *                      Number of threads using the barrier
 *//*-----------------------------------------------------------------*/

inline ThreadSpinBarrier::ThreadSpinBarrier(const int a_numThread)
  :
  m_count(a_numThread),
  m_clear(0),
  m_numThread(a_numThread)
{
  assert(a_numThread > 0);
}

/*--------------------------------------------------------------------*/
//  Set the number of threads (dangerous!).  Weak construction.
/** \param[in]  a_numThread
 *                      Number of threads at the barrier
 *  This is relatively safe as long as the thread calling define also
 *  waits at the barrier.  The code will abort if one tries to set
 *  a_numThread to a value less than the number of threads currently
 *  waiting.
 *//*-----------------------------------------------------------------*/

inline void ThreadSpinBarrier::define(const int a_numThread)
{
  assert(a_numThread > 0);
  while (m_lock.test_and_set(std::memory_order_acquire));
  if (m_clear != 0)
    {
      m_lock.clear(std::memory_order_release);
      while (m_clear)
        {
        }
      while (m_lock.test_and_set(std::memory_order_acquire));
    }
  const int threadsWaiting = m_numThread - m_count;
  m_count = a_numThread - threadsWaiting;
  assert(m_count > 0);
  m_numThread = a_numThread;
  m_lock.clear(std::memory_order_release);
}

/*--------------------------------------------------------------------*/
//  Wait at the barrier
/**
 *//*-----------------------------------------------------------------*/

inline void ThreadSpinBarrier::spinWait()
{
  if (m_numThread > 1)
    {
      while (m_clear)
        {
        }
      while (m_lock.test_and_set(std::memory_order_acquire));
      if (m_count.fetch_sub(1) == 1)
        {
          m_clear += m_numThread - 1;
        }
      else
        {
          if (m_count)
            {
              m_lock.clear(std::memory_order_release);
              while (m_count)
                {
                }
            }

          while (m_lock.test_and_set(std::memory_order_acquire));
          if (m_clear.fetch_sub(1) == 1)
            {
              m_count += m_numThread;
            }
        }
      m_lock.clear(std::memory_order_release);
    }
}
} // namespace Threads

#include "NamespaceFooter.H"

#endif /* ! defined _THREADSPINBARRIER_H_ */
