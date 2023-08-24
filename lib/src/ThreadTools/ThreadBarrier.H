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
 *  \file ThreadBarrier.H
 *
 *  \brief Barrier for threads
 *
 *//*+*************************************************************************/

#ifndef _THREADBARRIER_H_
#define _THREADBARRIER_H_

//----- Standard Library -----//

#include <atomic>
#include <mutex>
#include <condition_variable>
#include <cassert>
#include <thread>

#include "NamespaceHeader.H"

namespace ThreadTools
{

/*******************************************************************************
 */
/// Thread barrier
/**
 *  This barrier may be used in sequence without fear of conflicts.  Also,
 *  different threads can be safely used at subsequent barriers.
 *
 ******************************************************************************/

class ThreadBarrier
{


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Constructor
  ThreadBarrier(const int a_numThread);

  /// Destructor
  ~ThreadBarrier() = default;

//--Use synthesized copy, copy assignment, move, and move assignment


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Set the number of threads.  Weak construction.
  void define(const int a_numThread);

  /// Wait at the barrier
  void wait();

/*==============================================================================
 * Data members
 *============================================================================*/

private:

  int m_count;                        ///< Downwards count of threads at barrier
  int m_clear;                        ///< Downwards count as threads clear
                                      ///< barrier
  int m_numThread;                    ///< Number of threads in the barrier
  std::mutex m_cvLock;                ///< Lock
  std::condition_variable m_cv;       ///< Condition variable
};


/*******************************************************************************
 *
 * Class ThreadBarrier: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Constructor
/** \param[in]  a_numThread
 *                      Number of threads using the barrier
 *//*-----------------------------------------------------------------*/

inline
ThreadBarrier::ThreadBarrier(const int a_numThread)
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

inline void
ThreadBarrier::define(const int a_numThread)
{
  assert(a_numThread > 0);
  std::unique_lock<std::mutex> lock(m_cvLock);
  if (m_clear != 0)
    {
      m_cv.wait(lock,  // Wait until previous barrier cleared by all threads
                [&l_clear = m_clear]
                {
                  return l_clear == 0;
                });
    }
  const int threadsWaiting = m_numThread - m_count;
  m_count = a_numThread - threadsWaiting;
  assert(m_count > 0);
  m_numThread = a_numThread;
}

/*--------------------------------------------------------------------*/
//  Wait at the barrier
/**
 *//*-----------------------------------------------------------------*/

inline void
ThreadBarrier::wait()
{
  if (m_numThread > 1)
    {
      std::unique_lock<std::mutex> lock(m_cvLock);
      if (m_clear != 0)
        {
          m_cv.wait(lock, // Wait until previous barrier cleared by all threads
                    [&l_clear = m_clear]
                    {
                      return l_clear == 0;
                    });
        }
      --m_count;
      if (m_count == 0)
        {
          m_clear = m_numThread;
          m_cv.notify_all();
        }
      else
        {
          m_cv.wait(lock,  // Wait at this barrier
                    [&l_count = m_count]
                    {
                      return l_count == 0;
                    });
        }
      --m_clear;
      if (m_clear == 0)
        {
          m_count = m_numThread;
          m_cv.notify_all();
        }
    }
}
}

#include "NamespaceFooter.H"

#endif  /* ! defined _THREADBARRIER_H_ */
