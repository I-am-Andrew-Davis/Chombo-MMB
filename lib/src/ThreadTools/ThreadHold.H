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
 *  \file ThreadHold.H
 *
 *  \brief Simple CV for threads
 *
 *//*+*************************************************************************/

#ifndef _THREADHOLD_H_
#define _THREADHOLD_H_

//----- Standard Library -----//

#include <mutex>
#include <condition_variable>

#include "NamespaceHeader.H"

namespace ThreadTools
{

/*******************************************************************************
 */
/// Thread condition variable
/**
 *  Simplest of CVs but packages all three components: flag, lock, and CV.
 *  Use to hold a single or multiple threads and then release with command.
 *
 ******************************************************************************/

class ThreadHold
{


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Default constructor
  ThreadHold();

  /// Destructor
  ~ThreadHold() = default;

  ThreadHold(const ThreadHold&) = delete;
  ThreadHold& operator=(const ThreadHold&) = delete;

/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Reset the release flag
  void reset();

  /// Wait
  void wait();

  /// Release any waiting threads
  void release();

/*==============================================================================
 * Data members
 *============================================================================*/

private:

  bool m_release;                     ///< T - release
  std::mutex m_cvLock;                ///< Lock
  std::condition_variable m_cv;       ///< Condition variable
};


/*******************************************************************************
 *
 * Class ThreadHold: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Default constructor
/**
 *//*-----------------------------------------------------------------*/

inline
ThreadHold::ThreadHold()
  :
  m_release(false)
{ }

/*--------------------------------------------------------------------*/
//  Reset the release flag
/**
 *//*-----------------------------------------------------------------*/

inline void
ThreadHold::reset()
{
  m_cvLock.lock();
  m_release = false;
  m_cvLock.unlock();
}

/*--------------------------------------------------------------------*/
//  Wait
/**
 *//*-----------------------------------------------------------------*/

inline void
ThreadHold::wait()
{
  std::unique_lock<std::mutex> lock(m_cvLock);
  m_cv.wait(lock,
            [&l_release = m_release]
            {
              return l_release;
            });
}

/*--------------------------------------------------------------------*/
//  Release any waiting threads
/**
 *//*-----------------------------------------------------------------*/

inline void
ThreadHold::release()
{
  m_cvLock.lock();
  m_release = true;
  m_cvLock.unlock();
  m_cv.notify_all();
}
}

#include "NamespaceFooter.H"

#endif  /* ! defined _THREADHOLD_H_ */
