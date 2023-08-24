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
 *  \file ThreadBuilder.cpp
 *
 *  \brief Non-inline definitions for classes in ThreadBuilder.H
 *
 *//*+*************************************************************************/

//----- Standard Library -----//

#include <cstring>

//----- Internal -----//

#include "ThreadBarrier.H"
#include "ThreadBuilder.H"

#include <pthread.h>

#include "NamespaceHeader.H"

using namespace ThreadTools;

int ThreadBuilder::s_cnt{0};
thread_local int ThreadBuilder::tls_thrID = 0;

/*******************************************************************************
 *
 * Struct EnclosingArgs: information to pass to as argument to pthread
 *
 ******************************************************************************/

struct EnclosingArgs
{
  void* args;
  void* (*func)(void*);
  std::shared_ptr<const HardwareTopology> hwTopology;
  hwloc_obj_t hwlocBindObj;
  int policy;
  int prio;
  int thrId;
  std::shared_ptr<ThreadBarrier> barrier;
};

/*******************************************************************************
 *
 * Struct EnclosingArgsNoBind: information to pass to as argument to pthread
 *
 ******************************************************************************/

struct EnclosingArgsNoBind
{
  void* args;
  void* (*func)(void*);
  std::shared_ptr<const HardwareTopology> hwTopology;
  int policy;
  int prio;
  int thrId;
  std::shared_ptr<ThreadBarrier> barrier;
};

/*--------------------------------------------------------------------*/
//  Function executed by the thread that calls the provided function
/** \param[in]  args   void pointer to the EnclosingArgs to use to
 *                     setup the thread and execute the provided func
 *//*-----------------------------------------------------------------*/

void* enclosingFunc(void* args)
{
  EnclosingArgs* enArgs = static_cast<EnclosingArgs*>(args);
  ThreadBuilder::tls_thrID = ++ThreadBuilder::s_cnt;
  tls_threadId = enArgs->thrId;

  std::shared_ptr<ThreadBarrier> barrier(enArgs->barrier);

  if (enArgs->prio >= 0) // TODO cleanup this logic
    {
      ThreadSupport::setThreadPrio(enArgs->policy, enArgs->prio);
    }
  ThreadSupport::bindThreadToHW(enArgs->hwTopology.get(), enArgs->hwlocBindObj);

  barrier->wait();
  barrier.reset();
  return enArgs->func(enArgs->args);
}

/*--------------------------------------------------------------------*/
//  Function executed by the thread that calls the provided function
/** \param[in]  args   void pointer to the EnclosingArgs to use to
 *                     setup the thread and execute the provided func
 *//*-----------------------------------------------------------------*/

void* enclosingFuncNoBind(void* args)
{
  EnclosingArgsNoBind* enArgs = static_cast<EnclosingArgsNoBind*>(args);
  ThreadBuilder::tls_thrID = ++ThreadBuilder::s_cnt;
  tls_threadId = enArgs->thrId;

  std::shared_ptr<ThreadBarrier> barrier(enArgs->barrier);

  if (enArgs->prio >= 0) // TODO cleanup this logic
    {
      ThreadSupport::setThreadPrio(enArgs->policy, enArgs->prio);
    }

  barrier->wait();
  barrier.reset();
  return enArgs->func(enArgs->args);
}

/*******************************************************************************
 *
 * Class ThreadBuilder: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** \param[in]  a_hwTopology
 *                      HardwareTopology to use to create threads.
 *                      If null, make one
 *//*-----------------------------------------------------------------*/

ThreadBuilder::ThreadBuilder(
  std::shared_ptr<const HardwareTopology> a_hwTopology)
:
m_hwTopology{a_hwTopology ?
  a_hwTopology : std::make_shared<const HardwareTopology>()}
{
}

/*==============================================================================
 * Public member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Create a thread
/** \param[in]  a_thrId The thread ID
 *  \param[in]  a_func  Function executed by thread
 *  \param[in]  a_args  Args for function passed to thread
 *  \param[in]  a_status
 *                      Holds error codes of system calls
 *//*-----------------------------------------------------------------*/

std::thread ThreadBuilder::createPlainThread(const int    a_thrId,
                                             void*        a_func (void*),
                                             void*        a_args,
                                             int&         a_status)
{
  std::shared_ptr<ThreadBarrier> barrier = std::make_shared<ThreadBarrier>(2);

  auto enclosingFunc = [&, a_barrier = barrier]() {
    ThreadBuilder::tls_thrID = ++ThreadBuilder::s_cnt;
    tls_threadId = a_thrId;

    std::shared_ptr<ThreadBarrier> barrier(a_barrier);

    barrier->wait();
    barrier.reset();
    return a_func(a_args);
  };

  std::thread thr(enclosingFunc);

  barrier->wait();

  return thr;
}

/*--------------------------------------------------------------------*/
//  Create a thread and bind it to a hardware object
/** \param[in]  a_obj   Hardware object to bind to
 *  \param[in]  a_stackSize_B
 *                      Stacksize of the thread
 *  \param[in]  a_guardSize_B
 *                      Extra guard memory for stack overflowHardware object to
 *                      bind to
 *  \param[in]  a_policy
 *                      Scheduling policy
 *  \param[in]  a_prio  Priority of thread set with ThreadSupport
 *  \param[in]  a_thrId The thread ID
 *  \param[in]  a_func  Function executed by thread
 *  \param[in]  a_args  Args for function passed to thread
 *  \param[in]  a_status
 *                      Holds error codes of system calls
 *//*-----------------------------------------------------------------*/

std::thread ThreadBuilder::createThread(hwloc_obj_t  a_obj,
                                        const size_t a_stackSize_B,
                                        const size_t a_guardSize_B,
                                        const int    a_policy,
                                        const int    a_prio,
                                        const int    a_thrId,
                                        void*        a_func (void*),
                                        void*        a_args,
                                        int&         a_status)
{
  pthread_attr_t attr;
  pthread_attr_init(&attr); // Sets most defaults
  a_status = pthread_attr_setschedpolicy(&attr, SCHED_OTHER);
  if (a_status)
    {
      s_thrErr << "error in pthread_attr_setschedpolicy: " << a_status << std::endl;
      return std::thread();
    }
  struct sched_param schedParam;
  schedParam.sched_priority = 0; // Only value allowed for SCHED_OTHER
  a_status = pthread_attr_setschedparam(&attr, &schedParam);
  if (a_status)
    {
      s_thrErr << "error in pthread_attr_setschedparam: " << a_status << std::endl;
      return std::thread();
    }
  a_status = pthread_attr_setinheritsched(&attr, PTHREAD_EXPLICIT_SCHED);
  if (a_status)
    {
      s_thrErr << "allocating stack memory for pthread failed: " << a_status << std::endl;
      return std::thread();
    }

  void* addr = ThreadSupport::allocStackMemory(m_hwTopology.get(),
                                               a_obj,
                                               a_stackSize_B,
                                               a_guardSize_B,
                                               a_status);
  if (a_status)
    {
      s_thrErr << "allocating stack memory for pthread failed: " << a_status << std::endl;
      return std::thread();
    }
  a_status = pthread_attr_setstack(&attr, addr, a_stackSize_B);
  if (a_status)
    {
      s_thrErr << "pthread failed to set stack: " << a_status << std::endl;
      return std::thread();
    }

  std::shared_ptr<ThreadBarrier> barrier = std::make_shared<ThreadBarrier>(2);

  EnclosingArgs* enclosingArgs = new EnclosingArgs{a_args,
                                                   a_func,
                                                   m_hwTopology,
                                                   a_obj,
                                                   a_policy,
                                                   a_prio,
                                                   a_thrId,
                                                   barrier};

  pthread_t pthr;
  a_status = pthread_create(&pthr, &attr, &enclosingFunc, enclosingArgs);

  barrier->wait();

  // Highly non-portable, bad in every way, but only way to set stack
  // and it seems to work.
  std::thread thr;
  std::memcpy(&thr, &pthr, sizeof(std::thread));

  ThreadSupport::printPThreadAttr(pthr);
  pthread_attr_destroy(&attr);
  return thr;
}


/*--------------------------------------------------------------------*/
//  Create a thread and set its priority and bind it to a hardware object
/** \param[in]  a_policy
 *                      Scheduling policy
 *  \param[in]  a_prio  Priority of thread set with ThreadSupport
 *  \param[in]  a_thrId The thread ID
 *  \param[in]  a_func  Function executed by thread
 *  \param[in]  a_args  Args for function passed to thread
 *  \param[in]  a_status
 *                      Holds error codes of system calls
 *//*-----------------------------------------------------------------*/

std::thread ThreadBuilder::createThreadNoBind(const int    a_policy,
                                              const int    a_prio,
                                              const int    a_thrId,
                                              void*        a_func (void*),
                                              void*        a_args,
                                              int&         a_status)
{
  pthread_attr_t attr;
  pthread_attr_init(&attr); // Sets most defaults
  a_status = pthread_attr_setschedpolicy(&attr, SCHED_OTHER);
  if (a_status)
    {
      s_thrErr << "error in pthread_attr_setschedpolicy: " << a_status << std::endl;
      return std::thread();
    }
  struct sched_param schedParam;
  schedParam.sched_priority = 0; // Only value allowed for SCHED_OTHER
  a_status = pthread_attr_setschedparam(&attr, &schedParam);
  if (a_status)
    {
      s_thrErr << "error in pthread_attr_setschedparam: " << a_status << std::endl;
      return std::thread();
    }
  a_status = pthread_attr_setinheritsched(&attr, PTHREAD_EXPLICIT_SCHED);
  if (a_status)
    {
      s_thrErr << "pthread attr set inherit sched failed: " << a_status << std::endl;
      return std::thread();
    }

  std::shared_ptr<ThreadBarrier> barrier = std::make_shared<ThreadBarrier>(2);

  EnclosingArgsNoBind* enclosingArgs = new EnclosingArgsNoBind{a_args,
                                                               a_func,
                                                               m_hwTopology,
                                                               a_policy,
                                                               a_prio,
                                                               a_thrId,
                                                               barrier};

  pthread_t pthr;
  a_status = pthread_create(&pthr, &attr, &enclosingFuncNoBind, enclosingArgs);

  barrier->wait();

  // Highly non-portable, bad in every way, but only way to set stack
  // and it seems to work.
  std::thread thr;
  std::memcpy(&thr, &pthr, sizeof(std::thread));

  ThreadSupport::printPThreadAttr(pthr);
  pthread_attr_destroy(&attr);
  return thr;
}

#include "NamespaceFooter.H"
