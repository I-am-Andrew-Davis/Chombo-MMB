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
 *  \file ThreadBuilder.H
 *
 *  \brief Builds threads
 *
 *//*+*************************************************************************/

#ifndef _THREAD_BUILDER_H_
#define _THREAD_BUILDER_H_

//----- Standard Library -----//

#include <thread>

//----- Internal -----//

#include "ThreadSupport.H"
#include "HardwareTopology.H"

#include "NamespaceHeader.H"

namespace ThreadTools
{

/*******************************************************************************
 */
/// Builds threads
/**
 *  Builds threads and binds them to a hardware object
 *
 ******************************************************************************/

  class ThreadBuilder
  {


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  public:

    /// Constructor
    ThreadBuilder(std::shared_ptr<const HardwareTopology> hwTopology);
  
    /// Destructor
    ~ThreadBuilder() = default;

/*==============================================================================
 * Public member functions
 *============================================================================*/

    /// Create a thread and bind to a hwloc_obj_t
    std::thread createThread(hwloc_obj_t a_obj,
                             size_t      a_stackSize_B,
                             size_t      a_guardSize_B,
                             int         a_policy,
                             int         a_prio,
                             int         a_thrId,
                             void* a_func (void*),
                             void*       a_args,
                             int&        status);

    /// Create a thread and set its policy
    std::thread createThreadNoBind(int   a_policy,
                                   int   a_prio,
                                   int   a_thrId,
                                   void* a_func (void*),
                                   void* a_args,
                                   int&  status);

    /// Create a thread
    std::thread createPlainThread(int   a_thrId,
                                  void* a_func (void*),
                                  void* a_args,
                                  int&  status);

/*==============================================================================
 * Data members
 *============================================================================*/

  public:
    static int s_cnt;                 ///< Number of threads built
    static thread_local int tls_thrID;
                                      ///< Each thread has its thread id

  protected:
    std::shared_ptr<const HardwareTopology> m_hwTopology;
                                      ///< Hardware topology is needed to
                                      ///< setup each thread
  };

} // namespace Threads

#include "NamespaceFooter.H"

#endif  /* ! defined _THREAD_BUILDER_H_ */