
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
 *  \file ThreadSupport.H
 *
 *  \brief Definition of class ThreadSupport
 *
 *//*+*************************************************************************/

#ifndef _THREADSUPPORT_H_
#define _THREADSUPPORT_H_

//----- Standard Library -----//

#include <cstdint>

//----- System -----//

#include <hwloc.h>

#include "NamespaceHeader.H"

namespace ThreadTools
{

//-- Forward declarations
class HardwareTopology;

/*******************************************************************************
 */
/// Thread support
/**
 *  Handles topology, system, binding of threads to hardware, and priority
 *  of threads.
 *
 ******************************************************************************/

class ThreadSupport
{


/*==============================================================================
 * Constructors and destructors
 *============================================================================*/

  /// Prohibit construction of this class
  ThreadSupport() = delete;
  ~ThreadSupport() = delete;


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Get the hwloc api version string (X.Y.Z)
  static std::string hwlocVersion();

  /// Check alignment of some memory
  static bool isAligned(void *const          a_p,
                        const std::uintptr_t a_align) noexcept;

  /// Grow a chunk of memory so it is a multiple of a_align
  static constexpr std::size_t alignTo(std::size_t a_x,
                                       std::size_t a_align) noexcept;
  
  /// Name of hwloc membind policy
  static constexpr const char* membindPolicyName(
      const hwloc_membind_policy_t a_policy) noexcept;
  
  /// Name of a Linux schedule policy
  static constexpr const char* schedPolicyName(int a_policy) noexcept;
  
  /// Bind a thread to a hardware object
  static void bindThreadToHW(HardwareTopology const* const a_topology,
                             hwloc_obj_t                   a_obj) noexcept;

  /// Set thread priority
  static int setThreadPrio(int a_policy,
                           int a_prio) noexcept;

  /// Return if we can set realtime priority for threads
  static bool hasRTPrio() noexcept;

  /// Allocate stack memory for a thread
  static
  void* allocStackMemory(HardwareTopology const* const a_topology,
                          hwloc_obj_t                   a_obj,
                          size_t                        a_stackSize_B,
                          size_t                        a_guardSize_B,
                          int&                          a_status) noexcept;

  /// Free stack memory for a thread
  static void freeStackMemory(const hwloc_topology_t a_topology,
                              void*&                 a_addr,
                              size_t                 a_stackSize_B,
                              size_t                 a_guardSize_B) noexcept;

  /// Test if thread priority can be set
  static bool testRLimit() noexcept;

  /// Print thread CPU bindings
  static void printThreadCPUBinding(const hwloc_topology_t topology) noexcept;

  /// Print thread memory bindings
  static
  void printThreadMemoryBinding(const hwloc_topology_t topology) noexcept;

  /// Print thread memory bindings
  static void printPThreadAttr(pthread_t a_pthr) noexcept;

  /// Set the system realtime priorities
  static void setRealtimePriorities();


/*==============================================================================
 * Data members
 *============================================================================*/

private:

  static constexpr int c_strsz = 128; ///< Size of character arrays used to
                                      ///< retrieve names from hwloc
  static bool s_hasRTPrio;            ///< Has modifiable realtime priority
                                      ///< limits
  static bool s_hasTestedRTPrio;      ///< Has realtime priority limits been
                                      ///< tested for
};


/*******************************************************************************
 *
 * Class ThreadBarrier: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Check alignment of some memory
/** \param[in]  a_p     Memory address to check
 *  \param[in]  a_align Alignment value in bytes
 *  \return             T - aligned
 *                      F - not aligned
 *//*-----------------------------------------------------------------*/

inline bool  // Clang does not like this as constexpr
ThreadSupport::isAligned(void *const a_p, const std::uintptr_t a_align) noexcept
{
  return (reinterpret_cast<std::uintptr_t>(a_p) & (a_align - 1)) == 0;
}

/*--------------------------------------------------------------------*/
//  Grow a chunk of memory so it is a multiple of a_align
/** \param[in]  a_x     Chuck of memory in bytes
 *  \param[in]  a_align Alignment value in bytes
 *  \return             Chunk grow to be a multiple of a_align
 *//*-----------------------------------------------------------------*/

inline constexpr std::size_t
ThreadSupport::alignTo(std::uintptr_t a_x,
                       std::uintptr_t a_align) noexcept
{
  return (a_x + (a_align - 1)) & ~(a_align - 1);
}

/*--------------------------------------------------------------------*/
//  Return if we can set realtime priority for threads
/** \return             Whether realtime priority can be set
 *//*-----------------------------------------------------------------*/

inline
bool ThreadSupport::hasRTPrio() noexcept
{
  if (!s_hasTestedRTPrio)
  {
    testRLimit();
  }
  return s_hasRTPrio;
}

}

#include "NamespaceFooter.H"

#endif  /* ! defined _THREADSUPPORT_H_ */
