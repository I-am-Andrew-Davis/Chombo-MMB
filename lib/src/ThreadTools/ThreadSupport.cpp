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
 *  \file ThreadSupport.cpp
 *
 *  \brief Non-inline definitions for ThreadSupport
 *
 *//*+*************************************************************************/


//----- Standard Library -----//

#include <cstdio>
#include <sstream>

//----- System -----//

#include <unistd.h>
#include <sys/resource.h>
#include <sys/mman.h>
#include <pthread.h>

//----- Internal -----//

#include "ThreadSupport.H"
#include "ThreadParameters.H"
#include "HardwareTopology.H"

#include "NamespaceHeader.H"

using namespace ThreadTools;

//--Initialize static members

bool ThreadSupport::s_hasRTPrio = false;
bool ThreadSupport::s_hasTestedRTPrio = false;


/*==============================================================================
 * Public member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Return version string of hwloc API version
/** return              Version (X.Y.Z)
 *//*-----------------------------------------------------------------*/

std::string ThreadSupport::hwlocVersion()
{
  uint Z = HWLOC_API_VERSION % 256;
  uint Y = (HWLOC_API_VERSION & 0xFF00) >> 8;
  uint X = (HWLOC_API_VERSION & 0xFF0000) >> 16;

  std::stringstream ss;
  ss << X << '.' << Y << '.' << Z;
  return ss.str();
}

/*--------------------------------------------------------------------*/
//  Name of hwloc membind policy
/** \param[in]  a_policy
 *                      hwloc membind policy enumeration
 *  return              Name
 *//*-----------------------------------------------------------------*/

constexpr const char*
ThreadSupport::membindPolicyName(
  const hwloc_membind_policy_t a_policy) noexcept
{
#if HWLOC_API_VERSION >= 0x00020000
  switch (a_policy)
    {
    case HWLOC_MEMBIND_DEFAULT:
      return "default";
    case HWLOC_MEMBIND_FIRSTTOUCH:
      return "first touch";
    case HWLOC_MEMBIND_BIND:
      return "bind";
    case HWLOC_MEMBIND_INTERLEAVE:
      return "interleave";
    case HWLOC_MEMBIND_NEXTTOUCH:
      return "next touch";
    case HWLOC_MEMBIND_MIXED:
      return "mixed";
    default:
      return "unknown";
    }
#else
  switch (a_policy)
    {
    case HWLOC_MEMBIND_DEFAULT:
      return "default";
    case HWLOC_MEMBIND_FIRSTTOUCH:
      return "first touch";
    case HWLOC_MEMBIND_BIND:
      return "bind";
    case HWLOC_MEMBIND_INTERLEAVE:
      return "interleave";
    case HWLOC_MEMBIND_REPLICATE:
      return "replicate";
    case HWLOC_MEMBIND_NEXTTOUCH:
      return "next touch";
    case HWLOC_MEMBIND_MIXED:
      return "mixed";
    default:
      return "unknown";
    }
#endif

}

/*--------------------------------------------------------------------*/
//  Name of a Linux schedule policy
/** \param[in]  a_policy
 *                     Linux schedule policy
 *  return             Name
 *//*-----------------------------------------------------------------*/

constexpr const char*
ThreadSupport::schedPolicyName(const int a_policy) noexcept
{
  switch (a_policy)
    {
#ifdef SCHED_OTHER
    case SCHED_OTHER:
      return "SCHED_OTHER";
#endif
#ifdef SCHED_IDLE
    case SCHED_IDLE:
      return "SCHED_IDLE";
#endif
#ifdef SCHED_BATCH
    case SCHED_BATCH:
      return "SCHED_BATCH";
#endif
#ifdef SCHED_FIFO
    case SCHED_FIFO:
      return "SCHED_FIFO";
#endif
#ifdef SCHED_RR
    case SCHED_RR:
      return "SCHED_RR";
#endif
    default:
      return "unknown";
    }
}

/*--------------------------------------------------------------------*/
//  Bind a thread to a hardware object
/** \param[in]  a_topology
 *                     The hardware topology
 *  \param[in]  a_obj  The object to bind the thread to
 *//*-----------------------------------------------------------------*/

void
ThreadSupport::bindThreadToHW(const HardwareTopology* a_topology,
                              hwloc_obj_t a_obj) noexcept
{
  const int idThr = ThreadTools::tls_threadId;

  char objType[2][c_strsz];
  char listSet[c_strsz];
  hwloc_cpuset_t set = hwloc_bitmap_alloc();

  const bool hasNuma = a_topology->hasNUMA();
  const hwloc_topology_t hwlocTopology = a_topology->topology();

  std::stringstream ss;

  // Find object to bind to and print intended CPU bindings
  hwloc_obj_type_snprintf(objType[0], sizeof(objType[0]), a_obj, 0);
  unsigned osIndexCore = a_obj->os_index;
  // Children of core, presumably hardware support for threads
  if (a_obj->type == HWLOC_OBJ_CORE && a_obj->arity)
    {
      a_obj = a_obj->children[0];
      hwloc_obj_type_snprintf(objType[1], sizeof(objType[1]), a_obj, 0);
      hwloc_bitmap_list_snprintf(listSet, sizeof(listSet), a_obj->cpuset);
      if (s_verbosity >= Verb::v2)
        {
          ss << "Binding " << idThr << " to "
                 << objType[1] << a_obj->os_index << " of "
                 << objType[0] <<  osIndexCore << " with cpuset "
                 << listSet << "\n";
        }
    }
  else  // Bind to CPU or core
    {
      hwloc_bitmap_list_snprintf(listSet, sizeof(listSet), a_obj->cpuset);
      if (s_verbosity >= Verb::v2)
        {
          ss << "Binding " << idThr << " to "
                 << objType[0] <<  osIndexCore << " with cpuset "
                 << listSet << "\n";
        }
    }

  // Print intended memory bindings
  if (hasNuma)
    {
      hwloc_bitmap_list_snprintf(listSet, sizeof(listSet), a_obj->nodeset);
      if (strlen(listSet) == 2 && listSet[1] == '-')
        {
          listSet[1] = '\0';
        }
      if (s_verbosity >= Verb::v2)
        {
          ss << "Binding " << idThr << " memory to nodeset "
                 << listSet << "\n";
        }
    }

  // Print current CPU binding of this thread
  hwloc_get_cpubind(hwlocTopology, set, HWLOC_CPUBIND_THREAD);
  hwloc_bitmap_list_snprintf(listSet, sizeof(listSet), set);
  if (s_verbosity >= Verb::v3)
    {
      ss << "  Current CPU binding                 : (" << listSet << ")\n";
    }

  // Print current memory binding of this thread
  if (hasNuma)
    {
      hwloc_membind_policy_t policy;
#if HWLOC_API_VERSION >= 0x00020100
      hwloc_get_membind(hwlocTopology,
                        set,
                        &policy,
                        HWLOC_MEMBIND_THREAD | HWLOC_MEMBIND_BYNODESET);
#else
      hwloc_get_membind_nodeset(hwlocTopology,
                                set,
                                &policy,
                                HWLOC_MEMBIND_THREAD);
#endif
      hwloc_bitmap_list_snprintf(listSet, sizeof(listSet), set);
      if (strlen(listSet) == 2 && listSet[1] == '-')
        {
          listSet[1] = '\0';
        }
      if (s_verbosity >= Verb::v2)
        {
          ss << "  Current memory binding              : (" << listSet
                 << ") " << membindPolicyName(policy) << "\n";
        }
    }

  // This thread should be binded to a cpuset included in the cpuset of 
  // the process
  // FIXME might be using puset for hyperthreading
  // CH_assert(hwloc_bitmap_isincluded(a_obj->cpuset, m_cpuset));
  
  // Bind and print new CPU binding of this thread
  if (hwloc_set_cpubind(hwlocTopology, a_obj->cpuset,
                        HWLOC_CPUBIND_THREAD |
                        HWLOC_CPUBIND_STRICT))
    {
      s_thrErr << "CPU binding failed" << std::endl;
    }
  else
    {
      hwloc_get_cpubind(hwlocTopology, set, HWLOC_CPUBIND_THREAD);
      hwloc_bitmap_list_snprintf(listSet, sizeof(listSet), set);
      if (s_verbosity >= Verb::v1)
        {
          ss << "  New CPU binding                     : (" << listSet
                   << ')' << "\n";
        }
    }

  // Bind and print new memory binding of this thread
  if (hasNuma)
    {

#if HWLOC_API_VERSION >= 0x00020100
      const int memBindFailed = hwloc_set_membind(hwlocTopology,
                                                  a_obj->nodeset,
                                                  HWLOC_MEMBIND_BIND,
                                                  HWLOC_MEMBIND_THREAD |
                                                  HWLOC_MEMBIND_STRICT |
                                                  HWLOC_MEMBIND_BYNODESET);
#else
      const int memBindFailed =
        hwloc_set_membind_nodeset(hwlocTopology, a_obj->nodeset,
                                  HWLOC_MEMBIND_BIND,
                                  HWLOC_MEMBIND_THREAD |
                                  HWLOC_MEMBIND_STRICT);
#endif

      if (memBindFailed)
        {
          s_thrErr << "Memory binding failed" << std::endl;
        }
      else
        {
          hwloc_membind_policy_t policy;
#if HWLOC_API_VERSION >= 0x00020100
          hwloc_get_membind(hwlocTopology,
                            set,
                            &policy,
                            HWLOC_MEMBIND_THREAD | HWLOC_MEMBIND_BYNODESET);
#else
          hwloc_get_membind_nodeset(hwlocTopology,
                                    set,
                                    &policy,
                                    HWLOC_MEMBIND_THREAD);
#endif
          hwloc_bitmap_list_snprintf(listSet, sizeof(listSet), set);
          if (strlen(listSet) == 2 && listSet[1] == '-')
            {
              listSet[1] = '\0';
            }
          if (s_verbosity >= Verb::v3)
            {
              ss << "  New memory binding                  : (" << listSet
                       << ") " << membindPolicyName(policy) << "\n";
            }
        }
    }

  s_thrOut << ss.rdbuf() << std::flush;

  hwloc_bitmap_free(set);
}

/*--------------------------------------------------------------------*/
//  Set thread priority
/** \param[in]  a_policy
 *                      Scheduling policy.  Both SCHED_FIFO and
 *                      SCHED_OTHER are supported but only if
 *                      privileges allow.
 *  \param[in]  a_prio  Priority
 *                      SCHED_FIFO  : (low)  1 -- 5  (high)
 *                      SCHED_OTHER : (high) 0 -- 19 (low) 
 *  \return             0 : success
 *                      1 : error
 *                      2 : insufficient privileges
 *  This routine sets the priority of the calling thread.
 *//*-----------------------------------------------------------------*/

int
ThreadSupport::setThreadPrio(const int        a_policy,
                             const int        a_prio) noexcept
{
  if (!s_hasTestedRTPrio)
  {
    testRLimit();
  }

  std::stringstream ss;

  switch (a_policy)
    {
    case SCHED_FIFO:
      if (!s_hasRTPrio)
        {
          return 2;
        }
      if (a_prio < 1 || a_prio > 5)
        {
          s_thrErr << "Invalid priority" << std::endl;
          return 1;
        }
      if (sched_get_priority_min(SCHED_FIFO) > 1 ||
          sched_get_priority_max(SCHED_FIFO) < 5)
        {
          s_thrErr << "Incompatible system priorities" << std::endl;
          return 1;
        }
      break;
    case SCHED_OTHER:
      if (a_prio < 0 || a_prio > 19)
        {
          s_thrErr << "Invalid priority" << std::endl;
          return 1;
        }
      break;
    }

  pthread_t pthr = pthread_self();
  int schedPolicy;
  struct sched_param schedParam;
  if (pthread_getschedparam(pthr, &schedPolicy, &schedParam))
    {
      s_thrErr << "pthread_getschedparam" << std::endl;
      return 1;
    }
  if (s_verbosity >= Verb::v3)
    {
      ss << "  Current schedule                    : "
               << schedPolicyName(schedPolicy) << "\n";
      ss << "  Current priority                    : "
               << schedParam.sched_priority << "\n";
    }
  schedParam.sched_priority = a_prio;
  int err = pthread_setschedparam(pthr, a_policy, &schedParam);
  if (err)
    {
      s_thrErr << "pthread_setschedparam: ";
      switch (err)
        {
        case ESRCH:
          s_thrErr << "thread not found" << std::endl;
          return 1;
        case EINVAL:
          s_thrErr << "unrecognized policy: " << schedPolicyName(schedPolicy)
                   << " " << a_prio << std::endl;
          return 1;
        case EPERM:
          s_thrErr << "insufficient privileges for FIFO policy" << std::endl;
          return 2;
        default:
          s_thrErr << "unknown reasons" << std::endl;
          return 1;
        }
    }
  if (pthread_getschedparam(pthr, &schedPolicy, &schedParam))
    {
      s_thrErr << "pthread_getschedparam" << std::endl;
      return 1;
    }
    if (s_verbosity >= Verb::v3)
      {
        ss << "  New schedule                        : "
                 << schedPolicyName(schedPolicy) << "\n";
        ss << "  New priority                        : "
                 << schedParam.sched_priority << "\n";
      }
  s_thrOut << ss.rdbuf() << std::flush;
  return 0;
}

/*--------------------------------------------------------------------*/
//  Allocate stack memory for a thread
/** \param[in]  a_topology
 *                      The hardware topology
 *  \param[in]  a_obj   hwloc object describing the nodeset to bind to
 *  \param[in]  a_stackSize_B
 *                      The stacksize for the thread
 *  \param[in]  a_guardSize_B
 *                      The size of the guard for stack overflow
 *  \param[out] a_status
 *                      0 : success
 *                      1 : allocation failure
 *                      2 : binding failure
 *                      3 : failure to protect guard page
 *//*-----------------------------------------------------------------*/

void*
ThreadSupport::allocStackMemory(const HardwareTopology* const a_topology,
                                hwloc_obj_t                   a_obj,
                                const size_t                  a_stackSize_B,
                                const size_t                  a_guardSize_B,
                                int&                          a_status) noexcept
{
  const hwloc_topology_t hwlocTopology = a_topology->topology();

  std::stringstream ss;

#if HWLOC_API_VERSION >= 0x00020100
  void* addr = hwloc_alloc_membind(hwlocTopology,
                                   a_stackSize_B + a_guardSize_B,
                                   a_obj->nodeset,
                                   HWLOC_MEMBIND_BIND,
                                   HWLOC_MEMBIND_STRICT |
                                   HWLOC_MEMBIND_NOCPUBIND |
                                   HWLOC_MEMBIND_BYNODESET);
#else
  void* addr = hwloc_alloc_membind_nodeset(hwlocTopology,
                                           a_stackSize_B + a_guardSize_B,
                                           a_obj->nodeset,
                                           HWLOC_MEMBIND_BIND,
                                           HWLOC_MEMBIND_STRICT |
                                           HWLOC_MEMBIND_NOCPUBIND);
#endif

  if (addr == nullptr)
    {
      // Find object to bind to and print intended CPU bindings
      char objType[c_strsz];
      hwloc_obj_type_snprintf(objType, sizeof(objType), a_obj, 0);
      
      switch (errno)
        {
        case ENOSYS:
          s_thrErr << "Binding action or policy is not supported by OS"
                     << std::endl;
          a_status = 2;
          break;
        case EXDEV:
          s_thrErr << "Requested binding cannot be enforced"
                     << std::endl;
          a_status = 2;
          break;
        case ENOMEM:
          s_thrErr << "Allocation failed even before attempting binding"
                     << std::endl;
          a_status = 1;
          break;
        }
      s_thrErr << "on " << objType << a_obj->os_index << std::endl;

      return addr;
    }
  assert(isAligned(addr, a_topology->getPageSize()));
  if (s_verbosity >= Verb::v2)
    {
      ss << "Stack allocated from " << (void*)((char*)addr + a_guardSize_B)
               << " to "
               << (void*)((char*)addr + a_guardSize_B + a_stackSize_B - 1)
               << "\n";
    }
  // Set last page as a guard
  if (mprotect(addr,
               a_guardSize_B,
               PROT_NONE))
    {
      printf("Error: mprotect failed\n");
      hwloc_free(hwlocTopology, addr, a_stackSize_B + a_guardSize_B);
      a_status = 3;
      return nullptr;
    }
    if (s_verbosity >= Verb::v3)
      {
        ss << "Guard allocated from " << addr << " to "
                 << (void*)((char*)addr + a_guardSize_B - 1) << "\n";
      }
  s_thrOut << ss.rdbuf() << std::flush;
  return static_cast<char*>(addr) + a_guardSize_B;
}

/*--------------------------------------------------------------------*/
//  Free stack memory for a thread
/** \param[in]  a_topology
 *                     The hardware topology
 *  \param[out] a_addr Set to nullptr
 *  \param[in]  a_stackSize_B
 *                     The stack size
 *  \param[in]  a_guardSize_B 
 *                     Size of guard for stack overflow
 *//*-----------------------------------------------------------------*/

void
ThreadSupport::freeStackMemory(const hwloc_topology_t a_topology,
                               void*&                 a_addr,
                               const size_t           a_stackSize_B,
                               const size_t           a_guardSize_B) noexcept
{
  hwloc_free(a_topology, a_addr, a_stackSize_B + a_guardSize_B);
  a_addr = nullptr;
}

/*--------------------------------------------------------------------*/
//  Print thread CPU bindings
/** \param[in] a_topology
 *                     The hardware topology
 *//*-----------------------------------------------------------------*/

void
ThreadSupport::printThreadCPUBinding(const hwloc_topology_t a_topology) noexcept
{
  char listSet[c_strsz];
  hwloc_cpuset_t set = hwloc_bitmap_alloc();
  hwloc_get_cpubind(a_topology, set, HWLOC_CPUBIND_THREAD);
  hwloc_bitmap_list_snprintf(listSet, sizeof(listSet), set);
  std::stringstream ss;
  if (s_verbosity >= Verb::v2)
    {
      ss << "  CPU binding                         : (" << listSet
               << ')' << std::endl;
    }
  s_thrOut << ss.rdbuf() << std::flush;
  hwloc_bitmap_free(set);
}

/*--------------------------------------------------------------------*/
//  Print thread memory bindings
/** \param[in] a_topology
 *                     The hardware topology
 *//*-----------------------------------------------------------------*/

void
ThreadSupport::printThreadMemoryBinding(const hwloc_topology_t a_topology) noexcept
{
  char listSet[c_strsz];
  hwloc_nodeset_t set = hwloc_bitmap_alloc();
  hwloc_membind_policy_t policy;

#if HWLOC_API_VERSION >= 0x00020100
  hwloc_get_membind(a_topology, set, &policy, HWLOC_MEMBIND_THREAD |
                                              HWLOC_MEMBIND_BYNODESET);
#else
  hwloc_get_membind_nodeset(a_topology, set, &policy, HWLOC_MEMBIND_THREAD);
#endif

  hwloc_bitmap_list_snprintf(listSet, sizeof(listSet), set);
  if (strlen(listSet) == 2 && listSet[1] == '-')
    {
      listSet[1] = '\0';
    }
  std::stringstream ss;
  if (s_verbosity >= Verb::v2)
    {
      ss << "  Memory binding (NUMA set)           : (" << listSet
               << ") " << membindPolicyName(policy) << std::endl;
    }
  s_thrOut << ss.rdbuf() << std::flush;
  hwloc_bitmap_free(set);
}

/*--------------------------------------------------------------------*/
//  Print thread memory bindings
/** \param[in]  a_pthr A pthread handle
 *//*-----------------------------------------------------------------*/

void
ThreadSupport::printPThreadAttr(pthread_t a_pthr) noexcept
{
  pthread_attr_t attr;
  pthread_getattr_np(a_pthr, &attr);

  // This code is sourced from the man page for pthread_getattr_np
  int s;
  size_t stacksize;
  int policy;
  struct sched_param schedparam;
  int detachstate;
  // int inheritsched;

  std::stringstream ss;

  s = pthread_attr_getstacksize(&attr, &stacksize);
  if (s_verbosity >= Verb::v2)
    {
      ss << "  Stack size                          : "
               << stacksize/1024 << " kiB" << "\n";
    }
//  printf("  Stack size            : %zd\n", stacksize);

  s = pthread_attr_getschedpolicy(&attr, &policy);
  if (s_verbosity >= Verb::v2)
    {
      ss << "  Scheduling policy                   : "
               << schedPolicyName(policy) << "\n";
    }

  // printf("  Scheduling policy     : %s\n",
  //        (policy == SCHED_FIFO) ? "SCHED_FIFO" :
  //        (policy == SCHED_RR) ? "SCHED_RR" :
  //        (policy == SCHED_OTHER) ? "SCHED_OTHER" : "[unknown]");

  s = pthread_attr_getschedparam(&attr, &schedparam);
  if (s_verbosity >= Verb::v2)
    {
      ss << "  Scheduling priority                 : "
               << schedparam.sched_priority << "\n";
    }
  // printf("  Scheduling priority   : %d\n", schedparam.sched_priority);

  s = pthread_attr_getdetachstate(&attr, &detachstate);
  const char* jstate =
    (detachstate == PTHREAD_CREATE_DETACHED) ? "DETACHED" :
    (detachstate == PTHREAD_CREATE_JOINABLE) ? "JOINABLE" : "???";
  if (s_verbosity >= Verb::v2)
    {
      ss << "  Detach state                        : "
               << jstate << "\n";
    }
  // printf("  Detach state          : %s\n",
  //        (detachstate == PTHREAD_CREATE_DETACHED) ? "DETACHED" :
  //        (detachstate == PTHREAD_CREATE_JOINABLE) ? "JOINABLE" :
  //        "???");

  // s = pthread_attr_getinheritsched(&attr, &inheritsched);
  // printf("  Inherit scheduler     : %s\n",
  //        (inheritsched == PTHREAD_INHERIT_SCHED) ? "INHERIT" :
  //        (inheritsched == PTHREAD_EXPLICIT_SCHED) ? "EXPLICIT" :
  //        "???");
  (void)s;
  s_thrOut << ss.rdbuf() << std::flush;
  pthread_attr_destroy(&attr);
}

/*--------------------------------------------------------------------*/
//  Test if thread priority can be set
/** return             Whether the thread priority be set
 *//*-----------------------------------------------------------------*/

bool ThreadSupport::testRLimit() noexcept
{
  struct rlimit limits;

  if (getrlimit(RLIMIT_STACK, &limits))
  {
    s_thrErr << "getrlimit(RLIMIT_STACK)" << std::endl;
  }

  if (getrlimit(RLIMIT_RTPRIO, &limits))
  {
    s_thrErr << "getrlimit(RLIMIT_RTPRIO)" << std::endl;
  }

  if (s_verbosity >= Verb::v1)
  {
    s_thrOut << "Current realtime priority limits      : soft "
            << limits.rlim_cur << " hard " << limits.rlim_max
            << "\n";
  }

  if (limits.rlim_max < 3)
  {
    if (s_verbosity >= Verb::vW)
    {
      s_thrOut << "Unable to modify priority limits.  Add line\n"
        "*                hard    rtprio          5\n"
        "to /etc/security/limits.conf and reboot.  Continuing without "
        "optimal priorities.\n";
      s_thrOut << "Will decrease priority of worker threads as an alternative.\n";
    }
    s_hasRTPrio = false;
  }
  else
  {
    s_hasRTPrio = true;
  }
  s_hasTestedRTPrio = true;
  return s_hasRTPrio;
}

/*--------------------------------------------------------------------*/
//  Set system realtime priorities
//*-------------------------------------------------------------------*/

void ThreadSupport::setRealtimePriorities()
{
  struct rlimit limits;

  std::stringstream ss;

  if (getrlimit(RLIMIT_STACK, &limits))
  {
    s_thrErr << "getrlimit(RLIMIT_STACK)" << std::endl;
  }

  if (getrlimit(RLIMIT_RTPRIO, &limits))
  {
    s_thrErr << "getrlimit(RLIMIT_RTPRIO)" << std::endl;
  }

  if (s_verbosity >= Verb::v1)
  {
    ss << "Current realtime priority limits      : soft "
             << limits.rlim_cur << " hard " << limits.rlim_max
             << "\n";
  }
  
  if (limits.rlim_max < 3)
  {
    if (s_verbosity >= Verb::vW)
    {
      ss << "Unable to modify priority limits.  Add line\n"
        "*                hard    rtprio          5\n"
        "to /etc/security/limits.conf and reboot.  Continuing without "
        "optimal priorities.\n";
      ss << "Will decrease priority of worker threads as an alternative.\n";
    }
    s_hasRTPrio = false;
  }
  else
  {
    if (limits.rlim_cur < 3)
      {
        limits.rlim_cur = 3;
      }

    if (setrlimit(RLIMIT_RTPRIO, &limits))
      {
        s_thrErr << "setrlimit(RLIMIT_RTPRIO)" << std::endl;
      }

    if (getrlimit(RLIMIT_RTPRIO, &limits))
      {
        s_thrErr << "getrlimit(RLIMIT_RTPRIO)" << std::endl;
      }

    if (s_verbosity >= Verb::v1)
    {
      ss << "New realtime priority limits          : soft "
                << limits.rlim_cur << " hard " << limits.rlim_max
                << "\n";
    }
  }
  if (s_verbosity >= Verb::v1)
  {
    ss << "\n";
  }
  s_thrOut << ss.rdbuf() << std::flush;
}

#include "NamespaceFooter.H"
