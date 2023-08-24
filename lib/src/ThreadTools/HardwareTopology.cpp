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
 *  \file Topology.cpp
 *
 *  \brief Non-inline definitions for classes in Topology.H
 *
 *//*+*************************************************************************/

//----- Standard Library -----//

#include <cstdio>

//----- System -----//

#include <pthread.h>
#include <sys/mman.h>
#include <sys/resource.h>
#include <unistd.h>

//----- Internal -----//

#include "HardwareTopology.H"

#include "NamespaceHeader.H"

//--Compatible with older versions of hwloc

#if HWLOC_API_VERSION < 0x00010b00
#define HWLOC_OBJ_NUMANODE HWLOC_OBJ_NODE
#define HWLOC_OBJ_PACKAGE HWLOC_OBJ_SOCKET
#endif


using namespace ThreadTools;

/*******************************************************************************
 *
 * Class Topology: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** \param[in]  a_extraStackSize_MiB
 *                      Additional stacksize to soft resource limit
 *                      (or 8192 kiB if infinite).  At least 1 page
 *                      is added to the stack memory for use as
 *                      a guard page against stack overflow.
 *//*-----------------------------------------------------------------*/

HardwareTopology::HardwareTopology()
{
  queryHW();
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

HardwareTopology::~HardwareTopology()
{
  hwloc_bitmap_free(m_cpuset);
  hwloc_bitmap_free(m_coreset);
  hwloc_topology_destroy(m_topology);
}

/*==============================================================================
 * Public member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Get CPU object
/** \param[in]  a_idx   Index of CPU (0 <= a_idx < m_numCPU)
 *  \return             hwloc CPU (package) object
 *//*-----------------------------------------------------------------*/

hwloc_obj_t HardwareTopology::getCPUObj(const unsigned a_idx) const
{
  // This function is not inlined, so that hwloc backward compatibility
  // defines can be kept in the cpp file
  TH_assert(a_idx < (unsigned)m_numCPU);
  return hwloc_get_obj_by_type(m_topology, HWLOC_OBJ_PACKAGE, a_idx);
}

/*--------------------------------------------------------------------*/
//  Query the hw to set member variables
/**
 *//*-----------------------------------------------------------------*/

void HardwareTopology::queryHW()
{
  if (s_verbosity >= Verb::v4)
    {
      s_thrOut << "Using hwloc api version " << ThreadSupport::hwlocVersion()
               << " to query hardware.\n\n";
    }

  hwloc_topology_init(&m_topology);
  hwloc_topology_load(m_topology);

  if (s_verbosity >= Verb::v2)
    {
      //    ___012345678901234567890123456789012345678901234567890123456789
      s_thrOut
        << "              Hardware           CPU Set             NUMA Set\n"
        << "   ------------------------------------------------------------"
        << std::endl;
      reportChildren(hwloc_get_root_obj(m_topology), 0);
      s_thrOut << std::endl;
    }
  // Check for NUMA
  const int depthNUMA = hwloc_get_type_depth(m_topology, HWLOC_OBJ_NUMANODE);
  m_hasNUMA = (depthNUMA != HWLOC_TYPE_DEPTH_UNKNOWN);
  // Number of CPUs
  const int depthCPU = hwloc_get_type_depth(m_topology, HWLOC_OBJ_PACKAGE);
  if (depthCPU == HWLOC_TYPE_DEPTH_UNKNOWN)
    {
      if (s_verbosity >= Verb::vE)
        {
          s_thrOut << "The number of CPUs is unknown" << std::endl;
        }
    }
  m_numCPU = hwloc_get_nbobjs_by_depth(m_topology, depthCPU);
  if (s_verbosity >= Verb::v1)
    {
      s_thrOut << "Number of CPUs                        : " << m_numCPU
               << std::endl;
    }
  // Number of cores
  const int depthCore =
    hwloc_get_type_or_below_depth(m_topology, HWLOC_OBJ_CORE);
  if (depthCore == HWLOC_TYPE_DEPTH_UNKNOWN)
    {
      s_thrOut << "The number of Cores is unknown" << std::endl;
    }
  m_numCore = hwloc_get_nbobjs_by_depth(m_topology, depthCore);
  if (s_verbosity >= Verb::v1)
    {
      s_thrOut << "Number of cores                       : " << m_numCore
               << std::endl;
    }
  // Number of threads
  const int depthPU = hwloc_get_type_or_below_depth(m_topology, HWLOC_OBJ_PU);
  if (depthPU == HWLOC_TYPE_DEPTH_UNKNOWN)
    {
      s_thrOut << "The number of PUs is unknown" << std::endl;
    }
  m_numPU = hwloc_get_nbobjs_by_depth(m_topology, depthPU);
  if (s_verbosity >= Verb::v1)
    {
      s_thrOut << "Number of PUs                         : " << m_numPU
               << std::endl;
    }
  // HW this process is bounded to
  m_cpuset = hwloc_bitmap_alloc();
  m_coreset = hwloc_bitmap_alloc();
  hwloc_get_cpubind(m_topology, m_cpuset, HWLOC_CPUBIND_PROCESS);
  m_numLocalCore = hwloc_get_nbobjs_inside_cpuset_by_type(m_topology,
                                                          m_cpuset,
                                                          HWLOC_OBJ_CORE);
  if (s_verbosity >= Verb::v1)
    {
      s_thrOut << "Number of cores for this process      : " << m_numLocalCore
               << std::endl;
    }

  m_numLocalCPU = hwloc_get_nbobjs_inside_cpuset_by_type(m_topology,
                                                         m_cpuset,
                                                         HWLOC_OBJ_PACKAGE);

  m_numLocalPU = hwloc_get_nbobjs_inside_cpuset_by_type(m_topology,
                                                        m_cpuset,
                                                        HWLOC_OBJ_PU);

  {
    char listSet[c_strsz];
    hwloc_bitmap_list_snprintf(listSet, sizeof(listSet), m_cpuset);
    if (s_verbosity >= Verb::v1)
      {
        s_thrOut << "CPU set for this process              : (" << listSet
                 << ')' << std::endl;
      }
  }
  // Add each core this process is bound to to m_coreset
  hwloc_obj_t core = nullptr;
  for (int i = 0; i != m_numLocalCore; ++i)
    {
      core = hwloc_get_next_obj_inside_cpuset_by_type(m_topology,
                                                      m_cpuset,
                                                      HWLOC_OBJ_CORE,
                                                      core);
      char listSet[c_strsz];
      hwloc_bitmap_list_snprintf(listSet, sizeof(listSet), core->cpuset);
      char objType[c_strsz];
      hwloc_obj_type_snprintf(objType, sizeof(objType), core, 0);

      if (s_verbosity >= Verb::v2)
        {
          s_thrOut << objType << core->os_index << "-" << core->logical_index
                   << " cpuset                        : (" << listSet << ')'
                   << std::endl;
        }
      hwloc_bitmap_set(m_coreset, core->logical_index);
    }

  {
    char listSet[c_strsz];
    hwloc_bitmap_list_snprintf(listSet, sizeof(listSet), m_coreset);
    if (s_verbosity >= Verb::v1)
      {
        s_thrOut << "Core set for this process is          : (" << listSet
                 << ")\n";
      }
  }

  //--Get the page size

  auto pageSize_B = sysconf(_SC_PAGESIZE);

  if (pageSize_B <= 0)
    {
      s_thrOut << "Error(" << errno << ") while getting the page size"
               << std::endl;
    }

  m_pageSize_B = pageSize_B;
  if (s_verbosity >= Verb::v1)
    {
      s_thrOut << "Page size                             : " << m_pageSize_B
               << " B" << std::endl;
    }
}

/*--------------------------------------------------------------------*/
//  Report all hardware children of a hardware object
/** \param[in]  a_obj   hwloc hardware object
 *  \param[in]  a_depth Current depth in hardware levels
 *//*-----------------------------------------------------------------*/

void HardwareTopology::reportChildren(hwloc_obj_t a_obj,
                                      const int a_depth) noexcept
{
  char objType[c_strsz];
  char attrInfo[c_strsz];
  char listSet[c_strsz];
  constexpr int wsz = 2 * c_strsz;
  char wbuffer[wsz];
  int wc = 0;

  hwloc_obj_type_snprintf(objType, sizeof(objType), a_obj, 0);
  wc += snprintf(wbuffer + wc, wsz - wc, "   %*s%s", 2 * a_depth, "", objType);
  if (a_obj->os_index != (unsigned)-1)
    {
      wc += snprintf(wbuffer + wc, wsz - wc, "%u", a_obj->os_index);
    }
  hwloc_obj_attr_snprintf(attrInfo, sizeof(attrInfo), a_obj, " ", 0);
  if (*attrInfo)
    {
      wc += snprintf(wbuffer + wc, wsz - wc, "(%s)", attrInfo);
    }
  while (wc < 30)
    {
      wbuffer[wc++] = ' ';
    }
  hwloc_bitmap_list_snprintf(listSet, sizeof(listSet), a_obj->cpuset);
  wc += snprintf(wbuffer + wc, wsz - wc, "(%s)", listSet);
  while (wc < 50)
    {
      wbuffer[wc++] = ' ';
    }
  hwloc_bitmap_list_snprintf(listSet, sizeof(listSet), a_obj->nodeset);
  if (strlen(listSet) == 2 && listSet[1] == '-')
    {
      listSet[1] = '\0';
    }
  wc += snprintf(wbuffer + wc, wsz - wc, "(%s)", listSet);
  s_thrOut << wbuffer << std::endl;
  for (unsigned i = 0; i < a_obj->arity; ++i)
    {
      reportChildren(a_obj->children[i], a_depth + 1);
    }
}

#include "NamespaceFooter.H"
