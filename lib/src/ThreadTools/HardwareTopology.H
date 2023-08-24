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
 *  \file HardwareTopology.H
 *
 *  \brief Methods for getting the hardware topology
 *
 *//*+*************************************************************************/

#ifndef _HARDWARETOPOLOGY_H_
#define _HARDWARETOPOLOGY_H_

//----- Standard Library -----//

#include <cstdint>
#include <iostream>

//----- Third Party -----//

#include <hwloc.h>

//----- Internal -----//

#include "ThreadParameters.H"
#include "ThreadSupport.H"

#include "NamespaceHeader.H"

namespace ThreadTools
{

/*******************************************************************************
 */
/// Hardware Topology
/**
 *  Provides information on hardware and access to hardware objects
 *
 ******************************************************************************/

class HardwareTopology
{


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Constructor
  HardwareTopology();

  /// Destructor
  ~HardwareTopology();


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Report all hardware children of a hardware object
  void reportChildren(hwloc_obj_t a_obj, int a_depth) noexcept;

  /// Get CPU object
  hwloc_obj_t getCPUObj(unsigned a_idx) const;

  /// Get PU object
  hwloc_obj_t getPUObj(unsigned a_idx) const;

  /// Get core object
  hwloc_obj_t getCoreObj(unsigned a_idx) const;

  /// True if more than 1 NUMA node
  bool hasNUMA() const noexcept
  {
    return m_hasNUMA;
  }

  /// Return number of CPUs
  int numCPUs() const noexcept
  {
    return m_numCPU;
  }

  /// Return number of CPUs
  int numLocalCPUs() const noexcept
  {
    return m_numLocalCPU;
  }

  /// Return number of PUs
  int numPUs() const noexcept
  {
    return m_numPU;
  }

  /// Return number of cores
  int numCores() const noexcept
  {
    return m_numCore;
  }

  /// Return number of cores this process is bound to
  int numLocalCores() const noexcept
  {
    return m_numLocalCore;
  }

  /// Return number of cores this process is bound to
  int numLocalPUs() const noexcept
  {
    return m_numLocalPU;
  }

  /// Return whether the hardware has hyperthreads
  bool hasHyperThreads() const noexcept
  {
    return m_numLocalPU > m_numLocalCore;
  }

  /// Return the core set this process is bound to
  hwloc_const_bitmap_t coreSet() const noexcept
  {
    return m_coreset;
  }

  /// Return the cpu set
  hwloc_const_bitmap_t cpuSet() const noexcept
  {
    return m_cpuset;
  }

  /// Return the page size
  size_t getPageSize() const noexcept
  {
    return m_pageSize_B;
  }

  /// Return the hwloc_topology_t
  hwloc_topology_t topology() const
  {
    return m_topology;
  }

  /// Check alignment to page boundary
  bool isAlignedPage(void *const a_p) const noexcept
  {
    return ThreadSupport::isAligned(a_p, m_pageSize_B);
  }

  /// Grow a chunk of memory so it s a multiple of page size
  std::size_t alignToPage(std::size_t a_x) const noexcept
  {
    return ThreadSupport::alignTo(a_x, m_pageSize_B);
  }


/*==============================================================================
 * Private member functions
 *============================================================================*/

private:

  /// Queries the hardware to set member variables
  void queryHW();


/*==============================================================================
 * Data members
 *============================================================================*/

private:

  static constexpr int c_strsz = 128; ///< Size of character arrays used to
                                      ///< retrieve names from hwloc
  hwloc_topology_t m_topology;        ///< hwloc topology
  bool m_hasNUMA;                     ///< T - have more than 1 NUMA node
  int m_numCPU;                       ///< Number of CPUs (packages) in topo
  int m_numCore;                      ///< Number of cores in topology
  int m_numPU;                        ///< Number of PUs in topo
  int m_numLocalCPU;                  ///< Number of CPUs this process is
                                      ///< bound to
  int m_numLocalCore;                 ///< Number of cores this process is
                                      ///< bound to
  int m_numLocalPU;                   ///< Number of PUs this process is
                                      ///< bound to
  hwloc_cpuset_t m_cpuset;            ///< CPU set for this process
  hwloc_bitmap_t m_coreset;           ///< core set for this process
  int m_pageSize_B;                   ///< Page size
};

/*******************************************************************************
 *
 * Class ThreadBarrier: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Get PU object
/** \param[in]  a_idx   Index of PU (0 <= a_idx < m_numCPU)
 *  \return             hwloc PU (package) object
 *//*-----------------------------------------------------------------*/

inline hwloc_obj_t HardwareTopology::getPUObj(const unsigned a_idx) const
{
  TH_assert(a_idx < (unsigned)m_numPU);
  return hwloc_get_obj_by_type(m_topology, HWLOC_OBJ_PU, a_idx);
}

/*--------------------------------------------------------------------*/
//  Get core object
/** \param[in]  a_idx   Index of core (0 <= a_idx < m_numCore)
 *  \return             hwloc core object
 *//*-----------------------------------------------------------------*/

inline hwloc_obj_t HardwareTopology::getCoreObj(const unsigned a_idx) const
{
  TH_assert(a_idx < (unsigned)m_numCore);
  return hwloc_get_obj_by_type(m_topology, HWLOC_OBJ_CORE, a_idx);
}

} // namespace Threads

#include "NamespaceFooter.H"

#endif /* ! defined _HARDWARETOPOLOGY_H_ */
