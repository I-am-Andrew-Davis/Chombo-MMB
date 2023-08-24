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
 *  \file ThreadTopologyBuilder.H
 *
 *  \brief Definition of class ThreadTopologyBuilder
 *
 *//*+*************************************************************************/

#ifndef _THREADTOPOLOGYBUILDER_H_
#define _THREADTOPOLOGYBUILDER_H_

//----- Standard Library -----//

#include <memory>
#include <vector>
#include <map>

//----- Internal -----//

#include "ThreadSupport.H"
#include "HardwareTopology.H"

#include "NamespaceHeader.H"

namespace ThreadTools
{

  //-- Forward declarations
  class ThreadTopology;

/*******************************************************************************
 */
/// Builds ThreadTopology
/**
 *  From parameters and the HardwareTopology, builds a ThreadTopology
 *
 ******************************************************************************/

  class ThreadTopologyBuilder
  {


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

    public:

    /// Constructor
    ThreadTopologyBuilder(std::shared_ptr<const HardwareTopology> hwTopology);

    /// Destuctor
    ~ThreadTopologyBuilder() = default;


/*==============================================================================
 * Public member functions
 *============================================================================*/

  public:

    /// Build function returns a ThreadTopology
    ThreadTopology build(int numWorkerThreads,
                         bool useHyperThreads,
                         int useCommThread = 0,
                         int useGPUThread = 0);


/*==============================================================================
 * Private member functions
 *============================================================================*/

  private:

    /// Get the compute object by index. Core if not using hyper threads.
    /// Thread if using hyper threads.
    hwloc_obj_t getComputeObj(unsigned a_idxObj) const;

    /// Get the set of compute objects Cores if not using hyper threads.
    /// Threads if using hyper threads.
    hwloc_const_bitmap_t getComputeSet() const;

    /// Set the number of worker threads
    unsigned setNumberOfWorkers(int numWorkers) const;

    /// Set the number of process communication threads
    unsigned setNumberOfComms(int numCommThread) const;

    /// Set the number of GPU communication threads
    unsigned setNumberOfGpuComms(int numGPUThread) const;

    /// Return a list of indices in the compute set
    std::vector<unsigned> getIndices() const;

    /// Return the compute objects for worker threads
    std::vector<hwloc_obj_t> getComputeObjs() const;

    /// DO NOT USE. Get the bindings to use for each comm thread.
    std::vector<hwloc_const_bitmap_t>
      getCommBindings(unsigned numCommThread) const;

    /// DO NOT USE. Get the bindings to use for each GPU comm thread.
    std::vector<hwloc_const_bitmap_t>
      getGpuBindings(unsigned numGPUThread) const;


/*==============================================================================
 * Data members
 *============================================================================*/

private:

    std::shared_ptr<const HardwareTopology> m_hwTopology;
                                     ///< Used to get hardware resources
    bool m_useHyperThreads;          ///< Flag to use hyper threads
  };


/*******************************************************************************
 *
 * Class ThreadTopologyBuilder: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Get the compute object by index. Returns a PU if using hyper
//  threads. Returns a core if not using hyper threads.
/** \param[in]  a_idxObj Index of compute object
 *  \return              The compute object
 *//*-----------------------------------------------------------------*/

  inline
  hwloc_obj_t ThreadTopologyBuilder::getComputeObj(unsigned a_idxObj) const
    {
      hwloc_obj_t obj = m_useHyperThreads ?
        m_hwTopology->getPUObj(a_idxObj) :
        m_hwTopology->getCoreObj(a_idxObj);
      return obj;
    }

/*--------------------------------------------------------------------*/
//  Get the compute set depending on if hyperthreads are being used.
/** \return              The compute set
 *//*-----------------------------------------------------------------*/

  inline
  hwloc_const_bitmap_t ThreadTopologyBuilder::getComputeSet() const
  {
    return m_useHyperThreads ? m_hwTopology->cpuSet() : m_hwTopology->coreSet();
  }

}

#include "NamespaceFooter.H"

#endif  /* ! defined _THREADTOPOLOGYBUILDER_H_ */
