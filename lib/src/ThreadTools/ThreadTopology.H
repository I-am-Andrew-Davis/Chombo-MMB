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
 *  \file ThreadTopology.H
 *
 *  \brief Definition of class ThreadTopology
 *
 *//*+*************************************************************************/

#ifndef _THREADTOPOLOGY_H_
#define _THREADTOPOLOGY_H_

//----- Standard Library -----//

#include <memory>
#include <vector>
#include <map>

//----- Internal -----//

#include "ThreadSupport.H"

#include "NamespaceHeader.H"

namespace ThreadTools
{

/*******************************************************************************
 */
/// Describes sets and holds the topology of the threads on the hardware
/**
 *  Describes sets and holds the topology of the threads on the hardware.
 *  Must be constructed with a ThreadTopologyBuilder
 *
 ******************************************************************************/

  class ThreadTopology
  {
    friend class ThreadTopologyBuilder;


/*==============================================================================
 * Private constructors
 *============================================================================*/

private:

    /// Constructor with parameters
    ThreadTopology(unsigned numWorkers,
                   unsigned numCommThreads,
                   unsigned numGpuThreads);


/*==============================================================================
 * Public destructor
 *============================================================================*/

public:

    /// Destuctor
    ~ThreadTopology() = default;

/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

    /// Return the number of workers
    unsigned numWorkers() const;

    /// Return the number of process communication threads
    unsigned numProcComms() const;

    /// Return the number of GPU communication threads
    unsigned numGpuComms() const;

    /// Return the hardware object that a worker thread is binded to
    hwloc_obj_t getObjOfWorker(unsigned a_tid) const;

    ///
    hwloc_const_bitmap_t getBindingOfProcComm(unsigned a_tid) const;

    ///
    hwloc_const_bitmap_t getBindingOfGpuComm(unsigned a_tid) const;

    /// Return the compute set
    hwloc_const_bitmap_t getComputeSet() const;

    /// Return wether hyper threads are being used
    bool usingHyperThreads() const;
    
    /// Return if we can set realtime priority for threads
    bool hasRTPrio() const noexcept
      { return s_hasRTPrio; }


/*==============================================================================
 * Private member functions
 *============================================================================*/

private:

    /// Set the hardware objects 
    void setObjsOfWorkers(const std::vector<hwloc_obj_t>& a_objs);

    ///
    void setBindingsOfProcComms(
      const std::vector<hwloc_const_bitmap_t>& bindings);

    ///
    void setBindingsOfGpuComms(
      const std::vector<hwloc_const_bitmap_t>& bindings);


/*==============================================================================
 * Data members
 *============================================================================*/

protected:

    unsigned m_numWorkers {0};       ///< The number of worker threads
    unsigned m_numProcComms {0};     ///< The number of communicators to other
                                     ///< processes
    unsigned m_numGpuComms {0};      ///< The number of communicators to GPUs
    bool m_useHyperThreads;          ///< Are hyper threads used?
    std::vector<hwloc_obj_t> m_workerObjs;
                                     ///< Hardware objects the worker threads
                                     ///< are binded to
    std::vector<hwloc_const_bitmap_t> m_procCommBindings;
                                     ///< Hardware objects the proc comm
                                     ///< threads are binded to
    std::vector<hwloc_const_bitmap_t> m_gpuCommBindings;
                                     ///< Hardware objects the gpu comm
                                     ///< threads are binded to
   static bool s_hasRTPrio;          ///< T - can schedule threads with a
                                     ///<     realtime priority
    // m_coreset being used
    // m_cpuset being used
  };


/*******************************************************************************
 *
 * Class ThreadTopology: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Set the hardware objects that workers will bind to
/** \param[in]  a_objs  Hardware objects, PUs or cores
 *//*-----------------------------------------------------------------*/

  inline
  void ThreadTopology::setObjsOfWorkers(const std::vector<hwloc_obj_t>& a_objs)
  {
    m_workerObjs = a_objs;
  }

  inline
  void ThreadTopology::setBindingsOfProcComms(
    const std::vector<hwloc_const_bitmap_t>& a_bindings)
  {
    m_procCommBindings = a_bindings;
  }

  inline
  void ThreadTopology::setBindingsOfGpuComms(
    const std::vector<hwloc_const_bitmap_t>& a_bindings)
  {
    m_gpuCommBindings = a_bindings;
  }

/*--------------------------------------------------------------------*/
//  Return the number of worker threads
/** \return             The number of worker threads
 *//*-----------------------------------------------------------------*/

  inline
  unsigned ThreadTopology::numWorkers() const
  {
    return m_numWorkers;
  }

/*--------------------------------------------------------------------*/
//  Get hardware object that worker is binded to
/** \param[in]  a_id    Index of worker
 *  \return             hwloc object
 *//*-----------------------------------------------------------------*/

  inline
  hwloc_obj_t ThreadTopology::getObjOfWorker(unsigned a_id) const
  {
    return m_workerObjs[a_id];
  }

  inline
  hwloc_const_bitmap_t ThreadTopology::getBindingOfProcComm(unsigned a_id) const
  {
    return m_procCommBindings[a_id];
  }

  inline
  hwloc_const_bitmap_t ThreadTopology::getBindingOfGpuComm(unsigned a_id) const
  {
    return m_gpuCommBindings[a_id];
  }

/*--------------------------------------------------------------------*/
//  Return the number of process communication threads
/** \return             The number of process communication threads
 *//*-----------------------------------------------------------------*/

  inline
  unsigned ThreadTopology::numProcComms() const
  {
    return m_numProcComms;
  }

/*--------------------------------------------------------------------*/
//  Return the number of GPU communication threads
/** \return             The number of GPU communication threads
 *//*-----------------------------------------------------------------*/

  inline
  unsigned ThreadTopology::numGpuComms() const
  {
    return m_numGpuComms;
  }

/*--------------------------------------------------------------------*/
//  Return whether this is using hyper threads
/** \return             Whether hyper threads are turned on
 *//*-----------------------------------------------------------------*/

  inline
  bool ThreadTopology::usingHyperThreads() const
  {
    return m_useHyperThreads;
  }
}

#include "NamespaceFooter.H"

#endif  /* ! defined _THREADTOPOLOGY_H_ */
