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
 *  \file ThreadTopologyBuilder.cpp
 *
 *  \brief Non-inline definitions for classes in ThreadTopologyBuilder.H
 *
 *//*+*************************************************************************/


//----- Standard Library -----//

//----- Internal -----//

#include "ThreadTopologyBuilder.H"
#include "ThreadParameters.H"
#include "ThreadTopology.H"

#include "NamespaceHeader.H"

using namespace ThreadTools;

/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Construct ThreadTopologyBuilder
/** \param[in]  a_hwTopology
 *                      The hardware topology to build on
 *//*-----------------------------------------------------------------*/

ThreadTopologyBuilder::ThreadTopologyBuilder(
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
//  Build a ThreadTopology
/** \param[in]  a_numWorkers
 *                      The number of worker threads. Can not be 0.
 *                      If less than zero, the number of worker
 *                      threads will equal the number of cores this
 *                      process is bound to.
 *  \param[in]  a_useHyperThreads
 *                      Wether to use hyperthreads if present.
 *  \param[in]  a_numCommThread
 *                      The number of process communication threads to use
 *  \param[in]  a_numGPUThread
 *                      The number of GPU communication threads to use
 *//*-----------------------------------------------------------------*/

ThreadTopology ThreadTopologyBuilder::build(const int a_numWorkers,
                                            const bool a_useHyperThreads,
                                            const int a_numCommThread,
                                            const int a_numGPUThread)
{
  m_useHyperThreads = a_useHyperThreads & m_hwTopology->hasHyperThreads();
  assert(a_numWorkers != 0 &&
         "The number of worker threads can not be 0. If less\n"
         " than zero, the number of worker threads will equal\n"
         " the number of cores this process is bound to.");

  unsigned numWorkers = setNumberOfWorkers(a_numWorkers);
  unsigned numCommThread = setNumberOfComms(a_numCommThread);
  unsigned numGPUThread = setNumberOfGpuComms(a_numGPUThread);

  std::vector<hwloc_obj_t> workerObjs = getComputeObjs();
  assert(numWorkers <= workerObjs.size());
  std::vector<hwloc_const_bitmap_t> procCommBindings =
    getCommBindings(numCommThread);
  std::vector<hwloc_const_bitmap_t> gpuCommBindings =
    getGpuBindings(numGPUThread);

  ThreadTopology threadTopology(numWorkers, numCommThread, numGPUThread);
  threadTopology.m_useHyperThreads = m_useHyperThreads;
  threadTopology.setObjsOfWorkers(workerObjs);
  threadTopology.setBindingsOfProcComms(procCommBindings);
  threadTopology.setBindingsOfGpuComms(gpuCommBindings);
  return threadTopology;
}

/*--------------------------------------------------------------------*/
//  Set the number of workers
/** \param[in]  a_numWorkers
 *                      The number of worker threads. Can not be 0.
 *                      If less than zero, the number of worker
 *                      threads will equal the number of cores this
 *                      process is bound to or the number of threads
 *                      if using hyperthreads.
 *//*-----------------------------------------------------------------*/

unsigned ThreadTopologyBuilder::setNumberOfWorkers(const int a_numWorkers) const
{
  unsigned numWorkers;
  if (a_numWorkers < 0)
    {
      if (m_useHyperThreads)
        {
          numWorkers = m_hwTopology->numLocalPUs();
        }
      else
        {
          numWorkers = m_hwTopology->numLocalCores();
        }
    }
  else
    {
      numWorkers = a_numWorkers;
      TH_assert(
        !(m_useHyperThreads && a_numWorkers > m_hwTopology->numLocalPUs()) &&
        "The argument for the number of worker threads"
        "exceeds the number of local PUs.");
      TH_assert(
        !(!m_useHyperThreads && a_numWorkers > m_hwTopology->numLocalCores()) &&
        "The argument for the number of worker threads"
        "exceeds the number of local cores.");
    }
  return numWorkers;
}

/*--------------------------------------------------------------------*/
//  Return the number of communication threads based on an input
/** \param[in]  a_numCommThread
 *                      The arg for the number of communication
 *                      threads. If less than 0, the number is provided.
 *  \return             The number of communication threads.
 *//*-----------------------------------------------------------------*/

unsigned
ThreadTopologyBuilder::setNumberOfComms(const int a_numCommThread) const
{
  return a_numCommThread > 0 ? a_numCommThread : 0;
}

/*--------------------------------------------------------------------*/
//  Return the number of GPU communication threads based on an input
/** \param[in]  a_numGPUThread
 *                      The arg for the number of GPU communication
 *                      threads. If less than 0, the number is provided.
 *  \return             The number of GPU communication threads.
 *//*-----------------------------------------------------------------*/

unsigned
ThreadTopologyBuilder::setNumberOfGpuComms(const int a_numGPUThread) const
{
  return a_numGPUThread > 0 ? a_numGPUThread : 0;
}

/*--------------------------------------------------------------------*/
//  Return the indices in the compute set for worker threads
/** \return             The indices of the compute set.
 *//*-----------------------------------------------------------------*/

std::vector<unsigned> ThreadTopologyBuilder::getIndices() const
{
  std::vector<unsigned> indices;
  unsigned idx;
  hwloc_const_bitmap_t set = getComputeSet();

  hwloc_bitmap_foreach_begin(idx, set)
  {
    indices.push_back(idx);
  }
  hwloc_bitmap_foreach_end();

  return indices;
}

/*--------------------------------------------------------------------*/
//  Return the hardware objects in the compute set for worker threads
/** \return             The hardware objects of the compute set.
 *//*-----------------------------------------------------------------*/

std::vector<hwloc_obj_t> ThreadTopologyBuilder::getComputeObjs() const
{
  std::vector<hwloc_obj_t> objs;
  unsigned idx;
  hwloc_const_bitmap_t set = getComputeSet();

  hwloc_bitmap_foreach_begin(idx, set)
  {
    hwloc_obj_t obj = getComputeObj(idx);
    objs.push_back(obj);
  }
  hwloc_bitmap_foreach_end();

  return objs;
}

/*--------------------------------------------------------------------*/
//  DO NOT USE. Get the bindings to use for each comm thread.
/** \param[in]  a_numCommThreads
 *                      The number of communication threads
 *  \return             The bindings to use for each comm thread
 *//*-----------------------------------------------------------------*/

std::vector<hwloc_const_bitmap_t>
ThreadTopologyBuilder::getCommBindings(const unsigned a_numCommThreads) const
{
  std::vector<hwloc_const_bitmap_t> bindings;
  for (unsigned ctid = 0; ctid < a_numCommThreads; ++ctid)
    {
      bindings.push_back(m_hwTopology->cpuSet());
    }
  return bindings;
}

/*--------------------------------------------------------------------*/
//  DO NOT USE. Get the bindings to use for each GPU comm thread.
/** \param[in]  a_numGpuThreads
 *                      The number of communication threads
 *  \return             The bindings to use for each GPU comm thread
 *//*-----------------------------------------------------------------*/

std::vector<hwloc_const_bitmap_t>
ThreadTopologyBuilder::getGpuBindings(const unsigned a_numGpuThreads) const
{
  std::vector<hwloc_const_bitmap_t> bindings;
  for (unsigned ctid = 0; ctid < a_numGpuThreads; ++ctid)
    {
      bindings.push_back(m_hwTopology->cpuSet());
    }
  return bindings;
}

#include "NamespaceFooter.H"
