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
 *  \file ThreadTeamSubTasks.H
 *
 *  \brief Definition of class ThreadTeamSubTasks
 *
 *//*+*************************************************************************/

#ifndef _THREADTEAMSUBTASKS_H_
#define _THREADTEAMSUBTASKS_H_

//----- Standard Library -----//

#include <vector>
#include <functional>
#include <atomic>
#include <condition_variable>

//----- Internal -----//

#include "ThreadConfig.H"
#include "ThreadBarrier.H"
#include "Profiler.H"

#ifdef CH_MPI
#include <mpi.h>
#endif

#include "NamespaceHeader.H"

namespace ThreadTools
{

  template<typename ...Args>
  class ThreadTeamSubTasks;

  /*******************************************************************************
   */
  /// The argument for the pthread function for ThreadTeamSubTasks
  /**
   *
   ******************************************************************************/

  template<typename ...Args>
  struct ThreadTeamSubTasksArg
  {
    std::shared_ptr<ThreadTeamSubTasks<Args...>> team;
  };


/*******************************************************************************
 */
/// A thread team that has sub tasks for each member
/**
 *  Templated to have sub tasks that have a number of arguments that can be
 *  different from each other, but each sub task has the same set of arguments
 *
 ******************************************************************************/

  template<typename ...Args>
  class ThreadTeamSubTasks
  {


/*==============================================================================
 * Aliases
 *============================================================================*/

public:

    /// The container to hold the arguments of one subtask
    using TupleArgs = std::tuple<Args...>;
    /// Alias to the signature of the function for a subtask
    using SubTaskFunc_t = void(Args...);
    /// Container for subtasks
    using SubTasks_t = std::vector<std::function<SubTaskFunc_t>>;
    /// Argument passed to thread of team
    using SubTasksArg_t = ThreadTeamSubTasksArg<Args...>;

/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

    /// Constructor
    ThreadTeamSubTasks(int a_teamSize)
      : m_teamSize(a_teamSize),
        m_threadBarrier(m_teamSize),
        m_threadBarrierPlusOne(m_teamSize+1)
      {
      };

    /// Destructor
    ~ThreadTeamSubTasks() = default;

//--Use synthesized copy, copy assignment, move, and move assignment


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

    /// Add a new task to the team
    void addNewTask(const std::vector<std::function<SubTaskFunc_t>>& a_subTasks,
                    const std::vector<TupleArgs>&                    a_args);

    /// Get the function each thread of the team executes for its duration
    static void* threadTeamFunction(void* args);

    /// Tell the team they are done and wait for them to finish
    void setDoneAndWait()
      {
        m_done = true;
        m_threadBarrierPlusOne.wait();
      }

    /// Wait until the current task is complete
    void wait() const
      {
        while (!m_taskComplete);
      }

    /// Has the team been told they are completely done
    bool isDone() const
      {
        return m_done;
      }

  protected:

    /// Execute a subtask
    template<size_t... Is>
    void execSubTask(const size_t a_idx, std::index_sequence<Is...>);

    /// Setup a profiler for the calling thread
    static void setupProfiler();

/*==============================================================================
 * Data members
 *============================================================================*/

protected:

    int m_teamSize;                  ///< The size of the team
    ThreadBarrier m_threadBarrier;   ///< A thread barrier with the size of the
                                     ///< team.
    ThreadBarrier m_threadBarrierPlusOne;
                                     ///< A thread barrier with the size of the
                                     ///< team plus one. Used
    std::atomic<bool> m_taskComplete{ false };
                                     ///< Is the current task complete?
    bool m_done{ false };            ///< Is the thread team completely done
                                     ///< with work and won't be assigned new
                                     ///< work?
    size_t m_numSubTasks;            ///< How many subtasks are in this task?
    SubTasks_t m_subTasks;           ///< The subtasks in the current task
    std::vector<TupleArgs> m_args;   ///< The arguments for each subtask
    std::atomic<size_t> m_idxOfNextTask;
                                     ///< The index of the next subtask that
                                     ///< has not been assigned
    static thread_local std::unique_ptr<Profiling::Profiler> tls_profiler;
                                     ///< static thread local profiler
  };

/*******************************************************************************
 *
 * Class ThreadTeamSubTasks: member definitions
 *
 ******************************************************************************/


  template<typename ...Args>
  thread_local std::unique_ptr<Profiling::Profiler> ThreadTeamSubTasks<Args...>::tls_profiler;
                                     ///< static thread local profiler

/*--------------------------------------------------------------------*/
//  Add a new task
/** \param[in]  a_subTasks  The subtasks of the new task
 *  \return     a_args      Each tuple is for one subtask
 *//*-----------------------------------------------------------------*/

  template<typename ...Args>
  void ThreadTeamSubTasks<Args...>::addNewTask(
    const std::vector<std::function<SubTaskFunc_t>>& a_subTasks,
    const std::vector<TupleArgs>&                    a_args)
  {
    m_threadBarrierPlusOne.wait();
    m_taskComplete = false;
    m_subTasks = a_subTasks;
    m_numSubTasks = m_subTasks.size();
    m_args = a_args;
    m_idxOfNextTask = 0;
    m_threadBarrierPlusOne.wait();
  }

/*--------------------------------------------------------------------*/
//  Execute the sub task
/** \param[in]  a_args  ThreadTeamSubTasksArg<Args...>*
 *  \return     void*   nullptr
 *//*-----------------------------------------------------------------*/

  template<typename ...Args>
  void* ThreadTeamSubTasks<Args...>::threadTeamFunction(void* a_args)
  {
    ThreadTeamSubTasksArg<Args...>* args =
      (ThreadTeamSubTasksArg<Args...>*)a_args;
    std::shared_ptr<ThreadTeamSubTasks<Args...>> team = args->team;

    if (c_useProfiler)
      {
        setupProfiler();
      }

    team->m_threadBarrierPlusOne.wait();
    while (!team->isDone())
      {
        team->m_threadBarrierPlusOne.wait();
        for (size_t idxOfTask = team->m_idxOfNextTask++;
             idxOfTask < team->m_numSubTasks;
             idxOfTask = team->m_idxOfNextTask++)
          {
            team->execSubTask(idxOfTask, std::index_sequence_for<Args...>{});
          }
        // TODO: optimize the next line with taskComplete
        team->m_threadBarrier.wait();
        team->m_taskComplete = true;
        team->m_threadBarrierPlusOne.wait();
      }
    return nullptr;
  }

/*--------------------------------------------------------------------*/
//  Execute the sub task
/** \param[in]  a_idx   The index of the subtask
 *  \param[in]  std::index_sequence<Is...>
 *                      sequence to expand args used in the subtask
 *//*-----------------------------------------------------------------*/

  template<typename ...Args>
  template<size_t... Is>
  void ThreadTeamSubTasks<Args...>::execSubTask(const size_t a_idx,
                                                std::index_sequence<Is...>)
  {
    auto& subtask = m_subTasks[a_idx];
    TupleArgs& arg = m_args[a_idx];

    if (c_useProfiler)
      {
        tls_profiler->startEvent(0, a_idx, "", 0);
      }

    subtask(std::get<Is>(arg)...);

    if (c_useProfiler)
      {
        tls_profiler->endEvent();
      }

  }

/*--------------------------------------------------------------------*/
//  Setup the profiler for the calling thread
/** \param[in]  a_idx   The index of the subtask
 *  \param[in]  std::index_sequence<Is...>
 *                      sequence to expand args used in the subtask
 *//*-----------------------------------------------------------------*/

  template<typename ...Args>
  void ThreadTeamSubTasks<Args...>::setupProfiler()
    {

      int procId = 0;
      std::string profilerName;
      {
        std::string basefilename = "profile";
        std::stringstream thrId_ss;
        std::string identifier;
        thrId_ss << std::setw( 4 ) << std::setfill( '0' ) << tls_threadId;
        std::stringstream pid_ss;
#ifdef CH_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &procId);
#endif
        pid_ss << std::setw( 4 ) << std::setfill( '0' ) << procId;
        identifier =  pid_ss.str() + "." + thrId_ss.str();
        profilerName = basefilename + "." + identifier + ".bin";
      }

      tls_profiler = std::make_unique<Profiling::Profiler>(profilerName,
                                                           procId,
                                                           tls_threadId,
                                                           ThreadArchitect::s_startTime,
                                                           1024);
    }

}

#include "NamespaceFooter.H"

#endif  /* ! defined _THREADTEAMSUBTASKS_H_ */
