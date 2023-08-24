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
 *  \file ThreadTeamSubTaskArchitect.H
 *
 *  \brief Definition of class ThreadTeamSubTaskArchitect
 *
 *//*+*************************************************************************/

#ifndef _THREADTEAMSUBTASKARCHITECT_H_
#define _THREADTEAMSUBTASKARCHITECT_H_

//----- Standard Library -----//

//----- Internal -----//

#include "ThreadTeamArchitect.H"
#include "ThreadTeamSubTasks.H"

#include "NamespaceHeader.H"

namespace ThreadTools
{

/*******************************************************************************
 */
/// Architect for thread teams for Chord
/**
 *  Builds thread teams and binds threads to hardware.  A task is added as a
 *  collection of subtasks. All the subtasks have the same signature for their
 *  function. The next available threadteam gets the next task. The threads
 *  of the team iterate through the subtasks executing them. When all the
 *  subtasks of a task have been completed, the threadteam is available for
 *  another task.
 *
 ******************************************************************************/

template<typename T>
class ThreadTeamSubTaskArchitect : public ThreadTeamArchitect
{


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Constructor
  ThreadTeamSubTaskArchitect(
    std::shared_ptr<const HardwareTopology> a_hwTopology = nullptr)
  : ThreadTeamArchitect(a_hwTopology)
  {}

  /// Destructor
  virtual ~ThreadTeamSubTaskArchitect() = default;

//--Use synthesized copy, copy assignment, move, and move assignment


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Add a new team
  void addTeamSubTasks(std::shared_ptr<T> a_teamSubTasks);

  /// Add a new task
  void addNewTask(
    const std::vector<std::function<typename T::SubTaskFunc_t>> a_subTasks,
    const std::vector<typename T::TupleArgs>                    a_args);

  /// Wait for all the tasks to be completed
  void wait();

  /// Tell each team they are done
  void setToDone();

  /// Wait for the finish of all execution
  void waitToFinish() override;

  /// Return the team according to a team ID
  std::shared_ptr<T> teamSubTasks(unsigned teamID) const;

/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  std::vector<std::shared_ptr<T>> m_teamSubTasks;
                                     ///< The thread teams
  size_t m_nextTeamToGetTask{0};     ///< The next team that will pick up the
                                     ///< next task
  std::mutex m_mtx;                  ///< Lock used while testing if
                                     ///< m_nextTeamToGetTask is at the end of
                                     ///< the teams and should be reset to 0.
};


/*******************************************************************************
 *
 * Class ThreadTeamArchitectChord: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Get the subtasks for a team according to a team ID
/** \param[in]  a_teamID
 *                      The subtasks
 *//*-----------------------------------------------------------------*/

  template<typename T>
  std::shared_ptr<T> ThreadTeamSubTaskArchitect<T>::teamSubTasks(
    unsigned a_teamID) const
  {
    return m_teamSubTasks[a_teamID];
  }


/*--------------------------------------------------------------------*/
//  Add a new team
/** \param[in]  a_team
 *                      The new team to add
 *//*-----------------------------------------------------------------*/

  template<typename T>
  void ThreadTeamSubTaskArchitect<T>::addTeamSubTasks(
    std::shared_ptr<T> a_team)
  {
    TH_assert(a_team);
    m_teamSubTasks.push_back(a_team);
  }

/*--------------------------------------------------------------------*/
//  Add a new task
/** \param[in]  a_subTasks
 *                      Functions for each subtask
 *  \param[in]  a_args  Arguments for each subtask function
 *//*-----------------------------------------------------------------*/

  template<typename T>
  void ThreadTeamSubTaskArchitect<T>::addNewTask(
    const std::vector<std::function<typename T::SubTaskFunc_t>> a_subTasks,
    const std::vector<typename T::TupleArgs>                    a_args)
  {
    TH_assert(a_subTasks.size() == a_args.size());

    std::unique_lock<std::mutex> lock(m_mtx);
    auto nextTeam = m_nextTeamToGetTask++;
    if (nextTeam == m_teamSubTasks.size())
      {
        nextTeam = 0;
        m_nextTeamToGetTask = 1;
      }
    lock.unlock();
    teamSubTasks(nextTeam)->addNewTask(a_subTasks, a_args);
  }

/*--------------------------------------------------------------------*/
//  Wait for all the tasks to be completedoverride-----------------*/

  template<typename T>
  void ThreadTeamSubTaskArchitect<T>::wait()
    {
      for (auto& team : m_teamSubTasks)
        {
          team->wait();
        }
    }

/*--------------------------------------------------------------------*/
//  Tell each team they are done
/** 
 *//*-----------------------------------------------------------------*/

  template<typename T>
  void ThreadTeamSubTaskArchitect<T>::setToDone()
    {
      for (auto& team : m_teamSubTasks)
        {
          team->setDoneAndWait();
        }
    }

/*--------------------------------------------------------------------*/
//  Wait for the finish of all execution
/** 
 *//*-----------------------------------------------------------------*/

  template<typename T>
  void ThreadTeamSubTaskArchitect<T>::waitToFinish()
  {
    setToDone();
    ThreadArchitect::waitToFinish();
  }

}

#include "NamespaceFooter.H"

#endif  /* ! defined _THREADTEAMSUBTASKARCHITECT_H_ */
