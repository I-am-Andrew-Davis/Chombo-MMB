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
 *  \file ThreadTeamArchitect.H
 *
 *  \brief Definition of ThreadTeamArchitect
 *
 *//*+*************************************************************************/

#ifndef _THREADTEAMARCHITECT_H_
#define _THREADTEAMARCHITECT_H_

//----- Standard Library -----//

//----- Internal -----//

#include "ThreadArchitect.H"
#include "ThreadTeam.H"

#include "NamespaceHeader.H"

/*******************************************************************************
 */
/// Namespace to query hardware and bind threads
/**
 *  
 *
 ******************************************************************************/

namespace ThreadTools
{

/*******************************************************************************
 */
/// Architect for thread teams
/**
 *  Builds thread teams and binds threads to hardware
 *
 ******************************************************************************/

class ThreadTeamArchitect : public ThreadArchitect
{


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Constructor
  ThreadTeamArchitect(std::shared_ptr<const HardwareTopology> hwTopology = nullptr);

  /// Destructor
  virtual ~ThreadTeamArchitect() = default;


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Setup the thread topology
  void setupTopology(int numThreads,
                     int threadsPerTeam,
                     int numCommThreads,
                     int numGpuThreads,
                     size_t extraStackSize_MiB,
                     bool useHyperThreading,
                     bool bindThread = true);

  /// Return the number of thread teams
  unsigned numThreadTeams() const;

  /// Return the number of threads per team
  unsigned numThreadsPerTeam() const;

  /// Return a const reference to the team according to the team id
  const ThreadTeam& team(unsigned teamID) const;

  /// Return a reference to the team according to the team id
  ThreadTeam& team(unsigned teamID);

/*==============================================================================
 * Protected member functions
 *============================================================================*/

protected:

  /// Set the threads per team
  void setThreadsPerTeam(int  a_threadsPerTeam,
                         bool a_useHyperThreads);

  /// Set the threads per team when using hyper threads
  /// (Called by setThreadsPerTeam)
  void setThreadsPerTeamWithHyperThreads(int a_threadsPerTeam);

  /// Set the threads per team when using hyper threads
  /// (Called by setThreadsPerTeam)
  void setThreadsPerTeamNoHyperThreads(int a_threadsPerTeam);


/*==============================================================================
 * Data members
 *============================================================================*/

public:


protected:

  int m_threadsPerTeam {0};           ///< The number of threads per team
  std::vector<ThreadTeam> m_teams;    ///< The thread teams
};


/*******************************************************************************
 *
 * Class ThreadTeamArchitect: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Get a const reference to a thread team according to its team ID
/** \param[in]  a_teamID
 *                      The team ID
 *  \return             The thread team
 *//*-----------------------------------------------------------------*/

  inline
  const ThreadTeam& ThreadTeamArchitect::team(unsigned a_teamID) const
  {
    return m_teams[a_teamID];
  }

/*--------------------------------------------------------------------*/
//  Get a reference to a thread team according to its team ID
/** \param[in]  a_teamID
 *                      The team ID
 *  \return             The thread team
 *//*-----------------------------------------------------------------*/

  inline
  ThreadTeam& ThreadTeamArchitect::team(unsigned a_teamID)
  {
    return m_teams[a_teamID];
  }

/*--------------------------------------------------------------------*/
//  Get the number of thread teams
/** \return             The number of thread teams
 *//*-----------------------------------------------------------------*/

  inline
  unsigned ThreadTeamArchitect::numThreadTeams() const
  {
    return m_threadTopology ? m_threadTopology->numWorkers() / m_threadsPerTeam : 0;
  }

/*--------------------------------------------------------------------*/
//  Get the number of threads per team
/** \return             The number of threads per team
 *//*-----------------------------------------------------------------*/

  inline
  unsigned ThreadTeamArchitect::numThreadsPerTeam() const
  {
    return m_threadsPerTeam;
  }

}

#include "NamespaceFooter.H"

#endif  /* ! defined _THREADTEAMARCHITECT_H_ */
