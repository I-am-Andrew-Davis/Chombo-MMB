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
 * \file Profiler.H
 *
 * \brief Classes to create trace files
 *
 *//*+*************************************************************************/

#ifndef _PROFILER_H_
#define _PROFILER_H_

//----- Standard Library -----//

#include <chrono>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstring>
#include <vector>

//----- Internal -----//

#include "Event.H"
#include "ProfileHeader.H"

#include "NamespaceHeader.H"

namespace Profiling
{

/*******************************************************************************
 */
/// Create a trace of tasks
/**
*******************************************************************************/

class Profiler
{
  using steady_clock = std::chrono::steady_clock;
  using system_clock = std::chrono::system_clock;
  using time_point = std::chrono::time_point<steady_clock>;
  using duration = std::chrono::duration<double>;


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:
  
  /// Default constructor
  Profiler();

  /// Constructor
  Profiler(const std::string a_filename,
           const int a_processID,
           const int a_threadID,
           const time_point a_startTime,
           const int a_maxNumEvents);

  /// Destructor
  ~Profiler();

//--Use synthesized copy, copy assignment, move, and move assignment


/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Weak constructor
  void define(const std::string a_filename,
              const int a_processID,
              const int a_threadID,
              const time_point a_startTime,
              const int a_maxNumEvents);

  /// Allocate memory for the events
  void allocate();

  /// Deallocate memory for the events
  void deallocate();

  /// Write header to the ofstream
  void writeHeader();

  /// Start an event. Start time is at this call.
  void startEvent(const int a_eventtype,
                  const int a_eventID,
                  const std::string& a_eventName,
                  const int a_execNum);

  /// End the event
  void endEvent();

  /// Display the buffer for debugging
  void displayBuffer() const;

  /// Write the buffer to the ofstream
  void flush();

  /// Write the buffer and close the file
  void end();

  
/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  Event* m_events;                    ///< Array of recorded events
  int m_processID;                    ///< The process ID for this profiler
  int m_threadID;                     ///< The thread ID for this profiler
  time_point m_startTime;             ///< Universal start time.
                                      ///< Same across all profilers
  int m_maxNumEvents;                 ///< Max number of events in buffer 
                                      ///< before writing to ofstream
  int m_numEvents;                    ///< Current number of completed events 
                                      ///< in the buffer
  std::ofstream m_ofs;                ///< Binary file being written
  bool m_isOpenEvent;                 ///< Is there a started event 
                                      ///< that has not ended
};

/*--------------------------------------------------------------------*/
//  Start an event. Start time is at this call.
/** \param[in]  a_filename
 *  \param[in]  a_processID
 *  \param[in]  a_threadID
 *  \param[in]  a_startTime
 *  \param[in]  a_maxNumEvents
 *//*-----------------------------------------------------------------*/

inline
void Profiler::startEvent(const int          a_eventtype,
                          const int          a_eventID,
                          const std::string& a_eventName,
                          const int          a_execNum)
{
  assert(!m_isOpenEvent);
  m_isOpenEvent = true;
  m_events[m_numEvents].define(
    a_eventtype, 
    a_eventID, 
    m_processID, 
    m_threadID, 
    a_eventName, 
    a_execNum, 
    duration(steady_clock::now()-m_startTime).count());
}

/*--------------------------------------------------------------------*/
//  End the event
/**
 *//*-----------------------------------------------------------------*/

inline
void Profiler::endEvent()
{
  assert(m_isOpenEvent);
  m_events[m_numEvents].end(
    duration(steady_clock::now()-m_startTime).count());
  m_numEvents++;
  m_isOpenEvent = false;
  if (m_numEvents == m_maxNumEvents)
    {
      flush();
    }
}

/*--------------------------------------------------------------------*/
//  Write the buffer and close the file
/**
 *//*-----------------------------------------------------------------*/

inline
void Profiler::end()
{
  flush();
  m_ofs.close();
}

}

#define PROFILER_ADD(x,task)        \
  x->startEvent(task->attrType(),   \
                task->id(),         \
                task->name(),       \
                0);

#define PROFILER_POP(x)             \
  x->endEvent();

#include "NamespaceFooter.H"

#endif /* ! defined _PROFILER_H_ */
