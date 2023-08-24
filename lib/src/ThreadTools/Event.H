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
 * \file Event.H
 *
 * \brief Classes to create trace files
 *
 *//*+*************************************************************************/

#ifndef _EVENT_H_
#define _EVENT_H_

//----- Standard Library -----//

#include <chrono>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstring>

#include "NamespaceHeader.H"

namespace Profiling
{

/*******************************************************************************
 */
/// An event is a record of what the processor or thread is doing, e.g. a task
/**
*******************************************************************************/

class Event
{
  
public:

  Event()
    :
    m_type(0),
    m_id(-1)
    {
    }
  
  Event(const int a_eventtype, 
        const int a_eventID, 
        const int a_pid,
        const int a_tid,
        const std::string& a_eventName, 
        const int a_execNum, 
        const double a_startTime)
    :
    m_type(a_eventtype),
    m_id(a_eventID),
    m_pid(a_pid),
    m_tid(a_tid),
    m_execNum(a_execNum),
    m_startTime(a_startTime)
    {
      memcpy(m_name, a_eventName.c_str(), a_eventName.size()+1);
    }

  void define(const int a_eventtype, 
              const int a_eventID, 
              const int a_pid,
              const int a_tid,
              const std::string& a_eventName, 
              const int a_execNum, 
              const double a_startTime)
  {
    m_type       = a_eventtype;
    m_id         = a_eventID;
    m_pid        = a_pid;
    m_tid        = a_tid;
    m_execNum    = a_execNum;
    m_startTime  = a_startTime;
    memcpy(m_name, a_eventName.c_str(), a_eventName.size()+1);
  }

  void end(double a_endTime)
    {
      m_endTime = a_endTime;
    }

  int id() const
    {
      return m_id;
    }
  
  int execNum() const
    {
      return m_execNum;
    }

  const char* name() const
     {
       return m_name;
     }
  
  int procID() const
    {
      return m_pid;
    }
  
  int threadID() const
    {
      return m_tid;
    }

  double startTime() const
    {
      return m_startTime;
    }

  double endTime() const
    {
      return m_endTime;
    }

  friend std::ostream& operator<<(std::ostream& os, const Event a_event)
    {
      os << "{\"name\": \"" << a_event.m_name
         << "\", \"id\": " << a_event.m_id 
         << ", \"procID\":" << a_event.m_pid
         << ", \"threadID\":" << a_event.m_tid
         << ", \"type\": " << a_event.m_type
         << ", \"execution number\": " << a_event.m_execNum
         << std::setprecision(std::numeric_limits< long double>::digits10 + 1)
         << ", \"start\": " << a_event.m_startTime
         << ", \"end\": " << a_event.m_endTime << "}";
      return os;
    }

private:
  int         m_type;
  int         m_id;
  int         m_pid;
  int         m_tid;
  int         m_execNum;
  double      m_startTime;
  double      m_endTime;
  char        m_name[40];
};

class DependencyRecord
{
    friend std::ostream& operator<<(std::ostream& os, 
                                    const DependencyRecord a_dependency)
    {
      os << "{\"type\": \"" << a_dependency.m_type
         << "\", \"tailTask\": " << a_dependency.m_tailEventID 
         << ", \"headTask\":" << a_dependency.m_headEventID << "}";
      return os;
    }

private:
  int         m_type;
  int         m_tailEventID;
  int         m_headEventID;
};

}

#include "NamespaceFooter.H"

#endif  /* ! defined _EVENT_H_ */
