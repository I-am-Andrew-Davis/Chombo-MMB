#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//----- Standard Library -----//

#include <iostream>

//----- Internal -----//

#include "Profiler.H"

#include "NamespaceHeader.H"

using namespace Profiling;

/*******************************************************************************
 *
 * Class Profiler: member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Default constructor
/**
 *//*-----------------------------------------------------------------*/

Profiler::Profiler()
  :
  m_events(nullptr),
  m_processID(-1),
  m_threadID(-1),
  m_numEvents(0),
  m_isOpenEvent(false)
{
}

/*--------------------------------------------------------------------*/
//  Constructor
/** \param[in]  a_filename
 *                        The name of the file to create for the
 *                        profiler
 *  \param[in]  a_processID
 *                        The id of the process. Important for MPI
 *  \param[in]  a_threadID
 *                        The id of the thread that is being
 *                        profiled
 *  \param[in]  a_startTime
 *                        The time at the start of the profile
 *  \param[in]  a_maxNumEvents
 *                        The maximum number of events to hold in
 *                        the buffer
 *//*-----------------------------------------------------------------*/

Profiler::Profiler(const std::string a_filename,
                   const int         a_processID,
                   const int         a_threadID,
                   const time_point  a_startTime,
                   const int         a_maxNumEvents)
  :
  m_events(nullptr),
  m_processID(a_processID),
  m_threadID(a_threadID),
  m_startTime(a_startTime),
  m_maxNumEvents(a_maxNumEvents),
  m_numEvents(0),
  m_isOpenEvent(false)
{
  m_ofs.open(a_filename.c_str(), 
              std::ios::binary | std::ios::out | std::ios::trunc);
  writeHeader();
  allocate();
}

/*--------------------------------------------------------------------*/
//  Destructor
/**
 *//*-----------------------------------------------------------------*/

Profiler::~Profiler()
{
  end();
  deallocate();
}

/*--------------------------------------------------------------------*/
//  Weak constructor
/** \param[in]  a_filename
 *                        The name of the file to create for the
 *                        profiler
 *  \param[in]  a_processID
 *                        The id of the process. Important for MPI
 *  \param[in]  a_threadID
 *                        The id of the thread that is being
 *                        profiled
 *  \param[in]  a_startTime
 *                        The time at the start of the profile
 *  \param[in]  a_maxNumEvents
 *                        The maximum number of events to hold in
 *                        the buffer
 *//*-----------------------------------------------------------------*/

void Profiler::define(const std::string a_filename,
                      const int         a_processID,
                      const int         a_threadID,
                      const time_point  a_startTime,
                      const int         a_maxNumEvents)
{
  m_events = nullptr;
  m_processID = a_processID;
  m_threadID = a_threadID;
  m_startTime = a_startTime;
  m_maxNumEvents = a_maxNumEvents;
  m_numEvents = 0;
  m_ofs.open(a_filename.c_str(), 
              std::ios::binary | std::ios::out | std::ios::trunc);
  writeHeader();
  allocate();
}

/*--------------------------------------------------------------------*/
//  Allocate memory for the events
/**
 *//*-----------------------------------------------------------------*/

void Profiler::allocate()
{
  if(m_events == nullptr)
    {
      m_events = new Event[m_maxNumEvents];
    }
}

/*--------------------------------------------------------------------*/
//  Deallocate memory for the events
/**
 *//*-----------------------------------------------------------------*/

void Profiler::deallocate()
{
  if(m_events != nullptr)
    {
      delete[] m_events;
    }
}

/*--------------------------------------------------------------------*/
//  Write header to the ofstream
/**
 *//*-----------------------------------------------------------------*/

void Profiler::writeHeader()
{
  // need time_t to print time with std::ctime(time_t)
  // need a system_clock to convert to time_t
  std::time_t l_start_time = 
    system_clock::to_time_t(system_clock::now()
                            +(m_startTime-steady_clock::now()));
  
  ProfileHeader header;
  header.processID = m_processID;
  header.threadID = m_threadID;
  header.bufferSize = m_maxNumEvents;
  strcpy(header.startTime,std::ctime(&l_start_time));
  m_ofs.write((char*)&(header), sizeof(ProfileHeader));
  m_ofs.flush();
}

/*--------------------------------------------------------------------*/
//  Display the buffer for debugging
/**
 *//*-----------------------------------------------------------------*/

void Profiler::displayBuffer() const
{
  for (int i=0; i < m_numEvents; i++)
    {
      std::cout << m_events[i] << std::endl;
    }
}

/*--------------------------------------------------------------------*/
//  Write the buffer to the ofstream
/**
 *//*-----------------------------------------------------------------*/

inline
void Profiler::flush()
{
  if(m_isOpenEvent)
    {
      endEvent();
      // might want to print a warning
    }
  Event writeEvent(0,
                   0,
                   m_processID,
                   m_threadID,
                   "writeEvents",
                   0,
                   duration(
                     steady_clock::now()-m_startTime).count());
  m_ofs.write((char*)(m_events), m_numEvents*sizeof(Event));
  writeEvent.end(duration(steady_clock::now()-m_startTime).count());
  m_ofs.write((char*)(&writeEvent), sizeof(Event));
  m_numEvents = 0;
}

#include "NamespaceFooter.H"
