#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MEMTRACK_H_
#define _MEMTRACK_H_

#ifdef CH_USE_MEMORY_TRACKING

#include <iostream>
using std::ostream;

#include "REAL.H"
#include "CH_assert.H"

// NOTE:  These were added to avoid macro substitutions of "calloc(a,b)",
// "realloc(a,b)", and "free(a)" for member function(s) of the "Arena" class.
//
// These member function(s) are declared in "Arena.H" and used
// "BaseFabImplem.H" (which is included via "BaseFab.H").
#include "Arena.H"
//#include "BaseFab.H"

// Use these macros to include in memory tracking system
#define callocMT(a_a,a_b)   callocp (a_a, a_b, __FILE__, __LINE__)
#define reallocMT(a_a, a_b) reallocp(a_a, a_b, __FILE__, __LINE__)
#define mallocMT(a_a)       mallocp (a_a,      __FILE__, __LINE__)
#define freeMT(a_a)         freep   (a_a)

#include "BaseNamespaceHeader.H"

extern unsigned long long int ch_memcount;
///
void dumpmemorymap(int a_sig);

///
void dumpmemoryatexit();

///
void dumpmemoryabort();

///
int registerMemorySignals();

///
void ReportUnfreedMemory(ostream& a_os);

///
void ReportAllocatedMemory(ostream& a_os);

///
void ReportMemoryUsage(ostream& a_os);

/// calls ReportUnfreedMemory with pout()
void UnfreedMemory();

/// calls  ReportMemoryUsage with pout()
void MemoryUsage();

///
void memTrackingOn();

///
void memtrackingOff();

void overallMemoryUsage(long long& a_currentTotal,
                        long long& a_peak);

void overallMemoryUsage();

void memtrackStamp(Real& a_current,
                   Real& a_peak);

/// Memory tracking functions
/**
   Memory tracking functions
*/

class Memtrack
{
public:
  /// Produce a formated output onto os of memory usage.
  /**
     Memory is reported one a class-by-class basis.  When
     CH_USE_MEMORY_TRACKING is turned on you also get a report from
     the atexit() function.  This can be handy for spotting memory
     leaks.  The memory tracking functionality consumes very little
     runtime and does not impede optimized performance.
  */
  static void ReportUnfreedMemory(ostream& a_os);

  /// calls ReportUnfreedMemory with pout()
  static void UnfreedMemory();

  static void memTrackingOn();

  static void memtrackingOff();

  static void overallMemoryUsage(long long& a_currentTotal,
                                 long long& a_peak);
  static std::vector<void (*)()> s_memtrackAtExit;	
};

/// Memory tracking registry
/**
    A place to tally the bytes used.  The registry automatically records
    entries so that the results can be recorded at exit.  Allocators for
    array-type containers, such as BaseFab, use this to record the memory that
    has been allocated.  These concepts were originally implemented in the
    Arena class.
*/

class MemtrackRegistry;
typedef std::list<MemtrackRegistry*> MemtrackRegList;

const int g_MemtrackRegNameSize = 128;
                                      ///< Size of "name" character array

class MemtrackRegistry
{
public:

  enum class Device
  {
    CPU,
    GPU
  };
    
  MemtrackRegistry();
  MemtrackRegistry(const char*  a_name,
                   const Device a_device = Device::CPU);
  MemtrackRegistry(const std::string& a_name,
                   const Device       a_device = Device::CPU);
  ~MemtrackRegistry();

  /// Add bytes to the amount tracked (same to used and buffered)
  void addBytes(const size_t a_numBytes);

  /// Add bytes to the amount tracked (separate used and buffered)
  void addBytes(const size_t a_numBytes, const size_t a_numBytesBuffered);

  /// Subtract bytes to the amount tracked (same to used and buffered)
  void subtractBytes(const size_t a_numBytes);

  /// Subtract bytes to the amount tracked (separate used and buffered)
  void subtractBytes(const size_t a_numBytes, const size_t a_numBytesBuffered);

  /// Get total bytes (used and buffered)
  size_t bytes() const { return m_bytesBuffered; }

  /// Get  total MiB (1024^2 bytes) (used and buffered)
  Real mebibytes() const;

  /// Get peak bytes (used and buffered)
  size_t peakBytes() const { return m_peakBuffered; }

  /// Get peak MiB (1024^2 bytes) (used and buffered)
  Real peakMebibytes() const;

  /// Get name
  const char* name() const { return m_name; }

  /// Get device
  Device device() const { return m_device; }

  /// On CPU
  bool onCPU() const { return m_device == Device::CPU; }

  /// On GPU
  bool onGPU() const { return m_device == Device::GPU; }
  
private:

  /// Get a list of registry entries
  static MemtrackRegList& memtrackRegList();

  char m_name[g_MemtrackRegNameSize]; ///< Name associated with this registry
                                      ///< entry
  size_t m_bytes;                     ///< Number of bytes currently used
  size_t m_peak;                      ///< Peak number of bytes currently used
  size_t m_bytesBuffered;             ///< Number of bytes currently allocated
                                      ///< including buffers (e.g. for
                                      ///< alignment)
  size_t m_peakBuffered;              ///< Peak number of bytes currently 
                                      ///< allocated including buffers (e.g. for
                                      ///< alignment)
  Device m_device;                    ///< Processor device

  // Deletion of the list occurs herein
  friend void dumpmemoryatexit();
  // Access to the list
  friend void ReportUnfreedMemory(ostream& a_os);
  friend void ReportAllocatedMemory(ostream& a_os);
};

//  Add bytes to the amount tracked (same to used and buffered)
inline void
MemtrackRegistry::addBytes(const size_t a_numBytes)
{
  m_bytes += a_numBytes;
  m_bytesBuffered += a_numBytes;
  m_peak = std::max(m_peak, m_bytes);
  m_peakBuffered = std::max(m_peakBuffered, m_bytesBuffered);
}

//  Add bytes to the amount tracked (separate used and buffered)
/** \param[in]  a_numBytes
 *                      Number of bytes used
 *  \param[in]  a_numBytesBuffered
 *                      Total number of bytes used including used and
 *                      buffered
 */
inline void
MemtrackRegistry::addBytes(const size_t a_numBytes,
                           const size_t a_numBytesBuffered)
{
  CH_assert(a_numBytesBuffered >= a_numBytes);
  m_bytes += a_numBytes;
  m_peak = std::max(m_peak, m_bytes);
  m_bytesBuffered += a_numBytesBuffered;
  m_peakBuffered = std::max(m_peakBuffered, m_bytesBuffered);
}

//  Subtract bytes to the amount tracked (same to used and buffered)
inline void
MemtrackRegistry::subtractBytes(const size_t a_numBytes)
{
  m_bytes -= a_numBytes;
  m_bytesBuffered -= a_numBytes;
}

// Subtract bytes to the amount tracked (separate used and buffered)
/** \param[in]  a_numBytes
 *                      Number of bytes unused
 *  \param[in]  a_numBytesBuffered
 *                      Total number of bytes unused including used
 *                      and buffered
 */
inline void
MemtrackRegistry::subtractBytes(const size_t a_numBytes,
                                const size_t a_numBytesBuffered)
{
  CH_assert(a_numBytesBuffered >= a_numBytes);
  m_bytes -= a_numBytes;
  m_bytesBuffered -= a_numBytesBuffered;
}

//  Get total MiB (1024^2 bytes) (used and buffered)
inline Real
MemtrackRegistry::mebibytes() const
{
  return (Real)m_bytesBuffered/(Real)1048576;
}

//  Get peak MiB (1024^2 bytes) (used and buffered)
inline Real
MemtrackRegistry::peakMebibytes() const
{
  return (Real)m_peakBuffered/(Real)1048576;
}

// =========== end user functions.===================================

void AddTrack(void*       a_addr,
              size_t      a_asize,
              const char* a_fname,
              int         a_lnum,
              bool        a_malloc);

void RemoveTrack(void* a_addr,
                 bool  a_malloc);

// void* operator new (size_t      a_size,
//                     char const* a_file,
//                     int         a_line);

// void operator delete (void *a_p) throw();

// void* operator new[] (size_t      a_size,
//                       char const* a_file,
//                       int         a_line);

// void operator delete[] (void *a_p) throw();

// Don't call these directly -- use the callocMT/mallocMT/freeMT macros.
// Use these functions to track the source file and line number of particular
// memory alloc/frees.
void* mallocp(size_t      a_size,
              const char* a_file,
              int         a_line);

void* reallocp(void*       a_p,
               size_t      a_size,
               const char* a_file,
               int         a_line);

void* callocp(size_t      a_nelem,
              size_t      a_elsize,
              const char* a_file,
              int         a_line);

void freep(void* a_p);

#include "BaseNamespaceFooter.H"

#else // on CH_USE_MEMORY_TRACKING

// If not compiling with memory tracking, then just call the LIBC versions
#include <cstdlib>
#define callocMT(a_a,a_b)   calloc (a_a, a_b)
#define reallocMT(a_a, a_b) realloc(a_a, a_b)
#define mallocMT(a_a)       malloc (a_a)
#define freeMT(a_a)         free   (a_a)

inline void AddTrack(void*       a_addr,
              size_t      a_asize,
              const char* a_fname,
              int         a_lnum,
              bool        a_malloc)
{
  // null op
}

inline void RemoveTrack(void* a_addr,
                 bool  a_malloc)
{
  // null op
}

#endif // ifdef CH_USE_MEMORY_TRACKING

#endif // include guard
