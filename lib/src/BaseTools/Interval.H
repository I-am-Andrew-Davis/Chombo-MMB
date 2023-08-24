#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _INTERVAL_H_
#define _INTERVAL_H_
#include "BaseNamespaceHeader.H"

//
// Little struct for passing component ranges in code.
//

/// Structure for passing component ranges in code.
/**
   A class to encapsulate component ranges
*/
struct Interval
{
  constexpr Interval() noexcept
    :
    m_begin(0),
    m_end(-1)
  {}

  ///
  /**
      initialize with first and last component numbers
  */
  constexpr Interval(int a_firstComp, int a_lastComp) noexcept
    :
    m_begin(a_firstComp),
    m_end(a_lastComp)
  {}

  ///
  /**
   * Copy from another Interval
   */
  Interval& operator=(const Interval& p)
  {
    m_begin = p.m_begin;
    m_end = p.m_end;
    return *this;
  }

  void define(int a_firstComp,
              int a_lastComp)
  {
    m_begin = a_firstComp;
    m_end   = a_lastComp;
  }

#ifndef WRAPPER
  ///
  /**
     return first component number
  */
  constexpr int begin() const noexcept
    { return m_begin; }

  ///
  /**
     return last component number
  */
  constexpr int end() const noexcept
    { return m_end; }
#else
  int begin() const noexcept;
  int end() const noexcept;
#endif

  ///
  /**
     return last - first + 1
   */
  constexpr int size() const noexcept
  {
    return m_end-m_begin+1;
  }

  /// check if interval contains element
  constexpr bool contains(int a_val) const noexcept
  {
    return a_val >= m_begin && a_val <= m_end;
  }

  /// check if interval contains another interval
  constexpr bool contains(const Interval& a_intv) const noexcept
  {
    return a_intv.m_begin >= m_begin && a_intv.m_end <= m_end;
  }
  
  /// test for equality
  constexpr bool operator==(const Interval& p) const noexcept
  {
    return ( (m_begin == p.m_begin) && (m_end == p.m_end));
  }

  friend std::ostream& operator<< (std::ostream& os, const Interval& dit);

private:
  int m_begin, m_end;
};

#include "BaseNamespaceFooter.H"
#endif
