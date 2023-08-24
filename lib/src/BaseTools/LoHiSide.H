#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LOHISIDE_H_
#define _LOHISIDE_H_

#include <cassert>
#include <type_traits>

#include "CH_assert.H"

/*
  Notes:
    - tag __CUDACC__ means the file is processed with nvcc.  The class is
      modified whether on the cpu or gpu.
    - tag __CUDA_ARCH__ is defined for device code.  This only works in
      functions with the __device__ qualifier.  It cannot be used in the
      class declaration.
*/

#undef DEVICE
#ifdef __CUDACC__
#define DEVICE __device__
#else
#define DEVICE
#endif

#include "BaseNamespaceHeader.H"

//--Forward declarations

namespace CH_Hash
{
  template <typename T>
  struct isHashable;
}

/// Encapsulation of high or low side.
/**
   Encapsulation of high or low side.
*/
namespace Side
{
  ///
  /**
     Lo is the low side.<br>
     Hi is the high side.
  */
  enum LoHiSide : int
  {
    Invalid=-1,
    Lo=0,
    Hi,
    NUMSIDES
  };

  /// Returns the other side.  That is, swap Hi for Lo and vice versa.
  constexpr DEVICE LoHiSide
  flip(const LoHiSide a_side) noexcept
  {
    assert((a_side == Lo) || (a_side == Hi));
    return static_cast<LoHiSide>(!static_cast<int>(a_side));
  }

  /// Returns the sign of the side: -1 for Lo, 1 for Hi
  constexpr DEVICE int
  sign(const Side::LoHiSide& a_side) noexcept
  {
    assert((a_side == Side::Lo) || (a_side == Side::Hi));
    return 2*static_cast<int>(a_side) - 1;
  }
}

namespace CH_Hash
{

/// Verify that LoHiSide can be byte-hashed and give bytes to hash
template <>
struct isHashable<Side::LoHiSide> : public std::true_type
{
  // Assert LoHiSide is POD
  static_assert(sizeof(Side::LoHiSide) == sizeof(int),
                "Bad assumption for sizeof(LoHiSide)");
  static constexpr int c_hashSize = sizeof(Side::LoHiSide);
};

}

// ///
// /**
//    Returns the other side.  That is, swap Hi for Lo and vice versa.
// */
// inline DEVICE Side::LoHiSide
// flip(const Side::LoHiSide a_side)
// {
//   return Side::flip(a_side);
// }

// ///
// /**
//    Returns -1 for Lo Side, +1 for Hi Side.
// */
// constexpr DEVICE int
// sign(const Side::LoHiSide& a_side) noexcept
// {
//   return Side::sign(a_side);
// }

/// This class allows one to do 'for (auto side : EachSide)'
class _SideRange
{
public:
  constexpr _SideRange() noexcept
    :
    m_side(static_cast<int>(Side::LoHiSide::NUMSIDES))
    { }
  constexpr _SideRange(const Side::LoHiSide& a_side) noexcept
    :
    m_side(static_cast<int>(a_side))
    { }
  constexpr _SideRange begin() const noexcept
    {
      return _SideRange(Side::LoHiSide::Lo);
    }
  constexpr _SideRange end() const noexcept
    {
      return _SideRange(Side::LoHiSide::NUMSIDES);
    }
  constexpr Side::LoHiSide operator*() const noexcept
    { return static_cast<Side::LoHiSide>(m_side); }
  constexpr void operator++() noexcept
    { ++m_side; }
  constexpr bool operator!=(const _SideRange& a_right) const noexcept
    { return m_side != a_right.m_side; }
private:
  int m_side;
};
constexpr _SideRange EachSide = _SideRange{};

// NOTE: Use 'for (auto side : EachSide)' instead of the SideIterator
//       This is only for legacy code
/// Iterator for low and high side.
/**
   SideIterator iterates over the values Side::Lo and Side::Hi.
   There are no adjustable parameters whatsoever.
*/
class SideIterator
{
public:
  ///
  /**
     Default, and only, constructor.
  */
  DEVICE SideIterator()
    :
    m_current(0)
    { }

  /**
     \name Modification functions
  */
  /*@{*/

  ///
  /**
     Set the iterator to the first side.
  */
  DEVICE void begin()
    {
      m_current = 0;
    }

  ///
  /**
     Set the iterator to the first side.
  */
  DEVICE void reset()
  {
    begin();
  }

  ///
  /**
     Advance the iterator to the next face.
  */
  DEVICE void next()
    {
      ++m_current;
    }

  ///
  /**
     Advance the iterator to the next face.
  */
  DEVICE void operator++ ()
    {
      ++m_current;
    }

  /*@}*/

  /**
     \name Access functions
  */
  /*@{*/

  ///
  /**
     Returns the value of the side for the current iteration.
  */
  DEVICE Side::LoHiSide operator () () const
    {
      switch (m_current)
        {
        case 0:
          return Side::Lo;
        case 1:
          return Side::Hi;
        default:
          return Side::Invalid;
        }
    }  

  ///
  /**
     Returns <tt>true</tt> if there is a next side, <tt>false</tt> if not.
  */
  DEVICE bool ok() const
    {
      return ((m_current > -1) && (m_current < Side::NUMSIDES));
    }

  /*@}*/

private:
  int m_current;
};

#include "BaseNamespaceFooter.H"

#endif
