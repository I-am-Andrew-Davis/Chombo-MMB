#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MISC_H_
#define _MISC_H_

#undef HOSTDEVICE
#undef DEVICE
#ifdef __CUDACC__
#define HOSTDEVICE __host__ __device__
#define DEVICE __device__
#else
#define HOSTDEVICE
#define DEVICE
#endif  /* !CUDA */

#include "CH_assert.H"
#include "REAL.H"

#ifndef __CUDACC__
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <algorithm>
#include <bitset>

#include "BaseNamespaceHeader.H"
#endif  /* !CUDA */

/**
   Returns a copy of the minimum of the two values.  Will work with
   any type T that supplies a T::operator<(const T, const T), or an
   equivalent construct, as well as a copy constructor.   In
   particular, it will work with any integral or floating-point
   intrinsic type.
*/
template <class T>
HOSTDEVICE inline T
Min(const T& a_a, const T& a_b)
{
  return (a_a < a_b) ? a_a : a_b;
}

/**
   Returns a copy of the maximum of the two values.  Will work with
   any type T that supplies a T::operator>(const T, const T), or an
   equivalent construct, as well as a copy constructor.   In
   particular, it will work with any integral or floating-point
   intrinsic type.
*/
template <class T>
HOSTDEVICE inline T
Max (const T& a_a, const T& a_b)
{
  return (a_a > a_b) ? a_a : a_b;
}

/**
   Returns a copy of the absolute value of the value.  Will work with
   any type T that can be compared against zero with
   T::operator>(const T, const T), or an equivalent construct, as
   well as a copy constructor, and an T::operator-() or equivalent.
   In particular, it will work with any integral or floating-point
   intrinsic type.
*/
template <class T>
HOSTDEVICE inline T
Abs(const T& a_a)
{
  return (a_a > 0) ? a_a : -a_a;
}

/**
   Swaps the two values.  Type T must have a copy constructor and
   an assignment operator.
*/
template <class T>
HOSTDEVICE inline void
Swap(T& a_a, T& a_b)
{
  T tmp = a_a;
  a_a   = a_b;
  a_b   = tmp;
}

//inline int pow(int a, int b)
//{
//  CH_assert(a>=0);
//  CH_assert(b>=0);
//  int rtn = 1;
//  for (; b>0; b--)
//    rtn*=a;
//  return rtn;
//}
template <unsigned int P>
HOSTDEVICE constexpr inline int ipow(int M) { return M*ipow<P-1>(M);}
template<>
HOSTDEVICE constexpr inline int ipow<0>(int M) {return 1;}

HOSTDEVICE inline int
ipow(int a, int b)
{
  CH_assert(b>=0);
  int rtn = 1;
  for (; b>0; b--)
  {
     rtn *= a;
  }
  return rtn;
 // return std::pow(a,b);
}

HOSTDEVICE inline Real
ipow(const Real& a, const int& b)
{
  return std::pow(a, b);

//   double rtn = 1.0;
//   int r=Abs(b);
//   for (;r>0; r--)
//     rtn*=a;
//   return (b>0) ? rtn : 1.0/rtn;
}

/*------------------------------------------------------------------------------
 * The following rarely used functions are in namespace Misc to avoid conficts
 * and help identify their origin
 *----------------------------------------------------------------------------*/

namespace Misc
{

/// constexpr popcnt
/** Note: suitable for compile-time evaluation only.  See popcnt for the
 *  run-time version.  Will not work on a negative integer.
 */
template <typename T>
inline constexpr T
cpopcnt(const T a_i)
{
  static_assert(std::is_integral<T>::value, "popcnt on non-integer");
  return (a_i == (T)0) ? (T)0 : (T)(a_i & (T)1) + cpopcnt(a_i >> 1);
}

/*
  Portable popcnt (if it does not use popcnt instruction when available, get a
  decent compiler).  A test file is embedded here.  E.g., compile with
  g++ -march=native -O3 -S filename.cpp and look for popcntq in the assembly
*/
/*----------------------------------------------------------------------------*/
/*
#include <cstdio>
#include <bitset>
#include <limits>

#define HOSTDEVICE
*/
// <= 32 bit popcnt
template <typename T>
HOSTDEVICE inline auto static
popcntImpl(const T a_i, std::true_type)
{
  static_assert(sizeof(T) <= sizeof(unsigned),
                "Insufficient size of integer for type sent to 32-bit popcnt");
#ifdef __CUDA_ARCH__
  return popc(static_cast<unsigned>(a_i));
#else
  return std::bitset<std::numeric_limits<unsigned>::digits>(a_i).count();
#endif
}
// > 32 bit popcnt
template <typename T>
HOSTDEVICE inline auto static
popcntImpl(const T a_i, std::false_type)
{
  static_assert(sizeof(T) <= sizeof(unsigned long long),
    "Insufficient size of integer for type sent to popcnt");
#ifdef __CUDA_ARCH__
  return popcll(static_cast<unsigned long long>(a_i));
#else
  return std::bitset<std::numeric_limits<unsigned long long>::digits>(a_i)
    .count();
#endif
}
/// Popcnt
/** Dispatch based on size
 */
template <typename T>
HOSTDEVICE inline auto
popcnt(const T a_i)
{
  static_assert(std::is_integral<T>::value, "popcnt on non-integer");
  return popcntImpl(a_i, std::integral_constant<bool, (sizeof(T) <= 4)>{});
}
/*
int main()
{
  volatile unsigned u = 15;
  int i = popcnt(u);
  printf("%d\n", i);
}
*/
/*----------------------------------------------------------------------------*/

// Check for power of 2
template <typename T>
HOSTDEVICE inline constexpr bool
isPower2(const T a_i)
{
  static_assert(std::is_integral<T>::value, "isPower2 on non-integer");
  static_assert(sizeof(T) <= sizeof(unsigned long long),
    "Insufficient size of integer for type sent to isPower2");
  const unsigned long long u = a_i;
  return (a_i > 0) && !(u & (u-1));
}

/// Constexpr log2 (expects a_i is power of 2, otherwise this is the floor)
template <typename T>
HOSTDEVICE inline constexpr T
clog2(const T a_i)
{
  static_assert(std::is_integral<T>::value, "clog2 on non-integer");
  return (a_i == (T)1) ? (T)0 : (T)1 + clog2(a_i >> 1);
}

/// Return greater or equal size that is aligned
//  Consider using C++11 alignof and/or std::aligned_storage instead
HOSTDEVICE inline std::uintptr_t
align(std::uintptr_t a_s, const std::uintptr_t a_alignment)
{
  CH_assert(a_alignment > 0);
  return a_alignment*((a_s + a_alignment - 1)/a_alignment);
}

/// Check if an address is aligned (alignment must be a power of 2)
HOSTDEVICE inline bool
isAligned(const void* a_p, const std::uintptr_t a_alignment)
{
  CH_assert(isPower2(a_alignment));
  return (!(reinterpret_cast<std::uintptr_t>(a_p) & (a_alignment - 1)));
}

/// Check if an address has sufficient alignment for a type
template <typename T>
HOSTDEVICE inline bool
isAligned(const void* a_p) noexcept
{
  return isAligned(a_p, alignof(T));
}

// NOTE - As of C++11, use std::is_class<T>::value instead
/// Test for class type
/**
   Use the SFINAE principle as described in 15.2.2 "Determining
   Class Types" in Vandevoorde and Josuttis "C++ Templates" book to
   see if T is class type.  E.g.,
   TypeTr<int>::IsClass == false
   TypeTr<Box>::IsClass == true
*/
#ifndef __CUDACC__
template <typename T>
class TypeTr
{
private:
  typedef char One;
  typedef struct
  {
    char a[2];
  } Two;
  template <typename C> static One test(int C::*);
  template <typename C> static Two test(...);
public:
  enum
  {
    IsClass = sizeof(TypeTr<T>::template test<T>(0)) == 1
  };
};
// Use with above to test if false.  E.g.,
//   typedef typename TestFalse<TypeTr<T>::IsClass>::Check FalseValue;
template <bool B>
struct TestFalse;
template <>
struct TestFalse<false>
{
  typedef int Check;
};
#endif

/// Comparison with limit tol^2 as x and y -> 0.
/** \param[in]  x       First operand
 *  \param[in]  y       Second operand
 *  \param[in]  prec    Number of base 10 significant digits to compare
 *  \return             T - Not equal
 */
template <typename T>
HOSTDEVICE inline bool
compare(const T &x, const T &y, int prec)
{
  const T tol = std::pow(10., -std::abs(prec));
#ifdef __CUDACC__
  return std::fabs(x - y) >
    (::min(std::fabs(x), std::fabs(y)) + tol)*tol;
#else
  return std::fabs(x - y) >
    (std::min(std::fabs(x), std::fabs(y)) + tol)*tol;
#endif
}

/// Comparison with limit tol as x and y -> 0.
/** This is more useful if you expect x or y to be garbage as they approach 0
 *  \param[in]  x       First operand
 *  \param[in]  y       Second operand
 *  \param[in]  prec    Number of base 10 significant digits to compare
 *  \return             T - Not equal
 */
template <typename T>
HOSTDEVICE inline bool
compare1(const T &x, const T &y, int prec)
{
  const T tol = std::pow(10., -std::abs(prec));
#ifdef __CUDACC__
  return std::fabs(x - y) >
    ::min(std::fabs(x), std::fabs(y))*tol + tol;
#else
  return std::fabs(x - y) >
    std::min(std::fabs(x), std::fabs(y))*tol + tol;
#endif
}

}  // namespace Misc

#ifndef __CUDACC__
#include "BaseNamespaceFooter.H"
#endif  /* ! CUDA */

#endif  /* ! defined _MISC_H_ */
