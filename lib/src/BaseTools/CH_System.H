#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CH_SYSTEM_H_
#define _CH_SYSTEM_H_

/******************************************************************************/
/**
 * \file
 *
 * \brief System dependent functions
 *
 *//*+*************************************************************************/

#include <iostream>
#include <cstring>
#include <type_traits>

#include "CH_config.H"
#include "Misc.H"
#include "BaseNamespaceHeader.H"

namespace CH_System
{
  /// Check if a file exists
  int fileExists(const char *const a_filename);

  /// Allocate aligned memory
  int memalign(void **a_memptr, size_t a_alignment, size_t a_size);

  /// Get the path and name of the currently running executable
  int getProcessPath(char *const a_procPath, const int a_len);

/*--------------------------------------------------------------------*/
/// Coarsen by a power of 2
/**
 *  \param[in]  a_i     Integer
 *  \param[in]  a_nref  Refinement ratio (must be power of 2)
 *  Cannot be constexpr because of popcnt
 *//*-----------------------------------------------------------------*/

  template <typename T1, typename T2>
  inline T1 coarsen(const T1& a_i, const T2& a_nref)
  {
    static_assert(std::is_integral<T1>::value &&
                  std::is_integral<T2>::value,
                  "Coarsen can only be applied to integral vector types");
#if CHDEF_SYSTEM_CPPARS == 1
    CH_assert(Misc::isPower2(a_nref));
    return a_i >> Misc::popcnt(a_nref-1);
#else
    const int i_nref = (int)a_nref;
    return (a_i<0) ? (a_i+1)/i_nref - 1 : a_i/i_nref;
#endif
  }

/*--------------------------------------------------------------------*/
/// Coarsen by a power of 2 given by 2^a_x
/**
 *  \param[in]  a_i     Integer
 *  \param[in]  a_x     Refinement ratio given by 2^a_x (a_x must be
 *                      >= 0)
 *//*-----------------------------------------------------------------*/

  template <typename T1, typename T2>
  inline constexpr T1 coarsen_log2(const T1& a_i, const T2& a_x)
  {
    static_assert(std::is_integral<T1>::value &&
                  std::is_integral<T2>::value,
                  "Coarsen can only be applied to integral vector types");
#if CHDEF_SYSTEM_CPPARS == 1
    return a_i >> a_x;
#else
    return (a_i<0) ? (a_i+1)/((int)1 << a_x) - 1 : a_i>>a_x;
#endif
  }
}  // namespace CH_System

// BOOST_PP_STRINGIZE(expr) - Return expr as a string literal.
#define BOOST_PP_STRINGIZE_(expr) #expr
#define BOOST_PP_STRINGIZE(expr) BOOST_PP_STRINGIZE_(expr)

// BOOST_PRAGMA(args) - Emits a pragma.
#define BOOST_PRAGMA(args) _Pragma(BOOST_PP_STRINGIZE(args))
#if   defined(__INTEL_COMPILER)
    #define BOOST_BUILTIN_ASSUME(expr)                                         \
        __assume(expr)                                                         \
        /**/
    #define BOOST_BUILTIN_ASSUME_ALIGNED(ptr, align)                           \
        __assume_aligned(ptr, align)                                           \
        /**/
#elif defined(__clang__)
    #define BOOST_BUILTIN_ASSUME(expr)                                         \
        __builtin_assume(expr)                                                 \
        /**/
    #define BOOST_BUILTIN_ASSUME_ALIGNED(ptr, align)                           \
        {                                                                      \
            ptr = reinterpret_cast<decltype(ptr)>(                             \
                __builtin_assume_aligned(ptr, align)                           \
            );                                                                 \
        }                                                                      \
        /**/
#else // GCC
    #define BOOST_BUILTIN_ASSUME(expr)                                         \
        do { if (!(expr)) __builtin_unreachable(); } while (0)                 \
        /**/
    #define BOOST_BUILTIN_ASSUME_ALIGNED(ptr, align)                           \
        {                                                                      \
            ptr = reinterpret_cast<decltype(ptr)>(                             \
                __builtin_assume_aligned(ptr, align)                           \
            );                                                                 \
        }                                                                      \
        /**/
#endif
#if   defined(__INTEL_COMPILER)
    #define BOOST_DEMAND_VECTORIZATION                                         \
        BOOST_PRAGMA(simd)                                                     \
        /**/
#elif defined(__clang__)
    #define BOOST_DEMAND_VECTORIZATION                                         \
        BOOST_PRAGMA(clang loop vectorize(enable) interleave(enable))          \
        /**/
#else
    #define BOOST_DEMAND_VECTORIZATION                                         \
        BOOST_PRAGMA(GCC ivdep)                                                \
        /**/
#endif

#include "BaseNamespaceFooter.H"
#endif
