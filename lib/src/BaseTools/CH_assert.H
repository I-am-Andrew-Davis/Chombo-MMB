#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CH_ASSERT_H_
#define _CH_ASSERT_H_

#ifdef __CUDACC__
#include <cassert>
#endif

#ifdef CH_LANG_CC
#include "MayDay.H"
#endif

#include "BaseNamespaceHeader.H"

#ifdef CH_LANG_CC
/**
@ingroup macros
compile-time cpp macro that replaces 'assert()' and aborts
correctly when running in parallel. */
# ifndef NDEBUG
#  define CH_assert_str(s) # s
#  define CH_assert_xstr(s) CH_assert_str(s)
#  ifdef __CUDACC__
#   define CH_assert(cond) assert(cond)
#  else
#   define CH_assert(cond) if (!(cond))                                                                                      \
                           {                                                                                                 \
                             CH_XD::MayDay::Abort( __FILE__ ":" CH_assert_xstr(__LINE__) ": Assertion `" #cond "' failed."); \
                           }
#  endif  /* !CUDA */
# else
#  define CH_assert(cond) (void)0
#endif

/*  removing this form of static_assert since it does not follow the
std::static_assert specification and users can now use real static_assert
*/
/** Compile-time assertion */
/** (renamed this temp variable as static_assert will become a keyword -- and likely have same functionality) */
/*
template<bool> struct diy_static_assert;
template<> struct diy_static_assert<true>
{
};
*/
//#define STATIC_ASSERT( a ) diy_static_assert<a>()

#endif /* CH_LANG_CC */

#ifdef  CH_LANG_FORT
#ifndef NDEBUG
#  define CH_assert(cond) if (.NOT. cond) call MAYDAYERROR()
#else
#  define CH_assert(cond)
#endif
#endif /* CH_LANG_FORT */

#include "BaseNamespaceFooter.H"

#endif  /* ! defined _CH_ASSERT_H_ */
