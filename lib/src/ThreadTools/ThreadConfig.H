#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _THREADCONFIG_H_
#define _THREADCONFIG_H_

#include "NamespaceHeader.H"

namespace ThreadTools
{
#ifdef USE_PROFILER
  constexpr bool c_useProfiler = true;
#else
  constexpr bool c_useProfiler = false;
#endif
}

#include "NamespaceFooter.H"

#endif
