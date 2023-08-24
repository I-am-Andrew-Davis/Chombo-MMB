#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CH_Thread_H_
#define _CH_Thread_H_

#if CH_OPENMP==1
#include <omp.h>
#endif

#include "BaseNamespaceHeader.H"

//returns true if we are on thread 0 (or if not threaded)
extern bool onThread0();

//returns the value of OMP_NUM_THREADS (or 1 if not threaded)
extern int  getMaxThreads();

#include "BaseNamespaceFooter.H"

#endif  /* ! defined _CH_Thread_H_ */
