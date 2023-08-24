#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "CH_Thread.H"
#include "NamespaceHeader.H"

#ifdef CH_CTHR
namespace ThreadTools
{
extern thread_local int tls_threadId;
}
#endif

bool onThread0()
{
  bool retval = true;
#if CH_OPENMP==1
  int thread_num = omp_get_thread_num();
  retval = (thread_num== 0);
#elif defined CH_CTHR
  int thread_num = ThreadTools::tls_threadId;
  retval = (thread_num == 0);
#endif
  return retval;
}

int getMaxThreads()
{
  int retval = 1;
#if CH_OPENMP==1
  retval = omp_get_max_threads();
#endif
  return retval;
}

#include "NamespaceFooter.H"
