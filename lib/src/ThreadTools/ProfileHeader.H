#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

/******************************************************************************/
/**
 * \file ProfileHeader.H
 *
 * \brief
 *
 *//*+*************************************************************************/

#ifndef _PROFILERHEADER_H_
#define _PROFILERHEADER_H_

#include "NamespaceHeader.H"

namespace Profiling
{

struct ProfileHeader
{
  char hr_id[128];                    ///< 128 bytes to identify the
                                      ///< application "uniquely"
  int processID;                      ///< Rank of the process that generated
                                      ///< this profile
  int threadID;                       ///< The ID of the thread being profiled
  int numProcs;                       ///< The number of processes being
                                      ///< profiled
  int bufferSize;                     ///< Size of profile buffers used
  char startTime[40];                 ///< The start time of the profile
};

}

#include "NamespaceFooter.H"

#endif  /* ! defined _PROFILERHEADER_H_ */
