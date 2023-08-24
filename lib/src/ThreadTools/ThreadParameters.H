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
 *  \file ThreadParameters.H
 *
 *  \brief Utility functions and parameters for namespace Threads
 *
 *//*+*************************************************************************/

#ifndef _THREADPARAMETERS_H_
#define _THREADPARAMETERS_H_

//----- Standard Library -----//

#include <iostream>
#include <cassert>

#include "NamespaceHeader.H"

namespace ThreadTools
{

#ifdef NDEBUG
  #define TH_assert(cond) (void)0
#else
  #define TH_assert(cond) assert(cond)
#endif
  
  /// the verbosity enum
  enum class Verb
  {
    vE = -1,
    vW = 0,
    v0 = vW,
    v1 = 1,
    v2 = 2,
    v3 = 3,
    v4 = 4
  };

#ifdef FULLOPTIMIZEDRELEASE
  /// the verbosity: cannot be changed
  static constexpr Verb s_verbosity = {Verb::vW};
#else
  /// the verbosity: can be changed
  extern Verb s_verbosity;
#endif

/*--------------------------------------------------------------------*/
//  Set the verbosity if not FULLOPTIMIZEDRELEASE
/** \param[in]  a_verb    The verbosity
 *//*-----------------------------------------------------------------*/
 
  inline
  void setVerbosity(const Verb a_verb)
  {
#ifndef FULLOPTIMIZEDRELEASE
    s_verbosity = a_verb;
#endif
  }

/*--------------------------------------------------------------------*/
//  Get the verbosity
/** \return              The verbosity
 *//*-----------------------------------------------------------------*/
 
  inline
  Verb getVerbosity()
  {
    return s_verbosity;
  }

/*--------------------------------------------------------------------*/
//  Convert an int to a verbosity
/** \param[in]  a_val    The value of the int
 *  \return              The converted value
 *//*-----------------------------------------------------------------*/
 
  inline
  Verb intToVerb(const int a_val)
  {
    Verb verb = Verb::vE;
    if (a_val < 0)
      {
        verb = Verb::vE;
      }
    else if (a_val == 0)
      {
        verb = Verb::vW;
      }
    else if (a_val == 1)
      {
        verb = Verb::v1;
      }
    else if (a_val == 2)
      {
        verb = Verb::v2;
      }
    else if (a_val == 3)
      {
        verb = Verb::v3;
      }
    else if (a_val >= 4)
      {
        verb = Verb::v4;
      }
    return verb;
  }

  extern thread_local int tls_threadId;
  extern std::ostream s_thrOut;
  extern std::ostream s_thrErr;

/*--------------------------------------------------------------------*/
//  Set the output
/** \param[in]  a_sb    The streambuf of the output
 *//*-----------------------------------------------------------------*/
 
  inline
  void setThrOut(std::streambuf* a_sb)
  {
    s_thrOut.rdbuf(a_sb);
  }

/*--------------------------------------------------------------------*/
//  Set the error output
/** \param[in]  a_sb    The streambuf of the error output
 *//*-----------------------------------------------------------------*/

  inline
  void setThrErr(std::streambuf* a_sb)
  {
    s_thrErr.rdbuf(a_sb);
  }
}

#include "NamespaceFooter.H"

#endif  /* ! defined _THREADPARAMETERS_H_ */