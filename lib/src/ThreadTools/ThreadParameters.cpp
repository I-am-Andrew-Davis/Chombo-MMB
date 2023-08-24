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
 *  \file ThreadParameters.cpp
 *
 *  \brief Declare variables that are extern in the ThreadParameters.H
 *
 *//*+*************************************************************************/


//----- Standard Library -----//

//----- Internal -----//

#include "ThreadParameters.H"

#include "NamespaceHeader.H"

namespace ThreadTools
{
#ifndef FULLOPTIMIZEDRELEASE
Verb s_verbosity = {Verb::vW};
#endif

std::ostream s_thrOut{std::cout.rdbuf()};
std::ostream s_thrErr{std::cerr.rdbuf()};
thread_local int tls_threadId = 0;
} // namespace Threads

#include "NamespaceFooter.H"
