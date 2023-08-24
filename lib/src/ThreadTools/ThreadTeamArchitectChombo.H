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
 *  \file ThreadTeamArchitectChombo.H
 *
 *  \brief includes and aliases for Chombo
 *
 *//*+*************************************************************************/

#ifndef _THREADTEAMARCHITECTCHOMBO_H_
#define _THREADTEAMARCHITECTCHOMBO_H_

//----- Standard Library -----//

//----- Internal -----//

#include "ThreadTeamArchitect.H"
#include "ThreadTeamSubTasks.H"
#include "ThreadTeamSubTaskArchitect.H"

#include "NamespaceHeader.H"

//-- Forward declarations from Chombo for aliases

class DataIndex;
class Box;

namespace ThreadTools
{

//-- Public aliases

/// A container for the arguments of a subtask in Chombo
using ChomboSubTaskArgs_t = std::tuple<DataIndex, Box, Box>;

/// The thread team for subtasks in Chombo
using ChomboSubTasks_t = ThreadTeamSubTasks<DataIndex, Box, Box>;

/// The thread team architect in Chombo
using ThreadTeamArchitectChombo_t = ThreadTeamSubTaskArchitect<ChomboSubTasks_t>;

}

#include "NamespaceFooter.H"

#endif  /* ! defined _THREADTEAMARCHITECTCHOMBO_H_ */
