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
 * \file subTaskTeamsEx.cpp
 *
 * \brief This demo focuses on using sub task teams.
 *
 * This demo focuses on creating a ThreadTeamSubTaskArchitect, thread teams to
 * execute tasks made up of subtasks, and the tasks made up of subtasks.
 *
 *//*+*************************************************************************/

#include "ThreadTeamSubTaskArchitect.H"
#include "ThreadParameters.H"
#include <cstring>
#include <sstream>
#include <mutex>

const std::string exeName = "subTaskTeamsEx";

static const char *const usage =
  "Usage ./subTaskTeamsEx [OPTION]\n"
  "Demos the use of thread teams with subtasks\n"
  "\n"
  "  -h, --help                  print help\n"
  "  -p, --numThreads THREADS    number of threads.\n"
  "  -S, --teamSize TEAMSIZE     team size.\n"
  "  -t, --numTasks TASKS        number of tasks\n"
  "  -s, --numSubTasks SUBTASKS  number of subtasks\n"
  "  -ht, --hyperThreading       turn on hyper threading\n"
  "  -v, --verbosity             verbosity level\n"
  "\n";

using namespace ThreadTools;

//--Aliases
using ThreadTeamSubTasks_t = ThreadTeamSubTasks<size_t, size_t, size_t>;
using SubTaskFunc_t = ThreadTeamSubTasks_t::SubTaskFunc_t;
using SubTasksArg_t = ThreadTeamSubTasks_t::SubTasksArg_t;

//--Arguments to be parsed from argv
struct Args
{
  int numThreads = -1;
  int teamSize = 2;
  Verb verbosity = Verb::vW;
  size_t extraStackSize_MiB = 0;
  bool useHyperThreading = false;
  size_t numTasks = 1;
  size_t numSubtasksPerTask = 2;
};

//--Forward declarations

/// Parses the commandline arguments
Args parseArgs(int argc, const char* argv[]);

/// The subtask function
void SubTask(size_t task, size_t subtask, size_t cnt);

/*----------------------------------------------------------------------------*/

int main(int argc, const char* argv[])
{
  // Parse arguments
  if (argc > 1 && (std::strcmp(argv[1], "-h") == 0 ||
                   std::strcmp(argv[1], "--help") == 0))
    {
      std::cout << usage;
      return 0;
    }

  Args args = parseArgs(argc, argv);

  auto numTasks = args.numTasks;
  auto numSubtasksPerTask = args.numSubtasksPerTask;

  ThreadTools::setVerbosity(args.verbosity);

  // Create the architect that is used as the interface for building threads,
  // building teams, adding tasks made up of subtasks, and waiting for tasks
  // to complete
  ThreadTeamSubTaskArchitect<ThreadTeamSubTasks_t> architect;

  // Set the number of threads and teamsize. This will also figure out the
  //  hardware each thread is binded to
  architect.setupTopology(args.numThreads,
                          args.teamSize,
                          0,
                          0,
                          args.extraStackSize_MiB,
                          args.useHyperThreading,
                          true);

  // Build each thread team
  for (unsigned teamCnt = 0; teamCnt < architect.numThreadTeams(); teamCnt++)
    {
      // The team
      std::shared_ptr<ThreadTeamSubTasks_t> teamSubTasks =
        std::make_shared<ThreadTeamSubTasks_t>(architect.numThreadsPerTeam());
      // The argument used for the function for each thread of the team
      SubTasksArg_t arg{teamSubTasks};
      // Build each thread of the team
      for (unsigned threadCnt = 0; threadCnt < architect.numThreadsPerTeam();
            threadCnt++)
        {
          architect.buildNextWorker(
            ThreadTeamSubTasks_t::threadTeamFunction, &arg);
        }
      // Add the thread team to the architect
      architect.addTeamSubTasks(teamSubTasks);
    }

  const unsigned numThreads = architect.getThreadTopology()->numWorkers();
  const unsigned numThreadsPerTeam = architect.numThreadsPerTeam();
  const unsigned numThreadTeams = architect.numThreadTeams();

  std::cout << "Number of threads: " << numThreads << std::endl
            << "Threads per team:  " << numThreadsPerTeam << std::endl
            << "Number of teams:   " << numThreadTeams << std::endl
            << "Number of tasks:   " << numTasks << std::endl
            << "Subtasks per task: " << numSubtasksPerTask << std::endl;

  // Add the tasks
  unsigned cntSubtask = 0;
  for (unsigned task = 0; task < numTasks; ++task)
    {
      // Build the subtasks of a task
      std::vector<std::function<SubTaskFunc_t>>
        subTasks{numSubtasksPerTask, std::function<SubTaskFunc_t>(SubTask)};

      // Build a vector of tuples to hold the arguments of each subtask
      std::vector<ThreadTeamSubTasks_t::TupleArgs> subTaskArgs;
      for (size_t subtask = 0; subtask < numSubtasksPerTask; ++subtask)
        {
          subTaskArgs.emplace_back(task, subtask, cntSubtask++);
        }

      // Add the task to the architect, so a team can execute it
      architect.addNewTask(subTasks, subTaskArgs);
    }

  // Wait for all the tasks to be finished
  architect.waitToFinish();
}

/*--------------------------------------------------------------------*/
//  Parses the commandline arguments
/** \param[in]  argc    argc passed to main
 *  \param[in]  argc    argv passed to main
 *//*-----------------------------------------------------------------*/

Args parseArgs(int argc, const char* argv[])
{
  Args args;

  int iargc = 1;
  while (argc > iargc)
    {
      if (std::strcmp(argv[iargc], "-p") == 0)
        {
          args.numThreads = std::atoi(argv[iargc+1]);
          iargc += 2;
        }
      else if (std::strcmp(argv[iargc], "--teamSize") == 0 ||
          std::strcmp(argv[iargc], "-S") == 0)
        {
          args.teamSize = std::atoi(argv[iargc+1]);
          iargc += 2;
        }
      else if (std::strcmp(argv[iargc], "--numTasks") == 0 ||
          std::strcmp(argv[iargc], "-t") == 0)
        {
          args.numTasks = std::atoi(argv[iargc+1]);
          iargc += 2;
        }
      else if (std::strcmp(argv[iargc], "--numSubTasks")== 0 ||
          std::strcmp(argv[iargc], "-s") == 0)
        {
          args.numSubtasksPerTask = std::atoi(argv[iargc+1]);
          iargc += 2;
        }
      else if (std::strcmp(argv[iargc], "-v") == 0 ||
               std::strcmp(argv[iargc], "--verbosity") == 0)
        {
          args.verbosity = ThreadTools::intToVerb(std::atoi(argv[iargc+1]));
          iargc += 2;
        }
      else if (std::strcmp(argv[iargc], "-ht") == 0 ||
               std::strcmp(argv[iargc], "--hyperThreading") == 0)
        {
          args.useHyperThreading = true;
          ++iargc;
        }
      else
        {
          const std::string arg(argv[iargc]);
          throw std::invalid_argument("option is not recognized: " + arg);
        }
    }
  return args;
}

/*--------------------------------------------------------------------*/
//  The sub task function
/** \param[in]  a_task  The id of the task this subtask belongs to
 *  \param[in]  a_subTask
 *                      The id of this subtask within its task
 *  \param[in]  a_cnt   The total id of this subtask
 *//*-----------------------------------------------------------------*/

void SubTask(const size_t a_task, const size_t a_subTask, const size_t a_cnt)
{
  std::stringstream ss;
  ss << "In task:\n"
     << "  Thread id:  " << ThreadTools::tls_threadId << "\n"
     << "  Task:       " << a_task << "\n"
     << "  Subtask:    " << a_subTask << "\n"
     << "  Count:      " << a_cnt << "\n";

  using namespace std::chrono_literals;
  std::this_thread::sleep_for(200ms);

  std::mutex mtx;
  std::unique_lock<std::mutex> lck(mtx);
  std::cout << ss.rdbuf();
}
