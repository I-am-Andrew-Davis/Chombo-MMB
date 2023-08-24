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
 * \file threadsEx.cpp
 *
 * \brief This demo focuses on ThreadArchitect
 *
 * This demo focuses on creating a for cycle, closing the cycle, and waiting for
 * incomplete tasks when there is a cycle.
 *
 *//*+*************************************************************************/

#include "ThreadArchitect.H"
#include "ThreadParameters.H"
#include <cstring>
#include <sstream>
#include <mutex>
#ifdef USE_MPI
#include <mpi.h>
#endif

const std::string exeName = "threadsEx";

static const char *const usage =
  "Usage ./threadsEx [-nt x] [-v y]\n"
  "  x   : number of threads.\n"
  "  y   : verbosity\n\n";

using ThreadTools::Verb;
using ThreadTools::ThreadArchitect;

//--Arguments to be parsed from argv
struct Args
{
  int numThreads = -1;
  Verb verbosity = Verb::vW;
};

//--Forward declarations

/// Parses the commandline arguments
Args parseArgs(int argc, const char* argv[]);

/// The function executed by the worker threads
void* workerFunc(void* args);

/// The struct passed as an argument to the worker thread function
struct WorkerArgs
{
  int value = 0;
};

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

  ThreadTools::setVerbosity(args.verbosity);

#ifdef CH_MPI
  const int requested = MPI_THREAD_MULTIPLE;
  int provided;

  int numProcs = -1;
  int procID = -1;

  MPI_Init_thread(&argc, const_cast<char***>(&argv), requested, &provided);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &procID);

  const bool masterProc = (procID == 0);
  if (provided < requested && masterProc)
    {
      std::stringstream ss;
      ss << "Error: this test must be run with MPI thread level "
         << "support of at least " << requested << "!\n"
         << "MPI only provided " << provided << "\n"
         << "Check your mpi implementation.\n";
      if (ThreadTools::getVerbosity() >= ThreadTools::Verb::v0)
        {
          ThreadTools::s_thrOut << ss.rdbuf() << std::flush;
        }
      MPI_Finalize();
      return 1;
    }
#endif

  // Create the architect that is used as the interface for building worker
  // threads
  ThreadArchitect architect;

  // Set the number of threads. This will also figure out the hardware each
  // thread is binded to
  architect.setupTopology(args.numThreads,
                          0,
                          0,
                          0,
                          true,
                          true);

  // Get the number of workers from the architect
  const unsigned numThreads = architect.getThreadTopology()->numWorkers();

  // Build each worker with its function and argument
  std::vector<WorkerArgs> wargs(numThreads);
  for (unsigned tid = 0; tid < numThreads; ++tid)
  {
    wargs[tid].value = tid;
    architect.buildNextWorker(&workerFunc, &wargs[tid]);
  }

  // Wait for each worker to be done
  architect.waitToFinish();

#ifdef USE_MPI
  MPI_Finalize();
#endif
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
  while (argc > iargc && argv[iargc][0] == '-')
    {
      std::cout << argv[iargc] << std::endl;
      if (std::strcmp(argv[iargc], "-nt") == 0)
        {
          args.numThreads = std::atoi(argv[iargc+1]);
          iargc += 2;
        }
      else if (std::strcmp(argv[iargc], "-v") == 0)
        {
          args.verbosity = ThreadTools::intToVerb(std::atoi(argv[iargc+1]));
          iargc += 2;
        }
      else
        {
          ++iargc;
        }
    }
  return args;
}

/*--------------------------------------------------------------------*/
//  The function executed by the worker thread
/** \param[in]  a_args  Casted to a WorkerArgs*
 *//*-----------------------------------------------------------------*/

void* workerFunc(void* a_args)
{
  WorkerArgs* wargs = static_cast<WorkerArgs*>(a_args);

  std::stringstream ss;
  ss << "In task:\n"
     << "  value:  " << wargs->value << "\n";

  std::mutex mtx;
  std::unique_lock<std::mutex> lck(mtx);
  std::cout << ss.rdbuf();
  return nullptr;
}