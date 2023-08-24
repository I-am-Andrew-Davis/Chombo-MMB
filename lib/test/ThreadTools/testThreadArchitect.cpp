#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cstring>
#include <atomic>

#include "parstream.H"
#include "ThreadArchitect.H"
#include "ThreadBuilder.H"
#include "CH_Timer.H"
#ifdef CH_MPI
#include <mpi.h>
#endif
#include "UsingNamespace.H"

using namespace ThreadTools;





/// Prototypes:

void
parseTestOptions( int argc ,char* argv[] ) ;

int
testThreadArchitect();

/// Global variables for handling output:
static const char *pgmname = "testThreadArchitect" ;
static const char *indent = "   ", *indent2 = "      " ;
static bool verbose = true ;

/// Code:

int
main(int argc, char* argv[])
{
  CH_TIMERS("main()");
  CH_TIMER("wholething", t1);
  CH_START(t1);

#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  parseTestOptions( argc ,argv ) ;

    if ( verbose )
      pout() << indent2 << "Beginning " << pgmname << " ..." << endl ;

  ///
  // Run the tests
  ///
  int ret = testThreadArchitect() ;
#ifdef CH_MPI
  MPI_Finalize();
#endif

  CH_STOP(t1);
  CH_TIMER_REPORT();

  return ret;
}

void* func1(void* a_count)
{
  int* count = (int*) a_count;
  (*count)++;
  return nullptr;
}

void* waitForTeam(void* a_count)
{
  std::atomic<int>* count = (std::atomic<int>*) a_count;
  count->fetch_add(1);

  return nullptr;
}

int
testThreadArchitect()
{
  int return_code = 0;

  // build bound worker
  {
    ThreadArchitect architect;
    if (architect.getThreadTopology() != nullptr)
      {
        if (verbose)
          {
            pout() << indent2 << "ThreadTopology should not be defined until it is setup" << endl;
          }
        return_code = -1;
      }

    std::shared_ptr<const HardwareTopology> hwTopology = std::make_shared<const HardwareTopology>();

    unsigned numThreads = hwTopology->numCores() / 2;
    if (numThreads == 0) numThreads++;

    architect.setupTopology(numThreads,
                            0,
                            0,
                            0,
                            true,
                            false);

    if (architect.getThreadTopology()->numWorkers() != numThreads)
      {
        if (verbose)
          {
            pout() << indent2 << "The number of workers does not equal the number requested" << endl;
          }
        return_code = -1;
      }

    int count = 0;
    
    for (unsigned worker = 0; worker < architect.getThreadTopology()->numWorkers(); ++worker)
      {
        architect.buildNextWorker(&func1, &count);
      }
    architect.waitToFinish();

    if (count != static_cast<int>(architect.getThreadTopology()->numWorkers()))
      {
        if (verbose)
          {
            pout() << indent2 << "Bound worker threads were not successfully built and executed" << endl;
          }
        return_code = -1;
      }
  }

  // add non-bound worker
  {
    ThreadArchitect architect;
    std::shared_ptr<const HardwareTopology> hwTopology = architect.getHardwareTopology();

    unsigned numThreads = hwTopology->numCores() / 2;
    if (numThreads == 0)
      {
        numThreads++;
      }

    architect.setupTopology(numThreads,
                            0,
                            0,
                            0,
                            false,
                            false);

    ThreadBuilder builder(hwTopology);

    int count = 0;
    
    for (unsigned worker = 0; worker < architect.getThreadTopology()->numWorkers(); ++worker)
      {
        architect.buildNextWorker(&func1, &count);
      }
    architect.waitToFinish();
    if (count != static_cast<int>(architect.getThreadTopology()->numWorkers()))
      {
        if (verbose)
          {
            pout() << indent2 << "Non-bound worker threads were not successfully built and executed" << endl;
          }
        return_code = -1;
      }
  }

  if (return_code == 0)
    {
      pout() << indent << "hardware topology passed." << endl;
    }
  else
    {
      pout() << indent << "hardware topology FAILED!!!" << endl;
    }

  return return_code;
}

///
// Parse the standard test options (-v -q) out of the command line.
// Stop parsing when a non-option argument is found.
///
void
parseTestOptions( int argc ,char* argv[] )
{
  for ( int i = 1 ; i < argc ; ++i )
  {
    if ( argv[i][0] == '-' ) //if it is an option
    {
      // compare 3 chars to differentiate -x from -xx
      if ( strncmp( argv[i] ,"-v" ,3 ) == 0 )
      {
        verbose = true ;
        // argv[i] = "" ;
      }
      else if ( strncmp( argv[i] ,"-q" ,3 ) == 0 )
      {
        verbose = false ;
        // argv[i] = "" ;
      }
      else
      {
        break ;
      }
    }
  }
  return ;
}
