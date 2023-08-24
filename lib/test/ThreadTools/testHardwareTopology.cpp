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

#include "parstream.H"
#include "HardwareTopology.H"
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
testHardwareTopology();

/// Global variables for handling output:
static const char *pgmname = "testHardwareTopology" ;
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
  int ret = testHardwareTopology() ;
#ifdef CH_MPI
  MPI_Finalize();
#endif

  CH_STOP(t1);
  CH_TIMER_REPORT();

  return ret;
}


int
testHardwareTopology()
{
  int return_code = 0;

  HardwareTopology topology;

  topology.hasNUMA();
  topology.hasHyperThreads();

  if (topology.topology() == nullptr)
    {
      if (verbose)
        {
          pout() << indent2 << "Failure in creating the hwloc_topology_t in HardwareTopology" << endl;
        }
      return_code = -1;
    }

  if (topology.cpuSet() == nullptr)
    {
      if (verbose)
        {
          pout() << indent2 << "Failure in creating the cpuSet in HardwareTopology" << endl;
        }
      return_code = -1;
    }

  if (topology.coreSet() == nullptr)
    {
      if (verbose)
        {
          pout() << indent2 << "Failure in creating the cpuSet in HardwareTopology" << endl;
        }
      return_code = -1;
    }

  if (topology.numCPUs() <= 0)
    {
      if (verbose)
        {
          pout() << indent2 << "The number of CPUs must be greater than 0" << endl;
        }
      return_code = -1;
    }

  if (topology.numLocalCPUs() <= 0)
    {
      if (verbose)
        {
          pout() << indent2 << "The number of CPUs for this process must be greater than 0" << endl;
        }
      return_code = -1;
    }

  if (topology.numPUs() <= 0)
    {
      if (verbose)
        {
          pout() << indent2 << "The number of PUs must be greater than 0" << endl;
        }
      return_code = -1;
    }

  if (topology.numCores() <= 0)
    {
      if (verbose)
        {
          pout() << indent2 << "The number of cores must be greater than 0" << endl;
        }
      return_code = -1;
    }

  if (topology.numLocalCores() <= 0)
    {
      if (verbose)
        {
          pout() << indent2 << "The number of cores for this process must be greater than 0" << endl;
        }
      return_code = -1;
    }

  if (topology.numLocalPUs() <= 0)
    {
      if (verbose)
        {
          pout() << indent2 << "The number of PUs for this process must be greater than 0" << endl;
        }
      return_code = -1;
    }

  if (topology.getPageSize() <= 0)
    {
      if (verbose)
        {
          pout() << indent2 << "The page size must be greater than 0" << endl;
        }
      return_code = -1;
    }


  if (topology.numPUs() < topology.numLocalPUs())
    {
      if (verbose)
        {
          pout() << indent2 << "The total number of PUs must be equal to or greater than"
          " the number of PUs for this process" << endl;
        }
      return_code = -1;
    }

  if (topology.numCPUs() < topology.numLocalCPUs())
    {
      if (verbose)
        {
          pout() << indent2 << "The total number of CPUs must be equal to or greater than"
          " the number of CPUs for this process" << endl;
        }
      return_code = -1;
    }

  if (topology.numCores() < topology.numLocalCores())
    {
      if (verbose)
        {
          pout() << indent2 << "The total number of cores must be equal to or greater than"
          " the number of cores for this process" << endl;
        }
      return_code = -1;
    }


  if (topology.numPUs() < topology.numCPUs())
    {
      if (verbose)
        {
          pout() << indent2 << "The total number of PUs must be equal to or greater than"
          " the number of CPUs" << endl;
        }
      return_code = -1;
    }

  if (topology.numCores() < topology.numCPUs())
    {
      if (verbose)
        {
          pout() << indent2 << "The total number of cores must be equal to or greater than"
          " the number of CPUs" << endl;
        }
      return_code = -1;
    }

  if (topology.numPUs() < topology.numCores())
    {
      if (verbose)
        {
          pout() << indent2 << "The total number of PUs must be equal to or greater than"
          " the number of cores" << endl;
        }
      return_code = -1;
    }

  if (topology.getCPUObj(topology.numCPUs() - 1) == nullptr)
    {
      if (verbose)
        {
          pout() << indent2 << "The CPU object is not defined" << endl;
        }
      return_code = -1;
    }

  if (topology.getPUObj(topology.numPUs() - 1) == nullptr)
    {
      if (verbose)
        {
          pout() << indent2 << "The PU object is not defined" << endl;
        }
      return_code = -1;
    }

  if (topology.getCoreObj(topology.numCores() - 1) == nullptr)
    {
      if (verbose)
        {
          pout() << indent2 << "The core object is not defined" << endl;
        }
      return_code = -1;
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
