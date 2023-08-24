#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>
#include <cstring>
using std::endl;
#include "parstream.H"
#include "IntVect.H"
#include "ProblemDomain.H"
#include "DisjointBoxLayout.H"
#include "LoadBalance.H"
#include "LevelData.H"
#include "FluxBox.H"
#include "DebugOut.H"
#ifdef CH_MPI
#include <mpi.h>
#endif

#include "UsingNamespace.H"

/// Prototypes:

void
parseTestOptions( int argc ,char* argv[] ) ;

int
testSandbox();

/// Global variables for handling output:
static const char *pgmname = "testSandbox" ;
static const char *indent = "   ", *indent2 = "      " ;
static bool verbose = true ;

/// Code:

int
main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  parseTestOptions( argc ,argv ) ;

  if (verbose)
    pout() << indent2 << "Beginning " << pgmname << " ..." << endl ;

  ///
  // Run the tests
  ///
  int ret = testSandbox() ;

  if (ret == 0)
    {
      pout() << indent << pgmname << " passed all tests" << endl ;
    }
  else
    {
      pout() << indent << pgmname << " failed " << ret << " test(s)" << endl ;
    }

#ifdef CH_MPI
  MPI_Finalize();
#endif
  return ret;
}


/*******************************************************************************
 *
 * Routine testCHArray
 *
 ******************************************************************************/

int
testSandbox()
{
  int status = 0;

  bool periodic[SpaceDim] = { false, false };
  ProblemDomain domain(Box(IntVect_zero, 7*IntVect_unit), periodic);

  Vector<Box> newBoxes;
  newBoxes.push_back(Box(IntVect{0, 0}, IntVect{3, 3}));
  newBoxes.push_back(Box(IntVect{4, 0}, IntVect{7, 3}));
  newBoxes.push_back(Box(IntVect{0, 4}, IntVect{3, 7}));
  newBoxes.push_back(Box(IntVect{4, 4}, IntVect{7, 7}));
  Vector<int> procIDs;
  LoadBalance(procIDs, newBoxes);
  DisjointBoxLayout dbl(newBoxes, procIDs, domain);
  LevelData<FluxBox> lvl(dbl, 1, 2*IntVect_unit);
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      std::cout << "INTCODE: " << dit().intCode() << std::endl;
      const int dir = 0;
      // for (const int dir : EachDir)
        {
          FArrayBox& fab = lvl[dit][dir];
          fab.setVal(SpaceDim*procID() + dit().intCode() + 0.5*dir);
          std::cout << dir << ' ' << fab.box() << std::endl;
          dumpFAB2DSlicePretty(&fab, 0);
        }
    }
  lvl.exchange();
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
    {
      std::cout << "INTCODE: " << dit().intCode() << std::endl;
      const int dir = 0;
      // for (const int dir : EachDir)
        {
          FArrayBox& fab = lvl[dit][dir];
          std::cout << dir << ' ' << fab.box() << std::endl;
          dumpFAB2DSlicePretty(&fab, 0);
        }
    }
  return status;
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
      if (argv[i][0] == '-') //if it is an option
        {
          // compare 3 chars to differentiate -x from -xx
          if (strncmp( argv[i] ,"-v" ,3 ) == 0)
            {
              verbose = true ;
              // argv[i] = "" ;
            }
          else if (strncmp( argv[i] ,"-q" ,3 ) == 0)
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
