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
#include "HashIVSLayout.H"
#ifdef CH_MPI
#include <mpi.h>
#endif

#include "UsingNamespace.H"

/// Prototypes:

void
parseTestOptions( int argc ,char* argv[] ) ;

int
testHashIVSLayout1();

/// Global variables for handling output:
static const char *pgmname = "testHashIVSLayout" ;
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
  int ret = testHashIVSLayout1() ;

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
 * Routine testHashIVSLayout1
 *
 ******************************************************************************/

int
testHashIVSLayout1()
{
  int status = 0;

//--Experimentation

  using HashLayout =
    HashIVS::Layout<unsigned, int, unsigned,
                    CH_Hash::DynamicTable,
                    HashIVS::AuxInfoDefault,
                    8, 16, 16>;
  HashLayout ivl;
  ivl |= HashLayout::iixvec_type{ 1, 1, 1 };
  ivl |= HashLayout::iixvec_type{ 4, 1, 1 };
  ivl |= HashLayout::iixvec_type{ 0, 1, 1 };
  ivl |= HashLayout::iixvec_type{ -1, 0, 0 };
  ivl.forEachIV([]
                (const HashLayout::iixvec_type& a_iv)
                  {
                    std::cout << a_iv << std::endl;
                  });
  std::cout << "Sz : " << ivl.hashmap_size() << std::endl;
  std::cout << "BC : " << ivl.bucket_count() << std::endl;
  // std::cout << "LF : " << ivl.load_factor() << std::endl;
  // std::cout << "MLF: " << ivl.max_load_factor() << std::endl;

  // BS bs;
  // bs[0] = 1;
  // std::cout << bs.count() << std::endl;
  // std::cout << bs.rank << std::endl;
  // std::cout << bs.size() << std::endl;
  // std::cout << bs.NW << std::endl;
  // std::cout << bs.dimensions() << std::endl;
  // std::cout << "tallyb     : " << bs.c_tallyb << std::endl;
  // std::cout << "dimsW      : " << bs.c_dimsW << std::endl;
  // std::cout << "dimsWin0   : " << bs.c_dimsWin0 << std::endl;
  // std::cout << "strideWallb: " << bs.c_strideWallb << std::endl;
  // std::cout << "strideW    : " << bs.c_strideW << std::endl;
  // std::cout << "strideWin0 : " << bs.c_strideWin0 << std::endl;

//--Set 1

  // {
  // }

//--Set 2

  {
    int stat2 = 0;
    {
      using HashLayout =
        HashIVS::Layout<unsigned, int, unsigned,
                        CH_Hash::DynamicTable,
                        HashIVS::AuxInfoDefault,
                        16, 16, 16>;
      HashLayout ivl;
      Box box(-IntVect_unit, 3*IntVect_unit);
      ivl |= box;
      int checks = box.numPts();
      int cnt = 0;
      ivl.forEachIV([&box, &checks, &cnt]
                    (const HashLayout::iixvec_type& a_iv)
                      {
                        ++cnt;
                        const IntVect iv(a_iv);
                        if (box.contains(iv)) --checks;
                        // std::cout << a_iv << ' ' << iv << std::endl;
                      });
      stat2 += checks;
      if (cnt != box.numPts()) ++stat2;
    }
    if (stat2 != 0)
      {
        pout() << "Failure in set 2" << endl;
        status += stat2;
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
