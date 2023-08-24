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
#include "FArrayBox.H"
#include "CHArray.H"
#include "BaseFabMacros.H"
#ifdef CH_MPI
#include <mpi.h>
#endif
#include "UsingNamespace.H"

/// Prototypes:

void
parseTestOptions( int argc ,char* argv[] ) ;

int
testMDArrayMacros();

/// Global variables for handling output:
static const char *pgmname = "testMDArrayMacros" ;
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
  int ret = testMDArrayMacros() ;

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
 * Routine testMDArrayMacros
 *
 ******************************************************************************/

int
testMDArrayMacros()
{
  int errors = 0;

/*--------------------------------------------------------------------*
 * Parameters
 *--------------------------------------------------------------------*/

  const int n = 4;

/*--------------------------------------------------------------------*
 * Test MDArray access
 *--------------------------------------------------------------------*/

//--BaseFab

  {
    Box box(IntVect_unit, n*IntVect_unit);
    FArrayBox fab(box, 2);

    {  // Using BaseFab directly
      for (int iC = 0; iC != 2; ++iC)
        {
          MD_BOXLOOP(box, i)
            {
              if (&(fab[MD_IX(i, iC)]) !=
                  &(fab(IntVect(D_DECL6(i0, i1, i2, i3, i4, i5)), iC)))
                {
                  ++errors;
                }
            }
        }
    }

    {  // Using MD macros
      MD_ARRAY(arr, fab);
      for (int iC = 0; iC != 2; ++iC)
        {
          MD_BOXLOOP(box, i)
            {
              if (&(arr[MD_IX(i, 0)]) !=
                  &(fab(IntVect(D_DECL6(i0, i1, i2, i3, i4, i5)), 0)))
                {
                  ++errors;
                }
            }
        }
    }
  }

  {
    Box box(-n*IntVect_unit, -IntVect_unit);
    FArrayBox fab(box, 2);

    {  // Using BaseFab directly
      for (int iC = 0; iC != 2; ++iC)
        {
          MD_BOXLOOP(box, i)
            {
              if (&(fab[MD_IX(i, iC)]) !=
                  &(fab(IntVect(D_DECL6(i0, i1, i2, i3, i4, i5)), iC)))
                {
                  ++errors;
                }
            }
        }
    }

    {  // Using MD macros
      MD_ARRAY(arr, fab);
      for (int iC = 0; iC != 2; ++iC)
        {
          MD_BOXLOOP(box, i)
            {
              if (&(arr[MD_IX(i, 0)]) !=
                  &(fab(IntVect(D_DECL6(i0, i1, i2, i3, i4, i5)), 0)))
                {
                  ++errors;
                }
            }
        }
    }
  }

  {
    Box box(-IntVect_unit, (n-2)*IntVect_unit);
    FArrayBox fab(box, 2);

    {  // Using BaseFab directly
      for (int iC = 0; iC != 2; ++iC)
        {
          MD_BOXLOOP(box, i)
            {
              if (&(fab[MD_IX(i, iC)]) !=
                  &(fab(IntVect(D_DECL6(i0, i1, i2, i3, i4, i5)), iC)))
                {
                  ++errors;
                }
            }
        }
    }

    {  // Using MD macros
      MD_ARRAY(arr, fab);
      for (int iC = 0; iC != 2; ++iC)
        {
          MD_BOXLOOP(box, i)
            {
              if (&(arr[MD_IX(i, 0)]) !=
                  &(fab(IntVect(D_DECL6(i0, i1, i2, i3, i4, i5)), 0)))
                {
                  ++errors;
                }
            }
        }
    }
  }

//--CHArray

  {
    Box box(IntVect_unit, n*IntVect_unit);
    CHArray<Real, SpaceDim+1, ArRangeCol> charrayx(box, 2);

    {  // Regular
      MD_ARRAY_CHARRAY(arr, charrayx);
      for (int iC = 0; iC != 2; ++iC)
        {
          MD_BOXLOOP(box, i)
            {
              if (&(arr[MD_IX(i, iC)]) !=
                  &(charrayx(IntVect(D_DECL6(i0, i1, i2, i3, i4, i5)), iC)))
                {
                  ++errors;
                }
            }
        }
    }
  }

  {
    Box box(-n*IntVect_unit, -IntVect_unit);
    CHArray<Real, SpaceDim+1, ArRangeCol> charrayx(box, 2);

    {  // Regular
      MD_ARRAY_CHARRAY(arr, charrayx);
      for (int iC = 0; iC != 2; ++iC)
        {
          MD_BOXLOOP(box, i)
            {
              if (&(arr[MD_IX(i, iC)]) !=
                  &(charrayx(IntVect(D_DECL6(i0, i1, i2, i3, i4, i5)), iC)))
                {
                  ++errors;
                }
            }
        }
    }
  }

  {
    Box box(-IntVect_unit, (n-2)*IntVect_unit);
    CHArray<Real, SpaceDim+1, ArRangeCol> charrayx(box, 2);

    {  // Regular
      MD_ARRAY_CHARRAY(arr, charrayx);
      for (int iC = 0; iC != 2; ++iC)
        {
          MD_BOXLOOP(box, i)
            {
              if (&(arr[MD_IX(i, iC)]) !=
                  &(charrayx(IntVect(D_DECL6(i0, i1, i2, i3, i4, i5)), iC)))
                {
                  ++errors;
                }
            }
        }
    }
  }

  return errors;
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
