#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// Purpose:
//   Test the basic ideas of using a ProblemDomain per block.  This test only
//   examines the connected capability of ProblemDomain

//#include <cstring>
//#include <set>
//#include <map>
//#include <functional>
//#include <utility>
//#include <unordered_map>

#ifdef CH_MPI
#include <mpi.h>
#endif

#include "parstream.H"
#include "ProblemDomain.H"

#include "UsingNamespace.H"

/// Prototypes:

void
parseTestOptions(int argc, char* argv[]);

int
testBlockDomain();

/// Global variables for handling output:
static const char *pgmname = "testBlockDomain" ;
static const char *indent = "   ", *indent2 = "      ";
static bool verbose = true;

/// Code:

int
main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  parseTestOptions(argc, argv);

  if (verbose)
    {
      pout() << indent2 << "Beginning " << pgmname << " ..." << std::endl;
    }

  ///
  // Run the tests
  ///
  int ret = testBlockDomain();

  if (ret == 0)
    {
      pout() << indent << pgmname << " passed all tests" << std::endl;
    }
  else
    {
      pout() << indent << pgmname << " failed " << ret << " test(s)"
             << std::endl;
    }

#ifdef CH_MPI
  MPI_Finalize();
#endif
  return ret;
}


int testBlockDomain()
{
  int status = 0;

  stc::Vector<bool, SpaceDim> allPeriodic(true);

  Box box(IntVect_zero, 3*IntVect_unit);
  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      stc::Vector<bool, SpaceDim> tanPeriodic(allPeriodic);
      tanPeriodic[dir] = false;
      // a) Just connected on one side
      BlockDomain domainLo(box);
      domainLo.setAsConnected(dir, Side::Lo);
      BlockDomain domainHi(box);
      domainHi.setAsConnected(dir, Side::Hi);
      // b) Connected on one side and periodic in other directions
      BlockDomain domainPtanLo(box, tanPeriodic.dataPtr());
      domainPtanLo.setAsConnected(dir, Side::Lo);
      BlockDomain domainPtanHi(box, tanPeriodic.dataPtr());
      domainPtanHi.setAsConnected(dir, Side::Hi);
      // c) Connected on one side and periodic in all directions
      BlockDomain domainPallLo(box, allPeriodic.dataPtr());
      domainPallLo.setAsConnected(dir, Side::Lo);
      BlockDomain domainPallHi(box, allPeriodic.dataPtr());
      domainPallHi.setAsConnected(dir, Side::Hi);

      // For each set above, test:
// contains, intersects, &=, bdry{Lo|Hi},
      // adjCell{Lo|Hi}
      // 1) contains
      {
        int status1 = 0;

        // 1a)
        Box testBox(box);
        // Grow in various dir
        testBox.growLo(dir, 1);
        if (!domainLo.contains(testBox)) ++status1;
        testBox.growHi(dir, 1);
        if (domainLo.contains(testBox)) ++status1;
        if (domainHi.contains(testBox)) ++status1;
        testBox.growLo(dir, -1);
        if (!domainHi.contains(testBox)) ++status1;

        testBox = grow(box, 1);
        if (domainLo.contains(testBox)) ++status1;
        if (domainHi.contains(testBox)) ++status1;

        // 1b)
        testBox.growHi(dir, -1);
        if (!domainPtanLo.contains(testBox)) ++status1;
        testBox.growHi(dir, 1);
        if (domainPtanLo.contains(testBox)) ++status1;
        if (domainPtanHi.contains(testBox)) ++status1;
        testBox.growLo(dir, -1);
        if (!domainPtanHi.contains(testBox)) ++status1;

        testBox = grow(box, 1);
        if (domainPtanLo.contains(testBox)) ++status1;
        if (domainPtanHi.contains(testBox)) ++status1;

        // 1b) - IntVect
        if (!domainPtanLo.contains(IntVect_zero)) ++status1;
        if (!domainPtanHi.contains(IntVect_zero)) ++status1;
        if (SpaceDim > 2)
          {
            const int tanDir = !dir;
            IntVect iv(box.smallEnd());
            const int shiftBy = 2*box.size(tanDir);
            // Shift left and right in a periodic direction
            iv.shift(tanDir, -shiftBy);
            if (!domainPtanLo.contains(iv)) ++status1;
            if (!domainPtanHi.contains(iv)) ++status1;
            iv.shift(tanDir, 2*shiftBy);
            if (!domainPtanLo.contains(iv)) ++status1;
            if (!domainPtanHi.contains(iv)) ++status1;
            // Shift off the low end in dir
            iv.shift(dir, -1);
            if (!domainPtanLo.contains(iv)) ++status1;
            if ( domainPtanHi.contains(iv)) ++status1;
            iv = box.bigEnd();
            // Shift left in a periodic direction
            iv.shift(tanDir, -shiftBy);
            // Shift off the high end in dir
            iv.shift(dir, 1);
            if ( domainPtanLo.contains(iv)) ++status1;
            if (!domainPtanHi.contains(iv)) ++status1;
          }
        if (!domainPtanLo.contains(IntVect_zero)) ++ status1;

        // 1c)
        if (!domainPallLo.contains(testBox)) ++status1;
        if (!domainPallHi.contains(testBox)) ++status1;

        if (verbose)
          {
            pout() << "contains       (dir=" << dir << ") status: " << status1
                   << std::endl;
          }
        status += status1;
      }

      // 2) intersects
      {
        int status2 = 0;

        // 2a)
        Box testBox(box);
        testBox.grow(dir, 1);
        if (!domainLo.intersects(testBox)) ++status2;
        if (!domainHi.intersects(testBox)) ++status2;
        if ( domainLo.intersects(Box{})) ++status2;
        if ( domainHi.intersects(Box{})) ++status2;
        Box loTestBox = adjCellLo(testBox, dir, 1);
        Box hiTestBox = adjCellHi(testBox, dir, 1);
        if (!domainLo.intersects(loTestBox)) ++status2;
        if ( domainHi.intersects(loTestBox)) ++status2;
        if ( domainLo.intersects(hiTestBox)) ++status2;
        if (!domainHi.intersects(hiTestBox)) ++status2;
        if (!domainLo.intersectsNotEmpty(loTestBox)) ++status2;
        if ( domainHi.intersectsNotEmpty(loTestBox)) ++status2;
        if ( domainLo.intersectsNotEmpty(hiTestBox)) ++status2;
        if (!domainHi.intersectsNotEmpty(hiTestBox)) ++status2;

        // 2b) same tests for tangential periodic
        if (!domainPtanLo.intersects(testBox)) ++status2;
        if (!domainPtanHi.intersects(testBox)) ++status2;
        if ( domainPtanLo.intersects(Box{})) ++status2;
        if ( domainPtanHi.intersects(Box{})) ++status2;
        if (!domainPtanLo.intersects(loTestBox)) ++status2;
        if ( domainPtanHi.intersects(loTestBox)) ++status2;
        if ( domainPtanLo.intersects(hiTestBox)) ++status2;
        if (!domainPtanHi.intersects(hiTestBox)) ++status2;
        if (!domainPtanLo.intersectsNotEmpty(loTestBox)) ++status2;
        if ( domainPtanHi.intersectsNotEmpty(loTestBox)) ++status2;
        if ( domainPtanLo.intersectsNotEmpty(hiTestBox)) ++status2;
        if (!domainPtanHi.intersectsNotEmpty(hiTestBox)) ++status2;

        // 2c) all periodic should always intersect
        if (!domainPallLo.intersects(testBox)) ++status2;
        if (!domainPallHi.intersects(testBox)) ++status2;
        if ( domainPallLo.intersects(Box{})) ++status2;
        if ( domainPallHi.intersects(Box{})) ++status2;
        if (!domainPallLo.intersects(loTestBox)) ++status2;
        if (!domainPallHi.intersects(loTestBox)) ++status2;
        if (!domainPallLo.intersects(hiTestBox)) ++status2;
        if (!domainPallHi.intersects(hiTestBox)) ++status2;
        if (!domainPallLo.intersectsNotEmpty(loTestBox)) ++status2;
        if (!domainPallHi.intersectsNotEmpty(loTestBox)) ++status2;
        if (!domainPallLo.intersectsNotEmpty(hiTestBox)) ++status2;
        if (!domainPallHi.intersectsNotEmpty(hiTestBox)) ++status2;

        if (verbose)
          {
            pout() << "intersects     (dir=" << dir << ") status: " << status2
                   << std::endl;
          }
        status += status2;
      }

      // 3) &=
      {
        int status3 = 0;

        // 3a)
        Box testBox;
        Box growBox;
        testBox = box;
        testBox.growLo(dir, 1);
        growBox = grow(box, 1);
        growBox &= domainLo;
        if (growBox != testBox) ++status3;
        testBox = box;
        testBox.growHi(dir, 1);
        growBox = grow(box, 1);
        growBox &= domainHi;
        if (growBox != testBox) ++status3;

        // 3b)
        testBox = grow(box, 1);
        testBox.growHi(dir, -1);
        growBox = grow(box, 1);
        growBox &= domainPtanLo;
        if (growBox != testBox) ++status3;
        testBox = grow(box, 1);
        testBox.growLo(dir, -1);
        growBox = grow(box, 1);
        growBox &= domainPtanHi;
        if (growBox != testBox) ++status3;

        // 3c)
        testBox = grow(box, 1);
        growBox = grow(box, 1);
        growBox &= domainPallLo;
        if (growBox != testBox) ++status3;
        growBox &= domainPallHi;
        if (growBox != testBox) ++status3;

        if (verbose)
          {
            pout() << "&=             (dir=" << dir << ") status: " << status3
                   << std::endl;
          }
        status += status3;
      }

      // 4) bdry{Lo|Hi}
      {
        int status4 = 0;

        // 4a)
        if (!bdryLo(domainLo, dir, 1).isEmpty()) ++status4;
        if ( bdryHi(domainLo, dir, 1) != bdryHi(box, dir, 1)) ++status4;
        if ( bdryLo(domainHi, dir, 1) != bdryLo(box, dir, 1)) ++status4;
        if (!bdryHi(domainHi, dir, 1).isEmpty()) ++status4;

        // 4b)
        if (!bdryLo(domainPtanLo, dir, 1).isEmpty()) ++status4;
        if ( bdryHi(domainPtanLo, dir, 1) != bdryHi(box, dir, 1)) ++status4;
        if ( bdryLo(domainPtanHi, dir, 1) != bdryLo(box, dir, 1)) ++status4;
        if (!bdryHi(domainPtanHi, dir, 1).isEmpty()) ++status4;

        // 4c)
        if (!bdryLo(domainPallLo, dir, 1).isEmpty()) ++status4;
        if (!bdryHi(domainPallLo, dir, 1).isEmpty()) ++status4;
        if (!bdryLo(domainPallHi, dir, 1).isEmpty()) ++status4;
        if (!bdryHi(domainPallHi, dir, 1).isEmpty()) ++status4;

        if (verbose)
          {
            pout() << "bdry{Lo|Hi}    (dir=" << dir << ") status: " << status4
                   << std::endl;
          }
        status += status4;
      }

      // 5) adjCell{Lo|Hi}
      {
        int status5 = 0;

        // 5a)
        if (!adjCellLo(domainLo, dir, 1).isEmpty()) ++status5;
        if ( adjCellHi(domainLo, dir, 1) != adjCellHi(box, dir, 1)) ++status5;
        if ( adjCellLo(domainHi, dir, 1) != adjCellLo(box, dir, 1)) ++status5;
        if (!adjCellHi(domainHi, dir, 1).isEmpty()) ++status5;

        // 5b)
        if (!adjCellLo(domainPtanLo, dir, 1).isEmpty()) ++status5;
        if ( adjCellHi(domainPtanLo, dir, 1) != adjCellHi(box, dir, 1))
          ++status5;
        if ( adjCellLo(domainPtanHi, dir, 1) != adjCellLo(box, dir, 1))
          ++status5;
        if (!adjCellHi(domainPtanHi, dir, 1).isEmpty()) ++status5;

        // 5c)
        if (!adjCellLo(domainPallLo, dir, 1).isEmpty()) ++status5;
        if (!adjCellHi(domainPallLo, dir, 1).isEmpty()) ++status5;
        if (!adjCellLo(domainPallHi, dir, 1).isEmpty()) ++status5;
        if (!adjCellHi(domainPallHi, dir, 1).isEmpty()) ++status5;

        if (verbose)
          {
            pout() << "adjCell{Lo|Hi} (dir=" << dir << ") status: " << status5
                   << std::endl;
          }
        status += status5;
      }

    }  // Loop over directions

  return status;
}

///
// Parse the standard test options (-v -q) out of the command line.
// Stop parsing when a non-option argument is found.
///
void
parseTestOptions(int argc, char* argv[])
{
  for (int i = 1; i < argc; ++i)
    {
      if (argv[i][0] == '-') //if it is an option
        {
          // compare 3 chars to differentiate -x from -xx
          if (strncmp( argv[i], "-v", 3 ) == 0)
            {
              verbose = true;
              // argv[i] = "";
            }
          else if (strncmp( argv[i], "-q", 3 ) == 0)
            {
              verbose = false;
              // argv[i] = "";
            }
          else
            {
              break;
            }
        }
    }
  return;
}
