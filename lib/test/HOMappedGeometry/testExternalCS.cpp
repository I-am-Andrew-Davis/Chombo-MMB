#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <fstream>
#include <iostream>
#include <vector>
#include "newMappedGridIO.H"
using std::cout;
using std::endl;
using std::ofstream;


#include "parstream.H"
#include "ParmParse.H"
#include "ExternalCS.H"
#include "WarpedCS.H"
#include "ReadCGNS.H"
#include "FArrayBox.H"
#include "computeNorm.H"
#include "DebugOut.H"
#include "FABView.H"

#include "UsingNamespace.H"

/// Global variables for handling output:
static const char* pgmname = "testExternalCS" ;
static const char* indent = "   ";
static const char* indent2 = "      " ;
static bool verbose = true ;
static bool writePlotFiles = true;

#ifdef CH_USE_DOUBLE
static Real precision = 1.0e-15;
#else
static Real precision = 1.0e-7;
#endif

///
// Parse the standard test options (-v -q -h) and
// app-specific options (-S <domain_size>) out of the command line
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
              //argv[i] = "" ;
            }
          else if ( strncmp( argv[i] ,"-q" ,3 ) == 0 )
            {
              verbose = false ;
              //argv[i] = "" ;
            }
          else if ( strncmp( argv[i] ,"-h" ,3 ) == 0 )
            {
              pout() << "usage: " << pgmname << " [-hqv]" << std::endl ;
              exit( 99 ) ;
            }
        }
    }
  return ;
}

int testExternalCS();

int main(int argc ,char* argv[])
{
#ifdef CH_MPI
  MPI_Init (&argc, &argv);
#endif
  parseTestOptions( argc ,argv ) ;
  if ( verbose )
    pout () << indent2 << "Beginning " << pgmname << " ..." << endl ;

  int status = testExternalCS();

  if ( status == 0 )
    cout << indent << "External coordinate system test" << " passed." << endl ;
  else
    cout << indent << "External coordinate system test" << " failed with return code " << status << endl ;



  if ( status == 0 )
    cout << indent << pgmname << " passed." << endl ;
  else
    cout << indent << pgmname << " failed with return code " << status << endl ;


#ifdef CH_MPI
  CH_TIMER_REPORT();
  MPI_Finalize ();
#endif
  return status ;
}

int testExternalCS()
{
  int returnCode = 0;
#if CH_SPACEDIM == 2
#ifdef CH_USE_CGNS
  int numCells = 64;
  ProblemDomain baseDomain(IntVect::Zero, (numCells-1)*IntVect::Unit);
  int numLevels = 3;
  Vector<int> nRefVect(numLevels-1, 2); // all refinement by 2

  // for the warped grid
  Real physLength = 1.0;
  RealVect domainLength = physLength*RealVect::Unit;
  RealVect dxiVect_warp = (physLength/numCells)*RealVect::Unit;
  RealVect scale_warp = 0.1*RealVect::Unit;

  // for the external grid
  // dXi is is arbitrary here. It is chosen to match the warped grid but ordinarily when given no context it should be 1.
  // For multi-block it must always be 1.
  // RealVect dxiVect_extn = RealVect::Unit; safest
  RealVect dxiVect_extn = dxiVect_warp;
  
  // set numbers of ghost cells
  int numGhost = 1;
  IntVect numGhostVect = numGhost*IntVect::Unit;

  Real L1Crse, L2Crse, maxCrse;
  Real L1volCrse, L2volCrse, maxVolCrse;

  // simplest possible test
  ProblemDomain levelDomain = baseDomain;
  Vector<Box> boxes(1, baseDomain.domainBox());
  Vector<int> procAssign(1,0);

  // Make the coordinate systems
  std::string gridFileName = "warpedGrid64.cgns";
  ReadCGNS readcgns(gridFileName);
  CH_assert(readcgns.numBase() == 1);
  ReadCGNS::ZoneNameMapToBlkIdx zoneNames;
  readcgns.selectBase(1, zoneNames);
  std::string zoneName;
  IntVect zoneNumCell;
  readcgns.readZoneInfo(1, zoneName, zoneNumCell);
  // std::vector<DisjointBoxLayout> grids;
  // std::vector<DisjointBoxLayout> exactCS;
  // for (int level = 0; level<numLevels; level++)
  //   {
  //     DisjointBoxLayout levelGrids(boxes, procAssign, levelDomain);

  //     WarpedCS levelExactCS(dxVect, origin, stretch);
  //     ExternalCS levelExternalCS(dxVect, origin, stretch);
  //   }

  // Run test on the coordinate systems
  // This a simple series of test that serve to make sure the externalCS functionality is not blatantly wrong
  // FIXME :: Should add a convergence test to show the ExternalCS approaches the WarpedCS solution
  for (int level = 0; level<numLevels; level++)
    {
      DisjointBoxLayout levelGrids(boxes, procAssign, levelDomain);

      WarpedCS levelExactCS(dxiVect_warp,
                            scale_warp,
                            domainLength);
      ExternalCS levelExternCS(readcgns,
                               dxiVect_extn,
                               RealVect::Unit,
                               levelDomain,
                               gridFileName,
                               zoneName); // last argument is unused now

      // Check mapping locations, at nodes
      // When nodes line up with points specified from the CGNS file, mapping
      // should be exact, but elsewhere mapping error is expected
      {
        LevelData<FluxBox> exactLocs(levelGrids, SpaceDim, numGhostVect);
        LevelData<FluxBox> externLocs(levelGrids, SpaceDim, numGhostVect);
        LevelData<FluxBox> locError(levelGrids, SpaceDim, numGhostVect);

        Real maxErr = 0;
        DataIterator dit = levelGrids.dataIterator();
        for (dit.begin(); dit.ok(); ++dit)
          {
            FluxBox& thisError = locError[dit];
            thisError.setVal(0.0);
            for (const auto faceDir : EachDir)
              {
                FArrayBox& exactFab = exactLocs[dit][faceDir];
                FArrayBox& externFab = externLocs[dit][faceDir];
                FArrayBox& errorFab = locError[dit][faceDir];
                MD_BOXLOOP(exactFab.box(), i)
                  {
                    IntVect iv = MD_GETIV(i);
                    RealVect xi_exact = iv*dxiVect_warp;
                    RealVect xi_extrn = iv*dxiVect_extn;
                    auto xExact = levelExactCS.realCoord(xi_exact);
                    auto xExtern = levelExternCS.realCoord(xi_extrn);
                    for (int comp = 0; comp != SpaceDim; comp++)
                      {
                        exactFab[MD_IX(i, comp)] = xExact[comp];
                        externFab[MD_IX(i, comp)] = xExtern[comp];
                        errorFab[MD_IX(i, comp)] = xExact[comp] - xExtern[comp];
                        maxErr = std::max(maxErr,
                                          std::abs(errorFab[MD_IX(i, comp)]));
                      }
                  }
                //dumpFAB(&errorFab);
              }
            if (verbose)
              {
                pout() << "Locations: level " << level << " err = " << maxErr << std::endl;
              }
            if ((level == 0) && (maxErr > std::sqrt(precision)))
              {
                // fail
                returnCode += 1;
                if (verbose)
                  {
                    pout() << "Node locations at defined grid points do not match"
                           << " -- Max(error) for grid locations = "
                           << maxErr << endl;
                    //dumpLDFLoc(&externLocs);
                  }

              }

          }
      }

      // Check metrics, at nodes
      {
        int numMetrics = SpaceDim*SpaceDim;
        LevelData<FluxBox> exactMetrics(levelGrids, numMetrics, numGhostVect);
        LevelData<FluxBox> externMetrics(levelGrids, numMetrics, numGhostVect);
        LevelData<FluxBox> metricError(levelGrids, numMetrics, numGhostVect);

        Real maxErr = 0;
        DataIterator dit = levelGrids.dataIterator();
        for (dit.begin(); dit.ok(); ++dit)
          {
            FluxBox& thisError = metricError[dit];
            thisError.setVal(0.0);
            for (const auto faceDir : EachDir)
              {
                FArrayBox& exactFab = exactMetrics[dit][faceDir];
                FArrayBox& externFab = externMetrics[dit][faceDir];
                FArrayBox& errorFab = metricError[dit][faceDir];
                MD_BOXLOOP(exactFab.box(), i)
                  {
                    IntVect iv = MD_GETIV(i);
                    RealVect xi_exact = iv*dxiVect_warp;
                    RealVect xi_extrn = iv*dxiVect_extn;
                    for (int comp = 0; comp != numMetrics; comp++)
                      {
                        exactFab[MD_IX(i, comp)] =
                          levelExactCS.dXdXi(xi_exact,
                                             (comp % SpaceDim),
                                             (comp / SpaceDim));
                        externFab[MD_IX(i, comp)] =
                          levelExternCS.dXdXi(xi_extrn,
                                              (comp % SpaceDim),
                                              (comp / SpaceDim));
                        errorFab[MD_IX(i, comp)] =
                          exactFab[MD_IX(i, comp)] - externFab[MD_IX(i, comp)];
                        maxErr = std::max(maxErr,
                                          std::abs(errorFab[MD_IX(i, comp)]));
                      }
                  }
                //dumpFAB(&exactFab);
              }

            if (verbose)
              {
                pout() << "Metrics: level " << level << " err = " << maxErr << std::endl;
              }
            if ((level == 0) && (maxErr > std::sqrt(std::sqrt(precision))))
              {
                // fail
                returnCode += 1;
                if (verbose)
                  {
                    pout() << "Metrics at defined grid points do not match"
                           << " -- Max(error) for grid metrics = "
                           << maxErr << endl;
                    //dumpLDFLoc(&externLocs);
                  }

              }
          }
      }

      // Check the mapping inverse, at nodes


      // Do refinement
      if (level < nRefVect.size())
        {
          levelDomain.refine(nRefVect[level]);
          for (int i=0; i<boxes.size(); i++)
            {
              boxes[i].refine(nRefVect[level]);
            }
          dxiVect_warp /= nRefVect[level];
          dxiVect_extn /= nRefVect[level];
        }

    } // end loop over levels
#endif // CH_USE_CGNS
#endif // DIM == 2
  return returnCode;
}

