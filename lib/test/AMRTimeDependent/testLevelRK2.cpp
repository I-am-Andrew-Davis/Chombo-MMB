#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cmath>
#include <iostream>
using std::endl;

#include "BRMeshRefine.H"
#include "AMRLevel.H"
#include "CH_HDF5.H"
#include "LevelRK2.H"
#include "TimeInterpolatorRK2.H"
#include "parstream.H"

#include "AMRIO.H"
#include "DebugOut.H"

#include "FArrayBox.H"
#include "LevelData.H"
#include "LayoutIterator.H"
#include "AMRLevel.H"
#include "CoarseAverage.H"
#include "FineInterp.H"
#include "BoxIterator.H"
#include "LoadBalance.H"
#include "UsingNamespace.H"

/// Global variables for handling output:
static const char* pgmname = "testLevelRK2" ;
static const char* indent = "   ";
static const char* indent2 = "      " ;
static bool verbose = true ;

class EvalOp {
public:
  void evalRHS(LevelData<FArrayBox>& a_rhs,
               LevelData<FArrayBox>& a_soln,
               const int             a_stage,
               const Real            a_stageTime,
               const Real            a_stageWeight);

  void defineRHSData(LevelData<FArrayBox>&       a_newRHS,
                     const LevelData<FArrayBox>& a_oldRHS);

  void defineSolnData(LevelData<FArrayBox>&       a_newSoln,
                      const LevelData<FArrayBox>& a_oldSoln);

  void updateODE(LevelData<FArrayBox>& a_soln,
                 LevelData<FArrayBox>& a_rhs,
                 Real                  a_dt);

  void copySolnData(LevelData<FArrayBox>& a_dest,
                    LevelData<FArrayBox>& a_src);
};

/// Prototypes:
int
testLevelRK2();

void
parseTestOptions(int argc ,char* argv[]) ;

int
main(int argc ,char* argv[])
{
#ifdef CH_MPI
  MPI_Init (&argc, &argv);
#endif
  parseTestOptions( argc ,argv ) ;
  if ( verbose )
    pout () << indent2 << "Beginning " << pgmname << " ..." << endl ;

  int status = testLevelRK2();

  if ( status == 0 )
    pout() << indent << pgmname << " passed." << endl ;
  else
    pout() << indent << pgmname << " failed with return code " << status << endl ;

#ifdef CH_MPI
  MPI_Finalize ();
#endif
  return status ;
}

int
testLevelRK2 ()
{
  int refRatio = 2;
  int numGhosts = 0;

  Box prob_domain (IntVect::Zero,
                   IntVect::Unit);
  Box refined_domain = refine(prob_domain, refRatio);
  bool periodic[SpaceDim];
  D_TERM6(periodic[0] = true;,
          periodic[1] = true;,
          periodic[2] = true;,
          periodic[3] = true;,
          periodic[4] = true;,
          periodic[5] = true;)
  ProblemDomain fineDomain(refined_domain, periodic);

  Vector<Box> boxes;
  boxes.push_back(prob_domain);
  Vector<Box> refinedBoxes;
  refinedBoxes.push_back(refined_domain);

  Vector<int> procIds;
  LoadBalance(procIds, boxes);
  Vector<int> refinedProcIds;
  LoadBalance(refinedProcIds, refinedBoxes);

  int numComps = 1;
  DisjointBoxLayout dbl(boxes, procIds);
  LevelData<FArrayBox> data(dbl, numComps);
  LevelData<FArrayBox> oldData(dbl, numComps);

  DisjointBoxLayout refinedDbl(refinedBoxes, refinedProcIds);
  LevelData<FArrayBox> refinedData(refinedDbl, numComps);

  Real initialVal = 3.;

  // Initialize
  for(DataIterator dit = oldData.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& fab = oldData[dit];
      fab.setVal(initialVal);
    }
  for(DataIterator dit = refinedData.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& fab = refinedData[dit];
      fab.setVal(initialVal);
    }

  TimeInterpolatorRK2 interp;
  interp.defineCoarseParts(dbl, numComps, numGhosts);
  interp.defineFineParts(refinedDbl, fineDomain, refRatio, numGhosts);

  Real time = 0.;
  Real dt = 1.e-1;
  bool initializeNew = true;

  EvalOp op;

  RK2LevelAdvance<LevelData<FArrayBox>,
                  LevelData<FArrayBox>,
                  TimeInterpolatorRK2,
                  EvalOp>
    (data,
     oldData,
     &interp,
     time,
     dt,
     initializeNew,
     op);

  // Get some intermediate values to test
  LevelData<FArrayBox> interp0(dbl,numComps);
  interp.intermediate(interp0, 0., 0., 0, interp0.interval());

  LevelData<FArrayBox> interp1(dbl,numComps);
  interp.intermediate(interp1, 1., 0., 0, interp0.interval());

  LevelData<FArrayBox> interp2(dbl,numComps);
  interp.intermediate(interp2, 0.5, 0., 0, interp0.interval());

  LevelData<FArrayBox> interp3(dbl,numComps);
  interp.intermediate(interp3, 1./3., 0., 0, interp0.interval());

  LevelData<FArrayBox> interp4(dbl,numComps);
  interp.intermediate(interp4, 7./11., 0., 0, interp0.interval());

  // Now begin testing
  int status = 0;
  for(DataIterator dit = data.dataIterator(); dit.ok(); ++dit)
    {
      const FArrayBox& fab = data[dit];
      const FArrayBox& interp0Fab = interp0[dit];
      const FArrayBox& interp1Fab = interp1[dit];
      const FArrayBox& interp2Fab = interp2[dit];
      const FArrayBox& interp3Fab = interp3[dit];
      const FArrayBox& interp4Fab = interp4[dit];

      const Box& box = fab.box();


      const Real expected = 4.19025;
      const Real interp0Expected = initialVal;
      const Real interp1Expected = 14.9025;
      const Real interp2Expected = 8.225625;
      const Real interp3Expected = 6.3225;
      const Real interp4Expected = 9.902665289;


      MD_BOXLOOP(box, i)
        {
          const Real actual = fab[MD_IX(i, 0)];
          if((actual-expected)>1.e-14)
            {
              status += 1;
              pout() << "Failure in cell: " << MD_GETIV(i)
                     << " expected: " << expected
                     << " actual: " << actual
                     << std::endl;
            }

          Real interp0Actual = interp0Fab[MD_IX(i,0)];
          if((interp0Actual-interp0Expected)>1.e-14)
            {
              status += 1;
              pout() << "Interp0 failure in cell: " << MD_GETIV(i)
                     << " expected: " << interp0Expected
                     << " actual: " << interp0Actual
                     << " difference: " << (interp0Actual-interp0Expected)
                     << std::endl;
            }
          Real interp1Actual = interp1Fab[MD_IX(i,0)];
          if((interp1Actual-interp1Expected)>1.e-14)
            {
              status += 1;
              pout() << "Interp1 failure in cell: " << MD_GETIV(i)
                     << " expected: " << interp1Expected
                     << " actual: " << interp1Actual
                     << " difference: " << (interp1Actual-interp1Expected)
                     << std::endl;
            }
          Real interp2Actual = interp2Fab[MD_IX(i,0)];
          if((interp2Actual-interp2Expected)>1.e-14)
            {
              status += 1;
              pout() << "Interp2 failure in cell: " << MD_GETIV(i)
                     << " expected: " << interp2Expected
                     << " actual: " << interp2Actual
                     << " difference: " << (interp2Actual-interp2Expected)
                     << std::endl;
            }
          Real interp3Actual = interp3Fab[MD_IX(i,0)];
          if((interp3Actual-interp3Expected)>1.e-14)
            {
              status += 1;
              pout() << "Interp3 failure in cell: " << MD_GETIV(i)
                     << " expected: " << interp3Expected
                     << " actual: " << interp3Actual
                     << " difference: " << (interp3Actual-interp3Expected)
                     << std::endl;
            }
          Real interp4Actual = interp4Fab[MD_IX(i,0)];
          if((interp4Actual-interp4Expected)>3.e-10)
            {
              status += 1;
              pout() << "Interp4 failure in cell: " << MD_GETIV(i)
                     << " expected: " << interp4Expected
                     << " actual: " << interp4Actual
                     << " difference: " << (interp4Actual-interp4Expected)
                     << std::endl;
            }

        }
    }

  return status;
}

///
// Parse the standard test options (-v -q) out of the command line.
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
        }
    }
  return ;
}

void EvalOp::evalRHS(LevelData<FArrayBox>& a_rhs,
                     LevelData<FArrayBox>& a_soln,
                     const int             a_stage,
                     const Real            a_stageTime,
                     const Real            a_stageWeight)
{
  for(DataIterator dit = a_rhs.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& rhsFab = a_rhs[dit];
      FArrayBox& solnFab = a_soln[dit];

      rhsFab.copy(solnFab);
      rhsFab.mult(solnFab);
    }
}

void EvalOp::defineRHSData(LevelData<FArrayBox>&       a_newRHS,
                           const LevelData<FArrayBox>& a_oldRHS)
{
  a_newRHS.define(a_oldRHS.getBoxes(), a_oldRHS.nComp(), IntVect::Zero);
}

void EvalOp::defineSolnData(LevelData<FArrayBox>&       a_newSoln,
                            const LevelData<FArrayBox>& a_oldSoln)
{
  a_newSoln.define(a_oldSoln.getBoxes(),
                   a_oldSoln.nComp(),
                   a_oldSoln.ghostVect());
}

void EvalOp::updateODE(LevelData<FArrayBox>& a_soln,
                       LevelData<FArrayBox>& a_rhs,
                       Real                  a_dt)
{
  DisjointBoxLayout dbl = a_soln.getBoxes();
  for(DataIterator dit = a_rhs.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& rhsFab = a_rhs[dit];
      FArrayBox& solnFab = a_soln[dit];
      solnFab.plus(rhsFab, a_dt);
    }
}

void EvalOp::copySolnData(LevelData<FArrayBox>& a_dest,
                          LevelData<FArrayBox>& a_src)
{
  a_src.copyTo(a_dest);
}
