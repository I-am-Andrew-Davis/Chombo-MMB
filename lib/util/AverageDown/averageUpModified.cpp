#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// This code takes plotfiles of multilevel, cell-averaged data, 
// one finely resolved, one coarsely resolved, and averages 
// finely-resolved data down to coarsely-resolved data.  Averaging 
// is performed only on the valid regions of each AMR grid.  

// At present, the fine data and coarse data must have equivalent 
// ProblemDomains (i.e. they must occupy the same physical space), 
// and corresponding levels in the different data sets must also 
// occupy the same physical space.

// Perhaps someone else will later enable the code to average 
// down cell-centered, point-valued data.

#include <iostream>
#include <strstream>
using namespace std;

#include "AMRIO.H"
#include "ParmParse.H"
#include "BoxLayout.H"
#include "CoarseAverage.H"
#include "FourthOrderFineInterp.H"
#include "AMRLevel.H"
#include "MappedLevelData.H"
#include "LevelGridMetrics.H"
#include "MultiBlockCoordSys.H"
#include "MultiBlockFluxRegister.H"
#include "LoadBalance.H"
#include "computeSum.H"
#include "computeNorm.H"
#include "FourthOrderUtil.H"
#include "TimeInterpolatorRK4.H"
#include "LevelRK4_v2.H"
#include "MOLUtilities.H"
#include "NodeAMRIO.H"
#include "FABView.H"

void init(string&         a_aveDownRoot,
          string&         a_crseRoot,
          string&         a_fineRoot,
          int&            a_aveDownRatio,
          int&            a_intFieldSize,
          int&            a_numFineStart,
          int&            a_numFineFinish,
          int&            a_fineStep,
          int&            a_fineMult,
          int&            a_crseRef,
          Vector<string>& a_aveDownVars,
          bool&           a_isTimeDep,
          bool&           a_isCheckpoint);

void constructAveDownNames(Vector<string>&       a_aveDownNames,
                           const Vector<string>& a_aveDownVars);

void constructPlotFileName(ostrstream&   a_fileName,
                           const string& a_fileRoot,
                           const int     a_intFieldWidth,
                           const int     a_step);

// Function for MPI
void dumpmemoryatexit();

int main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
  // setChomboMPIErrorHandler();
  MPI_Barrier(Chombo_MPI::comm);  // Barrier #1
#endif
  
  // ------------------------------------------------
  // parse command line and input file
  // ------------------------------------------------
  // infile must be first
  if (argc < 2)
  {
    cerr << "  need inputs file" << endl;
    abort();
  }

  char* in_file = argv[1];
  ParmParse pp(argc-2, argv+2, NULL, in_file);

#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif

  // Set defaults and make declarations
  bool isTimeDep = false;
  
  string aveDownRoot, crseRoot, fineRoot;
  
  Vector<string> aveDownVars;
  
  int interpRatio;
  int intFieldSize, numFineStart, numFineFinish, fineStep, fineMult, crseRef;
  bool isCheckpoint;
  
  init(aveDownRoot, crseRoot, fineRoot, interpRatio, intFieldSize, 
       numFineStart, numFineFinish, fineStep, fineMult, crseRef, aveDownVars, 
       isTimeDep, isCheckpoint);
  
  int nStep;
  
  if (!isTimeDep)
    {
      fineStep = 1;
      numFineFinish = numFineStart;
    }
  
  for (nStep = numFineStart; nStep <= numFineFinish; nStep += fineStep)
    {
      ostrstream fineFile;
      ostrstream coarseFile;
      ostrstream aveDownFile;
      
      fineFile.fill('0');
      coarseFile.fill('0');
      aveDownFile.fill('0');

      if (isTimeDep)
        {
          constructPlotFileName(fineFile, fineRoot, intFieldSize, nStep);
          constructPlotFileName(coarseFile, crseRoot, intFieldSize, crseRef);
          constructPlotFileName(aveDownFile, aveDownRoot, intFieldSize, nStep);
        }
      else
        {
          // if not time dependent, file roots are really filenames
          fineFile << fineRoot << ends;
          coarseFile << crseRoot << ends;
          aveDownFile << aveDownRoot << ends;
        }
      
      //pout() << "fine Filename = " << fineFile.str() << endl;
      //pout() << "coarse Filename = " << crseFile.str() << endl;
      //pout() << "aveDown Filename = " << aveDownFile.str() << endl;

      // declare memory and read in fine data
      Vector<LevelData<FArrayBox>* > coarseSoln;
      Vector<string> coarseVars; // fine solution variable names
      Vector<DisjointBoxLayout> coarseGrids;
      Box coarseDomain;
      Real coarseDx, coarseDt, coarseTime;
      Vector<int> coarseRefRatio;
      int coarseNumLevels;
      IntVect ghostVect = IntVect::Zero;
      string coarseFileName(coarseFile.str());

      int max_level, cur_step;
      Real cur_time;
      Vector<int> regrid_intervals;
      Vector<int> steps_since_regrid;
      Vector<int> isPeriodic;

      // Get all the checkpoint data out of the file
      // ReadAMRHierarchyHDF5Checkpoint(
      //   coarseFileName,
      //   coarseGrids,
      //   coarseSoln,
      //   coarseVars,
      //   coarseDomain,
      //   coarseDx,
      //   coarseDt,
      //   coarseTime,
      //   coarseRefRatio,
      //   coarseNumLevels,
      //   max_level,
      //   cur_step,
      //   cur_time,
      //   regrid_intervals,
      //   steps_since_regrid,
      //   isPeriodic);

      ReadAMRHierarchyHDF5(coarseFileName,
                           coarseGrids,
                           coarseSoln,
                           coarseVars,
                           coarseDomain,
                           coarseDx,
                           coarseDt,
                           coarseTime,
                           coarseRefRatio,
                           coarseNumLevels);

      // declare memory and read in coarse data
      Vector<LevelData<FArrayBox>* > fineSoln;
      Vector<string> fineVars; // coarse soln variable names
      Vector<DisjointBoxLayout> fineGrids;
      Box fineDomain;
      Real fineDx, fineDt, fineTime;
      Vector<int> fineRefRatio;
      int fineNumLevels;
      string fineFileName(fineFile.str());
      int max_level_dummy, cur_step_dummy;
      Real cur_time_dummy;
      Vector<int> regrid_intervals_dummy;
      Vector<int> steps_since_regrid_dummy;
      Vector<int> isPeriodic_dummy;

      // ReadAMRHierarchyHDF5Checkpoint(
      //   fineFileName,
      //   fineGrids,
      //   fineSoln,
      //   fineVars,
      //   fineDomain,
      //   fineDx,
      //   fineDt,
      //   fineTime,
      //   fineRefRatio,
      //   fineNumLevels,
      //   max_level_dummy,
      //   cur_step_dummy,
      //   cur_time_dummy,
      //   regrid_intervals_dummy,
      //   steps_since_regrid_dummy,
      //   isPeriodic_dummy);

      ReadAMRHierarchyHDF5(fineFileName,
                           fineGrids,
                           fineSoln,
                           fineVars,
                           fineDomain,
                           fineDx,
                           fineDt,
                           fineTime,
                           fineRefRatio,
                           fineNumLevels);

      int numFine = fineVars.size();
      int numCrse = coarseVars.size();

      Vector<string> aveDownNames;
      aveDownNames.resize(numCrse);

      constructAveDownNames(aveDownNames, aveDownVars);

      IntVect coarseInterpGhostVect = 2*IntVect::Unit;

      Vector<LevelData<FArrayBox>* > aveDown(fineNumLevels);
      Vector<LevelData<FArrayBox>* > coarseWithGhosts(coarseNumLevels);

      {
        int level = 0;

        ProblemDomain probDomFine;
        ProblemDomain probDomCoarse;
        bool periodicity[SpaceDim];
        for (int i = 0; i != SpaceDim; ++i)
          {
            periodicity[i] = false;
            // if (!isPeriodic[i])
            //   {
            //     periodicity[i] = false;
            //   }
          }
        probDomFine = ProblemDomain(fineDomain, periodicity);
        probDomCoarse = ProblemDomain(coarseDomain, periodicity);

        HDF5Handle handleOne(fineFileName, HDF5Handle::OPEN_RDONLY);
        HDF5HeaderData headerOne;
        headerOne.readFromFile(handleOne);
        char levelStrOne[32];
        sprintf(levelStrOne,"%d",level);
        const std::string labelOne = std::string("level_") + levelStrOne;
        handleOne.setGroup(labelOne);
        Vector<Box> readBoxesOne;
        const int gridStatusOne = read(handleOne, readBoxesOne);
        Vector<int> procIDsOne;
        LoadBalance(procIDsOne, readBoxesOne);

        DisjointBoxLayout fineDBL;
        fineDBL = DisjointBoxLayout(readBoxesOne, procIDsOne, probDomFine);
        fineDBL.close();
        handleOne.close();

        fineGrids[0] = fineDBL;

        HDF5Handle handleTwo(coarseFileName, HDF5Handle::OPEN_RDONLY);
        HDF5HeaderData headerTwo;
        headerTwo.readFromFile(handleTwo);
        char levelStrTwo[32];
        sprintf(levelStrTwo,"%d",level);
        const std::string labelTwo = std::string("level_") + levelStrTwo;
        handleTwo.setGroup(labelTwo);
        Vector<Box> readBoxesTwo;
        const int gridStatusTwo = read(handleTwo, readBoxesTwo);
        Vector<int> procIDsTwo;
        LoadBalance(procIDsTwo, readBoxesTwo);

        DisjointBoxLayout coarseDBL;
        coarseDBL = DisjointBoxLayout(readBoxesTwo, procIDsTwo, probDomCoarse);
        coarseDBL.close();
        handleTwo.close();

        aveDown[level] = new LevelData<FArrayBox>(fineDBL,
                                                  numCrse, 
                                                  ghostVect);
        coarseWithGhosts[level] =
          new LevelData<FArrayBox>(coarseDBL,
                                   numCrse,
                                   coarseInterpGhostVect);
        (*(coarseSoln[level])).copyTo(*(coarseWithGhosts[level]));
        (*(coarseWithGhosts[level])).exchange();

        FourthOrderFineInterp interpOp;
        interpOp.define(fineDBL,
                        numCrse,
                        interpRatio,
                        probDomFine,
                        0);
        interpOp.interpToFine(*(aveDown[level]), *(coarseWithGhosts[level]));
      }

      // Output the plot file
      // We need to fix this so that fine data information is used
      // Fix this as follows:
      //   (1) crseDt is replaced by refineRatio*fineDt
      //   (2) crseTime is replaced by fineTime
      Real interpFineDt = coarseDt/interpRatio;

      // Write all the checkpoint data into the new file
      // WriteAMRHierarchyHDF5Checkpoint(
      //   aveDownFile.str(),
      //   fineGrids,
      //   aveDown,
      //   aveDownNames,
      //   fineDomain,
      //   fineDx,
      //   interpFineDt,
      //   coarseTime,
      //   fineRefRatio,
      //   fineNumLevels,
      //   max_level,
      //   cur_step,
      //   cur_time,
      //   regrid_intervals,
      //   steps_since_regrid,
      //   isPeriodic);

       // WriteAMRHierarchyHDF5(
       //      aveDownFile.str(),
       //      crseGrids,
       //      aveDown,
       //      aveDownNames,
       //      crseDomain,
       //      crseDx,
       //      coarsenedFineDt,
       //      fineTime,
       //      crseRefRatio,
       //      crseNumLevels);

       WriteAMRHierarchyHDF5(
            aveDownFile.str(),
            fineGrids,
            aveDown,
            aveDownNames,
            fineDomain,
            fineDx,
            interpFineDt,
            coarseTime,
            fineRefRatio,
            fineNumLevels);

      // clean up memory
      for (int level = 0; level < fineNumLevels; level++)
        {
          if (fineSoln[level] != NULL)
            {
              delete fineSoln[level];
              fineSoln[level] = NULL;
            }
        }
      
      for (int level = 0; level < coarseNumLevels; level++)
        {
          if (coarseSoln[level] != NULL)
            {
              delete coarseSoln[level];
              coarseSoln[level] = NULL;
            }
          
          if (aveDown[level] != NULL)
            {
              delete aveDown[level];
              aveDown[level] = NULL;
            }
        }
      
    } // end loop over time steps
  
#ifdef CH_MPI
  dumpmemoryatexit();
  MPI_Finalize();
#endif
} // end main


void init(string&         a_aveDownRoot,
          string&         a_crseRoot,
          string&         a_fineRoot,
          int&            a_interpRatio,
          int&            a_intFieldSize,
          int&            a_numFineStart,
          int&            a_numFineFinish,
          int&            a_fineStep,
          int&            a_fineMult,
          int&            a_crseRef,
          Vector<string>& a_aveDownVars,
          bool&           a_isTimeDep,
          bool&           a_isCheckpoint)
{
  ParmParse ppAveDown("avedown");
  
  ppAveDown.get("aveDownRoot", a_aveDownRoot);
  ppAveDown.get("crseRoot", a_crseRoot);
  ppAveDown.get("fineRoot", a_fineRoot);
  ppAveDown.get("interpRatio", a_interpRatio);
  
  int isTimeDepInt = a_isTimeDep;
  ppAveDown.query("isTimeDep", isTimeDepInt);
  a_isTimeDep = (isTimeDepInt == 1);

  int isCheckpoint = a_isCheckpoint;
  ppAveDown.query("isCheckpoint", isCheckpoint);
  a_isCheckpoint = (isCheckpoint == 1);
  
  a_numFineStart = 0;
  ppAveDown.query("numFineStart", a_numFineStart);
  
  if (a_isTimeDep)
    {
      ppAveDown.get("numFineFinish", a_numFineFinish);
      a_fineStep = 1;
      ppAveDown.query("fineStep", a_fineStep);
      ppAveDown.get("mult", a_fineMult);
      a_intFieldSize = 4;
      ppAveDown.query("intFieldSize", a_intFieldSize);
      ppAveDown.query("crseRef", a_crseRef);
    }

  int nVars = ppAveDown.countval("aveDownVars");
  if (nVars > 0)
    {
      a_aveDownVars.resize(nVars);
      ppAveDown.getarr("aveDownVars", a_aveDownVars, 0, nVars);
    }
}

void constructAveDownNames(Vector<string>&       a_aveDownNames,
                           const Vector<string>& a_aveDownVars)
{
  CH_assert(a_aveDownNames.size() == a_aveDownVars.size());

  // for now, don't do anything fancy -- just copy
  for (int i = 0; i < a_aveDownVars.size(); i++)
  {
    a_aveDownNames[i] = a_aveDownVars[i];
  }
}

void constructPlotFileName(ostrstream&   a_fileName,
                           const string& a_fileRoot,
                           const int     a_intFieldWidth,
                           const int     a_step)
{
  // this is kinda klugy, but what are ya gonna do?
  a_fileName << a_fileRoot
             << setw(a_intFieldWidth) << a_step
             << "."
             << setw(1) << CH_SPACEDIM
             << "d.hdf5" << ends;
}
