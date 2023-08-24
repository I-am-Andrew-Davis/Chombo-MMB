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
 * \file postProcess.cpp
 *
 * \brief Post processing utility
 *
 *//*+*************************************************************************/

//----- Standard Library -----//

#include <iostream>
#include <fstream>
#include <strstream>

//----- Chombo Library -----//

#include "ParmParse.H"
#include "AMR.H"
#include "AMRLevel.H"
#include "CONSTANTS.H"
#include "FABView.H"

//----- Internal -----//

#include "AllLevelsOp.H"
#include "SpectralUtil.H"


/*******************************************************************************
 *
 * Begin Post Processing
 *
 ******************************************************************************/

int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  MPI_Init(&a_argc,&a_argv);
#endif

  // Parallel setup
  // FIXME: currently setup for serial runs only

/*==============================================================================
 * User Input
 *============================================================================*/

  // Check for input file
  char* inFile = NULL;
  if (a_argc > 1)
    {
      inFile = a_argv[1];
    }
  else
    {
      std::cout << "Usage: <executable name> <inputfile>" << std::endl;
      std::cout << "No input file specified" << std::endl;
      return 1;
    }

  // Parse the command line and the input file
  ParmParse pp(a_argc-2, a_argv+2, NULL, inFile);

  ParmParse ppFile("postProcess");

  // This determines the amount of diagnositic output generated
  int verbosity = 0;
  ppFile.query("verbosity", verbosity);
  if (verbosity < 0)
    {
      std::cout << "\n" << "Warning:" << "\n" << std::endl;
      std::cout << "Input: 'verbosity' must be >= 0!" << std::endl;
      std::cout << "Setting verbosity to 0\n" << std::endl;
      verbosity = 0;
    }
  // Print a header when the program starts up
  std::cout << "\n\n" << std::endl;
  std::cout << "         "
            << "**** Post Processing for Turbulent Chemistry Simulations ****"
            << "\n\n\n" << std::endl;

  // -------------------------------- //
  //         File Information         //
  // -------------------------------- //

  std::string inputFilePath;
  ppFile.query("input_file_path", inputFilePath);
  if (verbosity > 0)
    {
      std::cout << "Input file path: ";
      std::cout << inputFilePath << std::endl;
    }

  std::string outputFilePath;
  ppFile.query("output_file_path", outputFilePath);
  if (verbosity > 0)
    {
      std::cout << "Output file path: ";
      std::cout << outputFilePath << std::endl;
    }

  std::string inputFileRoot;
  ppFile.query("input_file_root", inputFileRoot);
  if (verbosity > 0)
    {
      std::cout << "Input file root: ";
      std::cout << inputFileRoot << std::endl;
    }

  std::string outputCSVFile;
  ppFile.query("output_CSV_file", outputCSVFile);
  if (verbosity > 0)
    {
      std::cout << "Output CSV file root:";
      std::cout << outputCSVFile << std::endl;
    }

  std::string outputHDFFileRoot;
  ppFile.query("output_HDF_file_root", outputHDFFileRoot);
  if (verbosity > 0)
    {
      std::cout << "Output HDF5 file root:";
      std::cout << outputHDFFileRoot << std::endl;
    }

  bool computeESD = false;
  std::string outputESDFileRoot;
  if (ppFile.contains("output_ESD_file_root"))
    {
      computeESD = true;
      ppFile.get("output_ESD_file_root", outputESDFileRoot);
    }
  if (verbosity > 0)
    {
      std::cout << "Output ESD file root:";
      std::cout << outputESDFileRoot << std::endl;
    }

  int intFieldWidth = 6;
  ppFile.query("fileNumber_Field_Width", intFieldWidth);
  if (intFieldWidth >= 16)
    {
      std::cout << "\n" << "Warning:" << "\n" << std::endl;
      std::cout << "Input: Are you really sure you need "
                << "a quadrillion time steps?" << std::endl;
      std::cout << "Setting fileNumber_Field_Width = 15\n" << std::endl;
      intFieldWidth = 15;
    }
  if (verbosity > 0)
    {
      std::cout << "File number integer field width:";
      std::cout << intFieldWidth << std::endl;
    }

  int fieldCompare = 1;
  for (int i = 0; i != intFieldWidth; ++i)
    {
      fieldCompare *= 10;
    }
  if (verbosity > 1)
    {
      std::cout << "Field compare - debugging:";
                   std::cout << fieldCompare << std::endl;
    }

  int fileStart = 0;
  ppFile.query("file_Start", fileStart);
  if (fileStart < 0)
    {
      std::cout << "\n" << "Warning:" << "\n" << std::endl;
      std::cout << "Input: 'file_Start' must be >= 0" << std::endl;
      std::cout << "Setting file_Start to 0\n" << std::endl;
      fileStart = 0;
    }
  if (fileStart >= fieldCompare)
    {
      std::cout << "\n" << "Error:" << "\n" << std::endl;
      std::cout << "Input: 'file_Start' exceeds size allowed by "
                << "fileNumber_Field_Width" << std::endl;
      std::cout << "Correct the discrepancy and restart" << std::endl;
      return 1;
    }
  if (verbosity > 0)
    {
      std::cout << "File start number:";
      std::cout << fileStart << std::endl;
    }

  int fileStop = 0;
  ppFile.query("file_Stop", fileStop);
  if (fileStop < fileStart)
    {
      std::cout << "\n" << "Warning:" << "\n" << std::endl;
      std::cout << "Input: 'file_Stop' must be >= fileStart" << std::endl;
      std::cout << "Setting file_Stop = fileStart\n" << std::endl;
      fileStop = fileStart;
    }
  if (fileStop >= fieldCompare)
    {
      std::cout << "\n" << "Error:" << "\n" << std::endl;
      std::cout << "Input: 'file_Stop' exceeds size allowed by "
                << "fileNumber_Field_Width" << std::endl;
      std::cout << "Correct the discrepancy and restart" << std::endl;
      return 1;
    }
  if (verbosity > 0)
    {
      std::cout << "File stop number:";
      std::cout << fileStop << std::endl;
    }

  int fileStep = 0;
  ppFile.query("file_Step", fileStep);
  if (fileStep < 1)
    {
      std::cout << "\n" << "Warning:" << "\n" << std::endl;
      std::cout << "Input: 'file_Step' must be >= 1" << std::endl;
      std::cout << "Setting file_Step = 1\n" << std::endl;
      fileStep = 1;
    }
  if (fileStep >= fieldCompare)
    {
      std::cout << "\n" << "Error:" << "\n" << std::endl;
      std::cout << "Input: 'file_Step' exceeds size allowed by "
                << "fileNumber_Field_Width" << std::endl;
      std::cout << "Correct the discrepancy and restart" << std::endl;
      return 1;
    }
  if (verbosity > 0)
    {
      std::cout << "File step size:";
      std::cout << fileStep << std::endl;
    }
  int fileCompare = (fileStop - fileStart) % fileStep;
  if (verbosity > 1)
    {
      std::cout << "File compare - debugging:";
      std::cout << fileCompare << std::endl;
    }
  if (fileCompare != 0)
    {
      std::cout << "\n" << "Warning:" << "\n" << std::endl;
      std::cout << "Input: Note that (file_Stop - file_Start) is "
                << "not evenly divided by file_Step" << std::endl;
      std::cout << "The last file will not be processed\n" << std::endl;
    }

  int numOutGhostCells = 0;
  ppFile.query("num_ghost_cell_output", numOutGhostCells);
  if (numOutGhostCells < 0)
    {
      std::cout << "\n" << "Warning:" << "\n" << std::endl;
      std::cout << "Input: 'num_ghost_cell_output' must be >= 0" << std::endl;
      std::cout << "Setting num_ghost_cell_output = 0\n" << std::endl;
      numOutGhostCells = 0;
    }
  if (verbosity > 0)
    {
      std::cout << "Number of output ghost cells:";
      std::cout << numOutGhostCells << std::endl;
    }

  // Suggested warning: check for the files in a directory -- if they
  //                    don't exist, output a warning

  // -------------------------------- //
  //    Data Scaling/Manipulation     //
  // -------------------------------- //

  Real timeMult = 1.;
  ppFile.query("time_Mult_Scaling", timeMult);
  if (timeMult <= 0.)
    {
      std::cout << "\n" << "Warning:" << "\n" << std::endl;
      std::cout << "Input: 'time_Mult_Scaling' must be > 0" << std::endl;
      std::cout << "Setting time_Mult_Scaling = 1\n" << std::endl;
      timeMult = 1.;
    }
  if (verbosity > 0)
    {
      std::cout << "Time multiplication factor:";
      std::cout << timeMult << std::endl;
    }

  std::vector<Real> varScale;
  int numVarScale = ppFile.countval("var_scaling");
  varScale.resize(numVarScale);
  ppFile.getarr("var_scaling", varScale, 0, numVarScale);

  // Suggested warning: numVarScale should be connected
  //                    with number of outputs somehow

  // -------------------------------- //
  //   Mapping - FIXME: Implement!    //
  // -------------------------------- //

  bool mapped = false;
  ppFile.query("mapped_simulation", mapped);
  if (verbosity > 0)
    {
      std::cout << "Mapped simulation bool:";
      std::cout << mapped << std::endl;
    }

  std::string inputMappingFileRoot;
  ppFile.query("mapping_file_root", inputMappingFileRoot);
  if (verbosity > 0)
    {
      std::cout << "Mapping file root name:";
      std::cout << inputMappingFileRoot << std::endl;
    }

  // -------------------------------- //
  //        Domain Information        //
  // -------------------------------- //

  // (a) type of bc (for spatially dependent variables, e.g. y+)
  //     additionally, symmetry boundaries can have data "reflected"
  //     into ghost cells rather than exchanged (as in periodic)

  std::vector<int> periodicity(SpaceDim, 1);
  ppFile.getarr("periodicity", periodicity, 0, SpaceDim);

  std::vector<std::string> BCs(SpaceDim*2);
  int numBCs = ppFile.countval("boundary_conditions");
  if (numBCs > 0)
    {
      BCs.resize(numBCs);
      ppFile.getarr("boundary_conditions", BCs, 0, numBCs);
    }
  if (numBCs != (2*SpaceDim))
    {
      std::cout << "\n" << "Warning:" << "\n" << std::endl;
      std::cout << "Input: 'boundary_conditions' must have "
                << "2*(space dimensions) number of elements\n" << std::endl;
    }

  // Suggested warning: BCs must match a known, simple set of BCs
  //                    such as periodic (no boundary), symmetry
  //                    no-slip wall, slip wall, and inflow/outflow
  //                    (each list member has a unique property)

  // Suggested warning: number of BCs and number of periodicity
  //                    should match 2*SpaceDim and SpaceDim

  // -------------------------------- //
  //      Variables to Calculate      //
  // -------------------------------- //

  // (a) method of computing variables (points, interp, polynomial)

  int compOrderOfAccuracy = 4;
  ppFile.query("computing_order_of_accuracy", compOrderOfAccuracy);
  if (compOrderOfAccuracy <= 0)
    {
      std::cout << "\n" << "Warning:" << "\n" << std::endl;
      std::cout << "Input: 'computing_order_of_accuracy' must be >= 0" << std::endl;
      std::cout << "Setting computing_order_of_accuracy = 4\n" << std::endl;
      compOrderOfAccuracy = 4;
    }
  if (verbosity > 0)
    {
      std::cout << "Order of accuracy of derivative computations: ";
      std::cout << compOrderOfAccuracy << std::endl;
    }

  int deconvOrderOfAccuracy = 4;
  ppFile.query("computing_order_of_accuracy", deconvOrderOfAccuracy);
  if (deconvOrderOfAccuracy <= 0)
    {
      std::cout << "\n" << "Warning:" << "\n" << std::endl;
      std::cout << "Input: 'computing_order_of_accuracy' must be >= 0" << std::endl;
      std::cout << "Setting computing_order_of_accuracy = 4\n" << std::endl;
      deconvOrderOfAccuracy = 4;
    }
  if (verbosity > 0)
    {
      std::cout << "Order of accuracy of deconvolutions/convolutions: ";
      std::cout << deconvOrderOfAccuracy << std::endl;
    }

  // Suggested warning: number of outputVars should be limited

  int method = 0;
  ppFile.query("compute_variable_method", method);
  if (verbosity > 0)
    {
      std::cout << "Variable computation method: ";
      std::cout << method << std::endl;
    }

  // Suggested warning: if method is not one implemented, set to default

  std::vector<std::string> outputTimePointVars;
  int numCompVars = ppFile.countval("output_time_point_variables");
  if (numCompVars > 0)
    {
      outputTimePointVars.resize(numCompVars);
      ppFile.getarr("output_time_point_variables",
                    outputTimePointVars,
                    0,
                    numCompVars);
    }
  if (verbosity > 0)
    {
    }

  std::vector<std::string> outputTimeSumVars;
  int numCompTimeSumVars = ppFile.countval("output_time_domainSum_variables");
  if (numCompTimeSumVars > 0)
    {
      outputTimeSumVars.resize(numCompTimeSumVars);
      ppFile.getarr("output_time_domainSum_variables",
                    outputTimeSumVars,
                    0,
                    numCompTimeSumVars);
    }
  if (verbosity > 0)
    {
    }

  std::vector<std::string> outputGlobalVars;
  int numCompGlobalVars = ppFile.countval("output_global_variables");
  if (numCompGlobalVars > 0)
    {
      outputGlobalVars.resize(numCompGlobalVars);
      ppFile.getarr("output_global_variables",
                    outputGlobalVars,
                    0,
                    numCompGlobalVars);
    }
  if (verbosity > 0)
    {
    }

  // -------------------------------- //
  //       Nature of Simulation       //
  // -------------------------------- //

  bool compressible = true;
  ppFile.query("compressible", compressible);
  if (verbosity > 0)
    {
    }

  bool thermallyPerfect = false;
  ppFile.query("thermally_perfect", thermallyPerfect);
  if (verbosity > 0)
    {
    }

/*==============================================================================
 * End of Input: Begin Interpreting
 *============================================================================*/

  // Quantities to compute
  int derivative = 0;
  int numOutComps = 0;
  int derivOne = 0;
  for (int i = 0; i != numCompVars; ++i)
    {
      if (outputTimePointVars[i] == "enstrophy")
        {
          derivOne += 1;
        }
      else if (outputTimePointVars[i] == "KE")
        {
        }
    }
  int derivTwo = 0;
  for (int i = 0; i != numCompTimeSumVars; ++i)
    {
      if (outputTimePointVars[i] == "enstrophy")
        {
          derivTwo += 1;
        }
      else if (outputTimePointVars[i] == "KE")
        {
        }
    }
  int derivThree = 0;
  for (int i = 0; i != numCompGlobalVars; ++i)
    {
      if (outputTimePointVars[i] == "enstrophy")
        {
          derivThree += 1;
          numOutComps += 1;
        }
      else if (outputTimePointVars[i] == "KE")
        {
          numOutComps += 1;
        }
    }
  derivative = std::max(derivOne, std::max(derivTwo, derivThree));

  // What does incompressible mean for computations?
  if (!compressible)
    {
      // Density = 1
      // no energy equation
      // pressure?
    }

  // How large do individual stencil operations have to be?
  int derivStencil = (compOrderOfAccuracy+1)/2;
  int convStencil = ((deconvOrderOfAccuracy-1)/2);

  // What do various outputs mean for computations?
  //   --  have a pre-defined list of generic output variables
  //       - pressure, temp, vorticity, entropy, KE
  //   --  have a pre-defined list of more variables
  //       - enstrophy, palinstrophy, helicity
  //       - shear stress, pressure dilatation, dissipation
  //       - non-dimensional values (y+, z+, y*, etc.)
  //       - gradients, divergence, curl, cross-product
  //   --  have a simple scheme to interpret new input
  //       - using functional inputs
  //       - e.g. mag(curl(gradient(dot(divergence(velocity),pressure))))

  // How far out do we need the primitive point values?
  int varPrimOffset = derivStencil*derivative;
  // How far out do we need the cell averaged variables for output?
  int varCellAvg = numOutGhostCells;
  // How far out do we need the point variables for output?
  int varCellPnt = varCellAvg + convStencil;
  int primCellPnt = varCellPnt + varPrimOffset;
  int consCellAvg = primCellPnt + convStencil;

  // FIXME: numCompGhostCells may need to be larger depending on computation
  int numCompGhostCells = consCellAvg;

  // Note: don't do more work than necessary
  // FIXME: build in checks so that
  //       if a variable doesn't need exchanged or copied, don't
  //       if a variable doesn't need deconvolved, don't
  //       if a variable doesn't need consToPrim'ed, don't

/*==============================================================================
 * File Name Setup
 *============================================================================*/

  // Step interval assumed to remain constant, but last file may not be a
  // complete interval -- truncate the last file if necessary
  int numInterval = fileStop - fileStart;
  int numSteps = numInterval/fileStep;
  int numFiles = numSteps + 1;

  // Open the CSV file here
  std::ofstream outputCSV(outputCSVFile,
                          std::ofstream::out | std::ofstream::app);
  // Print the CSV header here -- x is used for time header
  outputCSV << "x";
  for (int i = 0; i != numCompTimeSumVars; ++i)
    {
      outputCSV << ", " << outputTimeSumVars[i];
    }
  outputCSV << "\n";

  // Energy spectral density
  EnergySpectralDensity esd;
  if (computeESD) esd.use();
#if 0  /* If enabling test, also define DEBUG_ESD in SpectralUtil.cpp */
  {
    esd.use();
    constexpr int n = 8;
    const Real L  = 2*Pi;
    const Real dL = L/n;
    Box domainBox(IntVect::Zero, (n-1)*IntVect::Unit);
    ProblemDomain testDomain(domainBox);
    esd.defineSerial(testDomain);
    FArrayBox testFab(domainBox, 3);
    for (BoxIterator bit(domainBox); bit.ok(); ++bit)
      {
        const IntVect& iv = bit();
        const Real x = iv[0]*dL;
        const Real y = iv[1]*dL;
        const Real z = iv[2]*dL;
        testFab(iv, 0) = (2.0) + (1.0*  std::cos(  x) +
                                  1.0*  std::cos(  y) +
                                  1.0*  std::cos(  x)*std::cos(  y) +  // +r
                                  1.0*  std::cos(2*x)*std::sin(  y) +  //  i
                                  1.0*  std::sin(  x)*std::sin(2*y) +  // -r
                                  0.5*  std::cos(2*x) +
                                  0.5*  std::cos(2*y) +
                                  0.5*  std::cos(2*x)*std::cos(2*y) +
                                  0.25* std::cos(3*x) +
                                  0.25* std::cos(3*y) +
                                  0.25* std::cos(3*x)*std::cos(3*y));
        testFab(iv, 1) = (2.0) + (1.0*  std::cos(  x) +
                                  1.0*  std::cos(  y) +
                                  1.0*  std::cos(  x)*std::cos(  y) +
                                  0.5*  std::cos(2*x) +
                                  0.5*  std::cos(2*y) +
                                  0.5*  std::cos(2*x)*std::cos(2*y) +
                                  0.25* std::cos(3*x) +
                                  0.25* std::cos(3*y) +
                                  0.25* std::cos(3*x)*std::cos(3*y));
        testFab(iv, 2) = (0.0);
      }
    esd.load(domainBox, testFab, Interval(0, 2));
    std::vector<Real> esdvec;
    RealVect physDomain(L*RealVect::Unit);
    esd.compute(esdvec, physDomain);
    std::cout << "esd:" << std::endl;
    for (int i = 0; i != n/2+1; ++i)
      {
        std::cout << i << ": " << esdvec[i] << std::endl;
      }
    std::cout << "surplus: " << esdvec[n/2+1] << std::endl;
  }
  exit(1);
#endif

  //**FIXME make sure paths exist

  // Begin the main computation loop
  for (int num = 0; num != numFiles; ++num)
    {
      // Input HDF5 file
      int currentNum = fileStart + (num*fileStep);
      std::ostringstream fileString;
      fileString.fill('0');
      fileString << inputFilePath << inputFileRoot
                 << setw(intFieldWidth) << currentNum << "."
                 << setw(1) << SpaceDim << "d.hdf5";
      std::string inFileName(fileString.str());
      // Output HDF5 file
      std::ostringstream outFileString;
      outFileString.fill('0');
      outFileString << outputFilePath << outputHDFFileRoot
                    << setw(intFieldWidth) << currentNum << "."
                    << setw(1) << SpaceDim << "d.hdf5";
      std::string outFileName(outFileString.str());
      // Output ESD CSV file
      outFileString.str("");
      outFileString << outputFilePath << outputESDFileRoot
                    << setw(intFieldWidth) << currentNum << "."
                    << setw(1) << SpaceDim << "d.csv";
      esd.outputFileName(outFileString.str());

/*==============================================================================
 * Compute Quantities
 *============================================================================*/

      // Note: there are several ways to do this
      // FIXME: only the first method is currently used
      //  (1) point-based method: <U>  -->  U  -->  W  -->  F(W)  -->  <F(W)>
      //  (2) interpolation: obtain finer cell-averaged values and proceed
      //      (a) assume smooth, polynomial basis for interpolation
      //      (b) pde constrained interpolation, avoid smoothness assumption
      //  (3) polynomial: construct in-cell polynomial and proceed analytically

      // FIXME: add functionality to
      //        (a) average-down plot files to a coarser solution using
      //            (i)   coarser plot files with comparable domains
      //            (ii)  a single coarser plot file with a comparable domain
      //            (iii) no comparable plot files
      //        (b) average-down checkpoint files to a coarser solution using
      //            (i)   pre-existing coarser checkpoint files
      //            (ii)  pre-existing coarser plot files
      //            (iii) no pre-existing coarser file
      //        (c) build checkpoint files from plot files and mapping
      //        (d) add levels to checkpoint files over a given region
      //        (e) interpolate plot files to a finer solution

      //**FIXME move this outside the file loop.
      // Start with the level operator
      AllLevelsOp levelObj(outputTimePointVars,
                           outputTimeSumVars,
                           outputGlobalVars,
                           BCs,
                           varScale,
                           timeMult,
                           method,
                           numOutComps,
                           deconvOrderOfAccuracy,
                           compOrderOfAccuracy,
                           numCompGhostCells,
                           numOutGhostCells,
                           verbosity,
                           compressible,
                           thermallyPerfect,
                           mapped);

      // FIXME: timeBaseLevel actually takes the time from the top level
      Real timeBaseLevel = 0;
      //**FIXME a member function name like loadFile would make more sense
      levelObj.define(inFileName, periodicity, timeBaseLevel);

      // FIXME: domSums only contains sums from the coarsest level
      std::vector<Real> domSums;
      domSums.resize(outputTimeSumVars.size(), 0);
      levelObj.runLevels(domSums, esd);

      // Write to HDF5 file
      levelObj.writeLevels(outFileName);

      // Write to CSV file
      Real time = timeBaseLevel*timeMult;
      outputCSV << time;
      for (int i = 0; i != domSums.size(); ++i)
        {
          outputCSV << ", " << domSums[i]*varScale[i];
        }
      outputCSV << "\n";

      // Write the ESD (performs internal check if used or not)
      esd.write();
    }
  outputCSV.close();
}
