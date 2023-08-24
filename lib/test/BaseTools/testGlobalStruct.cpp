#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//**FIXME testHashTable.cu is built and linked against libbasetools2d.Linux.64.g++.gfortran.DEBUG.OPT.GPU.cu.a before libbasetools2d.Linux.64.g++.gfortran.DEBUG.OPT.GPU.cu.a is constructed (change both testHashTable.cu and src/BaseTools/xxhash.cu

#include <iostream>
#include <cstring>
#include <unordered_map>

#include "parstream.H"
#include "GlobalStruct.H"
#ifdef CH_MPI
#include <mpi.h>
#endif
#ifdef CH_GPU
  #include "testGlobalStruct_CUX.H"
  #include "CudaDriver.H"
#endif

#include "UsingNamespace.H"

/*******************************************************************************
 *
 * Class GlobalStruct: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Static data member initialization
 *============================================================================*/

CH_Hash::StaticTable<int, GlobalStruct>* g_GlobalStruct_table = nullptr;
CH_Hash::DynamicTable<int, GlobalStruct>* g_GlobalStruct_builderTable = nullptr;


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/


void 
GlobalStruct::createInstance()
{
  g_GlobalStruct_builderTable = new CH_Hash::DynamicTable<int, GlobalStruct>;
}

void 
GlobalStruct::destroyInstance()
{
  if (g_GlobalStruct_builderTable)
    {
      delete g_GlobalStruct_builderTable;
    }
  if (g_GlobalStruct_table)
    {
      delete g_GlobalStruct_table;
#ifdef CH_GPU
      // Get the pointer on the device
      CH_Cuda::DevicePointer d_table;
      CH_Cuda::g_device.copySymbolToHostAsync(CH_Cuda::Application_lib,
                                              App_Symbol::g_table,
                                              d_table);
      checkCudaErrors(cuMemFree(d_table));
#endif
    }
}


/*******************************************************************************
 *
 * Standard test stuff
 *
 ******************************************************************************/

/// Prototypes:

void
parseTestOptions(int argc, char* argv[]) ;

int
testGlobalStruct1();

/// Global variables for handling output:
static const char *pgmname = "testGlobalStruct" ;
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
    pout() << indent2 << "Beginning " << pgmname << " ..." << std::endl ;

#ifdef CH_GPU
  CH_Cuda::g_device.initialize(2*(int)verbose);
#endif

  ///
  // Run the tests
  ///
  int ret = testGlobalStruct1();

  if (ret == 0)
    {
      pout() << indent << pgmname << " passed all tests" << std::endl ;
    }
  else
    {
      pout() << indent << pgmname << " failed " << ret << " test(s)"
             << std::endl ;
    }

#ifdef CH_GPU
  CH_Cuda::g_device.finalize();
#endif
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return ret;
}


/*******************************************************************************
 *
 * Routine testHashTable
 *
 ******************************************************************************/

int
testGlobalStruct1()
{
  int status = 0;

//--Load the module from the GPU (has functions for several tests)

#ifdef CH_GPU

  // Configure the functions and symbols for this application (test)
  CH_Cuda::g_device.configureChomboLibFunctions(
    CH_Cuda::Application_lib,
    []
    ()
      {
        std::vector<CH_Cuda::DeviceFunctionConfig> config;
        config.reserve(1);
        // App_Global 0
        config.emplace_back("testGlobalStructGet", 4);  // 4 threads
        return config;
      });

  CH_Cuda::g_device.configureChomboLibSymbols(
    CH_Cuda::Application_lib,
    []
    ()
      {
        std::vector<CH_Cuda::DeviceSymbolConfig> config;
        config.reserve(2);
        // App_Symbol 0
        config.emplace_back("g_verbose");
        // App_Symbol 1
        config.emplace_back("g_stat");
        // App_Symbol 2
        config.emplace_back("g_GlobalStruct_table");
        return config;
      });

  // Now load the module which automatically loads the functions and symbols
  CH_Cuda::g_device.loadApplicationModule("testGlobalStruct",
                                          testGlobalStruct);
#endif

//--Set 1 (test a HashGraphTable)

  {
    int lstat = 0;
    GlobalStruct::createInstance();
    {
      GlobalStruct& gs = GlobalStruct::get(1);
      gs.a_val = 1.1;
    }
    {
      GlobalStruct& gs = GlobalStruct::get(2);
      gs.a_val = 2.2;
    }
    GlobalStruct::finalize();
    {
      GlobalStruct& gs = GlobalStruct::get(1);
      if (gs.a_val != 1.1) ++lstat;
      if (verbose)
        {
          pout() << gs.a_val << std::endl;
        }
    }
    {
      GlobalStruct& gs = GlobalStruct::get(2);
      if (gs.a_val != 2.2) ++lstat;
      if (verbose)
        {
          pout() << gs.a_val << std::endl;
        }
    }
    if (lstat != 0)
      {
        pout() << "Failure in set 1" << std::endl;
        status += lstat;
      }

#ifdef CH_GPU
    lstat = 0;
    // Transfer to the device
    CH_Cuda::g_device.copySymbolToDeviceAsync(CH_Cuda::Application_lib,
                                              App_Symbol::g_verbose,
                                              verbose);
    CH_Cuda::g_device.copySymbolToDeviceAsync(CH_Cuda::Application_lib,
                                              App_Symbol::g_stat,
                                              lstat);

    // Launch the first kernel
    dim3 numBlock;
    dim3 numThread;
    numThread.x = 24;
    CH_Cuda::g_device.launch(CH_Cuda::Application_lib,
                             App_Global::testGlobalStructGet,
                             numBlock, numThread,
                             CH_Cuda::c_defaultStream);

    // Transfer to the host
    CH_Cuda::g_device.copySymbolToHostAsync(CH_Cuda::Application_lib,
                                            App_Symbol::g_stat,
                                            lstat);

    // Barrier
    checkCudaErrors(cuStreamSynchronize(CH_Cuda::c_defaultStream));

    if (lstat != 0)
      {
        pout() << "Failure in set 1 on GPU" << std::endl;
        status += lstat;
      }
#endif
    GlobalStruct::destroyInstance();
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
