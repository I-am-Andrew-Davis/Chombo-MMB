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
using std::endl;
#include "parstream.H"
#include "CH_HashTable.H"
#ifdef CH_MPI
#include <mpi.h>
#endif
#ifdef CH_GPU
  #include "testHashTable_CUX.H"
  #include "CudaDriver.H"
#endif

#include "UsingNamespace.H"

/// Prototypes:

void
parseTestOptions(int argc, char* argv[]) ;

int
testHashTable1();

/// Global variables for handling output:
static const char *pgmname = "testHashTable" ;
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

#ifdef CH_GPU
  CH_Cuda::g_device.initialize(2*(int)verbose);
#endif

  ///
  // Run the tests
  ///
  int ret = testHashTable1();

  if (ret == 0)
    {
      pout() << indent << pgmname << " passed all tests" << endl ;
    }
  else
    {
      pout() << indent << pgmname << " failed " << ret << " test(s)" << endl ;
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
testHashTable1()
{
  int status = 0;

//--Load the module from the GPU (has functions for several tests)

#ifdef CH_GPU
  /// Handles to module global functions for this application
  enum class App_Global
    {
      hashGraphFind = 0
    };

  /// Handles to module symbols for this application
  enum class App_Symbol
    {
      g_verbose   = 0,
      g_stat      = 1,
      g_devicePtr = 2
    };

  // Configure the functions and symbols for this application (test)
  CH_Cuda::g_device.configureChomboLibFunctions(
    CH_Cuda::Application_lib,
    []
    ()
      {
        std::vector<CH_Cuda::DeviceFunctionConfig> config;
        config.reserve(1);
        // App_Global 0
        config.emplace_back("hashGraphFind", 4);  // 4 threads
        return config;
      });

  CH_Cuda::g_device.configureChomboLibSymbols(
    CH_Cuda::Application_lib,
    []
    ()
      {
        std::vector<CH_Cuda::DeviceSymbolConfig> config;
        config.reserve(3);
        // App_Symbol 0
        config.emplace_back("g_verbose");
        config.emplace_back("g_stat");
        config.emplace_back("g_devicePtr");
        return config;
      });

  // Now load the module which automatically loads the functions and symbols
  CH_Cuda::g_device.loadApplicationModule("testHashTable", testHashTable);
#endif

//--Set 1 (test a HashGraphTable)

  {
    int stat1 = 0;
    using HashTable0 = std::unordered_map<int, double, CH_Hash::XXHash32<int>>;
    HashTable0 table0;
    table0.insert({  0, 0.0  });
    table0.insert({  9, 0.9  });
    table0.insert({ 19, 0.19 });
    table0.insert({  8, 0.8  });
    table0.insert({ 18, 0.18 });
    table0.insert({  1, 0.1  });
    using HashTable1 = CH_Hash::HashGraphTable<int, double>;
    HashTable1 table1(table0, -1);
    if (table1.find( 0)->second !=  0.0) ++stat1;
    if (table1.find( 1)->second !=  0.1) ++stat1;
    if (table1.find( 8)->second !=  0.8) ++stat1;
    if (table1.find( 9)->second !=  0.9) ++stat1;
    if (table1.find(18)->second != 0.18) ++stat1;
    if (table1.find(19)->second != 0.19) ++stat1;
    if (stat1 != 0)
      {
        pout() << "Failure in set 1.1" << endl;
        status += stat1;
      }

    // Test an alias
    {
      int stat2 = 0;
      HashTable1 table1Alias(CH_Hash::alias_format{}, table1);
      if (table1Alias.find( 0)->second !=  0.0) ++stat2;
      if (table1Alias.find( 1)->second !=  0.1) ++stat2;
      if (table1Alias.find( 8)->second !=  0.8) ++stat2;
      if (table1Alias.find( 9)->second !=  0.9) ++stat2;
      if (table1Alias.find(18)->second != 0.18) ++stat2;
      if (table1Alias.find(19)->second != 0.19) ++stat2;
      if (stat2 != 0)
        {
          pout() << "Failure in set 1.2" << endl;
          status += stat2;
        }
    }
    int stat3 = 0;
    if (table1.find( 0)->second !=  0.0) ++stat3;
    if (table1.find( 1)->second !=  0.1) ++stat3;
    if (table1.find( 8)->second !=  0.8) ++stat3;
    if (table1.find( 9)->second !=  0.9) ++stat3;
    if (table1.find(18)->second != 0.18) ++stat3;
    if (table1.find(19)->second != 0.19) ++stat3;
    if (stat3 != 0)
      {
        pout() << "Failure in set 1.3" << endl;
        status += stat3;
      }

#ifdef CH_GPU
    stat1 = 0;
    if (verbose)
      {
        pout() << "Table pointer locations:\n";
        static constexpr const char* ptrAttrLbl[] =
          { "unknown", "host  ", "device" };
        static_assert(CU_MEMORYTYPE_HOST == 1u, "Fix attribute mapping");
        static_assert(CU_MEMORYTYPE_DEVICE == 2u, "Fix attribute mapping");
        unsigned ptrAttr;
        checkCudaErrors(cuPointerGetAttribute(
                          &ptrAttr,
                          CU_POINTER_ATTRIBUTE_MEMORY_TYPE,
                          reinterpret_cast<CH_Cuda::DevicePointer>(
                            table1.m_elem.data())));
        ptrAttr = std::min(std::max(ptrAttr, 0u), 2u);
        pout() << "m_elem on " << ptrAttrLbl[ptrAttr] << "  : " <<
          table1.m_elem.data() << std::endl;
        checkCudaErrors(cuPointerGetAttribute(
                          &ptrAttr,
                          CU_POINTER_ATTRIBUTE_MEMORY_TYPE,
                          reinterpret_cast<CH_Cuda::DevicePointer>(
                            table1.m_elem.devicePtr())));
        ptrAttr = std::min(std::max(ptrAttr, 0u), 2u);
        pout() << "m_elem on " << ptrAttrLbl[ptrAttr] << "  : " <<
          table1.m_elem.devicePtr() << std::endl;
        checkCudaErrors(cuPointerGetAttribute(
                          &ptrAttr,
                          CU_POINTER_ATTRIBUTE_MEMORY_TYPE,
                          reinterpret_cast<CH_Cuda::DevicePointer>(
                            table1.m_offset.data())));
        ptrAttr = std::min(std::max(ptrAttr, 0u), 2u);
        pout() << "m_offset on " << ptrAttrLbl[ptrAttr] << ": " <<
          table1.m_offset.data() << std::endl;
        checkCudaErrors(cuPointerGetAttribute(
                          &ptrAttr,
                          CU_POINTER_ATTRIBUTE_MEMORY_TYPE,
                          reinterpret_cast<CH_Cuda::DevicePointer>(
                            table1.m_offset.devicePtr())));
        ptrAttr = std::min(std::max(ptrAttr, 0u), 2u);
        pout() << "m_offset on " << ptrAttrLbl[ptrAttr] << ": " <<
          table1.m_offset.devicePtr() << std::endl;
      }

    CH_Cuda::g_device.copySymbolToDeviceAsync(CH_Cuda::Application_lib,
                                              App_Symbol::g_verbose,
                                              verbose);
    CH_Cuda::g_device.copySymbolToDeviceAsync(CH_Cuda::Application_lib,
                                              App_Symbol::g_stat,
                                              stat1);
    uintptr_t l_devicePtr[2] = 
    {
      reinterpret_cast<uintptr_t>(table1.m_elem.devicePtr()),
      reinterpret_cast<uintptr_t>(table1.m_offset.devicePtr())
    };
    CH_Cuda::g_device.copySymbolToDeviceAsync(CH_Cuda::Application_lib,
                                              App_Symbol::g_devicePtr,
                                              l_devicePtr);
      
    table1.copyToDeviceAsync();

    // Launch the first kernel
    dim3 numBlock;
    dim3 numThread;
    numThread.x = 24;
    CH_Cuda::g_device.launch(CH_Cuda::Application_lib,
                             App_Global::hashGraphFind,
                             numBlock, numThread,
                             CH_Cuda::c_defaultStream,
                             table1);

    // Transfer to the host
    CH_Cuda::g_device.copySymbolToHostAsync(CH_Cuda::Application_lib,
                                            App_Symbol::g_stat,
                                            stat1);

    // Barrier
    checkCudaErrors(cuStreamSynchronize(CH_Cuda::c_defaultStream));

    if (stat1 != 0)
      {
        pout() << "Failure in set 1 on GPU" << endl;
        status += stat1;
      }
#endif
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
