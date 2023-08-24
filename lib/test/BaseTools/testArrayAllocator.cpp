#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// Define to see constructor calls.
// #define DEBUGARRAYALLOC

#include <iostream>
#include <limits>
#include <functional>
#include <cstring>
using std::endl;

#include "REAL.H"
#include "parstream.H"
#include "SPMD.H"
#include "ArrayAllocator.H"
#include "Misc.H"

#ifdef CH_GPU
#include "testArrayAllocator_CUX.H"
#include "CudaDriver.H"
#endif

#include "UsingBaseNamespace.H"

/// Prototypes:

Real BaseFabRealSetVal = BASEFAB_REAL_SETVAL;

void
parseTestOptions( int argc ,char* argv[] ) ;

int
testArrayAllocator1();

/// Global variables for handling output:
static const char *pgmname = "testArrayAllocator" ;
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
  int ret = testArrayAllocator1() ;

  if (ret == 0)
    {
      if (verbose)
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

template <typename T>
using BaseFabArrayImpl =
  Array_impl<T, DefaultArrayAlloc<T, ArrayClassIndex::BaseFab> >;

#ifdef CH_GPU
template <typename T>
using BaseFabArrayImplGPU =
  Array_impl<T, CUDAArrayAlloc<T, ArrayClassIndex::BaseFab> >;
#endif

// Dummy box class
class Box
{
public:
  Box(const int a_lo, const int a_hi)
    :
    m_lo(a_lo),
    m_hi(a_hi)
    { }
  int m_lo;
  int m_hi;
};

int testArrayAllocator1()
{
  int status = 0;

//--Load the module from the GPU (has functions for several tests)

#ifdef CH_GPU
  /// Handles to module global functions for this application
  enum class App_Global
    {
      arrayTestRW       = 0,
      arrayTestAliasRaw,
      arrayTestAlloc,
    };

  /// Handles to module symbols for this application
  enum class App_Symbol
    {
      g_verbose = 0,
      g_stat,
      g_devicePtr
    };

  // Configure the functions and symbols for this application (test)
  CH_Cuda::g_device.configureChomboLibFunctions(
    CH_Cuda::Application_lib,
    []
    ()
      {
        std::vector<CH_Cuda::DeviceFunctionConfig> config;
        // App_Global 0
        config.emplace_back("arrayTestRW", 4);  // 4 threads
        config.emplace_back("arrayTestAliasRaw", 16);  // 16 threads
        config.emplace_back("arrayTestAlloc", 16);  // 16 threads
        return config;
      });

  CH_Cuda::g_device.configureChomboLibSymbols(
    CH_Cuda::Application_lib,
    []
    ()
      {
        std::vector<CH_Cuda::DeviceSymbolConfig> config;
        // App_Symbol 0
        config.emplace_back("g_verbose");
        config.emplace_back("g_stat");
        config.emplace_back("g_devicePtr");
        return config;
      });

  // Now load the module which automatically loads the functions and symbols
  CH_Cuda::g_device.loadApplicationModule("testArrayAllocator",
                                          testArrayAllocator);
#endif

  if (verbose) pout() << "=> Testing DefaultArrayAlloc\n";

//--Set 1: Test construction/copy/move of default

  {
    if (verbose) pout() << "-> Construction:\n";
    DefaultArrayAlloc<Real, ArrayClassIndex::BaseFab> Aa;
    if (verbose) pout() << "-> Copy:\n";
    DefaultArrayAlloc<Real, ArrayClassIndex::BaseFab> Ab(Aa);
    if (verbose) pout() << "-> Move:\n";
    DefaultArrayAlloc<Real, ArrayClassIndex::BaseFab> Ac = std::move(Aa);
  }

//--Set 2: Test bits and size in size variable

  {
    int lstat = 0;
    BaseFabArrayImpl<Real> arrA;
    // Use hard numbers on RHS.  The tests should be reexamined if another
    // device is added.
    arrA.size(20, (unsigned)AllocOn::cpu, AllocBy::array);
    if (arrA.size() != 20) ++lstat;
    if (arrA.allocOn() != 4) ++lstat;
    if (static_cast<unsigned>(arrA.allocBy()) != 1) ++lstat;
    if (!arrA.allocOn_cpu()) ++lstat;
    if (arrA.allocOn_gpu()) ++lstat;

    arrA.size(0, (unsigned)AllocOn::gpu, AllocBy::alias);
    if (arrA.size() != 0) ++lstat;
    if (arrA.allocOn() != 8) ++lstat;
    if (static_cast<unsigned>(arrA.allocBy()) != 2) ++lstat;
    if (arrA.allocOn_cpu()) ++lstat;
    if (!arrA.allocOn_gpu()) ++lstat;

    arrA.size(1, (unsigned)(AllocOn::all), AllocBy::raw);
    if (arrA.size() != 1) ++lstat;
    if (arrA.allocOn() != 12) ++lstat;
    if (static_cast<unsigned>(arrA.allocBy()) != 3) ++lstat;
    if (!arrA.allocOn_cpu()) ++lstat;
    if (!arrA.allocOn_gpu()) ++lstat;

    arrA.size(BaseFabArrayImpl<Real>::c_max_size,
              (unsigned)(AllocOn::all),
              AllocBy::raw);
    if (arrA.size() != BaseFabArrayImpl<Real>::c_max_size) ++lstat;
    if (arrA.allocOn() != 12) ++lstat;
    if (static_cast<unsigned>(arrA.allocBy()) != 3) ++lstat;

    arrA.size(64, AllocBy::array);
    if (arrA.size() != 64) ++lstat;
    if (arrA.allocOn() != 12) ++lstat;
    if (static_cast<unsigned>(arrA.allocBy()) != 1) ++lstat;

    arrA.size(128);
    if (arrA.size() != 128) ++lstat;
    if (arrA.allocOn() != 12) ++lstat;
    if (static_cast<unsigned>(arrA.allocBy()) != 1) ++lstat;

    arrA.setAllocBy(AllocBy::raw);
    if (arrA.size() != 128) ++lstat;
    if (arrA.allocOn() != 12) ++lstat;
    if (static_cast<unsigned>(arrA.allocBy()) != 3) ++lstat;

    arrA.setAllocBy(AllocBy::alias);
    if (arrA.size() != 128) ++lstat;
    if (arrA.allocOn() != 12) ++lstat;
    if (static_cast<unsigned>(arrA.allocBy()) != 2) ++lstat;

    arrA.size(0, (unsigned)AllocOn::none, AllocBy::none);
    if (arrA.size() != 0) ++lstat;
    if (arrA.allocOn() != 0) ++lstat;
    if (static_cast<unsigned>(arrA.allocBy()) != 0) ++lstat;

    if (lstat != 0)
      {
        pout() << "Failure in set 2" << std::endl;
        status += lstat;
      }
  }

//--Set 3: Test allocation and initialization

  constexpr int prec = std::numeric_limits<Real>::digits10 - 2;

  {
    int lstat = 0;
    BaseFabArrayImpl<Real> arrayA;
    arrayA.define(4);
    if (verbose)
      {
        pout() << "Default setReal: " << BaseFabRealSetVal << std::endl;
      }
    if (arrayA.size() != 4) ++lstat;
    for (int i = 0; i != 4; ++i)
      {
        if (Misc::compare(arrayA[i], BaseFabRealSetVal, prec)) ++lstat;
      }
    arrayA.define(3, 1.5);
    if (arrayA.size() != 3) ++lstat;
    for (int i = 0; i != 3; ++i)
      {
        if (Misc::compare(arrayA[i], 1.5, prec)) ++lstat;
      }
    // Alias to arrayA
    BaseFabArrayImpl<Real> arrayA1;
    arrayA1.defineAlias(arrayA.dataPtr(), 2);
    if (arrayA1.size() != 2) ++lstat;
    for (int i = 0; i != 2; ++i)
      {
        if (Misc::compare(arrayA1[i], 1.5, prec)) ++lstat;
      }
    // Using defineRaw should reinitialize memory
    BaseFabArrayImpl<Real> arrayA2;
    arrayA2.defineRaw(arrayA.dataPtr(), 2);
    if (arrayA2.size() != 2) ++lstat;
    for (int i = 0; i != 2; ++i)
      {
        if (Misc::compare(arrayA2[i], BaseFabRealSetVal, prec)) ++lstat;
      }
    // Usine defneRaw should reinitialize memory
    BaseFabArrayImpl<Real> arrayA3;
    arrayA3.defineRaw(arrayA.dataPtr(), 2, -1.5);
    if (arrayA3.size() != 2) ++lstat;
    for (int i = 0; i != 2; ++i)
      {
        if (Misc::compare(arrayA[i], -1.5, prec)) ++lstat;
        if (Misc::compare(arrayA1[i], -1.5, prec)) ++lstat;
        if (Misc::compare(arrayA3[i], -1.5, prec)) ++lstat;
      }
    // Move of A should not affect aliases
    BaseFabArrayImpl<Real> arrayB = std::move(arrayA);
    if (arrayB.size() != 3) ++lstat;
    if (arrayA.dataPtr() != nullptr) ++lstat;
    if (arrayA.allocBy() != AllocBy::none) ++lstat;
    for (int i = 0; i != 2; ++i)
      {
        if (Misc::compare(arrayA1[i], -1.5, prec)) ++lstat;
        if (Misc::compare(arrayA2[i], -1.5, prec)) ++lstat;
        if (Misc::compare(arrayA3[i], -1.5, prec)) ++lstat;
        if (Misc::compare(arrayB[i], -1.5, prec)) ++lstat;
      }
    {
      if (Misc::compare(arrayA1[2], 1.5, prec)) ++lstat;
      if (Misc::compare(arrayA2[2], 1.5, prec)) ++lstat;
      if (Misc::compare(arrayA3[2], 1.5, prec)) ++lstat;
      if (Misc::compare(arrayB[2], 1.5, prec)) ++lstat;
    }
    // Test construction of class type
    Box* dataPtrA;
    BaseFabArrayImpl<Box> arrayBoxB;
    {
      BaseFabArrayImpl<Box> arrayBoxA;
      arrayBoxA.define(4, -2, -1);
      for (int i = 0; i != 4; ++i)
        {
          if (arrayBoxA[i].m_lo != -2) ++lstat;
          if (arrayBoxA[i].m_hi != -1) ++lstat;
        }
      dataPtrA = arrayBoxA.dataPtr();
      if (verbose) pout() << "-> Assignment move:\n";
      arrayBoxB = std::move(arrayBoxA);
    }  // arrayBoxA is now destroyed but arrayBoxB has assumed resources
    Box* dataPtrB = arrayBoxB.dataPtr();
    if (dataPtrA != dataPtrB) ++lstat;
    for (int i = 0; i != 4; ++i)
      {
        if (arrayBoxB[i].m_lo != -2) ++lstat;
        if (arrayBoxB[i].m_hi != -1) ++lstat;
      }

    if (lstat != 0)
      {
        pout() << "Failure in set 3" << std::endl;
        status += lstat;
      }
  }

//--Set 4: Test basic allocations on the GPU

#ifdef CH_GPU
  if (verbose) pout() << "=> Testing CudaArrayAlloc\n";
  {
    int lstat = 0;
    Real* hostPtrA = nullptr;
    Real* devicePtrA = nullptr;
    BaseFabArrayImplGPU<Real> arrayB;
    {
      BaseFabArrayImplGPU<Real> arrayA(AllocOn::all);
      arrayA.define(4);
      hostPtrA = arrayA.dataPtr();
      devicePtrA = arrayA.devicePtr();
      if (verbose) pout() << "-> Assignment move of CUDA alloc:\n";
      arrayB = std::move(arrayA);
    }  // arrayA is now destroyed but arrayB has assumed resources
    if (hostPtrA != arrayB.dataPtr()) ++lstat;
    if (devicePtrA != arrayB.devicePtr()) ++lstat;
    // And of course we should have proper initialization still
    for (int i = 0; i != 4; ++i)
      {
        if (Misc::compare(arrayB[i], BaseFabRealSetVal, prec)) ++lstat;
      }
    if (verbose)
      {
        ReportAllocatedMemory(pout());
      }

    // Test allocation only on CPU
    BaseFabArrayImplGPU<Real> arrayC(AllocOn::cpu);
    arrayC.define(4);
    if (verbose)
      {
        pout() << "-> Add allocation only on CPU\n";
        ReportAllocatedMemory(pout());
      }

    // Test allocation only on GPU
    BaseFabArrayImplGPU<Real> arrayD(AllocOn::gpu);
    arrayD.define(4);
    if (verbose)
      {
        pout() << "-> Add allocation only on GPU\n";
        ReportAllocatedMemory(pout());
      }

    if (lstat != 0)
      {
        pout() << "Failure in set 4" << std::endl;
        status += lstat;
      }
  }
#endif

//--Set 5: Test basic reading/writing on the GPU

#ifdef CH_GPU
  {
    int lstat = 0;

    BaseFabArrayImplGPU<Real> arrA;
    arrA.define(4);
    int c = 0;
    for (auto& x : arrA)
      {
        x = c++ + 0.5;
      }

    // Pointer locations
    if (verbose)
      {
        pout() << "Array pointer locations:\n";
        static constexpr const char* ptrAttrLbl[] =
          { "unknown", "host  ", "device" };
        static_assert(CU_MEMORYTYPE_HOST == 1u, "Fix attribute mapping");
        static_assert(CU_MEMORYTYPE_DEVICE == 2u, "Fix attribute mapping");
        unsigned ptrAttr;
        checkCudaErrors(cuPointerGetAttribute(
                          &ptrAttr,
                          CU_POINTER_ATTRIBUTE_MEMORY_TYPE,
                          reinterpret_cast<CH_Cuda::DevicePointer>(
                            arrA.data())));
        ptrAttr = std::min(std::max(ptrAttr, 0u), 2u);
        pout() << "arrA on " << ptrAttrLbl[ptrAttr] << "         : " <<
          arrA.data() << std::endl;
        checkCudaErrors(cuPointerGetAttribute(
                          &ptrAttr,
                          CU_POINTER_ATTRIBUTE_MEMORY_TYPE,
                          reinterpret_cast<CH_Cuda::DevicePointer>(
                            arrA.devicePtr())));
        ptrAttr = std::min(std::max(ptrAttr, 0u), 2u);
        pout() << "arrA on " << ptrAttrLbl[ptrAttr] << "         : " <<
          arrA.devicePtr() << std::endl;
      }

    // Transfer variables to the device
    CH_Cuda::g_device.copySymbolToDeviceAsync(CH_Cuda::Application_lib,
                                              App_Symbol::g_verbose,
                                              verbose);
    CH_Cuda::g_device.copySymbolToDeviceAsync(CH_Cuda::Application_lib,
                                              App_Symbol::g_stat,
                                              lstat);
    uintptr_t l_devicePtr = reinterpret_cast<uintptr_t>(arrA.devicePtr());
    CH_Cuda::g_device.copySymbolToDeviceAsync(CH_Cuda::Application_lib,
                                              App_Symbol::g_devicePtr,
                                              l_devicePtr);

    arrA.copyToDeviceAsync();
    // Launch the kernel
    dim3 numBlock;
    dim3 numThread(4);
    CH_Cuda::g_device.launch(CH_Cuda::Application_lib,
                             App_Global::arrayTestRW,
                             numBlock, numThread,
                             CH_Cuda::c_defaultStream,
                             arrA);
    arrA.copyToHostAsync();

    // Transfer status to the host
    CH_Cuda::g_device.copySymbolToHostAsync(CH_Cuda::Application_lib,
                                            App_Symbol::g_stat,
                                            lstat);

    // Barrier
    checkCudaErrors(cuStreamSynchronize(CH_Cuda::c_defaultStream));

    c = 0;
    for (const auto x : arrA)
      {
        if (x != (c++ + 0.75)) ++lstat;
      }

    if (lstat != 0)
      {
        pout() << "Failure in set 5" << std::endl;
        status += lstat;
      }
  }
#endif

//--Set 6: Test raw and aliased arrays on the GPU

#ifdef CH_GPU
  {
    int lstat = 0;
    // Transfer status to the device
    CH_Cuda::g_device.copySymbolToDeviceAsync(CH_Cuda::Application_lib,
                                              App_Symbol::g_stat,
                                              lstat);

    // Launch the kernel
    dim3 numBlock;
    dim3 numThread(16);
    CH_Cuda::g_device.launch(CH_Cuda::Application_lib,
                             App_Global::arrayTestAliasRaw,
                             numBlock, numThread,
                             CH_Cuda::c_defaultStream);

    // Transfer status to the host
    CH_Cuda::g_device.copySymbolToHostAsync(CH_Cuda::Application_lib,
                                            App_Symbol::g_stat,
                                            lstat);

    // Barrier
    checkCudaErrors(cuStreamSynchronize(CH_Cuda::c_defaultStream));

    if (lstat != 0)
      {
        pout() << "Failure in set 6" << std::endl;
        status += lstat;
      }
  }
#endif

//--Set 7: Test device allocations

#ifdef CH_GPU
  {
    int lstat = 0;
    // Transfer status to the device
    CH_Cuda::g_device.copySymbolToDeviceAsync(CH_Cuda::Application_lib,
                                              App_Symbol::g_stat,
                                              lstat);

    // Launch the kernel
    dim3 numBlock;
    dim3 numThread(16);
    CH_Cuda::g_device.launch(CH_Cuda::Application_lib,
                             App_Global::arrayTestAlloc,
                             numBlock, numThread,
                             CH_Cuda::c_defaultStream);

    // Transfer status to the host
    CH_Cuda::g_device.copySymbolToHostAsync(CH_Cuda::Application_lib,
                                            App_Symbol::g_stat,
                                            lstat);

    // Barrier
    checkCudaErrors(cuStreamSynchronize(CH_Cuda::c_defaultStream));

    if (lstat != 0)
      {
        pout() << "Failure in set 7" << std::endl;
        status += lstat;
      }
  }
#endif

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
