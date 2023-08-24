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

#include "parstream.H"
#include "SPMD.H"
#include "DynArray.H"
#include "CH_Cuda_shared_ptr.H"
#include "CH_Cuda_allocator.H"

#ifdef CH_GPU
  #include "testCudaSharedPtr_CUX.H"
  #include "CudaDriver.H"
#endif

#include "UsingBaseNamespace.H"

/// Prototypes:

void
parseTestOptions(int argc, char* argv[]);

int
testCudaSharedPtr1();

/// Global variables for handling output:
static const char *pgmname = "testCudaSharedPtr";
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

#ifdef CH_GPU
  CH_Cuda::g_device.initialize(2*(int)verbose);
#endif

  ///
  // Run the tests
  ///
  int ret = testCudaSharedPtr1();

  if (ret == 0)
    {
      if (verbose)
        {
          pout() << indent << pgmname << " passed all tests" << std::endl;
        }
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

// These macros ease writing of the nested loops.
#define BOUNDSLOOP3(lb, ub, x)                                          \
  for (int x##2 = (lb)[2],                                              \
         x##2_end = (ub)[2] + 1,                                        \
         x##1_end = (ub)[1] + 1,                                        \
         x##0_end = (ub)[0] + 1;                                        \
       x##2 != x##2_end; ++ x##2)                                       \
    for (int x##1 = (lb)[1]; x##1 != x##1_end; ++ x##1)                 \
      for (int x##0 = (lb)[0]; x##0 != x##0_end; ++ x##0)

int testCudaSharedPtr1()
{
  int status = 0;

//--Load the module from the GPU (has functions for several tests)

#ifdef CH_GPU
  /// Handles to module global functions for this application
  enum class App_Global
    {
      arrayTestRW = 0
    };

  /// Handles to module symbols for this application
  enum class App_Symbol
    {
      g_verbose = 0,
      g_stat,
      g_shrDevicePtr,
      g_arrDevicePtr
    };

  // Configure the functions and symbols for this application (test)
  CH_Cuda::g_device.configureChomboLibFunctions(
    CH_Cuda::Application_lib,
    []
    ()
      {
        std::vector<CH_Cuda::DeviceFunctionConfig> config;
        config.reserve(1);
        // All kernels use this block size
        dim3 numThread(4, 4, 4);
        config.emplace_back("arrayTestRW", numThread);
        return config;
      });

  CH_Cuda::g_device.configureChomboLibSymbols(
    CH_Cuda::Application_lib,
    []
    ()
      {
        std::vector<CH_Cuda::DeviceSymbolConfig> config;
        config.reserve(4);
        config.emplace_back("g_verbose");
        config.emplace_back("g_stat");
        config.emplace_back("g_shrDevicePtr");
        config.emplace_back("g_arrDevicePtr");
        return config;
      });

  // Now load the module which automatically loads the functions and symbols
  CH_Cuda::g_device.loadApplicationModule("testCudaSharedPtr",
                                          testCudaSharedPtr);
#endif

//--Set 1: Test construction (array format)

  {
    int lstat = 0;
    CH_Cuda::shared_ptr<dyn::Array<Real, 3>> arrA(
      CH_Cuda::Alloc<dyn::Array<Real, 3>>::newPair(4, 3, 2, 0.5));

    for (int k = 0; k != 2; ++k)
      for (int j = 0; j != 3; ++j)
        for (int i = 0; i != 4; ++i)
          {
            if ((*arrA)(i, j, k) != 0.5) ++lstat;
          }
    {
      auto arrB = arrA;
      if (arrA.use_count() != 2) ++lstat;
      if (arrB.use_count() != 2) ++lstat;
      // Test a reset to nullptr
      arrB.reset();
      if (arrA.use_count() != 1) ++lstat;
      if (arrB.use_count() != 0) ++lstat;
      // Test conversion to boolean
      if (!arrA) ++lstat;
      if (arrB) ++lstat;
    }
    if (arrA.use_count() != 1) ++lstat;

    if (lstat != 0)
      {
        pout() << "Failure in set 1" << std::endl;
        status += lstat;
      }
  }

//--Set 2: Test basic reading/writing of shared pointers on the GPU

#ifdef CH_GPU
  {
    int lstat = 0;

    stc::Vector<int, 3> lb{ -1, 0, 1 };
    stc::Vector<int, 3> ub{  2, 3, 4 };
    CH_Cuda::shared_ptr<dyn::Array<Real, 3>> arrA(
      CH_Cuda::Alloc<dyn::Array<Real, 3>>::newPair(lb, ub));
    BOUNDSLOOP3(lb, ub, i)
      {
        (*arrA)(i0, i1, i2) = i0 + i1 + i2;
      }
    auto arrB(arrA);
    // This copies the dynamic array.  Note that this cannot be part of an
    // asyncrhonous workflow.
    arrA.copyToDevice();

    // Pointer locations
    if (verbose)
      {
        static constexpr const char* ptrAttrLbl[] =
          { "unknown", "host  ", "device" };
        static_assert(CU_MEMORYTYPE_HOST == 1u, "Fix attribute mapping");
        static_assert(CU_MEMORYTYPE_DEVICE == 2u, "Fix attribute mapping");
        unsigned ptrAttr;
        pout() << "Shared pointer locations:\n";
        checkCudaErrors(cuPointerGetAttribute(
                          &ptrAttr,
                          CU_POINTER_ATTRIBUTE_MEMORY_TYPE,
                          reinterpret_cast<CH_Cuda::DevicePointer>(
                            arrA.get())));
        ptrAttr = std::min(std::max(ptrAttr, 0u), 2u);
        pout() << "arrA on " << ptrAttrLbl[ptrAttr] << "         : " <<
          arrA.get() << std::endl;
        checkCudaErrors(cuPointerGetAttribute(
                          &ptrAttr,
                          CU_POINTER_ATTRIBUTE_MEMORY_TYPE,
                          reinterpret_cast<CH_Cuda::DevicePointer>(
                            arrA.devicePtr())));
        ptrAttr = std::min(std::max(ptrAttr, 0u), 2u);
        pout() << "arrA on " << ptrAttrLbl[ptrAttr] << "         : " <<
          arrA.devicePtr() << std::endl;
        pout() << "Array data locations:\n";
        checkCudaErrors(cuPointerGetAttribute(
                          &ptrAttr,
                          CU_POINTER_ATTRIBUTE_MEMORY_TYPE,
                          reinterpret_cast<CH_Cuda::DevicePointer>(
                            arrA->data())));
        ptrAttr = std::min(std::max(ptrAttr, 0u), 2u);
        pout() << "arrA->data() on " << ptrAttrLbl[ptrAttr] << " : " <<
          arrA->data() << std::endl;
        checkCudaErrors(cuPointerGetAttribute(
                          &ptrAttr,
                          CU_POINTER_ATTRIBUTE_MEMORY_TYPE,
                          reinterpret_cast<CH_Cuda::DevicePointer>(
                            arrA->devicePtr())));
        ptrAttr = std::min(std::max(ptrAttr, 0u), 2u);
        pout() << "arrA->data() on " << ptrAttrLbl[ptrAttr] << " : " <<
          arrA->devicePtr() << std::endl;
      }

    // Transfer variables to the device
    CH_Cuda::g_device.copySymbolToDeviceAsync(CH_Cuda::Application_lib,
                                              App_Symbol::g_verbose,
                                              verbose);
    CH_Cuda::g_device.copySymbolToDeviceAsync(CH_Cuda::Application_lib,
                                              App_Symbol::g_stat,
                                              lstat);
    uintptr_t l_shrDevicePtr = reinterpret_cast<uintptr_t>(arrA.devicePtr());
    CH_Cuda::g_device.copySymbolToDeviceAsync(CH_Cuda::Application_lib,
                                              App_Symbol::g_shrDevicePtr,
                                              l_shrDevicePtr);
    if (reinterpret_cast<uintptr_t>(arrB.devicePtr()) != l_shrDevicePtr)
      ++lstat;
    uintptr_t l_arrDevicePtr = reinterpret_cast<uintptr_t>(arrA->devicePtr());
    CH_Cuda::g_device.copySymbolToDeviceAsync(CH_Cuda::Application_lib,
                                              App_Symbol::g_arrDevicePtr,
                                              l_arrDevicePtr);
    if (reinterpret_cast<uintptr_t>(arrB->devicePtr()) != l_arrDevicePtr)
      ++lstat;

    // This copies the data in the dynamic array
    arrA->copyToDeviceAsync();
    // Launch the kernel
    dim3 numBlock;
    CH_Cuda::g_device.launch(CH_Cuda::Application_lib,
                             App_Global::arrayTestRW,
                             numBlock,
                             CH_Cuda::c_defaultStream,
                             arrA, arrB);
    arrA->copyToHostAsync();

    // Transfer status to the host
    CH_Cuda::g_device.copySymbolToHostAsync(CH_Cuda::Application_lib,
                                            App_Symbol::g_stat,
                                            lstat);

    // Barrier
    checkCudaErrors(cuStreamSynchronize(CH_Cuda::c_defaultStream));
    
    BOUNDSLOOP3(lb, ub, i)
      {
        if ((*arrA)(i0, i1, i2) != i0 + i1 + i2 + 0.5) ++lstat;
      }

    if (lstat != 0)
      {
        pout() << "Failure in set 2" << std::endl;
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
