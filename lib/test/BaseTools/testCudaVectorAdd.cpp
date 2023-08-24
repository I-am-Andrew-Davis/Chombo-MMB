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
#include "SPMD.H"
#ifdef CH_GPU
  #include "testCudaVectorAdd_CUX.H"
  #include "CudaDriver.H"
#endif

#include "UsingBaseNamespace.H"

/// Prototypes:

void
parseTestOptions(int argc, char* argv[]);

int
testCUDAVectorAdd_libLevel();

int
testCUDAVectorAdd_lowLevel();

/// Global variables for handling output:
static const char *pgmname = "testCUDAVectorAdd" ;
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
    pout() << indent2 << "Beginning " << pgmname << " ..." << std::endl;

#ifdef CH_GPU
  CH_Cuda::g_device.initialize(2*(int)verbose);
#endif

  ///
  // Run the tests
  ///
  int ret = testCUDAVectorAdd_libLevel();
  ret +=    testCUDAVectorAdd_lowLevel();

  if (ret == 0)
    {
      pout() << indent << pgmname << " passed all tests" << std::endl;
    }
  else
    {
      pout() << indent << pgmname << " failed " << ret << " test(s)"
             << std::endl;
    }

#ifdef CH_GPU
  CH_Cuda::g_device.finalize();
#endif
#ifdef CH_MPI
  MPI_Finalize();
#endif
  return ret;
}

//--This is high level approach does not has built int support for using
//--deivce global functions and sybmols in Chombo libraries and in the
//--application (test in this case).  Recommended for most cases.

int testCUDAVectorAdd_libLevel()
{
  int status = 0;

#ifdef CH_GPU

  const int n = 64;
  const int elemPerBlock = 16;
  dim3 numBlock;
  dim3 numThread;
  numThread.x = elemPerBlock;
  numBlock.x = n/numThread.x;

  /// Handles to module global functions for this application
  enum class App_Global
    {
      vvAddTo = 0,
      vsAddTo = 1
    };

  // If there are functions or symbols from standard Chombo libraries, you
  // would mark them here using, e.g., BaseTools_lib
  // Here we mark the functions and symbols for this application (test).
  std::string namevv = "vvAddTo";  // Just for testing
  CH_Cuda::g_device.configureChomboLibFunctions(
    CH_Cuda::Application_lib,
    [&]
    ()
      {
        std::vector<CH_Cuda::DeviceFunctionConfig> config;
        config.reserve(2);
        // App_Global 0
        config.emplace_back(namevv, numThread);
        // App_Global 0
        config.emplace_back("vsAddTo", numThread);
        return config;
      });

  // Now load the module which automatically loads the functions and symbols

  CH_Cuda::g_device.loadApplicationModule("testCudaVectorAdd",
                                          testCudaVectorAdd);

  // Allocate host side memory
  CH_Cuda::SymbolPair<float> x;
  CH_Cuda::SymbolPair<float> y;
  checkCudaErrors(cuMemAllocHost(reinterpret_cast<void**>(&x.host),
                                 n*sizeof(float)));
  checkCudaErrors(cuMemAllocHost(reinterpret_cast<void**>(&y.host),
                                 n*sizeof(float)));

  for (int i = 0; i != n; ++i)
    {
      x.host[i] = i;
      y.host[i] = i;
    }

  // Allocate on device
  checkCudaErrors(cuMemAlloc(&x.device, n*sizeof(float)));
  checkCudaErrors(cuMemAlloc(&y.device, n*sizeof(float)));

  // Create the stream
  CUstream stream;
  checkCudaErrors(cuStreamCreate(&stream, 0));

  // Transfer to the device
  x.copyHtoDAsync(n, stream);
  y.copyHtoDAsync(n, stream);

  // Launch the first kernel
  CH_Cuda::g_device.launch(CH_Cuda::Application_lib,
                           App_Global::vvAddTo,
                           numBlock,
                           // Note: You don't have to provide numThread if the
                           // value configured for the function is desired
                           numThread,
                           stream,
                           x.device, y.device);

  // Launch the second kernel
  CH_Cuda::g_device.launch(CH_Cuda::Application_lib,
                           App_Global::vsAddTo,
                           numBlock,
                           // In this example, numThread from function
                           // configuration is used.
                           stream,
                           x.device, 2.0f);

  // Get the results
  x.copyDtoHAsync(n, stream);

  // Barrier
  checkCudaErrors(cuStreamSynchronize(stream));

  // Destroy the stream
  checkCudaErrors(cuStreamDestroy(stream));

  // Destroy resources on the device
  checkCudaErrors(cuMemFree(x.device));
  checkCudaErrors(cuMemFree(y.device));

  // Check the results on the host
  for (int i = 0; i != n; ++i)
    {
      if (x.host[i] != 2*y.host[i] + 2) ++status;
    }

  // Destroy host-side memory
  checkCudaErrors(cuMemFreeHost(x.host));
  checkCudaErrors(cuMemFreeHost(y.host));

  // Unload the modules.  Normally you don't have to do this but we may want to
  // call testCUDAVectorAdd_lowLevel() which loads the same module in a
  // different way.
  CH_Cuda::g_device.unloadModule0();
  
#endif  /* CH_GPU */

  return status;
}

//--This is low level approach does not have support for Chombo module functions
//--and symbols.  Recommended only for external modules.  Functions and symbols
//--must be manually loaded and tracked.

int testCUDAVectorAdd_lowLevel()
{
  int status = 0;

#ifdef CH_GPU

  const int n = 64;
  const int elemPerBlock = 16;
  dim3 numBlock;
  dim3 numThread;
  numThread.x = elemPerBlock;
  numBlock.x = n/numThread.x;

  class VecAddToModule : public CH_Cuda::DriverModule
  {
  public:
    CH_Cuda::DeviceFunction vvAddTo;
    CH_Cuda::DeviceFunction vsAddTo;
    CH_Cuda::SymbolPair<float> x;
    CH_Cuda::SymbolPair<float> y;
  };

  // Initialize module
  VecAddToModule *vecAddToMdl =
    CH_Cuda::g_device.loadModule<VecAddToModule>("testCudaVectorAdd",
                                                 testCudaVectorAdd);

  // Retrieve functions
  checkCudaErrors(cuModuleGetFunction(&vecAddToMdl->vvAddTo,
                                      vecAddToMdl->m_module,
                                      "vvAddTo"));
  checkCudaErrors(cuModuleGetFunction(&(vecAddToMdl->vsAddTo),
                                      vecAddToMdl->m_module,
                                      "vsAddTo"));

  // Allocate host side memory
  float *x;
  float *y;
  checkCudaErrors(cuMemAllocHost(reinterpret_cast<void**>(&x),
                                 n*sizeof(float)));
  checkCudaErrors(cuMemAllocHost(reinterpret_cast<void**>(&y),
                                 n*sizeof(float)));
  vecAddToMdl->x.host = x;
  vecAddToMdl->y.host = y;

  for (int i = 0; i != n; ++i)
    {
      x[i] = i;
      y[i] = i;
    }

  // Allocate on device
  checkCudaErrors(cuMemAlloc(&vecAddToMdl->x.device, n*sizeof(float)));
  checkCudaErrors(cuMemAlloc(&vecAddToMdl->y.device, n*sizeof(float)));

  // Create the stream
  CUstream stream;
  checkCudaErrors(cuStreamCreate(&stream, 0));

  // Transfer to the device
  vecAddToMdl->x.copyHtoDAsync(n, stream);
  vecAddToMdl->y.copyHtoDAsync(n, stream);

  // Launch the first kernel
  CH_Cuda::g_device.launch(vecAddToMdl->vvAddTo, numBlock, numThread, 0, stream,
                           vecAddToMdl->x.device, vecAddToMdl->y.device);

  // Launch the second kernel
  CH_Cuda::g_device.launch(vecAddToMdl->vsAddTo, numBlock, numThread, 0, stream,
                           vecAddToMdl->x.device, 2.0f);

  // Get the results
  vecAddToMdl->x.copyDtoHAsync(n, stream);

  // Barrier
  checkCudaErrors(cuStreamSynchronize(stream));

  // Destroy the stream
  checkCudaErrors(cuStreamDestroy(stream));

  // Destroy resources on the device
  checkCudaErrors(cuMemFree(vecAddToMdl->x.device));
  checkCudaErrors(cuMemFree(vecAddToMdl->y.device));

  // Check the results on the host
  for (int i = 0; i != n; ++i)
    {
      if (x[i] != 2*y[i] + 2) ++status;
    }

  // Destroy host-side memory
  checkCudaErrors(cuMemFreeHost(x));
  checkCudaErrors(cuMemFreeHost(y));

  // Unload the modules.  Normally you don't have to do this but we may want to
  // call testCUDAVectorAdd_libLevel() which loads the same module in a
  // different way.
  CH_Cuda::g_device.unloadModules();
  
#endif  /* CH_GPU */

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

