#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cstring>
#include <set>
#include <map>
#include <functional>
#include <utility>
#include <unordered_map>

#if CXXSTD>=14
#define D_DECL6(a,b,c,d,e,f) a,b,c
#include "StcVector.H"
#endif
#include "CH_Timer.H"
#include "parstream.H"
#include "CH_Hash.H"
#ifdef CH_MPI
#include <mpi.h>
#endif
#ifdef CH_GPU
  #include "testHash_CUX.H"
  #include "CudaDriver.H"
#endif

#include "UsingBaseNamespace.H"

/// Prototypes:

void
parseTestOptions(int argc, char* argv[]);

int
testCHHash();

/// Global variables for handling output:
static const char *pgmname = "testHash" ;
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
  int ret = testCHHash();

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

#if CXXSTD>=14
// Dummy IntVect class with SpaceDim=3
constexpr int c_SpaceDim = 3;
using IntVect = stc::IVec<c_SpaceDim>;
constexpr auto IntVect_zero = stc::make_IVec<c_SpaceDim>::zero();
constexpr auto IntVect_unit = stc::make_IVec<c_SpaceDim>::unit();

namespace CH_Hash
{

/// Remark that IntVect can be byte-hashed and give bytes to hash
template <>
struct isHashable<IntVect> : public std::true_type
{
  static_assert(sizeof(IntVect) == c_SpaceDim*sizeof(int),
                "Cannot guarantee that IntVect is not padded");
  static constexpr int c_hashSize = sizeof(IntVect);
};

/// Remark that int[3] can be byte-hashed and give bytes to hash
template <>
struct isHashable<int[3]> : public std::true_type
{
  static constexpr int c_hashSize = 3*sizeof(int);
};

}
#endif

int testCHHash()
{
  CH_TIMERS("testHash");
  int status = 0;
#if CXXSTD>=14

  CH_Hash::google_CityHash64<IntVect> hasher;
  IntVect iv{ 1, 2, 3 };
  if (verbose)
    {
      pout() << hasher(iv) << std::endl;
    }
  auto hash0A = CH_Hash::hash_google_CityHash64(IntVect_zero);
  int zero[3] = { 0, 0, 0 };
  auto hash0B = CH_Hash::hash_google_CityHash64(zero);
  if (hash0A != hash0B) ++status;

//--Set values of interest here

  constexpr size_t keyDim = 64;
  constexpr size_t numKey = keyDim*keyDim*keyDim;
  constexpr double loadFactor = 0.7;
  constexpr size_t setSize = (size_t)(numKey/loadFactor);
  // size_t divider = std::numeric_limits<size_t>::max()/setSize;
  // while (std::numeric_limits<size_t>::max()/divider > setSize - 1) ++divider;

  using KeyMap = std::map<uint64_t, int>;
  KeyMap keys;
  stc::nestedLoop(IntVect_zero, (keyDim-1)*IntVect_unit,
                  [&]
                  (const IntVect& a_iv)
                  {
                    std::pair<KeyMap::iterator, bool> ins =
                      // keys.insert(std::make_pair(hasher(a_iv) / divider, 1));
                      keys.insert(std::make_pair(hasher(a_iv) % setSize, 1));
                    if (!ins.second)  // Key already existed
                      {
                        ++(ins.first->second);
                      }
                  });
  int maxCollision = 0;
  for (auto& kv : keys)
    {
      maxCollision = std::max(maxCollision, kv.second);
    }
  if (verbose)
    {
      pout() << "Max number of collisions: " << maxCollision << std::endl;
    }
  std::vector<int> cntCollision(maxCollision+1, 0);
  for (auto& kv : keys)
    {
      ++cntCollision[kv.second];
    }
  if (verbose)
    {
      pout() << "Collision frequency:\n";
    }
  size_t totalk = 0;
  for (int i = 1; i <= maxCollision; ++i)
    {
      if (verbose)
        {
          pout() << i << ": " << cntCollision[i] << std::endl;
        }
      totalk += i*cntCollision[i];
    }
  if (numKey != totalk) ++status;
  double collisionlessRate = 100*((double)cntCollision[1])/numKey;
  if (collisionlessRate < 49.) ++status;
  if (verbose)
    {
      pout() << "Collisionless %: " << collisionlessRate << std::endl;
    }
  CH_TIMER("Hash:City64", t1);
  CH_START(t1);
  uint64_t sum = 0;
  stc::nestedLoop(IntVect_zero, 256*IntVect_unit,
                  [&]
                  (const IntVect& a_iv)
                  {
                    sum += hasher(a_iv);
                  });
  CH_STOP(t1);
  if (verbose)
    {
      pout() << "sum: " << sum << std::endl;
    }

  CH_Hash::google_CityHash32_r64<IntVect> hasher32;
  keys.clear();
  stc::nestedLoop(IntVect_zero, (keyDim-1)*IntVect_unit,
                  [&]
                  (const IntVect& a_iv)
                  {
                    std::pair<KeyMap::iterator, bool> ins =
                      // keys.insert(std::make_pair(hasher32(a_iv) / divider, 1));
                      keys.insert(std::make_pair(hasher32(a_iv) % setSize, 1));
                    if (!ins.second)  // Key already existed
                      {
                        ++(ins.first->second);
                      }
                  });
  maxCollision = 0;
  for (auto& kv : keys)
    {
      maxCollision = std::max(maxCollision, kv.second);
    }
  if (verbose)
    {
      pout() << "Max number of collisions: " << maxCollision << std::endl;
    }
  cntCollision.assign(maxCollision+1, 0);
  for (auto& kv : keys)
    {
      ++cntCollision[kv.second];
    }
  if (verbose)
    {
      pout() << "Collision frequency:\n";
    }
  totalk = 0;
  for (int i = 1; i <= maxCollision; ++i)
    {
      if (verbose)
        {
          pout() << i << ": " << cntCollision[i] << std::endl;
        }
      totalk += i*cntCollision[i];
    }
  if (numKey != totalk) ++status;
  collisionlessRate = 100*((double)cntCollision[1])/numKey;
  if (collisionlessRate < 49.) ++status;
  if (verbose)
    {
      pout() << "Collisionless %: " << collisionlessRate << std::endl;
    }
  CH_TIMER("Hash:City32", t2);
  CH_START(t2);
  sum = 0;
  stc::nestedLoop(IntVect_zero, 256*IntVect_unit,
                  [&]
                  (const IntVect& a_iv)
                  {
                    sum += hasher32(a_iv);
                  });
  CH_STOP(t2);
  if (verbose)
    {
      pout() << "sum: " << sum << std::endl;
    }
  
#endif  /* CXXSTD>=14 */

  {
    std::unordered_map<int, double, CH_Hash::google_CityHash<int>> hashmap{};
    hashmap[3]    = 1.5;
    hashmap[1045] = 2.5;
    if (hashmap[3   ] != 1.5) ++status;
    if (hashmap[1045] != 2.5) ++status;
    if (hashmap.find(8) != hashmap.end()) ++status;
  }

//--Test hashing on the GPU

#ifdef CH_GPU
  {
    unsigned cpuHash = CH_Hash::hash_google_CityHash32(42);

#if 0

//------------------------------------------------------------------------------
// This is the low-level approach also demonstrated in testCudaVectorAdd.cpp
// It should not be used for the application and Chombo

    class HashModule : public CH_Cuda::DriverModule
    {
    public:
      CH_Cuda::DeviceFunction func_hash;
      CH_Cuda::DevicePtr symb_hash;
    };

    // Initialize module
    HashModule *hashMdl =
      CH_Cuda::g_device.loadModule<HashModule>("testHash", testHash);

    // Retrieve functions (FIXME** should be made automatic)
    checkCudaErrors(cuModuleGetFunction(&hashMdl->func_hash,
                                        hashMdl->m_module,
                                        "hash"));

    // Retrieve symbols (FIXME** should be made automatic)
    size_t numBytes;
    checkCudaErrors(cuModuleGetGlobal(&hashMdl->symb_hash,
                                      &numBytes,
                                      hashMdl->m_module,
                                      "g_hash"));

    // Create the stream
    CUstream stream;
    checkCudaErrors(cuStreamCreate(&stream, 0));

    // Launch the first kernel
    dim3 numBlock;
    dim3 numThread;
    CH_Cuda::g_device.launch(hashMdl->func_hash, numBlock, numThread, 0, stream);

    // Transfer to the host
    unsigned gpuHash;
    checkCudaErrors(cuMemcpyDtoHAsync(&gpuHash,
                                      hashMdl->symb_hash,
                                      sizeof(unsigned),
                                      stream));

    // Barrier
    checkCudaErrors(cuStreamSynchronize(stream));

    // Destroy the stream
    checkCudaErrors(cuStreamDestroy(stream));
//------------------------------------------------------------------------------

#else

    /// Handles to module global functions for this application
    enum class App_Global
      {
        hash = 0
      };

    /// Handles to module symbols for this application
    enum class App_Symbol
      {
        g_hash = 0
      };

    // If there are functions or symbols from standard Chombo libraries, you
    // would mark them here using, e.g., BaseTools_lib
    // Here we mark the functions and symbols for this application (test).
    CH_Cuda::g_device.configureChomboLibFunctions(
      CH_Cuda::Application_lib,
      []
      ()
        {
          std::vector<CH_Cuda::DeviceFunctionConfig> config;
          config.reserve(1);
          // App_Global 0
          config.emplace_back("hash", 1);  // 1 thread
          return config;
        });

    CH_Cuda::g_device.configureChomboLibSymbols(
      CH_Cuda::Application_lib,
      []
      ()
        {
          std::vector<CH_Cuda::DeviceSymbolConfig> config;
          config.reserve(1);
          // App_Symbol 0
          config.emplace_back("g_hash");
          return config;
        });

    // Now load the module which automatically loads the functions and symbols
    CH_Cuda::g_device.loadApplicationModule("testHash", testHash);

    // Create the stream
    CUstream stream;
    checkCudaErrors(cuStreamCreate(&stream, 0));

    // Launch the first kernel (use threads specified for the function)
    dim3 numBlock;
    CH_Cuda::g_device.launch(CH_Cuda::Application_lib,
                             App_Global::hash,
                             numBlock,
                             stream);

    // Transfer to the host
    unsigned gpuHash;
    CH_Cuda::g_device.copySymbolToHostAsync(CH_Cuda::Application_lib,
                                            App_Symbol::g_hash,
                                            gpuHash,
                                            stream);

    // Barrier
    checkCudaErrors(cuStreamSynchronize(stream));

    // Destroy the stream
    checkCudaErrors(cuStreamDestroy(stream));
#endif

    // std::cout << "On CPU: " << cpuHash << std::endl;
    if (cpuHash != gpuHash) ++status;
  }
#endif

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
