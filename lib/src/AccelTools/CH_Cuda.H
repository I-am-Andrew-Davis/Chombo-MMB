#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _CH_CUDADRIVER_H_
#define _CH_CUDADRIVER_H_


/******************************************************************************/
/**
 * \file
 *
 * \brief Support macros, types, and routines for the Cuda driver framework
 *
 * It is safe to include this file even if not compiling for the GPU
 *
 *//*+*************************************************************************/

#include <iostream>
#include <cstdio>
#include <cstdint>

#ifdef CH_GPU
#include "cuda.h"
#include <builtin_types.h>  /* From Cuda toolkit */
#endif

#include "CH_assert.H"
#include "parstream.H"

/*
  Notes:
    - tag __CUDACC__ means the file is processed with nvcc.  The class can
      be modified whether on the cpu or gpu.
    - tag __CUDA_ARCH__ is defined for device code.  This only works in
      functions with the __device__ qualifier.  It cannot be used in the
      class declaration.
*/

#undef HOSTDEVICE
#undef DEVICE
#ifdef __CUDACC__
#define HOSTDEVICE __host__ __device__
#define DEVICE __device__
#else
#define HOSTDEVICE
#define DEVICE
#endif


/*******************************************************************************
 */
/// Error-checking MACROS (declarations)
/*
 ******************************************************************************/

// From CUDA SDK samples...
// This will output the proper CUDA error strings in the event that a CUDA host
// call returns an error
#ifndef checkCudaErrors
#define checkCudaErrors(err) __checkCudaErrors(err, __FILE__, __LINE__)
#endif

#ifdef CH_GPU
void __checkCudaErrors(CUresult err, const char *file, const int line);
#endif


/*******************************************************************************
 */
/// Diagnostic MACROS
/**
 *  This quiets the warning: "dynamic initialization is not supported for a
 *  function-scope static __shared__ variable".  E.g., change
 *    __shared__ dyn::Array<Real, 3> arrA;
 *  to
 *    UNINITIALIZED__shared__(dyn::Array<Real, 3> arrA);
 *
 *  The meaning of the warning is that no member data is initialized, so make
 *  sure it is properly initialized before use!
 *
 *//*+*************************************************************************/

#ifdef __NVCC_DIAG_PRAGMA_SUPPORT__
  #define UNINITIALIZED__shared__(...)                                  \
  _Pragma("nv_diagnostic push")                                         \
  _Pragma("nv_diag_suppress static_var_with_dynamic_init")              \
  __shared__ __VA_ARGS__;                                               \
  _Pragma("nv_diagnostic pop")
#else
  #define UNINITIALIZED__shared__(...)                                  \
  _Pragma("diag_suppress static_var_with_dynamic_init")                 \
  __shared__ __VA_ARGS__
#endif


/*******************************************************************************
 */
/// Types
/*
 ******************************************************************************/

namespace CH_Cuda
{

#ifdef CH_GPU
using DevicePointer = CUdeviceptr;
using DeviceFunction = CUfunction;

/*
  The default stream identifier is set to CU_STREAM_PER_THREAD which does not
  synchronize with other streams.  In many asynchronous functions, this defines
  the default argument to the stream identifier.  Alternatives are an
  explicitly-managed stream or the legacy stream, which syncrhonizes with all
  other streams, can be selected with CU_STREAM_LEGACY.
 */
const CUstream c_defaultStream = CU_STREAM_PER_THREAD;
#else
using DevicePointer = std::uintptr_t;
#endif

/// A handle for dual pointers to host and device memory
template <typename T>
struct SymbolPair
{
  /// Default constructor
  HOSTDEVICE SymbolPair()
    :
    host(nullptr)
#ifdef CH_GPU
    ,
    device(0)
#endif
    { }

  /// Construct with arguments
  HOSTDEVICE SymbolPair(T *const a_host, DevicePointer a_device = 0)
    :
    host(a_host)
#ifdef CH_GPU
    ,
    device(a_device)
#endif
    { }

  SymbolPair(const SymbolPair&) = default;
  SymbolPair(SymbolPair&&) = default;
  SymbolPair& operator=(const SymbolPair&) = default;
  SymbolPair& operator=(SymbolPair&&) = default;

#ifdef CH_GPU
  /// Synchronous copy from host to device
  /** \param[in]  a_num   Number of elements of type T to copy.
   */
  void copyHtoD(const size_t a_num) const
    {
      checkCudaErrors(cuMemcpyHtoD(device,
                                   host,
                                   a_num*sizeof(T)));
    }

  /// Asynchronous copy from host to device
  /** \param[in]  a_num   Number of elements of type T to copy.
   *  \param[in]  a_stream
   *                      Stream identifier
   */
  void copyHtoDAsync(const size_t a_num,
                     CUstream     a_stream = c_defaultStream) const
    {
      checkCudaErrors(cuMemcpyHtoDAsync(device,
                                        host,
                                        a_num*sizeof(T),
                                        a_stream));
    }

  /// Synchronous copy from device to host
  /** \param[in]  a_num   Number of elements of type T to copy.
   */
  void copyDtoH(const size_t a_num) const
    {
      checkCudaErrors(cuMemcpyDtoH(host,
                                   device,
                                   a_num*sizeof(T)));
    }

  /// Asynchronous copy from device to host
  /** \param[in]  a_num   Number of elements of type T to copy.
   *  \param[in]  a_stream
   *                      Stream identifier
   */
  void copyDtoHAsync(const size_t a_num,
                     CUstream     a_stream = c_defaultStream) const
    {
      checkCudaErrors(cuMemcpyDtoHAsync(host,
                                        device,
                                        a_num*sizeof(T),
                                        a_stream));
    }
#endif

//--Data members

  T *host;
#ifdef CH_GPU
  DevicePointer device;
#endif
};

using copy_to_device = std::integral_constant<unsigned, 0>;


/*******************************************************************************
 */
/// Configuration for a single function
/**
 *//*+*************************************************************************/

#ifdef CH_GPU
struct DeviceFunctionConfig
{
  /// Default constructor
  DeviceFunctionConfig()
    :
    m_cacheConfig(CU_FUNC_CACHE_PREFER_NONE),
    m_sharedMemBankSize(CU_SHARED_MEM_CONFIG_DEFAULT_BANK_SIZE),
    m_sharedMemBytes(0)
    {
      m_numThread.x = 0;
      m_numThread.y = 0;
      m_numThread.z = 0;
    }
  /// Constructor (string name)
  DeviceFunctionConfig(
    const std::string&   a_name,
    const dim3&          a_numThread,
    const CUfunc_cache   a_cacheConfig = CU_FUNC_CACHE_PREFER_NONE,
    const CUsharedconfig a_sharedMemBankSize =
    CU_SHARED_MEM_CONFIG_DEFAULT_BANK_SIZE,
    const unsigned a_sharedMemBytes = 0)
    :
    m_name(a_name),
    m_numThread(a_numThread),
    m_cacheConfig(a_cacheConfig),
    m_sharedMemBankSize(a_sharedMemBankSize),
    m_sharedMemBytes(a_sharedMemBytes)
    { }
  /// Constructor (string name) (threads in x-direction only)
  DeviceFunctionConfig(
    const std::string&   a_name,
    const int            a_numThread = 1,
    const CUfunc_cache   a_cacheConfig = CU_FUNC_CACHE_PREFER_NONE,
    const CUsharedconfig a_sharedMemBankSize =
    CU_SHARED_MEM_CONFIG_DEFAULT_BANK_SIZE,
    const unsigned a_sharedMemBytes = 0)
    :
    m_name(a_name),
    m_cacheConfig(a_cacheConfig),
    m_sharedMemBankSize(a_sharedMemBankSize),
    m_sharedMemBytes(a_sharedMemBytes)
    {
      m_numThread.x = a_numThread;
      m_numThread.y = 1;
      m_numThread.z = 1;
    }
  /// Constructor (char* name)
  DeviceFunctionConfig(
    const char *const    a_name,
    const dim3&          a_numThread,
    const CUfunc_cache   a_cacheConfig = CU_FUNC_CACHE_PREFER_NONE,
    const CUsharedconfig a_sharedMemBankSize =
    CU_SHARED_MEM_CONFIG_DEFAULT_BANK_SIZE,
    const unsigned a_sharedMemBytes = 0)
    :
    m_name(a_name),
    m_numThread(a_numThread),
    m_cacheConfig(a_cacheConfig),
    m_sharedMemBankSize(a_sharedMemBankSize),
    m_sharedMemBytes(a_sharedMemBytes)
    { }
  /// Constructor (char* name) (threads in x-direction only)
  DeviceFunctionConfig(
    const char *const    a_name,
    const int            a_numThread = 1,
    const CUfunc_cache   a_cacheConfig = CU_FUNC_CACHE_PREFER_NONE,
    const CUsharedconfig a_sharedMemBankSize =
    CU_SHARED_MEM_CONFIG_DEFAULT_BANK_SIZE,
    const unsigned a_sharedMemBytes = 0)
    :
    m_name(a_name),
    m_numThread(a_numThread),
    m_cacheConfig(a_cacheConfig),
    m_sharedMemBankSize(a_sharedMemBankSize),
    m_sharedMemBytes(a_sharedMemBytes)
    {
      m_numThread.x = a_numThread;
      m_numThread.y = 1;
      m_numThread.z = 1;
    }

//--Data

  std::string m_name;                 ///< Name of the function
  dim3 m_numThread;                   ///< Number of threads.  Note, launches
                                      ///< can be made with a different number
                                      ///< of threads.  This is particularly
                                      ///< useful if a function is designed to
                                      ///< use a specific number of threads
                                      ///< (which is actually pretty common).
                                      ///< Set to 0 to disable
  CUfunc_cache m_cacheConfig;         ///< Cache preferences for the kernel.
                                      ///< CU_FUNC_CACHE_PREFER_NONE: no
                                      ///<   preference for shared memory or L1
                                      ///< CU_FUNC_CACHE_PREFER_SHARED: prefer
                                      ///<   larger shared memory and smaller L1
                                      ///<   cache
                                      ///< CU_FUNC_CACHE_PREFER_L1: prefer
                                      ///<   larger L1 cache and smaller shared
                                      ///<   memory
                                      ///< CU_FUNC_CACHE_PREFER_EQUAL: prefer
                                      ///<   equal sized L1 cache and shared
                                      ///<   memory
  CUsharedconfig m_sharedMemBankSize; ///< Preferences for shared memory bank
                                      ///< size
                                      ///< CU_SHARED_MEM_CONFIG_DEFAULT_BANK_
                                      ///< SIZE: set bank width to the default
                                      ///<   initial setting (currently, four
                                      ///<   bytes).
                                      ///< CU_SHARED_MEM_CONFIG_FOUR_BYTE_BANK_
                                      ///< SIZE: set shared memory bank width to
                                      ///<   be natively four bytes.
                                      ///< CU_SHARED_MEM_CONFIG_EIGHT_BYTE_BANK_
                                      ///< SIZE: set shared memory bank width to
                                      ///<   be natively eight bytes.
  unsigned m_sharedMemBytes;          ///< Number of bytes to allocate in shared
                                      ///< memory
  DeviceFunction m_CUfunction;        ///< Function handle (do not set, used
                                      ///< internally)
};


/*******************************************************************************
 */
/// Configuration for a single symbol
/**
 *//*+*************************************************************************/

struct DeviceSymbolConfig
{
  /// Default constructor
  DeviceSymbolConfig()
    { }
  /// Constructor (string name)
  DeviceSymbolConfig(const std::string& a_name)
    :
    m_name(a_name)
    { }
  /// Constructor (char* name)
  DeviceSymbolConfig(const char *const a_name)
    :
    m_name(a_name)
    { }

//--Movement

  /// Synchronous copy from host to device
  /** \param[in]  a_host  Object on the host
   */
  template <typename T>
  void copyHtoD(const T& a_host) const noexcept
    {
      CH_assert(sizeof(T) == m_numBytes);
      checkCudaErrors(cuMemcpyHtoD(m_CUpointer,
                                   &a_host,
                                   sizeof(T)));
    }

  /// Asynchronous copy from host to device
  /** \param[in]  a_host  Object on the host
   *  \param[in]  a_stream
   *                      Stream identifier
   */
  template <typename T>
  void copyHtoDAsync(const T& a_host,
                     CUstream a_stream = c_defaultStream) const noexcept
    {
      CH_assert(sizeof(T) == m_numBytes);
      checkCudaErrors(cuMemcpyHtoDAsync(m_CUpointer,
                                        &a_host,
                                        sizeof(T),
                                        a_stream));
    }

  /// Synchronous copy from device to host
  /** \param[out] a_host  Object on the host
   */
  template <typename T>
  void copyDtoH(T& a_host) const noexcept
    {
      CH_assert(sizeof(T) == m_numBytes);
      checkCudaErrors(cuMemcpyDtoH(&a_host,
                                   m_CUpointer,
                                   sizeof(T)));
    }

  /// Asynchronous copy from device to host
  /** \param[out] a_host  Object on the host
   *  \param[in]  a_stream
   *                      Stream identifier
   */
  template <typename T>
  void copyDtoHAsync(T&       a_host,
                     CUstream a_stream = c_defaultStream) const noexcept
    {
      CH_assert(sizeof(T) == m_numBytes);
      checkCudaErrors(cuMemcpyDtoHAsync(&a_host,
                                        m_CUpointer,
                                        sizeof(T),
                                        a_stream));
    }

//--Data

  std::string m_name;                 ///< Name of the symbol
  size_t m_numBytes;                  ///< Number of bytes used on device
                                      ///< (do not set, set internally)
  DevicePointer m_CUpointer;          ///< Pointer to the symbol on the device
                                      ///< (do not set, set internally)
};


/*******************************************************************************
 *
 * Message output support
 *
 ******************************************************************************/

/// Write message
# define WRITE_MSG( msg, header )                                       \
  {                                                                     \
    pout() << header << msg << std::endl;                               \
  }
const char *const msgHeader1 = "=> ";
const char *const msgHeader2 = "   ";
const int msgLabelWidth = 40;
const int msgBufferSize = 512;
extern char msgBuffer[msgBufferSize];

/// Finds a better SI prefix for a count
const char* reduceUnitCount(const int a_unitType, float& a_count);
#endif  /* CH_GPU */

}  // End of namespace Cuda


/*******************************************************************************
 */
/// Error-checking MACROS (definitinos)
/*
 ******************************************************************************/

#ifdef CH_GPU

// From CUDA SDK samples...
// This will output the proper CUDA error strings in the event that a CUDA host
// call returns an error
#ifndef checkCudaErrors
#define checkCudaErrors(err) __checkCudaErrors(err, __FILE__, __LINE__)
#endif

inline void __checkCudaErrors(CUresult err, const char *file, const int line)
{
  if (CUDA_SUCCESS != err)
    {
      const char* errorName = nullptr;
      const char* errorStr = nullptr;
      cuGetErrorName(err, &errorName);
      cuGetErrorString(err, &errorStr);
#ifdef __CUDACC__
      printf("Cuda driver API error %04d, %s: \"%s\" from file <%s>, "
             "line %i", err, errorName, errorStr, file, line);
#else
      std::snprintf(CH_Cuda::msgBuffer, CH_Cuda::msgBufferSize,
                    "Cuda driver API error %04d, %s: \"%s\" from file <%s>, "
                    "line %i", err, errorName, errorStr, file, line);
      MayDay::Error(CH_Cuda::msgBuffer);
#endif
    }                                                          
}
#endif  /* CH_GPU */

#endif  /* ! defined _CH_CUDADRIVER_H_ */
