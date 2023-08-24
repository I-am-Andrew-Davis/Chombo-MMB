#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _CUDADRIVER_H_
#define _CUDADRIVER_H_


/******************************************************************************/
/**
 * \file
 *
 * \brief Device initialization and support for using the Cuda driver API
 *
 * This file should only be included under the protection of
 * #ifdef CH_GPU
 *   #include CudaDriver.H
 * #endif CH_GPU
 *
 *//*+*************************************************************************/

#include <vector>
#include <map>

#include "CH_Cuda.H"
#include "StcVector.H"

//--Forward declarations

// enum class ArrayClassIndex;

// template <typename T, ArrayClassIndex I>
// class CUDAArrayAlloc;

// namespace CH_Hash
// {
// template <typename T>
// struct google_CityHash;

// template<typename T>
// struct equal_to;

// template <typename Key,
//           typename Mapped,
//           typename Hash = google_CityHash<Key>,
//           typename KeyEqual = equal_to<Key>,
//           typename Alloc = CUDAArrayAlloc<std::pair<Key, Mapped>,
//                                           ArrayClassIndex::HashIVS>>
// class OpenAddressHashTable;
// }

namespace CH_Cuda
{

  // Forward declaration
  class Driver;

/// Handle to a global instance of the driver.
/** This supports 1 device for 1 or more (using Cuda MPS) processes.  For
 *  multiple devices on a process, it is probalby best to have a separate
 *  CudaDriver class for each device.  However, a better alternative would be to
 *  use more processes so one still has one or more processes per device.
 */
  extern Driver g_device;

//--Handles to the Chombo and application libraries in primary module 0

  // Basic
  static constexpr int BaseTools_lib             = 0;
  static constexpr int BoxTools_lib              = 1;
  static constexpr int AMRTools_lib              = 2;
  static constexpr int AMRTimeDependent_lib      = 3;
  static constexpr int AMRElliptic_lib           = 4;
  // High-order
  static constexpr int HOAMRTools_lib            = 5;
  static constexpr int HOAMRTimeDependent_lib    = 6;
  // Mapped geometry
  static constexpr int HOMappedGeometry_lib      = 7;
  static constexpr int HOMappedConservation_lib  = 8;
  static constexpr int Application_lib           = 9;
  static constexpr int c_numModule0Lib = Application_lib + 1;

/*******************************************************************************
 */
/// Help with alignments (use for adding parameters to functions)
/**
 *//*+*************************************************************************/

#if 0
template <typename T>
inline int
addParam(CUfunction& a_func,
         int         a_offset,
         T&          a_param)
{
  MayDay::Error("Unsupported parameter type");
  return 0;
}

// Specialized for int
template <>
inline int
addParam(CUfunction& a_func,
         int         a_offset,
         int&        a_param)
{
  a_offset = (a_offset + __alignof(a_param) - 1) & ~(__alignof(a_param) - 1);
  CU_SAFE_CALL(cuParamSeti(a_func, a_offset, a_param));
  a_offset += sizeof(a_param);
  return a_offset;
}

// Specialized for float
template <>
inline int
addParam(CUfunction& a_func,
         int         a_offset,
         float&      a_param)
{
  a_offset = (a_offset + __alignof(a_param) - 1) & ~(__alignof(a_param) - 1);
  CU_SAFE_CALL(cuParamSetf(a_func, a_offset, a_param));
  a_offset += sizeof(a_param);
  return a_offset;
}

// Specialized for device pointers
template <>
inline int
addParam(CUfunction&  a_func,
         int          a_offset,
         CUdeviceptr& a_param)
{
  void* ptr = (void*)(size_t)a_param;
  a_offset = (a_offset + __alignof(ptr) - 1) & ~(__alignof(ptr) - 1);
  CU_SAFE_CALL(cuParamSetv(a_func, a_offset, &ptr, sizeof(ptr)));
  a_offset += sizeof(ptr);
  return a_offset;
}
#endif


/*******************************************************************************
 */
/// Base class for Cuda driver module data
/** Any description of a module must be derived from this class.
 *
 *//*+*************************************************************************/

class DriverModule
{
public:
  CUmodule m_module;
  std::string m_name;
};


/*******************************************************************************
 */
/// Initializes generic CUDA components
/** Handles to devices, contexts, and modules are stored here.
 *
 *//*+*************************************************************************/

// Note: T must be deduced meaning you must use helpers converting_t and
//       converting_builder
template <typename T>
struct Converting
{
  using type = T;
  // Passthrough of constructor argument
  static T&& builder(T&& a_arg)
    {
      return std::forward<T>(a_arg);
    }
};

// template <typename Key,
//           typename Mapped,
//           typename Hash,
//           typename KeyEqual>
// struct Converting<CH_Hash::OpenAddressHashTable<
//                     Key, Mapped, Hash, KeyEqual,
//                     CUDAArrayAlloc<std::pair<Key, Mapped>,
//                                    ArrayClassIndex::HashIVS>>&>
// {
//   using type = CH_Hash::OpenAddressHashTable<Key, Mapped, Hash, KeyEqual,
//                                              CUDAArrayAlloc<std::pair<Key, Mapped>,
//                                                             ArrayClassIndex::HashIVS>>;
//   static type builder(type& a_arg)
//     {
//       type local(copy_to_device{}, a_arg);
//       return local;
//     }
// };

// template <>
// struct Converting<int&>
// {
//   using type = int;
//   static type builder(int& a_arg)
//     {
//       type local = a_arg + 2;
//       return local;
//     }
// };

template <typename T>
using converting_type = typename Converting<T>::type;

template <typename T>
typename Converting<T>::type converting_builder(T&& a_arg)
{
  return Converting<T>::builder(std::forward<T>(a_arg));
}

template <typename... Args>
struct LValArgs
{
  LValArgs(Args&&... a_args)
    :
    m_args(converting_builder(std::forward<Args>(a_args))...)
    {
      stc::forEachTElement<sizeof...(Args)>(*this);
    }

  template <stc::array_size_type I>
  void operator()(const int)
    {
      m_addrArgs[I] = &std::get<I>(m_args);
    }

  std::tuple<converting_type<Args>...> m_args;
  void* m_addrArgs[sizeof...(Args)];
};

/// Zero-sized
template <>
struct LValArgs<>
{
  LValArgs()
    :
    m_addrArgs(nullptr)
    { }

  void** m_addrArgs;
};

class Driver
{
private:

  /// Numbers (bit flags) for CUDA components
  enum class InitComponentFlag : unsigned
  {
    device     = (1 << 0),
    context    = (1 << 1),
    context0   = (1 << 2),
    anyContext = context + context0,
    module0    = (1 << 3),
    all        = (1 << 4) - 1
  };

public:

  /// Default constructor initializes all cuda components as undefined
  Driver();

  /// Destructor undefines defined cuda components
  ~Driver();

  /// Test if a cuda component has been defined by CudaDriver
  bool testInitCompFlag(const InitComponentFlag a_icf) const;

  /// Initialize the device
  void initialize(
    const int            a_verbose = 0,
    const CUfunc_cache   a_cachePreference = CU_FUNC_CACHE_PREFER_SHARED,
    const CUsharedconfig a_sharedMemBankSize =
      CU_SHARED_MEM_CONFIG_DEFAULT_BANK_SIZE);

  /// Unload all external modules
  void unloadModules();

  /// Unload primary module 0
  void unloadModule0();

  /// Finalize the state
  void finalize();

  template <typename T>
  void configureChomboLib(const int a_idxLib, const T& a_libData)
    {
      m_functionConfigs[a_idxLib] = a_libData.functions();
      m_symbolConfigs[a_idxLib] = a_libData.symbols();
    }

  template <typename T>
  void configureChomboLibFunctions(const int a_idxLib, const T& a_libData)
    {
      m_functionConfigs[a_idxLib] = a_libData();
    }

  template <typename T>
  void configureChomboLibSymbols(const int a_idxLib, const T& a_libData)
    {
      m_symbolConfigs[a_idxLib] = a_libData();
    }

  /// Load an application module
  void loadApplicationModule(const char *const a_moduleName,
                             const void*       a_embeddedBin,
                             const char *const a_moduleDir  = nullptr,
                             const char *const a_searchPath = nullptr);

  /// Load a module (GPU executable) - T is derived from class DriverModule
  /** This should only be used for 3rd-party libraries having a CUDA module.
   *  The standard way is to include a *_CUX.H file which contains the byte
   *  string for an embedded CUDA binary.  Alternatively, the module may be in
   *  a separate file.
   *  See member function baseLoadModule for more information on paths.
   *
   *  WARNING: you must load your functions and symbols separately!
   *
   *  \param[in]  a_moduleName
   *                      Name of the module.  If the module is loaded from
   *                      a file, this should be the name of the file.  If the
   *                      module is embedded, the name is only used as a
   *                      reference and can be anything.  It is suggested to use
   *                      a_moduleName = "a_embeddedBin"
   *  \param[in]  a_embeddedBin
   *                      Byte string for an embedded CUDA binary.  In the
   *                      Chombo build system, this is found in the *_CUX.H
   *                      files and has the same name as *.  It is suggested
   *                      to use a module name of "*" for parameter
   *                      a_moduleName.
   *  \param[in]  a_moduleDir
   *                      A directory (to append to a_searchPath) in
   *                      which the module might be found.  You can
   *                      specify the directory in a_searchPath and
   *                      leave this NULL -- it is only here to provide
   *                      more flexibility and for use with
   *                      'path_to_current_process'
   *  \param[in]  a_searchPath
   *                      Path to search for the file.
   */
  template <typename T>
  T* loadModule(const char *const a_moduleName,
                const void*       a_embeddedBin,
                const char *const a_moduleDir  = nullptr,
                const char *const a_searchPath = nullptr)
    {
      CH_assert(a_moduleName != nullptr);
      m_moduleVec.push_back(new T());
      m_moduleVec.back()->m_name = a_moduleName;
      int moduleStatus = 1;
      if (a_embeddedBin != nullptr)
        {
          moduleStatus = baseLoadEmbeddedModule(&m_moduleVec.back()->m_module,
                                                a_moduleName,
                                                a_embeddedBin);
        }
      if (moduleStatus)  // 0 = previous success
        {
          baseLoadModule(&m_moduleVec.back()->m_module,
                         a_moduleName,
                         a_moduleDir,
                         a_searchPath);
        }
      return static_cast<T*>(m_moduleVec.back());
    }


  /// Launch a kernel
  template <typename... Args>
  void launch(CUfunction   a_kernel,
              const dim3&  a_numBlock,
              const dim3&  a_numThread,
              unsigned     a_sharedMemBytes,
              CUstream     a_stream,
              Args&&...    a_args)
    {
      LValArgs<Args...> lvArgs(std::forward<Args>(a_args)...);
      // for (int i = 0; i != sizeof...(Args); ++i)
      //   {
      //     std::cout << "Addr " << i << ": " << lvArgs.m_addrArgs[i]
      //               << std::endl;
      //   }
      // std::cout << a_kernel << std::endl;
      // std::cout << a_numBlock.x << std::endl;
      // std::cout << a_numBlock.y << std::endl;
      // std::cout << a_numBlock.z << std::endl;
      // std::cout << a_numThread.x << std::endl;
      // std::cout << a_numThread.y << std::endl;
      // std::cout << a_numThread.z << std::endl;
      // std::cout << a_sharedMemBytes << std::endl;
      // std::cout << a_stream << std::endl;
      // std::cout << lvArgs.m_addrArgs << std::endl;
      // std::cout << lvArgs.m_addrArgs[0] << std::endl;
      // std::cout << lvArgs.m_addrArgs[1] << std::endl;
      checkCudaErrors(
        cuLaunchKernel(a_kernel,
                       a_numBlock.x, a_numBlock.y, a_numBlock.z,
                       a_numThread.x, a_numThread.y, a_numThread.z,
                       a_sharedMemBytes,
                       a_stream,
                       lvArgs.m_addrArgs,
                       nullptr));
    }

  /// Launch a kernel from a function in the primary module
  template <typename F, typename... Args>
  void launch(const int   a_idxLib,
              const F     a_idxFunc,
              const dim3& a_numBlock,
              const dim3& a_numThread,
              CUstream    a_stream,
              Args&&...   a_args)
    {
      CH_assert(testInitCompFlag(InitComponentFlag::module0));
      const DeviceFunctionConfig& funcConf =
        m_functionConfigs[a_idxLib][static_cast<int>(a_idxFunc)];
      launch(funcConf.m_CUfunction,
             a_numBlock,
             a_numThread,
             funcConf.m_sharedMemBytes,
             a_stream,
             std::forward<Args>(a_args)...);
    }

  /// Launch a kernel from a function in the primary module with preset block
  /// size
  template <typename F, typename... Args>
  void launch(const int   a_idxLib,
              const F     a_idxFunc,
              const dim3& a_numBlock,
              CUstream    a_stream,
              Args&&...   a_args)
    {
      CH_assert(testInitCompFlag(InitComponentFlag::module0));
      const DeviceFunctionConfig& funcConf =
        m_functionConfigs[a_idxLib][static_cast<int>(a_idxFunc)];
      launch(funcConf.m_CUfunction,
             a_numBlock,
             funcConf.m_numThread,
             funcConf.m_sharedMemBytes,
             a_stream,
             std::forward<Args>(a_args)...);
    }

  /// Get a symbol config from the primary module
  /** \tparam     S       Type of enum for symbol
   *  \param[in]  a_idxLib
   *                      Index of library in which to find symbol (e.g.,
   *                      c_BaseTools_lib or c_Application_lib
   *  \param[in]  a_idxSymb
   *                      Enum giving handle to the symbol
   */
  template <typename S>
  const DeviceSymbolConfig& getSymbol(const int a_idxLib,
                                      const S   a_idxSymb) noexcept
    {
      return m_symbolConfigs[a_idxLib][static_cast<int>(a_idxSymb)];
    }

//--Copy symbols to/from the device

  /// Synchronous copy from host to device
  /** \tparam     S       Type of enum for symbol
   *  \tparam     T       Type of object (should match symbol type in *.cu file)
   *  \param[in]  a_idxLib
   *                      Index of library in which to find symbol (e.g.,
   *                      c_BaseTools_lib or c_Application_lib
   *  \param[in]  a_idxSymb
   *                      Enum giving handle to the symbol
   *  \param[out] a_val   Object on the host
   */
  template <typename S, typename T>
  void copySymbolToHost(const int a_idxLib,
                        const S   a_idxSymb,
                        T&        a_val) noexcept
    {
      CH_assert(testInitCompFlag(InitComponentFlag::module0));
      getSymbol(a_idxLib, a_idxSymb).copyDtoH(a_val);
    }

  /// Asynchronous copy from host to device
  /** \tparam     S       Type of enum for symbol
   *  \tparam     T       Type of object (should match symbol type in *.cu file)
   *  \param[in]  a_idxLib
   *                      Index of library in which to find symbol (e.g.,
   *                      c_BaseTools_lib or c_Application_lib
   *  \param[in]  a_idxSymb
   *                      Enum giving handle to the symbol
   *  \param[out] a_val   Object on the host
   *  \param[in]  a_stream
   *                      Stream identifier
   */
  template <typename S, typename T>
  void copySymbolToHostAsync(const int a_idxLib,
                             const S   a_idxSymb,
                             T&        a_val,
                             CUstream  a_stream = c_defaultStream) noexcept
    {
      CH_assert(testInitCompFlag(InitComponentFlag::module0));
      getSymbol(a_idxLib, a_idxSymb).copyDtoHAsync(a_val, a_stream);
    }

  /// Synchronous copy from device to host
  /** \tparam     S       Type of enum for symbol
   *  \tparam     T       Type of object (should match symbol type in *.cu file)
   *  \param[in]  a_idxLib
   *                      Index of library in which to find symbol (e.g.,
   *                      c_BaseTools_lib or c_Application_lib
   *  \param[in]  a_idxSymb
   *                      Enum giving handle to the symbol
   *  \param[in]  a_val   Object on the host
   */
  template <typename S, typename T>
  void copySymbolToDevice(const int a_idxLib,
                          const S   a_idxSymb,
                          const T&  a_val) noexcept
    {
      CH_assert(testInitCompFlag(InitComponentFlag::module0));
      getSymbol(a_idxLib, a_idxSymb).copyHtoD(a_val);
    }

  /// Asynchronous copy from device to host
  /** \tparam     S       Type of enum for symbol
   *  \tparam     T       Type of object (should match symbol type in *.cu file)
   *  \param[in]  a_idxLib
   *                      Index of library in which to find symbol (e.g.,
   *                      c_BaseTools_lib or c_Application_lib
   *  \param[in]  a_idxSymb
   *                      Enum giving handle to the symbol
   *  \param[in]  a_val   Object on the host
   *  \param[in]  a_stream
   *                      Stream identifier
   */
  template <typename S, typename T>
  void copySymbolToDeviceAsync(const int a_idxLib,
                               const S   a_idxSymb,
                               const T&  a_val,
                               CUstream  a_stream = c_defaultStream) noexcept
    {
      CH_assert(testInitCompFlag(InitComponentFlag::module0));
      getSymbol(a_idxLib, a_idxSymb).copyHtoDAsync(a_val, a_stream);
    }

private:

  /// Set a cuda component as being defined
  void setInitCompFlag(const InitComponentFlag a_icf);

  /// Set a cuda component as undefined
  void unsetInitCompFlag(const InitComponentFlag a_icf);

  /// Get the device handle
  void getDevice(const int iDevice);

  /// Create a context
  int createContext(const CUfunc_cache   a_cachePreference,
                    const CUsharedconfig a_sharedMemBankSize);

  /// Use the primary context (allows use of runtime API)
  void retainPrimaryContext(const CUfunc_cache   a_cachePreference,
                            const CUsharedconfig a_sharedMemBankSize);

  /// Append a slash to a string if it's not the last character
  void appendSlash(std::string& a_str) const;

  /// Find and load an embedded module returning the base pointer
  int baseLoadEmbeddedModule(CUmodule*         a_module,
                             const char *const a_moduleName,
                             const void*       a_embeddedBin);

  /// Find and load an external module (GPU executable) returning base pointer
  void baseLoadModule(CUmodule*         a_module,
                      const char *const a_moduleName,
                      const char *const a_moduleDir,
                      const char *const a_searchPath);

//--Data members

public:

  CUdevice m_device;                  ///< Handle to the device

  // Primary context and module
  CUcontext m_context0;               ///< Handle to context 0 (this is the
                                      ///< the primary context for the device
                                      ///< and retrieved using
                                      ///< cuDevicePrimaryCtxRetain)
  CUmodule m_module0;                 ///< Handle to module 0 (for Chombo and
                                      ///< application functions and symbols
  std::string m_module0Name;          ///< Name of primary module
  stc::Vector<std::vector<DeviceFunctionConfig>, c_numModule0Lib>
  m_functionConfigs;                  ///< Global function configurations for
                                      ///< primary module
  stc::Vector<std::vector<DeviceSymbolConfig>, c_numModule0Lib> m_symbolConfigs;
                                      ///< Symbol configurations for primary
                                      ///< module

  // External or additional contexts and modules
  std::vector<CUcontext> m_contextVec;
                                      ///< Handles to the contexts
  std::vector<DriverModule*> m_moduleVec;
                                      ///< External module data

private:

  unsigned m_initCompFlag;            ///< Flags for which CUDA components have
                                      ///< been defined.
};

}  // End of namespace Cuda

#endif  /* ! defined _CUDADRIVER_H_ */
