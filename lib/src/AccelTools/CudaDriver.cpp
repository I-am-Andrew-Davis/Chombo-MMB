#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifdef CH_GPU


/******************************************************************************/
/**
 * \file
 *
 * \brief Device initialization and support for using the Cuda driver API
 *
 *//*+*************************************************************************/

#include "CH_assert.H"
#include "CH_System.H"
#include "CudaDriver.H"


/*******************************************************************************
 *
 * Definition of global device handle
 * Note: call initialize() before using and finalize() when finished.
 *
 ******************************************************************************/

CH_Cuda::Driver CH_Cuda::g_device;


/*******************************************************************************
 *
 * Class CudaDriver: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Default constructor initializes all cuda components as undefined
/*--------------------------------------------------------------------*/

CH_Cuda::Driver::Driver()
  : m_initCompFlag(0)
{ }

/*--------------------------------------------------------------------*/
//  Destructor undefines defined cuda components
/*--------------------------------------------------------------------*/

CH_Cuda::Driver::~Driver()
{
  finalize();
}


/*==============================================================================
 * Member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Test if a cuda component has been defined by CudaDriver
/** \param[in]  a_icf   Component to test
 *//*-----------------------------------------------------------------*/

bool
CH_Cuda::Driver::testInitCompFlag(const InitComponentFlag a_icf) const
{
  return (m_initCompFlag & (unsigned)a_icf);
}

/*--------------------------------------------------------------------*/
//  Initialize the device
/** \param[in]  a_verbose
 *                      0 - Only print errors (default)
 *                      1 - Some device data
 *                      2 - All device data
 *  \param[in]  a_cachePreference
 *                      Cache preferences for the context (these can
 *                      be adjusted for individual kernels)
 *                      CU_FUNC_CACHE_PREFER_NONE: no preference for
 *                        shared memory or L1
 *                      CU_FUNC_CACHE_PREFER_SHARED: prefer larger
 *                        shared memory and smaller L1 cache (default)
 *                      CU_FUNC_CACHE_PREFER_L1: prefer larger L1 cache
 *                        and smaller shared memory
 *                      CU_FUNC_CACHE_PREFER_EQUAL: prefer equal sized
 *                        L1 cache and shared memory
 *  \param[in]  a_sharedMemBankSize
 *                      Preferences for shared memory bank size
 *                      CU_SHARED_MEM_CONFIG_DEFAULT_BANK_SIZE: set 
 *                        bank width to the default initial setting
 *                        (currently, four bytes).  (default)
 *                      CU_SHARED_MEM_CONFIG_FOUR_BYTE_BANK_SIZE: set
 *                        shared memory bank width to be natively four
 *                        bytes.
 *                      CU_SHARED_MEM_CONFIG_EIGHT_BYTE_BANK_SIZE: set
 *                        shared memory bank width to be natively
 *                        eight bytes.
 *//*-----------------------------------------------------------------*/

void
CH_Cuda::Driver::initialize(const int            a_verbose,
                            const CUfunc_cache   a_cachePreference,
                            const CUsharedconfig a_sharedMemBankSize)
{

//--Initialize

  if (a_verbose > 0)
    {
      WRITE_MSG("Initializing CUDA GPU device.\n", msgHeader1);
    }
  checkCudaErrors(cuInit(0));

//--Find a device

  int numDevice;
  checkCudaErrors(cuDeviceGetCount(&numDevice));
  if (numDevice == 0)
    {
      MayDay::Error("No NVIDIA GPU device found");
    }

  int iDevice;
  int ccminor;
  int ccmajor;
  for (iDevice = 0; iDevice != numDevice; ++iDevice)
    {
      CUdevice device;
      checkCudaErrors(cuDeviceGet(&device, iDevice));
      checkCudaErrors(
        cuDeviceGetAttribute(&ccmajor,
                             CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR,
                             device));
      if (ccmajor >= 1) break;
    }
  if (iDevice == numDevice)
    {
      MayDay::Error("No GPU device found that supports CUDA");
    }
  getDevice(iDevice);

//--Report device name

   if (a_verbose > 0)
     {
       char deviceName[128];
       checkCudaErrors(cuDeviceGetName(deviceName, 128, m_device));
       std::sprintf(msgBuffer, "%-*s %s", msgLabelWidth, "Device name:",
                    deviceName);
       WRITE_MSG(msgBuffer, msgHeader2);
       checkCudaErrors(
         cuDeviceGetAttribute(&ccminor,
                              CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR,
                              m_device));
       std::sprintf(msgBuffer, "%-*s %d.%d", msgLabelWidth,
                    "Compute capability:", ccmajor, ccminor);
       WRITE_MSG(msgBuffer, msgHeader2);
     }

//--Report some attributes

   if (a_verbose > 1)
     {
       size_t globalMemByte;
       checkCudaErrors(cuDeviceTotalMem(&globalMemByte, m_device));
       float globalMem = globalMemByte;
       const char *memUnit = reduceUnitCount(0, globalMem);
       std::sprintf(msgBuffer, "%-*s %g %s", msgLabelWidth, "Global memory:",
                    globalMem, memUnit);
       WRITE_MSG(msgBuffer, msgHeader2);

       int clockRateKHz;
       checkCudaErrors(cuDeviceGetAttribute(&clockRateKHz,
                                         CU_DEVICE_ATTRIBUTE_CLOCK_RATE,
                                         m_device));
       float clockRate = 1000*clockRateKHz;
       memUnit = reduceUnitCount(1, clockRate);
       std::sprintf(msgBuffer, "%-*s %g %s\n", msgLabelWidth, "Clock rate:",
                    clockRate, memUnit);
       WRITE_MSG(msgBuffer, msgHeader2);
     }

//--Create a context

   /*
     To use the primary context, something must have used the runtime to
     "lazily" create it.  But we are probably created the first and only
     context so cannot use this
   */
   // retainPrimaryContext(a_cachePreference, a_sharedMemBankSize);
   createContext(a_cachePreference, a_sharedMemBankSize);
}

/*--------------------------------------------------------------------*/
//  Unload all external modules
/** Finalize will also accomplish this task.  Leaves primary module 0
 *  unchanged.
 *//*-----------------------------------------------------------------*/

void
CH_Cuda::Driver::unloadModules()
{
  for (auto ptrMdl : m_moduleVec)
    {
      checkCudaErrors(cuModuleUnload(ptrMdl->m_module));
      delete ptrMdl;
    }
  m_moduleVec.clear();
}

/*--------------------------------------------------------------------*/
//  Unload primary module 0
/** Finalize will also accomplish this task.  Leaves external modules
 *  unchanged.
 *//*-----------------------------------------------------------------*/

void
CH_Cuda::Driver::unloadModule0()
{
  if (testInitCompFlag(InitComponentFlag::module0))
    {
      checkCudaErrors(cuModuleUnload(m_module0));
      unsetInitCompFlag(InitComponentFlag::module0);
    }
}

/*--------------------------------------------------------------------*/
//  Fnalize the device
/** Same as destructor but required for object in global memory space
 *//*-----------------------------------------------------------------*/

void
CH_Cuda::Driver::finalize()
{
  // Primary module
  unloadModule0();
  // External modules
  unloadModules();
  // Primary context
  if (testInitCompFlag(InitComponentFlag::context0))
    {
      checkCudaErrors(cuDevicePrimaryCtxRelease(m_device));
    }
  // External contexts
  if (testInitCompFlag(InitComponentFlag::context))
    {
      for (auto ctx : m_contextVec)
        {
          checkCudaErrors(cuCtxDestroy(ctx));
        }
      m_contextVec.clear();
    }
  m_initCompFlag = 0;
}

/*--------------------------------------------------------------------*/
//  Load an application module
/** This is the primary means for loading functions and symbols from
 *  both the Chombo libraries and the application.  The standard way
 *  is to include a *_CUX.H file which contains the byte string for an
 *  embedded CUDA binary.  Alternatively, the module may be in a
 *  separate file.
 *  \param[in]  a_moduleName
 *                      Name of the module.  If the module is loaded
 *                      from a file, this should be the name of the
 *                      file.  If the module is embedded, the name is
 *                      only used as a reference and can be anything.
 *                      It is suggested to use a_moduleName =
 *                      "a_embeddedBin"
 *  \param[in]  a_embeddedBin
 *                      Byte string for an embedded CUDA binary.  In
 *                      the Chombo build system, this is found in the
 *                      *_CUX.H files and has the same name as *.  It
 *                      is suggested to use a module name of "*" for
 *                      parameter a_moduleName.
 *  \param[in]  a_moduleDir
 *                      A directory (to append to a_searchPath) in
 *                      which the module might be found.  You can
 *                      specify the directory in a_searchPath and
 *                      leave this NULL -- it is only here to provide
 *                      more flexibility and for use with
 *                      'path_to_current_process'
 *  \param[in]  a_searchPath
 *                      Path to search for the file.
 *//*-----------------------------------------------------------------*/
void
CH_Cuda::Driver::loadApplicationModule(const char *const a_moduleName,
                                       const void*       a_embeddedBin,
                                       const char *const a_moduleDir,
                                       const char *const a_searchPath)
{
  CH_assert(a_moduleName != nullptr);
  CH_assert(!testInitCompFlag(InitComponentFlag::module0));
  m_module0Name = a_moduleName;
  int moduleStatus = 1;
  if (a_embeddedBin != nullptr)
    {
      moduleStatus = baseLoadEmbeddedModule(&m_module0,
                                            a_moduleName,
                                            a_embeddedBin);
    }
  if (moduleStatus)  // 0 = previous success
    {
      baseLoadModule(&m_module0, a_moduleName, a_moduleDir, a_searchPath);
    }
  setInitCompFlag(InitComponentFlag::module0);

//--Load the functions and symbols

  for (int idxLib = 0; idxLib != c_numModule0Lib; ++idxLib)
    {
      // Functions
      for (auto& funcConfig : m_functionConfigs[idxLib])
        {
          checkCudaErrors(cuModuleGetFunction(&funcConfig.m_CUfunction,
                                              m_module0,
                                              funcConfig.m_name.c_str()));
          checkCudaErrors(cuFuncSetCacheConfig(funcConfig.m_CUfunction,
                                               funcConfig.m_cacheConfig));
          checkCudaErrors(cuFuncSetSharedMemConfig(
                            funcConfig.m_CUfunction,
                            funcConfig.m_sharedMemBankSize));
        }
      // Symbols
      for (auto& symbConfig : m_symbolConfigs[idxLib])
        {
          checkCudaErrors(cuModuleGetGlobal(&symbConfig.m_CUpointer,
                                            &symbConfig.m_numBytes,
                                            m_module0,
                                            symbConfig.m_name.c_str()));
        }
    }
}

/*--------------------------------------------------------------------*/
//  Find and load an embedded module returning the base pointer
/** Simply search the given bytes which should be an embedded CUDA
 *  cubin or fatbin
 *  \param[out] a_module
 *                      Handle for module
 *  \param[in]  a_moduleName
 *                      Name of the module
 *  \param[in]  a_embeddedBin
 *                      Byte string for an embedded CUDA binary
 *  \return             0 - success
 *                      1 - module not found
 *//*-----------------------------------------------------------------*/

int
CH_Cuda::Driver::baseLoadEmbeddedModule(CUmodule*         a_module,
                                        const char *const a_moduleName,
                                        const void*       a_embeddedBin)
{
  CH_assert(testInitCompFlag(InitComponentFlag::anyContext));
  CH_assert(a_moduleName != NULL);
#ifdef CH_FATBIN
  CUresult err = cuModuleLoadFatBinary(a_module, a_embeddedBin);
#else
  CUresult err = cuModuleLoadData(a_module, a_embeddedBin);
#endif
  if (CUDA_SUCCESS != err)
    {
      std::sprintf(CH_Cuda::msgBuffer, "Searching for embedded module: Cuda "
                   "driver error %d", err);
      MayDay::Warning(CH_Cuda::msgBuffer);
    }
  return (int)err;
}

/*--------------------------------------------------------------------*/
//  Find and load an external module (GPU executable) returning the
//  base pointer
/** Mostly methods for finding the module.  We search:
 *  <ol>
 *    <li> a_searchPath/a_moduleName
 *    <li> a_searchPath/a_moduleDir/a_moduleName
 *    <li> a_moduleName
 *    <li> a_moduleDir/a_moduleName
 *    <li> 'path_to_current_process'/a_moduleName
 *    <li> 'path_to_current_process'/a_moduleDir/a_moduleName
 *    <li> 'path_to_current_process'/cubin/
 *           {match/executable_name/.*([1-6]d\..*)\.ex/}/a_moduleName
 *  </ol>
 *  The last option above is the default for modules built with the
 *  Chombo make system.
 *  \param[out] a_module
 *                      Handle for module
 *  \param[in]  a_moduleName
 *                      Name of the module
 *  \param[in]  a_moduleDir
 *                      A directory (to append to a_searchPath) in
 *                      which the module might be found.  You can
 *                      specify the directory in a_searchPath and
 *                      leave this NULL -- it is only here to provide
 *                      more flexibility and for use with
 *                      'path_to_current_process'
 *  \param[in]  a_searchPath
 *                      Path to search for the file.
 *  \note
 *  <ul>
 *    <li> For arguments left NULL, paths involving those arguments
 *         are not searched.
 *    <li> 'path_to_current_process' is determined by system calls
 *         specified in System::getProcessPath().
 *  </ul>
 *//*-----------------------------------------------------------------*/

void
CH_Cuda::Driver::baseLoadModule(CUmodule*         a_module,
                                const char *const a_moduleName,
                                const char *const a_moduleDir,
                                const char *const a_searchPath)
{
  CH_assert(testInitCompFlag(InitComponentFlag::anyContext));
  CH_assert(a_moduleName != NULL);

  std::string modulePath;

//--In a_searchPath

  if (a_searchPath)
    {
      modulePath = a_searchPath;
      appendSlash(modulePath);
      modulePath += a_moduleName;
      if (CH_System::fileExists(modulePath.c_str()))
        {
          goto found_module;
        }
      if (a_moduleDir)
        {
          modulePath = a_searchPath;
          appendSlash(modulePath);
          modulePath += a_moduleDir;
          appendSlash(modulePath);
          modulePath += a_moduleName;
          if (CH_System::fileExists(modulePath.c_str()))
            {
              goto found_module;
            }
        }
    }

//--In ./

  modulePath = a_moduleName;
  if (CH_System::fileExists(modulePath.c_str()))
    {
      goto found_module;
    }
  if (a_moduleDir)
    {
      modulePath = a_moduleDir;
      appendSlash(modulePath);
      modulePath += a_moduleName;
      if (CH_System::fileExists(modulePath.c_str()))
        {
          goto found_module;
        }
    }

//--In path of executable

  {
    // Get the path to the executable
    int pathSize = 128;  // Start with length 256 (will be x2)
    int ier = -pathSize;
    while (ier == -pathSize && pathSize <= 4096)
      {
        pathSize *= 2;
        char *exePath = new char[pathSize];
        ier = CH_System::getProcessPath(exePath, pathSize);
        if (ier > 0)
          {
            modulePath = exePath;
          }
        delete[] exePath;
      }
    if (ier > 0)
      {
        // Strip the executable name
        size_t spos = modulePath.find_last_of('/');
        std::string execName;
        if (spos == std::string::npos)
          {
            modulePath = ".";
          }
        else
          {
            execName = modulePath.substr(spos+1);
            modulePath.erase(spos);
          }
        // Search for the module
        appendSlash(modulePath);
        std::string saveModulePath(modulePath);
        modulePath += a_moduleName;
        if (CH_System::fileExists(modulePath.c_str()))
          {
            goto found_module;
          }
        if (a_moduleDir)
          {
            modulePath = saveModulePath;
            modulePath += a_moduleDir;
            appendSlash(modulePath);
            modulePath += a_moduleName;
            if (CH_System::fileExists(modulePath.c_str()))
              {
                goto found_module;
              }
          }
        if (execName.size() > 0)
          {
            spos = execName.find_first_of('.');
            size_t len = execName.size();
            if (execName[--len] == 'x' &&
                execName[--len] == 'e' &&
                execName[--len] == '.' &&
                execName[--spos] == 'd')
              // This is a Chombo executable
              {
                --spos;
                --len;
                // Get the config string
                std::string config(execName.substr(spos, len-spos+1));
                modulePath = saveModulePath;
                modulePath += "cubin/";
                modulePath += config;
                modulePath += '/';
                modulePath += a_moduleName;
                if (CH_System::fileExists(modulePath.c_str()))
                  {
                    goto found_module;
                  }
              }
          }
      }
    else if (pathSize > 4096)
      {
        MayDay::Error("Size of path of executable is > 4096 characters?!?");
      }
  }

//--Could not find module

  std::sprintf(msgBuffer, "Could not find module %s.", a_moduleName);
  MayDay::Error(msgBuffer);

//--Found module

  found_module: ;
  std::sprintf(msgBuffer, "Found module at %s.\n", modulePath.c_str());
  WRITE_MSG(msgBuffer, msgHeader1);
  checkCudaErrors(cuModuleLoad(a_module, modulePath.c_str()));
}

/*--------------------------------------------------------------------*/
//  Set a cuda component as being defined
/** \param[in]  a_icf   Component
 *//*-----------------------------------------------------------------*/

void
CH_Cuda::Driver::setInitCompFlag(const InitComponentFlag a_icf)
{
  m_initCompFlag += (unsigned)a_icf;
}

/*--------------------------------------------------------------------*/
//  Set a cuda component as undefined
/** \param[in]  a_icf   Component
 *//*-----------------------------------------------------------------*/

void
CH_Cuda::Driver::unsetInitCompFlag(const InitComponentFlag a_icf)
{
  if (testInitCompFlag(a_icf))
    {
      m_initCompFlag -= (unsigned)a_icf;
    }
}

/*--------------------------------------------------------------------*/
//  Get the device handle
/** \param[in]  iDevice Device number to get handle for
 *//*-----------------------------------------------------------------*/

void
CH_Cuda::Driver::getDevice(const int iDevice)
{
  checkCudaErrors(cuDeviceGet(&m_device, iDevice));
  setInitCompFlag(InitComponentFlag::device);
}

/*--------------------------------------------------------------------*/
//  Create a context
/** This will create a new context (like a process) for the GPU.  It
 *  most cases, it is recommended to use the primary context.
 *  \param[in]  a_cachePreference
 *                      Cache preferences for the context (these can
 *                      be adjusted for individual kernels, see
 *                      initialize() for full description)
 *  \param[in]  a_sharedMemBankSize
 *                      Preferences for shared memory bank size (see
 *                      initialize() for full description)
 *  \return             The index of this context in the vector
 *//*-----------------------------------------------------------------*/

int
CH_Cuda::Driver::createContext(const CUfunc_cache   a_cachePreference,
                               const CUsharedconfig a_sharedMemBankSize)
{
  CH_assert(testInitCompFlag(InitComponentFlag::device));
  m_contextVec.push_back(CUcontext{});
  checkCudaErrors(cuCtxCreate(&m_contextVec.back(),
                              CU_CTX_SCHED_BLOCKING_SYNC | CU_CTX_MAP_HOST,
                              m_device));
  setInitCompFlag(InitComponentFlag::context);
  checkCudaErrors(cuCtxSetCacheConfig(a_cachePreference));
  checkCudaErrors(cuCtxSetSharedMemConfig(a_sharedMemBankSize));
  return m_contextVec.size() - 1;
}

/*--------------------------------------------------------------------*/
//  Use the primary context (allows use of runtime API)
/** To use this, the context must have been lazily created by a
 *  library or method using the runtime API.  Then we could reuse that
 *  context.
 *  \param[in]  a_cachePreference
 *                      Cache preferences for the context (these can
 *                      be adjusted for individual kernels, see
 *                      initialize() for full description)
 *  \param[in]  a_sharedMemBankSize
 *                      Preferences for shared memory bank size (see
 *                      initialize() for full description)
 *//*-----------------------------------------------------------------*/

void
CH_Cuda::Driver::retainPrimaryContext(const CUfunc_cache   a_cachePreference,
                                      const CUsharedconfig a_sharedMemBankSize)
{
  CH_assert(testInitCompFlag(InitComponentFlag::device));
  checkCudaErrors(cuDevicePrimaryCtxRetain(&m_context0, m_device));
  setInitCompFlag(InitComponentFlag::context0);
  checkCudaErrors(cuCtxSetCacheConfig(a_cachePreference));
  checkCudaErrors(cuCtxSetSharedMemConfig(a_sharedMemBankSize));
}

/*--------------------------------------------------------------------*/
//  Append a slash to a string if it's not the last character
/** \param[in]  str     String to modify
 *  \param[out] str     Slash as last character
 *//*-----------------------------------------------------------------*/

void
CH_Cuda::Driver::appendSlash(std::string& a_str) const
{
  if (a_str[a_str.length() - 1] != '/')
    {
      a_str += '/';
    }
}

#endif  /* CH_GPU */
