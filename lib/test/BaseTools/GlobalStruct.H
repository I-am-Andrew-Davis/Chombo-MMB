#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _GLOBALSTRUCT_H_
#define _GLOBALSTRUCT_H_

#include "CH_HashTable.H"
#ifdef CH_GPU
  #include "CudaDriver.H"
#endif

/*
  Notes:
    - tag __CUDACC__ means the file is process with nvcc.  The class is
      modified whether on the cpu or gpu.
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

/*
  For a global object with hash tables (as is explored here), it is considered
  best to only define once (so both tables must be defined).

  For a non-global object, it is probably better to have the hash table as a
  template argument and build a new static object from a dynamic object.
*/


/// Handles to module global functions for this application
enum class App_Global
{
  testGlobalStructGet = 0
};

/// Handles to module symbols for this application
enum class App_Symbol
{
  g_verbose = 0,
  g_stat,
  g_table
};

class GlobalStruct;

/*
  Namespaces are bad idea for anything that you want to be able to refer to
  in a non-mangled sense.  On the GPU we do need to find these symbols by name.

  Note: this could be made with C linkage using extern "C" here and extern "C" {
  CH_Hash::StaticTable<int, ::GlobalStruct>* g_GlobalStruct_table; } for the
  definition but is not required since it is just a pointer
*/
extern DEVICE CH_Hash::StaticTable<int, GlobalStruct>* g_GlobalStruct_table;
#ifndef __CUDACC__
  extern CH_Hash::DynamicTable<int, GlobalStruct>* g_GlobalStruct_builderTable;
#endif

class GlobalStruct
{

public:

  GlobalStruct()
    { }

#ifndef __CUDACC__
  static void createInstance();
  static void destroyInstance();
#endif

  /*
    This is still HOSTDEVICE but nvcc complains about using
    g_GlobalStruct_table from the host if we keep the host tag (and DEVICE is
    set to null when compiling for the CPU)
  */
  DEVICE static GlobalStruct& get(const int a_key)
    {
#ifndef __CUDACC__
      if (g_GlobalStruct_builderTable)
        {
          auto iter = g_GlobalStruct_builderTable->find(a_key);
          if (iter == g_GlobalStruct_builderTable->end())
            {
              auto ins = g_GlobalStruct_builderTable->insert(
                std::make_pair(a_key, GlobalStruct{}));
              iter = ins.first;
            }
          return iter->second;
        }
      else
#endif
        {
          auto iter = g_GlobalStruct_table->find(a_key);
          CH_assert(iter != g_GlobalStruct_table->end());
          return iter->second;
        }
    }

#ifndef __CUDACC__
  static void finalize()
    {
      // Change from using a dynamaic table to a static table.
      g_GlobalStruct_table =
        new CH_Hash::StaticTable<int, GlobalStruct>(
          *g_GlobalStruct_builderTable, -1);
      delete g_GlobalStruct_builderTable;
      g_GlobalStruct_builderTable = nullptr;
#ifdef CH_GPU
      // Copy table arrays to the device
      g_GlobalStruct_table->copyToDeviceAsync();
      // This copy adjusts pointers to memory on the device.
      CH_Hash::StaticTable<int, GlobalStruct> h_table(
        CH_Cuda::copy_to_device{}, *g_GlobalStruct_table);
      // New get that object on the device
      const size_t numBytes = sizeof(CH_Hash::StaticTable<int, GlobalStruct>);
      CH_Cuda::SymbolPair<CH_Hash::StaticTable<int, GlobalStruct>> pointers;
      pointers.host = &h_table;
      checkCudaErrors(cuMemAlloc(&pointers.device, numBytes));
      pointers.copyHtoD(1);  // Synchronous since h_table is local
      // Save the pointer to the object in a particular symbol in global memory
      // on the device
      CH_Cuda::g_device.copySymbolToDeviceAsync(CH_Cuda::Application_lib,
                                                App_Symbol::g_table,
                                                pointers.device);
#endif
    }
#endif

  double a_val;
};

#endif  /* _GLOBALSTRUCT_H_ */
