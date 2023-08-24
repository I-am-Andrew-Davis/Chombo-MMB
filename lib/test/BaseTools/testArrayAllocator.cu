#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif


/******************************************************************************/
/**
 * \file
 *
 * \brief Kernels for testArrayAllocator on the GPU
 *
 *//*+*************************************************************************/

#include <cstdint>

#include "ArrayAllocator.H"
#include "Misc.H"

__device__ bool g_verbose;

// Tracking of errors on the device
__device__ unsigned g_stat;
__device__ inline void
incrStat()
{
  atomicInc(&g_stat, (unsigned)-1);
}

// Record of device pointer from the CPU (to make sure the same ones are used
// here and that we are not simply accessing host memory)
__device__ uintptr_t g_devicePtr;

/// Test basic reading and writing to an argument array
extern "C"
__global__ void
arrayTestRW(Array_impl<Real,
                       CUDAArrayAlloc<Real, ArrayClassIndex::BaseFab>> a_arrA)
{
  // Check pointers
  if (threadIdx.x == 0)
    {
      if (g_verbose)
        {
          printf("Pointer location on GPU: %p\n", a_arrA.data());
        }
      if (reinterpret_cast<uintptr_t>(a_arrA.data()) != g_devicePtr) incrStat();
    }
  // The argument array aliases memory controlled by the array on the CPU
  if (a_arrA.allocBy() != AllocBy::alias) incrStat();
  if (a_arrA.size() != 4) incrStat();
  if (a_arrA[threadIdx.x] != threadIdx.x + 0.5) incrStat();
  a_arrA[threadIdx.x] += 0.25;
}

/// Test aliases and raw memory
extern "C"
__global__ void
arrayTestAliasRaw()
{
  __shared__ Real data[16];
  UNINITIALIZED__shared__(
    Array_impl<Real, CUDAArrayAlloc<Real, ArrayClassIndex::BaseFab>> raw);
  if (threadIdx.x == 0)
    {
      raw.INITIALIZEshared();
      raw.defineRaw(data, 16, 2.5);
      // alias.defineAlias(data, 16);
    }
  __syncthreads();
  Array_impl<Real, CUDAArrayAlloc<Real, ArrayClassIndex::BaseFab>> alias;
  alias.defineAlias(data, 16);
  alias[threadIdx.x] += 1.5;
  if (raw[threadIdx.x] != 4.0) incrStat();
}

/// Test device allocations
extern "C"
__global__ void
arrayTestAlloc()
{
  UNINITIALIZED__shared__(
    Array_impl<Real, CUDAArrayAlloc<Real, ArrayClassIndex::BaseFab>> arr);
  if (threadIdx.x == 0)
    {
      arr.INITIALIZEshared();
      arr.define(blockDim.x, 3.2);
    }
  __syncthreads();
  arr[threadIdx.x] += threadIdx.x + 0.2;
  if (Misc::compare(arr[threadIdx.x], threadIdx.x + 3.4, 6))
    {
      printf("Thread: %d\n", threadIdx.x);
      printf("Val: %lf %lf\n", arr[threadIdx.x], threadIdx.x + 3.4);
      incrStat();
    }
}
