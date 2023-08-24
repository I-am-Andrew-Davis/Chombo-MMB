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
 * \brief Kernels for testHashTable on the GPU
 *
 *//*+*************************************************************************/

#include <cstdint>

#include "DynArray.H"
#include "CH_Cuda_shared_ptr.H"

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
__device__ uintptr_t g_shrDevicePtr;
__device__ uintptr_t g_arrDevicePtr;

/// Test basic reading and writing to an argument array
extern "C"
__global__ void
arrayTestRW(CH_Cuda::shared_ptr<dyn::Array<Real, 3>> a_arrA,
            CH_Cuda::shared_ptr<dyn::Array<Real, 3>> a_arrB)
{
  // Check pointers
  if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0)
    {
      if (g_verbose)
        {
          printf("A: Shared pointer location on GPU: %p\n", a_arrA.get());
          printf("A: Array data location on GPU    : %p\n", a_arrA->data());
          printf("B: Shared pointer location on GPU: %p\n", a_arrB.get());
          printf("B: Array data location on GPU    : %p\n", a_arrB->data());
        }
      if (reinterpret_cast<uintptr_t>(a_arrA.get()) != g_shrDevicePtr)
        incrStat();
      if (reinterpret_cast<uintptr_t>(a_arrA->data()) != g_arrDevicePtr)
        incrStat();
      if (reinterpret_cast<uintptr_t>(a_arrB.get()) != g_shrDevicePtr)
        incrStat();
      if (reinterpret_cast<uintptr_t>(a_arrB->data()) != g_arrDevicePtr)
        incrStat();
    }
  // Test read and then perform write
  if ((*a_arrA)(threadIdx.x - 1, threadIdx.y, threadIdx.z + 1) !=
      threadIdx.x - 1.0 + threadIdx.y + threadIdx.z + 1.0) incrStat();
  // Writes are performed on B (tests on host will use A)
  (*a_arrB)[threadIdx.z + 1][threadIdx.y][threadIdx.x - 1] += 0.5;
}
