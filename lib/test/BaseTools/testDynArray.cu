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
arrayTestRW(dyn::Array<Real, 3> a_arrA)
{
  // Check pointers
  if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0)
    {
      if (g_verbose)
        {
          printf("Pointer location on GPU: %p\n", a_arrA.data());
        }
      if (reinterpret_cast<uintptr_t>(a_arrA.data()) != g_devicePtr) incrStat();
    }
  // Test read and then perform write
  if (a_arrA(threadIdx.x - 1, threadIdx.y, threadIdx.z + 1) !=
      threadIdx.x - 1.0 + threadIdx.y + threadIdx.z + 1.0) incrStat();
  a_arrA[threadIdx.z + 1][threadIdx.y][threadIdx.x - 1] += 0.5;
}

/// Test aliasing an argument array (aliased array per thread)
extern "C"
__global__ void
arrayTestAliasThread(stc::Vector<int, 3> a_lb,
                     dyn::Array<Real, 3> a_arrA)
{
  dyn::Array<Real, 3> arrB(dyn::alias_format{}, a_arrA);
  if (threadIdx.z == 0)
    {
      arrB(threadIdx.x + a_lb[0],
           threadIdx.y + a_lb[1],
           threadIdx.z + a_lb[2]) = 2.5;
    }

  // Another alias
  dyn::Array<Real, 3> arrC(dyn::alias_format{}, a_arrA);
  if (threadIdx.z == 1)
    {
      arrC(threadIdx.x + a_lb[0],
           threadIdx.y + a_lb[1],
           threadIdx.z + a_lb[2]) = 3.5;
    }

  // Another alias (we directly address memory and the LB on the alias is
  // therefore 0).
  arrC.defineAlias(a_arrA.data() + 2*a_arrA.stride(2),
                   a_arrA.size(0), a_arrA.size(1), 1);
  if (threadIdx.z == 2)
    {
      arrC(threadIdx.x, threadIdx.y, 0) = 4.5;
    }
}

/// Test aliasing an argument array (aliased array in shared memory)
/** Check the ptx output (OPT=OPTHIGH) and you will likely find that it takes
 *  more registers to keep the dyn::Array object in shared memory.  It is
 *  usually best to store the dyn::Array object per thread.  The previous
 *  kernel is preferred.  It likely has higher performance, better memory use,
 *  and is cleaner code.
 */
extern "C"
__global__ void
arrayTestAliasShared(stc::Vector<int, 3> a_lb,
                     dyn::Array<Real, 3> a_arrA)
{
  UNINITIALIZED__shared__(dyn::Array<Real, 3> arrB);
  if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0)
    {
      arrB.INITIALIZEshared();
      arrB.defineAlias(a_arrA);
    }
  __syncthreads();
  if (threadIdx.z == 0)
    {
      arrB(threadIdx.x + a_lb[0],
           threadIdx.y + a_lb[1],
           threadIdx.z + a_lb[2]) = 2.5;
    }

  // Another alias
  UNINITIALIZED__shared__(dyn::Array<Real, 3> arrC);
  if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0)
    {
      arrC.INITIALIZEshared();
      arrC.defineAlias(a_arrA);
    }
  __syncthreads();
  if (threadIdx.z == 1)
    {
      arrC(threadIdx.x + a_lb[0],
           threadIdx.y + a_lb[1],
           threadIdx.z + a_lb[2]) = 3.5;
    }

  // Another alias in shared memory (we directly address memory and the LB
  // on the alias is therefore 0).
  __syncthreads();
  if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0)
    {
      arrC.defineAlias(a_arrA.data() + 2*a_arrA.stride(2), a_arrA.size(0),
                       a_arrA.size(1), 1);
    }
  __syncthreads();
  if (threadIdx.z == 2)
    {
      arrC(threadIdx.x, threadIdx.y, 0) = 4.5;
    }
}

/// Test using raw memory (array stored per thread)
extern "C"
__global__ void
arrayTestRawThread()
{
  __shared__ Real data[64];
  dyn::Array<Real, 3> arrA(dyn::raw_format{}, data, 4, 4, 4);
  if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0)
    {
      arrA = 2.5;
    }
  __syncthreads();
  if (arrA(threadIdx.x, threadIdx.y, threadIdx.z) != 2.5) incrStat();
  __syncthreads();
  dyn::Array<Real, 3> arrB(dyn::raw_format{}, data, 2, 2, 2);
  if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0)
    {
      arrB = 3.5;
    }
  __syncthreads();
  if (threadIdx.x < 2 && threadIdx.y < 2 && threadIdx.z < 2)
    {
      if (arrB(threadIdx.x, threadIdx.y, threadIdx.z) != 3.5) incrStat();
    }
  Real val = 2.5;
  if (threadIdx.y < 2 && threadIdx.z == 0)  // First 8 elements in data
    {
      val = 3.5;
    }
  if (arrA(threadIdx.x, threadIdx.y, threadIdx.z) != val) incrStat();
}

/// Test using raw memory (array in shared memory)
/** Check the ptx output (OPT=OPTHIGH) and you will likely find that it takes
 *  more registers to keep the dyn::Array object in shared memory.  It is
 *  usually best to store the dyn::Array object per thread.  The previous
 *  kernel is preferred.  It likely has higher performance, better memory use,
 *  and is cleaner code.
 */
extern "C"
__global__ void
arrayTestRawShared()
{
  __shared__ Real data[64];
  UNINITIALIZED__shared__(dyn::Array<Real, 3> arrA);
  if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0)
    {
      arrA.INITIALIZEshared();
      arrA.defineRaw(data, 4, 4, 4);
      arrA = 2.5;
    }
  __syncthreads();
  if (arrA(threadIdx.x, threadIdx.y, threadIdx.z) != 2.5) incrStat();
  __syncthreads();
  UNINITIALIZED__shared__(dyn::Array<Real, 3> arrB);
  if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0)
    {
      arrB.INITIALIZEshared();
      arrB.defineRaw(data, 2, 2, 2);
      arrB = 3.5;
    }
  __syncthreads();
  if (threadIdx.x < 2 && threadIdx.y < 2 && threadIdx.z < 2)
    {
      if (arrB(threadIdx.x, threadIdx.y, threadIdx.z) != 3.5) incrStat();
    }
  Real val = 2.5;
  if (threadIdx.y < 2 && threadIdx.z == 0)  // First 8 elements in data
    {
      val = 3.5;
    }
  if (arrA(threadIdx.x, threadIdx.y, threadIdx.z) != val) incrStat();
}

/// Test device allocations
extern "C"
__global__ void
arrayTestAlloc()
{
  if (threadIdx.x == 0 && threadIdx.y == 0)  // 4 threads compute
    {
      dyn::Array<Real, 2> arrA(3, 2);
      arrA = 3.5;
      for (int j = 0; j != 2; ++j)
        for (int i = 0; i != 3; ++i)
          {
            if (arrA[j][i] != 3.5) incrStat();
          }
    }
  UNINITIALIZED__shared__(dyn::Array<Real, 3> arrB);
  if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0)
    {
      arrB.INITIALIZEshared();
      arrB.define(4, 4, 4);
      arrB = 2.5;
    }
  __syncthreads();
  if (arrB(threadIdx.x, threadIdx.y, threadIdx.z) != 2.5) incrStat();
  __syncthreads();
  if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0)
    {
      arrB.define(stc::Vector<int, 4>{ -1, 0, 1 },
                  stc::Vector<int, 4>{  2, 3, 4 });
      arrB = 1.5;
    }
  __syncthreads();
  if (arrB(threadIdx.x - 1, threadIdx.y, threadIdx.z + 1) != 1.5) incrStat();
}
