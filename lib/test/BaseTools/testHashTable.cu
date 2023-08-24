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

#include "CH_HashTable.H"

__device__ bool g_verbose;

// Tracking of errors on the device
__device__ unsigned g_stat;
__device__ void incrStat()
{
  atomicInc(&g_stat, (unsigned)-1);
}

// Record of device pointer from the CPU (to make sure the same ones are used
// here and that we are not simply accessing host memory)
__device__ uintptr_t g_devicePtr[2];

__device__ constexpr void
getValue(const int a_idx, CH_Hash::pair<int, double>& a_val)
{
  constexpr int key[6] = { 0, 9, 19, 8, 18, 1};
  constexpr double mapped[6] = { 0.0, 0.9, 0.19, 0.8, 0.18, 0.1};
  a_val = CH_Hash::make_pair(key[a_idx], mapped[a_idx]);
}

extern "C"
__global__ void
hashGraphFind(CH_Hash::HashGraphTable<int, double> a_table)
{
  __shared__ unsigned scratch[6];  // hash table index
  if (threadIdx.x == 0)
    {
      if (g_verbose)
        {
          printf("Table pointer locations on GPU:\n");
          printf("m_elem  : %p\n", a_table.m_elem.data());
          printf("m_offset: %p\n", a_table.m_offset.data());
        }
      if (reinterpret_cast<uintptr_t>(a_table.m_elem.data()) !=
          g_devicePtr[0]) incrStat();
      if (reinterpret_cast<uintptr_t>(a_table.m_offset.data()) !=
          g_devicePtr[1]) incrStat();
    }

  // 1 thread reads per query
  if (threadIdx.x < 6)
    {
      CH_Hash::pair<int, double> value;
      getValue(threadIdx.x, value);
      auto iter = a_table.find(value.first);
      if (iter == a_table.end()) incrStat();
      if (iter->second != value.second) incrStat();
    }
  // 4 threads read per query
  if (threadIdx.x < 4*6)
    {
      CH_Hash::pair<int, double> value;
      getValue(threadIdx.x/4, value);
      auto iter = a_table.find<4>(value.first,
                                  threadIdx.x%4,
                                  &scratch[threadIdx.x/4]);
      if (iter == a_table.end()) incrStat();
      if (iter->second != value.second) incrStat();
    }
}
