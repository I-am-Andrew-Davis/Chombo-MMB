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
 * \brief Kernels for testGlobalStruct on the GPU
 *
 *//*+*************************************************************************/

#include "GlobalStruct.H"

__device__ bool g_verbose;

// Definition of globals on the device
__device__ CH_Hash::StaticTable<int, GlobalStruct>* g_GlobalStruct_table;

// Tracking of errors on the device
__device__ unsigned g_stat;
__device__ void incrStat()
{
  atomicInc(&g_stat, (unsigned)-1);
}

extern "C"
__global__ void
testGlobalStructGet()
{
  GlobalStruct& a_gs = GlobalStruct::get(1);
  if (threadIdx.x < 2)
    {
      GlobalStruct& a_gs = GlobalStruct::get(threadIdx.x + 1);
      const double expected = 1.1*(threadIdx.x + 1);
      if (a_gs.a_val != expected) incrStat();
      if (g_verbose)
        {
          printf("%d: %lf\n", threadIdx.x, a_gs.a_val);
        }
    }
}
