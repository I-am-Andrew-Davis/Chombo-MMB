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
 * \brief Kernels for testHash on the GPU
 *
 *//*+*************************************************************************/

#include "CH_Hash.H"

// The result
__device__ unsigned g_hash;

extern "C"
__global__ void
hash()
{
  if (threadIdx.x == 0)
    {
      g_hash = CH_Hash::hash_google_CityHash32(42);
      // printf("On GPU: %d\n", g_hash);
    }
}
