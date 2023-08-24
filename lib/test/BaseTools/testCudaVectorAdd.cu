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
 * \brief Kernels for testing the GPU
 *
 *//*+*************************************************************************/

// __device__ unsigned
// reallySetVal(unsigned a_val);

extern "C"
__global__ void
vvAddTo(float *const x, const float *const y)
{
   const int i = threadIdx.x + blockIdx.x*blockDim.x;
   x[i] += y[i];
}

extern "C"
__global__ void
vsAddTo(float *const x, const float s)
{
   // float ss = reallySetVal(s);
   const int i = threadIdx.x + blockIdx.x*blockDim.x;
   x[i] += s;
}
