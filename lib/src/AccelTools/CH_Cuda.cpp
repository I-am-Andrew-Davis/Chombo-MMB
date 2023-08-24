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
 * \brief Support macros, types, and routines for the Cuda driver framework
 *
 *//*+*************************************************************************/

#include "CH_Cuda.H"


//--Definition of static variables in CudaDriverSupport.H

char CH_Cuda::msgBuffer[msgBufferSize];


/*============================================================================*/
//  Finds a better SI prefix for a count
/**
 *  Converts a count of bytes or hertz to to kilo, mega, ... if required such
 *  that it is in an easily readable format
 *  \param[in]  a_unitType
 *                      0 - bytes
 *                      1 - Hz
 *  \param[in]  a_count Number of bytes to reduce
 *  \param[out] a_count Number in units returned
 *   return             Character pointer describing units of reduced
 *                      count.
 *//*=========================================================================*/

const char*
CH_Cuda::reduceUnitCount(const int a_unitType, float& a_count)
{
  static const char* unitSI[2][4] =
    { 
      { "bytes", "kilobytes", "megabytes", "gigabytes" },
      { "Hz", "KHz", "MHz", "GHz" }
    };
  const char* memUnit = unitSI[a_unitType][0];
  if ( a_count > 10000.f ) {
    a_count /= 1000.f;
    memUnit = unitSI[a_unitType][1];
  }
  if ( a_count > 10000.f ) {
    a_count /= 1000.f;
    memUnit = unitSI[a_unitType][2];
  }
  if ( a_count > 10000.f ) {
    a_count /= 1000.f;
    memUnit = unitSI[a_unitType][3];
  }
  return memUnit;
}
