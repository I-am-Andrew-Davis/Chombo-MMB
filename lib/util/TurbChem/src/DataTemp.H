#ifdef CH_LANG_CC
/*
 *      _______               __
 *     / ___/ /  ___  __  ___/ /
 *    / /__/ _ \/ _ \/ _\/ _  /
 *    \___/_//_/\___/_/  \_._/
 *    Please refer to Copyright.txt, in Chord's root directory.
 */
#endif


/******************************************************************************/
/**
 * \file DataTemp.H
 *
 * \brief Abstractions for temporary data
 *
 *//*+*************************************************************************/

#ifndef _DATATEMP_H_
#define _DATATEMP_H_

//----- Chombo Library -----//

#include "FArrayBox.H"
#include "CHArray.H"

//--Temporary FAB

// #define USE_STACK
#ifdef USE_STACK

// Defines FAB on the stack
#define FABSTACKTEMP(fabname, box, ncomp)                               \
  Real fabname ## _smem[box.size().product()*ncomp];                    \
  FArrayBox fabname(box, ncomp, fabname ## _smem)


//--Fluxbox

#define FLUXBOXSTACKTEMP(flxname, box, ncomp)                           \
  D_TERM6(                                                              \
    Real flxname ## _smem0[surroundingNodes(box, 0).size().product()*ncomp];, \
    Real flxname ## _smem1[surroundingNodes(box, 1).size().product()*ncomp];, \
    Real flxname ## _smem2[surroundingNodes(box, 2).size().product()*ncomp];, \
    Real flxname ## _smem3[surroundingNodes(box, 3).size().product()*ncomp];, \
    Real flxname ## _smem4[surroundingNodes(box, 4).size().product()*ncomp];, \
    Real flxname ## _smem5[surroundingNodes(box, 5).size().product()*ncomp];) \
  FluxBox flxname(box, ncomp, D_DECL6(flxname ## _smem0,                \
                                      flxname ## _smem1,                \
                                      flxname ## _smem2,                \
                                      flxname ## _smem3,                \
                                      flxname ## _smem4,                \
                                      flxname ## _smem5))

//--CHArray

// Defines CHArray (BNN) on the stack
#define CHARRAYBNNSTACKTEMP(arrname, T, rank, config, box, dimc1, dimc0) \
  T arrname ## _smem[box.size().product()*(dimc1)*(dimc0)];             \
  CHArray<T, rank, config> arrname;                                     \
  arrname.define(arrname ## _smem, box, (dimc1), (dimc0))

// Defines CHArray (BN) on the stack
#define CHARRAYBNSTACKTEMP(arrname, T, rank, config, box, dimc0)        \
  T arrname ## _smem[box.size().product()*(dimc0)];                     \
  CHArray<T, rank, config> arrname;                                     \
  arrname.define(arrname ## _smem, box, (dimc0))

// Defines CHArray (NNB) on the stack
#define CHARRAYNNBSTACKTEMP(arrname, T, rank, config, dimc1, dimc0, box) \
  T arrname ## _smem[(dimc1)*(dimc0)*box.size().product()];             \
  CHArray<T, rank, config> arrname;                                     \
  arrname.define(arrname ## _smem, (dimc1), (dimc0), box)

// Defines CHArray (NB) on the stack
#define CHARRAYNBSTACKTEMP(arrname, T, rank, config, dimc0, box)        \
  T arrname ## _smem[(dimc0)*box.size().product()];                     \
  CHArray<T, rank, config> arrname;                                     \
  arrname.define(arrname ## _smem, (dimc0), box)

#else

// Allocate FAB on the heap
#define FABSTACKTEMP(fabname, box, ncomp)                               \
  FArrayBox fabname(box, ncomp)

//--Fluxbox

#define FLUXBOXSTACKTEMP(flxname, box, ncomp)                           \
  FluxBox flxname(box, ncomp)

//--CHArray

// Defines CHArray (BNN) on the stack
#define CHARRAYBNNSTACKTEMP(arrname, T, rank, config, box, dimc1, dimc0) \
  CHArray<T, rank, config> arrname(box, (dimc1), (dimc0))

// Defines CHArray (BN) on the stack
#define CHARRAYBNSTACKTEMP(arrname, T, rank, config, box, dimc0)        \
  CHArray<T, rank, config> arrname(box, (dimc0))

// Defines CHArray (NNB) on the stack
#define CHARRAYNNBSTACKTEMP(arrname, T, rank, config, dimc1, dimc0, box) \
  CHArray<T, rank, config> arrname((dimc1), (dimc0), box)

// Defines CHArray (NB) on the heap
#define CHARRAYNBSTACKTEMP(arrname, T, rank, config, dimc0, box)        \
  CHArray<T, rank, config> arrname((dimc0), box)

#endif


#endif /* ! defined _DATATEMP_H_ */
