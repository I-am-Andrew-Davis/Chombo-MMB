#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SINGLELEVELDIVERGENCE_H_
#define _SINGLELEVELDIVERGENCE_H_

#include "LevelData.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "LevelFluxRegister.H"
#include "QuadCFInterp.H"

#include "NamespaceHeader.H"

/// Class to encapsulate Divergence functions

class SingleLevelDivergence
{

public:

  /** computes cell-centered level-operator divergence of
      edge-centered vector field uEdge; assumes all coarse-fine BC's
      have already been set. */
  static void levelDivergenceMAC(
                                 ///
                                 LevelData<FArrayBox>& a_div,
                                 ///
                                 const LevelData<FluxBox>& a_uEdge,
                                 ///
                                 const Real a_dx);


  /// really basic single-fab MAC divergence
  static void simpleDivergenceMAC(///
                                  FArrayBox& a_div,
                                  ///
                                  const FluxBox& a_uEdge,
                                  ///
                                  const Real a_dx);
};

#include "NamespaceFooter.H"

#endif


