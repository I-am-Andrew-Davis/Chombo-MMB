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
 * \file PatchOps.H
 *
 * \brief Patch operators for post processing
 *
 *//*+*************************************************************************/

#ifndef _PATCHOPS_H_
#define _PATCHOPS_H_

//----- Standard Library -----//

//----- Chombo Library -----//

//----- Internal -----//

#include "DataTemp.H"

/*******************************************************************************
 */
///  Patch operators for post processing
/**
 *//*+*************************************************************************/

class PatchOps
{

/*====================================================================*
 * Types
 *====================================================================*/

public:

  using varStrVec = std::vector<std::string>;

/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Constructor
  PatchOps();

  /// Copy constructor not permitted
  PatchOps(const PatchOps&) = delete;

  /// Assignment constructor not permitted
  PatchOps& operator=(const PatchOps&) = delete;

  /// Destructor -- default destructor

/*==============================================================================
 * Public member functions
 *============================================================================*/

public:

  /// Weak construction
  void define();

  /// Convolve data -- fourth order accuracy
  void convolve4thOrd(FArrayBox&       a_convolvedState,
                      const FArrayBox& a_originalState,
                      const int        a_comp,
                      const Box&       a_box);

  /// Convolve data -- sixth order accuracy
  void convolve6thOrd(FArrayBox&       a_convolvedState,
                      const FArrayBox& a_originalState,
                      const int        a_comp,
                      const Box&       a_box);

  /// Deconvolve data -- fourth order accuracy
  void deconvolve4thOrd(FArrayBox&       a_deconvolvedState,
                        const FArrayBox& a_originalState,
                        const int        a_comp,
                        const Box&       a_box);

  /// Deconvolve data -- sixth order accuracy
  void deconvolve6thOrd(FArrayBox&       a_deconvolvedState,
                        const FArrayBox& a_originalState,
                        const int        a_comp,
                        const Box&       a_box);

  /// ConsToPrim -- density, velocity, pressure, temperature
  void consToPrim(FArrayBox&       a_primState,
                  const FArrayBox& a_consState,
                  const Box&       a_box);

  /// Gradient operator -- second order accuracy
  void gradOp2ndOrd(FArrayBox&      a_gradients,
                    FArrayBox&      a_inputState,
                    const RealVect& a_dx,
                    const Box&      a_box);

  /// Gradient operator -- fourth order accuracy
  void gradOp4thOrd(FArrayBox&      a_gradients,
                    FArrayBox&      a_inputState,
                    const RealVect& a_dx,
                    const Box&      a_box);

  /// Gradient operator -- sixth order accuracy
  void gradOp6thOrd(FArrayBox&      a_gradients,
                    FArrayBox&      a_inputState,
                    const RealVect& a_dx,
                    const Box&      a_box);

  /// Divergence operator -- derivatives available from gradient
  void divOpFromGrad(FArrayBox&       a_div,
                     const FArrayBox& a_gradients,
                     const Box&       a_box);

  /// Divergence operator -- calculate derivatives
  void divOp2ndOrd(FArrayBox&       a_div,
                     const FArrayBox& a_inputState,
                     const Box&       a_box);

  /// Divergence operator -- calculate derivatives
  void divOp4thOrd(FArrayBox&       a_div,
                     const FArrayBox& a_inputState,
                     const Box&       a_box);

  /// Divergence operator -- calculate derivatives
  void divOp6thOrd(FArrayBox&       a_div,
                     const FArrayBox& a_inputState,
                     const Box&       a_box);

  /// Curl operator -- derivatives availale from gradient
  void curlOpFromGrad(FArrayBox&       a_curl,
                      const FArrayBox& a_gradients,
                      const Box&       a_box);

  /// Curl operator -- calculate derivatives
  void curlOp2ndOrd(FArrayBox&       a_curl,
                      const FArrayBox& a_inputState,
                      const Box&       a_box);

  /// Curl operator -- calculate derivatives
  void curlOp4thOrd(FArrayBox&       a_curl,
                      const FArrayBox& a_inputState,
                      const Box&       a_box);

  /// Curl operator -- calculate derivatives
  void curlOp6thOrd(FArrayBox&       a_curl,
                      const FArrayBox& a_inputState,
                      const Box&       a_box);

  /// Laplacian operator -- compute divergence of gradients
  void lapOpFromGrad(FArrayBox&       a_laplacian,
                     const FArrayBox& a_gradients,
                     const Box&       a_box);

  /// Laplacian operator -- calculate derivatives from inputState
  void lapOp2ndOrd(FArrayBox&       a_laplacian,
                   const FArrayBox& a_inputState,
                   const Box&       a_box);

  /// Laplacian operator -- calculate derivatives from inputState
  void lapOp4thOrd(FArrayBox&       a_laplacian,
                   const FArrayBox& a_inputState,
                   const Box&       a_box);

  /// Laplacian operator -- calculate derivatives from inputState
  void lapOp6thOrd(FArrayBox&       a_laplacian,
                   const FArrayBox& a_inputState,
                   const Box&       a_box);

  /// Hessian matrix -- compute gradient of gradients
  void hessOpFromGrad(FArrayBox&       a_hessian,
                      const FArrayBox& a_gradients,
                      const Box&       a_box);

  /// Hessian matrix -- calculate derivatives from inputState
  void hessOp2ndOrd(FArrayBox&       a_hessian,
                    const FArrayBox& a_inputState,
                    const Box&       a_box);

  /// Hessian matrix -- calculate derivatives from inputState
  void hessOp4thOrd(FArrayBox&       a_hessian,
                    const FArrayBox& a_inputState,
                    const Box&       a_box);

  /// Hessian matrix -- calculate derivatives from inputState
  void hessOp6thOrd(FArrayBox&       a_hessian,
                    const FArrayBox& a_inputState,
                    const Box&       a_box);

  /// Matrix-Matrix, Matrix-Vector, Vector-Vector multiplication
  void matMatProdOp(FArrayBox&       a_product,
                    const FArrayBox& a_leftState,
                    const FArrayBox& a_rightState,
                    const Box&       a_box);

  /// Cross-product operator
  void crossProdOp(FArrayBox&       a_product,
                   const FArrayBox& a_leftState,
                   const FArrayBox& a_rightState,
                   const Box&       a_box);

  /// Domain wide summation
  void domSumOp(Real&            a_sum,
                const FArrayBox& a_inputState,
                const Box&       a_box);

  /// Magnitude of a vector quantity
  void magOp(FArrayBox&       a_magnitude,
             const FArrayBox& a_inputState,
             const Box&       a_box);

  /// Root mean square of a vector quantity
  void rmsOp(FArrayBox&       a_rms,
             const FArrayBox& a_inputState,
             const Box&       a_box);

/*==============================================================================
 * Protected member functions
 *============================================================================*/

protected:

  /// Any protected functions

/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  bool                 m_mapped;      ///< Simulation uses mapping
};

#endif  /* ! defined _LEVELPOSTPROCESSOP_H_ */
