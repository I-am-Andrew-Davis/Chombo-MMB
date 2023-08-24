#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


/*******************************************************************************
 * This is a rewrite of member functions in MOLAMRTimeDependent/MOLUtilities
 * that avoids a defined state, especially ProblemDomain which for MMB cannot
 * be generally defined.
 ******************************************************************************/

#ifndef _MOLUTILFUNC_H_
#define _MOLUTILFUNC_H_

#include "Box.H"
#include "IntVectSet.H"
#include "Vector.H"
#include "REAL.H"
#include "FArrayBox.H"
#include "FluxBox.H"

#include "NamespaceHeader.H"


/******************************************************************************/
/**
 * \file MOLUtilFunc.H
 *
 * \brief Utility functions for higher-order methods of lines
 *
 * slopes, parabolic interpolants, limiters.
 *
 *//*+*************************************************************************/

namespace MOLUtilFunc
{

  /// Compute the slope flattening coefficients
  void computeFlattening(FArrayBox&           a_flattening,
                         const FArrayBox&     a_W,
                         const Interval&      a_velInt,
                         const int&           a_presInd,
                         const Real&          a_smallPres,
                         const int&           a_bulkModulusInd,
                         const Box&           a_box,
                         const ProblemDomain& a_domain);

  /// Apply the flattening to slopes
  void applyFlattening(FArrayBox&       a_dW,
                       const FArrayBox& a_flat,
                       const Box&       a_box);

  /// PPM face-centered interpolant
  void PPMFaceValues(FArrayBox&           a_WFace,
                     const FArrayBox&     a_W,
                     const int&           a_numSlopes,
                     const int&           a_dir,
                     const Box&           a_box,
                     const ProblemDomain& a_domain);

  /// PPM Limiter
  void PPMLimiter(FArrayBox&           a_dWMinus,
                  FArrayBox&           a_dWPlus,
                  const FArrayBox&     a_W,
                  const int&           a_numSlopes,
                  const int&           a_dir,
                  const Box&           a_box,
                  const ProblemDomain& a_domain,
                  const int            a_useHOChecks,
                  const int            a_loBoundHOChecks,
                  const int            a_hiBoundHOChecks);

  /// PPM Limiter
  void PPMLimiter(FArrayBox&           a_dWMinus,
                  FArrayBox&           a_dWPlus,
                  const FArrayBox&     a_W,
                  const FArrayBox&     a_WfromU,
                  const int&           a_numSlopes,
                  const int&           a_dir,
                  const Box&           a_box,
                  const ProblemDomain& a_domain);


  /// Convolve or deconvolve on a cell
  void deconvolve(FArrayBox&           a_avgFab,
                  const FArrayBox&     a_cellFab,
                  const Box&           a_box,
                  const ProblemDomain& a_domain,
                  int                  a_sign = 1);

  /// Convolve or deconvolve on a face
  void deconvolveFace(FluxBox&             a_avgFlux,
                      const FluxBox&       a_cellFlux,
                      const Box&           a_box,
                      const ProblemDomain& a_domain,
                      int                  a_sign = 1);

  /// Convolve or deconvolve on a face described by an FarrayBox
  void deconvolveCenterFace(FArrayBox&           a_faceFab,
                            const FArrayBox&     a_faceExtFab,
                            const Box&           a_box,
                            const ProblemDomain& a_domain,
                            int                  a_dir,
                            int                  a_sign = 1);

}  // End of namespace MOLUtilFunc

#include "NamespaceFooter.H"

#endif  /* ! defined _MOLUTILFUNC_H_ */
