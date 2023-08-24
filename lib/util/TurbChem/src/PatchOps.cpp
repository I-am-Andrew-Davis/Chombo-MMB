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
 * \file PatchOps.cpp
 *
 * \brief Member functions for PatchOps
 *
 *//*+*************************************************************************/

//----- Standard Library -----//

//----- Chombo Library -----//

#include "AMR.H"
#include "AMRLevel.H"
#include "CONSTANTS.H"
#include "AMRIO.H"
#include "CH_HDF5.H"

//----- Internal -----//
#include "PatchOps.H"

/*******************************************************************************
 *
 * Class PatchOps: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** \param[in]  a_something 
 *                      Something
 *//*-----------------------------------------------------------------*/

PatchOps::PatchOps()
:
  m_mapped(false)
{
}


/*==============================================================================
 * Public member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Convolve the data (Taylor series of box filter for cell averaging)
/** \param[in]  a_level The level index
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *//*-----------------------------------------------------------------*/

void
PatchOps::define()
{
}

/*--------------------------------------------------------------------*/
//  Convolve the data (Taylor series of box filter for cell averaging)
/** \param[in]  a_level The level index
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *//*-----------------------------------------------------------------*/

void
PatchOps::convolve4thOrd(FArrayBox&       a_convolvedState,
                         const FArrayBox& a_originalState,
                         const int        a_comp,
                         const Box&       a_box)
{
  // Coefficient for 2nd-derivative in deconvolution expansion
  const Real delta = (1./24.);

  MD_ARRAY_RESTRICT(arrUavg, a_convolvedState);
  MD_ARRAY_RESTRICT(arrUpnt, a_originalState);
  MD_BOXLOOP(a_box, i)
    {
      Real Uj = arrUpnt[MD_IX(i, a_comp)];
      arrUavg[MD_IX(i,a_comp)] = Uj;
      for (int idir = 0; idir != SpaceDim; ++idir)
        {
          const int MD_ID(o, idir);
          Real Ujp = arrUpnt[MD_OFFSETIX(i,+,o,a_comp)];
          Real Ujm = arrUpnt[MD_OFFSETIX(i,-,o,a_comp)];
          Real d2x = Ujp - 2.*Uj + Ujm;
          arrUavg[MD_IX(i,a_comp)] += delta*d2x;
        }
    }
}

/*--------------------------------------------------------------------*/
//  Convolve the data (Taylor series of box filter for cell averaging)
/** \param[in]  a_level The level index
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *//*-----------------------------------------------------------------*/

void
PatchOps::convolve6thOrd(FArrayBox&       a_convolvedState,
                         const FArrayBox& a_originalState,
                         const int        a_comp,
                         const Box&       a_box)
{
}

/*--------------------------------------------------------------------*/
//  Deconvolve the data (inverse of convolve, i.e. de-averaging)
/** \param[in]  a_level The level index
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *//*-----------------------------------------------------------------*/

void
PatchOps::deconvolve4thOrd(FArrayBox&       a_deconvolvedState,
                           const FArrayBox& a_originalState,
                           const int        a_comp,
                           const Box&       a_box)
{
  // Coefficient for 2nd-derivative in deconvolution expansion
  const Real delta = (1./24.);

  MD_ARRAY_RESTRICT(arrUpnt, a_deconvolvedState);
  MD_ARRAY_RESTRICT(arrUavg, a_originalState);
  MD_BOXLOOP(a_box, i)
    {
      Real Uj = arrUavg[MD_IX(i, a_comp)];
      arrUpnt[MD_IX(i,a_comp)] = Uj;
      for (int idir = 0; idir != SpaceDim; ++idir)
        {
          const int MD_ID(o, idir);
          Real Ujp = arrUavg[MD_OFFSETIX(i,+,o,a_comp)];
          Real Ujm = arrUavg[MD_OFFSETIX(i,-,o,a_comp)];
          Real d2x = Ujp - 2.*Uj + Ujm;
          arrUpnt[MD_IX(i,a_comp)] -= delta*d2x;
        }
    }
}

/*--------------------------------------------------------------------*/
//  Deconvolve the data (inverse of convolve, i.e. de-averaging)
/** \param[in]  a_level The level index
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *//*-----------------------------------------------------------------*/

void
PatchOps::deconvolve6thOrd(FArrayBox&       a_deconvolvedState,
                           const FArrayBox& a_originalState,
                           const int        a_comp,
                           const Box&       a_box)
{
}

/*--------------------------------------------------------------------*/
//  Conservative to primitive variables
/** \param[in]  a_level The level index
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *//*-----------------------------------------------------------------*/

void
PatchOps::consToPrim(FArrayBox&       a_primState,
                     const FArrayBox& a_consState,
                     const Box&       a_box)
{
  // FIXME: create more general form of indexing to match
  //        what comes from the files (i.e. incompressible, etc.)
  // Indices
  const int rhoIndx = 0;
  const int velIndx = 1;
  const int pressIndx = velIndx + SpaceDim;
  const int tempIndx = pressIndx + 1;

  // Constants -- FIXME: should be read in from file
  //                     or set by user
  const Real gamma = 1.4;
  const Real R = 287.;

  // FIXME: account for multi-species flows
  MD_ARRAY_RESTRICT(arrW, a_primState);
  MD_ARRAY_RESTRICT(arrU, a_consState);
  MD_BOXLOOP(a_box, i)
    {
      Real rho  = arrU[MD_IX(i, rhoIndx)];
      D_TERM(Real u_mom = arrU[MD_IX(i, velIndx)];,
             Real v_mom = arrU[MD_IX(i, velIndx+1)];,
             Real w_mom = arrU[MD_IX(i, velIndx+2)];);
      Real eng = arrU[MD_IX(i, pressIndx)];

      // Copy the density
      arrW[MD_IX(i, rhoIndx)] = rho;

      // Compute the momentum
      D_TERM(
        arrW[MD_IX(i, velIndx)]   = u_mom/rho;,
        arrW[MD_IX(i, velIndx+1)] = v_mom/rho;,
        arrW[MD_IX(i, velIndx+2)] = w_mom/rho;);

      // Get the pressure
      arrW[MD_IX(i, pressIndx)] = (gamma - 1.)*(eng)
        - (gamma - 1.)*(1./2.)*((D_TERM(u_mom*u_mom,
                                      + v_mom*v_mom,
                                      + w_mom*w_mom))/rho);

      // Get the temperature
      arrW[MD_IX(i, tempIndx)] = arrW[MD_IX(i, pressIndx)]/(rho*R);
    }
}

/*--------------------------------------------------------------------*/
//  Compute the gradient of a scalar or vector
/** \param[in]  a_level The level index
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *//*-----------------------------------------------------------------*/

void
PatchOps::gradOp2ndOrd(FArrayBox&      a_gradients,
                       FArrayBox&      a_inputState,
                       const RealVect& a_dx,
                       const Box&      a_box)
{
  // Test the difference between
  //   <U>  -->  Grad(<U>)
  //   <U>  -->  U  -->  Grad(U)  -->  <Grad(U)>
  const int numInputComps = a_inputState.nComp();

  // Check to make sure that the input FArrayBox has enough components
  // to store the gradient of the input FArrayBox

  // Check to make sure that the two input FArrayBoxes have the same
  // size of box()

  MD_ARRAY_RESTRICT(arrGrad, a_gradients);
  MD_ARRAY_RESTRICT(arrIn, a_inputState);
  for (int j = 0; j != numInputComps; ++j)
    {
      for (int k = 0; k != SpaceDim; ++k)
        {
          // Gradients are stored in row vectors
          int gradIdx = j*SpaceDim + k;
          // Offset in direction of gradient evaluation
          const int MD_ID(ii, k);
          MD_BOXLOOP(a_box, i)
            {
              Real u_p  = arrIn[MD_OFFSETIX(i,+,ii,j)];
              Real u_m  = arrIn[MD_OFFSETIX(i,-,ii,j)];
              // 2nd-order gradient
              arrGrad[MD_IX(i, gradIdx)]
                = ((1./2.)*u_p - (1./2.)*u_m)/a_dx[k];
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Compute the gradient of a scalar or vector
/** \param[in]  a_level The level index
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *//*-----------------------------------------------------------------*/

void
PatchOps::gradOp4thOrd(FArrayBox&      a_gradients,
                       FArrayBox&      a_inputState,
                       const RealVect& a_dx,
                       const Box&      a_box)
{
  // Test the difference between
  //   <U>  -->  Grad(<U>)
  //   <U>  -->  U  -->  Grad(U)  -->  <Grad(U)>
  const int numInputComps = a_inputState.nComp();

  // Check to make sure that the input FArrayBox has enough components
  // to store the gradient of the input FArrayBox

  // Check to make sure that the two input FArrayBoxes have the same
  // size of box()

  MD_ARRAY_RESTRICT(arrGrad, a_gradients);
  MD_ARRAY_RESTRICT(arrIn, a_inputState);
  for (int j = 0; j != numInputComps; ++j)
    {
      for (int k = 0; k != SpaceDim; ++k)
        {
          // Gradients are stored in row vectors
          int gradIdx = j*SpaceDim + k;
          // Offset in direction of gradient evaluation
          const int MD_ID(ii, k);
          MD_BOXLOOP(a_box, i)
            {
              Real u_p  = arrIn[MD_OFFSETIX(i,+,ii,j)];
              Real u_pp = arrIn[MD_OFFSETIX(i,+,2*ii,j)];
              Real u_m  = arrIn[MD_OFFSETIX(i,-,ii,j)];
              Real u_mm = arrIn[MD_OFFSETIX(i,-,2*ii,j)];
              // 4th-order gradient
              arrGrad[MD_IX(i, gradIdx)]
                = ((1./12.)*u_mm - (2./3.)*u_m
                 + (2./3.)*u_p - (1./12.)*u_pp)/a_dx[k];
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Compute the gradient of a scalar or vector
/** \param[in]  a_level The level index
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *//*-----------------------------------------------------------------*/

void
PatchOps::gradOp6thOrd(FArrayBox&      a_gradients,
                       FArrayBox&      a_inputState,
                       const RealVect& a_dx,
                       const Box&      a_box)
{
  // Test the difference between
  //   <U>  -->  Grad(<U>)
  //   <U>  -->  U  -->  Grad(U)  -->  <Grad(U)>
  const int numInputComps = a_inputState.nComp();

  // Check to make sure that the input FArrayBox has enough components
  // to store the gradient of the input FArrayBox

  // Check to make sure that the two input FArrayBoxes have the same
  // size of box()

  MD_ARRAY_RESTRICT(arrGrad, a_gradients);
  MD_ARRAY_RESTRICT(arrIn, a_inputState);
  for (int j = 0; j != numInputComps; ++j)
    {
      for (int k = 0; k != SpaceDim; ++k)
        {
          // Gradients are stored in row vectors
          int gradIdx = j*SpaceDim + k;
          // Offset in direction of gradient evaluation
          const int MD_ID(ii, k);
          MD_BOXLOOP(a_box, i)
            {
              Real u_p   = arrIn[MD_OFFSETIX(i,+,ii,j)];
              Real u_pp  = arrIn[MD_OFFSETIX(i,+,2*ii,j)];
              Real u_ppp = arrIn[MD_OFFSETIX(i,+,3*ii,j)];
              Real u_m   = arrIn[MD_OFFSETIX(i,-,ii,j)];
              Real u_mm  = arrIn[MD_OFFSETIX(i,-,2*ii,j)];
              Real u_mmm = arrIn[MD_OFFSETIX(i,-,3*ii,j)];
              // 6th-order gradient
              arrGrad[MD_IX(i, gradIdx)]
                = ((1./60.)*u_ppp - (3./20.)*u_pp + (3./4.)*u_p
                   - (3./4.)*u_m + (3./20.)*u_mm - (1./60.)*u_mmm)/a_dx[k];
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Compute the divergence of a vector or matrix
/** \param[in]  a_level The level index
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *//*-----------------------------------------------------------------*/

void
PatchOps::divOpFromGrad(FArrayBox&       a_div,
                        const FArrayBox& a_gradients,
                        const Box&       a_box)
{
}

/*--------------------------------------------------------------------*/
//  Compute the divergence of a vector or matrix
/** \param[in]  a_level The level index
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *//*-----------------------------------------------------------------*/

void
PatchOps::divOp2ndOrd(FArrayBox&       a_div,
                      const FArrayBox& a_inputState,
                      const Box&       a_box)
{
}

/*--------------------------------------------------------------------*/
//  Compute the divergence of a vector or matrix
/** \param[in]  a_level The level index
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *//*-----------------------------------------------------------------*/

void
PatchOps::divOp4thOrd(FArrayBox&       a_div,
                      const FArrayBox& a_inputState,
                      const Box&       a_box)
{
}

/*--------------------------------------------------------------------*/
//  Compute the divergence of a vector or matrix
/** \param[in]  a_level The level index
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *//*-----------------------------------------------------------------*/

void
PatchOps::divOp6thOrd(FArrayBox&       a_div,
                      const FArrayBox& a_inputState,
                      const Box&       a_box)
{
}

/*--------------------------------------------------------------------*/
//  Define the level operator (weak construction)
/** \param[in]  a_level The level index
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *//*-----------------------------------------------------------------*/

void
PatchOps::curlOpFromGrad(FArrayBox&       a_curl,
                         const FArrayBox& a_gradients,
                         const Box&       a_box)
{
  // FIXME: this is set up only for 3D cases
  // FIXME: this requires the calculation of gradients, if gradients
  //        are not required elsewhere, check to see if only specific
  //        gradients can be calculated here from other input data

  MD_ARRAY_RESTRICT(arrCurl, a_curl);
  MD_ARRAY_RESTRICT(arrGrad, a_gradients);
  MD_BOXLOOP(a_box, i)
    {
      Real U_y = arrGrad[MD_IX(i, 1)];
      Real U_z = arrGrad[MD_IX(i, 2)];
      Real V_x = arrGrad[MD_IX(i, 3)];
      Real V_z = arrGrad[MD_IX(i, 5)];
      Real W_x = arrGrad[MD_IX(i, 6)];
      Real W_y = arrGrad[MD_IX(i, 7)];

      arrCurl[MD_IX(i, 0)] = W_y - V_z;
      arrCurl[MD_IX(i, 1)] = U_z - W_x;
      arrCurl[MD_IX(i, 2)] = V_x - U_y;
    }
}

/*--------------------------------------------------------------------*/
//  Define the level operator (weak construction)
/** \param[in]  a_level The level index
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *//*-----------------------------------------------------------------*/

void
PatchOps::curlOp2ndOrd(FArrayBox&       a_curl,
                       const FArrayBox& a_inputState,
                       const Box&       a_box)
{
}

/*--------------------------------------------------------------------*/
//  Define the level operator (weak construction)
/** \param[in]  a_level The level index
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *//*-----------------------------------------------------------------*/

void
PatchOps::curlOp4thOrd(FArrayBox&       a_curl,
                       const FArrayBox& a_inputState,
                       const Box&       a_box)
{
}

/*--------------------------------------------------------------------*/
//  Define the level operator (weak construction)
/** \param[in]  a_level The level index
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *//*-----------------------------------------------------------------*/

void
PatchOps::curlOp6thOrd(FArrayBox&       a_curl,
                       const FArrayBox& a_inputState,
                       const Box&       a_box)
{
}

/*--------------------------------------------------------------------*/
//  Define the level operator (weak construction)
/** \param[in]  a_level The level index
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *//*-----------------------------------------------------------------*/

void
PatchOps::lapOpFromGrad(FArrayBox&       a_laplacian,
                        const FArrayBox& a_gradients,
                        const Box&       a_box)
{
}

/*--------------------------------------------------------------------*/
//  Define the level operator (weak construction)
/** \param[in]  a_level The level index
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *//*-----------------------------------------------------------------*/

void
PatchOps::lapOp2ndOrd(FArrayBox&       a_laplacian,
                      const FArrayBox& a_inputState,
                      const Box&       a_box)
{
}

/*--------------------------------------------------------------------*/
//  Define the level operator (weak construction)
/** \param[in]  a_level The level index
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *//*-----------------------------------------------------------------*/

void
PatchOps::lapOp4thOrd(FArrayBox&       a_laplacian,
                      const FArrayBox& a_inputState,
                      const Box&       a_box)
{
}

/*--------------------------------------------------------------------*/
//  Define the level operator (weak construction)
/** \param[in]  a_level The level index
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *//*-----------------------------------------------------------------*/

void
PatchOps::lapOp6thOrd(FArrayBox&       a_laplacian,
                      const FArrayBox& a_inputState,
                      const Box&       a_box)
{
}

/*--------------------------------------------------------------------*/
//  Define the level operator (weak construction)
/** \param[in]  a_level The level index
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *//*-----------------------------------------------------------------*/

void
PatchOps::hessOpFromGrad(FArrayBox&       a_hessian,
                         const FArrayBox& a_gradients,
                         const Box&       a_box)
{
}

/*--------------------------------------------------------------------*/
//  Define the level operator (weak construction)
/** \param[in]  a_level The level index
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *//*-----------------------------------------------------------------*/

void
PatchOps::hessOp2ndOrd(FArrayBox&       a_hessian,
                       const FArrayBox& a_inputState,
                       const Box&       a_box)
{
}

/*--------------------------------------------------------------------*/
//  Define the level operator (weak construction)
/** \param[in]  a_level The level index
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *//*-----------------------------------------------------------------*/

void
PatchOps::hessOp4thOrd(FArrayBox&       a_hessian,
                       const FArrayBox& a_inputState,
                       const Box&       a_box)
{
}

/*--------------------------------------------------------------------*/
//  Define the level operator (weak construction)
/** \param[in]  a_level The level index
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *//*-----------------------------------------------------------------*/

void
PatchOps::hessOp6thOrd(FArrayBox&       a_hessian,
                       const FArrayBox& a_inputState,
                       const Box&       a_box)
{
}

/*--------------------------------------------------------------------*/
//  Define the level operator (weak construction)
/** \param[in]  a_level The level index
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *//*-----------------------------------------------------------------*/

void
PatchOps::matMatProdOp(FArrayBox&       a_product,
                       const FArrayBox& a_leftState,
                       const FArrayBox& a_rightState,
                       const Box&       a_box)
{
}

/*--------------------------------------------------------------------*/
//  Define the level operator (weak construction)
/** \param[in]  a_level The level index
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *//*-----------------------------------------------------------------*/

void
PatchOps::crossProdOp(FArrayBox&       a_product,
                      const FArrayBox& a_leftState,
                      const FArrayBox& a_rightState,
                      const Box&       a_box)
{
}

/*--------------------------------------------------------------------*/
//  Define the level operator (weak construction)
/** \param[in]  a_level The level index
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *//*-----------------------------------------------------------------*/

void
PatchOps::domSumOp(Real&            a_sum,
                   const FArrayBox& a_inputState,
                   const Box&       a_box)
{
  // Use the sum function from FArrayBox class
  int component = 0;
  int numComponents = 1;
  a_sum = a_inputState.sum(a_box,
                           component,
                           numComponents);
}

/*--------------------------------------------------------------------*/
//  Define the level operator (weak construction)
/** \param[in]  a_level The level index
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *//*-----------------------------------------------------------------*/

void
PatchOps::magOp(FArrayBox&       a_magnitude,
                const FArrayBox& a_inputState,
                const Box&       a_box)
{
}

/*--------------------------------------------------------------------*/
//  Define the level operator (weak construction)
/** \param[in]  a_level The level index
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *//*-----------------------------------------------------------------*/

void
PatchOps::rmsOp(FArrayBox&       a_rms,
                const FArrayBox& a_inputState,
                const Box&       a_box)
{
}
