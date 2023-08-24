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

#include "LoHiSide.H"
#include "LoHiCenter.H"
#include "MOLPhysics.H"
#include "GodunovUtilitiesF_F.H"
#include "DebugOut.H"
#include "MOLUtilFunc.H"

#include "NamespaceHeader.H"

/*--------------------------------------------------------------------*/
//  Compute the slope flattening coefficients
/** Compute the slope flattening coefficients, a_flattening, using the
 *  primitive variables, a_W, within a_box.
 *  \param[out] a_flattening
 *                      Flattening coeffs, 1 component on a_box
 *  \param[in]  a_W     Primitive variables, on a_box grown by 3 within
 *                      a_domain
 *  \param[in]  a_velInt
 *                      Interval of a_W with velocity; length SpaceDim
 *  \param[in]  a_presInd
 *                      Index of a_W with pressure
 *  \param[in]  a_smallPres
 *                      Minimum pressure to use in ratio
 *  \param[in]  a_bulkModulusInd
 *                      Index of a_W with bulk modulus
 *  \param[in]  a_box   Cell-centered box on which a_flattening lives
 *  \param[in]  a_domain
 *                      Problem or block domain
 *//*-----------------------------------------------------------------*/

void
MOLUtilFunc::computeFlattening(FArrayBox&           a_flattening,
                               const FArrayBox&     a_W,
                               const Interval&      a_velInt,
                               const int&           a_presInd,
                               const Real&          a_smallPres,
                               const int&           a_bulkModulusInd,
                               const Box&           a_box,
                               const ProblemDomain& a_domain)
{
  CH_TIME("MOLUtilFunc::computeFlattening");
  CH_assert(a_W.box().contains(a_box));

  // The current direction
  int idir;

  // The directional flattening coefficients
  FArrayBox zetaDir(a_box,SpaceDim);

  // The divergence of the velocity
  FArrayBox dVel(a_box,SpaceDim);

  // The interval of the primitive variables corresponding to the velocity
  Interval velInterval= a_velInt;
  int v0index = velInterval.begin();

  // Get the directional flattening coefficients in each direction
  for (idir = 0; idir < SpaceDim; idir++)
  {
    // A box one larger (in direction "idir") than the final result box
    Box box1 = a_box;
    box1.grow(idir,1);

    // A box two larger (in direction "idir") than the final result box
    Box box2 = a_box;
    box2.grow(idir,2);

    // A box three larger (in direction "idir") than the final result box
    Box box3 = a_box;
    box3.grow(idir,3);

    // The primitive variables need to be defined over the largest box.
    // This assert may fail when a_W has been defined on a box that's grown
    // by 3 but then intersected with the domain.
    // CH_assert(a_W.box().contains(box3));

    // *** Compute delta1p, the first differences in "pressure"

    // Compute where centered differences can be used and where one-sided
    // differences need to be used.  The data used for the differences is
    // defined on cell-centered "box3".
    Box loBox, hiBox, centerBox, entireBox;
    int hasLo, hasHi;

    loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
               box3, a_domain, idir);

    // Set delta1p to cell-centered difference of cell-centered pressure
    // in direction idir.  Result is on entireBox.
    // For i in centerBox, delta1p[i] = (pressure[i+1] - pressure[i-1])/2.
    // For i in loBox, delta1p[i] = pressure[i+1] - pressure[i].
    // For i in hiBox, delta1p[i] = pressure[i] - pressure[i-1].

    FArrayBox delta1p(entireBox,1);
    FORT_GETGRADF(CHF_FRA1(delta1p,0),
                  CHF_CONST_FRA1(a_W,a_presInd),
                  CHF_CONST_INT(idir),
                  CHF_BOX(loBox),
                  CHF_CONST_INT(hasLo),
                  CHF_BOX(hiBox),
                  CHF_CONST_INT(hasHi),
                  CHF_BOX(centerBox));

    // *** Compute delta2p, the second differences in "pressure"

    // Compute where centered differences can be used and where one-sided
    // differences need to be used.  The data used for the differences is
    // defined on cell-centered "box2"
    loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
               box2, a_domain, idir);

    // Set delta2p to twice the cell-centered average of cell-centered delta1p
    // in direction idir.  Result is on entireBox.
    // For i in centerBox, delta2p[i] = delta1p[i-1] + delta1p[i+1].
    // For i in loBox, delta2p[i] = delta1p[i] + delta1p[i+1].
    // For i in hiBox, delta2p[i] = delta1p[i-1] + delta1p[i].

    FArrayBox delta2p(entireBox,1);
    FORT_GETDPTWOF(CHF_FRA1(delta2p,0),
                   CHF_CONST_FRA1(delta1p,0),
                   CHF_CONST_INT(idir),
                   CHF_BOX(loBox),
                   CHF_CONST_INT(hasLo),
                   CHF_BOX(hiBox),
                   CHF_CONST_INT(hasHi),
                   CHF_BOX(centerBox));

    // *** Compute bulkMin, a 3-way minimum of the "bulk modulus"

    // Same loBox, centerBox, hiBox, and entireBox as delta2p.

    // Set bulkMin to minimum of cell-centered bulk in neighborhood.
    // Result is on entireBox.
    // For i in centerBox, bulkMin[i] = min{bulk[i-1], bulk[i], bulk[i+1]}.
    // For i in loBox, bulkMin[i] = min{bulk[i], bulk[i+1]}.
    // For i in hiBox, bulkMin[i] = min{bulk[i-1], bulk[i]}.

    FArrayBox bulkMin(entireBox,1);
    int bulkIndex = a_bulkModulusInd;
    FORT_MIN3PTSF(CHF_FRA1(bulkMin,0),
                  CHF_CONST_FRA1(a_W,bulkIndex),
                  CHF_CONST_INT(idir),
                  CHF_BOX(loBox),
                  CHF_CONST_INT(hasLo),
                  CHF_BOX(hiBox),
                  CHF_CONST_INT(hasHi),
                  CHF_BOX(centerBox));

    // *** Use the first and second differences normalized by the
    // 3-way minimum, bulkMin, to generate zetaTwiddleDir,
    // flattening coefficients in direction idir.

    // Same loBox, centerBox, hiBox, and entireBox as delta2p and bulkMin.

    // if abs(delta1p[i] / bulkMin[i]) < 0.33,
    //    zetaTwiddleDir[i] = 1.;
    // else
    //    ratio[i] = abs(delta1p[i]) / max{abs(delta2p[i]), smallp};
    //    if ratio[i] <= 0.75, zetaTwiddleDir[i] = 1.;
    //    elseif ratio[i] >= 0.85, zetaTwiddleDir[i] = 0.;
    //    else zetaTwiddleDir[i] = 1. - (ratio[i] - 0.75) / (0.85 - 0.75).

    FArrayBox zetaTwiddleDir(entireBox,1);
    FORT_GETFLATF(CHF_FRA1(zetaTwiddleDir,0),
                  CHF_CONST_FRA1(delta1p,0),
                  CHF_CONST_FRA1(delta2p,0),
                  CHF_CONST_REAL(a_smallPres),
                  CHF_CONST_FRA1(bulkMin,0),
                  CHF_BOX(entireBox));

    // *** Compute component idir of zetaDir, a 3-way minimum
    // of the directional flattening coefficients, zetaTwiddleDir.

    // Compute where centered differences can be used and where one sided
    // differences need to be used.  The data used for the differences is
    // defined on cell-centered "box1"
    loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
               box1, a_domain, idir);

    // Set zetaDir to minimum of cell-centered zetTwiddleDir in neighborhood.
    // Result is on entireBox.
    // For i in centerBox, zetaDir[i] = min{zTD[i-1], zTD[i], zTD[i+1]}.
    // For i in loBox, zetaDir[i] = min{zTD[i], zTD[i+1]}.
    // For i in hiBox, zetaDir[i] = min{zTD[i-1], zTD[i]}.
    FORT_MIN3PTSF(CHF_FRA1(zetaDir,idir),
                  CHF_CONST_FRA1(zetaTwiddleDir,0),
                  CHF_CONST_INT(idir),
                  CHF_BOX(loBox),
                  CHF_CONST_INT(hasLo),
                  CHF_BOX(hiBox),
                  CHF_CONST_INT(hasHi),
                  CHF_BOX(centerBox));

    // *** Compute component idir of dVel, the divergence of the velocity

    // Set dVel to cell-centered difference of cell-centered pressure
    // in direction idir.  Result is on entireBox.
    // For i in centerBox, dVel[i] = (velocity[i+1] - velocity[i-1])/2.
    // For i in loBox, dVel[i] = velocity[i+1] - velocity[i].
    // For i in hiBox, dVel[i] = velocity[i] - velocity[i-1].
    FORT_GETGRADF(CHF_FRA1(dVel,idir),
                  CHF_CONST_FRA1(a_W,v0index+idir),
                  CHF_CONST_INT(idir),
                  CHF_BOX(loBox),
                  CHF_CONST_INT(hasLo),
                  CHF_BOX(hiBox),
                  CHF_CONST_INT(hasHi),
                  CHF_BOX(centerBox));
  }

  // *** Where divergence of velocity is negative, set a_flattening
  // to minimum of directional flattening coefficients zetaDir at that point

  // sum{dVel[i]} is divergence of dVel[i] (multiplied by mesh spacing).
  //
  // if sum{dVel[i]} >= 0,
  //    a_flattening[i] = 1.;
  // else
  //    a_flattening[i] = min{zetaDir[i]}.

  // At each point, set the flattening coefficient to the minimum of all
  // the directional flattening coefficients if the divergence of the velocity
  // is negative, otherwise set it to 1 (no flattening).
  FORT_MINFLATF(CHF_FRA1(a_flattening,0),
                CHF_CONST_FRA(zetaDir),
                CHF_CONST_FRA(dVel),
                CHF_BOX(a_box));
}

/*--------------------------------------------------------------------*/
//  Apply the flattening to slopes
/** Multiply every component of a_dW by a_flat, on a_box.
 *  \param[in]  a_dW    slopes to be flattened, on at least a_box
 *  \param[out] a_dW    flattened slopes
 *  \param[in]  a_flat  Cell-centered flattening coefficients, on at
 *                      least a_box
 *  \param[in]  a_box   Cell-centered box on which to flatten
 *//*-----------------------------------------------------------------*/

void
MOLUtilFunc::applyFlattening(FArrayBox&       a_dW,
                             const FArrayBox& a_flat,
                             const Box&       a_box)
{
  CH_TIME("MOLUtilFunc::applyFlattening");
  CH_assert(a_dW.box().contains(a_box));
  CH_assert(a_flat.box().contains(a_box));

  int numSlopes = a_dW.nComp();

  for (int islope = 0;islope < numSlopes;islope++)
    {
      a_dW.mult(a_flat,a_box,0,islope);
    }
}

/*--------------------------------------------------------------------*/
//  PPM face-centered interpolant
/** Given the cell average a_W, compute fourth-order accurate face-
 *  centered values WFace on a_box by differentiating the indefinite
 *  integral. Limiting is performed in a separate pass.
 *//*-----------------------------------------------------------------*/

void
MOLUtilFunc::PPMFaceValues(FArrayBox&           a_WFace,
                           const FArrayBox&     a_W,
                           const int&           a_numSlopes,
                           const int&           a_dir,
                           const Box&           a_box,
                           const ProblemDomain& a_domain)
{
  CH_TIME("MOLUtilFunc::PPMFaceValues");
  // a_box is the face-centered box on which a_WFace is computed.
  CH_assert(a_WFace.box().contains(a_box));

  // petermc, 7 Jan 2010: changed from a_box to grown a_box.
  Box box1cells = grow(a_box, BASISV(a_dir));
  box1cells.enclosedCells();

  Box loFaces, nextLoFaces;
  Box hiFaces, nextHiFaces;
  Box centerFaces, innerCenterFaces, entireFaces;
  int hasLoFaces, hasHiFaces;
  loHiCenterFace4(loFaces, nextLoFaces, hasLoFaces,
                  hiFaces, nextHiFaces, hasHiFaces,
                  centerFaces, innerCenterFaces, entireFaces,
                  box1cells, a_domain, a_dir);

  // For i-e/2 in innerCenterFaces, set a_WFace[i-e/2] from a_W[i-2e:i+e].
  // If hasLoFaces:
  // For i-e/2 in loFaces, set a_WFace[i-e/2] from a_W[i:i+3e].
  //           in nextLoFaces,                from a_W[i-e:i+2e].
  // If hasHiFaces:
  // For i-e/2 in hiFaces, set a_WFace[i-e/2] from a_W[i-4e:i-e].
  //           in nextHiFaces,                from a_W[i-3e:i].
  FORT_FOURTHINTERPFACES(CHF_FRA(a_WFace),
                         CHF_CONST_FRA(a_W),
                         CHF_CONST_INT(a_numSlopes),
                         CHF_CONST_INT(a_dir),
                         CHF_BOX(loFaces),
                         CHF_BOX(nextLoFaces),
                         CHF_CONST_INT(hasLoFaces),
                         CHF_BOX(hiFaces),
                         CHF_BOX(nextHiFaces),
                         CHF_CONST_INT(hasHiFaces),
                         CHF_BOX(innerCenterFaces));
}

/*--------------------------------------------------------------------*/
//  PPM Limiter
/** a_dWMinus and a_dWPlus are the differences between the face values
 *  on the minus and plus sides of cells and the average in the cell.
 *  That is,
 *    a_dWMinus[i] = WFace[i - e/2] - a_W[i]
 *    a_dWPlus[i] = WFace[i + e/2] - a_W[i]
 *  where e is the unit vector in dimension a_dir.
 *  The PPM limiter is applied to these values to obtain a monotone
 *  interpolant in the cell.
 *  The function returns the limited a_dWMinus and a_dWPlus on a_box.
 *  petermc, 4 Sep 2008:  included a_W in argument list
 *
 *  Need a_dWMinus and a_dWPlus on a_box,
 *  and need a_W on on a_box grown by 3 in dimension a_dir.
 *  Returns limited a_dWMinus and a_dWPlus on a_box.
 *//*-----------------------------------------------------------------*/

void
MOLUtilFunc::PPMLimiter(FArrayBox&           a_dWMinus,
                        FArrayBox&           a_dWPlus,
                        const FArrayBox&     a_W,
                        const int&           a_numSlopes,
                        const int&           a_dir,
                        const Box&           a_box,
                        const ProblemDomain& a_domain,
                        const int            a_useHOChecks,
                        const int            a_loBoundHOChecks,
                        const int            a_hiBoundHOChecks)
{
  CH_TIME("MOLUtilFunc::PPMLimiter");
  // a_dWMinus[i] = WFace[i - e/2] - a_W[i]
  // a_dWPlus[i] = WFace[i + e/2] - a_W[i]
  // where e is unit vector in dimension a_dir.

  // this option added by petermc, 5 Sep 2008
  // Will need to recalculate some D^2's.

  // We calculate d2Wfcf on a_box,
  // and we need a_dWMinus on a_box,
  // and we need a_dWPlus on a_box.
  // d2Wfcf[i] = 6 * (a_dWMinus[i] + a_dWPlus[i])
  //           = 6 * (thisFaceWDir[i-e/2] - a_cellW[i] +
  //                  thisFaceWDir[i+e/2] - a_cellW[i])
  FArrayBox d2Wfcf(a_box, a_numSlopes);
  d2Wfcf.copy(a_dWMinus);
  d2Wfcf.plus(a_dWPlus, 0, 0, a_numSlopes);
  d2Wfcf *= 6.;

  // petermc, 21 Sep 2010:
  // In order to get a_dWMinus and a_dWPlus on a_box,
  // we need d2W on a_box grown by 3 in a_dir directions.
  Box box3 = a_box;
  box3.grow(a_dir, 3);

  Box loBox, hiBox, centerBox, entireBox;
  int hasLo, hasHi;
  loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
             box3, a_domain, a_dir);

  FArrayBox d2W(entireBox, a_numSlopes);
  // On centerBox, use 3-point stencil for d2W;
  // on loBox and hiBox, copy result from neighbor.
  // In the case of centerBox == entireBox == box1,
  // where box1 is a_box grown by 1 in dimension a_dir:
  // We calculate d2W on a_box grown by 2 (was 1) in dimension a_dir,
  // and we need a_W on a_box grown by 3 (was 2) in dimension a_dir.

  // petermc, 21 Sep 2010, changed layer of d2W from 2 to 1,
  // and of a_W from 3 to 2.
  FORT_GETSECONDDIFF(CHF_FRA(d2W),
                     CHF_CONST_FRA(a_W),
                     CHF_CONST_INT(a_numSlopes),
                     CHF_CONST_INT(a_dir),
                     CHF_BOX(loBox),
                     CHF_CONST_INT(hasLo),
                     CHF_BOX(hiBox),
                     CHF_CONST_INT(hasHi),
                     CHF_BOX(centerBox));

  Box box1 = a_box;
  box1.grow(a_dir, 1);
  Box nextLoBox, nextHiBox, innerCenterBox;
  loHiCenter5(loBox, nextLoBox, hasLo,
              hiBox, nextHiBox, hasHi,
              centerBox, innerCenterBox, entireBox,
              box1, a_domain, a_dir);
  CH_assert(entireBox == a_box);

  Real limitC = 1.25;
  Real eps = 1.0e-12;
  Real c3 = 0.1;
  FORT_CHECKCUBICLIMITERF(// <W>_(i-1/2) - <W>_i, on a_box
    CHF_FRA(a_dWMinus),
    // <W>_(i+1/2) - <W>_i, on a_box
    CHF_FRA(a_dWPlus),
    CHF_CONST_FRA(a_W),
    CHF_CONST_FRA(d2W),
    CHF_CONST_FRA(d2Wfcf),
    CHF_CONST_INT(a_numSlopes),
    CHF_CONST_INT(a_dir),
    CHF_BOX(loBox),
    CHF_BOX(nextLoBox),
    CHF_CONST_INT(hasLo),
    CHF_BOX(hiBox),
    CHF_BOX(nextHiBox),
    CHF_CONST_INT(hasHi),
    CHF_BOX(innerCenterBox),
    CHF_CONST_REAL(limitC),
    CHF_CONST_REAL(c3),
    CHF_CONST_REAL(eps),
    CHF_CONST_INT(a_useHOChecks),
    CHF_CONST_INT(a_loBoundHOChecks),
    CHF_CONST_INT(a_hiBoundHOChecks));
}

/*--------------------------------------------------------------------*/
//  PPM Limiter
/**
 *//*-----------------------------------------------------------------*/

void
MOLUtilFunc::PPMLimiter(FArrayBox&           a_dWMinus,
                        FArrayBox&           a_dWPlus,
                        const FArrayBox&     a_W,
                        const FArrayBox&     a_WfromU,
                        const int&           a_numSlopes,
                        const int&           a_dir,
                        const Box&           a_box,
                        const ProblemDomain& a_domain)
{
  CH_TIME("MOLUtilFunc::PPMLimiter");
  // a_dWMinus[i] = WFace[i - e/2] - a_W[i]
  // a_dWPlus[i] = WFace[i + e/2] - a_W[i]
  // where e is unit vector in dimension a_dir.

  // this option added by petermc, 5 Sep 2008
  // Will need to recalculate some D^2's.

  // We calculate d2Wfcf on a_box,
  // and we need a_dWMinus on a_box,
  // and we need a_dWPlus on a_box.
  // d2Wfcf[i] = 6 * (a_dWMinus[i] + a_dWPlus[i])
  //           = 6 * (thisFaceWDir[i-e/2] - a_cellW[i] +
  //                  thisFaceWDir[i+e/2] - a_cellW[i])
  FArrayBox d2Wfcf(a_box, a_numSlopes);
  d2Wfcf.copy(a_dWMinus);
  d2Wfcf.plus(a_dWPlus, 0, 0, a_numSlopes);
  d2Wfcf *= 6.;

  // petermc, 21 Sep 2010:
  // In order to get a_dWMinus and a_dWPlus on a_box,
  // we need d2W on a_box grown by 3 in a_dir directions.
  Box box3 = a_box;
  box3.grow(a_dir, 3);

  Box loBox, hiBox, centerBox, entireBox;
  int hasLo, hasHi;
  loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
             box3, a_domain, a_dir);

  FArrayBox d2W(entireBox, a_numSlopes);
  // On centerBox, use 3-point stencil for d2W;
  // on loBox and hiBox, copy result from neighbor.
  // In the case of centerBox == entireBox == box1,
  // where box1 is a_box grown by 1 in dimension a_dir:
  // We calculate d2W on a_box grown by 2 (was 1) in dimension a_dir,
  // and we need a_W on a_box grown by 3 (was 2) in dimension a_dir.

  // petermc, 21 Sep 2010, changed layer of d2W from 2 to 1,
  // and of a_W from 3 to 2.
  CH_assert(a_WfromU.box().contains(entireBox));
  FORT_GETSECONDDIFF(CHF_FRA(d2W),
                     CHF_CONST_FRA(a_WfromU),
                     CHF_CONST_INT(a_numSlopes),
                     CHF_CONST_INT(a_dir),
                     CHF_BOX(loBox),
                     CHF_CONST_INT(hasLo),
                     CHF_BOX(hiBox),
                     CHF_CONST_INT(hasHi),
                     CHF_BOX(centerBox));

  Box box1 = a_box;
  box1.grow(a_dir, 1);
  Box nextLoBox, nextHiBox, innerCenterBox;
  loHiCenter5(loBox, nextLoBox, hasLo,
              hiBox, nextHiBox, hasHi,
              centerBox, innerCenterBox, entireBox,
              box1, a_domain, a_dir);
  CH_assert(entireBox == a_box);
      
  Real limitC = 1.25;
  Real eps = 1.0e-12;
  Real c3 = 0.1;
  const int useHOChecks = 1;
  const int useLoBoundHOChecks = 1; // Always use 3rd-derivative checks
  const int useHiBoundHOChecks = 1; // near hi and lo boundaries
  FORT_CHECKCUBICLIMITERF(// <W>_(i-1/2) - <W>_i, on a_box
    CHF_FRA(a_dWMinus),
    // <W>_(i+1/2) - <W>_i, on a_box
    CHF_FRA(a_dWPlus),
    CHF_CONST_FRA(a_W),
    CHF_CONST_FRA(d2W),
    CHF_CONST_FRA(d2Wfcf),
    CHF_CONST_INT(a_numSlopes),
    CHF_CONST_INT(a_dir),
    CHF_BOX(loBox),
    CHF_BOX(nextLoBox),
    CHF_CONST_INT(hasLo),
    CHF_BOX(hiBox),
    CHF_BOX(nextHiBox),
    CHF_CONST_INT(hasHi),
    CHF_BOX(innerCenterBox),
    CHF_CONST_REAL(limitC),
    CHF_CONST_REAL(c3),
    CHF_CONST_REAL(eps),
    CHF_CONST_INT(useHOChecks),
    CHF_CONST_INT(useLoBoundHOChecks),
    CHF_CONST_INT(useHiBoundHOChecks));
}

/*--------------------------------------------------------------------*/
//  Convolve or deconvolve on a cell
/** a_avgFab += a_sign * laplacian(a_cellFab) * m_dx^2 / 24, on
 *  a_box
 *  \param[in]  a_sign   1 : convolve
 *                      -1 : deconvolve
 *//*-----------------------------------------------------------------*/

void
MOLUtilFunc::deconvolve(FArrayBox&           a_avgFab,
                        const FArrayBox&     a_cellFab,
                        const Box&           a_box,
                        const ProblemDomain& a_domain,
                        int                  a_sign)
{
  CH_TIME("MOLUtilFunc::deconvolve");
  int numSlopes = a_cellFab.nComp();
  FArrayBox d2Fab(a_box, numSlopes);
  Real scale = (a_sign * 1.) /24.;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      Box bx1(a_box);
      bx1.grow(idir, 1);
      Box loBox, hiBox, centerBox, entireBox;
      int hasLo, hasHi;
      // Generate the domain boundary boxes, loBox and hiBox, if there are
      // domain boundaries there
      loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                 bx1, a_domain, idir);

      // d2Fab[i] = a_cellFab[i-e] - 2*a_cellFab[i] + a_cellFab[i+e]
      FORT_GETSECONDDIFF(CHF_FRA(d2Fab),
                         CHF_CONST_FRA(a_cellFab),
                         CHF_CONST_INT(numSlopes),
                         CHF_CONST_INT(idir),
                         CHF_BOX(loBox),
                         CHF_CONST_INT(hasLo),
                         CHF_BOX(hiBox),
                         CHF_CONST_INT(hasHi),
                         CHF_BOX(centerBox));
      // We have started with a_avgFab set to a_cellFab.
      a_avgFab.plus(d2Fab, scale);
    }
}

/*--------------------------------------------------------------------*/
//  Convolve or deconvolve on a face
/** a_avgFlux += a_sign * laplacian(a_cellFlux) * m_dx^2 / 24, on
 *  faces of a_box
 *  \param[in]  a_sign   1 : convolve
 *                      -1 : deconvolve
 *//*-----------------------------------------------------------------*/

void
MOLUtilFunc::deconvolveFace(FluxBox&             a_avgFlux,
                            const FluxBox&       a_cellFlux,
                            const Box&           a_box,
                            const ProblemDomain& a_domain,
                            int                  a_sign)
{
  CH_TIME("MOLUtilFunc::deconvolveFace");
  int numSlopes = a_cellFlux.nComp();
  Real scale = (a_sign * 1.) /24.;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      FArrayBox& avgFabDir = a_avgFlux[idir];
      const FArrayBox& cellFabDir = a_cellFlux[idir];
      for (int tanDir = 0; tanDir < SpaceDim; tanDir++)
        {
          if (tanDir != idir)
            {
              Box bx1(a_box);
              bx1.grow(tanDir, 1);
              Box loBox, hiBox, centerBox, entireBox;
              int hasLo, hasHi;
              // Generate the domain boundary boxes, loBox and hiBox,
              // if there are domain boundaries there.
              // Don't use loHiCenterFace, because that's for 2-point stencils.
              loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                         bx1, a_domain, tanDir);
              // We'll set d2Fab on idir-faces (NOT tanDir-faces).
              centerBox.surroundingNodes(idir);
              if (hasLo) loBox.surroundingNodes(idir); // lowest layer
              if (hasHi) hiBox.surroundingNodes(idir); // highest layer

              Box bxFaces = a_box; // from loHiCenterFace
              bxFaces.surroundingNodes(idir);
              FArrayBox d2Fab(bxFaces, numSlopes);
              // for each face i in direction idir:
              // d2Fab[i] = a_cellFab[i-e] - 2*a_cellFab[i] + a_cellFab[i+e]
              FORT_GETSECONDDIFF(CHF_FRA(d2Fab),
                                 CHF_CONST_FRA(cellFabDir),
                                 CHF_CONST_INT(numSlopes),
                                 CHF_CONST_INT(tanDir),
                                 CHF_BOX(loBox),
                                 CHF_CONST_INT(hasLo),
                                 CHF_BOX(hiBox),
                                 CHF_CONST_INT(hasHi),
                                 CHF_BOX(centerBox));
              // We have started with avgFabDir set to cellFabDir.
              avgFabDir.plus(d2Fab, scale);
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Convolve or deconvolve on a face described by an FarrayBox
/** a_faceFab += a_sign * laplacian(a_faceExtFab) * m_dx^2 / 24, on
 *  faces of a_box, which is a_dir-face-centered
 *  Convert a_faceFab from face-centered to face-averaged (a_sign = 1)
 *  or vice versa (a_sign = -1) on a_dir-face-centered a_box
 *  by using the second derivative from a_faceExtFab, which will 
 *  typically extend farther than a_faceFab does.
 *  \param[in]  a_sign   1 : convolve
 *                      -1 : deconvolve
 *//*-----------------------------------------------------------------*/

void
MOLUtilFunc::deconvolveCenterFace(FArrayBox&           a_faceFab,
                                  const FArrayBox&     a_faceExtFab,
                                  const Box&           a_box,
                                  const ProblemDomain& a_domain,
                                  int                  a_dir,
                                  int                  a_sign)
{
  CH_TIME("MOLUtilFunc::deconvolveCenterFace");
  int nComp = a_faceFab.nComp();
  Real scale = (a_sign * 1.) /24.;
  for (int tanDir = 0; tanDir < SpaceDim; tanDir++)
    {
      if (tanDir != a_dir)
        {
          // We have data on bx1.
          Box bx1(a_box);
          bx1.grow(tanDir, 1);
          Box loBox, hiBox, centerBox, entireBox;
          int hasLo, hasHi;
          // Generate the domain boundary boxes, loBox and hiBox,
          // if there are domain boundaries there.
          // Don't use loHiCenterFace, because that's for 2-point stencils.
          loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                     bx1, a_domain, tanDir);
          // Note that we might not even need to fill data on loBox and hiBox.

          // We'll set d2Fab on a_dir-faces (NOT tanDir-faces).
          Box bxFaces = a_box; // from loHiCenterFace
          FArrayBox d2Fab(bxFaces, nComp);
          if (hasLo)
            {
              loBox &= bxFaces;
            }
          if (hasHi)
            {
              hiBox &= bxFaces;
            }
          // for each face i in direction a_dir:
          // d2Fab[i] =
          //   a_faceExtFab[i-e] - 2*a_faceExtFab[i] + a_faceExtFab[i+e]
          // where e is unit vector in direction tanDir.
          FORT_GETSECONDDIFF(CHF_FRA(d2Fab),
                             CHF_CONST_FRA(a_faceExtFab),
                             CHF_CONST_INT(nComp),
                             CHF_CONST_INT(tanDir),
                             CHF_BOX(loBox),
                             CHF_CONST_INT(hasLo),
                             CHF_BOX(hiBox),
                             CHF_CONST_INT(hasHi),
                             CHF_BOX(centerBox));
          // We have started with avgFabDir set to cellFabDir.
          a_faceFab.plus(d2Fab, scale);
        }
    }
}

#include "NamespaceFooter.H"
