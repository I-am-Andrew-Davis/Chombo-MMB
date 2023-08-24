#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cstdio>
#include <string>
using std::string;

#include "LoHiSide.H"

#include "MOLUtilities.H"
#include "MOLUtilFunc.H"
#include "LoHiCenter.H"
#include "MOLPhysics.H"
#include "GodunovUtilitiesF_F.H"
#include "DebugOut.H"

#include "NamespaceHeader.H"

//////////////////////////////////////////////////////////////////////////////
// Flag everything as not defined or set
MOLUtilities::MOLUtilities()
{
  m_highOrderLimiter = false;
  m_isDefined = false;
}

//////////////////////////////////////////////////////////////////////////////
MOLUtilities::~MOLUtilities()
{
}

//////////////////////////////////////////////////////////////////////////////
// Define this object and the boundary condition object
void MOLUtilities::define(const ProblemDomain& a_domain,
                              const Real&          a_dx)
{
  // Store the domain and grid spacing
  m_domain    = a_domain;
  m_dx        = a_dx;
  m_isDefined = true;
}

//////////////////////////////////////////////////////////////////////////////
// Compute the flattening coefficients from the primitive variables
void MOLUtilities::computeFlattening(FArrayBox&       a_flattening,
                                         const FArrayBox& a_W,
                                         const Interval&  a_velInt,
                                         const int&       a_presInd,
                                         const Real&      a_smallPres,
                                         const int&       a_bulkModulusInd,
                                         const Box&       a_box)
{
  CH_TIME("MOLUtilities::computeFlattening");
  CH_assert(m_isDefined);
  MOLUtilFunc::computeFlattening(a_flattening,
                                 a_W,
                                 a_velInt,
                                 a_presInd,
                                 a_smallPres,
                                 a_bulkModulusInd,
                                 a_box,
                                 m_domain);
}

//////////////////////////////////////////////////////////////////////////////
void MOLUtilities::applyFlattening(      FArrayBox& a_dW,
                                       const FArrayBox& a_flat,
                                       const Box&       a_box)
{
  CH_TIME("MOLUtilities::applyFlattening");
  MOLUtilFunc::applyFlattening(a_dW,
                               a_flat,
                               a_box);
}

//////////////////////////////////////////////////////////////////////////////
void MOLUtilities::vanLeerSlopes(FArrayBox&       a_dW,
                                     const FArrayBox& a_W,
                                     const int&       a_numSlopes,
                                     const bool&      a_useLimiting,
                                     const int&       a_dir,
                                     const Box&       a_box)
{
  CH_TIME("MOLUtilities::vanLeerSlopes");
  // A box one larger (in direction "a_dir") than the final result box
  // 2 Sep 2008:  For vanLeerSlopesExtPreserving, expand by 2?  17 Sep 2008
  Box box1 = a_box;
  int ghostbox1 = 1; // FIX, 19 sep 2008 (m_highOrderLimiter) ? 2 : 1;
  // int ghostbox1 = (m_highOrderLimiter) ? 2 : 1;
  box1.grow(a_dir, ghostbox1);

  // Compute where centered differences can be used and where one-sided
  // differences need to be used.
  Box loBox,hiBox,centerBox,entireBox;
  int hasLo,hasHi;

  loHiCenter(loBox,hasLo,hiBox,hasHi,centerBox,entireBox,
             box1,m_domain,a_dir);

  // Compute 2nd order slopes - including one-sided differences
  FArrayBox dWMinus(entireBox,a_numSlopes);
  FArrayBox dWPlus (entireBox,a_numSlopes);

  // We calculate a_dW and dWMinus and dWPlus on centerBox,
  // and we need a_W on centerBox grown by 1 in direction a_dir.
  slopes(a_dW, dWMinus, dWPlus, a_W, a_numSlopes, a_dir,
         loBox, hiBox, centerBox, entireBox, hasLo, hasHi);

  // Apply the slope limiter if requested
  if (a_useLimiting)
  {
    // Apply slopeLimiter only on centerBox; elsewhere, a_dW is unchanged.

    // 2 Sep 2008:  replace slopeLimiter with slopeLimiterExtPreserving

    if (m_highOrderLimiter)
      {
        // Modifies a_dW on centerBox,
        // and needs dWMinus on centerBox including shift down by 1 in direction a_dir,
        // and needs dWPlus on centerBox including shift up by 1 in direction a_dir.
        slopeLimiterExtPreserving(a_dW, dWMinus, dWPlus, a_numSlopes, centerBox, a_dir);
      }
    else
      {
        slopeLimiter(a_dW, dWMinus, dWPlus, a_numSlopes, centerBox);
      }
  }
}

//////////////////////////////////////////////////////////////////////////////
void MOLUtilities::fourthOrderSlopes(FArrayBox&       a_dW,
                                         const FArrayBox& a_W,
                                         const FArrayBox& a_dWvL,
                                         const int&       a_numSlopes,
                                         const int&       a_dir,
                                         const Box&       a_box)
{
  CH_TIME("MOLUtilities::fourthOrderSlopes");
  // Number of slopes to compute
  int numSlope = a_numSlopes;

  CH_assert(a_dW.nComp() == numSlope);
  CH_assert(a_W.nComp() >= numSlope);

  // A box one larger (in direction "a_dir") than the final result box
  Box box1 = a_box;
  box1.grow(a_dir,1);

  // Compute where centered differences can be used and where one sided
  // differences need to be used.
  Box loBox,hiBox,centerBox,entireBox;
  int hasLo,hasHi;

  loHiCenter(loBox,hasLo,hiBox,hasHi,centerBox,entireBox,
             box1,m_domain,a_dir);

  CH_assert(a_dW.box().contains(entireBox));

  FORT_FOURTHSLOPEDIFFSF(CHF_FRA(a_dW),
                         CHF_CONST_FRA(a_W),
                         CHF_CONST_FRA(a_dWvL),
                         CHF_CONST_INT(numSlope),
                         CHF_CONST_INT(a_dir),
                         CHF_BOX(loBox),
                         CHF_CONST_INT(hasLo),
                         CHF_BOX(hiBox),
                         CHF_CONST_INT(hasHi),
                         CHF_BOX(centerBox));
}

//////////////////////////////////////////////////////////////////////////////
void MOLUtilities::oneSidedDifferences(FArrayBox&       a_dWMinus,
                                           FArrayBox&       a_dWPlus,
                                           const FArrayBox& a_W,
                                           const int&       a_dir,
                                           const Box&       a_box)
{
  CH_TIME("MOLUtilities::oneSidedDifferences");
  int numSlopes = a_dWMinus.nComp();

  // A box one larger (in direction "a_dir") than the final result box
  Box box1 = a_box;
  box1.grow(a_dir,1);

  // Compute where centered differences can be used and where one sided
  // differences need to be used.
  Box loBox,hiBox,centerBox,entireBox;
  int hasLo,hasHi;

  loHiCenter(loBox,hasLo,hiBox,hasHi,centerBox,entireBox,
             box1,m_domain,a_dir);

  // Compute 2nd order slopes - including one sided differences
  FArrayBox deltaWC(entireBox, numSlopes);

  slopes(deltaWC, a_dWMinus, a_dWPlus, a_W, numSlopes, a_dir,
         loBox, hiBox, centerBox, entireBox, hasLo, hasHi);
}

//////////////////////////////////////////////////////////////////////////////
void MOLUtilities::slopes(FArrayBox&       a_dWCent,
                              FArrayBox&       a_dWMinus,
                              FArrayBox&       a_dWPlus,
                              const FArrayBox& a_W,
                              const int&       a_numSlopes,
                              const int&       a_dir,
                              const Box&       a_loBox,
                              const Box&       a_hiBox,
                              const Box&       a_centerBox,
                              const Box&       a_entireBox,
                              const int&       a_hasLo,
                              const int&       a_hasHi)

{
  CH_TIME("MOLUtilities::slopes");
  CH_assert(a_dWCent .nComp() == a_numSlopes);
  CH_assert(a_dWMinus.nComp() == a_numSlopes);
  CH_assert(a_dWPlus .nComp() == a_numSlopes);
  CH_assert(a_W.nComp() >= a_numSlopes);

  CH_assert(a_dWCent .box().contains(a_entireBox));
  CH_assert(a_dWMinus.box().contains(a_entireBox));
  CH_assert(a_dWPlus .box().contains(a_entireBox));
  CH_assert(a_W.box().contains( a_entireBox));

  CH_assert((a_dir >= 0 )&& (a_dir < SpaceDim));

  FORT_SECONDSLOPEDIFFSF(CHF_FRA(a_dWCent),
                         CHF_FRA(a_dWMinus),
                         CHF_FRA(a_dWPlus),
                         CHF_CONST_FRA(a_W),
                         CHF_CONST_INT(a_numSlopes),
                         CHF_CONST_INT(a_dir),
                         CHF_BOX(a_loBox),
                         CHF_CONST_INT(a_hasLo),
                         CHF_BOX(a_hiBox),
                         CHF_CONST_INT(a_hasHi),
                         CHF_BOX(a_centerBox));
}

//////////////////////////////////////////////////////////////////////////////
// Apply a van Leer limiter directly to the slopes.
void MOLUtilities::slopeLimiter(FArrayBox&       a_dW,
                                    const FArrayBox& a_dWLeft,
                                    const FArrayBox& a_dWRigh,
                                    const int&       a_numSlopes,
                                    const Box&       a_box)
{
  CH_TIME("MOLUtilities::slopeLimiter");
  CH_assert(m_isDefined);
  CH_assert(a_dW.nComp()     == a_numSlopes);
  CH_assert(a_dWLeft.nComp() == a_numSlopes);
  CH_assert(a_dWRigh.nComp() == a_numSlopes);
  CH_assert(a_dW.box().contains(a_box));
  CH_assert(a_dWLeft.box().contains(a_box));
  CH_assert(a_dWRigh.box().contains(a_box));

  FORT_VANLEERLIMITERF(CHF_FRA(a_dW),
                       CHF_CONST_FRA(a_dWLeft),
                       CHF_CONST_FRA(a_dWRigh),
                       CHF_CONST_INT(a_numSlopes),
                       CHF_BOX(a_box));
}

//////////////////////////////////////////////////////////////////////////////
// Apply an extremum-preserving van Leer limiter directly to the slopes.
void MOLUtilities::slopeLimiterExtPreserving(FArrayBox&       a_dW,
                                                 const FArrayBox& a_dWLeft,
                                                 const FArrayBox& a_dWRigh,
                                                 const int&       a_numSlopes,
                                                 const Box&       a_box,
                                                 const int&       a_dir)
{
  CH_TIME("MOLUtilities::slopeLimiterExtPreserving");
  CH_assert(m_isDefined);
  CH_assert(a_dW.nComp()     == a_numSlopes);
  CH_assert(a_dWLeft.nComp() == a_numSlopes);
  CH_assert(a_dWRigh.nComp() == a_numSlopes);
  CH_assert(a_dW.box().contains(a_box));
  CH_assert(a_dWLeft.box().contains(a_box));
  CH_assert(a_dWRigh.box().contains(a_box));

  // Compute where centered differences can be used and where one-sided
  // differences need to be used.
  Box loBox, hiBox, centerBox, entireBox;
  int hasLo, hasHi;

  loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
             a_box, m_domain, a_dir);

  if (hasLo)
    {
      FORT_VANLEERLIMITERF(CHF_FRA(a_dW),
                           CHF_CONST_FRA(a_dWLeft),
                           CHF_CONST_FRA(a_dWRigh),
                           CHF_CONST_INT(a_numSlopes),
                           CHF_BOX(loBox));
    }
  if (hasHi)
    {
      FORT_VANLEERLIMITERF(CHF_FRA(a_dW),
                           CHF_CONST_FRA(a_dWLeft),
                           CHF_CONST_FRA(a_dWRigh),
                           CHF_CONST_INT(a_numSlopes),
                           CHF_BOX(hiBox));
    }
  if (!centerBox.isEmpty())
    {
      // Modifies a_dW on centerBox,
      // and needs a_dWLeft on centerBox including shift down by 1 in direction a_dir,
      // and needs a_dWRigh on centerBox including shift up by 1 in direction a_dir.
      FORT_EXTPRESERVINGVANLEERLIMITERF(CHF_FRA(a_dW),
                                        CHF_CONST_FRA(a_dWLeft),
                                        CHF_CONST_FRA(a_dWRigh),
                                        CHF_CONST_INT(a_numSlopes),
                                        CHF_CONST_INT(a_dir),
                                        CHF_BOX(centerBox));
    }
  // dummy statement in order to get around gdb bug
  int dummy_unused = 0; dummy_unused = 0;
}

//////////////////////////////////////////////////////////////////////////////
void MOLUtilities::PPMFaceValues(FArrayBox&            a_WFace,
                                     const FArrayBox&      a_W,
                                     const int&            a_numSlopes,
                                     const bool&           a_useLimiting,
                                     const int&            a_dir,
                                     const Box&            a_box,
                                     const Real&           a_time,
                                     const MOLPhysics* a_physPtr)
{
  CH_TIME("MOLUtilities::PPMFaceValues");
  // a_box is the face-centered box on which a_WFace is computed.
  CH_assert(a_WFace.box().contains(a_box));

  if (m_highOrderLimiter)
    {
      MOLUtilFunc::PPMFaceValues(a_WFace,
                                 a_W,
                                 a_numSlopes,
                                 a_dir,
                                 a_box,
                                 m_domain);
    }
  else // !m_highOrderLimiter :  this is the old method
    {
      // A box one larger (in direction "a_dir") than the final result box
      // petermc, 14 Aug 2009:  first grow(), then enclosedCells(),
      // rather than reverse order, so you don't end up with empty box1cells.
      Box box1cells(a_box);
      int ghostbox1 = 1;
      box1cells.grow(a_dir, ghostbox1);
      box1cells.enclosedCells();

      FArrayBox dW(box1cells, a_numSlopes);
      vanLeerSlopes(dW, a_W, a_numSlopes, a_useLimiting, a_dir, box1cells);

      if (a_physPtr != NULL)
        {
          Real a_time = 0.0;
          a_physPtr->getPhysIBC()->setBdrySlopes(dW, a_W, a_dir, a_time);
        }

      Box loBox,hiBox,centerBox,entireBox;
      int hasLo,hasHi;

      loHiCenterFace(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                     box1cells, m_domain, a_dir);

      // a_Wface[i-e/2] = (a_W[i-e] + dW[i-e]/3)/2 + (a_W[i] - dW[i]/3)/2
      FORT_PPMFACEVALUESF(CHF_FRA(a_WFace),
                          CHF_CONST_FRA(a_W),
                          CHF_CONST_FRA(dW),
                          CHF_CONST_INT(a_numSlopes),
                          CHF_CONST_INT(a_dir),
                          CHF_BOX(loBox),
                          CHF_CONST_INT(hasLo),
                          CHF_BOX(hiBox),
                          CHF_CONST_INT(hasHi),
                          CHF_BOX(centerBox));
    }
}

//////////////////////////////////////////////////////////////////////////////
void MOLUtilities::PPMLimiter(FArrayBox& a_dWMinus,
                              FArrayBox& a_dWPlus,
                              const FArrayBox& a_W,
                              const int& a_numSlopes,
                              const int& a_dir,
                              const Box& a_box,
                              const int  a_useHOChecks,
                              const int  a_loBoundHOChecks,
                              const int  a_hiBoundHOChecks)
{
  CH_TIME("MOLUtilities::PPMLimiter");
  // a_dWMinus[i] = WFace[i - e/2] - a_W[i]
  // a_dWPlus[i] = WFace[i + e/2] - a_W[i]
  // where e is unit vector in dimension a_dir.

  if (m_highOrderLimiter)
    {
      MOLUtilFunc::PPMLimiter(a_dWMinus,
                              a_dWPlus,
                              a_W,
                              a_numSlopes,
                              a_dir,
                              a_box,
                              m_domain,
                              a_useHOChecks,
                              a_loBoundHOChecks,
                              a_hiBoundHOChecks);
    }
  else
    {
      FORT_PPMLIMITERF(// <W>_(i-1/2) - <W>_i, on a_box
                       CHF_FRA(a_dWMinus),
                       // <W>_(i+1/2) - <W>_i, on a_box
                       CHF_FRA(a_dWPlus),
                       CHF_CONST_INT(a_numSlopes),
                       CHF_BOX(a_box));
    }
}

//////////////////////////////////////////////////////////////////////////////
void MOLUtilities::PPMLimiter(FArrayBox& a_dWMinus,
                              FArrayBox& a_dWPlus,
                              const FArrayBox& a_W,
                              const FArrayBox& a_WfromU,
                              const int& a_numSlopes,
                              const int& a_dir,
                              const Box& a_box)
{
  CH_TIME("MOLUtilities::PPMLimiter");
  // a_dWMinus[i] = WFace[i - e/2] - a_W[i]
  // a_dWPlus[i] = WFace[i + e/2] - a_W[i]
  // where e is unit vector in dimension a_dir.

  if (m_highOrderLimiter)
    {
      MOLUtilFunc::PPMLimiter(a_dWMinus,
                              a_dWPlus,
                              a_W,
                              a_WfromU,
                              a_numSlopes,
                              a_dir,
                              a_box,
                              m_domain);
    }
  else
    {
      FORT_PPMLIMITERF(// <W>_(i-1/2) - <W>_i, on a_box
                       CHF_FRA(a_dWMinus),
                       // <W>_(i+1/2) - <W>_i, on a_box
                       CHF_FRA(a_dWPlus),
                       CHF_CONST_INT(a_numSlopes),
                       CHF_BOX(a_box));
    }
}

//////////////////////////////////////////////////////////////////////////////
// Compute a face centered divergence of the velocity
void MOLUtilities::divVel(FArrayBox&       a_divVel,
                              const FArrayBox& a_W,
                              const Interval&  a_velInt,
                              const int&       a_dir,
                              const Box&       a_box)
{
  CH_TIME("MOLUtilities::divVel");
  // Quick description of this function:
  // Let v[d] denote d-th component of velocity (from a_velInt within a_W).
  // Then this function returns, on face-centered a_box,
  // a_divVel[i + e[d]/2] := v[d][i+e[d]] - v[d][i] +
  //                 sum_{d' ~= d} ( Dv[d'][i+e[d]] + Dv[d'][i] ) / 2
  // where Dv[d'][i] = (v[d'][i+e[d']] - v[d'][i-e[d']])/2  in center
  //                or  v[d'][i+e[d']] - v[d'][i]  on low  d'-side
  //                or  v[d'][i] - v[d'][i-e[d']]  on high d'-side
  CH_assert((a_dir >= 0 )&& (a_dir < SpaceDim));
  CH_assert(a_divVel.box().contains(a_box));

  Box dveltanBox = a_box;
  dveltanBox.enclosedCells(a_dir);
  dveltanBox.grow(1);
  dveltanBox &= m_domain;

  // First, we need to calculate the directional derivatives of
  // the tangential components of velocity at the cell-centers.
  // want to make sure numTanComps is > 0 even in 1D, since
  // we pass it in to the divergence (although it isn't used there)
  int numTanComps = std::max(1, SpaceDim-1);
  FArrayBox dveltan(dveltanBox,numTanComps);

  // Get the interval of the primitive variables corresponding to the velocity
  int v0index = a_velInt.begin();

  // Go through the tangential directions
  for (int i = 0, dir = 0; dir < SpaceDim; ++dir)
  {
    if (dir != a_dir)
    {
      // This velocity component is tangential.  Build the box in which
      // d(v[dir])/d(x[dir]) is to be computed.
      Box primBox = a_box;
      primBox.enclosedCells(a_dir);
      primBox.grow(dir,1).grow(a_dir,1);

      // Compute where centered differences can be used and where one sided
      // differences need to be used.  The data used for the differences is
      // defined on "dveltanBox"
      Box gradBox,hiBox,loBox,centerBox;
      int hasLo,hasHi;

      loHiCenter(loBox,hasLo,hiBox,hasHi,centerBox,gradBox,
                 primBox,m_domain,dir);

      // Compute d(v[dir])/d(x[dir]).
      FORT_GETGRADF(CHF_FRA1(dveltan,i),
                    CHF_CONST_FRA1(a_W,v0index+dir),
                    CHF_CONST_INT(dir),
                    CHF_BOX(loBox),
                    CHF_CONST_INT(hasLo),
                    CHF_BOX(hiBox),
                    CHF_CONST_INT(hasHi),
                    CHF_BOX(centerBox));
      i++;
    }
  }

  // Now, we can calculate the divergence of the normal velocity
  // at the center normal-direction faces. To do this, we determine
  // the faces at which we have sufficient data to compute centered
  // estimates of h*(div(u)). At the remaining faces. i.e. those
  // corresponding to the physical boundaries, we use zeroth-order
  // extrapolation.

  Box divBox = a_box;
  divBox.enclosedCells(a_dir);
  divBox.grow(a_dir,1);

  Box loBox,hiBox,centerBox,entireBox;
  int hasLo,hasHi;

  loHiCenterFace(loBox,hasLo,hiBox,hasHi,centerBox,entireBox,
                 divBox,m_domain,a_dir);

  // All of the boxes computed above are shifted so as to be cell-centered,
  // with the index of the cell center being identified with the low face.
  // We then shift a_divVel to be compatible with that convention on input,
  // and undo the shift on output.  Basically, everything is made compatible
  // a_W (which is cell-centered).

  loBox.shiftHalf(a_dir,1);
  centerBox.shiftHalf(a_dir,1);
  hiBox.shiftHalf(a_dir,1);

  a_divVel.shiftHalf(a_dir,1);

  FORT_DIVUEDGEF(CHF_FRA1(a_divVel,0),
                 CHF_CONST_FRA1(a_W,v0index+a_dir),
                 CHF_CONST_FRA(dveltan),
                 CHF_CONST_INT(a_dir),
                 CHF_BOX(loBox),
                 CHF_CONST_INT(hasLo),
                 CHF_BOX(hiBox),
                 CHF_CONST_INT(hasHi),
                 CHF_BOX(centerBox));

  a_divVel.shiftHalf(a_dir,-1);
}

//////////////////////////////////////////////////////////////////////////////
// Increment fluxes with artificial viscosity.
void MOLUtilities::artificialViscosity(FArrayBox&       a_F,
                                           const FArrayBox& a_U,
                                           const FArrayBox& a_divVel,
                                           const Real&      a_coeff,
                                           const int&       a_dir,
                                           const Box&       a_box)
{
  CH_TIME("MOLUtilities::artificialViscosity");
  CH_assert((a_dir >= 0 )&& (a_dir < SpaceDim));
  CH_assert(a_F.box().contains(a_box));
  CH_assert(a_divVel.box().contains(a_box));
  CH_assert(a_coeff >= 0.0);

  // WHAT GETS CALCULATED:
  // with d = a_dir,
  // a_F[i+e[d]/2] += a_coeff*min{a_divVel[i+e[d]/2], 0}*(a_U[i+e[d]]-a_U[i]) .
  //
  // That is, if a_divVel[i] >= 0 then no change.
  // Otherwise,
  // a_F[i+e[d]/2] += a_coeff * a_divVel[i+e[d]/2] * (a_U[i+e[d]] - a_U[i]) .
  FORT_ARTVISCF(CHF_FRA(a_F),
                CHF_CONST_FRA(a_U),
                CHF_CONST_FRA1(a_divVel,0),
                CHF_CONST_REAL(a_coeff),
                CHF_CONST_INT(a_dir),
                CHF_BOX(a_box));
}

//////////////////////////////////////////////////////////////////////////////
void MOLUtilities::highOrderLimiter(bool a_highOrderLimiter)
{
  m_highOrderLimiter = a_highOrderLimiter;
}

//////////////////////////////////////////////////////////////////////////////

void MOLUtilities::deconvolve(FArrayBox&        a_avgFab,
                                  const FArrayBox&  a_cellFab,
                                  int               a_sign)
{
  CH_TIME("MOLUtilities::deconvolve");
  const Box& bx = a_cellFab.box();
  int numSlopes = a_cellFab.nComp();
  FArrayBox d2Fab(bx, numSlopes);
  Real scale = (a_sign * 1.) /24.;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      Box loBox, hiBox, centerBox, entireBox;
      int hasLo, hasHi;
      // Generate the domain boundary boxes, loBox and hiBox, if there are
      // domain boundaries there
      loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                 bx, m_domain, idir);

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
  // dummy statement in order to get around gdb bug
  int dummy_unused = 0; dummy_unused = 0;
}

//////////////////////////////////////////////////////////////////////////////

void MOLUtilities::deconvolve(FArrayBox&       a_avgFab,
                              const FArrayBox& a_cellFab,
                              const Box&       a_box,
                              int              a_sign)
{
  CH_TIME("MOLUtilities::deconvolve");
  MOLUtilFunc::deconvolve(a_avgFab,
                          a_cellFab,
                          a_box,
                          m_domain,
                          a_sign);
}

//////////////////////////////////////////////////////////////////////////////

void MOLUtilities::deconvolveFace(FluxBox&       a_avgFlux,
                                  const FluxBox& a_cellFlux,
                                  int            a_sign)
{
  CH_TIME("MOLUtilities::deconvolveFace");
  int numSlopes = a_cellFlux.nComp();
  Real scale = (a_sign * 1.) /24.;
  const Box& bxCells = a_cellFlux.box();
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      FArrayBox& avgFabDir = a_avgFlux[idir];
      const FArrayBox& cellFabDir = a_cellFlux[idir];
      for (int tanDir = 0; tanDir < SpaceDim; tanDir++)
        {
          if (tanDir != idir)
            {
              Box loBox, hiBox, centerBox, entireBox;
              int hasLo, hasHi;
              // Generate the domain boundary boxes, loBox and hiBox,
              // if there are domain boundaries there.
              // Don't use loHiCenterFace, because that's for 2-point stencils.
              loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                         bxCells, m_domain, tanDir);
              centerBox.surroundingNodes(idir);
              if (hasLo) loBox.surroundingNodes(idir); // lowest layer
              if (hasHi) hiBox.surroundingNodes(idir); // highest layer

              Box bxFaces = bxCells; // from loHiCenterFace
              bxFaces.surroundingNodes(idir);
              FArrayBox d2Fab(bxFaces, numSlopes);
              d2Fab.setVal(0.0); // FIXME: Necessary for now. -JNJ
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
  // dummy statement in order to get around gdb bug
  int dummy_unused = 0; dummy_unused = 0;
}

//////////////////////////////////////////////////////////////////////////////

void MOLUtilities::deconvolveFace(FluxBox&       a_avgFlux,
                                  const FluxBox& a_cellFlux,
                                  const Box&     a_box,
                                  int            a_sign)
{
  CH_TIME("MOLUtilities::deconvolveFace");
  MOLUtilFunc::deconvolveFace(a_avgFlux,
                              a_cellFlux,
                              a_box,
                              m_domain,
                              a_sign);
}

//////////////////////////////////////////////////////////////////////////////
// Compute a face-centered nonlinear function of the divergence suitable for use
// as an artificial viscosity coefficient for a fourth-order method.
// F_visc = -dx*K*dU/dx
// K = max(-dx*div(u)*min(1,h^2*|div(u)/(c^2)/a_M0sq|),0)

void MOLUtilities::divVelHO(FArrayBox&       a_divVel,
                            const FArrayBox& a_W,
                            const int&       a_dir,
                            const Box&       a_box,
                            MOLPhysics*      a_physPtr)
{
  CH_TIME("MOLUtilities::divVelHO");
  CH_assert((a_dir >= 0 )&& (a_dir < SpaceDim));
  CH_assert(a_divVel.box().contains(a_box));
  CH_assert(a_physPtr->fourthOrderArtificialViscosityIsDefined() == true);
  // Compute face-centered divergence.
  Interval velInt = a_physPtr->velocityInterval();
  // a_W on valid cells + 1 layer of ghosts;
  // a_box = a_dir-faces of valid cells
  divVel(a_divVel, a_W, velInt, a_dir, a_box);

  // box1cells = valid cells + 1 layer of ghosts;
  // boxcsq = box1cells intersected with domain.
  Box box1cells = a_box;
  box1cells.enclosedCells();
  int ghostbox1 = 1;
  box1cells.grow(ghostbox1);
  Box boxcsq = box1cells & m_domain;

  // Compute cell-centered (bulk modulus)/(density).
  int bulkIndex = a_physPtr->bulkModulusIndex();
  int densityIndex = a_physPtr->densityIndex();
  Real M0sq = a_physPtr->getFourthOrderArtificialViscosityParameter();
  FArrayBox csq1(boxcsq, 1);
  FArrayBox csq2(boxcsq, 1);
  // csq1 = a_W[bulkIndex] / a_W[densityIndex] on boxcsq,
  // which is valid cells + 1 ghost layer, intersected with domain.
  csq1.setVal(0.);
  csq2.setVal(0.);
  csq1.copy(a_W, boxcsq, bulkIndex, boxcsq, 0, 1);
  csq1.divide(a_W, boxcsq, boxcsq, densityIndex, 0, 1);
  Box hiBox,loBox,centerBox,entireBox;
  int hasLo,hasHi;
  FArrayBox* csqin_ptr = &csq1;
  FArrayBox* csqout_ptr = &csq2;

  // Compute min of csq in transverse direction(s).
  for (int dir = 0; dir < SpaceDim; ++dir)
    {
      if (dir != a_dir)
        {
          // WHAT GETS CALCULATED:
          // Pass on transverse direction dir
          // takes *csqin_ptr and stores result in *csqout_ptr.
          // csqout[i] =
          // min{csqin[i-e[dir]], csqin[i], csqin[i+e[dir]]}  in center;
          // min{csqin[i-e[dir]], csqin[i]}                   at high end;
          //                  min{csqin[i], csqin[i+e[dir]]}  at low end.
          FArrayBox& csqin = *csqin_ptr;
          FArrayBox& csqout = *csqout_ptr;
          loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                     box1cells, m_domain,dir);

          FORT_MIN3PTSF(CHF_FRA1(csqout,0),
                        CHF_CONST_FRA1(csqin,0),
                        CHF_CONST_INT(dir),
                        CHF_BOX(loBox),
                        CHF_CONST_INT(hasLo),
                        CHF_BOX(hiBox),
                        CHF_CONST_INT(hasHi),
                        CHF_BOX(centerBox));
          box1cells.grow(dir,-1);
          FArrayBox* csqtmp = csqin_ptr;
          csqin_ptr = csqout_ptr;
          csqout_ptr = csqtmp;
        }
    }
  // divBox = valid cells + 1 ghost layer in direction a_dir
  Box divBox = a_box;
  divBox.enclosedCells(a_dir);
  divBox.grow(a_dir, 1);
  loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
             divBox, m_domain, a_dir);

  // All of the boxes computed above are shifted so as to be cell-centered,
  // with the index of the cell center being identified with the low face.
  // We then shift a_divVel to be compatible with that convention on input,
  // and undo the shift on output.  Basically, everything is made compatible
  // with a_W (which is cell-centered).

  // WHAT GETS CALCULATED:
  // Here we multiply face-centered a_divVel by min(a_divVel**2/csqmin/M0sq, 1)
  // where csqmin is min of csq on the two neighboring cells
  // (or csq on the one neighboring cell if on a boundary).
  // That is, if a_divVel**2 >= csqmin*M0sq, then no change in a_divVel.
  // Otherwise, a_divVel *= a_divVel**2/csqmin/M0sq.

  // shift a_divVel from a_dir-faces of valid cells
  // to valid cells + 1 ghost layer in Hi direction a_dir
  a_divVel.shiftHalf(a_dir, 1);
  // if (!hasLo && !hasHi)
  //   centerBox = valid cells - 1 layer in Lo direction a_dir
  // centerBox.growDir(a_dir, Side::Lo, -1);
  // petermc, 10 Jul 2009, changed to this
  centerBox.growDir(a_dir, Side::Hi, +1);
  // csq on valid cells + 1 ghost layer, intersected with domain.
  // csq contains pressure / density, after passes along transverse
  // directions taking minimum with neighbors.
  FArrayBox& csq = *csqin_ptr;
  FORT_HODIVCOEF(CHF_FRA1(a_divVel,0), // set this on centerBox, loBox, hiBox
                 CHF_CONST_FRA1(csq,0),
                 CHF_CONST_INT(a_dir),
                 CHF_CONST_REAL(M0sq),
                 CHF_BOX(loBox),
                 CHF_CONST_INT(hasLo),
                 CHF_BOX(hiBox),
                 CHF_CONST_INT(hasHi),
                 CHF_BOX(centerBox));
  a_divVel.shiftHalf(a_dir,-1);

  // if (!hasLo && !hasHi)
  //   a_divVel has been set on a_dir-faces BETWEEN valid cells
}

//////////////////////////////////////////////////////////////////////////////

void MOLUtilities::deconvolveCenterFace(FArrayBox&        a_faceFab,
                                        const FArrayBox&  a_faceExtFab,
                                        const Box&        a_box,
                                        int               a_dir,
                                        int               a_sign)
{
  CH_TIME("MOLUtilities::deconvolveCenterFace");
  // Added this function because I don't have it elsewhere:
  // MOLUtilities::deconvolveFace() is for FluxBoxes,
  // and MOLUtilities::deconvolve() is for cell-centered FABs and takes
  // derivatives in all directions.
  // FourthOrderUtil functions do not take the additional "for grad" argument.
  MOLUtilFunc::deconvolveCenterFace(a_faceFab,
                                    a_faceExtFab,
                                    a_box,
                                    m_domain,
                                    a_dir,
                                    a_sign);
}

#include "NamespaceFooter.H"
