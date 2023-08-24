#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "TimeInterpolatorARK4.H"
#include "AMRIO.H"
#include "NamespaceHeader.H"

/*******************************************************************************
 *
 * \file TimeInterpolatorARK4.cpp
 *
 * \brief Time interpolation for 6-stage ARK4
 *
 ******************************************************************************/

/*==============================================================================
 * Constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
///  Default constructor
/**  The 'define' function must be called before using this class
 *//*-----------------------------------------------------------------*/

TimeInterpolatorARK4::TimeInterpolatorARK4()
{
  m_defined = false;
  resetData();
}

/*--------------------------------------------------------------------*/
///  Destructor
/*--------------------------------------------------------------------*/

TimeInterpolatorARK4::~TimeInterpolatorARK4()
{
}


/*==============================================================================
 * Member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
///  Define function allocates space
/**  \param[in]  a_thisDisjointBoxLayout
 *                      The DBL of this level (the finer level)
 *   \param[in]  a_coarserDisjointBoxLayout
 *                      The DBL at next coarser level
 *   \param[in]  a_domain
 *                      Problem domain on this level
 *   \param[in]  a_refineCoarse
 *                      Refinement ratio between this level and next coarser level
 *   \param[in]  a_numStates
 *                      Number of variables/components
 *   \param[in]  a_ghosts
 *                      Layers of ghost cells to be filled in on the coarsened layout at this level
 *//*-----------------------------------------------------------------*/

void TimeInterpolatorARK4::define(const DisjointBoxLayout&  a_thisDisjointBoxLayout,
                                  const DisjointBoxLayout&  a_coarserDisjointBoxLayout,
                                  const ProblemDomain&      a_domain,
                                  const int&                a_refineCoarse,
                                  const int&                a_numStates,
                                  const int&                a_ghosts)
{

  // Cache some data
  m_refineCoarse = a_refineCoarse;
  m_coarseDomain = coarsen(a_domain, m_refineCoarse);
  m_numStates = a_numStates;
  m_ghosts = a_ghosts;
  m_ghostVect = m_ghosts*IntVect::Unit;

  m_coarseLayout = a_coarserDisjointBoxLayout;
  m_thisCoarsenedLayout = DisjointBoxLayout(); // Prevents crash if define called again
  coarsen(m_thisCoarsenedLayout, a_thisDisjointBoxLayout, m_refineCoarse); // Creates coarsened-fine

  m_rhsCopy.define(m_thisCoarsenedLayout, m_numStates, m_ghostVect);
  m_denseCoeffs.define(m_thisCoarsenedLayout,
                       m_numCoeffs * m_numStates, // Number of comps
                       m_ghostVect);
  m_copier.define(m_coarseLayout, m_thisCoarsenedLayout,
                  m_coarseDomain, m_ghostVect);

  m_defined = true;
}


/*--------------------------------------------------------------------*/
///  Sets the timestep dt for the time domain being interpolated in
/**
 *   \param[in]  a_dt   The timestep
 *//*-----------------------------------------------------------------*/
void TimeInterpolatorARK4::setDt(const Real&  a_dt)
{
  resetData();
  m_dt = a_dt;
  m_gotDt = true;
}


/*--------------------------------------------------------------------*/
///  Saves the initial solution of the coarse grid
/**
 *   \param[in]  a_soln The initial solution
 *//*-----------------------------------------------------------------*/
void TimeInterpolatorARK4::saveInitialSoln(const LevelData<FArrayBox>&   a_soln)
{
  CH_assert(m_defined);
  CH_assert(m_gotDt);
  CH_assert(!m_gotInitialSoln);
  CH_assert(a_soln.nComp() == m_numStates);

  // Set initial coeffs to zero
  DataIterator dit = m_denseCoeffs.dataIterator();
  for(dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& denseFab = m_denseCoeffs[dit];
      denseFab.setVal(0.);
    }

  // a_soln is passed to us on the coarse layout,
  // and m_denseCoeffs is on the coarsened-fine layout
  // So we copy from coarse to coarsened-fine
  // So fill first m_numStates of m_denseCoeffs with a_soln
  const Interval& srcInt = a_soln.interval();
  a_soln.copyTo(srcInt, m_denseCoeffs, srcInt, m_copier);

  m_gotInitialSoln = true;
}

/*--------------------------------------------------------------------*/
///  Save a RHS stage from the time integration
/**
 *   \param[in]  a_rhs  The RHS to save
 *//*-----------------------------------------------------------------*/
void TimeInterpolatorARK4::saveRHS(const LevelData<FArrayBox>&   a_rhs)
{
  bool wrong_saveRHS_called = false;
  // If you get this assert, that means that you called the wrong saveRHS.
  // Call the saveRHS with two arguments
  // It doesn't make sense to have saveRHS with one argument for Additive RK,
  // but I need this function here to block any calls to the parent saveRHS.
  CH_assert(wrong_saveRHS_called);
}

/*--------------------------------------------------------------------*/
///  Save a RHS stage from the time integration
/**
 *   \param[in]  a_ExpRHS  The explicit RHS term to save
 *   \param[in]  a_ImpRHS  The implicit RHS term to save
 *//*-----------------------------------------------------------------*/
void TimeInterpolatorARK4::saveRHS(const LevelData<FArrayBox>& a_ExpRHS,
                                   const LevelData<FArrayBox>& a_ImpRHS)
{
  CH_assert(m_defined);
  CH_assert(m_gotDt);
  CH_assert(m_gotInitialSoln);
  CH_assert(!m_gotFullDensePoly);
  CH_assert(m_countRHS >= 0);
  CH_assert(m_countRHS < m_numStages);
  CH_assert(a_ExpRHS.nComp() == m_numStates);
  CH_assert(a_ImpRHS.nComp() == m_numStates);
  CH_assert(a_ExpRHS.disjointBoxLayout() == a_ImpRHS.disjointBoxLayout());

  Real tol= 1e-12;

  // a_rhs is on the coarse layout
  // m_rhsCopy is on the coarsened-fine layout
  a_ExpRHS.copyTo(m_rhsCopy, m_copier);

  // Since bstar is the same for both the explicit and implicit RHS,
  // we can add them together here
  a_ImpRHS.copyTo(m_rhsCopy, m_copier, LDaddOp<FArrayBox>());

  // Now that we have the RHS on the coarsened-fine layout
  // We can multiply it by the dense output coefficient and
  // add it to the denseCoeff data
  DataIterator dit = m_denseCoeffs.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const FArrayBox& rhsCopyFab = m_rhsCopy[dit];
      FArrayBox& denseFab = m_denseCoeffs[dit];

      for (int ind = 0; ind < 3; ++ind)
        {
          Real bstar = s_bstar[ind][m_countRHS];
          if(Abs(bstar) > tol) // If bstar is not zero
            {
              // Do denseFab += bstar*m_dt*rhsCopyFab
              // We do (ind+1) because ind goes from zero to 3 and we want 1 to 4 here
              // because ind=0 is reserved for the U_n solution
              denseFab.plus(rhsCopyFab,
                            bstar*m_dt, // multiply rhsCopyFab by this
                            0, // start of rhsCopyFab component
                            (ind+1)*m_numStates, // Start of denseFab components
                            m_numStates); // Number of components to copy
            }
        }
    }

  ++m_countRHS;
  if(m_countRHS == m_numStages)
    {
      m_denseCoeffs.exchange();
      m_gotFullDensePoly = true;
    }
}


/*--------------------------------------------------------------------*/
///  Interpolate the solution to some intermediate
/**
 *   \param[out] a_U      Interpolated solution on this level coarsened
 *   \param[in]  a_timeInterpCoeff
 *                        Time interpolation coefficient in rage [0:1]
 *   \param[in]  a_intvl  Interval of a_U to fill in
 *//*-----------------------------------------------------------------*/
void TimeInterpolatorARK4::interpolate(/// 
                                      LevelData<FArrayBox>&   a_U,
                                      /// time interpolation coefficient in range [0:1]
                                      const Real&             a_timeInterpCoeff,
                                      /// interval of a_U to fill in
                                      const Interval&         a_intvl)
{
  // I don't know where this function is used and it doesn't provide complete information
  // needed to interpolate, so kill it for now... 
  // You should use intermediate() instead.
  CH_assert(false);
  // May be able to foward this to intermediate() passing any stage value that is not 0 and
  // any value for dtRatio since we don't use dtRatio.
}


/*--------------------------------------------------------------------*/
///  Interpolate the solution to some intermediate
/**
 *   \param[out] a_U    Intermediate RK4 solution on this level coarsened
 *   \param[in]  a_timeInterpCoeff
 *                        Time interpolation coefficient in rage [0:1]
 *   \param[in]  a_dtRatio
 *                        Ratio of fine to coarse time steps. This is usually 1/nref unless extra subcycling occurs.
 *   \param[in]  a_stage  Which RK4 stage: 0,1,2,3,4,5
 *   \param[in]  a_intvl  Interval of a_U to fill in
 *//*-----------------------------------------------------------------*/
void TimeInterpolatorARK4::intermediate(LevelData<FArrayBox>&   a_U,
                                        const Real&             a_timeInterpCoeff,
                                        const Real              a_dtRatio,
                                        const int&              a_stage,
                                        const Interval&         a_intvl) const
{
  CH_assert(m_defined);
  CH_assert(m_gotFullDensePoly);
  CH_assert(a_U.nComp() == m_numStates);
  CH_assert(a_stage >= 0);
  CH_assert(a_stage < m_numStages);

  (void)a_dtRatio; // Don't need a_dtRatio for 6stage ARK4, so hide the unused warning

  // Okay, now we can interpolate (or "intermediate")...
  // denseCoeffs has the previous time solution U_n in the first m_numStates components
  // The second m_numStates components of denseCoeffs corresponds to the coefficients of a_timeInterpCoeff
  // The third m_numstates components of denseCoeffs corresponds to the coefficients of a_timeInterpCoeff^2 etc.

  // Size of interval to copy
  const int intvLength = a_intvl.size();

  // If first stage of first subcycle, its just the old U_n state
  if(a_timeInterpCoeff == 0. && a_stage == 0)
    {
      DataIterator dit = a_U.dataIterator();
      for(dit.begin(); dit.ok(); ++dit)
        {
          a_U[dit].copy(m_denseCoeffs[dit],
                        0, // Start of denseCoeffs comps
                        0, // start of a_U comps
                        intvLength); // number of comps
        }
      return;
    }

  // Compute the starting component for each coefficient of a_timeInterpCoeff
  int firstComp[m_numCoeffs];
  for (int i = 0; i != m_numCoeffs; i++)
    {
      firstComp[i] = a_intvl.begin() + i * m_numStates;
    }

  // Compute the powers of the time interp coefficient
  Real interpPower[4]; // This is four because that is the power of our polynomial
  interpPower[0] = 1.; // Not really used... but kept for posterity
  interpPower[1] = a_timeInterpCoeff;
  interpPower[2] = a_timeInterpCoeff*a_timeInterpCoeff;
  interpPower[3] = a_timeInterpCoeff*a_timeInterpCoeff*a_timeInterpCoeff;

  LevelData<FArrayBox> Ucomp;
  aliasLevelData(Ucomp, &a_U, a_intvl);

  DataIterator dit = Ucomp.dataIterator();
  for(dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& UFab = Ucomp[dit];
      const FArrayBox& denseFab = m_denseCoeffs[dit];

      // Copy the old state U_n first
      UFab.copy(denseFab, firstComp[0], 0, intvLength);

      for(int ind = 1; ind != m_numCoeffs; ++ind)
        {
          // Now add in the remaining coeffs, we add in
          // U += t^ind*denseCoeff[ind]
          // where t is the a_timeInterpCoeff
          // and denseCoeff[ind] is the dense coefficient
          // corresponding to the power of the time interp coeff
          UFab.plus(denseFab,
                    interpPower[ind], // Is multiplied by denseFab
                    firstComp[ind], // Start of denseCoeffs comps
                    0, // start of a_U comps
                    intvLength); // number of comps
        }
    }
}

/*--------------------------------------------------------------------*/
///  Resets the class for a new timestep
/**
 *//*-----------------------------------------------------------------*/
void TimeInterpolatorARK4::resetData()
{
  m_gotDt = false;
  m_gotInitialSoln = false;
  m_gotFullDensePoly = false;
  m_countRHS = 0;
}

/*--------------------------------------------------------------------*/
/// Coefficients for dense output, 4th-order interpolation
/**
 *//*-----------------------------------------------------------------*/
const Real TimeInterpolatorARK4::s_bstar[][TimeInterpolatorARK4::m_numStages] = {
  {0.961753400252887, 0., 0.787405595186356, -2.74544192086633, 3.70351728061223, -1.70723435518514},
  {-1.76418754019038, 0., -0.774504669155511, 9.64023584441292, -12.544886411271, 5.44334277620397},
  {0.960350435099165, 0., 0.173858014493155, -6.21422862823726, 8.56612859966376, -3.48610842101883}
};

#include "NamespaceFooter.H"
