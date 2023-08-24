#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// #include <cstdio>
#include <array>

#include "DebugOut.H"
#include "TimeInterpolatorRK2.H"
#include "NamespaceHeader.H"

//////////////////////////////////////////////////////////////////////////////
// Constructor - set up some defaults
TimeInterpolatorRK2::TimeInterpolatorRK2()
{
  m_defined = false;
  m_coarseDefined = false;
  m_fineDefined = false;
  resetData();
}


//////////////////////////////////////////////////////////////////////////////
// Destructor - free up storage
TimeInterpolatorRK2::~TimeInterpolatorRK2()
{
}


//////////////////////////////////////////////////////////////////////////////
// Define the object so that time stepping can begin
void TimeInterpolatorRK2::define(/// layout at this level
                                 const DisjointBoxLayout&  a_thisDisjointBoxLayout,
                                 /// layout at next coarser level
                                 const DisjointBoxLayout&  a_coarserDisjointBoxLayout,
                                 /// problem domain on this level
                                 const ProblemDomain&      a_domain,
                                 /// refinement ratio between this level and next coarser level
                                 const int&                a_refineCoarse,
                                 /// number of variables
                                 const int&                a_numStates,
                                 /// layers of ghost cells to be filled in on the coarsened layout at this level
                                 const int&                a_ghosts)
{
  // Cache data
  m_refineCoarse = a_refineCoarse;
  m_coarseDomain = coarsen(a_domain, m_refineCoarse);
  m_numStates = a_numStates;
  m_ghosts = a_ghosts;
  m_ghostVect = m_ghosts * IntVect::Unit;

  m_coarseLayout = a_coarserDisjointBoxLayout;
  // petermc, 19 Dec 2008:  prevents crash in case of calling define again
  m_thisCoarsenedLayout = DisjointBoxLayout();
  coarsen(m_thisCoarsenedLayout, a_thisDisjointBoxLayout, m_refineCoarse);

  m_initialSolution.define(m_thisCoarsenedLayout, m_numStates, m_ghostVect);
  m_finalSolution.define(m_thisCoarsenedLayout, m_numStates, m_ghostVect);

  // Even though all of the member LevelData are on the coarse layout, we still
  // create the copier because most of the calls to interpolate will be from
  // coarse layout to coarsened-fine layout
  m_copier.define(m_coarseLayout, m_thisCoarsenedLayout,
                  m_coarseDomain, m_ghostVect);
  // Everything is defined now.
  m_coarseDefined = true;
  m_fineDefined = true;
  m_defined = true;
}

bool TimeInterpolatorRK2::isDefined() const
{
  return m_defined;
}


//////////////////////////////////////////////////////////////////////////////
void TimeInterpolatorRK2::setDt(const Real&  a_dt)
{
  resetData();
  m_dt = a_dt;
  m_gotDt = true;

  // This is really a no-op for RK2
}

//////////////////////////////////////////////////////////////////////////////
void TimeInterpolatorRK2::saveInitialSoln(const LevelData<FArrayBox>&   a_soln)
{
  CH_assert(m_defined);
  CH_assert(m_coarseDefined);
  CH_assert(!m_gotInitialSoln);
  CH_assert(a_soln.nComp() == m_numStates);

  // First zero out member initial solution
  DataIterator dit = m_initialSolution.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& initialFab = m_initialSolution[dit];
      initialFab.setVal(0.);
    }

  // a_soln is on the coarse layout,
  // m_initialSolution is on the coarsened-fine layout
  const Interval& srcInt = a_soln.interval();
  a_soln.copyTo(srcInt, m_initialSolution, srcInt, m_copier);

  m_gotInitialSoln = true;
}

void TimeInterpolatorRK2::saveFinalSoln(const LevelData<FArrayBox>&   a_soln)
{
  CH_assert(m_defined);
  CH_assert(m_coarseDefined);
  CH_assert(m_gotInitialSoln);
  CH_assert(!m_gotFinalSoln);
  CH_assert(a_soln.nComp() == m_numStates);

  // First zero out member final solution
  DataIterator dit = m_finalSolution.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& finalFab = m_finalSolution[dit];
      finalFab.setVal(0.);
    }

  // a_soln is on the coarse layout,
  // m_finalSolution is on the coarsened-fine layout
  const Interval& srcInt = a_soln.interval();
  a_soln.copyTo(srcInt, m_finalSolution, srcInt, m_copier);

  m_gotFinalSoln = true;
}


//////////////////////////////////////////////////////////////////////////////
void TimeInterpolatorRK2::saveRHS(const LevelData<FArrayBox>&   a_rhs)
{
  // This is a null OP for RK2
  return;
}


//////////////////////////////////////////////////////////////////////////////
void TimeInterpolatorRK2::intermediate(/// intermediate RK2 solution on this level coarsened
                                       LevelData<FArrayBox>&   a_U,
                                       /// time interpolation coefficient in range [0:1]
                                       const Real&             a_timeInterpCoeff,
                                       /// Ratio of fine to coarse time steps.  This is usually 1/nref unless extra subcycling occurs.
                                       const Real              a_dtRatio,
                                       /// which RK2 stage:  0, 1, 2, 3
                                       const int&              a_stage,
                                       /// interval of a_U to fill in
                                       const Interval&         a_intvl) const
{
  CH_assert(m_defined);
  CH_assert(m_coarseDefined);
  CH_assert(a_U.nComp() == m_numStates);
  CH_assert(m_gotInitialSoln);
  CH_assert(m_gotFinalSoln);

  const int intvLength = a_intvl.size();

  // We are at the initial solution time
  if(a_timeInterpCoeff == 0.)
    {
      for (DataIterator dit = a_U.dataIterator(); dit.ok(); ++dit)
        {
          a_U[dit].copy(m_initialSolution[dit], 0, 0, intvLength);
        }
      return;
    }
  else if (a_timeInterpCoeff == 1.) // We are at the final solution time
    {
      for (DataIterator dit = a_U.dataIterator(); dit.ok(); ++dit)
        {
          a_U[dit].copy(m_finalSolution[dit], 0, 0, intvLength);
        }
      return;
    }
  else // We are somewhere in between, do interpolation
    {
      for(DataIterator dit = a_U.dataIterator(); dit.ok(); ++dit)
        {
          const FArrayBox& initialFab = m_initialSolution[dit];
          const FArrayBox& finalFab = m_finalSolution[dit];
          FArrayBox& interpToFab = a_U[dit];
          for(int comp = 0; comp != intvLength; ++comp)
          MD_BOXLOOP(initialFab.box(),cell)
            {
              // Compute U_0(1-coeff)+U_1(coeff)
              interpToFab[MD_IX(cell,comp)] =
                initialFab[MD_IX(cell,comp)]*(1-a_timeInterpCoeff)
                + finalFab[MD_IX(cell,comp)]*a_timeInterpCoeff;
            }
        }
    }


  // dummy statement in order to get around gdb bug
  int dummy_unused = 0; dummy_unused = 0;
}

//////////////////////////////////////////////////////////////////////////////
void TimeInterpolatorRK2::resetData()
{
  m_gotDt = false;
  m_gotInitialSoln = false;
  m_gotFullTaylorPoly = false;
  m_gotFinalSoln = false;
  m_countRHS = 0;
}

#include "NamespaceFooter.H"
