#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _TIMEINTERPOLATORRK2_H_
#define _TIMEINTERPOLATORRK2_H_

#include "FArrayBox.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "TimeInterpolatorRK4.H"
#include "NamespaceHeader.H"

/// Time interpolator class using 4th-order Runge-Kutta

/**
 */
class TimeInterpolatorRK2 : public TimeInterpolatorRK4
{
public:
  /// Default constructor
  /**
     Object requires define() to be called before all other functions.
   */
  TimeInterpolatorRK2();

  /// Destructor
  /**
     Destroys all objects created by define(). Passed in data references
     of define() are left alone.
   */
  virtual ~TimeInterpolatorRK2();

  /// Actual constructor.
  /**
     Set up object.
   */
  virtual void define(/// layout at this level
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
                      const int&                a_ghosts) override;

  /// Whether or not this time interpolator is defined
  bool isDefined() const;

  /// Set coarse timestep.
  virtual void setDt(const Real&  a_dt) override;

  /// Save the initial coarse solution time
  /**
     Save the initial coarse solution time for use in interpolation.
     Saves ghost cells, too.
   */
  virtual void saveInitialSoln(const LevelData<FArrayBox>&   a_soln) override;

  /// Save the final coarse solution time
  /**
     Save the final coarse solution time for use in interpolation.
     Saves ghost cells, too.
   */
  virtual void saveFinalSoln(const LevelData<FArrayBox>&   a_soln);

  /// Update Taylor polynomial coefficients with coarse right-hand side.
  /**
  * This is a null op for RK2
  */
  virtual void saveRHS(const LevelData<FArrayBox>&   a_rhs) override;

  /// Interpolate in time using the Taylor polynomial.
  /**
     Interpolate in time to a_U on interval a_intvl using the Taylor polynomial.
     Ghost cells, too.
  */

  /// Set RK2 intermediate values in time using the Taylor polynomial.
  /**
     Set RK4 intermediate values to a_U on interval a_intvl
     using the Taylor polynomial.
     Ghost cells, too.
  */
  virtual void intermediate(/// interpolated solution on this level coarsened
                            LevelData<FArrayBox>&   a_U,
                            /// time interpolation coefficient in range [0:1]
                            const Real&             a_timeInterpCoeff,
                            /// Ratio of fine to coarse time steps.  This is usually 1/nref unless extra subcycling occurs.
                            const Real              a_dtRatio,
                            /// which RK4 stage:  0, 1, 2, 3
                            const int&              a_stage,
                            /// interval of a_U to fill in
                            const Interval&         a_intvl) const override;

  /// Reset this object to use with new data
  void resetData();

protected:

  /// layout for this level, coarsened by m_refToCoarse
  DisjointBoxLayout m_thisCoarsenedLayout;

  /// layout for the coarse level
  DisjointBoxLayout m_coarseLayout;

  /// For copying from rhs on m_coarseLayout to m_rhsCopy on m_thisCoarsenedLayout
  Copier m_copier;

  /// layers of ghost cells around m_thisCoarsenedLayout for m_rhsCopy and m_taylorCoeffs
  int m_ghosts;

  /// ghost vector around m_thisCoarsenedLayout for m_rhsCopy and m_taylorCoeffs
  IntVect m_ghostVect;

  /// coarse timestep
  Real m_dt;

  /// Copy of rhs on m_thisCoarsenedLayout, to be used within saveRHS
  LevelData<FArrayBox> m_initialSolution;

  /// coefficients of the third-degree Taylor polynomial on m_thisCoarsenedLayout with ghost vector m_ghostVect; m_numCoeffs*m_numStates components
  LevelData<FArrayBox> m_finalSolution;

  /// whether m_dt has been set
  bool m_gotDt;

  /// whether initial solution has been saved
  bool m_gotInitialSoln;

  /// saveFinalSolution() has been called
  bool m_gotFinalSoln;

  /// Problem domain - index space for next coarser level
  ProblemDomain m_coarseDomain;

  /// Refinement ratio between this level and the next coarser
  int m_refineCoarse;

  /// Number of variables
  int m_numStates;

  /// define() has been called
  bool m_defined;

  /// defineCoarseParts() has been called
  bool m_coarseDefined;

  /// defineFineParts() has been called
  bool m_fineDefined;

private:

  // Disallowed for all the usual reasons
  void operator=(const TimeInterpolatorRK2&);
  TimeInterpolatorRK2(const TimeInterpolatorRK2&);
};

#include "NamespaceFooter.H"
#endif
