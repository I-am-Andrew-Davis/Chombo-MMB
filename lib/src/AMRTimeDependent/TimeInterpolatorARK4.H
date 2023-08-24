#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _TIMEINTERPOLATORARK4_H_
#define _TIMEINTERPOLATORARK4_H_

#include "FArrayBox.H"
#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "TimeInterpolatorRK4.H"
#include "NamespaceHeader.H"

/*******************************************************************************
 */
///  Time interpolator class using 6 stage 4th-order Additive Runge-Kutta
/**
 * Uses dense output to calculate a third order interpolation of the the solution
 * between two time steps. This works by saving the initial state and RHS of the
 * ODE during time integration (see LevelARK4.H). Those saved values are then used
 * with the Dense Output coefficients to calculate a third order accurate
 * intermediate time. Given a starting time of t and an ending time of t+dt,
 * the intermediate time is t+theta*dt where theta is between 0 and 1, inclusive.
 *
 * This class is derived from TimeInterpolatorRK4 and functions in the same way.
 * Since this class is derived and implements the same API, it can be used in
 * places where TimeInterpolatorRK4 are used, such as LevelGridMetrics. It is
 * meant to be a drop in replacement.
 *
 * Note that since this class is used for additive runge kutta, the saveRHS()
 * member function that accepts two arguments (RHS values for nonstiff and stiff
 * terms) must be used. Calling the single argument saveRHS() results in an error.
 *
 ******************************************************************************/

class TimeInterpolatorARK4 : public TimeInterpolatorRK4
{
public:
  /// Default constructor
  TimeInterpolatorARK4();

  /// Destructor
  ~TimeInterpolatorARK4();

  /// Full define function
  void define(const DisjointBoxLayout&  a_thisDisjointBoxLayout,
              const DisjointBoxLayout&  a_coarserDisjointBoxLayout,
              const ProblemDomain&      a_domain,
              const int&                a_refineCoarse,
              const int&                a_numStates,
              const int&                a_ghosts);

  /// Set coarse timestep.
  void setDt(const Real&  a_dt);

  /// Save the initial coarse solution
  void saveInitialSoln(const LevelData<FArrayBox>&   a_soln);

  /// Save coarse RHS stage value
  void saveRHS(const LevelData<FArrayBox>&   a_rhs);

  /// Save coarse implicit and explicit RHS stage value
  void saveRHS(const LevelData<FArrayBox>& a_ExpRHS,
               const LevelData<FArrayBox>& a_ImpRHS);

  /// Interpolate in time using the Dense Output polynomial.
  void interpolate(LevelData<FArrayBox>&   a_U,
                   const Real&             a_timeInterpCoeff,
                   const Interval&         a_intvl);

  /// Set RK4 intermediate values in time using the Dense Output polynomial.
  void intermediate(LevelData<FArrayBox>&   a_U,
                    const Real&             a_timeInterpCoeff,
                    const Real              a_dtRatio,
                    const int&              a_stage,
                    const Interval&         a_intvl) const;

protected:
  /// Reset this object to use with new data
  void resetData();

  static constexpr int m_numCoeffs = 4;
                                    ///< number of coefficients of Dense Output:  this is 4
                                    ///< 1st corresponds to U_n (previous state)
                                    ///< The remaining three correspond to the three powers of theta
  static constexpr int m_numStages = 6;
                                    ///< Number of stages for this Dense Output: this is always 6

  static const Real s_bstar[m_numCoeffs-1][m_numStages];
                                    ///< These are the coefficients for the dense output
  DisjointBoxLayout m_thisCoarsenedLayout;
                                    ///< Layout for this level, coarsened by m_refToCoarse
  DisjointBoxLayout m_coarseLayout; ///< Layout for the coarse level
  Copier m_copier;                  ///< For copying from rhs on m_coarseLayout to m_rhsCopy on m_thisCoarsenedLayout
  int m_ghosts;                     ///< Layers of ghost cells around m_thisCoarsenedLayout for m_rhsCopy and m_DenseCoeffs
  IntVect m_ghostVect;              ///< Ghost vector around m_thisCoarsenedLayout for m_rhsCopy and m_DenseCoeffs
  Real m_dt;                        ///< Coarse timestep
  LevelData<FArrayBox> m_rhsCopy;   ///< Copy of rhs on m_thisCoarsenedLayout, to be used within saveRHS
  LevelData<FArrayBox> m_denseCoeffs;
                                    ///< coefficients of the third-order Dense Output on m_thisCoarsenedLayout with ghost vector m_ghostVect; m_numCoeffs*m_numStates components
  bool m_gotDt;                     ///< Whether m_dt has been set
  bool m_gotInitialSoln;            ///< Whether initial solution has been saved
  bool m_gotFullDensePoly;          ///< Whether we have the full Dense Output polynomial
  int m_countRHS;                   ///< Number of times saveRHS function has been called
  ProblemDomain m_coarseDomain;     ///< Problem domain - index space for next coarser level
  int m_refineCoarse;               ///< Refinement ratio between this level and the next coarser
  int m_numStates;                  ///< Number of variables
  bool m_defined;                   ///< Whether define() has been called

private:

  // Disallowed for all the usual reasons
  void operator=(const TimeInterpolatorARK4&);
  TimeInterpolatorARK4(const TimeInterpolatorARK4&);
};

#include "NamespaceFooter.H"
#endif
