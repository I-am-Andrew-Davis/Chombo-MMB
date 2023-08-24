#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LEVELARK4_H_
#define _LEVELARK4_H_

// An enum to let the user select the method of stage value
// predictor/initial guess they want
namespace ARK4
{
  enum StageValuePredictorMethod
    {
     PreviousStageValues, ///< Uses the previous stage value as
                          ///  initial guess to next stage
     Extrapolate2ndOrder  ///< 2nd Order extrapolation of previous
                          ///  solution and stage function values.
                          ///  Will not work on first time step
                          ///  (and will revert to PreviousStageValues)
    };
}

/*******************************************************************************
 */
///  Templated implementation of ARK4 advance for a single AMRLevel, allowing interpolation in time
/**
 *   This implements the fourth-order Additive Runge-Kutta time advance of an ODE,
 *   and saves intermediate results for use in interpolation in time.
 *   ODE is d(Soln)/dt = RHS
 *
 *   Template types:
 *   \li TSoln is the datatype for the solution
 *   \li TRhs is the datatype for the RHS (could be same as solution...)
 *   \li TInterp is a storage class for intermediate values
 *   \li TOp is an object encapsulating the actual ODE operations
 *
 *   TInterp requirements: TInterp must have the following functions:
 *   \code
 *   // Set time step
 *   void setDt(Real dt)
 *
 *   // Save initial solution
 *   void saveInitialSoln(TSoln& soln)
 *
 *   // save the non-stiff and stiff terms of the RHS
 *   void saveRHS(TSoln& nonstiffRHS, TSoln& stiffRHS)
 *   \endcode
 *
 *   TOp requirements: TOp must have the following functions:
 *   \code
 *   // Evaluates non-stiff term
 *   void evalNonStiff(RHS&         FnonStiff,
 *                     TSoln&       stageSoln,
 *                     int          stage,
 *                     Real         stageTime,
 *                     Real         stageWeight,
 *                     TOpArgs&&... forwardedOpArgs)
 *
 *   // Evaluates stiff term
 *   void evalStiff(RHS&         Fstiff,
 *                  TSoln&       stageSoln,
 *                  int          stage,
 *                  Real         stageTime,
 *                  Real         stageWeight,
 *                  TOpArgs&&... forwardedOpArgs)
 *
 *   // Increments destRHS += val*evalRHS
 *   void increment(RHS&  destRHS,
 *                  RHS&  evalRHS,
 *                  Real& val)
 *
 *  // Increments soln += val*RHS
 *  void increment(TSoln& soln,
 *                 RHS&   rhs,
 *                 Real   val)
 *
 *  // Solves a nonlinear problem
 *  void solve(RHS&         stageSoln,
 *             TSoln&       initialGuess,
 *             TSoln&       previousTimeStep,
 *             RHS&         previousStagesRHS,
 *             int          currentStage,
 *             Real         stageTime,
 *             Real         stageWeight,
 *             TOpArgs&&... forwardedOpArgs)
 *
 *   // implements soln += dt*rhs
 *   void updateODE(TSoln& soln,
 *                  const TSoln& rhs,
 *                  Real dt)
 *
 *   // define data holder newSoln based on existingSoln,
 *   // including ghost cell specification
 *   void defineSolnData(TSoln& newSoln,
 *                       const TSoln& existingSoln)
 *
 *   // define data holder for RHS based on existingSoln
 *   // including ghost cell specification
 *   // (which in most cases is no ghost cells)
 *   void defineRHSData(TSoln& newRHS, const TSoln& existingSoln)
 *
 *   /// copy data in src into dest
 *   void copySolnData(TSoln& dest, const TSoln& src)
 *
 *   /// Allows the user to compute the error estimate.
 *   // Will only be called if LevelARK4::advance has the
 *   // parameter a_computeErrorEstimate set to true.
 *   // The user must return the global max norm between
 *   // the solution U and the 3rd order solution estimate
 *   // Uhat.
 *   Real computeErrorEstimate(TSoln& U, TSoln& Uhat)
 *   \endcode
 *
 ******************************************************************************/

template <typename TSoln, typename TRhs, typename TInterp, typename TOp>
class LevelARK4
{
public:
  /// Default constructor
  LevelARK4();

  /// Destructor
  ~LevelARK4();

  /// Constructs internal data structures, prepares for time stepping
  void define(const DisjointBoxLayout& a_dbl,
              int                      a_numComp,
              const IntVect&           a_ghostCells);


  /// Integrates in time using the 6 stage ARK4 method
  template <typename... TOpArgs>
  Real advance(TSoln&                                a_newSoln,
               const TSoln&                          a_oldSoln,
               Real                                  a_time,
               Real                                  a_dt,
               Real                                  a_prevStepDt,
               TInterp*                              a_timeInterp,
               TOp&                                  a_levelOp,
               const ARK4::StageValuePredictorMethod a_predictorMethod,
               bool                                  a_computeErrorEstimate,
               TOpArgs&&...                          a_opArgs);

protected:
  static const int  s_nStages = 6;  //< Number of stages
  static const Real s_gamma;        //< Butcher Tableau constant
  static const Real s_aE[s_nStages][s_nStages];
                                    //< Nonstiff method coefficients
  static const Real s_aI[s_nStages][s_nStages];
                                    //< Stiff method coefficients
  static const Real s_b[s_nStages]; //< Stage weight coefficients
  static const Real s_bhat[s_nStages];
                                    //< 3rd order embedded method coefficients
  static const Real s_bstar2ndOrder[2][s_nStages];
                                    //< 2nd order dense output coefficients
  static const Real s_c[s_nStages]; //< Stage time coefficients

  TSoln m_state[s_nStages];         //< Stage values
  TRhs m_FnonStiff[s_nStages];      //< Evaluations of the nonstiff term
  TRhs m_Fstiff[s_nStages];         //< Evaluations of the stiff term
};

/*==============================================================================
 * Constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
///  Default constructor
/**  Does nothing at this time
 *//*-----------------------------------------------------------------*/
template <typename TSoln, typename TRhs, typename TInterp, typename TOp>
LevelARK4<TSoln, TRhs, TInterp, TOp>::LevelARK4()
{
}

/*--------------------------------------------------------------------*/
///  Default destructor
/**  Does nothing at this time
 *//*-----------------------------------------------------------------*/
template <typename TSoln, typename TRhs, typename TInterp, typename TOp>
LevelARK4<TSoln, TRhs, TInterp, TOp>::~LevelARK4()
{
}

/*==============================================================================
 * Member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
///  Prepares this class for time stepping
/**  Constructs internal data structures needed for ARK
 *   \param[in]  a_dbl     Box layout for solution data that will be stepped
 *   \param[in]  a_numComp Number of components in solution data that will be stepped
 *//*-----------------------------------------------------------------*/
template <typename TSoln, typename TRhs, typename TInterp, typename TOp>
void LevelARK4<TSoln, TRhs, TInterp, TOp>::define(const DisjointBoxLayout& a_dbl,
                                                  int                      a_numComp,
                                                  const IntVect&           a_ghostCells)
{
  for (int stage = 0; stage < s_nStages; ++stage)
    {
      m_state[stage].define(a_dbl, a_numComp, a_ghostCells);
      m_FnonStiff[stage].define(a_dbl, a_numComp);
      m_Fstiff[stage].define(a_dbl, a_numComp);
    }
}

/*--------------------------------------------------------------------*/
///  Advances the solution forward in time
/**  Integrates the solution in time using the 6-stage ARK4 method,
 *   calling the evalNonStiff() and evalStiff() method from TOp.
 *
 *   \param[out] a_newSoln    Solution after time integration
 *   \param[in]  a_oldSoln    Solution before time integration
 *   \param[in]  a_time       Time before doing time integration
 *   \param[in]  a_dt         Timestep to be taken by this time integration
 *   \param[in]  a_prevStepDt The time step size from the previous time integration
 *   \param[out] a_timeInterp Time interpolator filled with data from this level
 *   \param[in]  a_levelOp    Object that encapsulates the ODE functionality
 *                            (must meet the requirements detailed above)
 *   \param[in]  a_predictorMethod
 *                            Which method should be used to provide the initial guess
 *                            for the nonlinear solver
 *   \param[in]  a_computeErrorEstimate
 *                            Compute the error estimate using the 3rd order embedded
 *                            method. The error estimate is the max norm of U-\hat{U}
 *                            where U is the solution associated with the 4th order
 *                            ARK and \hat{U} is the solution associated with the 3rd
 *                            order embedded method.
 *   \param[in]  a_opArgs     Additional arguments to forward to evalNonStiff() and evalStiff()
 *   \return     Returns the error estimate if a_computeErrorEstimate is true. If
 *               a_computeErrorEstimate is false, returns -1. If the error estimate is
 *               computed, the returned value will be the max norm from fall processors.
 *//*-----------------------------------------------------------------*/

template <typename TSoln, typename TRhs, typename TInterp, typename TOp>
template <typename... TOpArgs>
Real LevelARK4<TSoln, TRhs, TInterp, TOp>::advance(
     TSoln&                                a_newSoln,
     const TSoln&                          a_oldSoln,
     Real                                  a_time,
     Real                                  a_dt,
     Real                                  a_prevStepDt,
     TInterp*                              a_timeInterp,
     TOp&                                  a_levelOp,
     const ARK4::StageValuePredictorMethod a_predictorMethod,
     bool                                  a_computeErrorEstimate,
     TOpArgs&&...                          a_opArgs)
{
  CH_assert(m_state[0].isDefined());
  // Stage 0 always uses the previous solution
  a_levelOp.copySolnData(m_state[0], a_oldSoln);
  bool extrapStageValues = false;


  // Calculate stage predictor values if requested
  if(a_prevStepDt > 0. && a_predictorMethod == ARK4::Extrapolate2ndOrder)
    {
      extrapStageValues = true;
      Real ratio = a_dt/a_prevStepDt;
      // Computes U^i
      for (int stage = 1; stage < s_nStages; ++stage)
        {
          CH_assert(m_state[stage].isDefined());
          Real theta = 1 + ratio*s_c[stage];
          a_levelOp.copySolnData(m_state[stage], a_oldSoln);
          // Computes the dense output dt*sum b*_j(theta_i)(F^j_ns + F^j_s)
          for (int Sdense = 0; Sdense < s_nStages; ++Sdense)
            {
              CH_assert(m_FnonStiff[Sdense].isDefined());
              CH_assert(m_Fstiff[Sdense].isDefined());
              Real bstar = 0.;
              // Computes bstar
              for (int ind = 0; ind < 2; ++ind)
                {
                  bstar += s_bstar2ndOrder[ind][Sdense]*std::pow(theta,1+ind);
                }
              bstar *= a_dt;
              a_levelOp.increment(m_state[stage], m_FnonStiff[Sdense], bstar);
              a_levelOp.increment(m_state[stage], m_Fstiff[Sdense], bstar);
            }
        }
    }

  // Calculate stage times
  Real t[s_nStages];
  Real stageweight[s_nStages];
  for (int stage = 0; stage < s_nStages; ++stage)
    {
      stageweight[stage] = a_dt * s_b[stage];
      t[stage] = a_time + s_c[stage] * a_dt;
    }

  // Setup the time interpolator
  if(a_timeInterp != nullptr)
    {
      a_timeInterp->setDt(a_dt);
      a_timeInterp->saveInitialSoln(a_oldSoln);
    }

  // Compute non-stiff and stiff terms at first stage
  a_levelOp.evalNonStiff(m_FnonStiff[0], // All the zeros here are the stage number
                         m_state[0],
                         0,
                         t[0],
                         stageweight[0],
                         std::forward<TOpArgs>(a_opArgs)...);

  a_levelOp.evalStiff(m_Fstiff[0], // Same with the zeros here
                      m_state[0],
                      0,
                      t[0],
                      stageweight[0],
                      std::forward<TOpArgs>(a_opArgs)...);
  // Save the RHS for the first stage
  if(a_timeInterp != nullptr)
    {
      a_timeInterp->saveRHS(m_FnonStiff[0],m_Fstiff[0]);
    }

  for(int stage = 1; stage != s_nStages; ++stage)
    {
      TRhs X_vec;
      a_levelOp.defineRHSData(X_vec, a_oldSoln);

      // Compute X
      for(int j = 0; j != stage; ++j)
        {
          a_levelOp.increment(X_vec, m_FnonStiff[j], a_dt*s_aE[stage][j]);
          a_levelOp.increment(X_vec, m_Fstiff[j], a_dt*s_aI[stage][j]);
        }

      // Solve the Ax=b problem (A is computed numerically, x is this stage
      // stateb and b is X_vec plus a few other things added in by the solve
      // function)
      int predictorIdx = stage-1;
      if (extrapStageValues)
        {
          predictorIdx = stage;
        }
      a_levelOp.solve(m_state[stage],
                      m_state[predictorIdx], // Initial guess
                      m_state[0], // Solution from last time step
                      X_vec, // The RHS
                      stage, // The current stage
                      t[stage], // The time of this stage
                      stageweight[stage], // Stage weight (for flux corrections)
                      std::forward<TOpArgs>(a_opArgs)...);

      // X of this stage is based off Fnonstiff and Fstiff of the previous stage
      a_levelOp.evalNonStiff(m_FnonStiff[stage],
                             m_state[stage],
                             stage,
                             t[stage],
                             stageweight[stage],
                             std::forward<TOpArgs>(a_opArgs)...);

      a_levelOp.evalStiff(m_Fstiff[stage],
                          m_state[stage],
                          stage,
                          t[stage],
                          stageweight[stage],
                          std::forward<TOpArgs>(a_opArgs)...);

      // Save the RHS for the next stages
      if(a_timeInterp != nullptr)
        {
          a_timeInterp->saveRHS(m_FnonStiff[stage], m_Fstiff[stage]);
        }
    }

  // Computes U^{n+1}
  a_levelOp.copySolnData(a_newSoln, a_oldSoln);
  for(int j = 0; j != s_nStages; ++j)
    {
      // Standard U^{n+1}
      a_levelOp.increment(a_newSoln, m_FnonStiff[j], a_dt*s_b[j]);
      a_levelOp.increment(a_newSoln, m_Fstiff[j], a_dt*s_b[j]);
    }

  // Compute the error estimate if requested
  Real errorEstimate = -1.;
  if (a_computeErrorEstimate)
    {
      TSoln UhatTmp;
      a_levelOp.defineSolnData(UhatTmp, a_oldSoln);
      a_levelOp.copySolnData(UhatTmp, a_oldSoln);
      for(int j = 0; j != s_nStages; ++j)
        {
          // Embedded method \hat{U}^{n+1}
          a_levelOp.increment(UhatTmp, m_FnonStiff[j], a_dt*s_bhat[j]);
          a_levelOp.increment(UhatTmp, m_Fstiff[j], a_dt*s_bhat[j]);
        }
      errorEstimate = a_levelOp.computeErrorEstimate(a_newSoln, UhatTmp);
    }

  return errorEstimate;
}

/* Define the butcher tableaux values and other constants used by the ARK4 method */
template <typename TSoln, typename TRhs, typename TInterp, typename TOp>
const Real LevelARK4<TSoln, TRhs, TInterp, TOp>::s_gamma = 0.25;

/*--------------------------------------------------------------------*/
///  Stage time values
/**  A static variable that holds the stage time values for this ARK4
 *   method. Values are from the Butcher Tableau.
 *//*-----------------------------------------------------------------*/
template <typename TSoln, typename TRhs, typename TInterp, typename TOp>
const Real LevelARK4<TSoln, TRhs, TInterp, TOp>::s_c[] = { 0.0, 0.5, 0.332, 0.62, 0.85, 1.0 };

/*--------------------------------------------------------------------*/
///  Explicit Butcher Tableau coefficients
/**  A static variable that holds the explicit coefficients for this ARK4
 *   method. Values are from the Butcher Tableau.
 *//*-----------------------------------------------------------------*/
template <typename TSoln, typename TRhs, typename TInterp, typename TOp>
const Real LevelARK4<TSoln, TRhs, TInterp, TOp>::s_aE[][LevelARK4<TSoln, TRhs, TInterp, TOp>::s_nStages] = {
  {0., 0., 0., 0., 0., 0.},
  {0.5, 0., 0., 0., 0., 0.},
  {0.221776, 0.110224, 0., 0., 0., 0.},
  {-0.04884659515311857, -0.17772065232640102, 0.8465672474795197, 0., 0., 0.},
  {-0.15541685842491548, -0.3567050098221991, 1.0587258798684427, 0.30339598837867193, 0., 0.},
  { 0.2014243506726763, 0.008742057842904185, 0.15993995707168115, 0.4038290605220775, 0.22606457389066084, 0.}
};

/*--------------------------------------------------------------------*/
///  Implicit Butcher Tableau coefficients
/**  A static variable that holds the implicit coefficients for this ARK4
 *   method. Values are from the Butcher Tableau.
 *//*-----------------------------------------------------------------*/
template <typename TSoln, typename TRhs, typename TInterp, typename TOp>
const Real LevelARK4<TSoln, TRhs, TInterp, TOp>::s_aI[][LevelARK4<TSoln, TRhs, TInterp, TOp>::s_nStages] = {
  {0., 0., 0., 0., 0., 0.},
  {0.25, 0.25, 0., 0., 0., 0.},
  {0.137776, -0.055776, 0.25, 0., 0., 0.},
  {0.14463686602698217, -0.22393190761334475, 0.4492950415863626, 0.25, 0., 0.},
  {0.09825878328356477, -0.5915442428196704, 0.8101210538282996, 0.283164405707806, 0.25, 0.},
  {0.15791629516167136, 0., 0.18675894052400077, 0.6805652953093346, -0.27524053099500667, 0.25}
};

/*--------------------------------------------------------------------*/
///  Sets the stage weights for this ARK4 method
/**  A static variable that holds the stage weight values for this ARK4
 *   method. Values are from the Butcher Tableau. Each column is a stage.
 *//*-----------------------------------------------------------------*/
template <typename TSoln, typename TRhs, typename TInterp, typename TOp>
const Real LevelARK4<TSoln, TRhs, TInterp, TOp>::s_b[] =
  {0.15791629516167136, 0., 0.18675894052400077, 0.6805652953093346, -0.27524053099500667, 0.25};

/*--------------------------------------------------------------------*/
///  Embedded method coefficients
/**  A static variable that holds the stage weight values for this ARK4
 *   method. Values are from the Butcher Tableau. Each column is a stage.
 *//*-----------------------------------------------------------------*/
template <typename TSoln, typename TRhs, typename TInterp, typename TOp>
const Real LevelARK4<TSoln, TRhs, TInterp, TOp>::s_bhat[] =
  {0.154711800763212, 0., 0.189205191660680, 0.702045371228922, -0.319187399063579, 0.273225035410765};

/*--------------------------------------------------------------------*/
///  Second order dense output coefficients used for extrapolation of
///  next time step stage values
/**  A static variable that holds the 2nd order dense output coefficients
 *   for ARK4. Values are from the Butcher Tableau. Each column is a stage.
 *//*-----------------------------------------------------------------*/
template <typename TSoln, typename TRhs, typename TInterp, typename TOp>
constexpr Real LevelARK4<TSoln, TRhs, TInterp, TOp>::s_bstar2ndOrder[][LevelARK4<TSoln, TRhs, TInterp, TOp>::s_nStages] = {
  {0.924880883935497, 0., 0.738545394730690, -2.517224551339271, 3.588774145955500, -1.734975873282416},
  {-0.766964588773826, 0., -0.551786454206689, 3.197789846648605, -3.864014676950506, 1.984975873282416}
};

#endif // _LEVELARK4_H_
