#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRLEVELCONSFACTORY_H_
#define _AMRLEVELCONSFACTORY_H_

#include "AMRLevelFactory.H"
#include "AMRLevelCons.H"

#include "NamespaceHeader.H"

/// AMR Cons Equation factory
/**
 */
class AMRLevelConsFactory : public AMRLevelFactory
{
public:
  /// Null constructor
  /**
   */
  AMRLevelConsFactory();

  /// Virtual constructor
  /**
   */
  virtual AMRLevel* new_amrlevel() const;

  /// Destructor
  /**
   */
  virtual ~AMRLevelConsFactory();

  /// CFL number
  /**
   */
  virtual void CFL(Real a_cfl);

  /// spatial order of accuracy (default is 4)
  virtual void spaceOrder(int a_spaceOrder);

  /// sets whether to limit face values in advection (default is false)
  virtual void limitFaceValues(bool a_limitFaceValues);

  /// sets whether to use high-order limiter (default is true)
  virtual void highOrderLimiter(bool a_highOrderLimiter);

  /// sets whether initial data is average; default false
  virtual void initialAverage(bool a_initialAverage);

  /// sets whether to flatten extrapolations to faces
  virtual void useFlattening(bool a_useFlattening);

  /// sets whether to turn PPM off
  virtual void noPPM(bool a_noPPM);

  /// sets m_doDeconvolution
  virtual void doDeconvolution(bool a_doDeconvolution);

  /// sets m_doFaceDeconvolution
  virtual void doFaceDeconvolution(bool a_doFaceDeconvolution);

  /// sets whether to use artificial viscosity (from divergence)
  virtual void useArtificialViscosity(bool a_useArtificialViscosity);

  /// sets coefficient of artificial viscosity (from divergence)
  virtual void artificialViscosity(Real a_artificialViscosity);

  /// sets whether to use artificial viscosity
  virtual void useArtVisc(bool a_useArtVisc);

  /// sets coefficient of artificial viscosity
  virtual void ratioArtVisc(Real a_ratioArtVisc);

  /// sets whether forward Euler is used (instead of RK4)
  virtual void forwardEuler(bool a_forwardEuler);

  /// sets whether to enforce a min value in advection, along with valeu
  virtual void enforceMinVal(bool a_enforceMinVal, Real a_minVal);

  /// Physical dimension of the longest side of the domain
  /**
   */
  virtual void domainLength(Real a_domainLength);

  void x0(const Real& a_x0);

  ///
  /**
   */
  virtual void verbosity(const int& verbosity);

  /// Refinement threshold
  /**
   */
  virtual void refinementThreshold(Real a_refineThresh);

  /// Set whether refinement threshold is scaled with dx
  /**
   */
  virtual void refinementIsScaled(bool a_refinementIsScaled);

  /// Set whether to tag on pressure instead of on density
  /**
   */
  virtual void tagPressure(bool a_tagPressure);

  /// Set whether to tag on vorticity instead of on density
  /**
   */
  virtual void tagVorticity(bool a_tagVorticity);

  /// Tag buffer size
  /**
   */
  virtual void tagBufferSize(int a_tagBufferSize);

  /// MOLPhysics class
  /**
   */
  virtual void molPhysics(const MOLPhysics* const a_molPhysics);

  /// Initial dt multiplier
  /**
   */
  void initialDtMultiplier(Real a_initialDtMultiplier);

  /// Check that everything is defined
  /**
   */
  virtual bool isDefined() const;

protected:

  // Some default values
  virtual void setDefaultValues();

  virtual void transferSettings(AMRLevelCons* a_newPtr) const;

  int m_verbosity;

  // CFL number
  Real m_cfl;
  bool m_cflSet;
  //  PhysIBC* m_cons_ibc;

  /// spatial order of accuracy (default is 4)
  int m_spaceOrder;

  /// if true, limit face values in advection
  bool m_limitFaceValues;

  /// if true, use high-order limiter
  bool m_highOrderLimiter;

  // if true, initial values are average
  bool m_initialAverage;

  // if true, use flattening
  bool m_useFlattening;

  // Avoid PPM?  Default false.
  bool m_noPPM;

  // Do deconvolution?  Default yes.
  bool m_doDeconvolution;

  // Do face deconvolution?  Default yes.
  bool m_doFaceDeconvolution;

  // use artificial viscosity (from divergence)?
  bool m_useArtificialViscosity;

  // coefficient of artificial viscosity (from divergence)
  Real m_artificialViscosity;

  // if true, use artificial viscosity
  bool m_useArtVisc;

  // if m_useArtVisc, artificial viscosity coefficient
  Real m_ratioArtVisc;

  // use forward Euler if true, or RK4 if false
  bool m_forwardEuler;

  /// if true, enforce min value
  bool m_enforceMinVal;

  /// min value to enforce
  Real m_minVal;

  // Physical dimension of the longest side of the domain
  Real m_domainLength;
  bool m_domainLengthSet;

  // Refinement threshold for gradient
  Real m_refineThresh;
  bool m_refineThreshSet;

  // Whether refinement threshold is scaled with dx
  bool m_refinementIsScaled;

  // Tag on pressure instead of on density?
  bool m_tagPressure;

  // Tag on vorticity instead of on density?
  bool m_tagVorticity;

  // Tag buffer size
  int  m_tagBufferSize;
  bool m_tagBufferSizeSet;

  // Initial dt multiplier
  Real m_initialDtMultiplier;
  bool m_initialDtMultiplierSet;

  // Pointer to the object that supplies the physics needed by the integrator
  // (used as a factory here)
  const MOLPhysics* m_molPhysics;

private:

  // Disallowed for all the usual reasons
  void operator=(const AMRLevelConsFactory& a_input);

  // Disallowed for all the usual reasons
  AMRLevelConsFactory(const AMRLevelConsFactory& a_input);
};

#include "NamespaceFooter.H"

#endif
