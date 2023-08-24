#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MULTIBLOCKFACEREGISTER_H_
#define _MULTIBLOCKFACEREGISTER_H_

#include "MultiBlockCoordSys.H"
#include "MultiBlockRegions.H"
#include "NamespaceHeader.H"

//! \class MultiBlockFaceRegister
//! This class handles multi-block AMR refluxing.
//! The face register differs from the flux-register in that face-values
//! themselves are corrected rather than cell-averages based on a
//! divergence.
class MultiBlockFaceRegister
{
protected:

  /// Define status
  enum FaceRegDefineStatus : unsigned
  {
    FaceRegUndefined      =        0,
    FaceRegDefined        = (1 << 0), ///< Indicates define() was called
    FaceRegIncFine        = (1 << 1), ///< IncrementFine was called
  };

public:

  //! Creates a face register that handles multi-block refluxing.
  //! This capability requires a multi-block coordinate system, which
  //! handles block connectivity as well as coordinate mappings.
  MultiBlockFaceRegister(const MultiBlockCoordSys *const a_coordSys,
                         const MultiBlockRegions&        a_mbRegions);

  //! Destructor.
  ~MultiBlockFaceRegister() = default;

  // Generally prohibit all others
  MultiBlockFaceRegister(const MultiBlockFaceRegister&) = delete;
  MultiBlockFaceRegister(MultiBlockFaceRegister&&) = delete;
  MultiBlockFaceRegister& operator=(const MultiBlockFaceRegister&) = delete;
  MultiBlockFaceRegister& operator=(MultiBlockFaceRegister&&) = delete;

  /// Define everything (overrides base)
  void define(const DisjointBoxLayout& a_FnGrid,
              int                      a_nRefine,
              int                      a_nComp,
              bool                     a_scaleFineFaces = true);

  /// Define everything (overrides base)
  void define(const DisjointBoxLayout& a_FnGrid,
              IntVect                  a_nRefine,
              int                      a_nComp,
              bool                     a_scaleFineFaces = true);

  /// Increment with the fine fluxes
  void incrementFine(const LevelData<FluxBox>& a_fnFlux,
                     const Real                a_scale,
                     const Interval&           a_srcIntv,
                     const Interval&           a_dstIntv);

  // New reface method (analagous to reflux)
  void reface(LevelData<FluxBox>&       a_uCoarse,
              const Interval&           a_coarse_interval,
              const Interval&           a_face_interval,
              const Real                a_scale,
              const bool                a_directional = true,
              const LevelData<FluxBox>* a_N = nullptr,
              const RealVect&           a_dx = RealVect_zero);

  /// Has the class been defined?
  bool isDefined()
    { return (m_isDefined & FaceRegDefined); }

private:

  LevelData<FArrayBox> m_fineFlux;
  const MultiBlockCoordSys *const m_coordSys;
                                      ///< Coordinate system from the coarse
                                      ///< level
  const MultiBlockRegions& m_mbRegions;
                                      ///< Identifies intra- and inter-level
                                      ///< connectivity between blocks on a grid
                                      ///< hierarchy.  Typically expected to
                                      ///< be defined in LevelGridMetrics.
  IntVect m_nRefine;                  ///< Refinement ratio between levels.
  unsigned m_isDefined;               ///< Flags from FaceRegDefineStatus
  bool m_scaleFineFaces;              ///< Whether or not to scale fine fluxes
                                      ///< by the refinement ratio
};

#include "NamespaceFooter.H"
#endif
