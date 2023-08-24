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
 * \file MultiBlockRegions.H
 *
 * \brief Defines regions used by various multiblock registers
 *
 *//*+*************************************************************************/

#ifndef _MULTIBLOCKREGIONS_H_
#define _MULTIBLOCKREGIONS_H_

//----- Standard Library -----//

//----- Chombo Library -----//

#include "Vector.H"
#include "DisjointBoxLayout.H"
#include "LayoutData.H"
#include "RemoteCopier.H"
#include "CH_BoxKeys.H"
#include "CH_Hash.H"

#include "NamespaceHeader.H"

//--Forward declarations from Chombo

class MultiBlockCoordSys;
class LevelFluxRegister;


/*******************************************************************************
 */
///  MultiBlockRegions
/**
 *   This class defines data at block boundaries that connect to another block.
 *   It is intended to replace definitions in:
 *   <ul>
 *     <li> BlockBaseRegister      - flux averages (no AMR interface present)
 *     <li> MultiBlockFluxRegister - fine-to-coarse flux correction in cells
 *     <li> MultiBlockFaceRegister - fine-to-coarse flux replacement (for NtX)
 *   </ul>
 *
 *//*+*************************************************************************/

class MultiBlockRegions
{


/*==============================================================================
 * Public types
 *============================================================================*/

public:

  /// Identification of a face (co-dimension 1) in a layout
  struct FaceTag
  {

//--Types

    using key_type = CH_BoxKeys::DataIndexInt;

//--Constructors, etc.

    /// Default constructor
    FaceTag()
      :
      m_idxFace(-1),
      m_srcDidx()
      { }
    /// Constructor
    FaceTag(const int            a_dir,
            const Side::LoHiSide a_side,
            const DataIndex&     a_srcDidx)
      :
      m_idxFace(indexFace(a_dir, a_side)),
      m_srcDidx(a_srcDidx)
      { }
    /// Constructor
    FaceTag(const int        a_idxFace,
            const DataIndex& a_srcDidx)
      :
      m_idxFace(a_idxFace),
      m_srcDidx(a_srcDidx)
    { }
    // All other copy, move, assignment, destructor are synthesized
    FaceTag(const FaceTag&) = default;
    FaceTag(FaceTag&&) = default;
    FaceTag& operator=(const FaceTag&) = default;
    FaceTag& operator=(FaceTag&&) = default;

//--Member functions

    /// Return key for using this in a map
    key_type key() const
      { return key_type(m_srcDidx, m_idxFace); }
    /// Get dir
    int dir() const
      { return m_idxFace % SpaceDim; }
    // Get side
    Side::LoHiSide side() const
      { return static_cast<Side::LoHiSide>(m_idxFace/SpaceDim); }

//--Static member functions

    /// Get an index for the face
    static int indexFace(const int a_dir, const Side::LoHiSide a_side)
      { return a_side*SpaceDim + a_dir; }
    /// Get dir from a face index
    static int dir(const int a_face)
      { return a_face % SpaceDim; }
    /// Get side from a face index
    static Side::LoHiSide side(const int a_face)
      { return static_cast<Side::LoHiSide>(a_face/SpaceDim); }

//--Data

    int m_idxFace;                    ///< Encodes side and direction of face
                                      ///< Dir has unit stride
    DataIndex m_srcDidx;              ///< Corresponding data index in the
                                      ///< source layout
  };


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Default constructor
  MultiBlockRegions()
    { }

  /// Destructor
  ~MultiBlockRegions() = default;

  // Copy and move prohibited due to storing LayoutData
  MultiBlockRegions(const MultiBlockRegions&) = delete;
  MultiBlockRegions(MultiBlockRegions&&) = delete;
  MultiBlockRegions& operator=(const MultiBlockRegions&) = delete;
  MultiBlockRegions& operator=(MultiBlockRegions&&) = delete;

  /// Weak construction
  void define(const MultiBlockCoordSys&       a_coordSys,
              const DisjointBoxLayout&        a_CrGrid,
              const DisjointBoxLayout *const  a_FnGrid = nullptr,
              const int                       a_nRef = 0);

  /// Clear the data structures
  void clear();


/*==============================================================================
 * Public member functions
 *============================================================================*/

  /// Access locations for same level operations
  const LayoutData<stc::Vector<Vector<Box>, 2*SpaceDim>>& getLvlLocations()
    const
    { return m_lvlLocations; }

  /// Access locations for same level operations
  const stc::Vector<Vector<Box>, 2*SpaceDim>& getLvlLocations(
    const DataIndex& a_didx) const
    { return m_lvlLocations[a_didx]; }

  /// Access locations for same level operations
  const Vector<Box>& getLvlLocations(const DataIterator&  a_dit,
                                     const int            a_dir,
                                     const Side::LoHiSide a_side) const
    { return m_lvlLocations[a_dit][FaceTag::indexFace(a_dir, a_side)]; }

  /// Add coarse locations from block boundaries to the LevelFluxRegister
  void addCrLocationsToFluxReg(LevelFluxRegister& a_fluxReg) const;

  /// Access coarse locations for AMR interface operations
  const LayoutData<stc::Vector<Vector<Box>, 2*SpaceDim>>& getCrLocations() const
    { return m_crLocations; }

  /// Access coarse locations for AMR interface operations
  const Vector<Box>& getCrLocations(const DataIterator&  a_dit,
                                    const int            a_dir,
                                    const Side::LoHiSide a_side) const
    { return m_crLocations[a_dit][FaceTag::indexFace(a_dir, a_side)]; }

  /// Access fine locations for AMR interface operations
  const LayoutData<stc::Vector<Vector<Box>, 2*SpaceDim>>& getFnLocations() const
    { return m_fnLocations; }

  /// Access fine locations for AMR interface operations
  const Vector<Box>& getFnLocations(const DataIterator&  a_dit,
                                    const int            a_dir,
                                    const Side::LoHiSide a_side) const
    { return m_fnLocations[a_dit][FaceTag::indexFace(a_dir, a_side)]; }

  /// Access level copier (copier is mutable)
  RemoteCopier& lvlRemoteCopier() const
    { return m_lvlRemoteCopier; }

  /// Access AMR face copier (copier is mutable)
  RemoteCopier& amrFaceRemoteCopier() const
    { return m_amrFaceRemoteCopier; }

  /// Access AMR flux copier (copier is mutable)
  RemoteCopier& amrFluxRemoteCopier() const
    { return m_amrFluxRemoteCopier; }

  /// Access the level compact storage layout
  const BoxLayout& lvlStoLayout() const
    { return m_lvlStoLayout; }

  /// Access the mapping from compact to original layout
  const LayoutData<FaceTag>& lvlSto2Origin() const
    { return m_lvlSto2Origin; }

  /// Access the mapping from original to compact layout
  const BoxLayout::UMap_T_DataIndex<FaceTag::key_type>& lvlMapOrigin2Sto() const
    { return m_lvlMapOrigin2Sto; }

  /// Get a data index in the compact layout from location in original layout
  const DataIndex& getLvlStoDidx(const FaceTag& a_faceTag) const
    {
      auto iter = m_lvlMapOrigin2Sto.find(a_faceTag.key());
      CH_assert(iter != m_lvlMapOrigin2Sto.end());
      return iter->second;
    }


/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  LayoutData<stc::Vector<Vector<Box>, 2*SpaceDim>> m_lvlLocations;
                                      ///< Cells to adjust on the same level
                                      ///< Index using FaceTag::indxeFace
  LayoutData<stc::Vector<Vector<Box>, 2*SpaceDim>> m_crLocations;
                                      ///< Cells to adjust on the coarse level
                                      ///< (because of a flux correction where
                                      ///< refinement interfaces align with
                                      ///< block boundaries)
                                      ///< Index using FaceTag::indexFace
  LayoutData<stc::Vector<Vector<Box>, 2*SpaceDim>> m_fnLocations;
                                      ///< Cells on the fine level that are
                                      ///< adjacent to a refinement interface
                                      ///< overlapping a block boundary.  Note,
                                      ///< in classic flux correction, this
                                      ///< information is not used and all boxes
                                      ///< average fluxes on their boundaries
                                      ///< down to the coarsened-fine mesh.
                                      ///< Index using FaceTag::indexFace
  mutable RemoteCopier m_lvlRemoteCopier;
                                      ///< Copiers for moving flux information
                                      ///< across block boundaries on the same
                                      ///< level
  mutable RemoteCopier m_amrFaceRemoteCopier;
                                      ///< Copier for moving flux information
                                      ///< on faces across block boundaries
                                      ///< from a fine level to a coarser level.
                                      ///< For uniqueness, the copy is from
                                      ///< outer cells at source to outer cells
                                      ///< at destination
  mutable RemoteCopier m_amrFluxRemoteCopier;
                                      ///< Copier for moving flux corrections to
                                      ///< cells across block boundaries from a
                                      ///< fine level to a coarser level.  The
                                      ///< copy is from outer cells at source to
                                      ///< inner cells at destination.
  BoxLayout m_lvlStoLayout;           ///< A compact layout for storing data for
                                      ///< level exchanges.  As implied by the
                                      ///< type, this is most likely not
                                      ///< disjoint 
  LayoutData<FaceTag> m_lvlSto2Origin;
                                      ///< Gives the origin FaceTag from the
                                      ///< compact layout.  This layout is not
                                      ///< the same as those for the locations
                                      ///< above
  BoxLayout::UMap_T_DataIndex<FaceTag::key_type> m_lvlMapOrigin2Sto;
                                      ///< A map using a key from the origin
                                      ///< FaceTag to the DataIndex in the
                                      ///< compact layout.  If you use this,
                                      ///< remember that you need to extract the
                                      ///< key from the FaceTag.
  /* Not sure if we need this.  Given a destination data index and region, it
     is relatively easy to find the block boundary and construct the FaceTag.
     See getTransform in class MBAddFaceOp in MultiBlockFaceRegister.cpp for
     example.
  */
  // std::unordered_map<Box, FaceTag, CH_Hash::google_CityHash<Box>>
  // m_amrCrRcvToOrigin;                 ///< A map from the coarse receive outer
  //                                     ///< region of a message to the original
  //                                     ///< DataIndex and face.  Primarily used
  //                                     ///< to identify the side of the block
  //                                     ///< boundary in Copier operators.  This
  //                                     ///< is actually the reverse of
  //                                     ///< m_crLocations except that
  //                                     ///< m_crLocations points to the inner
  //                                     ///< cells.
};

#include "NamespaceFooter.H"

#endif  /* ! defined _MULTIBLOCKREGIONS_H_ */
