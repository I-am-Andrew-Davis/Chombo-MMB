#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EXTERNALMULTIBLOCKCS_H_
#define _EXTERNALMULTIBLOCKCS_H_

#include <string>
#include <unordered_map>

#include "CH_Hash.H"
#include "MultiBlockCoordSys.H"
#include "NamespaceHeader.H"

//--Forward declarations

class ReadCGNS;

/// MultiBlockCoordSys built from reading in external grids
class ExternalMultiBlockCS : public MultiBlockCoordSys
{


/*==============================================================================
 * Types
 *============================================================================*/

public:

  // A map from zone names to block indices
  using ZoneNameMapToBlkIdx =
    std::unordered_map<std::string, int, CH_Hash::google_CityHash<std::string>>;


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/

  ///
  /** Default constructor.
   */
  ExternalMultiBlockCS();

  ///
  /** Destructor.
   */
  virtual ~ExternalMultiBlockCS();

  virtual void
  define(const ProblemDomain& a_levelDomain,
         const RealVect&      a_dx);

  /*@}*/

  /// block mapping conversion function
  /** given xi_k in block k, returns which block (k_valid)
      contains that location in its valid region and its
      location in that block's computational coordinate space.
  */
  virtual void
  blockRemapping(RealVect&            a_xi_valid,
                 int&                 a_n_valid,
                 bool&                a_validExists,
                 RigidTransformation& a_extraDispl,
                 const RealVect&      a_xiSrc,
                 const IntVect&       a_iSrc,
                 int                  a_nSrc) const override;
    
  virtual Vector<RealVect>
  displacements(const Vector<RealVect>& a_dstCoords,
                const Vector<int>&      a_dstBlocks,
                const RealVect&         a_srcCoords,
                int                     a_srcBlock) const;

  // /// Have separate volume flux on each block.  Normally does nothing, but double-cartesian coordinate system needs it.
  // virtual void separateVolFlux(LevelData<FluxBox>& a_flux) const;

  /// set which grid to read from
  void setupForExternalGrid(const std::string a_file_name,
                            const RealVect a_scale);
  
protected:

  ZoneNameMapToBlkIdx m_zoneNameMap;  ///< Map of zone names to block indices
  int m_numBlocks;
  std::string m_fileName;             ///< Name of external file to read from
  RealVect m_scale;                   ///< Resizing of external grid
  bool m_gotExternalInfo;             ///< Flag if external source is read
//**Note: this class does not yet handle conforming connections that are not
//**      periodic.  Periodic directions must be conforming.
  std::array<bool, SpaceDim> m_periodicDir;
                                      ///< Flag directions which are periodic
  bool m_followGridLinesOffBC = true; ///< For multiblock ghost cells extending
                                      ///< in Dir. A, where on block there is a
                                      ///< boundary normal to A, where in the
                                      ///< adjacent block there is no boundary
                                      ///< normal to A, follow grids lines in
                                      ///< the adjacent block if the valid cell
                                      ///< is not already in the adjacent block.
                                      ///< This situation a convex corner with
                                      ///< three blocks (in 2-D).  Set to false
                                      ///< for 4th order MB exchange.
  int m_followGridLinesNumNormal = 1; ///< Number of ghost normal to the
                                      ///< boundary in which to follow grid
                                      ///< lines.
  int m_followGridLinesNumGhost = 2;  ///< Number of ghosts tangential to the
                                      ///< boundary in which to follow grid
                                      ///< lines (same as number of ghosts
                                      ///< outside the block with the boundary)

  // define m_boundaries
  void defineBoundaries(const ReadCGNS& a_cgns);

private:

  void operator = (const ExternalMultiBlockCS&)
    {
    }

  ExternalMultiBlockCS(const ExternalMultiBlockCS&)
    {
    }
};


class ExternalMultiBlockCSFactory : public MultiBlockCoordSysFactory
{
  
public:

  ExternalMultiBlockCSFactory(const std::string a_fileName,
                              const RealVect a_scale,
                              const int a_verbosity=0)
    :
    m_fileName(a_fileName),
    m_scale(a_scale),
    m_verbosity(a_verbosity)
    { }

  virtual ~ExternalMultiBlockCSFactory()
    { }
  
  virtual MultiBlockCoordSys* getCoordSys(const ProblemDomain& a_levelDomain,
                                          const RealVect& a_dx) const;
  protected:

  std::string m_fileName;
  RealVect m_scale;
  int m_verbosity;
};

#include "NamespaceFooter.H"

#endif
