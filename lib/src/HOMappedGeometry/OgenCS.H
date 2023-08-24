#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _OGENCS_H_
#define _OGENCS_H_

#include "NewFourthOrderCoordSys.H"
#include <functional>
#include <memory>
#include "ComputeCache.H"

#include "readOgen.H"

#include "NamespaceHeader.H"

/*******************************************************************************
 */
/// Retrieve a coordinate mapping from pregenerated source
/**
 *  Current file types are:
 *    CGNS       - read a single block in for discreate mapping
 *               - optional spline construction allows for AMR
 *    Ogen       - Interface with Overture which provides an
 *                 analytic grid generator
 *
 *  For the dircreate grids we have two options, either use them as is
 *  and then use FD to calculate grid metrics. AMR does not work with
 *  this method. The other option is to use a spline interpolation,
 *  which provides all grid metrics and allows AMR. The spline interpolation
 *  class relies on features of c++14
 ******************************************************************************/
class OgenCS : public NewFourthOrderCoordSys
{
public:
  /// Constructor
  OgenCS(const RealVect& a_dXi, // Most of Chombo says d_x when meaning d_xi
         const RealVect& a_scale,
         const ProblemDomain& a_domain,
         const std::string a_file_name,
         const int a_idxBlk = 0,
         const IntVect& a_evalGhost = IntVect::Zero);

  /// Destructor
  virtual ~OgenCS();
  
  /// Given coordinate in mapped space, return its location in real space
  // Cache the values for faster evaluation
  virtual RealVect realCoord(const RealVect& a_Xi) const;

  /// Given coordinate in real space, return its location in mapped space
  virtual RealVect mappedCoord(const RealVect& a_x) const;
  
  /// Given coordinate in real space, return its location in mapped space
  /// Using an initial guess makes this much faster.
  RealVect mappedCoord(const RealVect& a_x,
                       const RealVect& a_xiGuess) const;

  /// Calculate the derivative of each coordinate vector
  virtual Real dXdXi(const RealVect& a_Xi, int a_dirX, int a_dirXi) const;
  
  /// A coarse search for mapping inverse in a region
  // this is a brute force search and should be used with caution
  // does not adjust for periodic domains
  RealVect globalInverseSearch(const RealVect& a_x,
                               const RealVect& a_dXi,
                               const Box& a_search) const;
  
  /// A coarse search for mapping inverse in a region
  // A targeted search, starting in the region around a_xiEst
  // does not adjust for periodic domains
  RealVect inverseSearch(const RealVect& a_x,
                         const RealVect& a_dXi,
                         const Box& a_search,
                         const RealVect& a_XiEst) const;

  /// Check if a real point is in the mapping problem domain
  bool inBlockDomain(const RealVect& a_x,
                     RealVect& a_xiEst) const;

  /// enable 
  void enableCache();

  /// disable and clear
  void disableCache();
  
protected:
  std::string m_grid_file;              ///< Name of the file to read
  int m_idxZone;                        ///< Zone to read from
  ProblemDomain m_domain;               ///< The problem domain
  IntVect m_gridSize;                   ///< The number of points in a grid
  RealVect m_scale;                     ///< Uniformly scale the grid
  RealVect m_scaleXi;                   ///< Scale from grid points to domain size
  RealVect m_xiOrigin;                  ///< The location Xi starts from
  IntVect m_evalGhost;                  ///< Number of ghost cells to consider extrapolation to. No guarantees in ghost cells quality, but if you require ghost cells that is your own problem!
  std::vector<RealVect> m_domainOffset; ///< X length along Xi axes

  // if using Ogen
  readOgen ogen;

  /// Setup caching of mapping transformations, since they can be expensive
  bool m_useCache;

  /// setup cache type
  using coordCache = ComputeCacheUMap<
    IntVect,                                 // Key (index in cell+face space)
    RealVect,                                // Cached value
    std::function<RealVect(const IntVect&)> >;
  std::unique_ptr<coordCache> m_coordCache;

  
  /// Read and store grid data
  void readGridFile();

  /// Periodic adjustment
  IntVect coordAdj(RealVect& a_Xi) const;

  /// Given coordinate in mapped space, return its location in real space
  // evaluates the values, no caching
  RealVect realCoordRaw(const RealVect& a_Xi) const;
  
  /// Given coordinate in mapped space, return its location in real space
  // does not adjust for periodic domains
  RealVect realCoordNonAdjst(const RealVect& a_Xi) const;
};

/*******************************************************************************
 */
///  Factory object to create an external coordinate system
/**
 ******************************************************************************/
class
OgenCSFactory : public NewCoordSysFactory
{

public:
  /// Constructor
  OgenCSFactory(const std::string& a_grid_file,
                    const RealVect& a_scale);

  /// Return a pointer to the OgenCS object
  virtual NewCoordSys* getCoordSys(const ProblemDomain& a_domain,
                                   const RealVect& a_dx) const;

protected:
  std::string m_grid_file;              ///< name of cgns file
  RealVect m_scale;                     ///< grow or shrink to match specified size
};

#include "NamespaceFooter.H"

#endif
