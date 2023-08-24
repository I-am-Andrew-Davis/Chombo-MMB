#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EXTERNALCS_H_
#define _EXTERNALCS_H_

#include <functional>
#include <memory>

#include "NewFourthOrderCoordSys.H"
#include "CHArray.H"
#include "ComputeCache.H"
#include "nanoflann.H"

// This only works with splines, and by association StcVector and c++14
#define USE_SPLINE
#ifdef USE_SPLINE
#include "BSplineInterp.H"
#endif

#ifdef CH_USE_CGNS
#include "ReadCGNS.H"
#endif

#include "NamespaceHeader.H"

//--Forward declarations
class ReadCGNS;

//-- Type for use with nanoflann
struct CHArrayKDadaptor
{
  const CHArray<Real, SpaceDim+1 >& m_data; // A const ref to the data set origin

  /// The constructor that sets the data set source
  CHArrayKDadaptor(const CHArray<Real, SpaceDim+1 >& a_data)
    :
    m_data(a_data)
    { }

  /// CRTP helper method
  inline const CHArray<Real, SpaceDim+1 >& data() const
    {
      return m_data;
    }
  
  // Must return the number of data points
  inline size_t kdtree_get_point_count() const
    {
      auto sizeV(m_data.sizeVector());
      sizeV[SpaceDim] = 1; // don't count the number of components
      auto size = stc::product(sizeV);
      CH_assert(size >= 0);
      return size;
    }

  // Returns the dim'th component of the idx'th point in the class:
  inline Real kdtree_get_pt(const size_t idx, const size_t dim) const
    {
      CH_assert(dim >= 0);
      CH_assert(dim < SpaceDim);
      auto leng = kdtree_get_point_count();
      CH_assert(idx < leng);
      return m_data.begin()[idx + dim*leng];
    }

  // method for accessing elements, not required as part of adaptor
  inline RealVect get_pt(const size_t idx) const
    {
      RealVect sol;
      for (auto dir : EachDir)
        {
          sol[dir] = kdtree_get_pt(idx, dir);
        }
      return sol;
    }
  
  // Optional bounding-box computation: return false to determine to a standard bbox computation loop.
  template <class BBOX>
  bool kdtree_get_bbox(BBOX& /*bb*/) const { return false; }
};

/*******************************************************************************
 */
/// Retrieve a coordinate mapping from pregenerated source
/**
 *  Current file types are:
 *    CGNS       - read a single block in for discreate mapping
 *               - optional spline construction allows for AMR
 *
 *  Once the grid is read in, a spline interpolation is performed which provides
 *  all grid metrics at arbitrary point and allows AMR. The spline interpolation
 *  class relies on features of c++14
 ******************************************************************************/
class ExternalCS : public NewFourthOrderCoordSys
{
public:
  /// Constructor
  ExternalCS(const ReadCGNS& a_cgns,
             const RealVect& a_dXi, // Most of Chombo says d_x when meaning d_xi
             const RealVect& a_scale,
             const ProblemDomain& a_domain,
             const std::string a_file_name,
             const std::string& a_zoneName,
             const int a_idxBlk = 0,
             const IntVect& a_evalGhost = IntVect::Zero,
             const int a_startProc = 0,
             const int a_numProc = 1);

  /// Destructor
  virtual ~ExternalCS();
  
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

  /// A reasonable block name (here taken from the cgns file)
  virtual std::string name(const int a_idxBlk) const override
    { return m_zoneName; }

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

  /// A coarse search for mapping inverse in the problem domain
  // This is very cheap to evaluate, but does require some setup cost
  // does not adjust for periodic domains
  RealVect inverseSearch(const RealVect& a_x) const;

  /// Check if a real point is in the mapping problem domain
  bool inBlockDomain(const RealVect& a_x,
                     RealVect& a_xiEst) const;

  /// Get Xi of the cell center point nearest to X
  RealVect nearestDomainCell(const RealVect& a_x) const;

  /// enable 
  void enableCache();

  /// disable and clear
  void disableCache();

protected:
  std::string m_grid_file;              ///< Name of the file to read
  std::string m_zoneName;               ///< Zone name as given by the cgns file
  int m_idxZone;                        ///< Zone to read from
  ProblemDomain m_domain;               ///< The problem domain
  IntVect m_gridSize;                   ///< The number of points in a grid
  RealVect m_scale;                     ///< Uniformly scale the grid
  RealVect m_scaleXi;                   ///< Scale from grid points to domain size
  RealVect m_xiOrigin;                  ///< The location Xi starts from
  RealVect m_baseDx;                    ///< Computational grid size on the coarsest level
  int m_refRatio;                       ///< Refinement ratio from the coarsest level
  IntVect m_evalGhost;                  ///< Number of ghost cells to consider extrapolation to. No guarantees in ghost cells quality, but if you require ghost cells that is your own problem!
  int m_startProc;                      ///< Processor to start parallel work with
  int m_numProc;                        ///< Number of processors to use for parallel work

  std::vector<RealVect> m_domainOffset; ///< X length along Xi axes
  
  CHArray<Real, SpaceDim+1> m_gridNodes;
                                      ///< \vect{x} locations from external

#ifdef USE_SPLINE
  BSplineInterp<SpaceDim, 5> m_gridInterp;
                                      ///< solve for metric interpolation
#endif
  /// Setup caching of mapping transformations, since they can be expensive
  bool m_useCache;

  /// setup cache type
  using coordCache =
    ComputeCacheUMap<
    IntVect,                                 // Key (index in cell+face space)
    RealVect,                                // Cached value
    std::function<RealVect(const IntVect&)> >;
  
  std::unique_ptr<coordCache> m_coordCache;

  /// Setup KD-tree searching for fast mapping inverses
  // Type to be used for KD-trees from CHArray data
  using KDtree =
    nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<Real, CHArrayKDadaptor >,
    CHArrayKDadaptor,
    SpaceDim>;
  
  std::unique_ptr<CHArrayKDadaptor> m_KDadaptor;
  
  std::unique_ptr<KDtree> m_KDtree;

  /// Read and store grid data
  void readGridFile(const ReadCGNS& a_cgns);

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
ExternalCSFactory : public NewCoordSysFactory
{

public:
  /// Constructor
  ExternalCSFactory(const std::string& a_grid_file,
                    const RealVect& a_scale);

  /// Return a pointer to the ExternalCS object
  virtual NewCoordSys* getCoordSys(const ProblemDomain& a_domain,
                                   const RealVect& a_dx) const;

protected:
  std::string m_grid_file;              ///< name of cgns file
  RealVect m_scale;                     ///< grow or shrink to match specified size
};

#include "NamespaceFooter.H"
#endif //#ifndef _EXTERNALCS_
