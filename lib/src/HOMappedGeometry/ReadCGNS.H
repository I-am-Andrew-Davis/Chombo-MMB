#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _READCGNS_H_
#define _READCGNS_H_

#include <string>
#include <vector>
#include <unordered_map>

#ifdef CH_USE_CGNS
  #ifdef CH_MPI
    #include "pcgnslib.h"
  #else
    #include "cgnslib.h"
  #endif
#endif

#include "CH_Hash.H"
#include "LevelData.H"
#include "CHArray.H"

#include "NamespaceHeader.H"


/*******************************************************************************
 */
///  Interfaces with CGNS library to open a file and read data
/**
 *//*+*************************************************************************/

class ReadCGNS
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

  /// Default constructor
  ReadCGNS();

  /// Constructor
  ReadCGNS(std::string a_file_name,
           int a_verbosity=0);

//--Use synthesized copy, move, assignment, move assignment constructors

  ReadCGNS(const ReadCGNS&) = default;
  ReadCGNS(ReadCGNS&&) = default;
  ReadCGNS& operator=(const ReadCGNS&) = default;
  ReadCGNS& operator=(ReadCGNS&&) = default;

  /// Destructor
  ~ReadCGNS();


/*==============================================================================
 * Public member functions
 *============================================================================*/

#ifdef CH_USE_CGNS
  /// Open a CGNS file
  int openFile(std::string a_fileName);
  
  /// Close the cgns file
  int closeFile();

  /// Read base information
  int readBaseInfo(const int    a_idxBase,
                   std::string& a_baseName,
                   int&         a_cellDim,
                   int&         a_physDim) const;

  /// Get the number of base nodes in the file (no checking of components)
  int numBase() const { return m_numBase; }

  /// Select a CGNS base
  int selectBase(const int a_idxBase, ZoneNameMapToBlkIdx& a_zoneNameMap);

  /// Get the number of zones in the file (no checking of components)
  int numZone() const { return m_numZone; }

  /// Read the zone information
  int readZoneInfo(const int    a_idxZone,
                   std::string& a_zoneName,
                   IntVect&     a_numCell) const;

  /// Find all connections bewteen given zone and neighbors
  int readZoneConnections(const int                  a_idxZone,
                          const ZoneNameMapToBlkIdx& a_zoneNameMap,
                          std::vector<int>&          a_connectedIdxBlk,
                          std::vector<Box>&          a_rangesFrom,
                          std::vector<Box>&          a_rangesTo,
                          std::vector<IntVect>&      a_transformations) const;
    
  /// Read the physical location of grid points into an FArrayBox
  int readZoneCoords(const int a_idxZone, FArrayBox& a_grid) const;

  /// Read the physical locations of grid points into a CHArray
  int readZoneCoords(const int                  a_idxZone,
                     CHArray<Real, SpaceDim+1>& a_grid,
                     IntVect&                   a_numVert) const;


/*==============================================================================
 * Protected member functions
 *============================================================================*/

protected:

  /// Check for CGNS errors
  int errorHandler(const int cgerr) const;
#endif


/*==============================================================================
 * Data members
 *============================================================================*/

protected:

  std::string m_fileName;             ///< File name
  int m_idxFile;                      ///< CGNS index of file (-1 if no open
                                      ///< file)
  int m_idxBase;                      ///< CGNS index of base node (-1 if no
                                      ///< base selected or file not open)
  int m_numBase;                      ///< Number of base nodes (available after
                                      ///< opening a file)
  int m_numZone;                      ///< Number of zones (available after
                                      ///< selecting a base)
  int m_verbosity;                    ///< How loud to be about what is happening
};


/*******************************************************************************
 *
 * Class ReadCGNS: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
// Check for CGNS errors
/** \param[in]  a_cgerr CGNS error index
 *  \return             0 : success
 *                      1 : error
 *//*-----------------------------------------------------------------*/

#ifdef CH_USE_CGNS
inline int
ReadCGNS::errorHandler(const int cgerr) const
{
  if (cgerr != 0)
    {
      pout() << "cgerr number " << cgerr << " : " << "\n\t";
      cg_error_print();
      return 1;
    }
  else
    {
      return 0;
    }
}
#endif

#include "NamespaceFooter.H"

#endif  /* ! defined _READCGNS_H_ */
