#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "ReadCGNS.H"

#ifdef CH_USE_CGNS
  #ifdef CH_MPI
    #include "pcgnslib.h"
  #else
    #include "cgnslib.h"
  #endif
#endif

#include "NamespaceHeader.H"


/*******************************************************************************
 *
 * Class ReadCGNS: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Default constructor
/** 
 *//*-----------------------------------------------------------------*/

ReadCGNS::ReadCGNS()
  :
  m_idxFile(-1),
  m_idxBase(-1),
  m_numBase(0),
  m_numZone(0),
  m_verbosity(0)
{ }

/*--------------------------------------------------------------------*/
//  Constructor
/** \param[in] a_fileName
 *                      The name of the .cgns file to read
 *  \note
 *  <ul>
 *    <li> Software aborts if the file is not opened successfully.
 *         Use the default constructor and explicitly call openFile if
 *         this is not what you want.
 *  </ul>
 *//*-----------------------------------------------------------------*/

ReadCGNS::ReadCGNS(std::string a_fileName,
                   int a_verbosity)
  :
  m_idxFile(-1),
  m_idxBase(-1),
  m_numBase(0),
  m_numZone(0),
  m_verbosity(a_verbosity)
{
  m_verbosity = 3; //BDG HACK
#ifdef CH_USE_CGNS
  int status = 0;
  status += openFile(a_fileName);
  CH_assert(status == 0);
#endif
}

/*--------------------------------------------------------------------*/
//  Destructor
/** 
 *//*-----------------------------------------------------------------*/

ReadCGNS::~ReadCGNS()
{
#ifdef CH_USE_CGNS
  closeFile();
#endif
}


/*==============================================================================
 * Public member functions
 *============================================================================*/

//--All remaining functions are only defined if using CGNS

#ifdef CH_USE_CGNS

/*--------------------------------------------------------------------*/
//  Open a CGNS file
/** \param[in] a_fileName
 *                      The name of the .cgns file to open
 *  \return             0  : success
 *                      !0 : accumulation of cgerr
 *//*-----------------------------------------------------------------*/

int
ReadCGNS::openFile(std::string a_fileName)
{
  CH_TIME("ReadCGNS::openFile");
  int status = 0;
  int cgerr  = 0;
  if (m_idxFile >= 0)
    {
      if (a_fileName == m_fileName)
        {
          return 0;
        }
      status += closeFile();
    }
  m_fileName = a_fileName;
  // Open the CGNS file
#ifdef CH_MPI
  cgerr = cgp_open(m_fileName.c_str(), CG_MODE_READ, &m_idxFile);
#else
  cgerr =  cg_open(m_fileName.c_str(), CG_MODE_READ, &m_idxFile);
#endif
  status += errorHandler(cgerr);

  // Get number of bases
  cgerr = cg_nbases(m_idxFile, &m_numBase);
  status += errorHandler(cgerr);
  return status;
}

/*--------------------------------------------------------------------*/
//  Close the cgns file
/**
 *  \return             0  : success
 *                      !0 : accumulation of cgerr
 *//*-----------------------------------------------------------------*/

int ReadCGNS::closeFile()
{
  CH_TIME("ReadCGNS::closeFile");
  int status = 0;
  if (m_idxFile >= 0)
    {
      int cgerr;
#ifdef CH_MPI
      cgerr = cgp_close(m_idxFile);
#else
      cgerr =  cg_close(m_idxFile);
#endif
      status += errorHandler(cgerr);
    }
  m_fileName.clear();
  m_idxFile = -1;
  m_idxBase = -1;
  m_numBase = 0;
  m_numZone = 0;
  return status;
}

/*--------------------------------------------------------------------*/
//  Read base information
/** 
 *  \return             0  : success
 *                      !0 : accumulation of cgerr
 *//*-----------------------------------------------------------------*/

int
ReadCGNS::readBaseInfo(const int    a_idxBase,
                       std::string& a_baseName,
                       int&         a_cellDim,
                       int&         a_physDim) const
{
  int status = 0;
  int cgerr  = 0;
  CH_assert(m_idxFile > 0);
  CH_assert(a_idxBase > 0 && a_idxBase <= m_numBase);
  char name[33];
  cgerr = cg_base_read(m_idxFile, m_idxBase, name, &a_cellDim, &a_physDim);
  status += errorHandler(cgerr);
  a_baseName = name;
  return status;
}

/*--------------------------------------------------------------------*/
//  Select a CGNS base
/** \param[in]  a_idxBase
 *                      CGNS index of the base
 *  \param[out] a_zoneNameMap
 *                      Name-block index pairs added to the map
 *  \return             0  : success
 *                      !0 : accumulation of cgerr
 *
 *  \note
 *  <ul>
 *    <li> For the map, it is assumed that idxBlk = idxZone - 1
 *    <li> a_zoneNameMap should be cleared or otherwise conditioned
 *         before calling this routine
 *  </ul>
 *//*-----------------------------------------------------------------*/

int
ReadCGNS::selectBase(const int a_idxBase, ZoneNameMapToBlkIdx& a_zoneNameMap)
{
  int status = 0;
  int cgerr  = 0;
  CH_assert(m_idxFile >= 0);
  CH_assert(a_idxBase > 0 && a_idxBase <= m_numBase);
  m_idxBase = a_idxBase;

  std::string baseName;
  int cellDim;
  int physDim;
  status += readBaseInfo(m_idxBase, baseName, cellDim, physDim);
  if (cellDim != SpaceDim)
    {
      MayDay::Error("Base selected in CGNS file that has different dimension "
                    "than SpaceDim");
    }

  // Get number of zones
  cgerr = cg_nzones(m_idxFile, m_idxBase, &m_numZone);
  status += errorHandler(cgerr);

  // Hash each zone name in a lookup map.  This maps zone names to the block
  // index (not the CGNS zone index)
  a_zoneNameMap.reserve(m_numZone);
  std::string zoneName;
  IntVect numCell;
  for (int idxZone = 1; idxZone <= m_numZone; ++idxZone)
    {
      status += readZoneInfo(idxZone, zoneName, numCell);
      const int idxBlk = idxZone - 1;
      auto ins = a_zoneNameMap.insert({ zoneName, idxBlk });
      CH_assert(ins.second);  // Insertion must have succeeded
    }

  return status;
}

/*--------------------------------------------------------------------*/
//  Read the zone information
/** \param[in]  a_idxZone
 *                      CGNS index of the zone
 *  \param[out] a_zoneName
 *                      Name of the zone
 *  \param[out] a_numCell
 *                      Number of cells in the zone
 *  \return             0  : success
 *                      !0 : accumulation of cgerr
 *//*-----------------------------------------------------------------*/

int
ReadCGNS::readZoneInfo(const int    a_idxZone,
                       std::string& a_zoneName,
                       IntVect&     a_numCell) const
{
  int status = 0;
  int cgerr  = 0;
  CH_assert(m_idxFile >= 0);
  CH_assert(m_idxBase > 0 && m_idxBase <= m_numBase);
  CH_assert(a_idxZone > 0 && a_idxZone <= m_numZone);
  char name[33];
  cgsize_t size[3*SpaceDim];
  cgerr = cg_zone_read(m_idxFile, m_idxBase, a_idxZone, name, size);
  status += errorHandler(cgerr);
  a_zoneName = name;
  a_numCell = stc::Vector<cgsize_t, SpaceDim>(size + SpaceDim);
  return status;
}

/*--------------------------------------------------------------------*/
//  Find all connections bewteen given zone and neighbors
//  The ranges are based around a domain box with origin of zero
/** \param[in]  a_idxZone
 *                      CGNS index of the zone
 *  \param[in]  a_zoneNameMap
 *                      Map of zone-name to block indices
 *  \param[out] a_connectedIdxZone
 *                      CGNS zone index of neighbors
 *  \param[out] a_rangesFrom
 *                      Node centered box describing range (face) of
 *                      connectivity from a_idxZone
 *  \param[out] a_rangesTo
 *                      Node-centered box describing range (face) of
 *                      connetivity of neighbor zone
 *  \param[out] a_transformations
 *                      Transformation to get from a_idxZone to the
 *                      neighbor zone
 *  \return             0  : success
 *                      !0 : accumulation of cgerr
 *
 *  \note
 *  <ul>
 *    <li> From the map, it is assumed that idxZone = idxBlk + 1
 *  </ul>
 *//*-----------------------------------------------------------------*/

int
ReadCGNS::readZoneConnections(
  const int                  a_idxZone,
  const ZoneNameMapToBlkIdx& a_zoneNameMap,
  std::vector<int>&          a_connectedIdxBlk,
  std::vector<Box>&          a_rangesFrom,
  std::vector<Box>&          a_rangesTo,
  std::vector<IntVect>&      a_transformations) const
{
  int status = 0;
  int cgerr  = 0;
  CH_assert(m_idxFile >= 0);
  CH_assert(m_idxBase > 0 && m_idxBase <= m_numBase);
  CH_assert(a_idxZone > 0 && a_idxZone <= m_numZone);

  // Get number of connections
  int numConnect;
  cgerr = cg_n1to1(m_idxFile, m_idxBase, a_idxZone, &numConnect);
  status += errorHandler(cgerr);
  if (m_verbosity > 2) pout() << "Zone " << a_idxZone << " has " << numConnect << " connections" << std::endl;

  // Resize vectors
  a_connectedIdxBlk.resize(numConnect);
  a_rangesFrom.resize(numConnect);
  a_rangesTo.resize(numConnect);
  a_transformations.resize(numConnect);

  // loop over each connection
  stc::Vector<cgsize_t, 2*SpaceDim> zoneRange, donorRange;
  IntVect trfm;
  char connectName[33], donorName[33];
  for (int idxConnect = 1; idxConnect <= numConnect; ++idxConnect)
    {
      // Index to vectors
      auto idx = idxConnect - 1;
      // Read the connectivity
      cgerr = cg_1to1_read(m_idxFile, m_idxBase, a_idxZone, idxConnect,
                           connectName, donorName,
                           zoneRange.dataPtr(), donorRange.dataPtr(),
                           a_transformations[idx].dataPtr());
      status += errorHandler(cgerr);
      // Straightforward copy from 2*SpaceDim vectors to IntVects
      IntVect zr1 = zoneRange;
      IntVect zr2 = stc::Vector<cgsize_t, SpaceDim>(zoneRange, SpaceDim);
      IntVect dr1 = donorRange;
      IntVect dr2 = stc::Vector<cgsize_t, SpaceDim>(donorRange, SpaceDim);
      if (m_verbosity > 2)
        {
          pout() << "Zone " << a_idxZone << " conn " << idxConnect
                 << " : connections "
                 << zr1 << zr2
                 << dr1 << dr2 << std::endl;
        }
      // Setup the correct box ranges from the zones
      a_rangesFrom[idx].define(min(zr1, zr2) - IntVect_unit,
                               max(zr1, zr2) - IntVect_unit,
                               IndexType::TheNodeType());
      a_rangesTo[idx]  .define(min(dr1, dr2) - IntVect_unit,
                               max(dr1, dr2) - IntVect_unit,
                               IndexType::TheNodeType());
      if (m_verbosity > 2)
        {
          pout() << "Zone " << a_idxZone << " conn " << idxConnect
                 << " connections "
                 << "from: " <<  a_rangesFrom[idx]
                 << " to: " <<  a_rangesTo[idx]
                 << std::endl;
        }
      // Find the zone with a matching donorName
      auto iter = a_zoneNameMap.find(donorName);
      if (iter == a_zoneNameMap.end())
        {
          MayDay::Error("Connected donor zone not found in map of names");
        }
      a_connectedIdxBlk[idx] = iter->second;
      if (m_verbosity > 2)
        {
          pout() << "Zone " << a_idxZone << " is connected to zone " << a_connectedIdxBlk[idx]+1 << std::endl;
        }
    }

  return status;
}

/*--------------------------------------------------------------------*/
//  Read the physical location of grid points into an FArrayBox
/** \param[in]  a_idxZone
 *                      CGNS index of the zone
 *  \param[out] a_grid  Contains SpaceDim coordinates
 *  \return             0  : success
 *                      !0 : accumulation of cgerr
 *//*-----------------------------------------------------------------*/

int
ReadCGNS::readZoneCoords(const int a_idxZone, FArrayBox& a_grid) const
{
  CH_TIME("ReadCGNS::readZoneCoords(FArrayBox)");
  int status = 0;
  int cgerr  = 0;
  CH_assert(m_idxFile >= 0);
  CH_assert(m_idxBase > 0 && m_idxBase <= m_numBase);
  CH_assert(a_idxZone > 0 && a_idxZone <= m_numZone);

  std::string zoneName;
  IntVect numCell;
  status += readZoneInfo(a_idxZone, zoneName, numCell);
  const IntVect numVert = numCell + 1;

  // Set the FArrayBox to appopriate size
  Box gridBox(IntVect_zero, numCell, IndexType::TheNodeType());
  a_grid.resize(gridBox, SpaceDim);

  stc::Vector<cgsize_t, SpaceDim> rmin, rmax, memRmin, memRmax, memDims;
  rmin = 1;
  rmax = numVert;
  memRmin = 1;
  memRmax = numVert;
  memDims = numVert;

  static constexpr const char* coordName[] =
    { "CoordinateX", "CoordinateY", "CoordinateZ" };
  // Number of coordinates
  int numCoord;
  cgerr = cg_ncoords(m_idxFile, m_idxBase, a_idxZone, &numCoord);
  status += errorHandler(cgerr);
  CH_assert(numCoord >= SpaceDim);

  for (const int dir : EachDir)
    {
#ifdef CH_MPI
      // Get the coordinate index
      char name[33];
      CGNS_ENUMT(DataType_t) dataType;
      int idxCoord;
      for (idxCoord = 1; idxCoord <= numCoord; ++idxCoord)
        {
          cgerr = cg_coord_info(m_idxFile, m_idxBase, a_idxZone, idxCoord,
                                &dataType, name);
          status += errorHandler(cgerr);
          if (strncmp(name, coordName[dir], 33) == 0) break;
        }
      CH_assert(idxCoord <= numCoord);  // Coord must be found

      // Read the coordinate
      cgerr = cgp_coord_general_read_data(m_idxFile, m_idxBase, a_idxZone,
                                          idxCoord,
                                          rmin.dataPtr(), rmax.dataPtr(),
                                          CGNS_REAL, SpaceDim,
                                          memDims.dataPtr(),
                                          memRmin.dataPtr(), memRmax.dataPtr(),
                                          &a_grid(IntVect_zero, dir));
#else
      // Read the coordinate
      cgerr = cg_coord_general_read(m_idxFile, m_idxBase, a_idxZone,
                                    coordName[dir],
                                    rmin.dataPtr(), rmax.dataPtr(),
                                    CGNS_REAL, SpaceDim,
                                    memDims.dataPtr(),
                                    memRmin.dataPtr(), memRmax.dataPtr(),
                                    &a_grid(IntVect_zero, dir));
#endif
      status += errorHandler(cgerr);
    }
  return status;
}

/*--------------------------------------------------------------------*/
//  Read the physical locations of grid points into a CHArray
/** \param[in]  a_idxZone
 *                      CGNS index of the zone
 *  \param[out] a_grid  Contains SpaceDim coordinates
 *  \param[out] numVert Number of vertices read
 *  \return             0  : success
 *                      !0 : accumulation of cgerr
 *//*-----------------------------------------------------------------*/

int
ReadCGNS::readZoneCoords(const int                  a_idxZone,
                         CHArray<Real, SpaceDim+1>& a_grid,
                         IntVect&                   a_numVert) const
{
  CH_TIME("ReadCGNS::readZoneCoords(CHArray)");
  int status = 0;
  int cgerr  = 0;
  CH_assert(m_idxFile >= 0);
  CH_assert(m_idxBase > 0 && m_idxBase <= m_numBase);
  CH_assert(a_idxZone > 0 && a_idxZone <= m_numZone);

  std::string zoneName;
  IntVect numCell;
  status += readZoneInfo(a_idxZone, zoneName, numCell);
  a_numVert = numCell + 1;

  // Set the CHArray to appopriate size
  a_grid.define(a_numVert, SpaceDim);

  stc::Vector<cgsize_t, SpaceDim> rmin, rmax, memRmin, memRmax, memDims;
  rmin = 1;
  rmax = a_numVert;
  memRmin = 1;
  memRmax = a_numVert;
  memDims = a_numVert;

  static constexpr const char* coordName[] =
    { "CoordinateX", "CoordinateY", "CoordinateZ" };
  // Number of coordinates
  int numCoord;
  cgerr = cg_ncoords(m_idxFile, m_idxBase, a_idxZone, &numCoord);
  status += errorHandler(cgerr);
  CH_assert(numCoord >= SpaceDim);

  for (const int dir : EachDir)
    {
#ifdef CH_MPI
      // Get the coordinate index
      char name[33];
      CGNS_ENUMT(DataType_t) dataType;
      int idxCoord;
      for (idxCoord = 1; idxCoord <= numCoord; ++idxCoord)
        {
          cgerr = cg_coord_info(m_idxFile, m_idxBase, a_idxZone, idxCoord,
                                &dataType, name);
          status += errorHandler(cgerr);
          if (strncmp(name, coordName[dir], 33) == 0) break;
        }
      CH_assert(idxCoord <= numCoord);  // Coord must be found

      // Read the coordinate
      cgerr = cgp_coord_general_read_data(m_idxFile, m_idxBase, a_idxZone,
                                          idxCoord,
                                          rmin.dataPtr(), rmax.dataPtr(),
                                          CGNS_REAL, SpaceDim,
                                          memDims.dataPtr(),
                                          memRmin.dataPtr(), memRmax.dataPtr(),
                                          &a_grid(IntVect_zero, dir));
#else
      // Read the coordinate
      cgerr = cg_coord_general_read(m_idxFile, m_idxBase, a_idxZone,
                                    coordName[dir],
                                    rmin.dataPtr(), rmax.dataPtr(),
                                    CGNS_REAL, SpaceDim,
                                    memDims.dataPtr(),
                                    memRmin.dataPtr(), memRmax.dataPtr(),
                                    &a_grid(IntVect_zero, dir));
#endif
      status += errorHandler(cgerr);
    }
  return status;
}

#endif  /* defined CH_USE_CGNS */

#include "NamespaceFooter.H"
