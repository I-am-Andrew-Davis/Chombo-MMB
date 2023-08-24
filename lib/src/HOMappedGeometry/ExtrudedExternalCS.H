#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EXTRUDEDEXTERNALCS_H_
#define _EXTRUDEDEXTERNALCS_H_

#include "ExternalCS.H"

#include "NamespaceHeader.H"

#define USE_SPLINE

/// ExtrudedExternalCS takes a face from a block given by a .cgns file
/// and then extrudes that face into a separate block
/// This currently only has the infastructor to handle the Channel Separation case

class ExtrudedExternalCS : public ExternalCS
{
public:
          /**
             \name Constructors, destructors and defines
          */

          ///
          /** Default constructor.
           */
 ExtrudedExternalCS(const ReadCGNS& a_cgns,
                    const RealVect& a_dX,
                    const RealVect& a_scale,
                    const ProblemDomain& a_domain,
                    const std::string a_file_name,
                    const std::string& a_zoneName,
                    const int a_constIdx,
                    const IntVect& a_evalGhost,
                    const int a_startProc,
                    const int a_numProc,
                    const int a_blkIdx,
                    const RealVect a_extDx,
                    const IntVect a_extCell);

 ///
 /**
    Destructor */
 virtual ~ExtrudedExternalCS();
 //reads in a grid from a cgns file
 void readGridFile(const ReadCGNS& a_cgns);

protected:
 const int m_blkIdx; //correct number to loop through for loop to define grids
 const RealVect m_extDx; // contains the spacing and which direction to extrude
 const IntVect m_extCell; // the number of cells to extrude
 CHArray<Real, SpaceDim +1> m_extrudedGridNodes; //new gridnodes after extrusion
};

#include "NamespaceFooter.H"

#endif
