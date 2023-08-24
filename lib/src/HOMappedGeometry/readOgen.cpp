#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifdef CH_USE_OGEN
#include "readOgen.H"
#include "CH_Timer.H"

//----- Ogen lib -----//

#include "Overture.h"
#include <A++.h>


#include "NamespaceHeader.H"

/*--------------------------------------------------------------------*/
//  Constructor
/** \param[in] a_file_name
 *                      The name of the ogen *.hdf file to read, default empty
 *//*-----------------------------------------------------------------*/
readOgen::readOgen(std::string a_file_name)
  :
  m_file_name(a_file_name),
  m_fileOpen(false),
  m_gridIdx(0)
{
  m_cg = std::unique_ptr<CompositeGrid>(new CompositeGrid());
  if(!m_file_name.empty())
    {
      openFile(m_file_name);
    }
}

/*--------------------------------------------------------------------*/
//  Destructor
/*--------------------------------------------------------------------*/
readOgen::~readOgen()
{
  closeFile();
}

/*--------------------------------------------------------------------*/
// read the composite grid
/*--------------------------------------------------------------------*/
void readOgen::openFile(std::string a_file_name)
{
  CH_TIME("readOgen::openFile");
  // // Initialize Overture
  // int ac = 0;
  // char *av[0] = NULL;
  // Overture::start(ac, av);

  // set the file name
  m_file_name = a_file_name;
    
  // Open the Ogen file
  aString ov_str(m_file_name);
  getFromADataBase(*m_cg, ov_str);
  m_cg->update();

  // check Overture and Chombo dimension match
  CH_assert(SpaceDim == m_cg->numberOfDimensions());
  
  // file open flag
  m_fileOpen = true;
}

/*--------------------------------------------------------------------*/
// close the file
/*--------------------------------------------------------------------*/
void readOgen::closeFile()
{
  CH_TIME("readOgen::closeFile");
  if(m_fileOpen)
    {
      // Overture::finish();
      m_fileOpen = false;
    }
}

/*--------------------------------------------------------------------*/
// gets the dimensions of a given Ogen file
// zone must already be selected with readZone
/** \param[in] 
 *//*-----------------------------------------------------------------*/
void readOgen::getDim(std::vector<int>& a_num_cells,
                      std::vector<Real>& a_dom_len,
                      std::vector<Real>& a_dom_orgn) const
{
  // size the vectors
  a_num_cells.resize(SpaceDim);
  a_dom_len.resize(SpaceDim);
  a_dom_orgn.resize(SpaceDim);
  // get size from Overture
  MappedGrid& mg = m_cg->operator[](m_gridIdx);
  mg.changeToAllCellCentered();
  intArray size = mg.gridIndexRange();
  
  // number of cells
  D_TERM(a_num_cells[0] = size(1, 0) - size(0, 0);,
         a_num_cells[1] = size(1, 1) - size(0, 1);,
         a_num_cells[2] = size(1, 2) - size(0, 2););
  
  // computational domain length (the number of cells)
  D_TERM(a_dom_len[0] = (Real)a_num_cells[0];,
         a_dom_len[1] = (Real)a_num_cells[1];,
         a_dom_len[2] = (Real)a_num_cells[2];);

  // // origin (point at xi=0)
  // int ovDim = mg.numberOfDimensions();
  // realArray r_lo(1, ovDim);
  // realArray x_lo(1, ovDim);
  // r_lo = 0;
  // mg.mapping().map(r_lo, x_lo, Overture::nullRealDistributedArray());
  // D_TERM(a_dom_orgn[0] = x_lo(0, 0);,
  //        a_dom_orgn[1] = x_lo(0, 1);,
  //        a_dom_orgn[2] = x_lo(0, 2););
  // origin is zero, this provides no offset in physical space
  D_TERM(a_dom_orgn[0] = 0.0;,
         a_dom_orgn[1] = 0.0;,
         a_dom_orgn[2] = 0.0;);
}

/*--------------------------------------------------------------------*/
/// given coordinate in mapped space, return its location in real space
//  Xi must be given in unit domain!
/*--------------------------------------------------------------------*/
RealVect readOgen::realCoord(const RealVect& a_Xi) const
{
  // set up Overture types
  MappedGrid& mg = m_cg->operator[](m_gridIdx);
  int ovDim = mg.numberOfDimensions();
  realArray r(1, ovDim), x(1, ovDim);
  // populate Overture types from Chombo input
  D_TERM(r(0,0) = a_Xi[0];,
         r(0,1) = a_Xi[1];,
         r(0,2) = a_Xi[2];);
  // get the mapping, this should be const but isnt for some reason
  mg.mapping().map(r, x, Overture::nullRealDistributedArray());
  // convert back from Overture to Chombo
  RealVect loc;
  D_TERM(loc[0] = x(0,0);,
         loc[1] = x(0,1);,
         loc[2] = x(0,2););
  return loc;
}

/*--------------------------------------------------------------------*/
/// given coordinate in real space, return its location in the mapped space
/*--------------------------------------------------------------------*/
RealVect readOgen::mappedCoord(const RealVect& a_x) const
{
  // set up Overture types
  MappedGrid& mg = m_cg->operator[](m_gridIdx);
  int ovDim = mg.numberOfDimensions();
  realArray r(1, ovDim), x(1, ovDim);
  // populate Overture types from Chombo input
  D_TERM(x(0,0) = a_x[0];,
         x(0,1) = a_x[1];,
         x(0,2) = a_x[2];);
  // get the inverse mapping, this should be const but isnt for some reason
  mg.mapping().inverseMap(x, r, Overture::nullRealDistributedArray());
  // convert back from Overture to Chombo
  RealVect xi;
  D_TERM(xi[0] = x(0,0);,
         xi[1] = x(0,1);,
         xi[2] = x(0,2););
  return xi;
}

/*--------------------------------------------------------------------*/
// Calculate the derivative of each coordinate vector
//  Xi must be given in unit domain!
/*--------------------------------------------------------------------*/
Real readOgen::dXdXi(const RealVect& a_Xi, int a_dirX, int a_dirXi) const
{
  // set up Overture types
  MappedGrid& mg = m_cg->operator[](m_gridIdx);
  int ovDim = mg.numberOfDimensions();
  realArray r(1, ovDim), xr(1, ovDim, ovDim);
  // populate Overture types from Chombo input
  D_TERM(r(0,0) = a_Xi[0];,
         r(0,1) = a_Xi[1];,
         r(0,2) = a_Xi[2];);
  // get the mapping derivatives, this should be const but isnt for some reason
  mg.mapping().map(r, Overture::nullRealDistributedArray(), xr);
  // convert back from Overture to Chombo
  Real dXdXi = xr(0, a_dirX, a_dirXi);
  return dXdXi;
}

#include "NamespaceFooter.H"
#endif
