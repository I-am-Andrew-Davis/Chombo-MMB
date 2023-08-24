#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _READOGEN_H_
#define _READOGEN_H_

//#ifdef CH_USE_OGEN
#include <string>
#include <vector>
#include <memory>
#include "RealVect.H"

// forward declare from Overture
class CompositeGrid;

#include "NamespaceHeader.H"

class readOgen
{
public:
  // constructor
  readOgen(std::string a_file_name = "");

  // destructor
  ~readOgen();

  // sets initial file indecies
  void openFile(std::string a_file_name);
  
  // sets initial file indecies
  void closeFile();

  // get the grid dimensions
  void getDim(std::vector<int>& a_num_cells,
              std::vector<Real>& a_dom_len,
              std::vector<Real>& a_dom_orgn) const;

  // get the periodic conditions
  void getPeriodic(std::vector<int>& a_periodic) const;
  
  /// given coordinate in mapped space, return its location in real space
  RealVect realCoord(const RealVect& a_Xi) const;

  /// given coordinate in real space, return its location in the mapped space
  RealVect mappedCoord(const RealVect& a_x) const;

  // Calculate the derivative of each coordinate vector
  Real dXdXi(const RealVect& a_Xi, int a_dirX, int a_dirXi) const;

protected:
  std::string m_file_name;
  std::unique_ptr<CompositeGrid> m_cg;
  bool m_fileOpen;
  int m_gridIdx;
};


#include "NamespaceFooter.H"
//#endif
#endif
