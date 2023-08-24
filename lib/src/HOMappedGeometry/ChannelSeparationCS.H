#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// ChannelSeparationCS.H
// Andrew Davis, 15 Oct 2021

#ifndef _CHANNELSEPARATIONCS_H_
#define _CHANNELSEPARATIONCS_H_

#include "ExternalMultiBlockCS.H"
#include "ExtrudedExternalCS.H"
#include "NamespaceHeader.H"

//--Forward declarations
class ReadCGNS;

/// MultiBlockCoordSys consisting of two blocks in each dimension.
class ChannelSeparationCS : public ExternalMultiBlockCS
{
public:

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/

  ///
  /** Default constructor.
   */
 ChannelSeparationCS();

    ///
    /** Destructor.
     */
 virtual ~ChannelSeparationCS();

 virtual void define(const ProblemDomain& a_levelDomain,
                     const RealVect& a_dx,
                     const Real a_z_dx,
                     const Real a_l_dx,
                     const Real a_r_dx,
                     const int a_z_length,
                     const int a_l_length,
                     const int a_r_length,
                     const int a_blocksAreConformal,
                     const int a_numRampBlocks
 );

protected:

  int  NUMBLOCKS;
  int  m_blocksAreConformal;          ///< T -- block mappings are conformal
  int  m_numZCells;                   ///< Number of cells in z direction
  Real m_zDx;                         ///< Cell length in z direction

 // define m_boundaries
 void defineBoundaries();

private:

 void operator = (const ChannelSeparationCS&)
 {
 }

 ChannelSeparationCS(const ChannelSeparationCS&)
    {
    }
};


class ChannelSeparationCSFactory : public ExternalMultiBlockCSFactory
{

public:

 ChannelSeparationCSFactory(const std::string a_fileName,
                            const RealVect a_scale,
                            const int a_verbosity=0,
                            const Real a_z_dx=1.,
                            const Real a_l_dx=1.,
                            const Real a_r_dx=1.,
                            const int a_z_length=1,
                            const int a_l_length=1,
                            const int a_r_length=1,
                            const int a_blocksAreConformal=0,
                            const int a_numRampBlocks=3)
 :
      ExternalMultiBlockCSFactory(a_fileName,
                                  a_scale,
                                  a_verbosity),
      // these _dx are used to create the
      // spacing vector in the left, right, and z direction
      m_z_dx(a_z_dx),
      m_l_dx(a_l_dx),
      m_r_dx(a_r_dx),
      // these _lengths are used to create the
      // total number of cells in the left, right, and z directions
      m_z_length(a_z_length),
      m_l_length(a_l_length),
      m_r_length(a_r_length),
      m_blocksAreConformal(a_blocksAreConformal),
      m_numRampBlocks(a_numRampBlocks)

 {}
 
 virtual ~ChannelSeparationCSFactory()
 {}

 virtual ExternalMultiBlockCS* getCoordSys(const ProblemDomain& a_levelDomain,
                                           const RealVect& a_dx) const;

private:
 // these _dx are used to create the
 // spacing vector in the left, right, and z direction
 Real m_z_dx;
 Real m_l_dx;
 Real m_r_dx;
 // these _lengths are used to create the
 // total number of cells in the left, right, and z directions
 int m_z_length;
 int m_l_length;
 int m_r_length;
 int m_blocksAreConformal;
 int m_numRampBlocks;
};

#include "NamespaceFooter.H"

#endif
