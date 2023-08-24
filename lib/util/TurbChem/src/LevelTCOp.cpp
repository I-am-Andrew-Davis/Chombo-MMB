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
 * \file LevelTCOp.cpp
 *
 * \brief Member functions for LevelTCOp
 *
 *//*+*************************************************************************/

//----- Standard Library -----//

//----- Chombo Library -----//

#include "ParmParse.H"
#include "AMR.H"
#include "AMRLevel.H"
#include "CONSTANTS.H"
#include "AMRIO.H"
#include "CH_HDF5.H"

//----- Internal -----//

#include "LevelTCOp.H"
#include "SpectralUtil.H"


/*******************************************************************************
 *
 * Class LevelTCOp: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Constructor
/** \param[in]  a_something 
 *                      Something
 *//*-----------------------------------------------------------------*/

LevelTCOp::LevelTCOp(const int a_level)
:
  m_level(a_level)
{
}


/*==============================================================================
 * Public member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//  Define the level operator (weak construction)
/** \param[in]  a_level The level index
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *//*-----------------------------------------------------------------*/

void
LevelTCOp::define(HDF5Handle&            a_handle,
                  const std::string      a_fileName,
                  const int              a_numGhostCells,
                  const int              a_numOutGhostCells,
                  const int              a_numOutComps,
                  const std::vector<int> a_periodicity,
                  Real&                  a_timeBaseLevel)
{
  // Define the patch operator
  m_patchOp.define();
  HDF5HeaderData header;
  header.readFromFile(a_handle);
  // Shift data around due to issues with dbl and ghost cells
  LevelData<FArrayBox> inData;
  int numComps = header.m_int["num_components"];
  int ierr = readLevel<FArrayBox>(a_handle,
                                  m_level,
                                  inData,
                                  m_dx,
                                  m_dt,
                                  m_time,
                                  m_domainBox,
                                  m_refRatio);
  DisjointBoxLayout inDataDBL = inData.getBoxes();
  // Problem is, inDataDBL may not match periodicity
  bool periodicity[SpaceDim];
  for (int i = 0; i != SpaceDim; ++i)
    {
      periodicity[i] = true;
      if (!a_periodicity[i])
        {
          periodicity[i] = false;
        }
    }
  m_probDomain = ProblemDomain(m_domainBox, periodicity);
  // Read the header for this level -- change a_handle group to level
  char levelStr[32];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;
  a_handle.setGroup(label);
  // Try to get the correct DisjointBoxLayout
  Vector<Box> readBoxes;
  const int gridStatus = read(a_handle, readBoxes);
  ierr += gridStatus;
  Vector<int> procIDs;
  LoadBalance(procIDs, readBoxes);
  m_dbl = DisjointBoxLayout(readBoxes, procIDs, m_probDomain);
  m_dbl.close();

  IntVect ghostVect(a_numGhostCells*IntVect::Unit);
  // use weak construction for member level data
  m_dataOld.define(m_dbl, numComps, ghostVect);
  // use weak construction for output level data
  IntVect outputGhostVect(a_numOutGhostCells*IntVect::Unit);
  m_dataNew.define(m_dbl, a_numOutComps, outputGhostVect);
  // copy inData to m_dataOld
  inData.copyTo(m_dataOld);
  // Interval for exchanging ghost cells
  Interval compExchangeInterval(0,numComps-1);
  m_dataOld.exchange(compExchangeInterval);
  // Set time for CSV printing
  a_timeBaseLevel = m_time;
}

/*--------------------------------------------------------------------*/
//  Define the level operator (weak construction)
/** \param[in]  a_box   The disjoint box over which to post process
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *//*-----------------------------------------------------------------*/

void
LevelTCOp::postProcessCompute(const std::vector<std::string> a_outTimePntVars,
                              const std::vector<std::string> a_outTimeSumVars,
                              const std::vector<std::string> a_outGlobalVars,
                              std::vector<Real>&             a_domSums,
                              const int                      a_deconvOrder,
                              const int                      a_compOrder,
                              const int                      a_numReqGhost,
                              const int                      a_numComps,
                              const int                      a_numOutComps,
                              EnergySpectralDensity&         a_ESD)
{
  // FIXME: this whole function needs an overhaul
  // Define the energy spectral density class (note: these functions internally
  // check whether they should actually do something)
  a_ESD.defineSerial(m_probDomain);

  //int derivStencil = (a_compOrder + 1)/2;
  int convStencil = (a_deconvOrder - 1)/2;

  for (DataIterator dit = m_dbl.dataIterator(); dit.ok(); ++dit)
    {
      const Box& disjointBox = m_dbl[dit];
      // Create a box that covers all of the available ghost cells
      Box ghostBox   = grow(disjointBox, a_numReqGhost);
      // Create a box that represents the conservative point values
      Box consPntBox = grow(ghostBox, -convStencil);
      // Create a box that represents the primitive point values
      Box primPntBox = consPntBox;
      // Create a box that represents the primitive average values
      Box primAvgBox = grow(primPntBox, -convStencil);
      // Create a box that represents the post-processed point values
      Box zPntBox = grow(disjointBox, convStencil);
      // Make sure that the boxes only cover the physical domain
      ghostBox   &= m_probDomain;
      consPntBox &= m_probDomain;
      primPntBox &= m_probDomain;
      primAvgBox &= m_probDomain;
      zPntBox    &= m_probDomain;

      // Deconvolve the conservative variables as necessary
      FABSTACKTEMP(UpntFab, consPntBox, a_numComps);
      for (int i = 0; i != a_numComps; ++i)
        {
          m_patchOp.deconvolve4thOrd(UpntFab, m_dataOld[dit], i, consPntBox);
        }

      // Convert the conservative variables to primitive variables
      FABSTACKTEMP(WpntFab, primPntBox, (a_numComps+1));
      m_patchOp.consToPrim(WpntFab, UpntFab, primPntBox);

      // Load primitive state as required into spectral utility
      a_ESD.load(disjointBox, WpntFab, Interval(1, SpaceDim));

      // Compute the needed quantities -- testing with just enstrophy and KE
      // FIXME: the number of components in ZpntFab must depend on input
      FABSTACKTEMP(ZpntFab, zPntBox, a_numOutComps);
      MD_ARRAY_RESTRICT(arrZ, ZpntFab);

      // FIXME: each of the three input strings could have completely different
      //        variables from one another resulting in a larger number of
      //        necessary components than what we are accounting for here
      int numCompVars = std::max(
        a_outTimePntVars.size(), std::max(
          a_outTimeSumVars.size(), a_outGlobalVars.size()));
      for (int comps = 0; comps != numCompVars; ++comps)
        {
          if (a_outGlobalVars[comps] == "enstrophy")
            {
              // -- computing gradient of velocity
              Box gradVelBox = zPntBox;
              FABSTACKTEMP(GradVel, gradVelBox, SpaceDim*SpaceDim);
              Interval velIntv(1,SpaceDim);
              FArrayBox Velocity(velIntv, WpntFab);
              m_patchOp.gradOp4thOrd(GradVel, Velocity, m_dx, gradVelBox);

              // -- computing curl of velocity
              FABSTACKTEMP(CurlVel, gradVelBox, SpaceDim);
              m_patchOp.curlOpFromGrad(CurlVel, GradVel, gradVelBox);

              // -- computing the enstrophy
              FABSTACKTEMP(Enstrophy, gradVelBox, 1);
              MD_ARRAY_RESTRICT(arrE, Enstrophy);
              MD_ARRAY_RESTRICT(arrCurl, CurlVel);
              MD_BOXLOOP(gradVelBox, i)
                {
                  Real xCurl = arrCurl[MD_IX(i, 0)];
                  Real yCurl = arrCurl[MD_IX(i, 1)];
                  Real zCurl = arrCurl[MD_IX(i, 2)];
                  arrE[MD_IX(i, 0)]
                    = (1./2.)*(xCurl*xCurl + yCurl*yCurl + zCurl*zCurl);
                  // -- transfer enstrophy to output variables
                  arrZ[MD_IX(i, comps)] = arrE[MD_IX(i, 0)];
                }
              for (int sum = 0; sum != a_outTimeSumVars.size(); ++sum)
                {
                  if ((a_outTimeSumVars[sum] == "enstrophy") &&
                      (m_level == 0))
                    {
                      // sum up the enstrophy in the domain
                      Real domainSum = 0;
                      m_patchOp.domSumOp(domainSum,
                                         Enstrophy,
                                         disjointBox);
                      // weight it with the cell volume
                      a_domSums[sum] += domainSum*m_dx.product();
                    }
                }
            }
          else if (a_outGlobalVars[comps] == "KE")
            {
              // -- computing the kinetic energy
              FABSTACKTEMP(KE, zPntBox, 1);
              MD_ARRAY_RESTRICT(arrKE, KE);
              MD_ARRAY_RESTRICT(arrW, WpntFab);
              MD_BOXLOOP(zPntBox, i)
                {
                  Real rho  = arrW[MD_IX(i, 0)];
                  Real uVel = arrW[MD_IX(i, 1)];
                  Real vVel = arrW[MD_IX(i, 2)];
                  Real wVel = arrW[MD_IX(i, 3)];
                  arrKE[MD_IX(i, 0)]
                    = (1./2.)*rho*(uVel*uVel + vVel*vVel + wVel*wVel);
                }
              // -- transfer kinetic energy to output variables
              MD_BOXLOOP(zPntBox, i)
                {
                  arrZ[MD_IX(i, comps)] = arrKE[MD_IX(i, 0)];
                }
              for (int sum = 0; sum != a_outTimeSumVars.size(); ++sum)
                {
                  if ((a_outTimeSumVars[sum] == "KE") &&
                      (m_level == 0))
                    {
                      // sum up the kinetic energy in the domain
                      Real domainSum = 0;
                      m_patchOp.domSumOp(domainSum,
                                         KE,
                                         disjointBox);
                      // weight it with the cell volume
                      a_domSums[sum] += domainSum*m_dx.product();
                    }
                }
            }
          else if (a_outGlobalVars[comps] == "vorticity")
            {
              // compute vorticity here
            }
          else if (a_outGlobalVars[comps] == "palinstrophy")
            {
              // compute palinstrophy here
            }
          else if (a_outGlobalVars[comps] == "helicity")
            {
              // compute helicity here
            }
          else if (a_outGlobalVars[comps] == "entropy")
            {
              // compute entropy here
            }
          else if (a_outGlobalVars[comps] == "rms_vel")
            {
              // compute rms velocity here
            }
          else if (a_outGlobalVars[comps] == "grad_rms_vel")
            {
              // compute the gradient of rms velocity here
            }
        }

      // Average the needed quantities
      // FIXME: the number of components in ZavgFab must depend on ZpntFab
      FABSTACKTEMP(ZavgFab, disjointBox, a_numOutComps);
      for (int i = 0; i != a_numOutComps; ++i)
        {
          m_patchOp.convolve4thOrd(ZavgFab, ZpntFab, i, disjointBox);
        }

      // Copy the data for printing to HDF5 file
      m_dataNew[dit].copy(ZavgFab);
    }

  // Compute the ESD
  if (a_ESD.isUsed())
    {
      RealVect physDomain(m_probDomain.size());
      physDomain *= m_dx;
      a_ESD.compute(physDomain);
    }
}

/*--------------------------------------------------------------------*/
//  Define the level operator (weak construction)
/** \param[in]  a_level The level index
 *  \param[in]  a_dx    Vector of mesh spacing in computation space
 *//*-----------------------------------------------------------------*/

int
LevelTCOp::writeDataHDF5(HDF5Handle&                    a_handle,
                         const std::vector<std::string> a_outGlobalVars,
                         const int                      a_numOutComps,
                         const int                      a_numOutGhost)
{
  
  // FIXME: probably don't need a second header
  HDF5HeaderData header2;
  // FIXME: refRatio[0] only gives x-component refinement ratio
  //        change to allow anisotropic refinement
  header2.m_int["ref_ratio"]      = m_refRatio[0];
  header2.m_real["time"]          = m_time;
  header2.m_int["num_components"] = a_numOutComps;
  // FIXME: need general component name methodology
  for (int i = 0; i != a_numOutComps; ++i)
    {
      char compStr[30];
      sprintf(compStr,"component_%d", i);
      header2.m_string[compStr] = a_outGlobalVars[i];
    }

  header2.writeToFile(a_handle);

  const Interval numOutComps(0,(a_numOutComps-1));
  const IntVect numOutputGhost(a_numOutGhost*IntVect::Unit);

  // Write the data for this level
  int ierr = writeLevel<LevelData<FArrayBox>>(a_handle,
                                              m_level,
                                              m_dataNew,
                                              m_dx,
                                              m_dt,
                                              m_time,
                                              m_domainBox,
                                              m_refRatio,
                                              numOutputGhost,
                                              numOutComps);
  return ierr;
}
