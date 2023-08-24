#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "LogCS.H"
#include "NamespaceHeader.H"

/*--------------------------------------------------------------------*/
///  Constructor
/**
 *   \param[in]  a_dX   The computational cell size
 *   \param[in]  a_shift
 *                      The stretching factor to use for each direction.
 *                      Positive is compressed toward the high side, and
 *                      negative towards the low side. Stretching of 0 
 *                      results in a Cartesian grid
 *   \param[in]  a_domL
 *                      Length of computational domain
 *//*-----------------------------------------------------------------*/
LogCS::LogCS(const RealVect& a_dX,
             const RealVect& a_shift,
             const RealVect  a_domL)
  :
  m_shift(a_shift),
  m_eShift(a_shift),
  m_diffShift(a_shift),
  m_domL(a_domL)
{
  m_dx = a_dX;
  // scale determines the shift, positive near top and negative is bottom
  for(int d = 0; d != SpaceDim; ++d)
    {
      m_eShift[d] = exp(m_shift[d]) - 1;
      m_diffShift[d] = m_eShift[d]/m_shift[d];
    }
}

/*--------------------------------------------------------------------*/
/// Destructor
/**
 *//*-----------------------------------------------------------------*/
LogCS::~LogCS()
{
}

/*--------------------------------------------------------------------*/
/// Given coordinate in mapped space, return its location in real space
/**
 *   \param[in]  a_Xi   Location in computational space
 *   \return            Location in physical space
 *//*-----------------------------------------------------------------*/
RealVect
LogCS::realCoord(const RealVect& a_Xi) const
{
  RealVect realLoc;
  for (int d = 0; d != SpaceDim; ++d)
    {
      if (m_shift[d] == 0)
        {
          realLoc[d] = a_Xi[d];
        }
      else
        {
          realLoc[d] = log(a_Xi[d]/m_domL[d]*m_eShift[d]+1.)
            *m_domL[d]/m_shift[d];
        }
    }
  return realLoc;
}

/*--------------------------------------------------------------------*/
/// Given coordinate in real space, return its location in mapped space
/**
 *   \param[in]  a_x    Location in physical space
 *   \return            Location in mapped space
 *//*-----------------------------------------------------------------*/
RealVect
LogCS::mappedCoord(const RealVect& a_x) const
{
  RealVect mapLoc;
  for (int d = 0; d != SpaceDim; ++d)
    {
      if (m_shift[d] == 0.0)
        {
          mapLoc[d] = a_x[d];
        }
      else
        {
          mapLoc[d] = m_domL[d]*(exp(a_x[d]*m_shift[d]/m_domL[d])-1.)/
            m_eShift[d];
        }
    }
  return mapLoc;
}

/*--------------------------------------------------------------------*/
/// Calculate the derivative of each coordinate vector
/**
 *   \param[in]  a_Xi   Location in computational space
 *   \param[in]  a_dirX The component in physical space (x, y)
 *   \param[in]  a_dirXi
 *                      The component in computational space (xi, eta)  
 *   \return            derivative of dX/dXi
 *//*-----------------------------------------------------------------*/
Real
LogCS::dXdXi(const RealVect& a_Xi, int a_dirX, int a_dirXi) const
{
  Real value;
  if (a_dirX == a_dirXi)
    {
      if (m_shift[a_dirX] == 0)
        {
          value = 1.0;
        }
      else
        {
          value = m_diffShift[a_dirX]/
            (m_eShift[a_dirX]*a_Xi[a_dirX]/m_domL[a_dirX] + 1.);
        }
    }
  else
    {
      value = 0.0;
    }
  return value;
}


// -- begin factory implementations ---------------------------

/*--------------------------------------------------------------------*/
///  Factory Constructor
/**
 *   \param[in]  a_shift
 *                      The stretching factor to use for each direction.
 *                      Positive is compressed toward the high side, and
 *                      negative towards the low side. Stretching of 0 
 *                      results in a Cartesian grid
 *//*-----------------------------------------------------------------*/
LogCSFactory::LogCSFactory(const RealVect& a_shift)
{
  m_shift = a_shift;
  for(int d = 0; d != SpaceDim; ++d)
    {
      if(abs(m_shift[d]) > 3.5)
        {
          // large values cause problems
          MayDay::Error(
            "Bad value in LogCoordSys::shift, must be less than 3.5");
        }
    }
}

/*--------------------------------------------------------------------*/
///  Return the LogCS object
/**
 *   \param[in]  a_domain
 *                      Problem domain the coordinate system is in
 *   \return            Pointer to the coordinate system
 *//*-----------------------------------------------------------------*/
NewCoordSys*
LogCSFactory::getCoordSys(const ProblemDomain& a_domain,
                          const RealVect&      a_dx) const
{
  IntVect domSize = a_domain.size();
  RealVect domL;
  for(int dir = 0; dir != SpaceDim; ++dir)
    {
      domL[dir] = domSize[dir]*a_dx[dir];
    }
  LogCS* newCSPtr = new LogCS(a_dx, m_shift, domL);
  return static_cast< NewCoordSys* >(newCSPtr);
}

#include "NamespaceFooter.H"
