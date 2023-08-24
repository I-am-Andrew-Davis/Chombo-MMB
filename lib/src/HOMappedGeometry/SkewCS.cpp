#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "SkewCS.H"
#include "NamespaceHeader.H"

/*--------------------------------------------------------------------*/
///  Constructor
/**
 *   \param[in]  a_dX   The computational cell size
 *   \param[in]  a_angle
 *                      The angle to skew the grid by
 *   \param[in]  a_physLength
 *                      The length and height of the parallelogram
 *   \param[in]  a_domainLen
 *                      Length of computational domain
 *//*-----------------------------------------------------------------*/
SkewCS::SkewCS(const RealVect& a_dX,
               const Real&     a_angle,
               const RealVect& a_physLength,
               const RealVect& a_domainLength)
  :
  m_angle(a_angle*4.0*std::atan(1.0)/180.),
  m_physLeng(a_physLength),
  m_domainLeng(a_domainLength)
{
  m_dx = a_dX;
}

/*--------------------------------------------------------------------*/
/// Destructor
/**
 *//*-----------------------------------------------------------------*/
SkewCS::~SkewCS()
{
}

/*--------------------------------------------------------------------*/
/// Given coordinate in mapped space, return its location in real space
/**
 *   \param[in]  a_Xi   Location in computational space
 *   \return            Location in physical space
 *//*-----------------------------------------------------------------*/
RealVect
SkewCS::realCoord(const RealVect& a_Xi) const
{
  RealVect realLoc;
  // adjust the scale
  realLoc = a_Xi*m_physLeng/m_domainLeng;
  // account for the angle
  D_TERM(,
         realLoc[0] += realLoc[1]/std::tan(m_angle);,
         realLoc[0] += realLoc[2]/std::tan(m_angle);
         realLoc[1] += realLoc[2]/std::tan(m_angle);)
  return realLoc;
}

/*--------------------------------------------------------------------*/
/// Given coordinate in real space, return its location in mapped space
/**
 *   \param[in]  a_x    Location in physical space
 *   \return            Location in mapped space
 *//*-----------------------------------------------------------------*/
RealVect
SkewCS::mappedCoord(const RealVect& a_x) const
{
  RealVect mapLoc;
  // adjust the scale
  mapLoc = a_x*m_domainLeng/m_physLeng;
  // account for the angle
  D_TERM(,
         mapLoc[0] -= a_x[1]*m_domainLeng[1]/(m_physLeng[1]*std::tan(m_angle));,
         mapLoc[0] -= a_x[2]*m_domainLeng[2]/(m_physLeng[2]*std::tan(m_angle));
         mapLoc[1] -= a_x[2]*m_domainLeng[2]/(m_physLeng[2]*std::tan(m_angle)););
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
SkewCS::dXdXi(const RealVect& a_Xi, int a_dirX, int a_dirXi) const
{
  Real value;
  if (a_dirX == a_dirXi)
    {
      value = m_physLeng[a_dirXi]/m_domainLeng[a_dirXi];
    }
  else if (a_dirXi > a_dirX)
    {
      value = m_physLeng[a_dirXi]/(m_domainLeng[a_dirXi]*std::tan(m_angle));
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
 *   \param[in]  a_angle
 *                      The angle to skew the grid by
 *   \param[in]  a_physLength
 *                      The length and height of the parallelogram
 *//*-----------------------------------------------------------------*/
SkewCSFactory::SkewCSFactory(const Real&     a_angle,
                             const RealVect& a_physicalLength)
{
  m_angle = a_angle;
  m_physLeng = a_physicalLength;
  if ((a_angle <=0) || (a_angle>=90))
    {
      MayDay::Error(
        "Bad value in SkewCS, angle be between 0 and 90");
    }
}

/*--------------------------------------------------------------------*/
///  Return the SkewCS object
/**
 *   \param[in]  a_domain
 *                      Problem domain the coordinate system is in
 *   \return            Pointer to the coordinate system
 *//*-----------------------------------------------------------------*/
NewCoordSys*
SkewCSFactory::getCoordSys(const ProblemDomain& a_domain,
                           const RealVect&      a_dx) const
{
  IntVect domSize = a_domain.size();
  RealVect domL = domSize*a_dx;
  SkewCS* newCSPtr = new SkewCS(a_dx, m_angle, m_physLeng, domL);
  return static_cast< NewCoordSys* >(newCSPtr);
}

#include "NamespaceFooter.H"


