#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "AnnulusCS.H"
#include "NamespaceHeader.H"

/*--------------------------------------------------------------------*/
///  Constructor
/**
 *   \param[in]  a_dX   The computational cell size
 *   \param[in]  a_inRadius
 *                      Radius of the inner cylinder
 *   \param[in]  a_outRadius
 *                      Radius of the outer cylinder
 *   \param[in]  a_domainLen
 *                      Length of computational domain
 *   \param[in]  a_Cs   Radial stretching factor
 *//*-----------------------------------------------------------------*/
AnnulusCS::AnnulusCS(const RealVect& a_dX,
                     const Real&     a_inRadius,
                     const Real&     a_outRadius,
                     const RealVect& a_domainLen,
                     const Real&     a_Cs)
  : m_inRadius(a_inRadius),
    m_outRadius(a_outRadius),
    m_domainLen(a_domainLen),
    m_Cs(a_Cs),
    m_twoPi(8.0*atan(1.0))
{
  CH_assert(m_Cs != 0.0);
  m_dx = a_dX;
}

/*--------------------------------------------------------------------*/
/// Destructor
/**
 *//*-----------------------------------------------------------------*/
AnnulusCS::~AnnulusCS()
{
}

/*--------------------------------------------------------------------*/
/// Given coordinate in mapped space, return its location in real space
/**
 *   \param[in]  a_Xi   Location in computational space
 *   \return            Location in physical space
 *//*-----------------------------------------------------------------*/
RealVect
AnnulusCS::realCoord(const RealVect& a_Xi) const
{
  // only x, y mapped. z is un changed
  RealVect x_loc(a_Xi);
  D_TERM(
    Real r = m_inRadius + (m_outRadius - m_inRadius)
    *(exp(m_Cs*a_Xi[0]/m_domainLen[0]) - 1.)/(exp(m_Cs)-1.); ,
    Real theta = m_twoPi*(a_Xi[1]/m_domainLen[1]); ,
    Real z = a_Xi[2]; )
#if CH_SPACEDIM > 1
  D_TERM(x_loc[0] = r*cos(theta); ,
         x_loc[1] = r*sin(theta); ,
         x_loc[2] = z; )
#else
    //**FIXME Workaround for a compiler error in 1-D, has
    //not been verified that this definition is appropriate
    x_loc[0] = r;
#endif
  
  return x_loc;
}

/*--------------------------------------------------------------------*/
/// Given coordinate in real space, return its location in mapped space
/**
 *   \param[in]  a_x    Location in physical space
 *   \return            Location in mapped space
 *//*-----------------------------------------------------------------*/
RealVect
AnnulusCS::mappedCoord(const RealVect& a_x) const
{
  RealVect mappedXi(a_x);
  D_TERM(
    mappedXi[0] = m_Cs*m_domainLen[0]*log((m_outRadius - m_Cs*m_inRadius)/
                                          (m_outRadius - m_inRadius));,
    mappedXi[1] = atan(a_x[1]/a_x[0])*m_domainLen[1]/m_twoPi;,
    mappedXi[2] = a_x[2]; )
  return mappedXi;
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
AnnulusCS::dXdXi(const RealVect& a_Xi, int a_dirX, int a_dirXi) const
{
  Real r, theta, f;
  D_TERM(r = m_inRadius + (m_outRadius - m_inRadius)
         *(exp(m_Cs*a_Xi[0]/m_domainLen[0]) - 1.)/(exp(m_Cs)-1.); ,
         theta = m_twoPi*(a_Xi[1]/m_domainLen[1]); ,
    )
    
  // Dx/Dxi or Dx/Dr
  if ((a_dirXi == 0) && (a_dirX == 0))
    {
      f = cos(theta);
      r = (m_outRadius - m_inRadius)*m_Cs*exp(m_Cs*a_Xi[0]/m_domainLen[0])/
        (m_domainLen[0]*(exp(m_Cs) - 1.));
    }
  // Dy/Dxi or Dy/Dr
  else if ((a_dirXi == 0) && (a_dirX == 1))
    {
      f = sin(theta);
      r = (m_outRadius - m_inRadius)*m_Cs*exp(m_Cs*a_Xi[0]/m_domainLen[0])/
        (m_domainLen[0]*(exp(m_Cs) - 1.));
    }
  // Dx/Deta or Dx/Dtheta
  else if ((a_dirXi == 1) && (a_dirX == 0))
    {
      // r remains unchanged
      f = -sin(theta)*m_twoPi/m_domainLen[1];
    }
  // Dy/Deta or Dy/Dtheta
  else if ((a_dirXi == 1) && (a_dirX == 1))
    {
      // r remains unchanged
      f = cos(theta)*m_twoPi/m_domainLen[1];
    }
  else if ((a_dirXi == 2) && (a_dirX == 2))
    {
      return 1.0;
    }
  else
    {
      return 0.0;
    }
  
  return r*f;
}


// -- begin factory implementations ---------------------------

/*--------------------------------------------------------------------*/
///  Factory Constructor
/**
 *   \param[in]  a_inRadius
 *                      Radius of the inner cylinder
 *   \param[in]  a_outRadius
 *                      Radius of the outer cylinder
 *   \param[in]  a_radialStretch
 *                      Radial stretching factor
 *//*-----------------------------------------------------------------*/
AnnulusCSFactory::AnnulusCSFactory(const Real& a_inRadius,
                                   const Real& a_outRadius,
                                   const Real& a_radialStretch)
  : m_inRadius(a_inRadius),
    m_outRadius(a_outRadius),
    m_radialStretch(a_radialStretch)
{
    if (a_radialStretch == 0)
    {
      MayDay::Error("Bad value in AnnulusCoordSys::radialStretch");
    }
}

/*--------------------------------------------------------------------*/
///  Return the AnnulusCS object
/**
 *   \param[in]  a_domain
 *                      Problem domain the coordinate system is in
 *   \return            Pointer to the coordinate system
 *//*-----------------------------------------------------------------*/
NewCoordSys*
AnnulusCSFactory::getCoordSys(const ProblemDomain& a_domain,
                              const RealVect& a_dx) const
{
  IntVect domSize = a_domain.size();
  RealVect domL;
  for (int dir = 0; dir != SpaceDim; ++dir)
    {
      domL[dir] = domSize[dir]*a_dx[dir];
    }
  AnnulusCS* newCSPtr = new AnnulusCS(a_dx, 
                                      m_inRadius,
                                      m_outRadius,
                                      domL,
                                      m_radialStretch);
  return static_cast< NewCoordSys* >(newCSPtr);
}

#include "NamespaceFooter.H"


