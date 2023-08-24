#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "WarpedCS.H"
#include "NamespaceHeader.H"

/*--------------------------------------------------------------------*/
///  Constructor
/**
 *   \param[in]  a_dX   The computational cell size
 *   \param[in]  a_scale
 *                      scale for the amount of warping displacement
 *   \param[in]  a_domainLength
 *                      Length of computational domain
 *   \param[in]  a_rtol Numerical tolerance
 *   \param[in]  a_atol Numerical tolerance
 *   \param[in]  a_imax The maximum number of iterations
 *//*-----------------------------------------------------------------*/

WarpedCS::WarpedCS(const RealVect& a_dX,
                   const RealVect& a_scale,
                   const RealVect a_domainLength,
                   const Real     a_rtol,
                   const Real     a_atol,
                   const int      a_imax)
  :
  m_staticScale(a_scale*a_domainLength),
  m_twoPi(8.0*atan(1.0)),
  m_domainLength(a_domainLength)
{
  m_dx = a_dX;
  if ( a_rtol < 0 )
    {
      MayDay::Error("Bad value in WarpedCoordSys::relative_tolerance");
    }

  m_rtol = a_rtol;
  if ( a_atol < 0 )
    {
      MayDay::Error("Bad value in WarpedCoordSys::absolute_tolerance");
    }
  m_atol = a_atol;

  if ( a_imax < 0 )
    {
      MayDay::Error("Bad value in WarpedCoordSys::maximum_iterations");
    }
  m_imax = a_imax;
  RealVect scaleR;
  for(int dir = 0; dir != SpaceDim; ++dir)
    {
      scaleR[dir] = a_scale[dir]*m_domainLength[dir];
    }
  scale(scaleR);
}

/*--------------------------------------------------------------------*/
/// Destructor
/**
 *//*-----------------------------------------------------------------*/
WarpedCS::~WarpedCS()
{
}

/*--------------------------------------------------------------------*/
/// Given coordinate in mapped space, return its location in real space
/**
 *   \param[in]  a_Xi   Location in computational space
 *   \return            Location in physical space
 *//*-----------------------------------------------------------------*/
RealVect
WarpedCS::realCoord(const RealVect& a_Xi) const
{
  Real sine_product = 1.0;
  for (int d=0; d<SpaceDim; d++)
    sine_product *= sin( m_twoPi/m_domainLength[d] * a_Xi[d] );
  RealVect xyzLoc;
  for (int d=0; d<SpaceDim; d++)
    {
      xyzLoc[d] = a_Xi[d] + m_scale[d] * sine_product;
    }
  return xyzLoc;
}

/*--------------------------------------------------------------------*/
/// Given coordinate in real space, return its location in mapped space
/**
 *   \param[in]  a_x    Location in physical space
 *   \return            Location in mapped space
 *//*-----------------------------------------------------------------*/
RealVect
WarpedCS::mappedCoord(const RealVect& a_x) const
{
   RealVect mappedXi;
   for (int d=0; d<m_rootDir; d++)
      mappedXi[d] = a_x[d];

   if (m_rootDir<SpaceDim)
   {

      // find root by fixed point iteration
      Real root = a_x[m_rootDir];
      Real residual = g( root, a_x ) - root;
      Real bound = m_rtol * fabs(residual) + m_atol;
      int count = 0;
      while ( fabs(residual) > bound && count < m_imax )
      {
         root += residual;
         residual = g( root, a_x ) - root;
         count++;
      }
      if (count==m_imax)
      {
         MayDay::Error("Convergence failure in WarpedCoordSys::mappedCoord iteration!");
      }

      mappedXi[m_rootDir] = root;
      Real delta = (mappedXi[m_rootDir] - a_x[m_rootDir]) / m_scale[m_rootDir];
      for (int d=m_rootDir+1; d<SpaceDim; d++)
         mappedXi[d] = a_x[d] + m_scale[d] * delta;
   }
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
WarpedCS::dXdXi(const RealVect& a_Xi, int a_dirX, int a_dirXi) const
{
  Real product = 1.0;
  for (int d=0; d<SpaceDim; d++)
    {
      Real arg = m_twoPi/m_domainLength[d] * a_Xi[d];
      product *= ( d == a_dirXi ) ? cos( arg ) : sin( arg );
    }

  Real value = m_twoPi/m_domainLength[a_dirXi] * m_scale[a_dirX] * product;
   if (a_dirX == a_dirXi)
     {
       value += 1.0;
     }
   return value;
}

/*--------------------------------------------------------------------*/
//  Sets the time for moving grids and adjusts m_scale
/** m_scale is varied by cos(2 Pi a_time).  For directions < 3, the
 *  each cosine is offset by Pi/4 to ensure m_scale is never 0 for
 *  all directions.  But this only affects the inverse mapping which
 *  does not seem to be used.
 *//*-----------------------------------------------------------------*/

void WarpedCS::setTime(const Real& a_time) const
{
  // FIXME : This function has not been tested with non-unit based domains
  // Set base
  NewFourthOrderCoordSys::setTime(a_time);
  // Rescale
  m_scale = RealVect::Zero;
  m_rootDir = SpaceDim;
  Real maxScale = 0.;
  for (int d = 0; d != std::min(3, SpaceDim); ++d)
    {
      // Each direction < 3 is offset by Pi/4 to ensure a component of scale is
      // always > 0
      m_scale[d] = m_staticScale[d]*std::cos(m_twoPi*a_time/m_domainLength[d] + 
                                             0.125*d*m_twoPi/m_domainLength[d]);
      if (std::fabs(m_scale[d]) > maxScale)
        {
          m_rootDir = d;
          maxScale = std::fabs(m_scale[d]);
        }
    }
  // To avoid problems with precision, require a rootDir (but note that this
  // only affects the inverse mapping)
  CH_assert(m_rootDir != SpaceDim);
}

void
WarpedCS::scale(const RealVect& a_scale)
{
  Real value = 1.0 / m_twoPi;
  RealVect bound = value * m_domainLength;
  if ( a_scale > bound || a_scale < RealVect::Zero )
    {
      MayDay::Error("Bad scale in WarpedCoordSys::scale");
    }
  m_scale = a_scale;

  m_rootDir = SpaceDim;
  for (int d=0; d<SpaceDim; d++)
    {
      if (m_scale[d]>0)
        {
          m_rootDir = d;
          break;
        }
    }
}

// RHS of fixed-point iteration for inverse mapping
Real
WarpedCS::g( const Real xi, const RealVect& x ) const
{
  const Real delta = (xi - x[m_rootDir]) / m_scale[m_rootDir];
   Real sine_product = 1.0;
   for (int d=0; d<SpaceDim; d++)
     {
       Real eta = x[d];
       if (d!=m_rootDir)
         eta += m_scale[d] * delta;
       else
         eta = xi;
       sine_product *= sin( m_twoPi/m_domainLength[d] * eta );
     }
   return x[m_rootDir] - m_scale[m_rootDir] * sine_product;
}


// -- begin factory implementations ---------------------------

/*--------------------------------------------------------------------*/
///  Factory Constructor
/**
 *   \param[in]  a_scale
 *                      scale for the amount of warping displacement
 *   \param[in]  a_rtol Numerical tolerance
 *   \param[in]  a_atol Numerical tolerance
 *   \param[in]  a_imax The maximum number of iterations
 *//*-----------------------------------------------------------------*/
WarpedCSFactory::WarpedCSFactory(const RealVect& a_scale,
                                 const Real a_rtol,
                                 const Real a_atol,
                                 const int a_imax)
  :
  m_scale(a_scale),
  m_rtol(a_rtol),
  m_atol(a_atol),
  m_imax(a_imax)
{
}

/*--------------------------------------------------------------------*/
///  Return the WarpedCS object
/**
 *   \param[in]  a_domain
 *                      Problem domain the coordinate system is in
 *   \return            Pointer to the coordinate system
 *//*-----------------------------------------------------------------*/
NewCoordSys*
WarpedCSFactory::getCoordSys(const ProblemDomain& a_domain,
                             const RealVect&      a_dx) const
{
  IntVect domSize = a_domain.size();
  RealVect domainLength;
  for(int dir = 0; dir != SpaceDim; ++dir)
    {
      domainLength[dir] = domSize[dir]*a_dx[dir];
    }
  WarpedCS* newCSPtr = new WarpedCS(a_dx, 
                                    m_scale, 
                                    domainLength, 
                                    m_rtol, 
                                    m_atol, 
                                    m_imax);
  return static_cast< NewCoordSys* >(newCSPtr);
}

#include "NamespaceFooter.H"


