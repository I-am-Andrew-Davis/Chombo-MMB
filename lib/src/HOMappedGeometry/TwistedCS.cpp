#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "TwistedCS.H"

#include "NamespaceHeader.H"

// helper functions
inline Real ipow_new( const Real x, const int n )
{
   Real val( x );
   for (int i=1; i<n; i++) val *= x;
   return val;
}

inline Real f_new( const Real r, const Real a )
{
   double val = (1.0 + cos( a * r ));
   return 0.5 * ipow_new(val,3);
}

#if 1
Real sinc_new( const Real x )
{
   Real val = 1.0;
   if ( fabs(x) > 1.0e-16 )
      val = sin(x) / x;
   else
   {
      static const Real A[] =
      {
        -1./6.,1./20., -1./42., 1./72., -1./110.
      };
      Real xsq = ipow_new(x,2);
      for (int i=4; i>=0; i--) val *= (1 - A[i] * xsq);
   }
   return val;
}

inline Real dfdronr_new( const Real r, const Real a )
{
   Real val = (1.0 + cos( a * r ));
   return (-1.5 * sinc_new( a * r ) * ipow_new(a*val,2) );
}
#else
inline Real dfdronr_new( const Real r, const Real a )
{
   Real val = (1.0 + cos( a * r ));
   return (-1.5 * a * sin( a * r ) * ipow_new(val,2) / r );
}
#endif



/// Calhoun-Leveque "twisted" Coordinate mapping (constant Jacobian)


/*--------------------------------------------------------------------*/
///  Constructor
/**
 *   \param[in]  a_dX   The computational cell size
 *   \param[in]  a_R    Radius of twist
 *   \param[in]  a_twist
 *                      The amount to twist by
 *   \param[in]  a_domainLength
 *                      Length of computational domain
 *//*-----------------------------------------------------------------*/
TwistedCS::TwistedCS(const RealVect& a_dX,
                     const Real& a_R,
                     const Real& a_twist,
                     const RealVect& a_domainLength)
{
  m_dx = a_dX;
  m_Pi = 4.0*atan(1.0);
  m_R = a_R;
  m_scale = m_Pi / m_R;
  m_theta = a_twist;
  m_domainLength = a_domainLength;
}

/*--------------------------------------------------------------------*/
/// Destructor
/**
 *//*-----------------------------------------------------------------*/
TwistedCS::~TwistedCS()
{
}

/*--------------------------------------------------------------------*/
/// Given coordinate in mapped space, return its location in real space
/**
 *   \param[in]  a_Xi   Location in computational space
 *   \return            Location in physical space
 *//*-----------------------------------------------------------------*/
RealVect
TwistedCS::realCoord(const RealVect& a_Xi) const
{
   RealVect xi(a_Xi);
   // Periodic offset for a domain from 0 to 1
   RealVect periodicOffset;
   for (int idir = 0; idir < SpaceDim; idir++)
     {
       periodicOffset[idir] =
         (int)(a_Xi[idir] < 0.) - (int)(
           a_Xi[idir] > m_domainLength[idir]) - 0.5*m_domainLength[idir];
     }
   xi += periodicOffset;
   Real r = sqrt( ipow_new(xi[0],2) + ipow_new(xi[1],2) );
   Real alpha = (r<m_R) ? f_new( r, m_scale) : 0;
   Real beta = alpha * m_theta;
   Real cosb = cos( beta );
   Real sinb = sin( beta );
   RealVect xyzLoc;
   D_TERM6(xyzLoc[0] = cosb * xi[0] + sinb * xi[1];,
           xyzLoc[1] = cosb * xi[1] - sinb * xi[0];,
           xyzLoc[2] = xi[2];,
           xyzLoc[3] = xi[3];,
           xyzLoc[4] = xi[4];,
           xyzLoc[5] = xi[5];)
   xyzLoc -= periodicOffset;
   return xyzLoc;
}

/*--------------------------------------------------------------------*/
/// Given coordinate in real space, return its location in mapped space
/**
 *   \param[in]  a_x    Location in physical space
 *   \return            Location in mapped space
 *//*-----------------------------------------------------------------*/
RealVect
TwistedCS::mappedCoord(const RealVect& a_x) const
{
   RealVect x(a_x);
   x -= 0.5*m_domainLength;
   Real r = sqrt( ipow_new(x[0],2) + ipow_new(x[1],2) );
   Real alpha = (r<m_R) ? f_new( r, m_scale) : 0;
   Real beta = alpha * m_theta;
   Real cosb = cos( beta );
   Real sinb = sin( beta );
   RealVect mappedXi;
   D_TERM6(mappedXi[0] = cosb * x[0] - sinb * x[1];,
           mappedXi[1] = cosb * x[1] + sinb * x[0];,
           mappedXi[2] = x[2];,
           mappedXi[3] = x[3];,
           mappedXi[4] = x[4];,
           mappedXi[5] = x[5];)
   mappedXi += 0.5*m_domainLength;
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
TwistedCS::dXdXi(const RealVect& a_Xi, int a_dirX, int a_dirXi) const
{
  RealVect xi(a_Xi);
  // Periodic offset for a domain from 0 to 1
   RealVect periodicOffset;
   for (int idir = 0; idir < SpaceDim; idir++)
     {
       periodicOffset[idir] =
         (int)(a_Xi[idir] < 0.) - (int)(
           a_Xi[idir] > m_domainLength[idir]) - 0.5*m_domainLength[idir];
     }
  xi += periodicOffset;
  Real r = sqrt( ipow_new(xi[0],2) + ipow_new(xi[1],2) );
  Real alpha = (r<m_R) ? f_new( r, m_scale) : 0;
  Real beta = alpha * m_theta;
  Real cosb = cos( beta );
  Real sinb = sin( beta );
  Real factor = (r<m_R) ? m_theta * dfdronr_new( r, m_scale ) : 0;

  Real value = 0.0;
  if (a_dirX == 0)
    {
      Real y = cosb * xi[1] - sinb * xi[0];
      if (a_dirXi == 0)
        {
          value =  cosb + xi[0] * y * factor;
        }
      else if (a_dirXi == 1)
        {
          value =  sinb + xi[1] * y * factor;
        }
      // else value stays 0
    }
  else if (a_dirX == 1)
    {
      Real x = cosb * xi[0] + sinb * xi[1];
      if (a_dirXi == 0)
        {
          value = -sinb - xi[0] * x * factor;
        }
      else if (a_dirXi == 1)
        {
          value =  cosb - xi[1] * x * factor;
        }
      // else value stays 0
    }
  else if (a_dirX == 2)
    {
      if (a_dirXi == 2)
        {
          value = 1.0;
        }
    }
  else
    {
      MayDay::Error("Bad dirX in TwistedCoordSys::dXdXi");
    }

  return value;
}


// -- begin factory implementations ---------------------------

/*--------------------------------------------------------------------*/
///  Factory Constructor
/**
 *   \param[in]  a_R    Radius of twist
 *   \param[in]  a_twist
 *                      The amount to twist by
 *//*-----------------------------------------------------------------*/
TwistedCSFactory::TwistedCSFactory(const Real& a_R,
                                   const Real& a_twist)
{
  m_R = a_R;
  m_twist = a_twist;
}

/*--------------------------------------------------------------------*/
///  Return the TwistedCS object
/**
 *   \param[in]  a_domain
 *                      Problem domain the coordinate system is in
 *   \return            Pointer to the coordinate system
 *//*-----------------------------------------------------------------*/
NewCoordSys*
TwistedCSFactory::getCoordSys(const ProblemDomain& a_domain,
                              const RealVect& a_dx) const
{
  RealVect domainLength = a_domain.size() * a_dx;
  TwistedCS* newCSPtr = new TwistedCS(a_dx, m_R, m_twist, domainLength);
  return static_cast< NewCoordSys* >(newCSPtr);
}

#include "NamespaceFooter.H"
