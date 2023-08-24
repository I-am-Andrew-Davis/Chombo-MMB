#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <limits>

#include "LogSchwarzChristoffelRampCS.H"
#include "parstream.H"

#include "NamespaceHeader.H"


/*******************************************************************************
 *
 * Floating-point utilities
 *
 ******************************************************************************/

namespace FP
{
  const int defaultPrec = std::numeric_limits<Real>::digits10 - 1;
                                      ///< Default number of base 10 significant
                                      ///< digits to compare (1 less than
                                      ///< represented)
  const Real defaultTol = std::pow(10., -std::abs(defaultPrec));
                                      ///< Default tolerance

/*--------------------------------------------------------------------*/
/// Floating-point comparison with limit tol^2 as x and y -> 0.
/** \param[in] x        First real
 *  \param[in] y        Second real
 *  \param[in] prec     Number of base 10 significant digits to
 *                      compare
 *  \return             T - Not equal within tolerance
 *//*-----------------------------------------------------------------*/

  template <typename T>
  inline bool compare(const T &x, const T &y, int prec)
  {
    const T tol = std::pow(10., -std::abs(prec));
    return std::fabs(x - y) >
      (std::min(std::fabs(x), std::fabs(y)) + tol)*tol;
  }

  template <typename T>
  inline bool compare(const T &x, const T &y)
  {
    const T tol = (T)defaultTol;
    return std::fabs(x - y) >
      (std::min(std::fabs(x), std::fabs(y)) + tol)*tol;
  }

/*--------------------------------------------------------------------*/
/// Floating-point comparison with limit tol as x and y -> 0.
/** \param[in] x        First real
 *  \param[in] y        Second real
 *  \param[in] prec     Number of base 10 significant digits to
 *                      compare
 *  \return             T - Not equal within tolerance
 *//*-----------------------------------------------------------------*/

  template <typename T>
  inline bool compare1(const T &x, const T &y, int prec)
  {
    const T tol = std::pow(10., -std::abs(prec));
    return std::fabs(x - y) >
      std::min(std::fabs(x), std::fabs(y))*tol + tol;
  }

  template <typename T>
  inline bool compare1(const T &x, const T &y)
  {
    const T tol = (T)defaultTol;
    return std::fabs(x - y) >
      std::min(std::fabs(x), std::fabs(y))*tol + tol;
  }
}  // namespace FP


const double zeroTol = FP::defaultTol;


/*******************************************************************************
 *
 * Class LogSchwarzChristoffelRampCS: member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Constructor
/** \param[in]  a_dx    Computational grid spacing
 *  \param[in]  a_nCellBaseDir0
 *                      Number of x-Dir. cells in the domain for the
 *                      base mesh
 *  \param[in]  a_alpha The ramp angle  (0 <= alpha < Pi)
 *                      (default (Pi/6))
 *  \param[in]  a_XLeadLength
 *                      The length of the wall ahead of the ramp
 *                      measured in physical space.  This will likely
 *                      be increased slightly to ensure a grid line
 *                      starts at the corner (default 0.5)
 *  \param[in]  a_XRampLength
 *                      The length of the ramp measured in physical
 *                      space along the ramp (default 4.0)
 *//*-----------------------------------------------------------------*/

LogSchwarzChristoffelRampCS::LogSchwarzChristoffelRampCS(
  const RealVect& a_dx,
  const int       a_nCellBaseDir0,
  const Real      a_alpha,
  const Real      a_XLeadLength,
  const Real      a_XRampLength,
  const RealVect& a_shift,
  const RealVect& a_domL)
  :
  m_alpha(a_alpha),
  m_c(1. - a_alpha/Pi),
  m_XiCorner(D_DECL6(0., 0., 0., 0., 0., 0.)),
  m_shift(a_shift),
  m_domL(a_domL)
{
  CH_assert(a_nCellBaseDir0 > 2);
  CH_assert(a_alpha >= 0.);
  CH_assert(a_XLeadLength >= 0.);
  CH_assert(a_XRampLength > 0.);
  CH_assert(m_c > zeroTol);

//--Base class

  m_dx = a_dx;

//--This class

  // Initial estimate of lead distance in computation space
  const Real XiLeadFraction =
    1./(1. + std::pow(a_XRampLength/a_XLeadLength, 1./m_c));
  // But we have to adjust this to ensure a grid point lies at Xi = 0
  m_XiCorner[0] = -std::ceil(XiLeadFraction*a_nCellBaseDir0)/a_nCellBaseDir0;
  // Now recompute k (to keep the ramp length fixed)
  m_k = (a_XLeadLength + a_XRampLength)/
    (std::pow(std::fabs(m_XiCorner[0]), m_c) +
     std::pow(1. - std::fabs(m_XiCorner[0]), m_c));
  // scale determines the shift, positive near top and negative is bottom
  for (int d = 0; d != SpaceDim; ++d)
    {
      m_eShift[d] = exp(m_shift[d]) - 1;
      m_diffShift[d] = m_eShift[d]/m_shift[d];
    }
}

/*--------------------------------------------------------------------*/
//  Destructor
/*--------------------------------------------------------------------*/

LogSchwarzChristoffelRampCS::~LogSchwarzChristoffelRampCS()
{
}

/*--------------------------------------------------------------------*/
//  Given coordinate in mapped space, return its location in real
//  space
/** \note
 *  <ul>
 *    <li> Quadratures along edges (in 3D) for cells adjacent to the
 *         line (xi0 = 0, xi1 < 0) would not be quite right.  Theta
 *         on the line is adjusted to \f$\pi + \alpha)/(-2*c)\f$ and
 *         this adjustment should be applied linearly through those
 *         cells.  But, not a problem in 2D.
 *  </ul>
 *//*-----------------------------------------------------------------*/

RealVect
LogSchwarzChristoffelRampCS::realCoord(const RealVect& a_Xi) const
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
          realLoc[d] = std::log(a_Xi[d]/m_domL[d]*m_eShift[d] + 1.)
            *m_domL[d]/m_shift[d];
        }
    }
  // Translation so corner at Xi=(0,0,0)
  const RealVect Xi = realLoc + m_XiCorner;

  const RealVect XiPolar = getPolar(Xi);
  const Real radc = m_k*std::pow(XiPolar[0], m_c);

  const Real beta = m_c*XiPolar[1] + m_alpha;
  RealVect z;
  D_TERM6(z[0] = radc*(cos(beta));,
          z[1] = radc*(sin(beta));,
          z[2] = m_k*Xi[2];,
          z[3] = 0.;,
          z[4] = 0.;,
          z[5] = 0.;)

// Debug -- test reversal
#if 0
  RealVect XiRev = mappedCoord(z);
  if (D_TERM(FP::compare1(a_Xi[0], XiRev[0]),
          || FP::compare1(a_Xi[1], XiRev[1]),
          || FP::compare1(a_Xi[2], XiRev[2])))
    {
      pout() << "Who: "
        D_TERM(<< FP::compare1(a_Xi[0], XiRev[0]),
               << ' ' << FP::compare1(a_Xi[1], XiRev[1]),
               << ' ' << FP::compare1(a_Xi[2], XiRev[2]))
             << std::endl;
      pout() << "X   : " << z << std::endl;
      pout() << "Orig: " << a_Xi << std::endl;
      pout() << "Calc: " << XiRev << std::endl;
      MayDay::Abort("Coordinate mappings not within tolerances");
    }
#endif

  return z;
}

/*--------------------------------------------------------------------*/
// Given coordinate in real space, return its location in the mapped
// space
/** \note
 *  <ul>
 *    <li> Quadratures along edges (in 3D) for cells adjacent to the
 *         line (xi0 = 0, xi1 < 0) would not be quite right.  Theta
 *         on the line is adjusted to \f$\pi + \alpha)/(-2*c)\f$ and
 *         this adjustment should be applied linearly through those
 *         cells.  But, not a problem in 2D.
 *  </ul>
 *//*-----------------------------------------------------------------*/

RealVect
LogSchwarzChristoffelRampCS::mappedCoord(const RealVect& a_X) const
{
  RealVect mapLoc;
  for (int d = 0; d != SpaceDim; ++d)
    {
      if (m_shift[d] == 0.0)
        {
          mapLoc[d] = a_X[d];
        }
      else
        {
          mapLoc[d] = m_domL[d]*(std::exp(a_X[d]*m_shift[d]/m_domL[d])-1.)/
            m_eShift[d];
        }
    }
  // Radius in mapped space
  const Real rXi = std::pow(sqrt(std::pow(mapLoc[0], 2) +
                                 std::pow(mapLoc[1], 2))/m_k,
                            1./m_c);
  // Angle in physical space
  Real beta;
  if (std::fabs(mapLoc[0]) < zeroTol)
    {
      if (mapLoc[1] < 0.)
        {
          beta = -0.5*Pi;
        }
      else
        {
          beta = 0.5*Pi;
        }
    }
  else if (std::fabs(mapLoc[1]) < zeroTol)
    {
      if (mapLoc[0] < 0)
        {
          beta = Pi;
        }
      else
        {
          beta = 0.;
        }
    }
  else
    {
      beta = atan(mapLoc[1]/mapLoc[0]);
      if (mapLoc[0] < 0.)
        {
          beta += Pi;
        }
    }
  // Angle in mapped space
  Real thetaXi = beta;
  // Remap beta = (alpha - Pi)/2 to -Pi/2 to handle extension of corner
  if (!FP::compare(beta, (Real)(0.5*(m_alpha - Pi))))
    {
      thetaXi = -0.5*Pi;
    }
  // We scale to theta only if the radius is > 0.  Otherwise, the meaning of
  // scaling is undefined and thetaXi is already = Pi/2, the assumed angle at
  // this point.
  else if (rXi > zeroTol)
    {
      thetaXi = (thetaXi - m_alpha)/m_c;
    }
  const RealVect Xi(D_DECL6(rXi*cos(thetaXi), rXi*sin(thetaXi), mapLoc[2]/m_k,
                            0., 0., 0.));
  return (Xi - m_XiCorner);
}

/*--------------------------------------------------------------------*/
//  Calculate the derivative of each coordinate vector
/** \note
 *  <ul>
 *    <li> Fails when \f$r=0\f$ where the derivatives go to plus or
 *         minus infinity.  Derivatives set to 1.
 *    <li> Quadratures along edges (in 3D) for cells adjacent to the
 *         line (xi0 = 0, xi1 < 0) would not be quite right.  Theta
 *         on the line is adjusted to \f$\pi + \alpha)/(-2*c)\f$ and
 *         this adjustment should be applied linearly through those
 *         cells.  But, not a problem in 2D.
 *  </ul>
 *//*-----------------------------------------------------------------*/

Real LogSchwarzChristoffelRampCS::dXdXi(const RealVect& a_Xi,
                                     int             a_dirX,
                                     int             a_dirXi) const
{
  CH_assert(a_dirX  >= 0 && a_dirX  < SpaceDim);
  CH_assert(a_dirXi >= 0 && a_dirXi < SpaceDim);

  // Use chain rule to fine dX/dXi_L * dXi_L/dXi

  // Find d x/d xi_L, where xi_L is the log stretched mapping
  Real dXdXiL;
  if (a_dirX == a_dirXi)
    {
      if (m_shift[a_dirX] == 0)
        {
          dXdXiL = 1.0;
        }
      else
        {
          dXdXiL = m_diffShift[a_dirX]/
            (m_eShift[a_dirX]*a_Xi[a_dirX]/m_domL[a_dirX] + 1.);
        }
    }
  else
    {
      dXdXiL = 0.0;
    }

  // Transform Xi to Xi_L which is the location in the log stretched mapping
  RealVect XiL;
  for (int d = 0; d != SpaceDim; ++d)
    {
      if (m_shift[d] == 0)
        {
          XiL[d] = a_Xi[d];
        }
      else
        {
          XiL[d] = log(a_Xi[d]/m_domL[d]*m_eShift[d]+1.)
            *m_domL[d]/m_shift[d];
        }
    }
  const RealVect Xi = XiL + m_XiCorner;

  switch (3*a_dirX + a_dirXi)
    {
    case 3*0 + 0:  // dX0/dXi0
    case 3*1 + 1:  // dX1/dXi1
    {
      const RealVect XiPolar = getPolar(Xi);
      if (std::fabs(XiPolar[0]) < zeroTol)
        {
          return 1.*dXdXiL;
        }
      else
        {
          const Real beta = m_c*XiPolar[1] + m_alpha;
          return m_k*std::pow(XiPolar[0], m_c - 2.)*m_c*
            (Xi[0]*cos(beta) + Xi[1]*sin(beta))*dXdXiL;
        }
      break;
    }
    case 3*0 + 1:  // dX0/dXi1
    case 3*1 + 0:  // dX1/dXi0
    {
      const RealVect XiPolar = getPolar(Xi);
      if (std::fabs(XiPolar[0]) < zeroTol)
        {
          return -1.*dXdXiL;
        }
      else
        {
          const Real beta = m_c*XiPolar[1] + m_alpha;
          const Real retVal = m_k*std::pow(XiPolar[0], m_c - 2.)*m_c*
            (Xi[1]*cos(beta) - Xi[0]*sin(beta));
          return dXdXiL*((a_dirX == 0) ? retVal : -retVal);
        }
      break;
    }
    case 3*2 + 2:  // dX2/dXi2
      return m_k*dXdXiL;
      break;
    }
  return 0.;
}

/*--------------------------------------------------------------------*/
//  Get a complex in polar coordinates
/** \param[in]  a_Xi    Location in computational space.  A
 *                      translation should have already been applied
 *                      so that a_Xi = (0,0) is the corner
 *  \return             a_Xi in polar coordinates (a_Xi[2] unchanged
 *                      if present)
 *  \note
 *  <ul>
 *    <li> theta is set to Pi/2 at (0,0)
 *  </ul>
 *//*-----------------------------------------------------------------*/

RealVect
LogSchwarzChristoffelRampCS::getPolar(const RealVect& a_Xi) const
{
  // Theta
  Real theta;
  if (std::fabs(a_Xi[0]) < zeroTol)
    {
      if (a_Xi[1] < 0.)
        {
          // At the extension of the corner, there are two solutions,
          // approaching from -theta or +theta.  Instead, extend the line
          // intersecting the corner with angle -(Pi - alpha)/2.
          theta = (Pi + m_alpha)/(-2*m_c);
        }
      else
        {
          theta = 0.5*Pi;
        }
    }
  else if (std::fabs(a_Xi[1]) < zeroTol)
    {
      if (a_Xi[0] < 0.)
        {
          theta = Pi;
        }
      else
        {
          theta = 0.;
        }
    }
  else
    {
      theta = atan(a_Xi[1]/a_Xi[0]);
      if (a_Xi[0] < 0.)
        {
          theta = Pi + theta;
        }
    }

  // Radius
  const Real r = sqrt(std::pow(a_Xi[0], 2) + std::pow(a_Xi[1], 2));

  return RealVect(D_DECL6(r, theta, a_Xi[2], 0., 0., 0.));
}


/*******************************************************************************
 *
 * Class LogSchwarzChristoffelRampCSFactory: member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
// Constructor
/** \param[in]  a_nCellBaseDir0
 *                      Number of x-Dir. cells in the domain for the
 *                      base mesh
 *  \param[in]  a_alpha The ramp angle  (0 <= alpha < Pi)
 *                      (default (Pi/6))
 *  \param[in]  a_XLeadLength
 *                      The length of the wall ahead of the ramp
 *                      measured in physical space.  This will likely
 *                      be increased slightly to ensure a grid line
 *                      starts at the corner (default 0.5)
 *  \param[in]  a_XRampLength
 *                      The length of the ramp measured in physical
 *                      space along the ramp (default 4.0)
 *//*-----------------------------------------------------------------*/

LogSchwarzChristoffelRampCSFactory::LogSchwarzChristoffelRampCSFactory(
  const int       a_nCellBaseDir0,
  const Real      a_alpha,
  const Real      a_XLeadLength,
  const Real      a_XRampLength,
  const RealVect& a_shift)
  :
  m_nCellBaseDir0(a_nCellBaseDir0),
  m_alpha(a_alpha),
  m_XLeadLength(a_XLeadLength),
  m_XRampLength(a_XRampLength),
  m_shift(a_shift)
{
  for(int d = 0; d != SpaceDim; ++d)
    {
      if(abs(m_shift[d]) > 3.5)
        {
          // large values cause problems
          MayDay::Error(
            "Bad value in LogSchwarzChristoffelRampCSCoordSys::shift, must be less than 3.5");
        }
    }
}

/*--------------------------------------------------------------------*/
// Return a new coordinate system
/*--------------------------------------------------------------------*/

NewCoordSys*
LogSchwarzChristoffelRampCSFactory::getCoordSys(const ProblemDomain& a_domain,
                                                const RealVect&      a_dx) const
{
  IntVect domSize = a_domain.size();
  RealVect domL;
  for(int dir = 0; dir != SpaceDim; ++dir)
    {
      domL[dir] = domSize[dir]*a_dx[dir];
    }
  LogSchwarzChristoffelRampCS* newCSPtr = new LogSchwarzChristoffelRampCS(
    a_dx, m_nCellBaseDir0, m_alpha, m_XLeadLength, m_XRampLength, m_shift,
    domL);
  return static_cast<NewCoordSys*>(newCSPtr);
}

#include "NamespaceFooter.H"
