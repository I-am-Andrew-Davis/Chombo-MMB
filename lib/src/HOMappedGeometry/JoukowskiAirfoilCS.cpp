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

//----- Standard Library -----//

#include <cmath>
#include <complex>

//----- Chombo -----//

#include "JoukowskiAirfoilCS.H"
#include "parstream.H"
#include "StcMatrix.H"
#include "Misc.H"

#include "NamespaceHeader.H"


/*******************************************************************************
 *
 * Class JoukowskiAirfoilCS: member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Constructor
/** \param[in]  a_dx    Computational grid spacing
 *  \param[in]  a_chord Airfoil chord
 *  \param[in]  a_span  Airfoil span (if 3-D)
 *  \param[in]  a_domainRatio
 *                      There is a mapping from the zeta plane, where
 *                      both the airfoil and domain exterior are
 *                      circles.  This is a ratio of the outer domain
 *                      circle to the inner airfoil circle.
 *  \param[in]  a_camberRatio
 *                      Ratio of maximum camber to chord (approximate)
 *  \param[in]  a_thicknessRatio
 *                      Ratio of maximum thickness to chord
 *                      (approximate)
 *  \param[in]  a_alpha Angle of attack (in degrees!)
 *  \param[in]  a_cellRatioRadial
 *                      Cells are stretched using a tanh function in
 *                      the radial direction to have thin layers in
 *                      the boundary layer.  This is the ratio of the
 *                      smallest to largest cell.  Set to 1 for no
 *                      stretching.
 *  \param[in]  a_cellRatioAzimuth
 *                      Cells are stretched using a cosine function
 *                      in the azimuthal direction.  This can be used
 *                      to focus cells at the leading edge to resolve
 *                      a bow shock.  This is the ratio of the
 *                      smallest to largest cell.  Set to 1 for no
 *                      stretching.
 *//*-----------------------------------------------------------------*/

JoukowskiAirfoilCS::JoukowskiAirfoilCS(
  const RealVect& a_dx,
  const Real      a_chord,
  const Real      a_span,
  const Real      a_domainRatio,
  const Real      a_camberRatio,
  const Real      a_thicknessRatio,
  const Real      a_alpha,
  const Real      a_cellRatioRadial,
  const Real      a_cellRatioAzimuth)
  :
  m_chord(a_chord),
  m_span(a_span),
  m_camberRatio(a_camberRatio),
  m_thicknessRatio(a_thicknessRatio),
  m_alpha(a_alpha*Pi/180.),
  m_c(0.25*a_chord),
  m_kappa(2*a_camberRatio),
  m_epsilon(4*a_thicknessRatio/(3*sqrt(3.))),
  m_r0(m_c*sqrt(std::pow(1. + m_epsilon, 2) + std::pow(m_kappa, 2))),
  m_r1(a_domainRatio*m_r0),
  m_mu(-m_epsilon*m_c, m_kappa*m_c),
  m_stretchCosX(a_cellRatioAzimuth),
  m_stretchTanhY(a_cellRatioRadial, 0.7, 0.2, 0.8)
{

//--Base class

  m_dx = a_dx;

//--This class

  RealVect Xi(0.);
  m_XTE = realCoord(Xi);
  m_mappedGuessXi.reserve(13);
  m_mappedGuessX.reserve(13);
  Xi = { 0.01, 0.01, 0.5 };
  m_mappedGuessXi.emplace_back(Xi);
  m_mappedGuessX.emplace_back(realCoord(Xi));
  Xi = { 0.99, 0.01, 0.5 };
  m_mappedGuessXi.emplace_back(Xi);
  m_mappedGuessX.emplace_back(realCoord(Xi));
  Xi = { 0.01, 0.2, 0.5 };
  m_mappedGuessXi.emplace_back(Xi);
  m_mappedGuessX.emplace_back(realCoord(Xi));
  Xi = { 0.99, 0.2, 0.5 };
  m_mappedGuessXi.emplace_back(Xi);
  m_mappedGuessX.emplace_back(realCoord(Xi));
  for (int i = 0; i != 9; ++i)
    {
      RealVect Xi(0., 0.1, 0.5);
      Xi[0] = i*0.1 + 0.1;
      m_mappedGuessXi.emplace_back(Xi);
      m_mappedGuessX.emplace_back(realCoord(Xi));
    }
}

/*--------------------------------------------------------------------*/
//  Destructor
/*--------------------------------------------------------------------*/

JoukowskiAirfoilCS::~JoukowskiAirfoilCS()
{
}

/*--------------------------------------------------------------------*/
//  Given coordinate in mapped space, return its location in real
//  space
/** \param[in]  a_Xi    Location in computational space
 *  \return             Location in physical space
 *//*-----------------------------------------------------------------*/

RealVect
JoukowskiAirfoilCS::realCoord(const RealVect& a_Xi) const
{
  // Stretch
  const Real xi0Hat = m_stretchCosX.x(a_Xi[0]);
  const Real xi1Hat = m_stretchTanhY.x(a_Xi[1]);
  // To polar
  const Real r = xi1Hat*(m_r1 - m_r0) + m_r0;
  const Real th = (1. - xi0Hat)*(2*Pi);
  // To complex
  const std::complex<Real> zetaHat = std::polar(r, th);
  // Camber and thickness
  const std::complex<Real> zeta = zetaHat + m_mu;
  // Map to airfoil
  std::complex<Real> z = zeta + (m_c*m_c)/zeta;
  // Add AoA
  z *= std::polar<Real>(1., -m_alpha);
  return RealVect(D_DECL(z.real(), z.imag(), a_Xi[2]*m_span));
}

/*--------------------------------------------------------------------*/
//  Given coordinate in real space, return its location in the mapped
//  space
/** \param[in]  a_X     Location in physical space
 *  \return             Location in computational space
 *//*-----------------------------------------------------------------*/

RealVect
JoukowskiAirfoilCS::mappedCoord(const RealVect& a_X) const
{
  CH_assert(false);  // Needs cminpack to work.  Also would very much like
                     // to depreciate this routine.

  // Check if at singularity at TE
  constexpr int prec = 14;
  if (!Misc::compare1(a_X[0], m_XTE[0], prec) &&
      !Misc::compare1(a_X[1], m_XTE[1], prec))
    {
      return RealVect(D_DECL(0., 0., a_X[2]/m_span));
    }

  // Otherwise find a guess.
  int idxGuess[2] = { -1, -1 };
  Real len0 = std::numeric_limits<Real>::max();
  Real len1 = std::numeric_limits<Real>::max();
  for (int i = 0, i_end = m_mappedGuessX.size(); i != i_end; ++i)
    {
      auto diff = a_X - m_mappedGuessX[i];
      Real lenNew = stc::dot(diff, diff);
      if (lenNew < len0)
        {
          idxGuess[1] = idxGuess[0];
          idxGuess[0] = i;
          len1 = len0;
          len0 = lenNew;
        }
      else if (lenNew < len1)
        {
          idxGuess[1] = i;
          len1 = lenNew;
        }
    }

  m_X = a_X;
  RealVect fvec;
  stc::RVec<SpaceDim*SpaceDim> fjac;
  // Generally, this needs to be the same precision used for
  // compare1(Xi, mappedCoord(realCoord(Xi)))
  const Real xtol = 1.E-12;//sqrt(__cminpack_func__(dpmpar)(1));
  const int maxfev = SpaceDim*100;
  RealVect diag;
  const int mode = 1;
  const Real factor = 0.1;
  const int nprint = 0;
  int nfev;
  int njev;
  constexpr int lr = (SpaceDim*(SpaceDim + 1))/2;
  stc::RVec<lr> r;
  RealVect qtf;
  RealVect wa1, wa2, wa3, wa4;

  RealVect Xi;
  for (int iTrial = 0; iTrial != 2; ++iTrial)
    {
      Xi = m_mappedGuessXi[idxGuess[iTrial]];
// Needs cminpack as a dependency
#if 0
      const int info = __cminpack_func__(hybrj)(
        fcn, this, SpaceDim, Xi.data(), fvec.data(), fjac.data(), SpaceDim,
        xtol, maxfev, diag.data(),  mode, factor, nprint, &nfev, &njev,
        r.data(), lr, qtf.data(),
        wa1.data(), wa2.data(), wa3.data(), wa4.data());
#else
      (void)xtol;
      (void)maxfev;
      (void)mode;
      (void)nprint;
      int info = 1;
#endif
      // The actual Xi[1]
      constexpr Real Xi1min = -0.001;
      constexpr Real Xi1max = 1.001;
      constexpr Real Xi1mid = 0.5*(Xi1min + Xi1max);
      constexpr Real Xi1wth = 0.5*(Xi1max - Xi1min);
      const Real Xi1tanh = std::tanh((Xi[1] - Xi1mid)/Xi1wth);
      Xi[1] = Xi1mid + Xi1tanh*Xi1wth;
      if (info != 1 && info != 3)
        {
          if (iTrial == 1)
            {
              pout() << "Mapped coord failed multidim root solve\n";
              pout() << "  Guess : " << m_mappedGuessXi[idxGuess[iTrial]]
                     << std::endl;
              pout() << "  Xi    : " << Xi << std::endl;
              pout() << "  Factor: " << factor << std::endl;
              pout() << "  Params: " << nfev << ' ' << njev << ' ' << info
                     << std::endl;
            }
        }
      else
        {
          break;
        }
    }

  return Xi;
}

/*--------------------------------------------------------------------*/
//  Calculate the derivative of each coordinate vector in a Jacobian
/** \param[in]  a_Xi    Point in computational space
 *  \param[out] a_J     Jacobian matrix dXdXi with column major
 *                      ordering
 *  \param[in]  a_ldJ   The leading dimension of matrix a_J in memory
 *                      (>= a_Xi.size())
 *  The Jacobian matrix has column-major ordering with leading
 *  dimension a_ldJ
 *//*-----------------------------------------------------------------*/

void
JoukowskiAirfoilCS::dXdXiJacobian(const RealVect& a_Xi,
                                  Real*           a_J,
                                  int             a_ldJ) const
{
  // Stretch
  const Real xi0Hat = m_stretchCosX.x(a_Xi[0]);
  const Real xi1Hat = m_stretchTanhY.x(a_Xi[1]);
  // To polar
  const Real r = xi1Hat*(m_r1 - m_r0) + m_r0;
  const Real th = (1. - xi0Hat)*(2*Pi);
  // To complex
  const std::complex<Real> zetaHat = std::polar(r, th);
  // Camber and thickness
  const std::complex<Real> zeta = zetaHat + m_mu;

  // dz/d\zeta (embeds AoA)
  const std::complex<Real> dZdZeta =
    std::polar<Real>(1., -m_alpha)*(1. - m_c*m_c/(zeta*zeta));
  const stc::Matrix<Real, 2> JdZdZeta
    {
      dZdZeta.real(), -dZdZeta.imag(),
      dZdZeta.imag(),  dZdZeta.real()
    };

  // d\zeta/d\hat{\zeta} is identity

  // Switch to polar coordinates
  // d\hat{\zeta}/d\hat{\xi}_p
  const Real cosTheta = std::cos(th);
  const Real sinTheta = std::sin(th);
  const stc::Matrix<Real, 2> JdZetaHatdXip
    {
      cosTheta, -r*sinTheta,
      sinTheta,  r*cosTheta
    };

  // Then to stretched Cartesian
  // d\hat{\xi}_p/d\hat{\xi}
  const stc::Matrix<Real, 2> JdXipdXiHat
    {
      0.,    m_r1 - m_r0,
      -2*Pi, 0.
    };

  // Stetching
  // d\hat{\xi}/d\xi
  const stc::Matrix<Real, 2> JdXiHatdXi
    {
      m_stretchCosX.dxdxi(a_Xi[0]),  0.,
      0.,                            m_stretchTanhY.dxdxi(a_Xi[1])
    };

  const stc::Matrix<Real, 2> J =
    JdZdZeta*JdZetaHatdXip*JdXipdXiHat*JdXiHatdXi;

  //  col       row    r  c
  a_J[0*a_ldJ + 0] = J(0, 0);
  a_J[1*a_ldJ + 0] = J(0, 1);
  a_J[0*a_ldJ + 1] = J(1, 0);
  a_J[1*a_ldJ + 1] = J(1, 1);
  
  if (SpaceDim > 2)
    {
      a_J[0*a_ldJ + 2] = 0.;
      a_J[1*a_ldJ + 2] = 0.;
      a_J[2*a_ldJ + 0] = 0.;
      a_J[2*a_ldJ + 1] = 0.;
      a_J[2*a_ldJ + 2] = m_span;
    }
}

/*--------------------------------------------------------------------*/
//  Calculate the derivative of each coordinate vector
/** 
 *//*-----------------------------------------------------------------*/

Real
JoukowskiAirfoilCS::dXdXi(const RealVect& a_Xi,
                          int             a_dirX,
                          int             a_dirXi) const
{
  CH_assert(a_dirX  >= 0 && a_dirX  < SpaceDim);
  CH_assert(a_dirXi >= 0 && a_dirXi < SpaceDim);

  switch (3*a_dirX + a_dirXi)
    {
    case 3*0 + 0:
    case 3*0 + 1:
    case 3*1 + 0:
    case 3*1 + 1:
    {
      stc::RVec<SpaceDim*SpaceDim> J;
      dXdXiJacobian(a_Xi, J.data(), SpaceDim);
      //       col                row
      return J[a_dirXi*SpaceDim + a_dirX];
    }  
    case 3*2 + 2:  // dX2/dXi2
      return m_span;
    }
  return 0.;
}


/*******************************************************************************
 *
 * Class JoukowskiAirfoilCSFactory: member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Constructor
/** \param[in]  a_chord Airfoil chord
 *  \param[in]  a_span  Airfoil span (if 3-D)
 *  \param[in]  a_domainRatio
 *                      There is a mapping from the zeta plane, where
 *                      both the airfoil and domain exterior are
 *                      circles.  This is a ratio of the outer domain
 *                      circle to the inner airfoil circle.
 *  \param[in]  a_camberRatio
 *                      Ratio of maximum camber to chord (approximate)
 *  \param[in]  a_thicknessRatio
 *                      Ratio of maximum thickness to chord
 *                      (approximate)
 *  \param[in]  a_alpha Angle of attack (in degrees!)
 *  \param[in]  a_cellRatioRadial
 *                      Cells are stretched using a tanh function in
 *                      the radial direction to have thin layers in
 *                      the boundary layer.  This is the ratio of the
 *                      smallest to largest cell.  Set to 1 for no
 *                      stretching.
 *  \param[in]  a_cellRatioAzimuth
 *                      Cells are stretched using a cosine function
 *                      in the azimuthal direction.  This can be used
 *                      to focus cells at the leading edge to resolve
 *                      a bow shock.  This is the ratio of the
 *                      smallest to largest cell.  Set to 1 for no
 *                      stretching.
 *//*-----------------------------------------------------------------*/

JoukowskiAirfoilCSFactory::JoukowskiAirfoilCSFactory(
  const Real a_chord,
  const Real a_span,
  const Real a_domainRatio,
  const Real a_camberRatio,
  const Real a_thicknessRatio,
  const Real a_alpha,
  const Real a_cellRatioRadial,
  const Real a_cellRatioAzimuth)
  :
  m_chord(a_chord),
  m_span(a_span),
  m_domainRatio(a_domainRatio),
  m_camberRatio(a_camberRatio),
  m_thicknessRatio(a_thicknessRatio),
  m_alpha(a_alpha),
  m_cellRatioRadial(a_cellRatioRadial),
  m_cellRatioAzimuth(a_cellRatioAzimuth)
{ }

/*--------------------------------------------------------------------*/
//  Return a new coordinate system
/** \param[in]  a_domain
 *                      Unused
 *  \param[in]  a_dx    Domain spacing (note that the domain is 0:1
 *                      for this problem)
 *//*-----------------------------------------------------------------*/

NewCoordSys*
JoukowskiAirfoilCSFactory::getCoordSys(const ProblemDomain& a_domain,
                                       const RealVect&      a_dx) const
{
  JoukowskiAirfoilCS* newCSPtr = new JoukowskiAirfoilCS(
    a_dx,
    m_chord,
    m_span,
    m_domainRatio,
    m_camberRatio,
    m_thicknessRatio,
    m_alpha,
    m_cellRatioRadial,
    m_cellRatioAzimuth);
    
  return static_cast<NewCoordSys*>(newCSPtr);
}

#include "NamespaceFooter.H"
