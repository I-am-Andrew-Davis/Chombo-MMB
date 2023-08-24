#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "CubedSphereShellPanel3DCS.H"
#include "BoxIterator.H"
#include <cmath>
#include "MBMiscUtil.H"
#include "csignal"
#include "RootSolver.H"

#include "CubedSphereShellF_F.H"
#include "FourthOrderUtil.H"
#include "FourthOrderUtilF_F.H"
#include "FCDivergenceF_F.H"
#include "AdvectOpF_F.H"

#include "NamespaceHeader.H"
// #include "UsingNamespace.H"

CubedSphereShellPanel3DCS::CubedSphereShellPanel3DCS(
  const RealVect&      a_dX,
  const ProblemDomain& a_domain,
  const int            a_blkIdx,
  const IntVect&       a_evalGhost,
  IntVect&             a_ix,
  const Real           a_radiusInner,
  const Real           a_radiusOuter,
  const Real           a_alpha,
  const Real           a_beta)
  :
  m_domain(a_domain),
  m_evalGhost(a_evalGhost),
  m_blkIdx(a_blkIdx),
  m_radiusInner(a_radiusInner),
  m_radiusOuter(a_radiusOuter),
  m_alpha(a_alpha),
  m_beta(a_beta),
  m_ix(a_ix)
{
  CH_assert((a_blkIdx >= 0) && (a_blkIdx < SpaceDim*2));

  m_dx = a_dX;
  // This is for AMR
  m_xiOrigin = RealVect{m_dx*m_domain.domainBox().smallEnd()};
}

CubedSphereShellPanel3DCS::~CubedSphereShellPanel3DCS()
{
}

RealVect
CubedSphereShellPanel3DCS::realCoord(const RealVect& a_Xi) const
{
  CH_TIME("CubedSphereShellPanel3DCS::realCoord");

  Real c_0 = 0.25*M_PI;

  // Set the output RealVect to nan -- this could help for debugging
  RealVect xVect(D_DECL(0./0., 0./0., 0./0.));

#if CH_SPACEDIM == 3
  RealVect halfBlockLength = 0.5*m_dx*(m_domain.domainBox().size());
  RealVect blockLength = m_dx*(m_domain.domainBox().size());
  Real xi_hat = c_0*(a_Xi[0] - m_xiOrigin[0] - halfBlockLength[0])/(
    halfBlockLength[0]);
  Real eta_hat = a_Xi[1]/blockLength[1];
  Real zeta_hat = c_0*(a_Xi[2] - m_xiOrigin[2] - halfBlockLength[2])/(
    halfBlockLength[2]);

  Real a_eta = (m_radiusOuter - m_radiusInner) - m_alpha*std::tanh(m_beta);
  Real f_eta = m_radiusInner + a_eta*eta_hat
    + m_alpha*std::tanh(m_beta*eta_hat);
  Real dGXi = std::tan(xi_hat);
  Real dGZeta = std::tan(zeta_hat);
  Real invBeta = 1./std::sqrt(1. + dGXi*dGXi + dGZeta*dGZeta);

  if (m_blkIdx == 0)
    {
      D_TERM(xVect[0] =  f_eta*invBeta;,
             xVect[1] = -f_eta*invBeta*dGZeta;,
             xVect[2] = -f_eta*invBeta*dGXi;);
    }
  else if (m_blkIdx == 1)
    {
      D_TERM(xVect[0] = -f_eta*invBeta*dGXi;,
             xVect[1] = -f_eta*invBeta*dGZeta;,
             xVect[2] = -f_eta*invBeta;);
    }
  else if (m_blkIdx == 2)
    {
      D_TERM(xVect[0] = -f_eta*invBeta;,
             xVect[1] = -f_eta*invBeta*dGZeta;,
             xVect[2] =  f_eta*invBeta*dGXi;);
    }
  else if (m_blkIdx == 3)
    {
      D_TERM(xVect[0] =  f_eta*invBeta*dGXi;,
             xVect[1] = -f_eta*invBeta*dGZeta;,
             xVect[2] =  f_eta*invBeta;);
    }
  else if (m_blkIdx == 4) // North polar block
    {
      D_TERM(xVect[0] = f_eta*invBeta*dGXi;,
             xVect[1] = f_eta*invBeta;,
             xVect[2] = f_eta*invBeta*dGZeta;);
    }
  else if (m_blkIdx == 5)
    {
      D_TERM(xVect[0] =  f_eta*invBeta*dGZeta;,
             xVect[1] = -f_eta*invBeta;,
             xVect[2] =  f_eta*invBeta*dGXi;);
    }
  else
    {
      MayDay::Error("CubedSphere block index must be 0 -- 5");
    }
#elif CH_SPACEDIM == 2
  RealVect halfBlockLength = 0.5*m_dx*(m_domain.domainBox().size());
  RealVect blockLength = m_dx*(m_domain.domainBox().size());
  Real xi_hat = c_0*(a_Xi[0] - m_xiOrigin[0] - halfBlockLength[0])/(
    halfBlockLength[0]);
  Real eta_hat = a_Xi[1]/blockLength[1];

  Real a_eta = (m_radiusOuter - m_radiusInner) - m_alpha*std::tanh(m_beta);
  Real f_eta = m_radiusInner + a_eta*eta_hat
    + m_alpha*std::tanh(m_beta*eta_hat);
  Real dGXi = std::tan(xi_hat);
  Real invBeta = 1./std::sqrt(1. + dGXi*dGXi);

  if (m_blkIdx == 0)
    {
      D_TERM(xVect[0] =  f_eta*invBeta;,
             xVect[1] = -f_eta*invBeta*dGXi;,);
    }
  else if (m_blkIdx == 1)
    {
      D_TERM(xVect[0] = -f_eta*invBeta*dGXi;,
             xVect[1] = -f_eta*invBeta;,);
    }
  else if (m_blkIdx == 2)
    {
      D_TERM(xVect[0] = -f_eta*invBeta;,
             xVect[1] =  f_eta*invBeta*dGXi;,);
    }
  else if (m_blkIdx == 3)
    {
      D_TERM(xVect[0] =  f_eta*invBeta*dGXi;,
             xVect[1] =  f_eta*invBeta;,);
    }
  else
    {
      MayDay::Error("CubedSphere block index must be 0 -- 3");
    }
#endif

  return xVect;
}

// given coordinate in real space, return its location in the mapped space
RealVect
CubedSphereShellPanel3DCS::mappedCoord(const RealVect& a_x) const
{
  CH_TIME("CubedSphereShellPanel3DCS::mappedCoord");
  // No matter the block orientation, the r component can always be
  // computed from x_0, x_1, and x_2
  const Real r = a_x.vectorLength();

  // Set the sphere-panel angle extent -- pi/4 from center either direction
  Real c_0 = 0.25*M_PI;

  RealVect xiVect = RealVect::Zero;

#if CH_SPACEDIM == 3
  Real xi = 0.;
  Real eta = 0.;
  Real zeta = 0.;
  switch (m_blkIdx)
  {
    case 0:
      xi   = std::atan2(-a_x[2], a_x[0]);
      zeta = std::atan2(-a_x[1], a_x[0]);
      break;

    case 1:
      xi   = std::atan2(-a_x[0], -a_x[2]);
      zeta = std::atan2(-a_x[1], -a_x[2]);
      break;

    case 2:
      xi   = std::atan2(a_x[2], -a_x[0]);
      zeta = std::atan2(-a_x[1], -a_x[0]);
      break;

    case 3:
      xi   = std::atan2(a_x[0], a_x[2]);
      zeta = std::atan2(-a_x[1], a_x[2]);
      break;

    case 4:
      xi   = std::atan2(a_x[0], a_x[1]);
      zeta = std::atan2(a_x[2], a_x[1]);
      break;

    case 5:
      xi   = std::atan2(a_x[2], -a_x[1]);
      zeta = std::atan2(a_x[0], -a_x[1]);
      break;
  }

  IntVect intSpaceBlockOrigin = m_domain.domainBox().smallEnd();
  RealVect realSpaceBlockOrigin = intSpaceBlockOrigin*m_dx;
  IntVect intSpaceBlockSize = m_domain.domainBox().size();
  RealVect realSpaceBlockSize = intSpaceBlockSize*m_dx;

  // First, compute the radial coordinate in the 0 - 1 space
  // Assume that we have perhaps 12 ghost cells (at most) in the radial dir
  int numRadGhost = 12;
  Real ghostLength = numRadGhost*m_dx[1];
  Real etaMin = 0. - ghostLength; // min eta_hat (without ghosts) is 0
  Real etaMax = 1. + ghostLength; // max eta_hat (without ghosts) is 1 
  int iterBrent  = 0;
  int errorBrent = 0;
  const RadiusFunc& f =
    RadiusFunc(r, m_radiusInner, m_radiusOuter, m_alpha, m_beta);
  Real etaSoln =
    RootSolver::BrentER(iterBrent, errorBrent, f, etaMin, etaMax);
  if (errorBrent != 0 || etaSoln != etaSoln)
    {
      MayDay::Error("CubedSphereShellPanel3DCS::mappedCoord bad etaSoln");
    }
  eta = etaSoln;
  // Next, scale it and shift it
  Real eta_scaled =
    eta*realSpaceBlockSize[1] + realSpaceBlockOrigin[1];

  // Scale xi
  Real xi_scaled = 0.5*realSpaceBlockSize[0]*(1./c_0)*xi
    + realSpaceBlockOrigin[0] + 0.5*realSpaceBlockSize[0];

  // Scale zeta
  Real zeta_scaled = 0.5*realSpaceBlockSize[2]*(1./c_0)*zeta
    + realSpaceBlockOrigin[2] + 0.5*realSpaceBlockSize[2];

  D_TERM(xiVect[0] = xi_scaled;,
         xiVect[1] = eta_scaled;,
         xiVect[2] = zeta_scaled;);
#elif CH_SPACEDIM == 2
  Real xi = 0.;
  Real eta = 0.;
  switch (m_blkIdx)
  {
    case 0:
      xi = std::atan2(-a_x[1], a_x[0]);
      break;
    case 1:
      xi = std::atan2(-a_x[0], -a_x[1]);
      break;
    case 2:
      xi = std::atan2(a_x[1], -a_x[0]);
      break;
    case 3:
      xi = std::atan2(a_x[0], a_x[21]);
      break;
  }

  IntVect intSpaceBlockOrigin = m_domain.domainBox().smallEnd();
  RealVect realSpaceBlockOrigin = intSpaceBlockOrigin*m_dx;
  IntVect intSpaceBlockSize = m_domain.domainBox().size();
  RealVect realSpaceBlockSize = intSpaceBlockSize*m_dx;

  // First, compute the radial coordinate in the 0 - 1 space
  // Assume that we have perhaps 12 ghost cells (at most) in the radial dir
  int numRadGhost = 12;
  Real ghostLength = numRadGhost*m_dx[1];
  Real etaMin = 0. - ghostLength; // min eta_hat (without ghosts) is 0
  Real etaMax = 1. + ghostLength; // max eta_hat (without ghosts) is 1 
  int iterBrent  = 0;
  int errorBrent = 0;
  const RadiusFunc& f =
    RadiusFunc(r, m_radiusInner, m_radiusOuter, m_alpha, m_beta);
  Real etaSoln =
    RootSolver::BrentER(iterBrent, errorBrent, f, etaMin, etaMax);
  if (errorBrent != 0 || etaSoln != etaSoln)
    {
      MayDay::Error("CubedSphereShellPanel3DCS::mappedCoord bad etaSoln");
    }
  eta = etaSoln;
  // Next, scale it and shift it
  Real eta_scaled =
    eta*realSpaceBlockSize[1] + realSpaceBlockOrigin[1];

  // Scale xi
  Real xi_scaled = 0.5*realSpaceBlockSize[0]*(1./c_0)*xi
    + realSpaceBlockOrigin[0] + 0.5*realSpaceBlockSize[0];

  D_TERM(xiVect[0] = xi_scaled;,
         xiVect[1] = eta_scaled;,);
#endif

  return xiVect;
}

Real
CubedSphereShellPanel3DCS::dXdXi(const RealVect& a_Xi,
                                 int             a_dirX,
                                 int             a_dirXi) const
{
  CH_TIME("CubedSphereShellPanel3DCS::dXdXi");
  CH_assert((m_blkIdx >= 0) && (m_blkIdx <= 5));

  Real value = 0.0;

  // First, we need to transform a_Xi into the shifted, scaled mapped coords
  IntVect blockSize = m_domain.domainBox().size();
  RealVect blockLength = blockSize*m_dx;
  RealVect halfBlockLength = 0.5*blockSize*m_dx;
  // Make it safe to divide by halfBlockTanLength (doesn't change eta)
  halfBlockLength[1] = 1;
  // Now shift a_Xi so (0, 0, 0) is at the center of the block
  RealVect shiftedXi = a_Xi - m_xiOrigin - halfBlockLength;
  // Scale xi so that it goes from -1 to +1
  RealVect scaledShiftedXi = (0.25*M_PI)*shiftedXi/halfBlockLength;
  // Set the sphere-normal coordinate back to what it was originally
  // except for the scaling
  scaledShiftedXi[1] = a_Xi[1]/blockLength[1];

#if CH_SPACEDIM == 3

  // Now, we can essentially get the gnomonic coordinates
  Real dGXi   = std::tan(scaledShiftedXi[0]);
  Real dGZeta = std::tan(scaledShiftedXi[2]);
  Real c_x = dGXi*dGXi;
  Real c_z = dGZeta*dGZeta;
  Real dInvDelta  = 1./std::sqrt(1. + c_x + c_z);
  Real dInvDelta3 = dInvDelta*dInvDelta*dInvDelta;

  // Set up some derivatives for using chain rule
  Real dXiScaled_dXi = (0.25*M_PI)/(halfBlockLength[0]);
  Real dEtaScaled_dEta = 1./blockLength[1];
  Real dZetaScaled_dZeta = (0.25*M_PI)/(halfBlockLength[2]);

  // Set up some other functions that can be precomputed independent of
  // which block is being considered
  Real a_eta = (m_radiusOuter - m_radiusInner) - m_alpha*std::tanh(m_beta);
  Real f_eta = m_radiusInner + a_eta*scaledShiftedXi[1]
    + m_alpha*std::tanh(m_beta*scaledShiftedXi[1]);
  // Note that the following is d(f(eta))/d(eta) not some really tasty food
  Real tanh_eta = std::tanh(m_beta*scaledShiftedXi[1]);
  Real dFeta_dEta = a_eta + m_alpha*m_beta*(1. - tanh_eta*tanh_eta);

  // Derivatives based on block 4
  Real dxdxi   = (1. + c_x)*(1. + c_z)*dInvDelta3*f_eta*dXiScaled_dXi;
  Real dxdeta  = dFeta_dEta*dInvDelta*dEtaScaled_dEta*dGXi;
  Real dxdzeta = -dGZeta*(1. + c_z)*dInvDelta3*f_eta*dGXi*dZetaScaled_dZeta;
  Real dydxi   = -dGXi*(1. + c_x)*dInvDelta3*f_eta*dXiScaled_dXi;
  Real dydeta  = dFeta_dEta*dInvDelta*dEtaScaled_dEta;
  Real dydzeta = -dGZeta*(1. + c_z)*dInvDelta3*f_eta*dZetaScaled_dZeta;
  Real dzdxi   = -dGXi*(1. + c_x)*dInvDelta3*f_eta*dGZeta*dXiScaled_dXi;
  Real dzdeta  = dFeta_dEta*dInvDelta*dEtaScaled_dEta*dGZeta;
  Real dzdzeta = (1. + c_x)*(1. + c_z)*dInvDelta3*f_eta*dZetaScaled_dZeta;

  if (a_dirX == 0)
    {
      if (a_dirXi == 0)
        {
          switch (m_blkIdx)
            {
            case 0:
              value = dydxi;
              break;

            case 1:
              value = -dxdxi;
              break;

            case 2:
              value = -dydxi;
              break;

            case 3:
              value = dxdxi;
              break;

            case 4:
              value = dxdxi;
              break;

            case 5:
              value = dzdxi;
              break;
            }
        }
      else if (a_dirXi == 1)
        {
          switch (m_blkIdx)
            {
            case 0:
              value = dydeta;
              break;

            case 1:
              value = -dxdeta;
              break;

            case 2:
              value = -dydeta;
              break;

            case 3:
              value = dxdeta;
              break;

            case 4:
              value = dxdeta;
              break;

            case 5:
              value = dzdeta;
              break;
            }
        }
      else if (a_dirXi == 2)
        {
          switch (m_blkIdx)
            {
            case 0:
              value = dydzeta;
              break;

            case 1:
              value = -dxdzeta;
              break;

            case 2:
              value = -dydzeta;
              break;

            case 3:
              value = dxdzeta;
              break;

            case 4:
              value = dxdzeta;
              break;

            case 5:
              value = dzdzeta;
              break;
            }
        }
    }
  else if (a_dirX == 1)
    {
      if (a_dirXi == 0)
        {
          switch (m_blkIdx)
            {
            case 0:
              value = -dzdxi;
              break;

            case 1:
              value = -dzdxi;
              break;

            case 2:
              value = -dzdxi;
              break;

            case 3:
              value = -dzdxi;
              break;

            case 4:
              value = dydxi;
              break;

            case 5:
              value = -dydxi;
              break;
            }
        }
      else if (a_dirXi == 1)
        {
          switch (m_blkIdx)
            {
            case 0:
              value = -dzdeta;
              break;

            case 1:
              value = -dzdeta;
              break;

            case 2:
              value = -dzdeta;
              break;

            case 3:
              value = -dzdeta;
              break;

            case 4:
              value = dydeta;
              break;

            case 5:
              value = -dydeta;
              break;
            }
        }
      else if (a_dirXi == 2)
        {
          switch (m_blkIdx)
            {
            case 0:
              value = -dzdzeta;
              break;

            case 1:
              value = -dzdzeta;
              break;

            case 2:
              value = -dzdzeta;
              break;

            case 3:
              value = -dzdzeta;
              break;

            case 4:
              value = dydzeta;
              break;

            case 5:
              value = -dydzeta;
              break;
            }
        }
    }
  else if (a_dirX == 2)
    {
      if (a_dirXi == 0)
        {
          switch (m_blkIdx)
            {
            case 0:
              value = -dxdxi;
              break;

            case 1:
              value = -dydxi;
              break;

            case 2:
              value = dxdxi;
              break;

            case 3:
              value = dydxi;
              break;

            case 4:
              value = dzdxi;
              break;

            case 5:
              value = dxdxi;
              break;
            }
        }
      else if (a_dirXi == 1)
        {
          switch (m_blkIdx)
            {
            case 0:
              value = -dxdeta;
              break;

            case 1:
              value = -dydeta;
              break;

            case 2:
              value = dxdeta;
              break;

            case 3:
              value = dydeta;
              break;

            case 4:
              value = dzdeta;
              break;

            case 5:
              value = dxdeta;
              break;
            }
        }
      else if (a_dirXi == 2)
        {
          switch (m_blkIdx)
            {
            case 0:
              value = -dxdzeta;
              break;

            case 1:
              value = -dydzeta;
              break;

            case 2:
              value = dxdzeta;
              break;

            case 3:
              value = dydzeta;
              break;

            case 4:
              value = dzdzeta;
              break;

            case 5:
              value = dxdzeta;
              break;
            }
        }
    }
#elif CH_SPACEDIM == 2
  // Now, we can essentially get the gnomonic coordinates
  Real dGXi   = std::tan(scaledShiftedXi[0]);
  Real c_x = dGXi*dGXi;
  Real dInvDelta  = 1./std::sqrt(1. + c_x);
  Real dInvDelta3 = dInvDelta*dInvDelta*dInvDelta;

  // Set up some derivatives for using chain rule
  Real dXiScaled_dXi = (0.25*M_PI)/(halfBlockLength[0]);
  Real dEtaScaled_dEta = 1./blockLength[1];

  // Set up some other functions that can be precomputed independent of
  // which block is being considered
  Real a_eta = (m_radiusOuter - m_radiusInner) - m_alpha*std::tanh(m_beta);
  Real f_eta = m_radiusInner + a_eta*scaledShiftedXi[1]
    + m_alpha*std::tanh(m_beta*scaledShiftedXi[1]);
  Real tanh_eta = std::tanh(m_beta*scaledShiftedXi[1]);
  Real dFeta_dEta = a_eta + m_alpha*m_beta*(1. - tanh_eta*tanh_eta);

  // Derivatives based on block 4
  Real dxdxi   = (1. + c_x)*dInvDelta3*f_eta*dXiScaled_dXi;
  Real dxdeta  = dFeta_dEta*dInvDelta*dEtaScaled_dEta*dGXi;
  Real dydxi   = -dGXi*(1. + c_x)*dInvDelta3*f_eta*dXiScaled_dXi;
  Real dydeta  = dFeta_dEta*dInvDelta*dEtaScaled_dEta;

  if (a_dirX == 0)
    {
      if (a_dirXi == 0)
        {
          switch (m_blkIdx)
            {
            case 0:
              value = dydxi;
              break;

            case 1:
              value = -dxdxi;
              break;

            case 2:
              value = -dydxi;
              break;

            case 3:
              value = dxdxi;
              break;
            }
        }
      else if (a_dirXi == 1)
        {
          switch (m_blkIdx)
            {
            case 0:
              value = dydeta;
              break;

            case 1:
              value = -dxdeta;
              break;

            case 2:
              value = -dydeta;
              break;

            case 3:
              value = dxdeta;
              break;
            }
        }
    }
  else if (a_dirX == 1)
    {
      if (a_dirXi == 0)
        {
          switch (m_blkIdx)
            {
            case 0:
              value = -dxdxi;
              break;

            case 1:
              value = -dydxi;
              break;

            case 2:
              value = dxdxi;
              break;

            case 3:
              value = dydxi;
              break;
            }
        }
      else if (a_dirXi == 1)
        {
          switch (m_blkIdx)
            {
            case 0:
              value = -dxdeta;
              break;

            case 1:
              value = -dydeta;
              break;

            case 2:
              value = dxdeta;
              break;

            case 3:
              value = dydeta;
              break;
            }
        }
    }
#endif

  return value;
}

#include "NamespaceFooter.H"
