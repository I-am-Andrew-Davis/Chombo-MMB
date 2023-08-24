#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <iostream>

#include "RigidTransformation.H"
#include "BoxIterator.H"
// #include "UGIO.H"
// #include "CH_HDF5.H"

#include "NamespaceHeader.H"

// ---------------------------------------------------------
RigidTransformation::RigidTransformation()
{
  defineAsIdentity();
}

// ---------------------------------------------------------
RigidTransformation::RigidTransformation(const RealVect& a_rotation,
                                         const RealVect& a_translation)
{
  define(a_rotation, a_translation);
}


// ---------------------------------------------------------
void
RigidTransformation::define(const RealVect& a_rotation,
                            const RealVect& a_translation)
{
  m_rotation = a_rotation;
  m_translation = a_translation;
}


// ---------------------------------------------------------
void
RigidTransformation::defineAsIdentity()
{
  m_rotation = RealVect::Zero;
  m_translation = RealVect::Zero;
}


// ---------------------------------------------------------
bool
RigidTransformation::operator==(const RigidTransformation& a_itOther) const
{
  return ((m_rotation == a_itOther.m_rotation) &&
          (m_translation == a_itOther.m_translation));
}


// ---------------------------------------------------------
RealVect
RigidTransformation::transformFwd(const RealVect& a_ivOld) const
{
  return transformVectorFwd(a_ivOld) + m_translation;
}


// ---------------------------------------------------------
RealVect
RigidTransformation::transformBack(const RealVect& a_ivNew) const
{
  return transformVectorBack(a_ivNew) - m_translation;
}


// ---------------------------------------------------------
RealVect
RigidTransformation::transform(const RealVect& a_iv,
                               bool a_forward) const
{
  if (a_forward)
    return transformFwd(a_iv);
  else
    return transformBack(a_iv);
}


// ---------------------------------------------------------
RealVect
RigidTransformation::transformVectorFwd(const RealVect& a_vecOld) const
{
  // Rotation only, no translation
  CH_assert(m_rotation == RealVect::Zero); // FIXME
  //CH_assert(SpaceDim == 2 || SpaceDim == 3);
  return a_vecOld;
}


// ---------------------------------------------------------
RealVect
RigidTransformation::transformVectorBack(const RealVect& a_vecNew) const
{
  // Inverse rotation only
  CH_assert(m_rotation == RealVect::Zero); // FIXME
  //CH_assert(SpaceDim == 2 || SpaceDim == 3);
  return a_vecNew;
}


// ---------------------------------------------------------
RealVect
RigidTransformation::transformVector(const RealVect& a_vec,
                                     bool a_forward) const
{
  if (a_forward)
    return transformVectorFwd(a_vec);
  else
    return transformVectorBack(a_vec);
}


// ---------------------------------------------------------
RigidTransformation
RigidTransformation::inverse() const
{
  return RigidTransformation(-m_rotation, -m_translation);
}


// ---------------------------------------------------------
RigidTransformation
RigidTransformation::compose(const RigidTransformation& a_next) const
{
  return RigidTransformation(m_rotation + a_next.m_rotation,
                             m_translation + a_next.m_translation);
}

/// Output the transformation
std::ostream&
operator<<(std::ostream& a_os, const RigidTransformation& a_trfm)
{
  a_os << "rot" << a_trfm.getRotation() << ",trn" << a_trfm.getTranslation();
  return a_os;
}

#include "NamespaceFooter.H"
