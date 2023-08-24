#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _RIGIDTRANSFORMATION_H_
#define _RIGIDTRANSFORMATION_H_

#include <iosfwd>

#include "RealVect.H"
#include "RealVect.H"
#include "Box.H"

#include "NamespaceHeader.H"

//**FIXME Currently this only supports translation.  For hops accross periodic
//**boundaries, we do not need rotation but do need translation in cylindrical
//**and spherical coordinates.
// ---------------------------------------------------------
/// Class to describe transformation of indices from one block to another.
class RigidTransformation
{
public:

  ///
  /*
   * Default is Identity transformation
   */
  RigidTransformation();

  ///
  /**
     Constructor.

     RealVect pOld:  indices of a cell in old index space
     RealVect pNew:  indices of the same cell in new index space

     Then for each direction idir:
     pNew[idir] = pOld[a_rotation[idir]] + a_translation[idir].
  */
  RigidTransformation(const RealVect& a_rotation,
                      const RealVect& a_translation);

  ///
  /*
    Same as constructor.
   */
  void define(const RealVect& a_rotation,
              const RealVect& a_translation);

  ///
  /*
   * Define as an identity transformation, for identical regions
   */
  void defineAsIdentity();

  bool operator==(const RigidTransformation& a_itOther) const;

  ///
  /**
     RealVect pOld:  indices of a cell in old index space
     RealVect pNew:  indices of the same cell in new index space

     Then for each direction idir:
     pNew[idir] = pOld[m_rotation[idir]] + m_translation[idir]
     and hence:
     pOld[m_rotation[idir]] = pNew[idir]
                               - m_translation[idir]
  */
  RealVect transformFwd(const RealVect& a_ivOld) const;

  RealVect transformBack(const RealVect& a_ivNew) const;

  RealVect transform(const RealVect& a_iv, bool a_forward = true) const;
  
  ///
  /**
     Transform a vector, not an index.
   */
  RealVect transformVectorFwd(const RealVect& a_vecOld) const;

  RealVect transformVectorBack(const RealVect& a_vecNew) const;

  RealVect transformVector(const RealVect& a_vec, bool a_forward = true) const;

  ///
  /**
     Return the inverse transformation.
   */
  RigidTransformation inverse() const;

  ///
  /**
     Return the composite transformation:
     that is, this transformation followed by a_next.
   */
  RigidTransformation compose(const RigidTransformation& a_next) const;

  ///
  /**
     Modify this to be the composite transformation:
     that is, this transformation followed by a_next.
   */
  RigidTransformation& operator+=(const RigidTransformation& a_next)
    {
      m_rotation += a_next.m_rotation;
      m_translation += a_next.m_translation;
      return *this;
    }

  /// Reset the transformation to zero
  void reset()
    {
      m_rotation = RealVect_zero;
      m_translation = RealVect_zero;
    }

  RealVect getRotation() const
  {
    return m_rotation;
  }
  RealVect getTranslation() const
  {
    return m_translation;
  }


protected:

  RealVect m_rotation;
  RealVect m_translation;
};

/// Composite, transform a_x followed by a_y
/** Right now, order does not matter but could be an issue if rotations are
 *  considered.
 */
inline
RigidTransformation operator+(const RigidTransformation& a_x,
                              const RigidTransformation& a_y)
{
  return RigidTransformation(a_x.getRotation() + a_y.getRotation(),
                             a_x.getTranslation() + a_y.getTranslation());
}

/// Remove transform a_y from a_x
/** Right now, order does not matter but could be an issue if rotations are
 *  considered.
 */
inline
RigidTransformation operator-(const RigidTransformation& a_x,
                              const RigidTransformation& a_y)
{
  return RigidTransformation(a_x.getRotation() - a_y.getRotation(),
                             a_x.getTranslation() - a_y.getTranslation());
}

/// Output the transformation
std::ostream&
operator<<(std::ostream& a_os, const RigidTransformation& a_trfm);

#include "NamespaceFooter.H"

#endif // include guard
