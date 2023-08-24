#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BOX_H_
#define _BOX_H_

#ifdef __CUDACC__
#define DEVICE __device__
#else
#define DEVICE
#endif  /* CUDA */

#ifdef NODEV
#undef NODEV
#endif

#include "SPACE.H"

#ifndef WRAPPER
#ifndef __CUDACC__
#include <iostream>
#endif  /* !CUDA */

#include "IntVect.H"
#include "Misc.H"
#include "LoHiSide.H"

#ifndef __CUDACC__
#include "Vector.H"
#include "SPMD.H"
#endif  /* !CUDA */

#include "BaseNamespaceHeader.H"
// Forward declaration
namespace CH_Hash
{
  template <typename T>
  struct isHashable;
}
#include "BaseNamespaceFooter.H"

#include "NamespaceHeader.H"

#ifdef USE_PROTO
#include "Proto.H"
#endif

struct SliceSpec;

/// Cell-Based or Node-Based Indices
/**

  The class IndexType defines an index as being cell based or node (face)
  based in each of the CH_SPACEDIM directions.  This class defines an
  enumerated type CellIndex to be either CELL or NODE; i.e. each of the
  CH_SPACEDIM dimensions must be either CELL or NODE.
*/

class IndexType
{
public:
  ///
  /**
     The cell index type: one of CELL or NODE.

  */
  enum CellIndex
  {
    CELL = 0,
    NODE = 1
  };

  ///
  /**
     The default constructor.

  */
  DEVICE constexpr IndexType () noexcept
    :
    itype(0)
    { }

  // Let the compiler synthesize some defaults
  IndexType(const IndexType&) = default;
  IndexType(IndexType&&) = default;
  IndexType& operator=(const IndexType&) = default;
  IndexType& operator=(IndexType&&) = default;
  ~IndexType() = default;

  ///
  /**
     Constructs an IndexType identical to an IntVect.

  */
  DEVICE explicit IndexType (const IntVect& iv);

  ///
  /**
     Another assignment operator.

  */
  DEVICE IndexType& operator= (const IntVect& iv);

  ///
  /**
      Constructs an IndexType given an explicit CellIndex for each
      direction.  D_DECL6 is a macro that sets the constructor to take
      BL_SPACEDIM arguments.

  */
  DEVICE IndexType (D_DECL6(CellIndex i, CellIndex j, CellIndex k,
                            CellIndex l, CellIndex m, CellIndex n));

  ///
  /**
     Sets this IndexType to be NODE based in direction dir.

  */
  DEVICE void set (int dir);

  ///
  /**
     Sets this IndexType to be CELL based in direction dir.

  */
  DEVICE void unset (int dir);

  ///
  /**
     True if this IndexType is NODE based in direction dir.

  */
  DEVICE bool test (int dir) const;

  ///
  /**
     Set this IndexType to be NODE based in all directions.

  */
  DEVICE void setall ();

  ///
  /**
     Set this IndexType to be CELL based in all directions.

  */
  DEVICE void clear ();

  ///
  /**
     True if this IndexType is NODE based in any direction.

  */
  DEVICE bool any () const;

  ///
  /**
     True if this IndexType is valid.

  */
  DEVICE bool ok () const;

  ///
  /**
     Change this IndexType from CELL to NODE or NODE to CELL
     in direction dir.

  */
  DEVICE void flip (int i);

  ///
  /**
     True if IndexTypes are identical.

  */
  DEVICE bool operator== (const IndexType& t) const;

  ///
  /**
    IndexTypes comparison.
  */
  DEVICE bool operator<(const IndexType& t) const;

  ///
  /**
     True if IndexTypes are not identical.

  */
  DEVICE bool operator!= (const IndexType& t) const;

  ///
  /**
     True if this IndexType is CELL based in all directions.

  */
  DEVICE bool cellCentered () const;

  ///
  /**
     True if this IndexType is NODE based in all directions.

  */
  DEVICE bool nodeCentered () const;

  ///
  /**
     Set this IndexType to CellIndex type t in direction dir.

  */
  DEVICE void setType (int dir, CellIndex t);

  ///
  /**
     Returns  the CellIndex of this IndexType in direction dir.

  */
  DEVICE CellIndex ixType (int dir) const;

  ///
  /**
     Return an integer representing this IndexType in direction dir.

  */
  DEVICE int operator[] (int dir) const;

  ///
  /**
     Convert this IndexType to an IntVect via component-wise
     conversion of enum CellIndex to int.

  */
  DEVICE IntVect ixType () const;

#ifndef __CUDACC__
  ///
  /**
      This static member function returns an IndexType object all of
      whose components are of value IndexType::CELL.  It is provided as
      a convenience to our users when defining a Box all of whose faces
      should be of type IndexType::CELL.

  */
  static IndexType TheCellType ();

  ///
  /**
     This static member function returns an IndexType object all of
     whose components are of value IndexType::NODE.  It is provided as
     a convenience to our users when defining a Box all of whose faces
     should be of type IndexType::NODE.

  */
  static IndexType TheNodeType ();

  ///
  /**
     Write an IndexType to an ostream in ASCII.

  */
  friend std::ostream& operator<< (std::ostream&         os,
                                   const IndexType& itype);

  ///
  /**
     Read an IndexType from an istream.

  */
  friend std::istream& operator>> (std::istream&   is,
                                   IndexType& itype);
#endif  /* !CUDA */

private:
  //
  // Returns 1<<k.
  //
  //  static int mask (int k);
  DEVICE static unsigned char mask (int k);

  //
  // An integer holding the CellIndex in bits 0 - BL_SPACEDIM-1.
  //
  //unsigned int itype;
  unsigned char itype;
};


//
// Inlines.
//

inline DEVICE unsigned char
IndexType::mask (int k)
{
  return 1<<k;
}

inline DEVICE
IndexType::IndexType (const IntVect& iv)
{
  itype = D_TERM6((iv[0]?1:0), | ((iv[1]?1:0)<<1), | ((iv[2]?1:0)<<2),
                  | ((iv[3]?1:0)<<3), | ((iv[4]?1:0)<<4), | ((iv[5]?1:0)<<5));
}

inline DEVICE
IndexType& IndexType::operator= (const IntVect& iv)
{
  itype = D_TERM6((iv[0]?1:0), | ((iv[1]?1:0)<<1), | ((iv[2]?1:0)<<2),
                  | ((iv[3]?1:0)<<3), | ((iv[4]?1:0)<<4), | ((iv[5]?1:0)<<5));
  return *this;
}

inline DEVICE
IndexType::IndexType (D_DECL6(CellIndex i, CellIndex j, CellIndex k,
                              CellIndex l, CellIndex m, CellIndex n))
{
  itype = D_TERM6(i, | (j<<1), | (k<<2),
           | (l<<3), | (m<<4), | (n<<5));
}

inline DEVICE void
IndexType::set (int dir)
{
  itype |= mask(dir);
}

inline DEVICE void
IndexType::unset (int dir)
{
  itype &= ~mask(dir);
}

inline DEVICE bool
IndexType::test (int dir) const
{
  return (itype & mask(dir)) != 0;
}

inline DEVICE void
IndexType::setall ()
{
  itype = (1 << CH_SPACEDIM) - 1;
}

inline DEVICE void
IndexType::clear ()
{
  itype = 0;
}

inline DEVICE bool
IndexType::any () const
{
  return itype != 0;
}

inline DEVICE bool
IndexType::ok () const
{
  return itype < (1 << CH_SPACEDIM);
}

inline DEVICE void
IndexType::flip (int i)
{
  itype ^= mask(i);
}

inline DEVICE bool
IndexType::operator== (const IndexType& t) const
{
  return t.itype == itype;
}

inline DEVICE bool
IndexType::operator< (const IndexType& t) const
{
  return itype < t.itype;
}

inline DEVICE bool
IndexType::operator!= (const IndexType& t) const
{
  return t.itype != itype;
}

inline DEVICE bool
IndexType::cellCentered () const
{
  return itype == 0;
}

inline DEVICE bool
IndexType::nodeCentered () const
{
  return itype == (1<<CH_SPACEDIM)-1;
}

inline DEVICE void
IndexType::setType (int       dir,
                    CellIndex t)
{
  t == CELL ? unset(dir) : set(dir);
}

inline DEVICE IndexType::CellIndex
IndexType::ixType (int dir) const
{
  return (CellIndex) ((itype & (1<<dir)) >> dir);
}

inline DEVICE int
IndexType::operator[] (int dir) const
{
  return test(dir);
}

inline DEVICE IntVect
IndexType::ixType () const
{
  return IntVect(D_DECL6(itype&1, (itype>>1)&1, (itype>>2)&1,
                         (itype>>3)&1, (itype>>4)&1, (itype>>5)&1));
}

#endif /* WRAPPER */

//
/// A Rectangular Domain on an Integer Lattice
/**

  A Box is an abstraction for defining discrete rectangular regions of
  SpaceDim-dimensioned indexing space.  Boxes have an IndexType, which
  defines IndexType::CELL or IndexType::NODE based points for each
  direction and a low and high IntVect which defines the lower and
  upper corners of the Box.  Boxes can exist in positive and negative
  indexing space.

  Box is a dimension dependent class, so SpaceDim must be
  defined as either 1, 2, or 3 when compiling.

  @ingroup chombo
*/

class Box
{
public:

  ///

  ///
  /**
     The default constructor.  The constructed Box is empty.

  */
  DEVICE Box ();

  // Let the compiler synthesize some defaults
  Box(const Box&) = default;
  Box(Box&&) = default;
  Box& operator=(const Box&) = default;
  Box& operator=(Box&&) = default;
  ~Box() = default;

  ///
  /**
     Construct cell-centered type Box.  It is an error if small is
     greater than big.

  */
  DEVICE constexpr Box(const IntVect& small, const IntVect& big) noexcept
    :
    smallend(small),
    bigend(big)
    {
      assert(small <= big);
    }

  ///
  /**
     Modify cell-centered type Box.  It is an error if small is
     greater than big.

  */
  DEVICE void define(const IntVect& small, const IntVect& big);

  ///
  /**
     Construct Box with specified lengths.  It is an error
     if the lengths are negative.

  */
  DEVICE Box (const IntVect& small, const int* vec_len);

  ///
  /**
      Construct Box with given type.  small and big are expected to be
      consistent with given type.  It is an error if small is greater
      than big.  It is an error if typ is not convertible to an
      IndexType.

  */
  DEVICE Box (const IntVect& small,
              const IntVect& big,
              const IntVect& typ);

  ///
  /**
      Modify Box with given type.  small and big are expected to be
      consistent with given type.  It is an error if small is greater
      than big.  It is an error if typ is not convertible to an
      IndexType.

  */
  DEVICE void define(const IntVect& small,
                     const IntVect& big,
                     const IntVect& typ);

  ///
  /**
      Construct Box with given type.  small and big are expected to be
      consistent with given type.  It is an error if small is greater
      than big.

  */
  DEVICE Box (const IntVect&   small,
              const IntVect&   big,
              const IndexType& t);

  ///
  /**
      Modify Box with given type.  small and big are expected to be
      consistent with given type.  It is an error if small is greater
      than big.

  */
  DEVICE void define(const IntVect&   small,
                     const IntVect&   big,
                     const IndexType& t);

  DEVICE void define(const Box& b);

#ifndef __CUDACC__
  Box copy() const
  {
    return *this;
  }
#endif  /* !CUDA */

  /// {\bf Accessors}

  ///
  /**
      Returns the lower corner of this Box.

  */
  DEVICE const IntVect& smallEnd () const;

  ///
  /**
     return smallend if side is low,
     bigend if side is high.
  */
  DEVICE IntVect sideEnd(Side::LoHiSide a_side) const;

  ///
  /**
     Returns the coordinate of the low end of this Box in the given
     direction.  Directions are zero-based.  It is an error if not 0
     <= dir < SpaceDim.

  */
  DEVICE int smallEnd (int dir) const;

  ///
  /**
     Returns the upper corner of this Box.

  */
  DEVICE const IntVect& bigEnd () const;

  ///
  /**
     Returns the coordinate of the high end of this Box in the given
     direction.  Directions are zero-based.  It is an error if not 0
     <= dir < SpaceDim.

  */
  DEVICE int bigEnd (int dir) const;

  ///
  /**
     Returns a constant pointer to the array of low end coordinates.
     Useful for calls to Fortran.  It should not be used in any other
     circumstances.

  */
  DEVICE const int* loVect () const;

  ///
  /**
     Returns a constant pointer the array of high end coordinates.
     Useful for calls to Fortran.  It should not be used in any other
     circumstances.

  */
  DEVICE const int* hiVect () const;

  ///
  /**
     Returns a constant pointer to the array of coordinates in the
     Box.  Useful for calls to Fortran, but otherwise too dangerous
     for use.

  */
  DEVICE const int* getVect () const;

  ///
  /**
     Returns offset of point from smallend; i.e.  index(smallend) ->
     0, bigend would return volume()-1.  Is used in accessing
     FArrayBox.

  */
  DEVICE long index (const IntVect& v) const;

  /** \name index() optimized for GPU.
   *
   *  index() caches the dimensions of the box.  If each thread did this, there
   *  would be significant register pressure.  The following routines allow
   *  for caching the stride (ideally to shared memory), and then each thread
   *  can use that to compute linear or vector indices as required.
   *
   *  Boxes are not zero-based (the smallEnd is not necessarily zero).
   *  Converting vector indices to linear indices must account for this.
   *  index() does this directly by subtracting the smallEnd.  We call this
   *  a box-based linear index and it can be used as (ignoring components) \code
fab.dataPtr()[lin]                                                      \endcode
   *
   *  In the following routines, the displacement of the smallEnd is separated
   *  into an offset.  The corresponding index routines, vecToLin, returns
   *  a zero-based index.  They can be used as \code
fab.dataPtr()[lin0 + offset]                                            \endcode
   *  However, it is probably more convenient to store \code
Real* p = fab.dataPtr() + offset
   *  and then use \code
p[lin0]                                                                 \endcode
   *  to get the data.
   *
   *  WhyTF???  All this avoids is excessive math for considering the smallEnd.
   *  Is it worth it or will people go insane?
   *
   *  IntVects and arrays representing an IntVect are stored as 32-bit integers.
   *  However, strides and linear indices are multiples of these values and
   *  higher precision may be needed.  Associated routines are templated to
   *  allow for flexible precision (but you cannot store anything in an
   *  IntVect which is always 32-bit int).
   */
  ///@{

  /// Get the strides
  /** \return             Strides in an IntVect
   */
  DEVICE IntVect strides () const
    {
      IntVect iv;
      strides<int>(iv.dataPtr());
      return iv;
    }

  /// Get the strides
  /** \tparam     I       Resolution of linear index
   *  \param[out] a_stride
   *                      Strides in an array
   */
  template <typename I>
  DEVICE void strides (I* const a_stride) const
    {
      D_TERM6(a_stride[0] = (I)1;,
              a_stride[1] = a_stride[0]*(I)(bigEnd(0) - smallEnd(0) + 1);,
              a_stride[2] = a_stride[1]*(I)(bigEnd(1) - smallEnd(1) + 1);,
              a_stride[3] = a_stride[2]*(I)(bigEnd(2) - smallEnd(2) + 1);,
              a_stride[4] = a_stride[3]*(I)(bigEnd(3) - smallEnd(3) + 1);,
              a_stride[5] = a_stride[4]*(I)(bigEnd(4) - smallEnd(4) + 1);)
    }

  /// Get the offset
  /** Add the offset to the pointer to get how far the pointer needs to be
   *  offset to account for a non-zero smallEnd.  If smallEnd is zero, this
   *  is a zero-based array and no offset is necessary.
   *  \param[in]  a_stride
   *                      Strides in an IntVect
   *  \return             Offset
   */
  DEVICE int offset (const IntVect& a_stride)
    {
      return vecToLin0<int>(IntVect(-smallEnd()).dataPtr(), a_stride.dataPtr());
    }

  /// Get the offset
  /** Add the offset to the pointer to get how far the pointer needs to be
   *  offset to account for a non-zero smallEnd.  If smallEnd is zero, this
   *  is a zero-based array and no offset is necessary.
   *  \tparam     I       Resolution of linear index
   *  \param[in]  a_stride
   *                      Strides in an array
   *  \return             Offset
   */
  template <typename I>
  DEVICE I offset (const I* const a_stride)
    {
      return vecToLin0<I>(IntVect(-smallEnd()).dataPtr(), a_stride);
    }

  /// Convert a vector index (in the box) to a zero-based linear index
  /** No consideration is made for a non-zero smallEnd.  This must instead be
   *  accounted for by the offset.  An offset data pointer
   * (fab.dataPtr() + offset) plus lin0 yields a valid index into the data.
   *  Alternatively a direct data pointer indexed with [lin + offset] yields a
   *  valid index into the data (fab.dataPtr() + lin0 + offset).
   *  \param[in]  a_iv    Vector index (contained in the box) represented as
   *                      an IntVect
   *  \param[in]  a_stride
   *                      Strides in an IntVect
   *  \return             Linear index (used with offset data pointer)
   */
  DEVICE int vecToLin0 (const IntVect& a_iv, const IntVect& a_stride)
    {
      return vecToLin0<int>(a_iv.dataPtr(), a_stride.dataPtr());
    }

  /// Convert a vector index (in the box) to a zero-based linear index
  /** No consideration is made for a non-zero smallEnd.  This must instead be
   *  accounted for by the offset.  An offset data pointer
   *  (fab.dataPtr() + offset) plus lin0 yields a valid index into the data.
   *  Alternatively a direct data pointer indexed with [lin + offset] yields a
   *  valid index into the data (fab.dataPtr() + lin0 + offset).
   *  \tparam     I       Resolution of linear index
   *  \param[in]  a_vec   Vector index (contained in the box) represented as
   *                      an array
   *  \param[in]  a_stride
   *                      Strides in an array
   *  \return             Linear index (used with offset data pointer)
   */
  template <typename I>
  DEVICE I vecToLin0 (const int* a_vec, const I* const a_stride)
    {
      return D_TERM6(  (I)a_vec[0],
                     + (I)a_vec[1]*a_stride[1],
                     + (I)a_vec[2]*a_stride[2],
                     + (I)a_vec[3]*a_stride[3],
                     + (I)a_vec[4]*a_stride[4],
                     + (I)a_vec[5]*a_stride[5]);
    }

  /// Convert a box-based linear index (in the box) to a vector index
  /** If you have a zero-based linear index (see vecToLin0), you must add the
   *  offset before calling this routine.  E.g.
   *  linToVec(lin0 + offset, strides, vec).
   *  \tparam     I       Resolution of linear index
   *  \param[in]  a_lin   A box-based linear index (add offset if you have
   *                      a zero-based index)
   *  \param[in]  a_stride
   *                      Strides in an IntVect
   *  \param[out] a_vec   Vector index (contained in the box) represented as
   *                      an IntVect
   */
  DEVICE void linToVec (int a_lin, const IntVect& a_stride, IntVect& a_vec)
    {
      linToVec<int>(a_lin, a_stride.dataPtr(), a_vec.dataPtr());
    }

  /// Convert a box-based linear index (in the box) to a vector index
  /** If you have a zero-based linear index (see vecToLin0), you must add the
   *  offset before calling this routine.  E.g.
   *  linToVec(lin0 + offset, strides, vec).
   *  \tparam     I       Resolution of linear index
   *  \param[in]  a_lin   A box-based linear index (add offset if you have
   *                      a zero-based index)
   *  \param[in]  a_stride
   *                      Strides in an array
   *  \param[out] a_vec   Vector index (contained in the box) represented as
   *                      an array
   */
  template <typename I>
  DEVICE void linToVec (I a_lin, const I* const a_stride, int* a_vec)
    {
      //**FIXME Workaround to avoid compiler warnings in 1-D
      // is there a better solution? (SEP)
#if CH_SPACEDIM > 1
      I tmp;
#endif
      D_INVTERM6(a_vec[0] = a_lin + smallEnd(0);,

                 tmp = a_lin/a_stride[1];
                 a_vec[1] = tmp + smallEnd(1);
                 a_lin -= tmp*a_stride[1];,

                 tmp = a_lin/a_stride[2];
                 a_vec[2] = tmp + smallEnd(2);
                 a_lin -= tmp*a_stride[2];,

                 tmp = a_lin/a_stride[3];
                 a_vec[3] = tmp + smallEnd(3);
                 a_lin -= tmp*a_stride[3];,

                 tmp = a_lin/a_stride[4];
                 a_vec[4] = tmp + smallEnd(4);
                 a_lin -= tmp*a_stride[4];,

                 tmp = a_lin/a_stride[5];
                 a_vec[5] = tmp + smallEnd(5);
                 a_lin -= tmp*a_stride[5];)
    }
  
  ///@}
  
  /// centering

  ///
  /**
     Return the indexing type of this Box.

  */
  DEVICE IndexType ixType () const;

  ///
  /**
     Return the indexing type of this Box.

  */
  DEVICE IntVect type () const;

  ///
  /**
     Return the indexing type of this Box in the specified direction.
     Directions are zero-based.  It is an error if not 0 <= dir <
     SpaceDim.

  */
  DEVICE IndexType::CellIndex type (int dir) const;

  /// size functions

  ///
  /**
     Return an IntVect containing the size of this Box in each
     coordinate direction.

  */
  DEVICE IntVect size () const;

  ///
  /**
     Return the size of this Box in the specified coordinate
     direction.  Directions are zero-based.  It is an error
     if not 0 <= dir < SpaceDim.

  */
  DEVICE int size (int dir) const;

 
  ///
  /**
     Return the volume, in indexing space, of region enclosed by this
     Box.  This is identical to volume().  

  */
  DEVICE size_t numPts () const;

 
  ///
  /**
     Return the volume, in indexing space, of region enclosed by this
     Box.  This is identical to numPts().  

  */
  DEVICE size_t volume () const;

  ///
  /**
     Returns length of the longest side of this Box.  The argument dir
     is modified to give direction with longest side:
     0...SpaceDim-1. Ignores type.

  */
  DEVICE int longside (int& dir) const;

  ///
  /**
     Returns length of the longest side of this Box.  Ignores type.

  */
  DEVICE int longside () const;

  ///
  /**
     Returns length of shortest side of this Box.  The argument dir is
     modified to give direction with shortest side: 0...SpaceDim-1.
     Ignores type.

  */
  DEVICE int shortside (int& dir) const;

  ///
  /**
     Returns length of the shortest side of this Box.  Ignores type.

  */
  DEVICE int shortside () const;

  /// {\bf Comparison Functions}

  ///
  /**
     Returns true if this Box is empty.

  */
  DEVICE bool isEmpty () const;

  ///
  /**
     Returns true if this Box is not empty and all
     sizes are right
  */
  DEVICE bool ok() const;

  ///
  /**
     Returns true if argument is contained within this Box.  An empty
     Box does not contain and is not contained by any Box, including
     itself.

  */
  DEVICE bool contains (const IntVect& p) const;

  ///
  /**
     Returns true if argument is contained within this Box.  It is an
     error if the Boxes have different types.  An empty Box does not
     contain any IntVect.

  */
  DEVICE bool contains (const Box& b) const;

  ///
  /**
     Returns true if this Box and the argument have non-null
     intersections.  It is an error if the Boxes have different types.
     An empty Box does not intersect any Box, including itself.

  */
  DEVICE bool intersects (const Box& b) const;

  ///
  /**
     Returns true if this Box and the argument have non-null
     intersections.  It is an error if the Boxes have different types.
     This routine does not perform the check to see if *this or b are
     empty boxes.  It is the callers responsibility to ensure that
     this never happens.  If you are unsure, the use the .intersects(..) routine.

  */
  DEVICE bool intersectsNotEmpty (const Box& b) const;

  ///
  /**
     Returns true if this Box and the argument are the same size, ie
     translates of each other.  It is an error if they have different
     types.
  */
  DEVICE bool sameSize (const Box& b) const;

  ///
  /**
     Returns true if this Box and the argument have same type.

  */
  DEVICE bool sameType (const Box &b) const;

  ///
  /**
     Returns true if this Box and the argument are identical,
     including type.

  */
  DEVICE bool operator== (const Box& b) const;

  DEVICE bool eq(const Box& b) const;
  ///
  /**
     Returns true if this Box and the argument differ, including type.

  */
  DEVICE bool operator!= (const Box& b) const;

  DEVICE bool neq(const Box& b) const;
  ///
  /**
     Returns true if this Box is cell-centered in all indexing
     directions.

  */
  DEVICE bool cellCentered () const;

  // following operators added Sept. 14, 1999.  bvs
  ///
  /**
     Returns true if this Box is lexigraphically less than rhs box.
     All comparison is based on lower box corner.  In the name of
     coding efficiency, we do not handle the case where two boxes have
     the same lower left corner. In DEBUG mode, this is checked for
     with an assert.

  */
  DEVICE bool operator < (const Box& rhs) const;

  DEVICE bool lt(const Box& rhs) const;

  /// {\bf Modification Functions}

  ///
  /**
     Redefine the lower corner of this Box.  It is an error
     if the specified corner is greater than the exisiting upper
     corner of this Box.

  */
  DEVICE Box& setSmall (const IntVect& sm);

  ///
  /**
     Redefines the lower end of this Box in the specified coordinate
     direction.  It is an error if the specified value is greater than
     the exisiting upper end of this Box.  Directions are zero-based.
     It is an error if not 0 <= dir < SpaceDim.

  */
  DEVICE Box& setSmall (int dir, int sm_index);

  ///
  /**
     Redefines the upper corner of this Box.  It is an error if the
     specified corner is less than the exisiting lower corner of this
     Box.

  */
  DEVICE Box& setBig (const IntVect& bg);

  ///
  /**
     Redefines the upper end of this Box in the specified coordinate
     direction.  It is an error if the specified value is less than the
     exisiting lower end of this Box.  Directions are zero-based.  It
     is an error if not 0 <= dir < SpaceDim.

  */
  DEVICE Box& setBig (int dir, int bg_index);

  ///
  /**
     Set the entire range in a given direction, starting at sm_index
     with length n_cells.  It is an error if n_cells <= 0.

  */
  DEVICE Box& setRange (int dir,
                        int sm_index,
                        int n_cells = 1);

  /// centering type conversion functions

  ///
  /**
     Modifies this Box by converting from the current type into the
     argument type.  This may change the Box coordinates:\\

     type CELL -> NODE : increase coordinate by one on high end.\\

     type NODE -> CELL : reduce coordinate by one on high end.\\

     other type mappings make no change.

  */
  DEVICE Box& convert (IndexType typ);

  ///
  /**
     Modifies this Box by converting from the current type into the
     argument type.  This may change the Box coordinates:\\ type CELL
     -> NODE : increase coordinate by one on high end.\\ type NODE ->
     CELL : reduce coordinate by one on high end.\\ other type
     mappings make no change.  It is an error if typ is not
     convertible to an IndexType.

  */
  DEVICE Box& convert (const IntVect& typ);

  ///
  /**
     Modifies this Box by converting from the current type into the
     argument type.  This may change the Box coordinates:\\

     type CELL -> NODE : increase coordinate by one on high end.\\

     type NODE -> CELL : reduce coordinate by one on high end.\\

     Other type mappings make no change.  Directions are zero-based.
     It is an error if not 0 <= dir < SpaceDim.

  */
  DEVICE Box& convert (int dir, IndexType::CellIndex typ);

  ///
  /**
     Modifies this Box by converting it to NODE type in all
     directions.  This increases all components of the upper corner by
     one.  The Empty Box is not modified by this function except for
     changing of type.

  */
  DEVICE Box& surroundingNodes ();

  ///
  /**
     Modifies this Box by converting it to NODE type in given
     direction.  This increases the component of the upper corner in
     the given direction by one.  Directions are zero-based.  It is an
     error if not 0 <= dir < SpaceDim.  The Empty Box is not modified
     by this function except for changing of type.

  */
  DEVICE Box& surroundingNodes (int dir);

  DEVICE Box& surroundingNodes_int(int dir);

  ///
  /**
     Modifies this Box by converting it to CELL type in all directions.
     This decreases all components of the upper corner by one.  The
     Empty Box is not modified by this function except for changing of
     type.

  */
  DEVICE Box& enclosedCells ();

  ///
  /**
     Modifies this Box by converting it to CELL type in given
     direction.  This decreases the component of the upper corner in
     the given direction by one.  Directions are zero-based.  It is an
     error if not 0 <= dir < SpaceDim.  The Empty Box is not modified
     by this function except for changing of type.

  */
  DEVICE Box& enclosedCells (int dir);

  DEVICE Box& enclosedCells_int (int dir);

  /// shift functions

  ///
  /**
     Modifies this Box by shifting it nzones indexing positions in
     coordinate direction dir.  Directions are zero-based.  It is an
     error if not 0 <= dir < SpaceDim.  The Empty Box is not modified
     by this function.

  */
  DEVICE Box& shift (int dir, int nzones);

  ///
  /**
     Modifies this Box by shifting.  Equivalent to
     b.shift(0,iv[0]).shift(1,iv[1]) ....  The Empty Box is not
     modified by this function.

  */
  DEVICE Box& shift (const IntVect& iv);

  DEVICE Box& shift_intvect (const IntVect& iv);

  ///
  /**
     This member modifies this Box by shifting by "half" indices,
     thereby converting the Box from type CELL to NODE or vice-versa.
     b.shiftHalf(0,1) shifts b to the right by 1/2 cells.
     b.shiftHalf(1,-3) shifts b in the -j direction by 3/2 cells.
     NOTE: If num_halfs is EVEN the shift is num_halfs/2 full zones
     and hence will not change the type.  This is: b.shifthalf(4) ==
     b.shift(2).  Directions are zero-based.  It is an error if not 0
     <= dir < SpaceDim.  The Empty Box is not modified by this
     function except for changing of type.

  */
  DEVICE Box& shiftHalf (int dir, int num_halfs);

#ifndef __CUDACC__
  ///multiblock stuff.
  void convertOldToNew(const IntVect& a_permutation,
                       const IntVect& a_sign,
                       const IntVect& a_translation);


  ///multiblock stuff
  void convertNewToOld(const IntVect& a_permutation,
                       const IntVect& a_sign,
                       const IntVect& a_translation);
#endif

  ///
  /**
     Modifies this Box by shifting by half indices.  Equivalent to
     b.shiftHalf(0,iv[0]).shiftHalf(1,iv[1]) ...  The Empty Box is not
     modified by this function except for changing of type.

  */
  DEVICE Box& shiftHalf (const IntVect& iv);
  DEVICE Box& shiftHalf_intvect (const IntVect& iv);

  ///
  /**
     Modifies this Box by shifting it by given IntVect.  The Empty Box
     is not modified by this function.

  */
  DEVICE Box& operator+= (const IntVect& v);

  ///
  /**
     Returns a Box that is this Box shifted by the given IntVect.  The
     Empty Box is not modified by this function.

  */
  DEVICE Box operator+ (const IntVect& v) const;

  ///
  /**
     Modifies this Box by shifting it by given IntVect.  The Empty Box
     is not modified by this function.

  */
  DEVICE Box& operator-= (const IntVect& v);

  ///
  /**
     Returns a Box that is this Box shifted by the given IntVect.  The
     Empty Box is not modified by this function.

  */
  DEVICE Box operator- (const IntVect& v) const;

  /// neighbor box functions

  ///
  friend Box bdryBox(const Box& b,
                     int        dir,
                     Side::LoHiSide a_sd,
                     int        len);

  ///
  /**
     Returns the edge-centered Box (in direction dir) defining the low
     side of the argument Box.  The output Box will have the given
     length in the given direction.  Directions are zero-based.  It is
     an error if not 0 <= dir < SpaceDim.  The neighbor of an Empty Box
     is an Empty Box of the appropriate type.

  */
  friend Box bdryLo (const Box& b,
                     int        dir,
                     int        len);

  ///
  /**
     Returns the edge-centered Box (in direction dir) defining the
     high side of the argument Box.  The return Box will have the
     given length in the given direction.  Directions are zero-based.
     It is an error if not 0 <= dir < SpaceDim.  The neighbor of an
     Empty Box is an Empty Box of the appropriate type.

  */
  friend Box bdryHi (const Box& b,
                     int        dir,
                     int        len);

  ///
  /**
     Returns the cell centered Box of the given length adjacent to the
     argument Box on the low end along the given coordinate direction.
     The return Box is identical to the argument Box in the other
     directions.  The return Box and the argument Box have an empty
     intersection.  \\

     NOTE: len != 0. A negative length results in cells _inside_ the original box. \\

     NOTE: Box retval = adjCellLo(b,dir,len) is equivalent to the
     following set of operations: \\

     Box retval(b); \\

     retval.convert(dir,Box::CELL); \\

     retval.setrange(dir,retval.smallEnd(dir)-len,len); \\

     Directions are zero-based.  It is an error if not 0 <= dir <
     SpaceDim.  The neighbor of an Empty Box is an Empty Box of the
     appropriate type.

  */
  friend Box adjCellLo (const Box& b,
                        int        dir,
                        int        len);

  ///
  /**
     Returns the cell centered Box of the given length adjacent to the
     argument Box on the high end along the given coordinate
     direction.  The return Box is identical to the argument Box in
     the other directions.  The return Box and the argument Box have
     an empty intersection.  \\

     NOTE: len != 0. A Negative length results in cells _inside_ the original box. \\

     NOTE: Box retval = adjCellHi(b,dir,len) is equivalent to the
     following set of operations: \\

     Box retval(b); \\

     retval.convert(dir,Box::CELL); \\

     retval.setrange(dir,retval.bigEnd(dir)+1,len);\\

     Directions are zero-based.  It is an error if not 0 <= dir <
     SpaceDim.  The neighbor of an Empty Box is an Empty Box of the
     appropriate type.

  */
  friend Box adjCellHi (const Box& b,
                        int        dir,
                        int        len);

  /// intersection functions

  ///
  /**
     Returns the Box that is intersection of this Box and the argument
     Box.  Intersection is commutative.  The Boxes MUST be of same
     type.  The intersection of the Empty Box and any Box is the Empty
     Box.

  */
  DEVICE Box operator& (const Box&) const;

  ///
  /**
     Modifies this Box by intersection with the argument Box.  The
     Boxes MUST be of the same type.  The intersection of the Empty
     Box and any Box is the Empty Box.

  */
  DEVICE Box& operator&= (const Box&);

  ///
  friend Box adjCellBox (const Box& b,
                         int        dir,
                         Side::LoHiSide a_side,
                         int        len);
  ///
  /**
     Modifies this Box to that of the minimum Box containing both this
     Box and the argument Box.  Both Boxes must have identical type.

  */
  DEVICE Box& minBox (const Box& b);

  ///
  /**
     Returns a Box that is the minimum Box containing both the
     argument Boxes.  Both Boxes must have identical type.

  */
  DEVICE friend Box minBox (const Box& b1,
                            const Box& b2);

  /// grow functions

  ///
  /**
     Modifies this Box by growing it in all directions by given
     amount.  The Empty Box is not modified by this function.  NOTE:
     n_cell negative shrinks the Box by that number of cells.  If
     shrinking produces an empty Box, the Box is transformed into the
     canonical Empty Box.

  */
  DEVICE Box& grow (int i);

  ///
  /**
     Returns a Box that is the argument Box grown in all directions by
     given amount.  The Empty Box is not modified by this function.
     NOTE: n_cell negative shrinks the Box by that number of cells.
     If shrinking produces an empty Box, the Box is transformed into
     the canonical Empty Box.

  */
  friend Box grow (const Box& b,
                   int        i);

  ///
  /**
     Modifies this Box by growing it in each direction by specified
     amount.  The Empty Box is not modified by this function.  NOTE:
     components of iv may be negative, which would shrink this Box.  If
     shrinking produces an empty Box, the Box is transformed into the
     canonical Empty Box.

  */
  DEVICE Box& grow (const IntVect& v);

  ///
  /**
     Returns a Box that is the argument Box grown in each
     direction by specified amount.  The Empty Box is not modified by
     this function.  NOTE: components of iv may be negative, which
     would return a shrunken Box.  If shrinking produces an empty
     Box, the Box is transformed into the canonical Empty Box.

  */
  friend  Box grow (const Box&     b,
                    const IntVect& v);

  ///
  /**
     Modifies this Box by growing it on the low and high end by n_cell
     cells in direction idir.  The Empty Box is not modified by this
     function.  NOTE: n_cell negative shrinks this Box by that number
     of cells.  If shrinking produces an empty Box, the Box is
     transformed into the canonical Empty Box.  Directions are
     zero-based.  It is an error if not 0 <= dir < SpaceDim.

  */
  DEVICE Box& grow (int idir, int n_cell);

  ///
  /**
     Modifies this Box by growing it on the low end by n_cell cells
     in direction idir.  The Empty Box is not modified by this
     function.  NOTE: n_cell negative shrinks this Box by that number
     of cells.  If shrinking produces an empty Box, the Box is
     transformed into the canonical Empty Box.  Directions are
     zero-based.  It is an error if not 0 <= dir < SpaceDim.

  */
  DEVICE Box& growLo (int idir, int n_cell=1);

  ///
  /**
     growLo if a_sd == Side::Lo \\
     growHi if a_sd == Side::Hi
   */
  DEVICE Box& growDir (int a_idir,
                       const Side::LoHiSide& a_sd,
                       int a_cell);

  ///
  /**
     Modifies this Box by growing it on the high end by n_cell cells
     in direction idir.  The Empty Box is not modified by this
     function.  NOTE: n_cell negative shrinks the Box by that number
     of cells.  If shrinking produces an empty Box, the Box is
     transformed into the canonical Empty Box.  Directions are
     zero-based.  It is an error if not 0 <= dir < SpaceDim.

  */
  DEVICE Box& growHi (int idir, int n_cell=1);

  /// refinement

  ///
  /**
     Modifies this Box by refining it by given (positive) refinement
     ratio.  The Empty Box is not modified by this function.  \\

     NOTE: if type(dir) = CELL centered: lo <- lo*ratio and hi <-
     (hi+1)*ratio - 1. \\

     NOTE: if type(dir) = NODE centered: lo <- lo*ratio and hi <-
     hi*ratio.

  */
  DEVICE Box& refine (int refinement_ratio);

  ///
  /**
     Returns a Box that is the argument Box refined by given
     (positive) refinement ratio.  The Empty Box is not modified by
     this function.  \\

     NOTE: if type(dir) = CELL centered: lo <- lo*ratio and hi <-
     (hi+1)*ratio - 1.\\

     NOTE: if type(dir) = NODE centered: lo <- lo*ratio and hi <-
     hi*ratio.

  */
  friend Box refine (const Box& b,
                     int        refinement_ratio);

  ///
  /**
     Modifies this Box by refining it by given (positive) refinement
     ratio.  The Empty Box is not modified by this function.  \\

     NOTE: if type(dir) = CELL centered: lo <- lo*ratio and hi <-
     (hi+1)*ratio - 1.\\

     NOTE: if type(dir) = NODE centered: lo <- lo*ratio and hi <-
     hi*ratio.

  */
  DEVICE Box& refine (const IntVect& refinement_ratio);

  ///
  /**
     Returns a Box that is the argument Box refined by given
     (positive) refinement ratio.  The Empty Box is not modified by
     this function.  \\

     NOTE: if type(dir) = CELL centered: lo <- lo*ratio and hi <-
     (hi+1)*ratio - 1.\\

     NOTE: if type(dir) = NODE centered: lo <- lo*ratio and hi <-
     hi*ratio.

  */
  friend Box refine (const Box&     b,
                     const IntVect& refinement_ratio);

  /// coarsening

  /*
      True if the if (refine(coarsen(int))) gets the same box back
  */
  DEVICE bool coarsenable(int refrat) const
    {
      Box test = *this;
      test.coarsen(refrat);
      test.refine(refrat);
      return (test==*this);
    }

  ///
  /**
     Modifies this Box by coarsening it by given (positive) refinement
     ratio.  The Empty Box is not modified by this function.  \\

     NOTE: if type(dir) = CELL centered: lo <- lo/ratio and hi <-
     hi/ratio.\\

     NOTE: if type(dir) = NODE centered: lo <- lo/ratio and hi <-
     hi/ratio + ((hi%ratio)==0 ? 0 : 1).\\

     That is, refinement of coarsened Box must contain the original
     Box.

  */

  DEVICE Box& coarsen (int refinement_ratio);

  ///
  /**
     Returns a Box that is the argument Box coarsened by given
     (positive) refinement ratio.  The Empty Box is not modified by
     this function.  \\

     NOTE: if type(dir) = CELL centered: lo <- lo/ratio and hi <-
     hi/ratio.\\

     NOTE: if type(dir) = NODE centered: lo <- lo/ratio and hi <-
     hi/ratio + ((hi%ratio)==0 ? 0 : 1).\\

     That is, refinement of coarsened Box must contain the original
     Box.

  */
  friend Box coarsen (const Box& b,
                      int        refinement_ratio);

  ///
  /**
     Modifies this Box by coarsening by given (positive) refinement
     ratio.  The Empty Box is not modified by this function.  \\

     NOTE: if type(dir) = CELL centered: lo <- lo/ratio and hi <-
     hi/ratio.\\

     NOTE: if type(dir) = NODE centered: lo <- lo/ratio and hi <-
     hi/ratio + ((hi%ratio)==0 ? 0 : 1).\\

     That is, refinement of coarsened Box must contain the original
     Box.

  */
  DEVICE Box& coarsen (const IntVect& refinement_ratio);

  ///
  /**
     Returns a Box that is the argument Box coarsened by given
     (positive) refinement ratio.  The Empty Box is not modified by
     this function.  \\

     NOTE: if type(dir) = CELL centered: lo <- lo/ratio and hi <-
     hi/ratio.\\

     NOTE: if type(dir) = NODE centered: lo <- lo/ratio and hi <-
     hi/ratio + ((hi%ratio)==0 ? 0 : 1).\\

     That is, refinement of coarsened Box must contain the original
     Box.

  */
  friend Box coarsen (const Box&     b,
                      const IntVect& refinement_ratio);

  // next(...) is out of favor.  use BoxIterator.
  /*
    Step through the rectangle.  It is a runtime error to give
    a point not inside rectangle.  Iteration may not be efficient.
  */
#ifndef __CUDACC__
  DEVICE void next (IntVect &) const;
#endif  /* !CUDA */

  /*
    Scan argument IntVect over object second arg is
    increment vector.  Runtime error if IntVect is not
    contained in object Box.  Iteration may not be efficient.
  */
#ifndef __CUDACC__
  DEVICE void next (IntVect& p, const int* shv) const;
#endif  /* !CUDA */

  /// misc

  ///
  /**
     Chops the Box at the chop_pnt in the dir direction returns one
     Box, modifies the object Box.  The union of the two is the
     original Box.  The modified Box is the low end, the returned Box
     is the high end.  If type(dir) = CELL, the Boxes are disjoint with
     the chop_pnt included in the high end (new Box).  It is an ERROR
     if chop_pnt is the low end of the orig Box.  If type(dir) = NODE,
     the chop_pnt is included in both Boxes but is the only point in
     common.  It is also an error if the chop_pnt is an end node of
     the Box.  Directions are zero-based.  It is an error if not 0 <=
     dir < SpaceDim.

  */
#ifndef __CUDACC__
  DEVICE Box chop (int dir, int chop_pnt);
#endif  /* !CUDA */

  /**
    Makes a_to a box that's the same as *this except that it's just one cell
    thick in the a_sliceSpec.direction-th direction, and its (one) coordinate
    in that direction is a_sliceSpec.position.

    If a_sliceSpec.position is outside the range of *this, and a_outofbounds
    isn't null, then we set *a_outofbounds to true, otherwise to false.

    It would be more natural to return a Box, than to set the value of a Box
    reference, but we do it this way to be consistent with BaseFab::degenerate()
    where we can't return an object because of the policy of disallowing copy
    constructors for heavyweight classes.
  */
#ifndef __CUDACC__
  DEVICE void degenerate( Box&             a_to,
                          const SliceSpec& a_sliceSpec,
                          bool*            a_outofbounds=0 ) const;
#endif  /* !CUDA */

  ///{\bf I/O Functions}

  ///
  /**
     Write an ASCII representation to the ostream.

  */
#ifndef __CUDACC__
  friend std::ostream& operator<< (std::ostream&   os,
                                   const Box& bx);
#endif  /* !CUDA */

  ///
  /**
     Read from istream.

  */
#ifndef __CUDACC__
  friend std::istream& operator>> (std::istream& os,
                                   Box&     bx);
#endif  /* !CUDA */

  ///
  /**
     print to pout()
  */
#ifndef __CUDACC__
  void p() const;
#endif  /* !CUDA */

  ///
  /**
     Gives more detail than printOn.  Useful for exiting due to an
     error.

  */
#ifndef __CUDACC__
  void dumpOn (std::ostream& strm) const;
#endif  /* !CUDA */

  /// {\bf Box Constants}

  ///
  /**
     This is a canonical empty Box of cell-centered type.

  */
  //static const Box Empty;

  // TheUnitBox is out of favor.
  /*
    This static member function returns a constant reference to
    an object of type Box representing the unit box in
    BL_SPACEDIM-dimensional space.
  */

  //static const Box& TheUnitBox ();

  //
  // Sets the 'len' element of the Box.  Aborts on integer overflow.
  //
  DEVICE void computeBoxLen ();
  DEVICE void computeBoxLenNotEmpty();

#ifndef __CUDACC__
  /** Doesn't print the IndexType, and uses square outer brackets. */
  static void setTempestOutputFormat( bool );
#endif  /* !CUDA */

#ifdef USE_PROTO
    Box(Proto::Box a_box)
    {
        define(a_box.low(), a_box.high());
    }
    
    operator Proto::Box() const
    {
        return Proto::Box(smallend, bigend); //calls IntVect -> Proto::Point converter
    }
#endif

//--Data members

protected:
  friend class HDF5Handle;

  IntVect   smallend;
  IntVect   bigend;
  //IntVect   len;
  IndexType btype;

#ifndef __CUDACC__
  static bool s_tempestOutputFormat;
#endif  /* !CUDA */
};

#ifndef __CUDACC__
#include "NamespaceFooter.H"

#include "BaseNamespaceHeader.H"

//**FIXME I think you can add a scope and include UsingNamespace.H instead of
//**      NamespaceVar.H.  At least, IntVect.H gets away with this.
#include "NamespaceVar.H"

//Box specialization of linearSize
template < >
int linearSize(const CH_XDIR::Box& a_input);

//Box specialization of linearIn
template < >
void linearIn(CH_XDIR::Box& a_outputT, const void* const a_inBuf);

//Box specialization of linearOut
template < >
void linearOut(void* const a_outBuf, const CH_XDIR::Box& a_inputT);


//Vector<Box>  specialization
template < >
int linearSize(const Vector<CH_XDIR::Box>& a_input);
template < >
void linearIn(Vector<CH_XDIR::Box>& a_outputT, const void* const inBuf);
template < >
void linearOut(void* const a_outBuf, const Vector<CH_XDIR::Box>& a_inputT);

//Vector<Vector<Box> > specialization
template < >
int linearSize(const Vector<Vector<CH_XDIR::Box> >& a_input);
template < >
void linearIn(Vector<Vector<CH_XDIR::Box> >& a_outputT, const void* const inBuf);
template < >
void linearOut(void* const a_outBuf, const Vector<Vector<CH_XDIR::Box> >& a_inputT);

namespace CH_Hash
{

/// Remark that Box can be byte-hashed and give bytes to hash
template <>
struct isHashable<CH_XDIR::Box> : public std::true_type
{
  // Assert IndexType is POD
  static_assert(sizeof(CH_XDIR::IndexType) == sizeof(unsigned char),
                "Cannot guarantee that IndexType is not padded");
  static constexpr int c_hashSize =
    2*isHashable<CH_XDIR::IntVect>::c_hashSize +
    sizeof(CH_XDIR::IndexType);
  // It is also runtime-asserted that Box is not padded in Box::define
};

}

#include "BaseNamespaceFooter.H"
#include "NamespaceHeader.H"
#endif  /* !CUDA */

// global function prototypes
DEVICE Box surroundingNodes (const Box& b,
                            int        dir);
DEVICE Box surroundingNodes (const Box& b);
DEVICE Box enclosedCells (const Box& b,
                          int        dir);
DEVICE Box enclosedCells (const Box& b);
DEVICE Box bdryBox (const Box& b,
                    int        dir,
                    Side::LoHiSide a_sd,
                    int        len=1);
DEVICE Box bdryLo (const Box& b,
                   int        dir,
                   int        len=1);
DEVICE Box bdryHi (const Box& b,
                   int        dir,
                   int        len=1);
DEVICE Box adjCellLo (const Box& b,
                      int        dir,
                      int        len=1);
DEVICE Box adjCellHi (const Box& b,
                      int        dir,
                      int        len=1);
DEVICE Box minBox (const Box& b1,
                   const Box& b2);
DEVICE Box coarsen (const Box&     b,
                    const IntVect& refinement_ratio);
DEVICE Box coarsen (const Box&     b,
                    int refinement_ratio);
DEVICE Box refine (const Box&     b,
                   const IntVect& refinement_ratio);
DEVICE Box refine (const Box&     b,
                   int refinement_ratio);

//
// Inlines.
//

#ifndef WRAPPER

//--For CUDA, constructors are inline.  For CPU, constructors are in Box.cpp

#ifdef __CUDACC__
inline DEVICE
Box::Box ()
  :
  smallend(D_DECL6(1, 1, 1, 1, 1, 1)),
  bigend(D_DECL6(0, 0, 0, 0, 0, 0)),
  btype()
{
}

inline DEVICE
Box::Box (const IntVect& small,
          const IntVect& big)
  :
  smallend(small),
  bigend(big)
{
  CH_assert (small <= big);
  // computeBoxLen();
}

inline DEVICE void
Box::define (const IntVect& small,
             const IntVect& big)
{
  smallend = small;
  bigend = big;
  CH_assert (small <= big);
  // computeBoxLen();
  // Ensure no padding in Box (supports byte hashing).  This assert can be
  // anywhere but we chose to place it here since define cannot be constexpr
  CH_assert(reinterpret_cast<char*>(this) +
            CH_Hash::isHashable<Box>::c_hashSize ==
            reinterpret_cast<char*>(&btype) + sizeof(IndexType));
}

inline DEVICE
Box::Box (const IntVect& small,
          const int*     vec_len)
  :
  smallend(small),
  bigend(D_DECL6(small[0]+vec_len[0]-1,
                 small[1]+vec_len[1]-1,
                 small[2]+vec_len[2]-1,
                 small[3]+vec_len[3]-1,
                 small[4]+vec_len[4]-1,
                 small[5]+vec_len[5]-1))
{
  CH_assert(D_TERM6(vec_len[0] >= 0, && vec_len[1] >= 0, && vec_len[2] >= 0, &&
                    vec_len[3] >= 0, && vec_len[4] >= 0, && vec_len[5] >= 0));
}

inline DEVICE
Box::Box (const IntVect& small,
          const IntVect& big,
          const IntVect& typ)
  :
  smallend(small),
  bigend(big),
  btype(typ)
{
  CH_assert (small <= big);
  CH_assert(typ >= IntVect(D_DECL6(0, 0, 0, 0, 0, 0)) &&
            typ <= IntVect(D_DECL6(1, 1, 1, 1, 1, 1)));
  // computeBoxLen();
}

inline DEVICE void
Box::define(const IntVect& small,
            const IntVect& big,
            const IntVect& typ)
{
  CH_assert (small <= big);
  CH_assert(typ >= IntVect(D_DECL6(0, 0, 0, 0, 0, 0)) &&
            typ <= IntVect(D_DECL6(1, 1, 1, 1, 1, 1)));
  smallend = small;
  bigend   = big;
  btype    = typ;
  // computeBoxLen();
}

inline DEVICE
Box::Box (const IntVect&   small,
          const IntVect&   big,
          const IndexType& t)
  :
  smallend(small),
  bigend(big),
  btype(t)
{
  CH_assert (small <= big);
  // computeBoxLen();
}

inline DEVICE void
Box::define(const IntVect&   small,
            const IntVect&   big,
            const IndexType& t)
{
  CH_assert (small <= big);
  smallend = small;
  bigend   = big;
  btype    = t;
  // computeBoxLen();
}

// note that if b is undefined, then this will be undefined as well
// (but is not an error)
inline DEVICE void
Box::define(const Box& b)
{
  smallend = b.smallend;
  bigend = b.bigend;
  btype = b.btype;
}
#endif  /* CUDA */

//inline
//Box&
//Box::operator= (const Box& b)
//{
//  smallend = b.smallend;
//  bigend = b.bigend;
//  btype = b.btype;
//   D_EXPR6(len[0] = b.len[0],
//          len[1] = b.len[1],
//          len[2] = b.len[2]);
//  return *this;
//}

inline DEVICE
IntVect
Box::sideEnd(Side::LoHiSide a_side) const
{
  IntVect retval;
  if (a_side == Side::Lo)
    retval = smallEnd();
  else
    retval = bigEnd();
  return retval;
}

inline DEVICE
const IntVect&
Box::smallEnd () const
{
  return smallend;
}

inline DEVICE
int
Box::smallEnd (int dir) const
{
  return smallend[dir];
}

inline DEVICE
const IntVect&
Box::bigEnd () const
{
  return bigend;
}

inline DEVICE
int
Box::bigEnd (int dir) const
{
  return bigend[dir];
}

inline DEVICE
IndexType
Box::ixType () const
{
  return btype;
}

inline DEVICE
IntVect
Box::type () const
{
  return btype.ixType();
}

inline DEVICE
IndexType::CellIndex
Box::type (int dir) const
{
  return btype.ixType(dir);
}

inline DEVICE
IntVect
Box::size () const
{
#ifdef USE_STCVECTOR
  return IntVect(bigend - smallend + 1);
#else
  return IntVect(D_DECL6(bigend[0]-smallend[0]+1,
                         bigend[1]-smallend[1]+1,
                         bigend[2]-smallend[2]+1,
                         bigend[3]-smallend[3]+1,
                         bigend[4]-smallend[4]+1,
                         bigend[5]-smallend[5]+1));
#endif
}

namespace BoxHelper
{

//--Helper counters

  template <size_t DIM>
  struct Count
  {
    static const size_t dimsize(const IntVect& lo, const IntVect& hi)
      {
        return (hi[DIM-1]-lo[DIM-1]+1)*Count<DIM-1>::dimsize(lo, hi);
      }
  };
  // Termination
  template <>
  struct Count<1>
  {
    static const size_t dimsize(const IntVect& lo, const IntVect& hi)
      {
        return hi[0]-lo[0]+1;
      }
  };
}  // namespace BoxHelper

inline DEVICE size_t
Box::numPts () const
{
#ifdef USE_STCVECTOR
  stc::VectorImpl<size_t, SpaceDim> dims(bigend - smallend + 1);
  return stc::product(dims);
#else  
  return BoxHelper::Count<CH_SPACEDIM>::dimsize(smallend, bigend);
#endif
}

inline DEVICE size_t
Box::volume () const
{
#ifdef USE_STCVECTOR
  stc::VectorImpl<size_t, SpaceDim> dims(bigend - smallend + 1);
  return stc::product(dims);
#else
  return BoxHelper::Count<CH_SPACEDIM>::dimsize(smallend, bigend);
#endif
}

inline DEVICE int
Box::longside () const
{
#ifdef USE_STCVECTOR
  return IntVect(bigend - smallend + 1).max();
#else
  int maxlen = bigend[0]-smallend[0]+1;
  for (int i = 1; i < CH_SPACEDIM; i++)
    {
      int len=bigend[i]-smallend[i]+1;
        if (len > maxlen)
          maxlen = len;
    }
  return maxlen;
#endif
}

inline DEVICE int
Box::longside (int& dir) const
{
#ifdef USE_STCVECTOR
  IntVect dims = bigend - smallend + 1;
  dir = stc::indexMaxElem(dims);
  return dims[dir];
#else
  int maxlen = bigend[0]-smallend[0]+1;
  dir=0;
  for (int i = 1; i < CH_SPACEDIM; i++)
    {
      int len=bigend[i]-smallend[i]+1;
        if (len > maxlen)
          {
            maxlen = len;
            dir=i;
          }
    }
  return maxlen;
#endif
}

inline DEVICE int
Box::shortside () const
{
#ifdef USE_STCVECTOR
  return IntVect(bigend - smallend + 1).min();
#else
  int dir;
  return shortside(dir);
#endif
}

inline DEVICE int
Box::shortside (int& dir) const
{
#ifdef USE_STCVECTOR
  IntVect dims = bigend - smallend + 1;
  dir = stc::indexMinElem(dims);
  return dims[dir];
#else
  int minlen = bigend[0]-smallend[0]+1;
  dir=0;
  for (int i = 1; i < CH_SPACEDIM; i++)
    {
      int len=bigend[i]-smallend[i]+1;
        if (len < minlen)
          {
            minlen = len;
            dir=i;
          }
    }
  return minlen;
#endif
}

inline DEVICE
int
Box::size (int dir) const
{
  // return len[dir];
  return bigend[dir]-smallend[dir]+1;
}

inline DEVICE
const int*
Box::loVect () const
{
  return smallend.getVect();
}

inline DEVICE
const int*
Box::hiVect () const
{
  return bigend.getVect();
}

inline DEVICE
const int*
Box::getVect () const
{
  return smallend.getVect();
}

inline DEVICE
bool
Box::isEmpty () const
{
#ifdef USE_STCVECTOR
  stc::BinaryOpLessThan_E<IntVect, IntVect> f(bigend, smallend);
  return stc::orEachElement<SpaceDim>(f);
#else
  return D_TERM6(   (bigend[0] < smallend[0]),
                 || (bigend[1] < smallend[1]),
                 || (bigend[2] < smallend[2]),
                 || (bigend[3] < smallend[3]),
                 || (bigend[4] < smallend[4]),
                 || (bigend[5] < smallend[5]) );
#endif
}

inline DEVICE
bool
Box::ok () const
{
  return bigend >= smallend;
}

inline DEVICE
bool
Box::contains (const IntVect& p) const
{
  return p >= smallend && p <= bigend;
}

inline DEVICE
bool
Box::sameType (const Box &b) const
{
  return btype == b.btype;
}

inline DEVICE
bool
Box::contains (const Box& b) const
{
  CH_assert(sameType(b));
  return b.isEmpty() ||
    (!isEmpty() && b.smallend >= smallend && b.bigend <= bigend);
}

//
// Boolean functions.
//

inline DEVICE bool
Box::intersects (const Box& b) const
{
  if (b.isEmpty() || isEmpty()) return false;
  CH_assert(sameType(b));
  IntVect low(smallend);
  IntVect hi(bigend);
  low.max(b.smallend);
  hi.min(b.bigend);
  return low <= hi;
}

inline DEVICE bool
Box::intersectsNotEmpty(const Box& b) const
{
  IntVect low(smallend);
  IntVect hi(bigend);
  low.max(b.smallend);
  hi.min(b.bigend);
  return low <= hi;
}

inline DEVICE
bool
Box::sameSize (const Box& b) const
{
  CH_assert(sameType(b));
  return D_TERM6(bigend[0] -smallend[0] == b.bigend[0]-b.smallend[0],
                 && bigend[1] -smallend[1] == b.bigend[1]-b.smallend[1],
                 && bigend[2] -smallend[2] == b.bigend[2]-b.smallend[2],
                 && bigend[3] -smallend[3] == b.bigend[3]-b.smallend[3],
                 && bigend[4] -smallend[4] == b.bigend[4]-b.smallend[4],
                 && bigend[5] -smallend[5] == b.bigend[5]-b.smallend[5]);
}

inline DEVICE bool
Box::operator== (const Box& b) const
{
  return smallend == b.smallend && bigend == b.bigend && b.btype == btype;
}

inline DEVICE bool
Box::eq(const Box& b) const
{
  return *this == b;
}

inline DEVICE bool
Box::operator!= (const Box& b) const
{
  return !operator==(b);
}

inline DEVICE bool
Box::neq(const Box& b) const
{
  return *this != b;
}

inline DEVICE
bool
Box::cellCentered () const
{
  return !btype.any();
}

inline DEVICE
long
Box::index (const IntVect& v) const
{
  long result = v[0]-smallend[0];
#if CH_SPACEDIM == 1
  // nothing more goes here
#elif   CH_SPACEDIM==2
  int len0=(bigend[0]-smallend[0]+1);
  result += len0*(v[1]-smallend[1]);
#elif CH_SPACEDIM==3
  long int len0=(bigend[0]-smallend[0]+1);
  long int len1=(bigend[1]-smallend[1]+1);
  long int dif1= v[1]-smallend[1];
  long int dif2= v[2]-smallend[2];
  result += len0*(dif1
                + dif2*len1);
#elif CH_SPACEDIM==4
  int len0=(bigend[0]-smallend[0]+1);
  int len1=(bigend[1]-smallend[1]+1);
  int len2=(bigend[2]-smallend[2]+1);
  result += len0*((v[1]-smallend[1])
                  +len1*((v[2]-smallend[2])
                         +len2*(v[3]-smallend[3])));
#elif CH_SPACEDIM==5
  int len0=(bigend[0]-smallend[0]+1);
  int len1=(bigend[1]-smallend[1]+1);
  int len2=(bigend[2]-smallend[2]+1);
  int len3=(bigend[3]-smallend[3]+1);
  result += len0*(v[1]-smallend[1]
                  +len1*(v[2]-smallend[2]
                         +len2*(v[3]-smallend[3]
                                +len3*(v[4]-smallend[4]))));
#elif CH_SPACEDIM==6
  int len0=(bigend[0]-smallend[0]+1);
  int len1=(bigend[1]-smallend[1]+1);
  int len2=(bigend[2]-smallend[2]+1);
  int len3=(bigend[3]-smallend[3]+1);
  int len4=(bigend[4]-smallend[4]+1);
  result += len0*(v[1]-smallend[1]
                  +len1*(v[2]-smallend[2]
                         +len2*(v[3]-smallend[3]
                                +len3*(v[4]-smallend[4]
                                       +len4*(v[5]-smallend[5])))));

#else
  DIM > 6 undefined!
#endif
  return result;
}

inline DEVICE
void
Box::computeBoxLen ()
{
//   if (isEmpty())
//     {
//       len = IntVect::Zero;
//     }
//   else
//     {
//       D_EXPR6(len[0] = bigend[0]-smallend[0] + 1,
//              len[1] = bigend[1]-smallend[1] + 1,
//              len[2] = bigend[2]-smallend[2] + 1);
//     }
}

inline DEVICE
void
Box::computeBoxLenNotEmpty()
{
//   D_EXPR6(len[0] = bigend[0]-smallend[0] + 1,
//          len[1] = bigend[1]-smallend[1] + 1,
//          len[2] = bigend[2]-smallend[2] + 1);
}

inline DEVICE
Box&
Box::setSmall (const IntVect& sm)
{
  CH_assert (sm <= bigend);

  smallend = sm;
  computeBoxLen();
  return *this;
}

inline DEVICE
Box&
Box::setSmall (int dir,
               int sm_index)
{
  CH_assert (sm_index <= bigend[dir]);

  smallend.setVal(dir,sm_index);
  computeBoxLen();
  return *this;
}

inline DEVICE
Box&
Box::setBig (const IntVect& bg)
{
  CH_assert (bg >= smallend);

  bigend = bg;
  computeBoxLen();
  return *this;
}

inline DEVICE
Box&
Box::setBig (int dir,
             int bg_index)
{
  CH_assert (bg_index >= smallend[dir]);

  bigend.setVal(dir,bg_index);
  computeBoxLen();
  return *this;
}

inline DEVICE
Box&
Box::setRange (int dir,
               int sm_index,
               int n_cells)
{
  CH_assert (n_cells > 0);

  smallend.setVal(dir,sm_index);
  bigend.setVal(dir,sm_index+n_cells-1);
  computeBoxLen();
  return *this;
}

inline DEVICE
Box&
Box::shift (int dir,
            int nzones)
{
  if (!isEmpty())
    {
      smallend.shift(dir,nzones);
      bigend.shift(dir,nzones);
    }
  return *this;
}

inline DEVICE Box&
Box::shift (const IntVect& iv)
{
  if (!isEmpty())
    {
      smallend.shift(iv);
      bigend.shift(iv);
    }
  return *this;
}

inline DEVICE Box&
Box::shift_intvect (const IntVect& iv)
{
  return shift(iv);
}

inline DEVICE Box&
Box::shiftHalf (int dir,
                int nzones)
{
  if (!isEmpty())
  {
    int nbit = (nzones<0 ? -nzones : nzones)%2;
    int nshift = nzones/2;
    unsigned int bit_dir = btype[dir];
    //
    // Toggle btyp bit if nzones is odd.
    //
    if (nbit)
        btype.flip(dir);
    if (nzones < 0)
        nshift -= (bit_dir ? nbit : 0);
    else
        nshift += (bit_dir ? 0 : nbit);
    smallend.shift(dir,nshift);
    bigend.shift(dir,nshift);
  }
  return *this;
}

inline DEVICE Box&
Box::shiftHalf (const IntVect& nz)
{
  if (!isEmpty())
  {
   for (int dir = 0; dir < CH_SPACEDIM; dir++)
   {
       int nzones = nz[dir];
       int nbit = (nzones<0 ? -nzones : nzones)%2;
       int nshift = nzones/2;
       unsigned int bit_dir = btype[dir];
       //
       // Toggle btype bit if nzones is odd.
       //
       if (nbit)
           btype.flip(dir);
       if (nzones < 0)
           nshift -= (bit_dir ? nbit : 0);
       else
           nshift += (bit_dir ? 0 : nbit);
       smallend.shift(dir,nshift);
       bigend.shift(dir,nshift);
   }
  }
   return *this;
}

inline DEVICE Box&
Box::shiftHalf_intvect (const IntVect& iv)
{
  return shiftHalf(iv);
}

inline DEVICE Box&
Box::convert (const IntVect& typ)
{
  CH_assert(typ >= IntVect::Zero && typ <= IntVect::Unit);
  if (!isEmpty())
    {
      IntVect shft(typ - btype.ixType());
      bigend += shft;
    }
  btype = IndexType(typ);
  computeBoxLen();
  return *this;
}

inline DEVICE Box&
Box::convert (int                  dir,
              IndexType::CellIndex typ)
{
  if (!isEmpty())
  {
   unsigned int bitval = btype[dir];
   int off = typ - bitval;
   bigend.shift(dir,off);
   if (off != 0)
      computeBoxLen();
  }
   //
   // Set dir'th bit to typ.
   //
   btype.setType(dir,typ);
   return *this;
}

inline DEVICE Box&
Box::convert (IndexType t)
{
   for (int dir = 0; dir < CH_SPACEDIM; dir++)
   {
      unsigned int typ = t[dir];
      if (!isEmpty())
      {
        unsigned int bitval = btype[dir];
        int off = typ - bitval;
        bigend.shift(dir,off);
      }
      btype.setType(dir, (IndexType::CellIndex) typ);
   }
   computeBoxLen();
   return *this;
}

inline DEVICE Box&
Box::surroundingNodes ()
{
  if (!isEmpty())
  {
    for (int i = 0; i < CH_SPACEDIM; ++i)
        if ((btype[i] == 0))
            bigend.shift(i,1);
  }
    btype.setall();
    computeBoxLen();
    return *this;
}

inline DEVICE Box&
Box::surroundingNodes (int dir)
{
  if (!(btype[dir]))
    {
      if (!isEmpty())
        {
          bigend.shift(dir,1);
        }
      //
      // Set dir'th bit to 1 = IndexType::NODE.
      //
      btype.set(dir);
      computeBoxLen();
    }
  return *this;
}

inline DEVICE Box&
Box::surroundingNodes_int(int dir)
{
  return surroundingNodes(dir);
}

inline DEVICE Box&
Box::enclosedCells ()
{
  if (!isEmpty())
  {
    for (int i = 0 ; i < CH_SPACEDIM; ++i)
        if (btype[i])
            bigend.shift(i,-1);
  }
    btype.clear();
    computeBoxLen();
    return *this;
}

inline DEVICE Box&
Box::enclosedCells (int dir)
{
  if (btype[dir])
    {
      if (!isEmpty())
        {
          bigend.shift(dir,-1);
        }
      //
      // Set dir'th bit to 0 = IndexType::CELL.
      //
      btype.unset(dir);
      computeBoxLen();
    }
  return *this;
}

inline DEVICE Box&
Box::enclosedCells_int (int dir)
{
  return enclosedCells(dir);
}

///
/**
   Returns a Box with NODE based coordinates in direction dir that
   encloses the argument Box.  NOTE: equivalent to
   b.convert(dir,NODE).  NOTE: error if b.type(dir) == NODE.
   Directions are zero-based.  It is an error if not 0 <= dir <
   SpaceDim.  The Empty Box is not modified by this function except
   for changing of type.

*/
inline DEVICE
Box
surroundingNodes (const Box& b,
                  int        dir)
{
  Box bx(b);
  return bx.surroundingNodes(dir);
}

///
/**
   Returns a Box with NODE based coordinates in all directions that
   encloses argument Box.  The upper corner of the return Box will
   be one more in all components than the upper corner of the
   argument Box.  The Empty Box is not modified by this function
   except for changing of type.

*/
inline DEVICE
Box
surroundingNodes (const Box& b)
{
  Box bx(b);
  return bx.surroundingNodes();
}

///
/**
   Returns a Box with CELL based coordinates in direction dir that
   is enclosed by argument Box.  NOTE: equivalent to
   b.convert(dir,CELL).  NOTE: error if b.type(dir) == CELL.
   Directions are zero-based.  It is an error if not 0 <= dir <
   SpaceDim.  The Empty Box is not modified by this function except
   for changing of type.

*/
inline DEVICE
Box
enclosedCells (const Box& b,
               int        dir)
{
  Box bx(b);
  return bx.enclosedCells(dir);
}

///
/**
   Returns a Box with CELL based coordinates in all directions that
   is enclosed by argument Box.  The upper corner of the return Box
   will be one less in all components than the upper corner of the
   argument Box. The Empty Box is not modified by this function
   except for changing of type.

*/
inline DEVICE
Box
enclosedCells (const Box& b)
{
  Box bx(b);
  return bx.enclosedCells();
}

inline DEVICE
Box&
Box::operator+= (const IntVect& v)
{
  if (!isEmpty())
    {
      smallend += v;
      bigend += v;
    }
  return *this;
}

inline DEVICE
Box
Box::operator+  (const IntVect& v) const
{
  if (isEmpty())
    {
      return(Box().convert(btype));
    }
  else
    {
      IntVect small(smallend);
      small += v;
      IntVect big(bigend);
      big += v;
      return Box(small,big,btype);
    }
}

inline DEVICE
Box&
Box::operator-= (const IntVect& v)
{
  if (!isEmpty())
    {
      smallend -= v;
      bigend -= v;
    }
  return *this;
}

inline DEVICE bool
Box::operator < (const Box& rhs) const
{
  if (smallend == rhs.smallend) return bigend.lexLT(rhs.bigend);
  return(!isEmpty() && (rhs.isEmpty() || smallend.lexLT(rhs.smallend)));
}

inline DEVICE bool
Box::lt(const Box& rhs) const
{
  return *this < rhs;
}

inline DEVICE
Box
Box::operator- (const IntVect& v) const
{
  if (isEmpty())
    {
      return(Box().convert(btype));
    }
  else
    {
      IntVect small = smallend;
      small -= v;
      IntVect big = bigend;
      big -= v;
      return Box(small,big,btype);
    }
}

inline DEVICE Box
bdryBox (const Box& b,
         int        dir,
         Side::LoHiSide a_side,
         int        len)
{
  Box retval;
  if (a_side == Side::Lo)
    retval = bdryLo(b, dir, len);
  else
    retval = bdryHi(b, dir, len);

  return retval;
}

inline DEVICE Box
bdryLo (const Box& b,
        int        dir,
        int        len)
{
  //
  // set dir'th bit to 1 = IndexType::NODE.
  //
  IndexType typ(b.btype);
  typ.set(dir);
  if (b.isEmpty())
  {
    return(Box().convert(typ));
  }
  else
  {
    IntVect low(b.smallend);
    IntVect hi(b.bigend);
    int sm = low[dir];
    hi.setVal(dir,sm+len-1);
    return Box(low,hi,typ);
  }
}

inline DEVICE Box
bdryHi (const Box& b,
        int        dir,
        int        len)
{
  //
  // Set dir'th bit to 1 = IndexType::NODE.
  //
  IndexType typ(b.btype);
  typ.set(dir);
  if (b.isEmpty())
  {
    return(Box().convert(typ));
  }
  else
  {
    IntVect low(b.smallend);
    IntVect hi(b.bigend);
    unsigned int bitval = b.btype.ixType(dir);
    int bg              = hi[dir]  + 1 - bitval%2;
    low.setVal(dir,bg);
    hi.setVal(dir,bg+len-1);
    return Box(low,hi,typ);
  }
}

inline DEVICE Box
adjCellLo (const Box& b,
           int        dir,
           int        len)
{
  CH_assert(len != 0);
  //
  // Set dir'th bit to 0 = IndexType::CELL.
  //
  IndexType typ(b.btype);
  typ.unset(dir);
  if (b.isEmpty())
  {
    return(Box().convert(typ));
  }
  else
  {
    IntVect low(b.smallend);
    IntVect hi(b.bigend);
    int sm = low[dir];
    if (len > 0)
      {
        low.setVal(dir,sm - len);
        hi.setVal(dir,sm - 1);
      }
    else if (len < 0)
      {
        hi.setVal(dir, sm - len - 1);
      }
    return Box(low,hi,typ);
  }
}

inline DEVICE Box
adjCellHi (const Box& b,
           int        dir,
           int        len)
{
  CH_assert(len != 0);
  //
  // Set dir'th bit to 0 = IndexType::CELL.
  //
  IndexType typ(b.btype);
  typ.unset(dir);
  if (b.isEmpty())
  {
    return(Box().convert(typ));
  }
  else
  {
    IntVect low(b.smallend);
    IntVect hi(b.bigend);
    unsigned int bitval = b.btype.ixType(dir);
    int bg              = hi[dir]  + 1 - bitval%2;
    if (len > 0)
      {
        low.setVal(dir,bg);
        hi.setVal(dir,bg + len - 1);
      }
    else if (len < 0)
      {
        hi.setVal(dir,bg-1);
        low.setVal(dir, bg + len);
      }
    return Box(low,hi,typ);
  }
}

inline DEVICE Box
Box::operator& (const Box& rhs) const
{
  Box lhs(*this);
  return lhs &= rhs;
}

inline DEVICE Box&
Box::operator&= (const Box& b)
{
  CH_assert(sameType(b));
  if (isEmpty() || b.isEmpty())
  {
    *this = Box().convert(b.btype);
  }
  else
  {
    smallend.max(b.smallend);
    bigend.min(b.bigend);
    if (!(bigend >= smallend)) *this = Box().convert(b.btype);
    else computeBoxLen(); // I'd like to stop doing this....bvs
  }
  return *this;
}

inline DEVICE Box
adjCellBox (const Box& b,
            int        dir,
            Side::LoHiSide a_side,
            int        len)
{
  Box retval;
  if (a_side == Side::Lo)
    retval = adjCellLo(b, dir, len);
  else
    retval = adjCellHi(b, dir, len);

  return retval;
}

inline DEVICE Box&
Box::minBox (const Box &b)
{
  CH_assert(sameType(b));
  if (isEmpty())
  {
    *this = b;
  }
  else if (!b.isEmpty())
  {
    smallend.min(b.smallend);
    bigend.max(b.bigend);
    computeBoxLen();
  }
  return *this;
}

inline DEVICE Box
minBox (const Box& b,
        const Box& o)
{
    CH_assert(o.sameType(b));
    if (b.isEmpty())
    {
      return(o);
    }
    else if (o.isEmpty())
    {
      return(b);
    }
    else
    {
      IntVect small = b.smallend;
      IntVect big = b.bigend;
      small.min(o.smallend);
      big.max(o.bigend);
      return Box(small,big,b.btype);
    }
}

inline DEVICE Box&
Box::grow (int i)
{
  if (!isEmpty())
    {
      smallend.diagShift(-i);
      bigend.diagShift(i);
      if (!(bigend >= smallend)) *this = Box().convert(btype);
      computeBoxLen();
    }
  return *this;
}

///
/**
   Returns a Box that is the argument Box grown in all directions by
   given amount.  The Empty Box is not modified by this function.
   NOTE: n_cell negative shrinks the Box by that number of cells.
   If shrinking produces an empty Box, the Box is transformed into
   the canonical Empty Box.

*/
inline DEVICE
Box
grow (const Box& b, int i)
{
  if (b.isEmpty())
    {
      return Box().convert(b.btype);
    }

  IntVect small = diagShift(b.smallend,-i);
  IntVect big   = diagShift(b.bigend,i);
  if (!(big >= small))
    {
      return Box().convert(b.btype);
    }
  return Box(small,big,b.btype);
}

inline DEVICE
Box&
Box::grow (const IntVect& v)
{
  if (!isEmpty())
    {
      smallend -= v;
      bigend   += v;
      if (!(bigend >= smallend)) *this = Box().convert(btype);
      computeBoxLen();
    }
  return *this;
}

inline DEVICE
Box
grow (const Box&     b,
      const IntVect& v)
{
  if (b.isEmpty())
    {
      return Box().convert(b.btype);
    }
  else
    {
      IntVect small = b.smallend - v;
      IntVect big   = b.bigend   + v;
      if (!(big >= small))
        {
          return Box().convert(b.btype);
        }
      else
        {
          return Box(small,big,b.btype);
        }
    }
}

inline DEVICE
Box&
Box::grow (int idir,
           int n_cell)
{
  if (!isEmpty())
    {
      smallend.shift(idir, -n_cell);
      bigend.shift(idir, n_cell);
      if (!(bigend >= smallend)) *this = Box().convert(btype);
      computeBoxLen();
    }
  return *this;
}

inline DEVICE
Box&
Box::growLo (int idir,
             int n_cell)
{
  if (!isEmpty())
    {
      smallend.shift(idir, -n_cell);
      if (!(bigend >= smallend)) *this = Box().convert(btype);
      computeBoxLen();
    }
  return *this;
}

inline DEVICE
Box&
Box::growDir (int idir,
              const Side::LoHiSide& a_sd,
              int n_cell)
{
  if (a_sd == Side::Lo)
    {
      growLo(idir, n_cell);
    }
  else
    {
      growHi(idir, n_cell);
    }
  return *this;
}

inline DEVICE
Box&
Box::growHi (int idir,
             int n_cell)
{
  if (!isEmpty())
    {
      bigend.shift(idir,n_cell);
      if (!(bigend >= smallend)) *this = Box().convert(btype);
      computeBoxLen();
    }
  return *this;
}

inline DEVICE Box&
Box::refine (int refinement_ratio)
{
  CH_assert(refinement_ratio > 0);
  if (!isEmpty())
  {
    IntVect shft(IntVect::Unit);
    shft -= btype.ixType();
    smallend.scale(refinement_ratio);
    bigend += shft; // Bigend does more than just multiply.
    bigend.scale(refinement_ratio);
    bigend -= shft;
    computeBoxLen();
  }
  return *this;
}

inline DEVICE Box
refine (const Box& b,
        int        refinement_ratio)
{
  CH_assert(refinement_ratio > 0);
  if (b.isEmpty()) return(b);

  IntVect small(b.smallend);
  small.scale(refinement_ratio);
  IntVect big(b.bigend);
  IntVect shft(IntVect::Unit);
  shft -= b.btype.ixType();
  big += shft;        // Large end is not just multiply.
  big.scale(refinement_ratio);
  big -= shft;
  return Box(small,big,b.btype);
}

inline DEVICE Box&
Box::refine (const IntVect& refinement_ratio)
{
  if (!isEmpty())
  {
    IntVect shft(IntVect::Unit);
    shft -= btype.ixType();
    smallend *= refinement_ratio;
    bigend += shft;
    bigend   *= refinement_ratio;
    bigend -= shft;
    computeBoxLen();
  }
  return *this;
}

inline DEVICE Box
refine (const Box&     b,
        const IntVect& refinement_ratio)
{
  if (b.isEmpty()) return(b);

  IntVect small(b.smallend);
  small *= refinement_ratio;
  IntVect big(b.bigend);
  IntVect shft(IntVect::Unit);
  shft -= b.btype.ixType();
  big += shft;
  big *= refinement_ratio;
  big -= shft;
  return Box(small,big,b.btype);
}

inline DEVICE Box&
Box::coarsen (const IntVect& refinement_ratio)
{
  if (!isEmpty())
  {
    smallend.coarsen(refinement_ratio);

    if (btype.any())
    {
        IntVect off(IntVect::Zero);
        for (int dir = 0; dir < CH_SPACEDIM; dir++)
        {
            if (btype[dir])
            {
                int b = bigend[dir];
                int r = refinement_ratio[dir];
                if (b%r)
                    off.setVal(dir,1);
            }
        }
        bigend.coarsen(refinement_ratio);
        bigend += off;
    }
    else
        bigend.coarsen(refinement_ratio);

    computeBoxLen();
  }
    return *this;
}

inline DEVICE Box
coarsen (const Box& b,
         int  refinement_ratio)
{
  if (b.isEmpty()) return(b);

    IntVect small(b.smallend);
    small.coarsen(refinement_ratio);
    IntVect big(b.bigend);

    if (b.btype.any())
    {
        IntVect off(IntVect::Zero);
        for (int dir = 0; dir < CH_SPACEDIM; dir++)
        {
            if (b.btype[dir])
                if (big[dir]%refinement_ratio)
                    off.setVal(dir,1);
        }
        big.coarsen(refinement_ratio);
        big += off;
    }
    else
        big.coarsen(refinement_ratio);

    return Box(small,big,b.btype);
}

inline DEVICE Box&
Box::coarsen (int refinement_ratio)
{
  if (!isEmpty())
  {
    smallend.coarsen(refinement_ratio);
    if (btype.any())
    {
        IntVect off(IntVect::Zero);
        for (int dir = 0; dir < CH_SPACEDIM; dir++)
        {
            if (btype[dir])
                if (bigend[dir]%refinement_ratio)
                    off.setVal(dir,1);
        }
        bigend.coarsen(refinement_ratio);
        bigend += off;
    }
    else
        bigend.coarsen(refinement_ratio);

    computeBoxLen();
  }
    return *this;
}

inline DEVICE Box
coarsen (const Box&     b,
         const IntVect& refinement_ratio)
{
  if (b.isEmpty()) return(b);

    IntVect small(b.smallend);
    small.coarsen(refinement_ratio);
    IntVect big(b.bigend);

    if (b.btype.any())
    {
        IntVect off(IntVect::Zero);
        for (int dir = 0; dir < CH_SPACEDIM; dir++)
        {
            if (b.btype[dir])
                if (big[dir]%refinement_ratio[dir])
                    off.setVal(dir,1);
        }
        big.coarsen(refinement_ratio);
        big += off;
    }
    else
        big.coarsen(refinement_ratio);

    return Box(small,big,b.btype);
}

#endif  /* WRAPPER */

#include "NamespaceFooter.H"

#endif  /* ! defined _BOX_H_ */
