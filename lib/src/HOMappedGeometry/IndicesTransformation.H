#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _INDICESTRANSFORMATION_H_
#define _INDICESTRANSFORMATION_H_

#include "IntVect.H"
#include "RealVect.H"
#include "FArrayBox.H"
#include "Box.H"

#include "NamespaceHeader.H"

// ---------------------------------------------------------
/// Class to describe transformation of indices from one block to another.
class IndicesTransformation
{
public:

  /// null constructor leaves object in undefined state.
  IndicesTransformation();

  ///
  /**
     Constructor.

     IntVect pOld:  indices of a cell in old index space
     IntVect pNew:  indices of the same cell in new index space

     Then for each direction idir:
     pNew[idir] = a_sign[idir]*pOld[a_permutation[idir]] + a_translation[idir].
  */
  IndicesTransformation(const IntVect& a_permutation,
                        const IntVect& a_sign,
                        const IntVect& a_translation);

  ///
  /*
    Same as constructor.
   */
  void define(const IntVect& a_permutation,
              const IntVect& a_sign,
              const IntVect& a_translation);

  ///
  /*
    Define an IndicesTransformation
    which is just a translation.
   */
  void defineFromTranslation(const IntVect& a_translation);

  ///
  /*
    Define an IndicesTransformation
    such that the NODE a_pivotOld maps to the NODE a_pivotNew
    and for each idir, [idir] maps to [a_permutation[idir]]
    with sign a_sign[idir].
   */
  void defineFromPivot(const IntVect& a_pivotOld,
                       const IntVect& a_pivotNew,
                       const IntVect& a_permutation,
                       const IntVect& a_sign);

  ///
  /*
    Define an IndicesTransformation
    such that a_srcBox's face in dimension a_srcDim side a_srcSide
    maps to a_dstBox's face in dimension a_dstDim side a_dstSide
    and for each idir, a_sign[idir] is -1 or +1 according to whether
    source direction idir is flipped or not in the opposite block.
   */
  void defineFromFaces(const Box& a_srcBox,
                       int a_srcDim,
                       Side::LoHiSide a_srcSide,
                       const Box& a_dstBox,
                       int a_dstDim,
                       Side::LoHiSide a_dstSide,
                       const IntVect& a_sign = IntVect::Unit);

  ///
  /*
    Is this transformation defined?
   */
  bool isDefined() const;

  bool operator==(const IndicesTransformation& a_itOther) const;

  ///
  /**
     IntVect pOld:  indices of a cell in old index space
     IntVect pNew:  indices of the same cell in new index space

     Then for each direction idir:
     pNew[idir] = m_sign[idir]*pOld[m_permutation[idir]] + m_translation[idir]
     and hence:
     pOld[m_permutation[idir]] = m_sign[idir]*pNew[idir]
                               - m_sign[idir]*m_translation[idir]
  */
  IntVect transformFwd(const IntVect& a_ivOld) const;

  IntVect transformBack(const IntVect& a_ivNew) const;

  IntVect transform(const IntVect& a_iv, bool a_forward = true) const;

  ///
  /**
     Transform NODE indices.
   */
  IntVect transformNode(const IntVect& a_iv) const;

  ///
  /**
     Transform coordinates in mapped space, not an index.
   */
  RealVect transformMapped(const RealVect& a_pointOld,
                           const RealVect& a_dxOld,
                           const RealVect& a_dxNew) const;

  ///
  /**
     Transform a vector, not an index.
   */
  IntVect transformVectorFwd(const IntVect& a_vecOld) const;

  IntVect transformVectorBack(const IntVect& a_vecNew) const;

  IntVect transformVector(const IntVect& a_vec, bool a_forward = true) const;

  ///
  /**
     Transform indices of either cells or nodes or faces or edges.
   */
  IntVect transformWithType(const IntVect& a_iv,
                            const IntVect& a_tp,
                            bool a_forward = true) const;

  ///
  /**
     Transform the type of Box:
     no change if a_tp is the cell type or the node type,
     but if a face or edge type, then the transformed type may differ.
   */
  IntVect transformType(const IntVect& a_tp, bool a_forward = true) const;

  ///
  /**
     Transform a whole box.
   */
  Box transformFwd(const Box& a_bxOld) const;

  Box transformBack(const Box& a_bxNew) const;

  Box transform(const Box& a_bx, bool a_forward = true) const;

  ///
  /**
     Transform a FAB.
   */
  void transformFwd(FArrayBox& a_transformedFab,
                    const FArrayBox& a_originalFab,
                    const Box& a_origBox,
                    const Interval& a_transIntvl,
                    const Interval& a_origIntvl) const;

  void transformBack(FArrayBox& a_transformedFab,
                     const FArrayBox& a_originalFab,
                     const Box& a_origBox,
                     const Interval& a_transIntvl,
                     const Interval& a_origIntvl) const;

  void transform(FArrayBox& a_transformedFab,
                 const FArrayBox& a_originalFab,
                 const Box& a_origBox,
                 const Interval& a_transIntvl,
                 const Interval& a_origIntvl,
                 bool a_forward = true) const;

  ///
  /**
     Return the inverse transformation.
   */
  IndicesTransformation inverse() const;

  ///
  /**
     Return the composite transformation:
     that is, this transformation followed by a_next.
   */
  IndicesTransformation compose(const IndicesTransformation& a_next) const;

  ///
  /**
     Return this transformation with index spaces coarsened by a_refRatio.
   */
  IndicesTransformation coarsen(int a_refRatio) const;

  ///
  /**
     Return this transformation with index spaces refined by a_refRatio.
   */
  IndicesTransformation refine(int a_refRatio) const;

  ///
  /** The identity transformation, which has no permutation,
      no sign change, and no translation.
   */
  static const IndicesTransformation Identity;

  ///
  /** The undefined transformation.
   */
  static const IndicesTransformation Undefined;

  IntVect getPermutation() const
  {
    return m_permutation;
  }
  IntVect getSign() const
  {
    return m_sign;
  }
  IntVect getTranslation() const
  {
    return m_translation;
  }

  ///
  /** Initializes IndicesTransformation::Identity .
   */
  static int InitStatics();

protected:

  // Zero if undefined
  IntVect m_permutation;

  // Zero if undefined
  IntVect m_sign;

  // Zero if undefined
  IntVect m_translation;

  // Return inverse permutation of a_permutation.
  IntVect inversePermutation(const IntVect& a_permutation);

  // Make sure the defined transformation is valid
  void checkValid();
};

//
//// Static initialization.  Gotta make sure there are no static object
//// definitions above here (except possibly stuff in headers).  Otherwise,
//// the danger is that some static object's constructor calls
//// IndicesTransformation::Identity, the very things the following
//// definition is supposed to initialize.
//// (I got this from IntVect.H)
static int s_dummyForIndicesTransformationH = IndicesTransformation::InitStatics();

#include "NamespaceFooter.H"

#endif // include guard
