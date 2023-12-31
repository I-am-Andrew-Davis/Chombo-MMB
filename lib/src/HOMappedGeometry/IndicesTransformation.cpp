#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "IndicesTransformation.H"
#include "BoxIterator.H"
// #include "UGIO.H"
// #include "CH_HDF5.H"

#include "NamespaceHeader.H"

// ---------------------------------------------------------
IndicesTransformation::IndicesTransformation()
{
  m_permutation = IntVect::Zero;
  m_sign = IntVect::Zero;
  m_translation = IntVect::Zero;
}


// ---------------------------------------------------------
IndicesTransformation::IndicesTransformation(const IntVect& a_permutation,
                                             const IntVect& a_sign,
                                             const IntVect& a_translation)
{
  define(a_permutation, a_sign, a_translation);
}


// ---------------------------------------------------------
void
IndicesTransformation::define(const IntVect& a_permutation,
                              const IntVect& a_sign,
                              const IntVect& a_translation)
{
  m_permutation = a_permutation;
  m_sign = a_sign;
  m_translation = a_translation;
  checkValid();
}


// ---------------------------------------------------------
void
IndicesTransformation::defineFromTranslation(const IntVect& a_translation)
{
  m_permutation = IntVect(D_DECL6(0, 1, 2, 3, 4, 5));
  m_sign = IntVect::Unit;
  m_translation = a_translation;
  checkValid();
}


// ---------------------------------------------------------
void
IndicesTransformation::defineFromPivot(const IntVect& a_pivotOld,
                                       const IntVect& a_pivotNew,
                                       const IntVect& a_permutation,
                                       const IntVect& a_sign)
{
  m_permutation = inversePermutation(a_permutation);
  m_sign = IntVect::Zero;
  m_translation = IntVect::Zero;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      int idirPerm = m_permutation[idir];
      m_sign[idir] = a_sign[idirPerm];
      m_translation[idir] =
        a_pivotNew[idir] - m_sign[idir] * a_pivotOld[idirPerm];
      if (m_sign[idir] == -1) // because we're transforming CELLs
        m_translation[idir]--;
    }
  checkValid();
}


// ---------------------------------------------------------
void
IndicesTransformation::defineFromFaces(const Box& a_srcBox,
                                       int a_srcDim,
                                       Side::LoHiSide a_srcSide,
                                       const Box& a_dstBox,
                                       int a_dstDim,
                                       Side::LoHiSide a_dstSide,
                                       const IntVect& a_sign)
{
  Box srcFaceNodes = surroundingNodes(a_srcBox);
  const IntVect& srcFaceEnd = srcFaceNodes.sideEnd(a_srcSide);
  srcFaceNodes.setRange(a_srcDim, srcFaceEnd[a_srcDim]);

  Box dstFaceNodes = surroundingNodes(a_dstBox);
  const IntVect& dstFaceEnd = dstFaceNodes.sideEnd(a_dstSide);
  dstFaceNodes.setRange(a_dstDim, dstFaceEnd[a_dstDim]);

  IntVect perm = IntVect(D_DECL6(0, 1, 2, 3, 4, 5));
  if (a_srcDim != a_dstDim) // if srcDim == a_dstDim then perm does not change
    {
      perm[a_srcDim] = a_dstDim;
      perm[a_dstDim] = a_srcDim;
    }

  IntVect srcFaceCorner = srcFaceNodes.smallEnd();

  // dstFaceCorner is point where srcFaceCorner ends up on dest face
  IntVect dstFaceCorner = dstFaceNodes.smallEnd();
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (a_sign[perm[idir]] == -1)
        dstFaceCorner[idir] = dstFaceNodes.bigEnd(idir);
    }

  defineFromPivot(srcFaceCorner, dstFaceCorner, perm, a_sign);
  checkValid();
}


// ---------------------------------------------------------
bool
IndicesTransformation::isDefined() const
{
  return (m_sign != IntVect::Zero);
}


// ---------------------------------------------------------
bool
IndicesTransformation::operator==(const IndicesTransformation& a_itOther) const
{
  CH_assert(isDefined());
  CH_assert(a_itOther.isDefined());
  return ((m_permutation == a_itOther.m_permutation) &&
          (m_sign == a_itOther.m_sign) &&
          (m_translation == a_itOther.m_translation));
}


// ---------------------------------------------------------
IntVect
IndicesTransformation::transformFwd(const IntVect& a_ivOld) const
{
  CH_assert(isDefined());
  IntVect ivNew;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      ivNew[idir] =
        m_sign[idir]*a_ivOld[m_permutation[idir]] + m_translation[idir];
    }
  return ivNew;
}


// ---------------------------------------------------------
IntVect
IndicesTransformation::transformBack(const IntVect& a_ivNew) const
{
  CH_assert(isDefined());
  IntVect ivOld;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      ivOld[m_permutation[idir]] = m_sign[idir] *
        (a_ivNew[idir] - m_translation[idir]);
    }
  return ivOld;
}


// ---------------------------------------------------------
IntVect
IndicesTransformation::transform(const IntVect& a_iv,
                                 bool a_forward) const
{
  if (a_forward)
    return transformFwd(a_iv);
  else
    return transformBack(a_iv);
}


// ---------------------------------------------------------
IntVect
IndicesTransformation::transformNode(const IntVect& a_iv) const
{
  CH_assert(isDefined());
  // first transform as cell indices
  IntVect ivNew = transformFwd(a_iv);
  // then correct for reflections
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (m_sign[idir] == -1) ivNew[idir]++;
    }
  return ivNew;
}


// ---------------------------------------------------------
IntVect
IndicesTransformation::transformWithType(const IntVect& a_iv,
                                         const IntVect& a_tp,
                                         bool a_forward) const
{
  CH_assert(isDefined());
  // first transform as cell indices
  IntVect ivNew = transform(a_iv, a_forward);
  // then correct for reflections
  IndicesTransformation it = (a_forward) ? *this : inverse();
  IntVect tpNew = it.transformType(a_tp);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      int dirTransformed = it.m_permutation[idir];
      if (tpNew[dirTransformed] == IndexType::NODE)
        {
          if (it.m_sign[dirTransformed] == -1) ivNew[dirTransformed]++;
        }
    }
  return ivNew;
}


// ---------------------------------------------------------
IntVect
IndicesTransformation::transformType(const IntVect& a_tp,
                                     bool a_forward) const
{
  CH_assert(isDefined());
  IntVect tpNew;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      int dirTransformed = m_permutation[idir];
      if (a_forward)
        tpNew[idir] = a_tp[dirTransformed];
      else
        tpNew[dirTransformed] = a_tp[idir];
    }
  return tpNew;
}


// ---------------------------------------------------------
RealVect
IndicesTransformation::transformMapped(const RealVect& a_pointOld,
                                       const RealVect& a_dxOld,
                                       const RealVect& a_dxNew) const
{
  RealVect pointNew;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      int translationDir = m_translation[idir];
      if (m_sign[idir] == -1) translationDir++;
      int permDir = m_permutation[idir];
      pointNew[idir] = a_dxNew[idir] *
        (m_sign[idir] * a_pointOld[permDir] / a_dxOld[permDir] +
         translationDir);
    }
  return pointNew;
}


// ---------------------------------------------------------
IntVect
IndicesTransformation::transformVectorFwd(const IntVect& a_vecOld) const
{
  CH_assert(isDefined());
  // Index transformation is
  // pNew[idir] = m_sign[idir]*pOld[m_permutation[idir]] + m_translation[idir]
  // so in direction idir, the transformation of pOld + vec is
  // m_sign[idir]*(pOld[m_permutation[idir]] + vec[m_permutation[idir]]) + m_translation[idir]
  // and hence the difference in direction idir is
  // m_sign[idir] * vec[m_permutation[idir]].
  IntVect vecNew;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      vecNew[idir] = m_sign[idir] * a_vecOld[m_permutation[idir]];
    }
  return vecNew;
}


// ---------------------------------------------------------
IntVect
IndicesTransformation::transformVectorBack(const IntVect& a_vecNew) const
{
  // Index transformation is
  // pOld[m_permutation[idir]] = m_sign[idir]*pNew[idir]
  //                           - m_sign[idir]*m_translation[idir]
  // so in direction m_permutation[idir], the back-transformation of pNew + vec is
  // m_sign[idir]*(pNew[idir] + vec[idir]) - m_sign[idir]*m_translation[idir]
  // and hence the difference in direction m_permutation[idir] is
  // m_sign[idir]*vec[idir].
  CH_assert(isDefined());
  IntVect vecOld;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      vecOld[m_permutation[idir]] = m_sign[idir] * a_vecNew[idir];
    }
  return vecOld;
}


// ---------------------------------------------------------
IntVect
IndicesTransformation::transformVector(const IntVect& a_vec,
                                       bool a_forward) const
{
  if (a_forward)
    return transformVectorFwd(a_vec);
  else
    return transformVectorBack(a_vec);
}


// ---------------------------------------------------------
Box
IndicesTransformation::transformFwd(const Box& a_bx) const
{
  return transform(a_bx, true);
}


// ---------------------------------------------------------
Box
IndicesTransformation::transformBack(const Box& a_bx) const
{
  return transform(a_bx, false);
}


// ---------------------------------------------------------
Box
IndicesTransformation::transform(const Box& a_bx,
                                 bool a_forward) const
{
  if (a_bx.isEmpty())
    {
      return Box();
    }
  IntVect tp = a_bx.type();
  IntVect tpTransformed = transformType(tp, a_forward);

  IntVect bxTransformedLo =
    transformWithType(a_bx.smallEnd(), tp, a_forward);
  IntVect bxTransformedHi =
    transformWithType(a_bx.bigEnd(), tp, a_forward);

  IntVect bxTransformedMin = min(bxTransformedLo, bxTransformedHi);
  IntVect bxTransformedMax = max(bxTransformedLo, bxTransformedHi);

  Box bxTransformed = Box(bxTransformedMin, bxTransformedMax, tpTransformed);
  return bxTransformed;
}


// ---------------------------------------------------------
void
IndicesTransformation::transformFwd(FArrayBox& a_transformedFab,
                                    const FArrayBox& a_originalFab,
                                    const Box& a_origBox,
                                    const Interval& a_transIntvl,
                                    const Interval& a_origIntvl) const
{
  transform(a_transformedFab, a_originalFab, a_origBox,
            a_transIntvl, a_origIntvl, true);
}


// ---------------------------------------------------------
void
IndicesTransformation::transformBack(FArrayBox& a_transformedFab,
                                     const FArrayBox& a_originalFab,
                                     const Box& a_origBox,
                                     const Interval& a_transIntvl,
                                     const Interval& a_origIntvl) const
{
  transform(a_transformedFab, a_originalFab, a_origBox,
            a_transIntvl, a_origIntvl, false);
}


// ---------------------------------------------------------
void
IndicesTransformation::transform(FArrayBox& a_transformedFab,
                                 const FArrayBox& a_originalFab,
                                 const Box& a_origBox,
                                 const Interval& a_transIntvl,
                                 const Interval& a_origIntvl,
                                 bool a_forward) const
{
  CH_assert(a_originalFab.box().contains(a_origBox));
  Box transformedBox = transform(a_origBox, a_forward);
  CH_assert(a_transformedFab.box().contains(transformedBox));
  IntVect tp = a_origBox.type();
  int transComp = a_transIntvl.begin();
  int origComp = a_origIntvl.begin();
  int origCompHi = a_origIntvl.end();
  for (; origComp <= origCompHi; ++origComp, ++transComp)
    {
      for (BoxIterator bit(a_origBox); bit.ok(); ++bit)
        {
          IntVect ivOrig = bit();
          IntVect ivTrans = transformWithType(ivOrig, tp, a_forward);
          a_transformedFab(ivTrans, transComp) =
            a_originalFab(ivOrig, origComp);
        }
    }
}


// ---------------------------------------------------------
IndicesTransformation
IndicesTransformation::inverse() const
{
  CH_assert(isDefined());
  IntVect permutationInv, signInv, translationInv;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      int dirPermuted = m_permutation[idir];
      permutationInv[dirPermuted] = idir;
      signInv[dirPermuted] = m_sign[idir];
      translationInv[dirPermuted] = -m_sign[idir] * m_translation[idir];
    }
  return IndicesTransformation(permutationInv, signInv, translationInv);
}


// ---------------------------------------------------------
IndicesTransformation
IndicesTransformation::compose(const IndicesTransformation& a_next) const
{
  CH_assert(isDefined());
  CH_assert(a_next.isDefined());
  IntVect permutationComp, signComp, translationComp;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      int dirPermutedNext = a_next.m_permutation[idir];
      permutationComp[idir] = m_permutation[dirPermutedNext];
      signComp[idir] = a_next.m_sign[idir] * m_sign[dirPermutedNext];
      translationComp[idir] =
        a_next.m_sign[idir] * m_translation[dirPermutedNext] +
        a_next.m_translation[idir];
    }
  return IndicesTransformation(permutationComp, signComp, translationComp);
}


// ---------------------------------------------------------
IndicesTransformation
IndicesTransformation::coarsen(int a_refRatio) const
{
  CH_assert(isDefined());
  IntVect translationCoarsened;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (m_sign[idir] == +1)
        translationCoarsened[idir] = m_translation[idir] / a_refRatio;
      else // m_sign[idir] == -1
        // Think of it this way:
        // the transformation
        // interval [i, i+1] => interval [-i-1 + k, -i + k]
        // maps cell i to cell -i-1 + k,
        // so it has sign -1 and translation k-1.
        // Coarsening this transformation:
        // interval [i/r, i/r+1] => interval [-i/r-1 + k/r, -i/r + k/r]
        // maps cell i/r to cell -i/r-1 + k/r
        // so it has sign -1 and translation k/r-1.
        // Hence, in general, the refined translation is the original
        // translation with 1 added (to get k) and the result divided
        // by r and then subtracted by 1.
        translationCoarsened[idir] = (m_translation[idir] + 1) / a_refRatio - 1;
    }
  return IndicesTransformation(m_permutation, m_sign, translationCoarsened);
}


// ---------------------------------------------------------
IndicesTransformation
IndicesTransformation::refine(int a_refRatio) const
{
  CH_assert(isDefined());
  IntVect translationRefined;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (m_sign[idir] == +1)
        translationRefined[idir] = a_refRatio * m_translation[idir];
      else // m_sign[idir] == -1
        // Think of it this way:
        // the transformation
        // interval [i, i+1] => interval [-i-1 + k, -i + k]
        // maps cell i to cell -i-1 + k,
        // so it has sign -1 and translation k-1.
        // Refining this transformation:
        // interval [r*i, r*i+1] => interval [-r*i-1 + r*k, -r*i + r*k]
        // maps cell r*i to cell -r*i-1 + r*k
        // so it has sign -1 and translation r*k-1.
        // Hence, in general, the refined translation is the original
        // translation with 1 added (to get k) and the result
        // multiplied by r and then subtracted by 1.
        translationRefined[idir] = a_refRatio * (m_translation[idir] + 1) - 1;
    }
  return IndicesTransformation(m_permutation, m_sign, translationRefined);
}


// ---------------------------------------------------------
IntVect
IndicesTransformation::inversePermutation(const IntVect& a_permutation)
{
  IntVect permInverse;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      permInverse[a_permutation[idir]] = idir;
    }
  return permInverse;
}

void
IndicesTransformation::checkValid()
{
  // Are the values in valid bounds
  for (const auto dir : EachDir)
    {
      CH_assert((m_permutation[dir] >= 0) && (m_permutation[dir] < SpaceDim));
      CH_assert((m_sign[dir] == 1) || (m_sign[dir] == -1));
    }
  // Is the permutation a valid one. Must contain 0 to SpaceDim in any order
  for (const auto p : EachDir)
    {
      bool p_in_permutation = false;
      for (const auto dir : EachDir)
        {
          p_in_permutation = (p_in_permutation || (m_permutation[dir] == p));
        }
      CH_assert(p_in_permutation);
    }
}

// ---------------------------------------------------------
int
IndicesTransformation::InitStatics()
{
  // taken from IntVect.cpp
  IndicesTransformation* ptrID =
    const_cast<IndicesTransformation*>( &IndicesTransformation::Identity );
  IntVect perm(D_DECL6(0, 1, 2, 3, 4, 5));
  *ptrID = IndicesTransformation(perm, IntVect::Unit, IntVect::Zero);

  IndicesTransformation* ptrUD =
    const_cast<IndicesTransformation*>( &IndicesTransformation::Undefined );
  // *ptrUD = IndicesTransformation(IntVect::Zero, IntVect::Zero, IntVect::Zero);
  *ptrUD = IndicesTransformation();

  return 0; // arbitrary
}

const IndicesTransformation IndicesTransformation::Identity;

const IndicesTransformation IndicesTransformation::Undefined;

static int s_dummyForIndicesTransformationCpp( IndicesTransformation::InitStatics() );

#include "NamespaceFooter.H"
