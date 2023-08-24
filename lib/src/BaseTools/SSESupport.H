#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _SSESUPPORT_H_
#define _SSESUPPORT_H_


/******************************************************************************/
/**
 * \file
 *
 * \brief Support routines for SSE/AVX operations
 *
 *//*+*************************************************************************/

#include <stdint.h>

#include "SSETypes.H"

// Requires at least SSE2
#if (CHDEF_SYSTEM_X86VECEXT_COMPILER_BITS & CHDEF_BIT_SSE2)

/*--------------------------------------------------------------------*
 * Constants for SSE operations
 *--------------------------------------------------------------------*/

const __mvi minusone_vi = _mm_i32(set1)(-1);
const __mvi      one_vi = _mm_i32(set1)( 1);

const __mvr  minusone_vr = _mm_vr(set1)(-1.);
const __mvr       one_vr = _mm_vr(set1)( 1.);

// Flags with all bits set or no bits set
const uint32_t i32true  = uint32_t(0) - uint32_t(1);
const uint32_t i32false = uint32_t(0);

const uint64_t i64true  = uint64_t(0) - uint64_t(1);
const uint64_t i64false = uint64_t(0);

#ifdef CH_USE_FLOAT
#define itrue i32true
#define ifalse i32false
#else
#define itrue i64true
#define ifalse i64false
#endif

// For accumulation of test conditional into bits of an int, this integer
// marks all bits as true
const int allVecTrue = (1 << (VecSz_r + 1)) - 1;

/*--------------------------------------------------------------------*/
/// mult_vr_si
/** Multiplication of a FP vector x by an integer scalar -1 <= y <= 4
 *//*-----------------------------------------------------------------*/

template <int y>
inline __mvr mult_vr_si(const __mvr& a_x);
template<> inline __mvr mult_vr_si<-1>(const __mvr& a_x)
{ return minusone_vr*a_x; }
template<> inline __mvr mult_vr_si< 0>(const __mvr& a_x)
{ return _mm_vr(setzero)(); }
template<> inline __mvr mult_vr_si< 1>(const __mvr& a_x)
{ return a_x; }
template<> inline __mvr mult_vr_si< 2>(const __mvr& a_x)
{ return a_x + a_x; }
template<> inline __mvr mult_vr_si< 3>(const __mvr& a_x)
{ return a_x + a_x + a_x; }
template<> inline __mvr mult_vr_si< 4>(const __mvr& a_x)
{
  const __mvr tmp = a_x + a_x;
  return tmp + tmp;
}

/*--------------------------------------------------------------------*/
/// pow_vr_si
/** Raise a FP vector x to an integer power -1 <= y <= 5
 *//*-----------------------------------------------------------------*/

template <int y>
inline __mvr pow_vr_si(const __mvr& a_x);
template<> inline __mvr pow_vr_si<-1>(const __mvr& a_x)
{ return one_vr/a_x; }
template<> inline __mvr pow_vr_si< 0>(const __mvr& a_x)
{ return one_vr; }
template<> inline __mvr pow_vr_si< 1>(const __mvr& a_x)
{ return a_x; }
template<> inline __mvr pow_vr_si< 2>(const __mvr& a_x)
{ return a_x*a_x; }
template<> inline __mvr pow_vr_si< 3>(const __mvr& a_x)
{ return a_x*a_x*a_x; }
template<> inline __mvr pow_vr_si< 4>(const __mvr& a_x)
{
  const __mvr tmp = a_x*a_x;
  return tmp*tmp;
}
template<> inline __mvr pow_vr_si< 5>(const __mvr& a_x)
{
  const __mvr tmp = a_x*a_x;
  return tmp*tmp*a_x;
}

inline __mvr pow_vr_si(const __mvr& a_x, const int a_y)
{
  CH_assert(a_y >= -1 && a_y <= 5);
  switch (a_y)
    {
    case -1:
      return pow_vr_si<-1>(a_x);
    case 0:
      return pow_vr_si< 0>(a_x);
    case 1:
      return pow_vr_si< 1>(a_x);
    case 2:
      return pow_vr_si< 2>(a_x);
    case 3:
      return pow_vr_si< 3>(a_x);
    case 4:
      return pow_vr_si< 4>(a_x);
    case 5:
      return pow_vr_si< 5>(a_x);
    }
  return _mm_vr(setzero)();  // Avoid compiler warning
}

/*--------------------------------------------------------------------*/
/// abs_vr
/** Absolute value of floats
 *  We only support output for 128 bit or 256 bit, not both due to
 *  name mangling issues.  This can be fixed with certain compiler
 *  flags but I'm not sure it's worth it for I/O.
 *//*-----------------------------------------------------------------*/

//**FIXME see if we can use broadcast
#if (CHDEF_SYSTEM_X86VECEXT_COMPILER_BITS & CHDEF_BIT_AVX)

// 256 bit sf
inline __m256 abs_vr(const __m256& a_x)
{
  union mask256_t
  {
    int32_t i[8];
    float f[8];
  };
  const mask256_t mask __attribute__((aligned (32))) =
    {{ 0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF,
       0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF }};
  return _mm256_and_ps(a_x, _mm256_load_ps(mask.f));
//   const int32_t mask = 0x7FFFFFFF;
//   return _mm256_and_ps(
//     a_x, _mm256_broadcast_ss(reinterpret_cast<const float*>(&mask)));
}

// 256 bit df
inline __m256d abs_vr(const __m256d& a_x)
{
  union mask256d_t
  {
    int64_t i[4];
    double f[4];
  };
  const mask256d_t mask __attribute__((aligned (32))) =
    {{ 0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF,
       0x7FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF }};
  return _mm256_and_pd(a_x, _mm256_load_pd(mask.f));
//   const int64_t mask = 0x7FFFFFFFFFFFFFFF;
//   return _mm256_and_pd(
//     a_x, _mm256_broadcast_sd(reinterpret_cast<const double*>(&mask)));
}

#else

// 128 bit sf
inline __m128 abs_vr(const __m128& a_x)
{
  union mask128_t
  {
    int32_t i[4];
    float   f[4];
  };
  const mask128_t mask __attribute__((aligned (16))) =
    {{ 0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF, 0x7FFFFFFF }};
  return _mm_and_ps(a_x, _mm_load_ps(mask.f));
}

// 128 bit df
inline __m128d abs_vr(const __m128d& a_x)
{
  union mask128d_t
  {
    int64_t i[2];
    double  f[2];
  };
  const mask128d_t mask __attribute__((aligned (16))) =
    {{ 0x7FFFFFFFFFFFFFFFL, 0x7FFFFFFFFFFFFFFFL }};
  return _mm_and_pd(a_x, _mm_load_pd(mask.f));
}

#endif

/*--------------------------------------------------------------------*/
/**  \name Floating point comparisons.  If using AVX, these are
 *   non-signaling and return false if either argument is NAN.
 *//*-----------------------------------------------------------------*/
//@{

/// (x == y)
inline __mvr cmpeq_vr(const __mvr& a_x, const __mvr& a_y)
{
#if (CHDEF_SYSTEM_X86VECEXT_COMPILER_BITS & CHDEF_BIT_AVX)
  return _mm_vr(cmp)(a_x, a_y, _CMP_EQ_OQ);
#else
  return _mm_vr(cmpeq)(a_x, a_y);
#endif
}

/// (x != y)
inline __mvr cmpneq_vr(const __mvr& a_x, const __mvr& a_y)
{
#if (CHDEF_SYSTEM_X86VECEXT_COMPILER_BITS & CHDEF_BIT_AVX)
  return _mm_vr(cmp)(a_x, a_y, _CMP_NEQ_OQ);
#else
  return _mm_vr(cmpneq)(a_x, a_y);
#endif
}

/// (x > y)
inline __mvr cmpgt_vr(const __mvr& a_x, const __mvr& a_y)
{
#if (CHDEF_SYSTEM_X86VECEXT_COMPILER_BITS & CHDEF_BIT_AVX)
  return _mm_vr(cmp)(a_x, a_y, _CMP_GT_OQ);
#else
  return _mm_vr(cmpgt)(a_x, a_y);
#endif
}

/// (x < y)
inline __mvr cmplt_vr(const __mvr& a_x, const __mvr& a_y)
{
#if (CHDEF_SYSTEM_X86VECEXT_COMPILER_BITS & CHDEF_BIT_AVX)
  return _mm_vr(cmp)(a_x, a_y, _CMP_LT_OQ);
#else
  return _mm_vr(cmplt)(a_x, a_y);
#endif
}

/// (x >= y)
inline __mvr cmpge_vr(const __mvr& a_x, const __mvr& a_y)
{
#if (CHDEF_SYSTEM_X86VECEXT_COMPILER_BITS & CHDEF_BIT_AVX)
  return _mm_vr(cmp)(a_x, a_y, _CMP_GE_OQ);
#else
  return _mm_vr(cmpge)(a_x, a_y);
#endif
}

/// (x <= y)
inline __mvr cmple_vr(const __mvr& a_x, const __mvr& a_y)
{
#if (CHDEF_SYSTEM_X86VECEXT_COMPILER_BITS & CHDEF_BIT_AVX)
  return _mm_vr(cmp)(a_x, a_y, _CMP_LE_OQ);
#else
  return _mm_vr(cmple)(a_x, a_y);
#endif
}
//@}

/*--------------------------------------------------------------------*/
/// Swap upper 128 to lower for 'a' and lower 128 to upper for 'b'
/** \param[in]  a_xup2lo
 *                      a = (a7, a6, a5, a4, a3, a2, a1, a0)
 *  \param[in]  a_ylo2up
 *                      b = (b7, b6, b5, b4, b3, b2, b1, b0)
 *  \return             (b3, b2, b1, b0, a7, a6, a5, a4)
 *//*-----------------------------------------------------------------*/

inline __mvr swapuplo_vr(const __mvr& a_xup2lo, const __mvr& a_ylo2up)
{
#if (CHDEF_SYSTEM_X86VECEXT_COMPILER_BITS & CHDEF_BIT_AVX)
  return _mm_vr(permute2f128)(a_xup2lo, a_ylo2up, _MM_SHUFFLE(0, 2, 0, 1));
#else
#if (VecSz_r == 4)
  return _mm_vr(shuffle)(a_xup2lo, a_ylo2up, _MM_SHUFFLE(1, 0, 3, 2));
#else
  return _mm_vr(shuffle)(a_xup2lo, a_ylo2up, 1);
#endif
#endif
}

/*--------------------------------------------------------------------*/
/// Circular shift 128 bit sequences in direction 'lo'
/** \tparam     S       Amount to shift by
 *  \param[in]  a_x     a = (a7, a6, a5, a4, a3, a2, a1, a0)
 *  \return             S = 0: (a7, a6, a5, a4, a3, a2, a1, a0)
 *                          1: (a4, a7, a6, a5, a0, a3, a2, a1)
 *                          2: (a5, a4, a7, a6, a1, a0, a3, a2)
 *                          3: (a6, a5, a4, a7, a2, a1, a0, a3)
 *  \note
 *  <ul>
 *    <li> If vr = pd, shift by 64 bits (0 <= S <= 1)
 *    <li> If S is out of range, (S % VecSz_r) is used for shifting
 *  </ul>
 *//*-----------------------------------------------------------------*/

template <int S>
inline __mvr shiftlo_vr(const __mvr& a_x)
{
#if (CHDEF_SYSTEM_X86VECEXT_COMPILER_BITS & CHDEF_BIT_AVX)
  const int ncf = VecSz_r/2;  // Number of control fields
  const int fw  = VecSz_r/4;  // Field width in bits
#else
  const int ncf = VecSz_r;    // Number of control fields
  const int fw  = VecSz_r/2;  // Field width in bits
#endif
  const int control = 
    (3 % ncf) << fw*((7 - (S % 4)) % 4) |
    (2 % ncf) << fw*((6 - (S % 4)) % 4) |
    (1 % ncf) << fw*((5 - (S % 4)) % 4) |
    (0 % ncf) << fw*((4 - (S % 4)) % 4);
#if (CHDEF_SYSTEM_X86VECEXT_COMPILER_BITS & CHDEF_BIT_AVX)
  return _mm_vr(permute)(a_x, control);
#else
  return _mm_vr(shuffle)(a_x, a_x, control);
#endif
}

/*--------------------------------------------------------------------*/
/// Circular shift 128 bit sequences in direction 'up'
/** \tparam     S       Amount to shift by
 *  \param[in]  a_x     a = (a7, a6, a5, a4, a3, a2, a1, a0)
 *  \return             S = 0: (a7, a6, a5, a4, a3, a2, a1, a0)
 *                          1: (a6, a5, a4, a7, a2, a1, a0, a3)
 *                          2: (a5, a4, a7, a6, a1, a0, a3, a2)
 *                          3: (a4, a7, a6, a5, a0, a3, a2, a1)
 *  \note
 *  <ul>
 *    <li> If vr = pd, shift by 64 bits (0 <= S <= 1)
 *    <li> If S is out of range, (S % VecSz_r) is used for shifting
 *  </ul>
 *//*-----------------------------------------------------------------*/

template <int S>
inline __mvr shiftup_vr(const __mvr& a_x)
{
#if (CHDEF_SYSTEM_X86VECEXT_COMPILER_BITS & CHDEF_BIT_AVX)
  const int ncf = VecSz_r/2;  // Number of control fields
  const int fw  = VecSz_r/4;  // Field width in bits
#else
  const int ncf = VecSz_r;    // Number of control fields
  const int fw  = VecSz_r/2;  // Field width in bits
#endif
  const int control = 
    (3 % ncf) << fw*((3 + S) % 4) |
    (2 % ncf) << fw*((2 + S) % 4) |
    (1 % ncf) << fw*((1 + S) % 4) |
    (0 % ncf) << fw*((0 + S) % 4);
#if (CHDEF_SYSTEM_X86VECEXT_COMPILER_BITS & CHDEF_BIT_AVX)
  return _mm_vr(permute)(a_x, control);
#else
  return _mm_vr(shuffle)(a_x, a_x, control);
#endif
}

/*--------------------------------------------------------------------*/
/// Test if an address is aligned to 128 for SSE and 256 for AVX
/*  \param[in]  a_p     Memory address
 *  \return             T - Aligned
 *//*-----------------------------------------------------------------*/

inline bool isVecAligned(void* a_p)
{
  return (!(reinterpret_cast<size_t>(a_p) & (CH_VECLS_ALIGN - 1)));
}

/*--------------------------------------------------------------------*/
/// Return the alignment of an address.
/*  \param[in]  a_p     Memory address
 *  \return             Largest power of 2 alignment boundary
 *//*-----------------------------------------------------------------*/

//**FIXME change size_t to uintptr_t once c++11 is commonplace
size_t alignofaddr(void* a_p)
{
  size_t x = reinterpret_cast<size_t>(a_p);
  return x & ~(x-1);
}

/*--------------------------------------------------------------------*/
/// Printing of vectors
/** We only support output for 128 bit or 256 bit, not both due to
 *  name mangling issues.  This can be fixed with certain compiler
 *  flags but I'm not sure it's worth it for I/O.
 *//*-----------------------------------------------------------------*/

#if (CHDEF_SYSTEM_X86VECEXT_COMPILER_BITS & CHDEF_BIT_AVX)

//--256i (prints signed values)

inline  std::ostream& operator<<(std::ostream& a_os, const CH256i_t& a_vec)
{
  a_os << '('
     << a_vec.i[0] << ", " << a_vec.i[1] << ", "
     << a_vec.i[2] << ", " << a_vec.i[3] << ", "
     << a_vec.i[4] << ", " << a_vec.i[5] << ", "
     << a_vec.i[6] << ", " << a_vec.i[7] << ')';
  return a_os;
}
inline  std::ostream& operator<<(std::ostream& a_os, const __m256i& a_vec)
{
  const CH256i_t tmp = { a_vec };
  a_os << tmp;
  return a_os;
}

//--256sf

inline  std::ostream& operator<<(std::ostream& a_os, const CH256sf_t& a_vec)
{
  a_os << '('
     << a_vec.f[0] << ", " << a_vec.f[1] << ", "
     << a_vec.f[2] << ", " << a_vec.f[3] << ", "
     << a_vec.f[4] << ", " << a_vec.f[5] << ", "
     << a_vec.f[6] << ", " << a_vec.f[7] << ')';
  return a_os;
}
inline  std::ostream& operator<<(std::ostream& a_os, const __m256& a_vec)
{
  const CH256sf_t tmp = { a_vec };
  a_os << tmp;
  return a_os;
}

//--256df

inline  std::ostream& operator<<(std::ostream& a_os, const CH256df_t& a_vec)
{
  a_os << '('
     << a_vec.f[0] << ", " << a_vec.f[1] << ", "
     << a_vec.f[2] << ", " << a_vec.f[3] << ')';
  return a_os;
}
inline  std::ostream& operator<<(std::ostream& a_os, const __m256d& a_vec)
{
  const CH256df_t tmp = { a_vec };
  a_os << tmp;
  return a_os;
}

#else

//--128i (prints signed values)

inline  std::ostream& operator<<(std::ostream& a_os, const CH128i_t& a_vec)
{
  a_os << '('
     << a_vec.i[0] << ", " << a_vec.i[1] << ", "
     << a_vec.i[2] << ", " << a_vec.i[3] << ')';
  return a_os;
}
inline  std::ostream& operator<<(std::ostream& a_os, const __m128i& a_vec)
{
  const CH128i_t tmp = { a_vec };
  a_os << tmp;
  return a_os;
}

//--128sf

inline  std::ostream& operator<<(std::ostream& a_os, const CH128sf_t& a_vec)
{
  a_os << '('
     << a_vec.f[0] << ", " << a_vec.f[1] << ", "
     << a_vec.f[2] << ", " << a_vec.f[3] << ')';
  return a_os;
}
inline  std::ostream& operator<<(std::ostream& a_os, const __m128& a_vec)
{
  const CH128sf_t tmp = { a_vec };
  a_os << tmp;
  return a_os;
}

//--128df

inline  std::ostream& operator<<(std::ostream& a_os, const CH128df_t& a_vec)
{
  a_os << '('
     << a_vec.f[0] << ", " << a_vec.f[1] << ')';
  return a_os;
}
inline  std::ostream& operator<<(std::ostream& a_os, const __m128d& a_vec)
{
  const CH128df_t tmp = { a_vec };
  a_os << tmp;
  return a_os;
}

#endif

#endif
#endif
