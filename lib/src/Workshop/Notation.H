#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _NOTATION_H_
#define _NOTATION_H_

#include <map>

#include "SPACE.H"
#include "IndexTM.H"

#include "NamespaceHeader.H"

#ifdef CH_REALM

#if CH_SPACEDIM != 2
void THIS_IS_AN_ERROR_MESSAGE(void)
{
  REALM_WILL_ONLY_COMPILE_WHEN_CH_SPACEDIM_IS_2;
}
#else
#define GLOBALDIM              3
#define USING_TOP_FACE_MOMENTS 1
#endif

#else

#define GLOBALDIM              SpaceDim
#define USING_TOP_FACE_MOMENTS 0

#endif

#define RECURSIVE_GEOMETRY_GENERATION 0

#if RECURSIVE_GEOMETRY_GENERATION == 0

// Order of the Taylor expansion of the normal used to compute the EB
// geometry.  Only 0 or 1 are valid currently.
#define GEOMETRY_ORDER 1

// The degree of the moments computed on the embedded boundary
#define GEOMETRY_DEGREEP 2

#else

// Order of the Taylor expansion of the normal used to compute the EB
// geometry for the maximum degree moments
#define GEOMETRY_ORDERPMAX  1

// The maximum degree of the moments computed on the embedded boundary;
// the maximum degree of the volume moments is one less than this
#define GEOMETRY_DEGREEPMAX 2

#endif

// true -> use constaints; false -> use clipping
#define GEOMETRY_CONSTRAINTS true

// Value 1 = refining in Z direction; 0 = not refining
#define REFINEMENT_IN_ZDIRECTION 0

// Value 1 = refining based on residuals value; 0 = not refining
#define REFINEMENT_WITH_RESIDUALS 0

// Value 1 = making an EBIS; 0 = not making an EBIS
#define MAKING_EB_GEOM 1

// General constants
#define LARGEINTVAL -999999999
#define LARGEREALVAL 123456789.0
#define TOLERANCE            1.0e-10
#define MACHINEPRECISION     1.0e-14
#define REFINEMENTTHRESHOLD  1.0e-14

// These terms are possible values for cornerSigns
#define OUT -1
#define ON   0
#define IN   1

// These are for indexing into boundary maps
#define BDID_DIR   0
#define BDID_HILO  1

// Ordering for map
#ifndef USE_STCVECTOR
template<typename T> struct LexLT
{
  bool operator()(const T& a_i1,
                  const T& a_i2) const
  {
    return a_i1.lexLT(a_i2);
  }
};
#else
template<typename Vec>
using LexLT = stc::CompareLT<typename Vec::value_type, Vec::size()>;
#endif

enum EBorVol
{
  EBMoment  = 0,
  VolMoment = 1
};

// N-tuples of integers and reals
typedef IndexTM<int,2>            Iv2;
typedef IndexTM<int,GLOBALDIM>    IvgDim;
typedef IndexTM<int,GLOBALDIM-1 > IvgLess1;
typedef IndexTM<Real,GLOBALDIM>   RvgDim;

#include "NamespaceFooter.H"

#endif
