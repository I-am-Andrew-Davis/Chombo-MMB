#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _INDEXTM_H_
#define _INDEXTM_H_

//
// This code supports the chombotomization-of-Tempest project.
//
// What we have here (and in IndexTMI.H) specifically is a templatized
// unification of IntVect and RealVect.
//
// Author: Ted
//

#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include "REAL.H"
#include "Misc.H"

#define USE_STCVECTOR
#ifdef USE_STCVECTOR
#include "StcVector.H"
#else
#include "GenericArithmetic.H"
#endif

#include "BaseNamespaceHeader.H"

//class HDF5Handle;

#ifdef USE_STCVECTOR

template<typename T, stc::array_size_type N=CH_SPACEDIM>
using IndexTM = stc::Vector<T, N>;
template<typename T, stc::array_size_type N=CH_SPACEDIM>
using IndexTMAlias = stc::VectorAlias<T, N>;
template<typename T, stc::array_size_type N=CH_SPACEDIM>
using IndexTMConstAlias = stc::VectorConstAlias<T, N>;
template<typename T, stc::array_size_type N>
constexpr auto IndexTM_zero = stc::make_vector<T, N>::zero();
template<typename T, stc::array_size_type N>
constexpr auto IndexTM_unit = stc::make_vector<T, N>::unit();
template<typename T, stc::array_size_type N>
constexpr auto BASISV_TM(stc::array_size_type a_dir) noexcept
{
  return stc::make_vector<T, N>::basis(a_dir);
}

template<typename T, stc::array_size_type N, typename Op>
constexpr auto IndexTM_make_array(const Op& a_op) noexcept
{
  return stc::Make_array<T, N, Op>::with(a_op);
}

#include "BaseNamespaceFooter.H"

#else  /* ! USE_STCVECTOR  */

template<typename T, int N> class IndexTM;

template<typename T, int N> IndexTM<T,N> min(const IndexTM<T,N> & a_p1,
                                             const IndexTM<T,N> & a_p2);

template<typename T, int N> IndexTM<T,N> max(const IndexTM<T,N> & a_p1,
                                             const IndexTM<T,N> & a_p2);

template<typename T, int N> IndexTM<T,N> scale(const IndexTM<T,N> & a_p,
                                               T                    a_s);

//template<typename T, int N> IndexTM<T,N> scale(T                    a_s,
//                                               const IndexTM<T,N> & a_p);

template<typename T, int N> IndexTM<T,N> reflect(const IndexTM<T,N> & a_a,
                                                 T                    a_refIx,
                                                 int                  a_idir);

template<typename T, int N> IndexTM<T,N> diagShift(const IndexTM<T,N> & a_p,
                                                   T                    a_s);

template<typename T, int N> IndexTM<T,N> coarsen(const IndexTM<T,N> & a_p,
                                                 T                    a_s);

template<typename T, int N> IndexTM<T,N> coarsen(const IndexTM<T,N> & a_p1,
                                                 const IndexTM<T,N> & a_p2);

template<typename T, int N> std::ostream& operator<<(std::ostream       & a_os,
                                                     const IndexTM<T,N> & a_iv);

template<typename T, int N> std::istream& operator>>(std::istream & a_os,
                                                     IndexTM<T,N> & a_iv);

///
/**
   Returns a basis vector in the given coordinate direction.<br>
   In 3-D:
   BASISV(0) == (1,0,0); BASISV(1) == (0,1,0); BASISV(2) == (0,0,1).<br>
   In 2-D:
   BASISV(0) == (1,0); BASISV(1) == (0,1).<br>
   Note that the coordinate directions are based at zero.
*/
template<typename T, int N> IndexTM<T,N> BASISV_TM(int a_dir);

template<typename T, int N=CH_SPACEDIM> class IndexTM: public GenericArithmeticable< T,IndexTM<T,N> >
{
public:
  /**
     \name Constructors and Accessors
    See derived class IndexTM for more constructors.
  */
  /*@{*/

  ///
  /**
     Construct an IndexTM whose components are uninitialized.
  */
  IndexTM() : GenericArithmeticable< T, IndexTM<T,N> >(this)
  {
  }

  ///
  /**
     Destructor.
  */
  ~IndexTM()
  {
  }

  // Each of these, in its implementation, has a STATIC_ASSERT to check that
  // N is the same as the number of arguments:
  explicit
  IndexTM(T a_i);
  IndexTM(T a_i, T a_j);
  IndexTM(T a_i, T a_j, T a_k);
  IndexTM(T a_i, T a_j, T a_k, T a_l);
  IndexTM(T a_i, T a_j, T a_k, T a_l, T a_m);
  IndexTM(T a_i, T a_j, T a_k, T a_l, T a_m, T a_n);

  ///
  /**
     Construct an IndexTM setting the coordinates to the corresponding
     values in the integer array <i>a_a</i>.
  */
  explicit IndexTM(const T* a_a);

  ///
  /**
     The copy constructor.
  */
  IndexTM(const IndexTM & a_rhs);

  IndexTM copy() const
  {
    return *this;
  }

//operator IndexTM<Real,N>();

  ///
  /**
     The assignment operator.
  */
  IndexTM& operator=(const IndexTM & a_rhs);

  ///
  /**
     Returns a modifiable lvalue reference to the <i>i</i>'th coordinate of the
     IndexTM.
  */
  inline T& operator[](int a_i);

  ///
  /**
     Returns the <i>i</i>'th coordinate of the IndexTM.
  */
  inline T operator[](int a_i) const;

  ///
  /**
     Set <i>i</i>'th coordinate of IndexTM to <i>val</i>.
  */
  void setVal(int a_i,
              T   a_val);

  ///
  /**
     set all values to <i>val</i>
  */
  void setAll(T a_val);

  /*@}*/

  /**
     \name Data pointer functions
  */
  /*@{*/

  ///
  /**
     Returns a const pointer to an array of coordinates of the IndexTM.
     Useful for arguments to FORTRAN calls.
  */
  const T* getVect() const;

  ///
  /**
     Only for sending to Fortran
   */
  const T* dataPtr() const;

  ///
  /**
     Only for sending to Fortran
   */
  T* dataPtr();

  /*@}*/

  /**
     \name Comparison Operators
  */
  /*@{*/

  ///
  /**
     Returns true if this IndexTM is equivalent to argument IndexTM.  All
     comparisons between analogous components must be satisfied.
  */
  bool operator==(const IndexTM & a_p) const;

  ///
  /**
     Returns true if this IndexTM is different from argument IndexTM.
     All comparisons between analogous components must be satisfied.
  */
  bool operator!=(const IndexTM & a_p) const;

  ///
  /**
     Returns true if this IndexTM is lexically less than the argument.
     An IndexTM MUST BE either lexically less than, lexically greater
     than, or equal to another IndexTM.

     iv1 is lexically less than iv2 if:

     in 2-D:<br>
     (iv1[0] < iv2[0]) || ((iv1[0] == iv2[0]) && (iv1[1] < iv2[1]));

     in 3-D:<br>
     (iv1[0] < iv2[0]) || (iv1[0]==iv2[0] && ((iv1[1] < iv2[1] || ((iv1[1] == iv2[1]) && (iv1[2] < iv2[2])))));
  */
  bool lexLT(const IndexTM & a_s) const;

  ///
  /**
     Returns true if this IndexTM is lexically less than the argument.
     An IndexTM MUST BE either lexically less than, lexically greater
     than, or equal to another IndexTM.

     iv1 is lexically less than iv2 if:

     in 2-D:<br>
     (iv1[0] < iv2[0]) || ((iv1[0] == iv2[0]) && (iv1[1] < iv2[1]));

     in 3-D:<br>
     (iv1[0] < iv2[0]) || (iv1[0]==iv2[0] && ((iv1[1] < iv2[1] || ((iv1[1] == iv2[1]) && (iv1[2] < iv2[2])))));
  */
  bool operator<(const IndexTM & a_s) const
  {
    return lexLT(a_s);
  }

  ///
  /**
     Does operator <= the same way that intvect does it
     Intvect code looks like this:
     return D_TERM6(vect[0] >= p[0], && vect[1] >= p[1], && vect[2] >= p[2],
     && vect[3] >= p[3], && vect[4] >= p[4], && vect[5] >= p[5]);

   */
  bool componentwiseLE(const IndexTM& a_s)
  {
    bool retval = true;
    for (int idir = 0; idir < N; idir++)
      {
        retval = (retval && ((*this)[idir] <= a_s[idir]));
      }
    return retval;
  }
  ///
  /**
     Returns true if this IndexTM is lexically greater than the
     argument.  An IndexTM MUST BE either lexically less than,
     lexically greater than, or equal to another IndexTM.

     iv1 is lexically less than iv2 if:

     in 2-D:<br>
     (iv1[0] > iv2[0]) || ((iv1[0] == iv2[0]) && (iv1[1] > iv2[1]));

     in 3-D:<br>
     (iv1[0] > iv2[0]) || (iv1[0]==iv2[0] && ((iv1[1] > iv2[1] || ((iv1[1] == iv2[1]) && (iv1[2] > iv2[2])))));
  */
  bool lexGT(const IndexTM & a_s) const;

  /*@}*/

  /**
     \name Unary operators
  */
  /*@{*/

  ///
  /**
     Unary plus -- for completeness.
  */
  IndexTM operator+() const;

  ///
  /**
     Unary minus -- negates all components of this IndexTM.
  */
  IndexTM operator-() const;

  ///
  T dotProduct(const IndexTM & a_rhs) const;

  ///
  /**
     Sum of all components of this IndexTM.
  */
  T sum() const;

  ///
  /**
     Product of all components of this IndexTM.
  */
  T product() const;

  /*@}*/

  /** Used by GenericArithmeticable to implement all pointwise arithmetic
   *  and comparison functions.
  */
  template<typename OP> bool operatorCompare(const IndexTM<T,N> & a_p,
                                             const OP           & a_op) const;

  template<typename OP> IndexTM<T,N>& operatorOpEquals(const IndexTM<T,N> & a_p,
                                                       const OP           & a_op);

  template<typename OP> IndexTM<T,N>& operatorOpEquals(const T  & a_p,
                                                       const OP & a_op);

  ///
  /**
     Modifies this IndexTM by division of each component into T(1).
  */
  IndexTM& reciprocal();

  ///
  /**
     Component with the minimum value of this Index(returns 0 if they are all the same).
     a_doAbs : if true then take the absolute value before comparing
  */
  int minDir(bool a_doAbs) const;

  ///
  /**
     Component with the maximum value of this Index (returns 0 if they are all the same).
     a_doAbs : if true then take the absolute value before comparing.
  */
  int maxDir(bool a_doAbs) const;

  /**
     \name Other arithmetic operators
  */
  /*@{*/

  ///
  /**
     Modifies this IndexTM by taking component-wise min with IndexTM
     argument.
  */
  IndexTM& min(const IndexTM & a_p);

  ///
  /**
     Modifies this IndexTM by taking component-wise max with IndexTM
     argument.
  */
  IndexTM& max(const IndexTM & a_p);

  ///
  /**
     Modifies this IndexTM by multiplying each component by a scalar.
  */
  IndexTM& scale(T a_s);

  ///
  /**
     Modifies IndexTM by reflecting it in the plane defined by the
     index <i>ref_ix</i> and with normal in the direction of <i>idir</i>.
     Directions are based at zero.
  */
  IndexTM& reflect(T   a_refIx,
                   int a_idir);

  ///
  /**
     Modifies this IndexTM by adding <i>s</i> to component in given coordinate
     direction.
  */
  IndexTM& shift(int a_coord,
                 T   a_s);

  ///
  /**
     Modifies this IndexTM by component-wise addition with IndexTM
     argument.
  */
  IndexTM& shift(const IndexTM & a_iv);

  ///
  /**
     Modifies this IndexTM by adding a scalar <i>s</i> to each component.

  */
  IndexTM& diagShift(T a_s);

  ///
  /**
     Modify IndexTM by component-wise integer projection.
  */
  IndexTM& coarsen(const IndexTM & a_p);

  ///
  /**
     Modify IndexTM by component-wise integer projection.
  */
  IndexTM& coarsen(T a_p);

  /*@}*/

  /**
     \name I/O Functions
  */
  /*@{*/

  ///
  /**
     Print an IndexTM to the ostream.
  */
  void printOn(std::ostream & a_os) const;

  ///
  /**
     Print an IndexTM to the pout().
  */
  void p() const;

  ///
  /**
     Print an IndexTM to the ostream a bit more verbosely.
  */
  void dumpOn(std::ostream & a_os) const;

  ///
  /**
     Print the IndexTM to given output stream in ASCII.
  */
  friend std::ostream& operator<< <>(std::ostream  & a_os,
                                     const IndexTM & a_iv);

  ///
  /**
     Read next IndexTM from given input stream.
  */
  friend std::istream& operator>> <> (std::istream & a_os,
                                      IndexTM      & a_iv);

  /*@}*/

  /**
     \name Constants
  */
  /*@{*/

  /**
     This is an IndexTM all of whose components are equal to zero.
  */
  static const IndexTM Zero;

  /**
     This is an IndexTM all of whose components are equal to one.
  */
  static const IndexTM Unit;

  /**
     Initializes Zero and Unit.
  */
  static int InitStatics();

  /**
    Low-level data copy.
  */
  void linearIn(const void* a_inBuf);

  void linearOut(void* a_outBuf) const;

  /*@}*/

  ///
  T integer_factorial(int n) const
  {
    CH_assert(n >= 0);
    T retval;
    if(n== 1 || n== 0)
      {
        retval = 1;
      }
    else
      {
        retval = integer_factorial(n-1)*n;
      }
    return retval;
  }
  
  ///
  T factorial() const
  {
    int retval = 1;
    for(int idir = 0; idir < N; idir++)
      {
        retval *= integer_factorial(m_vect[idir]);
      }
    return retval;
  }
protected:
  IndexTM(const char*);// used to build Zero and Unit static members

  //  friend class HDF5Handle;

  /**
     The individual components of this IndexTM.
  */
  T m_vect[N];
};

/**
    Useful for preventing division by IntVects, e.g.
        CH_assert( ! IndexTraits<T>::m_isInt );
    or better yet,
        STATIC_ASSERT( ! IndexTraits<T>::m_isInt );
*/
template<typename T> struct IndexTraits
{
    static bool const m_isInt=false;
    static bool const m_isReal=false;
};
template<> struct IndexTraits<int>
{
    static bool const m_isInt=true;
    static bool const m_isReal=false;
};
template<> struct IndexTraits<Real>
{
    static bool const m_isInt=false;
    static bool const m_isReal=true;
};

//
// Static initialization.  Gotta make sure there are no static object
// definitions above here (except possibly stuff in headers).  Otherwise,
// the danger is that some static object's constructor calls IndexTM::Zero or
// IndexTM::Unit -- the very things the following definition is supposed to
// initialize.
//
/*
static int s_dummyForIntVectH1 (IndexTM<int ,1>::InitStatics());
static int s_dummyForIntVectH2 (IndexTM<int ,2>::InitStatics());
static int s_dummyForIntVectH3 (IndexTM<int ,3>::InitStatics());
static int s_dummyForIntVectH4 (IndexTM<int ,4>::InitStatics());
static int s_dummyForIntVectH5 (IndexTM<int ,5>::InitStatics());
static int s_dummyForIntVectH6 (IndexTM<int ,6>::InitStatics());
static int s_dummyForRealVectH1(IndexTM<Real,1>::InitStatics());
static int s_dummyForRealVectH2(IndexTM<Real,2>::InitStatics());
static int s_dummyForRealVectH3(IndexTM<Real,3>::InitStatics());
static int s_dummyForRealVectH4(IndexTM<Real,4>::InitStatics());
static int s_dummyForRealVectH5(IndexTM<Real,5>::InitStatics());
static int s_dummyForRealVectH6(IndexTM<Real,6>::InitStatics());
*/

template <typename T, int N> const IndexTM<T,N> IndexTM<T,N>::Zero("0");
template <typename T, int N> const IndexTM<T,N> IndexTM<T,N>::Unit("1");

#include "BaseNamespaceFooter.H"

#include "IndexTMI.H"

#endif  /* ! USE_STCVECTOR */

#endif  /* ! defined _INDEXTM_H_ */
