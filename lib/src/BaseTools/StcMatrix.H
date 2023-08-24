
#ifndef _STCMATRIX_H_
#define _STCMATRIX_H_


/******************************************************************************/
/**
 * \file StcMatrix.H
 *
 * \brief Static square matrix class
 *
 * \note
 * <ul>
 *   <li> The term 'rank', with symbol R, is used herein to mean the number
 *        of rows or columns in the matrix.  The term 'size', with symbol N,
 *        refers to the total number of elements and is N=R*R
 *   <li> In the code, matrix implementations are actual matrices based on
 *        either std::array or an alias and having element values directly of
 *        type T.  Matrix expression involve some type of expression for
 *        the element values.  A matrix representation can be either an
 *        implementation or an expression.
 *   <li> Function objects that are used for immediate evaluation have the
 *        suffix _E.  One could use lambdas except that constexpr is not
 *        supported until C++17.  Also, there is sometimes reuse.
 * </ul>
 *
 *//*+*************************************************************************/

#include <iomanip>

#include "StcVector.H"

/*--------------------------------------------------------------------*
 * Select the storage order.  Default is column major.  Comment out
 * define on next line for row major ordering.
 *--------------------------------------------------------------------*/

#if !defined (STC_COLMAJOR) && !defined (STC_ROWMAJOR)
  #define STC_COLMAJOR
#endif
#ifdef STC_COLMAJOR
  #define STC_RO(i, stride) (i)           // To linear
  #define STC_CO(j, stride) (j)*(stride)
  #define STC_RI(idx, R) (idx)%(R)        // From linear
  #define STC_CI(idx, R) (idx)/(R)
#else
  #define STC_RO(i, stride) (i)*(stride)  // To linear
  #define STC_CO(j, stride) (j)
  #define STC_RI(idx, R) (idx)/(R)        // From linear
  #define STC_CI(idx, R) (idx)%(R)
#endif

namespace stc
{

/*--------------------------------------------------------------------*
 * Makers for matrices
 *
 * A bit of wizardry to initialize the array.  The approach was taken
 * from the way make_index_sequence itself is implemented in
 * gcc-6.1.0.
 *
 * Make_array is defined in StcVector.H where a_op() takes one
 * argument.  For Make_array2, a_op() takes two arguments (row and
 * column of matrix).
 *--------------------------------------------------------------------*/

/*
 *  Default maker's only purpose is to initialize last template parameter
 *  \tparam     T       Type of array
 *  \tparam     R       Rank of square matrix
 *  \tparam     Op      Function object to initialize element of array
 *  \tparam     ISeq    Generated integer sequence 0..N-1 for  array
 */
template <typename T, array_size_type R, typename Op,
          typename ISeq = std::make_integer_sequence<array_size_type, R*R>>
struct Make_array2;

/*
 *  Specialization matches generated integer sequence
 *  The generated sequence has type std::index_sequence<Is...> and
 *  now this specialization matches
 */
template <typename T, array_size_type R, typename Op, array_size_type... Is>
struct Make_array2<T, R, Op, std::integer_sequence<array_size_type, Is...>>
{
  /*
   *  \param[in]  a_op    Function applied to initialize element of
   *                      array
   */
  static constexpr DefaultImpl<T, R*R> with(const Op& a_op)
    {
      // op is applied to each integer in the sequence
      return { a_op(STC_RI(Is, R), STC_CI(Is, R))... };
    }
};

/*
 *  Interface to Make_array2 where Op can be deduced
 *  \tparam     T       Type of array
 *  \tparam     R       Rank of square matrix
 *  \tparam     Op      Object type to initialize element of array
 *  \param[in]  a_op    Function applied to initialize element of
 *                      array
 */
template <typename T, array_size_type R, typename Op>
HOSTDEVICE constexpr auto
make_array2(const Op& a_op)
{
  return Make_array2<T, R, Op>::with(a_op);
}


/*******************************************************************************
 *
 * Class MatrixImpl
 *
 * Usually built around a std::array component (or similar) but aliases are also
 * supported.  These static arrays can be any size and of any type.  Expressions
 * are also provided.
 *
 ******************************************************************************/

//--Alias implementations.  These take the same parameters as std::array and
//--discard irrelevant.

/// Matrix alias implementation of pointer
template <typename T, array_size_type>
using MatrixAliasImpl = T*;

/// Matrix alias implementation of constant pointer
template <typename T, array_size_type>
using MatrixConstAliasImpl = const T*;


/*============================================================================*/
/// Static matrix
/** Normally, this class is not used directly.  There are numerous typedefs
 *  that create more accessible objects.
 *  \tparam     T       Type of array
 *  \tparam     R       Rank of square matrix
 *  \tparam     Impl    Underlying implementation.  Default is DefaultImpl
 *                      (std::array or similar) but aliases are also possible
 *                      (use a type from above).
 *//*=========================================================================*/

template <typename T,
          array_size_type R,
          template<typename, array_size_type> class Impl = DefaultImpl>
class MatrixImpl : StcOp
{

/*--------------------------------------------------------------------*
 * Types
 *--------------------------------------------------------------------*/

public:

  using value_type = T;
  using is_scalar = std::false_type;
  using is_vector = std::false_type;
  using is_rowvec = std::false_type;
  using is_matrix = std::true_type;
  using is_string = std::false_type;
  static constexpr array_size_type N = R*R;

/*--------------------------------------------------------------------*
 * Constructors and destructors
 *--------------------------------------------------------------------*/

  /// Default construct default initializes (uninitialized data)
  HOSTDEVICE MatrixImpl()
    { }

  /// Construct with a constant T
  // This cannot be used in implicit conversions
  HOSTDEVICE explicit MatrixImpl(const T a_x)
    :
    m_data(make_array<T, N>(Uniform_E(a_x)))
    { }

  /// Constructor that only builds from rvalues.
  /** Use with make_matrix helper functions
   */
  HOSTDEVICE constexpr MatrixImpl(DefaultImpl<T, N>&& a_data) noexcept
    :
    m_data(std::move(a_data))
    { }

  /// Constructor taking pointer (use at your own risk)
  /** This does not make an alias.  Data is taken from the array at a_ptr.
   */
  HOSTDEVICE explicit MatrixImpl(const T* a_ptr) noexcept
    :
    m_data(make_array<T, N>(
             [a_ptr]
             (const array_size_type a_idx)
             { return a_ptr[a_idx]; }))
    { }

  // Helper function object for initializer_list constructor.  Used instead of
  // lambda to allow constexpr.
  struct IL2_E
  {
    HOSTDEVICE constexpr IL2_E(std::initializer_list<T> a_il) noexcept
      :
      m_il(a_il)
      { }
    HOSTDEVICE constexpr T operator()(const array_size_type a_i,
                                      const array_size_type a_j)
      const noexcept
      {
        // Note: the IL is always in row-major order
        return (a_i*R + a_j < m_il.size()) ? m_il.begin()[a_i*R + a_j] : T(0);
      }
    const std::initializer_list<T> m_il;
  };

  /// Constructor taking initializer... and forwarded to Make_array2
  HOSTDEVICE constexpr MatrixImpl(std::initializer_list<T> a_il) noexcept
    :
    m_data(make_array2<T, R>(IL2_E(a_il)))
    { }

  /// Copy constructor
  MatrixImpl(const MatrixImpl&) = default;
  /// Move constructor
  MatrixImpl(MatrixImpl&&) noexcept = default;
  /// Assignment constructor
  MatrixImpl& operator=(const MatrixImpl&) = default;
  /// Move assignment constructor
  MatrixImpl& operator=(MatrixImpl&&) noexcept = default;
  /// Destructor
  ~MatrixImpl() = default;

//--Alias (via constructor): if source is matrix and destination is alias

  /// Alias a const matrix
  /** \tparam     S       Type of source matrix, needed for enable_if
  *                       but otherwise S must equal T
  *   \param[in]  a_mat   Source matrix
  */
  template <typename S,
            typename = std::enable_if_t<
              std::is_same<Impl<S, N>, MatrixConstAliasImpl<S, N> >::value>
            >
  HOSTDEVICE constexpr MatrixImpl(
    const MatrixImpl<S, R, DefaultImpl>& a_mat) noexcept
    :
    m_data(a_mat.m_data.data())
    {
      static_assert(std::is_same<S, T>::value, "Matrix alias must have same "
                    "type");
    }

  /// Alias a non-const vector
  /** \tparam     S       Type of source matrix, needed for enable_if
  *                       but otherwise S must equal T
  *   \param[in]  a_mat   Source matrix
  */
  template <typename S,
            typename = std::enable_if_t<
              std::is_same<Impl<S, N>, MatrixAliasImpl<S, N> >::value>
            >
  HOSTDEVICE MatrixImpl(MatrixImpl<S, R, DefaultImpl>& a_mat) noexcept
    :
    m_data(a_mat.m_data.data())
    {
      static_assert(std::is_same<S, T>::value, "Matrix alias must have same "
                    "type");
    }

//--Copy: if source is alias and destination is matrix

  /// Source is a const alias
  /** \tparam     S       Type of source matrix, needed for enable_if
  *                       but otherwise S must equal T
  *   \param[in]  a_mat   Source matrix
  */
  template <typename S,
            typename = std::enable_if_t<
              std::is_same<Impl<S, N>, DefaultImpl<S, N> >::value>
            >
  HOSTDEVICE constexpr MatrixImpl(
    const MatrixImpl<S, R, MatrixConstAliasImpl>& a_alias) noexcept
    :
    m_data(make_array<T, N>(
             typename Vector<T, N>::template ArrayCopy_E<
               MatrixConstAliasImpl<S, N>, N
             >(a_alias.m_data)))
    {
      static_assert(std::is_same<S, T>::value, "Matrix alias must have same "
                    "type");
    }

  /// Source is a non-const alias
  /** \tparam     S       Type of source matrix, needed for enable_if
  *                       but otherwise S must equal T
  *   \param[in]  a_mat   Source matrix
  */
  template <typename S,
            typename = std::enable_if_t<
              std::is_same<Impl<S, N>, DefaultImpl<S, N> >::value>
            >
  HOSTDEVICE constexpr MatrixImpl(
    const MatrixImpl<S, R, MatrixAliasImpl>& a_alias) noexcept
    :
    m_data(make_array<T, N>(
             typename Vector<T, N>::template ArrayCopy_E<
               MatrixAliasImpl<S, N>, N
             >(a_alias.m_data)))
    {
      static_assert(std::is_same<S, T>::value, "Matrix alias must have same "
                    "type");
    }

  /// Construct from an expression template given by REP
  /** This is never called if RHS is any matrix implementation
   *  \tparam     REP     Representation of expression
   *  \param[in]  a_rhs   Right-hand-side expression
   */
  template <typename REP,
            typename = std::enable_if_t<
              !std::is_base_of<MatrixImpl<T, R, DefaultImpl>,
                               std::decay_t<REP> >::value
              &&
              !std::is_base_of<MatrixImpl<T, R, MatrixConstAliasImpl>,
                               std::decay_t<REP> >::value
              &&
              !std::is_base_of<MatrixImpl<T, R, MatrixAliasImpl>,
                               std::decay_t<REP> >::value
              &&
              std::is_base_of<StcOp, std::decay_t<REP> >::value>
            >
  HOSTDEVICE constexpr MatrixImpl(REP&& a_rhs) noexcept
    :
    m_data(make_array<T, N>(
             typename Vector<T, N>::template ArrayCopy_E<
               std::decay_t<REP>,
               cmin(MatrixImpl::size(), std::decay_t<REP>::size())
             >(a_rhs)))
    { }

/*--------------------------------------------------------------------*
 * Member functions
 *--------------------------------------------------------------------*/

  /// Total number of elements in the matrix
  HOSTDEVICE static constexpr array_size_type size() noexcept
    {
      return N;
    }

  /// Rank of the matrix
  HOSTDEVICE static constexpr array_size_type rank() noexcept
    {
      return R;
    }

  /// Assignment that evaluates an expression given by REP (=)
  /** This is never called if RHS is any matrix implementation
   *  \tparam     REP     Representation of expression
   *  \param[in]  a_rhs   Right-hand-side expression
   */
  template <typename REP,
            typename = std::enable_if_t<
              !std::is_base_of<MatrixImpl<T, R, DefaultImpl>,
                               std::decay_t<REP> >::value
              &&
              !std::is_base_of<MatrixImpl<T, R, MatrixConstAliasImpl>,
                               std::decay_t<REP> >::value
              &&
              !std::is_base_of<MatrixImpl<T, R, MatrixAliasImpl>,
                               std::decay_t<REP> >::value
              &&
              std::is_base_of<StcOp, std::decay_t<REP> >::value>
            >
  HOSTDEVICE MatrixImpl& operator=(REP&& a_rhs)
    {
      typename Vector<T, N>::template ArrayEq_E<std::decay_t<REP>>
        f(m_data, a_rhs);
      constexpr array_size_type I = cmin(MatrixImpl::size(),
                                         std::decay_t<REP>::size());
      // If not for nvcc, here we could use a lambda since LHS cannot be
      // constexpr, e.g.:
      // auto f = [&](const array_size_type a_idx)
      //   {
      //     m_data[a_idx] = a_rhs[a_idx];
      //   };
      forEachElement<I>(f);
      return *this;
    }

  /// Evaluate an expression template given by Rep (+=)
  template <typename REP,
            typename = std::enable_if_t<
              std::is_base_of<StcOp, std::decay_t<REP> >::value>
            >
  HOSTDEVICE MatrixImpl& operator+=(REP&& a_rhs) noexcept
    {
      this->operator=(*this + static_cast<const REP&>(a_rhs));
      return *this;
    }

  /// Evaluate an expression template given by Rep (-=)
  template <typename REP,
            typename = std::enable_if_t<
              std::is_base_of<StcOp, std::decay_t<REP> >::value>
            >
  HOSTDEVICE MatrixImpl& operator-=(REP&& a_rhs) noexcept
    {
      this->operator=(*this - static_cast<const REP&>(a_rhs));
      return *this;
    }

//--x= scalars

  /// Assign constant
  HOSTDEVICE MatrixImpl& operator=(const T& a_val) noexcept
    {
      this->operator=(stc::ScalarRep<T>(a_val));
      return *this;
    }

  /// Add scalar to each element
  HOSTDEVICE MatrixImpl& operator+=(const T& a_val) noexcept
    {
      this->operator=(*this + a_val);
      return *this;
    }

  /// Subtract scalar from each element
  HOSTDEVICE MatrixImpl& operator-=(const T& a_val) noexcept
    {
      this->operator=(*this - a_val);
      return *this;
    }

  /// Multiply scalar with each element
  HOSTDEVICE MatrixImpl& operator*=(const T& a_val) noexcept
    {
      this->operator=(*this * a_val);
      return *this;
    }

  /// Divide each element by scalar
  HOSTDEVICE MatrixImpl& operator/=(const T& a_val) noexcept
    {
      this->operator=(*this / a_val);
      return *this;
    }

//--Access

  /// Constant linear indexing
  HOSTDEVICE constexpr const T& operator[](
    const array_size_type a_idx) const noexcept
    {
      return m_data[a_idx];
    }

  /// Linear indexing
  HOSTDEVICE T& operator[](const array_size_type a_idx) noexcept
    {
      return m_data[a_idx];
    }

  /// Constant indexing (function/"lambda-like")
  HOSTDEVICE constexpr const T& operator()(
    const array_size_type a_idx) const noexcept
    {
      return m_data[a_idx];
    }

  /// Indexing (function/"lambda-like")
  HOSTDEVICE constexpr T& operator()(const array_size_type a_idx) noexcept
    {
      return m_data[a_idx];
    }

  /// Constant matrix indexing
  HOSTDEVICE constexpr const T& operator()(
    const array_size_type a_i,
    const array_size_type a_j) const noexcept
    {
      return m_data[STC_RO(a_i, R) + STC_CO(a_j, R)];
    }

  /// Matrix indexing
  HOSTDEVICE T& operator()(const array_size_type a_i,
                           const array_size_type a_j) noexcept
    {
      return m_data[STC_RO(a_i, R) + STC_CO(a_j, R)];
    }

  /// Constant data pointer
  HOSTDEVICE const T* dataPtr() const noexcept
    {
      return &m_data[0];
    }

  /// Data pointer
  HOSTDEVICE T* dataPtr() noexcept
    {
      return &m_data[0];
    }

/*--------------------------------------------------------------------*
 * Member data
 *--------------------------------------------------------------------*/

//--Friend declarations so different implementations can access data

  template <typename T1, array_size_type R1,
            template <typename, array_size_type> class Impl1>
  friend class MatrixImpl;

protected:

  Impl<T, N> m_data;                  ///< Underlying implementation giving data
                                      ///< for the matrix
};


/*==============================================================================
 * Types of matrices that are more accessible
 *
 * Matrix<T, R>
 * MatrixAlias<T, R>
 * MatrixConstAlias<T, R>
 *
 * Integer
 * IMat<R>
 * IMatAlias<R>
 * IMatConstAlias<R>
 *
 * Real
 * RMat<R>
 * RMatAlias<R>
 * RMatConstAlias<R>
 *
 * Spatial (Real, R=3) Real
 * SMat
 * SMatAlias
 * SMatConstAlias
 *
 * Transformation (Real, R=4)
 * TMat
 * TMatAlias
 * TMatConstAlias
 *
 *============================================================================*/

//--Matrix of any type and size

template <typename T, array_size_type R>
using Matrix = MatrixImpl<T, R, DefaultImpl>;
template <typename T, array_size_type R>
using MatrixAlias = MatrixImpl<T, R, MatrixAliasImpl>;
template <typename T, array_size_type R>
using MatrixConstAlias = MatrixImpl<T, R, MatrixConstAliasImpl>;

//--Integer matrix of any size

template <array_size_type R>
using IMat = MatrixImpl<int, R, DefaultImpl>;
template <array_size_type R>
using IMatAlias = MatrixImpl<int, R, MatrixAliasImpl>;
template <array_size_type R>
using IMatConstAlias = MatrixImpl<int, R, MatrixConstAliasImpl>;

//--Real matrix of any size

template <array_size_type R>
using RMat = MatrixImpl<Real, R, DefaultImpl>;
template <array_size_type R>
using RMatAlias = MatrixImpl<Real, R, MatrixAliasImpl>;
template <array_size_type R>
using RMatConstAlias = MatrixImpl<Real, R, MatrixConstAliasImpl>;

//--Spatial matrix of type Real and size spatial dimension

using SMat = MatrixImpl<Real, SPACEDIM, DefaultImpl>;
using SMatAlias = MatrixImpl<Real, SPACEDIM, MatrixAliasImpl>;
using SMatConstAlias = MatrixImpl<Real, SPACEDIM, MatrixConstAliasImpl>;

//--Transformation matrix of type Real and size spatial dimension + 1

using TMat = MatrixImpl<Real, SPACEDIM+1, DefaultImpl>;
using TMatAlias = MatrixImpl<Real, SPACEDIM+1, MatrixAliasImpl>;
using TMatConstAlias = MatrixImpl<Real, SPACEDIM+1, MatrixConstAliasImpl>;


/*******************************************************************************
 */
/// Makers for matrices
/**
 *  Each maker constructs a special type of Matrix, e.g, identity matrix.  Each
 *  maker function uses a function object.  Lambdas are avoided since they
 *  cannot be made constexpr until C++17.
 *
 *  \tparam     T       Type of array
 *  \tparam     R       Rank of square matrix
 *
 *//*+*************************************************************************/

/// Maker class supporting Matrix Impl
template <typename T, array_size_type R>
struct make_matrix
{
  static constexpr array_size_type N = R*R;
  /*
   * Using a lambda is ideal but we cannot make the constexpr until C++17
   */
  // static auto zero()
  //   {
  //     return MatrixImpl<T, R, DefaultImpl>(
  //       make_array<T, N>(
  //         [](const array_size_type a_idx)
  //         {
  //           return (T)0;
  //         }));
  //   }
  /*
   * So instead, our function objects are declared before each function
   */

  struct Zero_E
  {
    constexpr T operator()(const array_size_type a_idx) const noexcept
      {
        return (T)0;
      }
  };
  /// Set all elements to 0
  static constexpr auto zero() noexcept
    {
      return Matrix<T, R>(make_array<T, N>(Zero_E{}));
    }

  struct Unit_E
  {
    constexpr T operator()(const array_size_type a_idx) const noexcept
      {
        return (T)1;
      }
  };
  /// Set all elements to 1
  static constexpr auto unit() noexcept
    {
      return Matrix<T, R>(make_array<T, N>(Unit_E{}));
    }

  struct Identity_E
  {
    constexpr T operator()(const array_size_type a_i,
                           const array_size_type a_j) const noexcept
      {
        return (T)(a_i == a_j);
      }
  };
  /// Set identity matrix
  static constexpr auto identity() noexcept
    {
      return Matrix<T, R>(make_array2<T, R>(Identity_E{}));
    }

  struct FillDiagonal_E
  {
    constexpr FillDiagonal_E(const T& a_val) noexcept
      :
      m_val(a_val)
      { }
    constexpr T operator()(const array_size_type a_i,
                           const array_size_type a_j) const noexcept
      {
        return (T)((a_i == a_j)*m_val);
      }
    T m_val;
  };
  /// Fill diagonal matrix
  static constexpr auto fillDiagonal(const T& a_val) noexcept
    {
      return Matrix<T, R>(make_array2<T, R>(FillDiagonal_E(a_val)));
    }

  struct FillUpper_E
  {
    constexpr FillUpper_E(const T& a_val) noexcept
      :
      m_val(a_val)
      { }
    constexpr T operator()(const array_size_type a_i,
                           const array_size_type a_j) const noexcept
      {
        return (T)((R - a_i + a_j > R)*m_val);
      }
    T m_val;
  };
  /// Fill upper part of matrix
  static constexpr auto fillUpper(const T& a_val) noexcept
    {
      return Matrix<T, R>(make_array2<T, R>(FillUpper_E(a_val)));
    }

  struct FillLower_E
  {
    constexpr FillLower_E(const T& a_val) noexcept
      :
      m_val(a_val)
      { }
    constexpr T operator()(const array_size_type a_i,
                           const array_size_type a_j) const noexcept
      {
        return (T)((R - a_i + a_j < R)*m_val);
      }
    T m_val;
  };
  /// Fill lower part of matrix
  static constexpr auto fillLower(const T& a_val) noexcept
    {
      return Matrix<T, R>(make_array2<T, R>(FillLower_E(a_val)));
    }

  struct FillLDU_E
  {
    constexpr FillLDU_E(const T& a_valL,
                        const T& a_valD,
                        const T& a_valU) noexcept
      :
      m_valL(a_valL),
      m_valD(a_valD),
      m_valU(a_valU)
      { }
    constexpr T operator()(const array_size_type a_i,
                           const array_size_type a_j) const noexcept
      {
        return (T)((R - a_i + a_j < R)*m_valL +
                   (a_i == a_j)       *m_valD +
                   (R - a_i + a_j > R)*m_valU);
      }
    T m_valL;
    T m_valD;
    T m_valU;
  };
  /// Fill LDU sections with different values
  static constexpr auto fillLDU(const T& a_valL,
                                const T& a_valD,
                                const T& a_valU) noexcept
    {
      return Matrix<T, R>(
        make_array2<T, R>(FillLDU_E(a_valL, a_valD, a_valU)));
    }

  struct Fill_E
  {
    constexpr Fill_E(const T& a_val) noexcept
      :
      m_val(a_val)
      { }
    constexpr T operator()(const array_size_type a_idx) const noexcept
      {
        return m_val;
      }
    T m_val;
  };
  /// Fill all values
  static constexpr auto fill(const T& a_val) noexcept
    {
      return Matrix<T, R>(make_array<T, N>(Fill_E(a_val)));
    }
  /// Antisymmetric (requires R=3)
  // static constexpr auto antiSymmetric(const Vector<T, 3>& a_vec)
  //   {
  //     static_assert(R == 3);
  //   }
};

//--More accessible definitions

template <array_size_type R>
using make_IMat = make_matrix<int, R>;
template <array_size_type R>
using make_RMat = make_matrix<Real, R>;
using make_SMat = make_matrix<Real, SPACEDIM>;
using make_TMat = make_matrix<Real, 4>;

//--Handy definitions

constexpr auto SMat_zero = make_SMat::zero();
constexpr auto TMat_zero = make_TMat::zero();
constexpr auto SMat_identity = make_SMat::identity();
constexpr auto TMat_identity = make_TMat::identity();


/*******************************************************************************
 *
 * Operator Representations
 *
 * There is a single representation for operator types, e.g., binary and unary.
 * Each operator is defined by a function that creates the appropriate
 * representations with a function object.  In simple cases, lambdas could be
 * used except we make use of constexpr and that is not supported until C++17.
 * You will find the function objects before each operator function.
 *
 * We are quite careful about moving rvalues (temporary objects) which are
 * stored in class, and copying lvalue references, which are stored as
 * references.  In otherwords, anything that started as an lvalue is referred to
 * throughout the chain of representations; anything that started as an rvalue
 * is stored in the representation class and moved when required (the original
 * ceases to exist).
 *
 * Enable_if is used in the operators.  The operands must derive from StcOp.
 * Otherwise, the operators bind to absolutely everything...
 *
 * Some operators are complex (e.g, matrix multiply) and are immediately
 * evaluated (they do not created expressions).  Because the operands are
 * used repeatedly, the are evaluated if they are expressions.  Type
 * dispatch may be used to avoid storing temporary intermediate operanads if
 * the operands are initially both matrix implementations.
 *
 * Operations may return matrix implementations or expressions.  As such
 * one should avoid code such as
 *   auto matC = matA + matB;
 * since matC may actually be an expression type.  Instead do
 *   Matrix<Real, 3> matC = matA + matB;
 *
 ******************************************************************************/


/*============================================================================*/
/// Representation of a binary operation between matrices
/** \tparam REP1        Left operand
 *  \tparam REP2        Right operand
 *  \tparam BinaryOp    Function object performing a desired binary operation
 *
 *  \note
 *  <ul>
 *    <li> The function object operates on the underyling types which must be
 *         available through the representations
 *  </ul>
 *
 *//*=========================================================================*/

//--For basic Ops, the representations are evaluated at the index before
//--invoking the UnaryOp operator.

template <typename REP1, typename REP2, typename BinaryOp>
class MatrixBinaryOpRep : public BinaryOp, private StcOp
{
public:
  using value_type = typename BinaryOp::value_type;
  using is_scalar = std::false_type;
  using is_vector = std::false_type;
  using is_rowvec = std::false_type;
  using is_matrix = std::true_type;
  using is_string = std::false_type;
private:
  REP1 m_rep1;
  REP2 m_rep2;
public:
  template <typename... OpArgs>
  HOSTDEVICE constexpr MatrixBinaryOpRep(REP1&&      a_rep1,
                                         REP2&&      a_rep2,
                                         OpArgs&&... a_opArgs) noexcept
    :
    BinaryOp(std::forward<OpArgs>(a_opArgs)...),
    m_rep1(std::forward<REP1>(a_rep1)),
    m_rep2(std::forward<REP2>(a_rep2))
    { }
  HOSTDEVICE constexpr MatrixBinaryOpRep(
    const MatrixBinaryOpRep& a_other) noexcept
    :
    BinaryOp(a_other),
    m_rep1(a_other.m_rep1),
    m_rep2(a_other.m_rep2)
    { }
  HOSTDEVICE constexpr MatrixBinaryOpRep(MatrixBinaryOpRep&& a_other) noexcept
    :
    BinaryOp(std::move(a_other)),
    m_rep1(std::forward<REP1>(a_other.m_rep1)),
    m_rep2(std::forward<REP2>(a_other.m_rep2))
    { }
  /// Assignment constructor
  MatrixBinaryOpRep& operator=(const MatrixBinaryOpRep&) = delete;
  /// Move assignment constructor
  MatrixBinaryOpRep& operator=(MatrixBinaryOpRep&&) noexcept = delete;
  /// Size (number of elements) of the matrix representation (minimum of Reps)
  HOSTDEVICE static constexpr array_size_type size() noexcept
    {
      return cmin(std::decay_t<REP1>::size(), std::decay_t<REP2>::size());
    }
  /// Rank of the matrix representation (minimum of Reps)
  HOSTDEVICE static constexpr array_size_type rank() noexcept
    {
      return cmin(std::decay_t<REP1>::rank(), std::decay_t<REP2>::rank());
    }

//--Access

  /// Array indexing
  HOSTDEVICE constexpr auto operator[](
    const array_size_type a_idx) const noexcept
    {
      return this->eval(m_rep1[a_idx], m_rep2[a_idx]);
    }
  /// Function/"lambda-like"
  HOSTDEVICE constexpr auto operator()(
    const array_size_type a_idx) const noexcept
    {
      return this->eval(m_rep1[a_idx], m_rep2[a_idx]);
    }
  /// Matrix indexing
  HOSTDEVICE constexpr auto operator()(const array_size_type a_i,
                                       const array_size_type a_j) const noexcept
    {
      return this->eval(m_rep1[STC_RO(a_i, std::decay_t<REP1>::rank()) +
                               STC_CO(a_j, std::decay_t<REP1>::rank())],
                        m_rep2[STC_RO(a_i, std::decay_t<REP2>::rank()) +
                               STC_CO(a_j, std::decay_t<REP2>::rank())]);
    }
};

/*============================================================================*/
/// Representation of a unary operation on a matrix
/** \tparam REP         Operand
 *  \tparam UnaryOp     Function object performing a desired unary operation
 *
 *  \note
 *  <ul>
 *    <li> The function object operates on the underyling type which must be
 *         available through the representation
 *  </ul>
 *
 *//*=========================================================================*/

//--For basic Ops, the representations are evaluated at the index before
//--invoking the UnaryOp operator.

template <typename REP, typename UnaryOp>
class MatrixUnaryOpRep : public UnaryOp, private StcOp
{
public:
  using value_type = typename UnaryOp::value_type;
  using is_scalar = std::false_type;
  using is_vector = std::false_type;
  using is_rowvec = std::false_type;
  using is_matrix = std::true_type;
  using is_string = std::false_type;
private:
  REP m_rep;
public:
  template <typename... OpArgs>
  HOSTDEVICE constexpr MatrixUnaryOpRep(REP&&       a_rep,
                                        OpArgs&&... a_opArgs) noexcept
    :
    UnaryOp(std::forward<OpArgs>(a_opArgs)...),
    m_rep(std::forward<REP>(a_rep))
    { }
  HOSTDEVICE constexpr MatrixUnaryOpRep(
    const MatrixUnaryOpRep& a_other) noexcept
    :
    UnaryOp(a_other),
    m_rep(a_other.m_rep)
    { }
  HOSTDEVICE constexpr MatrixUnaryOpRep(MatrixUnaryOpRep&& a_other) noexcept
    :
    UnaryOp(std::move(a_other)),
    m_rep(std::forward<REP>(a_other.m_rep))
    { }
  /// Assignment constructor
  MatrixUnaryOpRep& operator=(const MatrixUnaryOpRep&) = delete;
  /// Move assignment constructor
  MatrixUnaryOpRep& operator=(MatrixUnaryOpRep&&) noexcept = delete;
  /// Size (number of elements) of the matrix representation
  HOSTDEVICE static constexpr array_size_type size() noexcept
    {
      return std::decay_t<REP>::size();
    }
  /// Rank of the matrix representation
  HOSTDEVICE static constexpr array_size_type rank() noexcept
    {
      return std::decay_t<REP>::rank();
    }

//--Access

  /// Array indexing
  HOSTDEVICE constexpr auto operator[](
    const array_size_type a_idx) const noexcept
    {
      return this->eval(m_rep[a_idx]);
    }
  /// Function/"lambda-like"
  HOSTDEVICE constexpr auto operator()(
    const array_size_type a_idx) const noexcept
    {
      return this->eval(m_rep[a_idx]);
    }
  /// Matrix indexing
  HOSTDEVICE constexpr auto operator()(const array_size_type a_i,
                                       const array_size_type a_j) const noexcept
    {
      return this->eval(m_rep[STC_RO(a_i, std::decay_t<REP>::rank()) +
                              STC_CO(a_j, std::decay_t<REP>::rank())]);
    }
};

//--For ROps, the representation is not evaluated at the index before
//--passing to the UnaryROp operator.  This is often used for more complicated
//--operations where multiple indices of the Rep must be accessed to compute
//--a single value (e.g., matrix-vector multiplication), but we still want to
//--compose an expression.

template <typename REP, typename UnaryROp>
class MatrixUnaryROpRep : public UnaryROp, private StcOp
{
public:
  using value_type = typename UnaryROp::value_type;
  using is_scalar = std::false_type;
  using is_vector = std::false_type;
  using is_rowvec = std::false_type;
  using is_matrix = std::true_type;
  using is_string = std::false_type;
private:
  REP m_rep;
public:
  template <typename... OpArgs>
  HOSTDEVICE constexpr MatrixUnaryROpRep(REP&&       a_rep,
                                         OpArgs&&... a_opArgs) noexcept
    :
    UnaryROp(std::forward<OpArgs>(a_opArgs)...),
    m_rep(std::forward<REP>(a_rep))
    { }
  HOSTDEVICE constexpr MatrixUnaryROpRep(
    const MatrixUnaryROpRep& a_other) noexcept
    :
    UnaryROp(a_other),
    m_rep(a_other.m_rep)
    { }
  HOSTDEVICE constexpr MatrixUnaryROpRep(MatrixUnaryROpRep&& a_other) noexcept
    :
    UnaryROp(std::move(a_other)),
    m_rep(std::forward<REP>(a_other.m_rep))
    { }
  /// Assignment constructor
  MatrixUnaryROpRep& operator=(const MatrixUnaryROpRep&) = delete;
  /// Move assignment constructor
  MatrixUnaryROpRep& operator=(MatrixUnaryROpRep&&) noexcept = delete;
  /// Size (number of elements) of the matrix representation
  HOSTDEVICE static constexpr array_size_type size() noexcept
    {
      return std::decay_t<REP>::size();
    }
  /// Rank of the matrix representation
  HOSTDEVICE static constexpr array_size_type rank() noexcept
    {
      return std::decay_t<REP>::rank();
    }

//--Access

  /// Array indexing
  HOSTDEVICE constexpr auto operator[](
    const array_size_type a_idx) const noexcept
    {
      return this->eval(m_rep, a_idx);
    }
  /// Function/"lambda-like"
  HOSTDEVICE constexpr auto operator()(
    const array_size_type a_idx) const noexcept
    {
      return this->eval(m_rep, a_idx);
    }
  /// Matrix indexing
  HOSTDEVICE constexpr auto operator()(const array_size_type a_i,
                                       const array_size_type a_j) const noexcept
    {
      return this->eval(m_rep,
                        STC_RO(a_i, std::decay_t<REP>::rank()) +
                        STC_CO(a_j, std::decay_t<REP>::rank()));
    }
};


/*******************************************************************************
 *
 * Operators
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*
 * Addition
 *--------------------------------------------------------------------*/

//--Use BinaryOpAdd from StcVector.H

//--Operator between two matrix representations

template <typename REP1, typename REP2,
          std::enable_if_t<
            std::decay_t<REP1>::is_matrix::value
            &&
            std::decay_t<REP2>::is_matrix::value>* = nullptr
          >
HOSTDEVICE constexpr auto
operator+(REP1&& a_rep1, REP2&& a_rep2) noexcept
{
  return stc::MatrixBinaryOpRep<REP1, REP2,
                                stc::BinaryOpAdd<
                                  typename std::decay_t<REP1>::value_type,
                                  typename std::decay_t<REP2>::value_type>>(
    std::forward<REP1>(a_rep1),
    std::forward<REP2>(a_rep2));
}

//--Operator between left scalar and right matrix representation

template <typename T1, typename REP2,
          std::enable_if_t<
            !std::is_base_of<stc::StcOp, std::decay_t<T1> >::value
            &&
            std::decay_t<REP2>::is_matrix::value>* = nullptr
          >
HOSTDEVICE constexpr auto
operator+(const T1& a_scalar, REP2&& a_rep2) noexcept
{
  return stc::MatrixBinaryOpRep<stc::ScalarRep<T1>, REP2,
                                stc::BinaryOpAdd<
                                  T1,
                                  typename std::decay_t<REP2>::value_type>>(
    stc::ScalarRep<T1>(a_scalar),
    std::forward<REP2>(a_rep2));
}

//--Operator between left matrix representation and right scalar

template <typename REP1, typename T2,
          std::enable_if_t<
            std::decay_t<REP1>::is_matrix::value
            &&
            !std::is_base_of<stc::StcOp, std::decay_t<T2> >::value>* = nullptr
          >
HOSTDEVICE constexpr auto
operator+(REP1&& a_rep1, const T2& a_scalar) noexcept
{
  return stc::MatrixBinaryOpRep<REP1, stc::ScalarRep<T2>,
                                stc::BinaryOpAdd<
                                  typename std::decay_t<REP1>::value_type,
                                  T2>>(
    std::forward<REP1>(a_rep1),
    stc::ScalarRep<T2>(a_scalar));
}

/*--------------------------------------------------------------------*
 * Subtraction
 *--------------------------------------------------------------------*/

//--Use BinaryOpSubtract from StcVector.H

//--Operator between two matrix representations

template <typename REP1, typename REP2,
          std::enable_if_t<
            std::decay_t<REP1>::is_matrix::value
            &&
            std::decay_t<REP2>::is_matrix::value>* = nullptr
          >
HOSTDEVICE constexpr auto
operator-(REP1&& a_rep1, REP2&& a_rep2) noexcept
{
  return stc::MatrixBinaryOpRep<REP1, REP2,
                                stc::BinaryOpSubtract<
                                  typename std::decay_t<REP1>::value_type,
                                  typename std::decay_t<REP2>::value_type>>(
    std::forward<REP1>(a_rep1),
    std::forward<REP2>(a_rep2));
}

//--Operator between left scalar and right matrix representation

template <typename T1, typename REP2,
          std::enable_if_t<
            !std::is_base_of<stc::StcOp, std::decay_t<T1> >::value
            &&
            std::decay_t<REP2>::is_matrix::value>* = nullptr
          >
HOSTDEVICE constexpr auto
operator-(const T1& a_scalar, REP2&& a_rep2) noexcept
{
  return stc::MatrixBinaryOpRep<stc::ScalarRep<T1>, REP2,
                                stc::BinaryOpSubtract<
                                  T1,
                                  typename std::decay_t<REP2>::value_type>>(
    stc::ScalarRep<T1>(a_scalar),
    std::forward<REP2>(a_rep2));
}

//--Operator between left matrix representation and right scalar

template <typename REP1, typename T2,
          std::enable_if_t<
            std::decay_t<REP1>::is_matrix::value
            &&
            !std::is_base_of<stc::StcOp, std::decay_t<T2> >::value>* = nullptr
          >
HOSTDEVICE constexpr auto
operator-(REP1&& a_rep1, const T2& a_scalar) noexcept
{
  return stc::MatrixBinaryOpRep<REP1, stc::ScalarRep<T2>,
                                stc::BinaryOpSubtract<
                                  typename std::decay_t<REP1>::value_type,
                                  T2>>(
    std::forward<REP1>(a_rep1),
    stc::ScalarRep<T2>(a_scalar));
}

/*--------------------------------------------------------------------*
 * Multiplication
 * Only supports multiplication by a scalar
 *--------------------------------------------------------------------*/

//--Use BinaryOpMultiply from StcVector.H

//--Operator between left scalar and right matrix representation

template <typename T1, typename REP2,
          std::enable_if_t<
            !std::is_base_of<stc::StcOp, std::decay_t<T1> >::value
            &&
            std::decay_t<REP2>::is_matrix::value>* = nullptr
          >
HOSTDEVICE constexpr auto
operator*(const T1& a_scalar, REP2&& a_rep2) noexcept
{
  return stc::MatrixBinaryOpRep<stc::ScalarRep<T1>, REP2,
                                stc::BinaryOpMultiply<
                                  T1,
                                  typename std::decay_t<REP2>::value_type>>(
    stc::ScalarRep<T1>(a_scalar),
    std::forward<REP2>(a_rep2));
}

//--Operator between left matrix representation and right scalar

template <typename REP1, typename T2,
          std::enable_if_t<
            std::decay_t<REP1>::is_matrix::value
            &&
            !std::is_base_of<stc::StcOp, std::decay_t<T2> >::value>* = nullptr
          >
HOSTDEVICE constexpr auto
operator*(REP1&& a_rep1, const T2& a_scalar) noexcept
{
  return stc::MatrixBinaryOpRep<REP1, stc::ScalarRep<T2>,
                                stc::BinaryOpMultiply<
                                  typename std::decay_t<REP1>::value_type,
                                  T2>>(
    std::forward<REP1>(a_rep1),
    stc::ScalarRep<T2>(a_scalar));
}

/*--------------------------------------------------------------------*
 * Division
 * Only supports division by a scalar
 *--------------------------------------------------------------------*/

//--Use BinaryOpDivision from StcVector.H

//--Operator between left matrix representation and right scalar

template <typename REP1, typename T2,
          std::enable_if_t<
            std::decay_t<REP1>::is_matrix::value
            &&
            !std::is_base_of<stc::StcOp, std::decay_t<T2> >::value>* = nullptr
          >
HOSTDEVICE constexpr auto
operator/(REP1&& a_rep1, const T2& a_scalar) noexcept
{
  return stc::MatrixBinaryOpRep<REP1, stc::ScalarRep<T2>,
                                stc::BinaryOpDivision<
                                  typename std::decay_t<REP1>::value_type,
                                  T2>>(
    std::forward<REP1>(a_rep1),
    stc::ScalarRep<T2>(a_scalar));
}

/*--------------------------------------------------------------------*
 * Transpose
 *--------------------------------------------------------------------*/

template <typename REP>
struct UnaryROpTranspose
{
  using value_type = typename REP::value_type;
  HOSTDEVICE constexpr auto eval(const REP&            a_rep,
                                 const array_size_type a_idx) const noexcept
    {
      constexpr array_size_type R = REP::rank();
      const array_size_type ir = STC_RI(a_idx, R);
      const array_size_type ic = STC_CI(a_idx, R);
      return a_rep[STC_RO(ic, R) + STC_CO(ir, R)];
    }
};

/// Transpose Operator
template <typename REP,
          std::enable_if_t<std::decay_t<REP>::is_matrix::value>* = nullptr
          >
HOSTDEVICE constexpr auto
operator^(REP&& a_rep, stc::TransposeType) noexcept
{
  using F = stc::UnaryROpTranspose<std::decay_t<REP>>;
  return stc::MatrixUnaryROpRep<REP, F>(std::forward<REP>(a_rep));
}

/*--------------------------------------------------------------------*
 * Negation
 *--------------------------------------------------------------------*/

// Use UnaryOpNegation from StcVector.H

/// Negation operator
template <typename REP,
          std::enable_if_t<std::decay_t<REP>::is_matrix::value>* = nullptr
          >
HOSTDEVICE constexpr auto
operator-(REP&& a_rep) noexcept
{
  return stc::MatrixUnaryOpRep<REP,
                               stc::UnaryOpNegation<
                                 typename std::decay_t<REP>::value_type>>(
    std::forward<REP>(a_rep));
}


/*******************************************************************************
 *
 * Other operators
 *
 * These have complexity usually requiring immediate evaluation.
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*
 * Equivalence
 *--------------------------------------------------------------------*/

// Use BinaryOpEquivalence_E from StcVector.H

/// Equivalence operator immediately evaluates and returns
/** The operands must have the same type, size, and values for true equivalence
 */
template <typename REP1, typename REP2,
          std::enable_if_t<
            std::decay_t<REP1>::is_matrix::value
            &&
            std::decay_t<REP2>::is_matrix::value>* = nullptr
          >
HOSTDEVICE constexpr bool
operator==(REP1&& a_rep1, REP2&& a_rep2) noexcept
{
  using T1 = typename std::decay_t<REP1>::value_type;
  using T2 = typename std::decay_t<REP2>::value_type;
  constexpr stc::array_size_type R = std::decay_t<REP1>::rank();
  static_assert(std::is_same<T1, T2>::value && R == std::decay_t<REP2>::rank(),
                "Invalid comparison between different matrix types");
  auto f = stc::MatrixBinaryOpRep<REP1, REP2,
                                  stc::BinaryOpEquivalence<T1, T2>>(
    std::forward<REP1>(a_rep1),
    std::forward<REP2>(a_rep2));
  return stc::andEachElement<R*R>(f);
}

/*--------------------------------------------------------------------*
 * Non-equivalence
 *--------------------------------------------------------------------*/

// Use BinaryOpEquivalence_E from StcVector.H

/// Non-equivalence operator immediately evaluates and returns
/** The operands must have the same type, size, and values for true equivalence
 */
template <typename REP1, typename REP2,
          std::enable_if_t<
            std::decay_t<REP1>::is_matrix::value
            &&
            std::decay_t<REP2>::is_matrix::value>* = nullptr
          >
HOSTDEVICE constexpr bool
operator!=(REP1&& a_rep1, REP2&& a_rep2) noexcept
{
  using T1 = typename std::decay_t<REP1>::value_type;
  using T2 = typename std::decay_t<REP2>::value_type;
  constexpr stc::array_size_type N = std::decay_t<REP1>::size();
  if (!std::is_same<T1, T2>::value ||
      std::decay_t<REP1>::rank() != std::decay_t<REP2>::rank())
    {
      return true;
    }
  const auto f = stc::MatrixBinaryOpRep<REP1, REP2,
                                        stc::BinaryOpEquivalence<T1, T2>>(
    std::forward<REP1>(a_rep1),
    std::forward<REP2>(a_rep2));
  return !stc::andEachElement<N>(f);
}

/*--------------------------------------------------------------------*
 * Matrix-vector multiply
 *--------------------------------------------------------------------*/

// Evaluates a row in matrix-vector multiply.
/* Evaluated immediately
 */
template <typename REP1, typename REP2>
struct BinaryOpMVMul_Row_E
{
  using value_type = decltype(typename REP1::value_type()*
                              typename REP2::value_type());
  HOSTDEVICE constexpr BinaryOpMVMul_Row_E(const REP1&           a_rep1,
                                           const REP2&           a_rep2,
                                           const array_size_type a_ir) noexcept
    :
    m_rep1(a_rep1),
    m_rep2(a_rep2),
    m_ro(STC_RO(a_ir, REP1::rank()))
    { }
  HOSTDEVICE constexpr auto operator()(
    const array_size_type a_idx) const noexcept
    {
      constexpr array_size_type R = REP1::rank(); (void)R;
      return m_rep1[m_ro + STC_CO(a_idx, R)]*m_rep2[a_idx];
    }
  const REP1& m_rep1;
  const REP2& m_rep2;
  const array_size_type m_ro;
};

// Evaluates all in matrix-vector multiply.
/* Evaluated immediately
 */
template <typename REP1, typename REP2>
struct BinaryOpMVMul_E
{
  using value_type = decltype(typename REP1::value_type()*
                              typename REP2::value_type());
  HOSTDEVICE constexpr BinaryOpMVMul_E(const REP1& a_rep1,
                                       const REP2& a_rep2) noexcept
    :
    m_rep1(a_rep1),
    m_rep2(a_rep2)
    { }
  HOSTDEVICE constexpr auto operator()(
    const array_size_type a_idx) const noexcept
    {
      constexpr array_size_type R = std::decay_t<REP1>::rank();
      using F = BinaryOpMVMul_Row_E<REP1, REP2>;
      auto f = F(m_rep1, m_rep2, a_idx);
      return sumEachElement<R>(f);
    }
  const REP1& m_rep1;
  const REP2& m_rep2;
};

// Represents all evaluations in matrix-vector multiply
/* Passed to representation
 */
template <typename REP1, typename REP2>
struct BinaryROpMVMul
{
  using value_type = decltype(typename REP1::value_type{}*
                              typename REP2::value_type{});
  HOSTDEVICE constexpr auto eval(const REP1&           a_rep1,
                                 const REP2&           a_rep2,
                                 const array_size_type a_idx) const noexcept
    {
      constexpr array_size_type R = REP1::rank();
      using F = BinaryOpMVMul_Row_E<REP1, REP2>;
      auto f = F(a_rep1, a_rep2, a_idx);
      return sumEachElement<R>(f);
    }
};

// Next two functions are for dispatch

//--Implementation of matrix-vector multiply
/*  This returns a vector expression since evaluating Rep2 (which must be done
 *  for each row) is trivial.
 *  \param[in]  a_rep1  Matrix
 *  \param[in]  a_rep2  A vector implementation
 */
template <typename REP1, typename REP2>
HOSTDEVICE constexpr auto
mvmulImpl(REP1&& a_rep1, REP2&& a_rep2, std::true_type) noexcept
{
  using F = BinaryROpMVMul<std::decay_t<REP1>, std::decay_t<REP2>>;
  return VectorBinaryROpRep<REP1, REP2, F>(std::forward<REP1>(a_rep1),
                                           std::forward<REP2>(a_rep2));
}

//--Implementation of matrix-vector multiply
/*  Since we evaluate Rep2 in advance, there is no point in returning an
 *  expression.  This returns a vector implementation.
 *  \param[in]  a_rep1  Matrix
 *  \param[in]  a_rep2  A vector expression
 */
template <typename REP1, typename REP2>
HOSTDEVICE constexpr auto
mvmulImpl(REP1&& a_rep1, REP2&& a_rep2, std::false_type) noexcept
{
  constexpr array_size_type R = std::decay_t<REP1>::rank();
  using T = typename std::decay_t<REP1>::value_type;
  Vector<T, R> vec2 = a_rep2;  // Evaluates a_rep2
  using F = BinaryOpMVMul_E<std::decay_t<REP1>, decltype(vec2)>;
  auto f = F(a_rep1, vec2);
  return Vector<T, R>(make_array<T, R>(f));
}

/// Matrix-vector multiply operator
/** \param[in]  a_rep1  Matrix
 *  \param[in]  a_rep2  A vector implementation or expression
 *
 *  Type dispatch is used to call a specific function based on whether REP2
 *  (the vector) is an implementation or an expression.  Depending on the
 *  dispatch, this function might return an expression or vector
 *  implementation.
 */
template <typename REP1, typename REP2,
          // We need more than a default template argument to differentiate
          // the signature of arguments.  Here enable_if is on the LHS.
          std::enable_if_t<
            std::decay_t<REP1>::is_matrix::value
            &&
            std::decay_t<REP2>::is_vector::value>* = nullptr
          >
HOSTDEVICE constexpr auto
operator*(REP1&& a_rep1, REP2&& a_rep2) noexcept
{
  constexpr array_size_type R = std::decay_t<REP1>::rank();
  static_assert(R == std::decay_t<REP2>::size(),
                "Vectors must have same size as matrix rank for matrix-vector "
                "multiply");
  using T = typename std::decay_t<REP1>::value_type;
  static_assert(std::is_same<T, typename std::decay_t<REP2>::value_type>::value,
                "Operands must have same type for matrix-vector multiply");
  return stc::mvmulImpl(
    std::forward<REP1>(a_rep1),
    std::forward<REP2>(a_rep2),
    std::integral_constant<
      bool,
      std::is_base_of<stc::VectorImpl<T, R, DefaultImpl>,
                      std::decay_t<REP2>>::value
      ||
      std::is_base_of<stc::VectorImpl<T, R, stc::VectorConstAliasImpl>,
                      std::decay_t<REP2>>::value
      ||
      std::is_base_of<stc::VectorImpl<T, R, stc::VectorAliasImpl>,
                      std::decay_t<REP2>>::value
    >{});
}

/*--------------------------------------------------------------------*
 * Matrix multiplication
 *--------------------------------------------------------------------*/

// Evaluates the matrix multiplication at a single source index
/* Evaluated immediately
 */
template <typename REP1, typename REP2>
struct BinaryOpMatmul_Rank_E
{
  using value_type = decltype(typename REP1::value_type()*
                              typename REP2::value_type());
  HOSTDEVICE constexpr BinaryOpMatmul_Rank_E(
    const REP1&           a_rep1,
    const REP2&           a_rep2,
    const array_size_type a_ir,
    const array_size_type a_ic) noexcept
    :
    m_rep1(a_rep1),
    m_rep2(a_rep2),
    m_ro(STC_RO(a_ir, REP1::rank())),
    m_co(STC_CO(a_ic, REP1::rank()))
    { }
  HOSTDEVICE constexpr auto operator()(
    const array_size_type a_idx) const noexcept
    {
      constexpr array_size_type R = REP1::rank();
      return m_rep1[m_ro + STC_CO(a_idx, R)]*m_rep2[STC_RO(a_idx, R) + m_co];
    }
  const REP1& m_rep1;
  const REP2& m_rep2;
  const array_size_type m_ro;
  const array_size_type m_co;
};

// Evaluates the matrix multiplication at a given destination index
/* Evaluated immediately
 */
template <typename REP1, typename REP2>
struct BinaryOpMatmul_E
{
  using value_type = decltype(typename REP1::value_type()*
                              typename REP2::value_type());
  HOSTDEVICE constexpr BinaryOpMatmul_E(const REP1& a_rep1,
                                        const REP2& a_rep2) noexcept
    :
    m_rep1(a_rep1),
    m_rep2(a_rep2)
    { }
  HOSTDEVICE constexpr auto operator()(
    const array_size_type a_ir,
    const array_size_type a_ic) const noexcept
    {
      constexpr array_size_type R = std::decay_t<REP1>::rank();
      using F = BinaryOpMatmul_Rank_E<REP1, REP2>;
      auto f = F(m_rep1, m_rep2, a_ir, a_ic);
      return sumEachElement<R>(f);
    }
  const REP1& m_rep1;
  const REP2& m_rep2;
};

// Next two functions are for dispatch

//--Implementation of matrix multiply from matrix implementations
/*  There is no need to evaluate the matrices since the element values are not
 *  expressions.
 *  \param[in]  a_rep1  Matrix implementation
 *  \param[in]  a_rep2  Matrix implementation
 */
template <typename REP1, typename REP2>
HOSTDEVICE constexpr auto
matmulImpl(REP1&& a_rep1, REP2&& a_rep2, std::true_type) noexcept
{
  constexpr array_size_type R = std::decay_t<REP1>::rank();
  using T = typename std::decay_t<REP1>::value_type;
  using F = BinaryOpMatmul_E<std::decay_t<REP1>, std::decay_t<REP2>>;
  auto f = F(a_rep1, a_rep2);
  return Matrix<T, R>(make_array2<T, R>(f));
}

//--Implementation of matrix multiply from one or more matrix expressions
/*  At least one of the operands is an expression and needs to be evaluated
 *  into a (temporary) implementation first.
 *  \param[in]  a_rep1  Matrix
 *  \param[in]  a_rep2  Matrix
 */
template <typename REP1, typename REP2>
HOSTDEVICE constexpr auto
matmulImpl(REP1&& a_rep1, REP2&& a_rep2, std::false_type) noexcept
{
  constexpr array_size_type R = std::decay_t<REP1>::rank();
  using T = typename std::decay_t<REP1>::value_type;
  Matrix<T, R> mat1 = a_rep1;  // Evaluates a_rep1
  Matrix<T, R> mat2 = a_rep2;  // Evaluates a_rep2
  using F = BinaryOpMatmul_E<decltype(mat1), decltype(mat2)>;
  auto f = F(mat1, mat2);
  return Matrix<T, R>(make_array2<T, R>(f));
}

/// Matrix-matrix multiply operator
/** \param[in]  a_rep1  Matrix
 *  \param[in]  a_rep2  Matrix
 *
 *  Type dispatch is used to call an optimized function if both REP1 and REP2
 *  are matrix implementations (versus expressions).  Otherwise the matrix
 *  representations are evaluated since each index is accessed multiple times.
 */
template <typename REP1, typename REP2,
          // We need more than a default template argument to differentiate
          // the signature of arguments.  Here enable_if is on the LHS.
          std::enable_if_t<
            std::decay_t<REP1>::is_matrix::value
            &&
            std::decay_t<REP2>::is_matrix::value>* = nullptr
          >
HOSTDEVICE constexpr auto
operator*(REP1&& a_rep1, REP2&& a_rep2) noexcept
{
  constexpr array_size_type R = std::decay_t<REP1>::rank();
  static_assert(R == std::decay_t<REP2>::rank(),
                "Matrices must have same rank for matrix multiply");
  using T = typename std::decay_t<REP1>::value_type;
  static_assert(std::is_same<T, typename std::decay_t<REP2>::value_type>::value,
                "Operands must have same type for matrix multiply");
  return stc::matmulImpl(
    std::forward<REP1>(a_rep1),
    std::forward<REP2>(a_rep2),
    std::integral_constant<
      bool,
      (
        std::is_base_of<stc::MatrixImpl<T, R, DefaultImpl>,
                        std::decay_t<REP1>>::value
        ||
        std::is_base_of<stc::MatrixImpl<T, R, stc::MatrixConstAliasImpl>,
                        std::decay_t<REP1>>::value
        ||
        std::is_base_of<stc::MatrixImpl<T, R, stc::MatrixAliasImpl>,
                        std::decay_t<REP1>>::value
      )
      &&
      (
        std::is_base_of<stc::MatrixImpl<T, R, DefaultImpl>,
                        std::decay_t<REP2>>::value
        ||
        std::is_base_of<stc::MatrixImpl<T, R, stc::MatrixConstAliasImpl>,
                        std::decay_t<REP2>>::value
        ||
        std::is_base_of<stc::MatrixImpl<T, R, stc::MatrixAliasImpl>,
                        std::decay_t<REP2>>::value
      )
    >{});
}

/*--------------------------------------------------------------------*
 * Outer product of vectors
 *--------------------------------------------------------------------*/

//  Evaluates the outer product at a given destination index
/*  Evaluated immediately
 */
template <typename REP1, typename REP2>
struct BinaryOpOuter_E
{
  using value_type = decltype(typename REP1::value_type()*
                              typename REP2::value_type());
  HOSTDEVICE constexpr BinaryOpOuter_E(const REP1& a_rep1,
                                       const REP2& a_rep2) noexcept
    :
    m_rep1(a_rep1),
    m_rep2(a_rep2)
    { }
  HOSTDEVICE constexpr auto operator()(
    const array_size_type a_ir,
    const array_size_type a_ic) const noexcept
    {
      return m_rep1[a_ir]*m_rep2[a_ic];
    }
  const REP1& m_rep1;
  const REP2& m_rep2;
};

// Next two functions are for dispatch

//--Implementation of outer product from vector implementations
/*  There is no need to evaluate the vectors since the element values are not
 *  expressions.
 *  \param[in]  a_rep1  Vector implementation
 *  \param[in]  a_rep2  Vector implementation
 */
template <typename REP1, typename REP2>
HOSTDEVICE constexpr auto
outerImpl(REP1&& a_rep1, REP2&& a_rep2, std::true_type) noexcept
{
  constexpr array_size_type R = std::decay_t<REP1>::size();
  using T1 = typename std::decay_t<REP1>::value_type;
  using T2 = typename std::decay_t<REP2>::value_type;
  using T = decltype(T1()*T2());
  using F = BinaryOpOuter_E<std::decay_t<REP1>, std::decay_t<REP2>>;
  auto f = F(a_rep1, a_rep2);
  return Matrix<T, R>(make_array2<T, R>(f));
}

//--Implementation of outer product from one or more vector expressions
/*  At least one of the operands is an expression and needs to be evaluated
 *  into a (temporary) implementation first.
 *  \param[in]  a_rep1  Vector implementation
 *  \param[in]  a_rep2  Vector implementation
 */
template <typename REP1, typename REP2>
HOSTDEVICE constexpr auto
outerImpl(REP1&& a_rep1, REP2&& a_rep2, std::false_type) noexcept
{
  constexpr array_size_type R = std::decay_t<REP1>::size();
  using T1 = typename std::decay_t<REP1>::value_type;
  using T2 = typename std::decay_t<REP2>::value_type;
  using T = decltype(T1()*T2());
  Vector<T, R> vec1 = a_rep1;  // Evaluates a_rep1
  Vector<T, R> vec2 = a_rep2;  // Evaluates a_rep2
  using F = BinaryOpOuter_E<decltype(vec1), decltype(vec2)>;
  auto f = F(vec1, vec2);
  return Matrix<T, R>(make_array2<T, R>(f));
}

/// Outer product operator
/** \param[in]  a_rep1  Vector implementation
 *  \param[in]  a_rep2  Vector implementation
 *
 *  Type dispatch is used to call an optimized function if both REP1 and REP2
 *  are vector implementations (versus expressions).  Otherwise the vector
 *  representations are evaluated since each index is accessed multiple times.
 *  This operator immediately evaluates the matrix implementation.
 */
template <typename REP1, typename REP2,
          std::enable_if_t<
            std::decay_t<REP1>::is_vector::value
            &&
            std::decay_t<REP2>::is_vector::value>* = nullptr
          >
HOSTDEVICE constexpr auto
outer(REP1&& a_rep1, REP2&& a_rep2) noexcept
{
  constexpr array_size_type I = std::decay_t<REP1>::size();
  static_assert(I == std::decay_t<REP2>::size(),
                "Vectors must have same size for outer product");
  using T1 = typename std::decay_t<REP1>::value_type;
  using T2 = typename std::decay_t<REP2>::value_type;
  using T = decltype(T1()*T2());
  return stc::outerImpl(
    std::forward<REP1>(a_rep1),
    std::forward<REP2>(a_rep2),
    std::integral_constant<
      bool,
      (
        std::is_base_of<stc::VectorImpl<T, I, stc::DefaultImpl>,
                        std::decay_t<REP1>>::value
        ||
        std::is_base_of<stc::VectorImpl<T, I, stc::VectorConstAliasImpl>,
                        std::decay_t<REP1>>::value
        ||
        std::is_base_of<stc::VectorImpl<T, I, stc::VectorAliasImpl>,
                        std::decay_t<REP1>>::value
      )
      &&
      (
        std::is_base_of<stc::VectorImpl<T, I, stc::DefaultImpl>,
                        std::decay_t<REP2>>::value
        ||
        std::is_base_of<stc::VectorImpl<T, I, stc::VectorConstAliasImpl>,
                        std::decay_t<REP2>>::value
        ||
        std::is_base_of<stc::VectorImpl<T, I, stc::VectorAliasImpl>,
                        std::decay_t<REP2>>::value
      )
    >{});
}

/// Alternative method to invoke outer product
/** Operator* between column and row vector (e.g., m = u*(u^t))
 */
template <typename REP1, typename REP2,
          // We need more than a default template argument to differentiate
          // the signature of arguments.  Here enable_if is on the LHS.
          std::enable_if_t<
            std::decay_t<REP1>::is_vector::value
            &&
            std::decay_t<REP2>::is_rowvec::value>* = nullptr
          >
HOSTDEVICE constexpr auto
operator*(REP1&& a_rep1, REP2&& a_rep2) noexcept
{
  return outer(std::forward<REP1>(a_rep1), a_rep2.get());
}

/*--------------------------------------------------------------------*
 * Plane projection tensor
 *--------------------------------------------------------------------*/

template <typename REP,
          std::enable_if_t<std::decay_t<REP>::is_vector::value>* = nullptr
          >
HOSTDEVICE constexpr auto
planeProjectionTensor(REP&& a_rep) noexcept
{
  constexpr array_size_type I = std::decay_t<REP>::size();
  using T = typename std::decay_t<REP>::value_type;
  return stc::make_matrix<T, I>::identity() -
    outer(std::forward<REP>(a_rep), std::forward<REP>(a_rep));
}


/*******************************************************************************
 *
 * Class MatrixImpl: external related functions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Output the matrix
/** \param[in] a_os     The output stream
 *  \param[in] a_mat    Matrix to output
 *  \return             The output stream
 *//*-----------------------------------------------------------------*/

template <typename T,
          array_size_type R,
          template<typename, stc::array_size_type> class Impl>
std::ostream &operator<<(std::ostream&                      a_os,
                         const stc::MatrixImpl<T, R, Impl>& a_mat)
{
  a_os << "\n";
  for (int i = 0; i != R; ++i)
    {
      a_os << std::setw(4) << a_mat(i, 0);
      for (int j = 1; j != R; ++j)
        {
          a_os << " " << std::setw(4) << a_mat(i, j);
        }
      a_os << '\n';
    }
  return a_os;
}

}  // namespace stc

#endif  /* ! defined _STCMATRIX_H_ */
