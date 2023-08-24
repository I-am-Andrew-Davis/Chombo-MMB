/*
  TODO
  - comments
*/
/*
  Notes
  (1) unsigned is used instead of std::size_t in case the latter creeps into
      32-bit GPU code.
  (2) Strides are stored instead of dimension sizes because it seemed to work
      better on GPUs.
  (3) Storing the UB is an option.  One can either store the UB or compute it
      from the LB and the dimension size.  But if strides are stored, computing
      the size involves a division.  Mostly this affects the performance of
      range checking.  In GCC, it appears better to compute the UB.  In clang++
      it appears better to store the UB.  On GPUs it seems to make little
      difference either way but range checking is very expensive.
  (4) Class organization
      ArrayBase --------------+-->ArrayTally
      {                       |
        BareImpl---------+    +-->ArrayBracket-->array
        {                |
          T*   m_data    +-->BoundsImpl-------->RangeCheckImpl
          RVec m_stride      {
        };                     RVec m_LB
      };                       Rvec m_UB (Opt)
                             };
 */
#ifndef _SHAPEARRAY_H_
#define _SHAPEARRAY_H_

/* Forces all predefined implementations to be range checked */
#ifndef NDEBUG
#define SHAPE_ARRAY_DEBUG
#endif

/* If errors, use throw instead of assert.  Note that assert is used if file is
   processed with nvcc, ignoring this define */
// #define SHAPE_ARRAY_THROW

/* Do not check for runtime errors */
// #define SHAPE_ARRAY_NDEBUG

/* Store the upper bound instead of calculating it (recommended for clang++) */
// #define SHAPE_ARRAY_STOUB

#include <utility>
#include <type_traits>

/*
  Notes:
    - tag __CUDACC__ means the file is processed with nvcc.  The class can
      be modified whether on the cpu or gpu.
    - tag __CUDA_ARCH__ is defined for device code.  This only works in
      functions with the __device__ qualifier.  It cannot be used in the
      class declaration.
*/

#undef HOSTDEVICE
#undef DEVICE
#ifdef __CUDACC__
#define HOSTDEVICE __host__ __device__
#define DEVICE __device__
#else
#define HOSTDEVICE
#define DEVICE
#endif

#ifdef SHAPE_ARRAY_NDEBUG
  #ifdef SHAPE_ARRAY_DEBUG
    #error "Both SHAPE_ARRAY_DEBUG and SHAPE_ARRAY_NDEBUG cannot be defined!"
  #endif
  /* Applied to check for a valid dimension */
  #define SHAPE_ARRAY_DIMCHECK(cond) ((void)0)
  /* Applied to check that array has been defined */
  #define SHAPE_ARRAY_DEFCHECK(cond) ((void)0)
  /* Applied to check that an index is in the range of a dimension */
  #define SHAPE_ARRAY_RNGCHECK(cond) ((void)0)
#else
  #if defined(SHAPE_ARRAY_THROW) && !defined(__CUDAAC__)
    #include <stdexcept>
    #define SHAPE_ARRAY_DIMCHECK(cond) \
      if (!(cond)) throw std::out_of_range("invalid dimension")
    #define SHAPE_ARRAY_DEFCHECK(cond) \
      if (!(cond)) throw std::logic_error("array not defined")
    #define SHAPE_ARRAY_RNGCHECK(cond) \
      if (!(cond)) throw std::out_of_range("index outside dimension range")
  #else
    #include <cassert>
    #define SHAPE_ARRAY_DIMCHECK(cond) assert(cond)
    #define SHAPE_ARRAY_DEFCHECK(cond) assert(cond)
    #define SHAPE_ARRAY_RNGCHECK(cond) assert(cond)
  #endif
#endif

#include "BaseNamespaceHeader.H"

namespace shape
{

using row_ordered        = std::true_type;   ///< Row-ordered tag
using column_ordered     = std::false_type;  ///< Column-ordered tag
using full_size          = std::integral_constant<unsigned, 0>;
using min_size           = std::integral_constant<unsigned, 1>;
using bare_config        = std::integral_constant<unsigned, 2>;
using bounds_config      = std::integral_constant<unsigned, 3>;
using range_check_config = std::integral_constant<unsigned, 4>;


/*******************************************************************************
 *
 * An aggregate for storing IIx[R].  Used instead of static_array to support
 * CUDA.
 *
 ******************************************************************************/

template <unsigned R, typename IIx>
struct RVec
{
  // This is supposed to be an aggregate and use aggregate initialization so no
  // constructors are allowed
  using index_type = IIx;
  HOSTDEVICE constexpr IIx operator[](const IIx a_idx) const noexcept
    { return m_RVec[a_idx]; }
  HOSTDEVICE constexpr IIx& operator[](const IIx a_idx) noexcept
    { return m_RVec[a_idx]; }
  HOSTDEVICE constexpr IIx operator()(const IIx a_idx) const noexcept
    { return m_RVec[a_idx]; }
  HOSTDEVICE constexpr IIx& operator()(const IIx a_idx) noexcept
    { return m_RVec[a_idx]; }
  IIx m_RVec[R];
};

// Zero-size specialization
template <typename IIx>
struct RVec<0, IIx>
{
  using index_type = IIx;
  HOSTDEVICE constexpr IIx operator[](const IIx a_idx) const noexcept
    { return *static_cast<IIx*>(nullptr); }
  HOSTDEVICE constexpr IIx& operator[](const IIx a_idx) noexcept
    { return *static_cast<IIx*>(nullptr); }
  HOSTDEVICE constexpr IIx operator()(const IIx a_idx) const noexcept
    { return *static_cast<IIx*>(nullptr); }
  HOSTDEVICE constexpr IIx& operator()(const IIx a_idx) noexcept
    { return *static_cast<IIx*>(nullptr); }
};


/*******************************************************************************
 *
 * Metaprogramming for unrolling loops
 *
 ******************************************************************************/

//--Standard unrolls

template <unsigned N, typename F>
struct UnrollLoop
{
  /// Apply f to each element
  HOSTDEVICE static constexpr void eval(F&& f)
    {
      UnrollLoop<N-1, F>::eval(std::forward<F>(f));
      f(N-1);
    }
  /// Apply f to each element in reverse order
  HOSTDEVICE static constexpr void reval(F&& f)
    {
      f(N-1);
      UnrollLoop<N-1, F>::reval(std::forward<F>(f));
    }
  /// Sum result of f
  HOSTDEVICE static constexpr auto sum(F&& f)
    {
      return UnrollLoop<N-1, F>::sum(std::forward<F>(f)) + f(N-1);
    }
  /// Product result of f
  HOSTDEVICE static constexpr auto product(F&& f)
    {
      return UnrollLoop<N-1, F>::product(std::forward<F>(f))*f(N-1);
    }
};

template <typename F>
struct UnrollLoop<1, F>
{
  HOSTDEVICE static constexpr void eval(F&& f)
    {
      f(0);
    }
  HOSTDEVICE static constexpr void reval(F&& f)
    {
      f(0);
    }
  HOSTDEVICE static constexpr auto sum(F&& f)
    {
      return f(0);
    }
  HOSTDEVICE static constexpr auto product(F&& f)
    {
      return f(0);
    }
};

// In case N = 0, don't do anything
template <typename F>
struct UnrollLoop<0, F>
{
  HOSTDEVICE static constexpr auto eval(F&& f)
    { }
  HOSTDEVICE static constexpr auto reval(F&& f)
    { }
  HOSTDEVICE static constexpr auto sum(F&& f)
    {
      return static_cast<decltype(f(0))>(0);
    }
  HOSTDEVICE static constexpr auto product(F&& f)
    {
      return static_cast<decltype(f(0))>(1);
    }
};

//--Linear indices

template <unsigned N, typename Strd, typename Idxs>
struct LinearIndex
{
  using IIx = typename std::decay_t<Idxs>::index_type;
  /// Index math from each element
  HOSTDEVICE static constexpr IIx sum(Strd&& strd,
                                      Idxs&& idxs)
    {
      return idxs(N-1)*strd(N-2) +
        LinearIndex<N-1, Strd, Idxs>::sum(std::forward<Strd>(strd),
                                          std::forward<Idxs>(idxs));
    }
};

template <typename Strd, typename Idxs>
struct LinearIndex<1, Strd, Idxs>
{
  using IIx = typename std::decay_t<Idxs>::index_type;
  HOSTDEVICE static constexpr IIx sum(Strd&& strd, 
                                      Idxs&& idxs)
    {
      return idxs(0);
    }
};

// In case N = 0, don't do anything
template <typename Strd, typename Idxs>
struct LinearIndex<0, Strd, Idxs>
{
  using IIx = typename std::decay_t<Idxs>::index_type;
  HOSTDEVICE static constexpr IIx sum(Strd&& strd,
                                      Idxs&& idxs)
    {
      return (IIx)0;
    }
};

//--Linear indices with range checking.  The unit stride index is not included
//--in the sum but is range checked

template <unsigned N, typename Impl, typename Idxs>
struct LinearIndexFrom1RC
{
  using IIx = typename std::decay_t<Idxs>::index_type;
  /// Index math from each element
  HOSTDEVICE static constexpr IIx sum(const Impl& impl, Idxs&& idxs)
    {
      impl.checkRange(N-1, idxs(N-1));
      return idxs(N-1)*impl.m_stride(N-2) +
        LinearIndexFrom1RC<N-1, Impl, Idxs>::sum(impl,
                                                 std::forward<Idxs>(idxs));
    }
};

// In case N = 2, terminate cascade
template <typename Impl, typename Idxs>
struct LinearIndexFrom1RC<2, Impl, Idxs>
{
  using IIx = typename std::decay_t<Idxs>::index_type;
  HOSTDEVICE static constexpr IIx sum(const Impl& impl, Idxs&& idxs)
    {
      impl.checkRange(1, idxs(1));
      return idxs(1)*impl.m_stride(0);
    }
};

// In case N = 1, only range check
template <typename Impl, typename Idxs>
struct LinearIndexFrom1RC<1, Impl, Idxs>
{
  using IIx = typename std::decay_t<Idxs>::index_type;
  HOSTDEVICE static constexpr IIx sum(const Impl& impl, Idxs&& idxs)
    {
      // impl.checkRange0(idxs(0));
      return (IIx)0;
    }
};

// In case N = 0, don't do anything
template <typename Impl, typename Idxs>
struct LinearIndexFrom1RC<0, Impl, Idxs>
{
  using IIx = typename std::decay_t<Idxs>::index_type;
  HOSTDEVICE static constexpr IIx sum(const Impl& impl, Idxs&& idxs)
    {
      return (IIx)0;
    }
};

//--Functions for deducing types other than N

/// Apply function object 'f' to each element 
template <unsigned N, typename F>
HOSTDEVICE constexpr void
forEachElement(F&& f)
{
  UnrollLoop<N, F>::eval(std::forward<F>(f));
}

/// Apply function object 'f' to each element passing index as template
template <unsigned N, typename F>
HOSTDEVICE constexpr void
forEachElementT(F&& f)
{
  UnrollLoop<N, F>::teval(std::forward<F>(f));
}

/// Apply function object 'f' to each element in reverse order
template <unsigned N, typename F>
HOSTDEVICE constexpr void
forEachReverseElement(F&& f)
{
  UnrollLoop<N, F>::reval(std::forward<F>(f));
}

/// Apply function object 'f' to each element and sum results
template <unsigned N, typename F>
HOSTDEVICE constexpr auto
sumEachElement(F&& f)
{
  return UnrollLoop<N, F>::sum(std::forward<F>(f));
}

/// Apply function object 'f' to each element and multiply results
template <unsigned N, typename F>
HOSTDEVICE constexpr auto
prodEachElement(F&& f)
{
  return UnrollLoop<N, F>::product(std::forward<F>(f));
}

/// Linear index from strides and indices
template <unsigned N, typename Strd, typename Idxs>
HOSTDEVICE constexpr auto
linearIndexEachElement(Strd&& a_strd, Idxs&& a_idxs)
{
  return LinearIndex<N, Strd, Idxs>::sum(std::forward<Strd>(a_strd),
                                         std::forward<Idxs>(a_idxs));
}

/// Linear index from implementation and indices with range checking
/** Linear indices with range checking.  The unit stride index is not included
 *  in the sum but is range checked
 */
template <unsigned N, typename Impl, typename Idxs>
HOSTDEVICE constexpr auto
linearIndexEachElementFrom1RC(const Impl& a_impl, Idxs&& a_idxs)
{
  return LinearIndexFrom1RC<N, Impl, Idxs>::sum(
    a_impl, std::forward<Idxs>(a_idxs));
}


/*******************************************************************************
 *
 * Makers for RVec
 *
 * A bit of wizardry to initialize a static array.  The approach was taken from
 * the way make_index_sequence itself is implemented in gcc-6.1.0.
 *
 ******************************************************************************/

// The default parameter here generates the integer sequence
template <unsigned R, typename IIx, typename Op,
          typename ISeq = std::make_integer_sequence<unsigned, R>>
struct Make_RVec;

// The sequence has type std::index_sequence<Is...> and now this specialization
// matches
template <unsigned R, typename IIx, typename Op, unsigned... Is>
struct Make_RVec<R, IIx, Op, std::integer_sequence<unsigned, Is...>>
{
  HOSTDEVICE static constexpr RVec<R, IIx> with(Op&& a_op) noexcept
    {
      return { std::forward<Op>(a_op)(Is)... };  // op is applied to each
                                                 // integer in the sequence
    }
};

// Create a function so Op can be deduced
template <unsigned R, typename IIx, typename Op>
HOSTDEVICE constexpr auto
make_RVec(Op&& a_op) noexcept
{
  return Make_RVec<R, IIx, Op>::with(std::forward<Op>(a_op));
}


/*******************************************************************************
 */
/// Returns an index or size from list of indices
/**
 *  This is where row-ordered versus column-ordered happens
 *  \tparam IIx        Indexing type
 *  \tparam Ord        Row or column ordering
 *  \tparam N          Number of indices
 *  
 *//*+*************************************************************************/

template <typename IIx, typename Ord, unsigned N>
struct GetIdx;

template <typename IIx, unsigned N>
struct GetIdx<IIx, column_ordered, N>
{
  using index_type = IIx;
  using order_type = column_ordered;
  template <typename... Idxs>
  HOSTDEVICE constexpr GetIdx(const Idxs... a_idxs) noexcept
    :
    m_idxs({ (IIx)a_idxs... })
    { }
  HOSTDEVICE constexpr IIx operator()(const IIx a_idx) const noexcept
    {
      return m_idxs[a_idx];
    }
private:
  const RVec<N, IIx> m_idxs;
};

template <typename IIx, unsigned N>
struct GetIdx<IIx, row_ordered, N>
{
  using index_type = IIx;
  using order_type = row_ordered;
  template <typename... Idxs>
  HOSTDEVICE constexpr GetIdx(const Idxs... a_idxs) noexcept
    :
    m_idxs({ (IIx)a_idxs... })
    { }
  HOSTDEVICE constexpr IIx operator()(const IIx a_idx) const noexcept
    {
      return m_idxs[N-1 - a_idx];
    }
private:
  const RVec<N, IIx> m_idxs;
};


/*******************************************************************************
 */
/// Returns an index or size from dimension D of a vector-like object
/**
 *  This is where row-ordered versus column-ordered happens
 *  \tparam Vec        Vector-like type
 *  \tparam IIx        Indexing type
 *  \tparam Ord        Row or column ordering
 *  \tparam N          Dimension of vector
 *  
 *//*+*************************************************************************/

template <typename Vec, typename IIx, typename Ord, unsigned N>
struct GetIdxVec;

template <typename Vec, typename IIx, unsigned N>
struct GetIdxVec<Vec, IIx, column_ordered, N>
{
  using index_type = IIx;
  using order_type = column_ordered;
  HOSTDEVICE constexpr GetIdxVec(const Vec& a_vec) noexcept
    :
    m_vec(a_vec)
    { }
  HOSTDEVICE constexpr IIx operator()(const IIx a_idx) const noexcept
    {
      return m_vec[a_idx];
    }
private:
  const Vec& m_vec;
};

template <typename Vec, typename IIx, unsigned N>
struct GetIdxVec<Vec, IIx, row_ordered, N>
{
  using index_type = IIx;
  using order_type = row_ordered;
  HOSTDEVICE constexpr GetIdxVec(const Vec& a_vec) noexcept
    :
    m_vec(a_vec)
    { }
  HOSTDEVICE constexpr IIx operator()(const IIx a_idx) const noexcept
    {
      return m_vec[N-1 - a_idx];
    }
private:
  const Vec& m_vec;
};


/*******************************************************************************
 */
/// Data implementations
/**
 *//*+*************************************************************************/

// Operator for zeroing all indices
template <typename IIx>
struct UniformIdx
{
  HOSTDEVICE constexpr UniformIdx(const IIx a_val = (IIx)0) noexcept
    :
    m_val(a_val)
    { }
  HOSTDEVICE constexpr IIx operator()(const IIx a_idx) const noexcept
    {
      return m_val;
    }
  const IIx m_val;
};

// Operator for stride
template <typename Dims>
struct StrideIdx
{
  using IIx = typename std::decay_t<Dims>::index_type;
  HOSTDEVICE constexpr StrideIdx(Dims&& a_dims) noexcept
    :
    m_dims(std::forward<Dims>(a_dims)),
    m_cache((IIx)1)
    { }
  HOSTDEVICE constexpr IIx operator()(const IIx a_idx) noexcept
    {
      m_cache *= m_dims(a_idx);
      return m_cache;
    }
  Dims&& m_dims;
  IIx m_cache;
};

// Operator to make UB from dimensions (assumes LB = 0)
template <typename Dims>
struct UBIdxDim
{
  using IIx = typename std::decay_t<Dims>::index_type;
  HOSTDEVICE constexpr UBIdxDim(Dims&& a_dims) noexcept
    :
    m_dims(std::forward<Dims>(a_dims))
    { }
  HOSTDEVICE constexpr IIx operator()(const IIx a_idx) noexcept
    {
      return m_dims(a_idx) - (IIx)1;
    }
  Dims&& m_dims;
};

/// Implementation with compile-time dimensions
/** Use this to save space, there is little to no benefit in terms of
 *  performance
 */
template <typename T, unsigned R, typename Ord, typename IIx, IIx... Dims>
struct CBareImpl
{
private:
  using Sz = full_size;
  using BoundsDefSz = void;
  using RangeCheckDefSz = void;
  template <typename Impl>
  friend class ArrayBase;

public:
  using value_type = T;
  using index_type = IIx;
  using order_type = Ord;
  using stosize_type = Sz;
  static constexpr unsigned rank = R;
  static constexpr RVec<R-Sz::value, IIx> m_stride =
    make_RVec<R-Sz::value, IIx>(
      StrideIdx<GetIdx<IIx, Ord, sizeof...(Dims)>>(
        GetIdx<IIx, Ord, sizeof...(Dims)>{ Dims... }));

  HOSTDEVICE CBareImpl() noexcept
    :
    m_data(nullptr)
    { }

  template <typename... NullDims>
  HOSTDEVICE CBareImpl(T *const a_data, const NullDims... a_dims) noexcept
    :
    m_data(a_data)
    { }

  CBareImpl(const CBareImpl&) = default;
  CBareImpl(CBareImpl&&) = default;
  CBareImpl& operator=(const CBareImpl&) = default;
  CBareImpl& operator=(CBareImpl&&) = default;
  ~CBareImpl() = default;

  template <typename... NullDims>
  HOSTDEVICE void define(T *const a_data, const NullDims... a_dims) noexcept
    {
      m_data = a_data;
    }

  HOSTDEVICE void clear() noexcept
    {
      m_data = nullptr;
    }

  template <typename... Idxs>
  HOSTDEVICE void setAbsoluteLBIndices(const Idxs... a_absLB) noexcept
    { }

  template <typename... Idxs>
  HOSTDEVICE void setRelativeLBIndices(const Idxs... a_relLB) noexcept
    { }

  HOSTDEVICE IIx size(const unsigned a_dim) const noexcept
    {
      const IIx lower = (a_dim == 0) ? (IIx)1 : m_stride[a_dim - 1];
      return m_stride[a_dim]/lower;
    }

  HOSTDEVICE constexpr void checkRange(const unsigned a_dim,
                                       const IIx      a_idx) const noexcept
    { }

  HOSTDEVICE constexpr void checkRange0(const IIx a_idx) const noexcept
    { }

  T* m_data;
};

//--Definitions (concerns about ODR use)

template <typename T, unsigned R, typename Ord, typename IIx, IIx... Dims>
constexpr RVec<R-full_size::value, IIx>
CBareImpl<T, R, Ord, IIx, Dims...>::m_stride;

/*
  As a general policy, when an RVec type is constructed or assigned or intended
  to support constexpr, make_RVec is used.  Otherwise a lambda is used with a
  metaprogramming unroll such as forEachElement.
 */
template <typename T, unsigned R, typename Ord, typename IIx, typename Sz>
struct BareImpl
{
private:
  using BoundsDefSz = void;
  using RangeCheckDefSz = void;
  template <typename Impl>
  friend class ArrayBase;

public:
  using value_type = T;
  using index_type = IIx;
  using order_type = Ord;
  using stosize_type = Sz;
  static constexpr unsigned rank = R;

  HOSTDEVICE BareImpl() noexcept
    :
    m_data(nullptr),
    m_stride(make_RVec<R-Sz::value, IIx>(UniformIdx<IIx>((IIx)0)))
    { }

  template <typename... Dims>
  HOSTDEVICE BareImpl(T *const a_data, const Dims... a_dims) noexcept
    :
    m_data(a_data),
    m_stride(make_RVec<R-Sz::value, IIx>(
               StrideIdx<GetIdx<IIx, Ord, sizeof...(Dims)>>(
                 GetIdx<IIx, Ord, sizeof...(Dims)>(a_dims...))))
    { }

  BareImpl(const BareImpl&) = default;
  BareImpl(BareImpl&&) = default;
  BareImpl& operator=(const BareImpl&) = default;
  BareImpl& operator=(BareImpl&&) = default;
  ~BareImpl() = default;

  template <typename... Dims>
  HOSTDEVICE void define(T *const a_data, const Dims... a_dims) noexcept
    {
      m_data = a_data;
      m_stride = make_RVec<R-Sz::value, IIx>(
        StrideIdx<GetIdx<IIx, Ord, sizeof...(Dims)>>(
          GetIdx<IIx, Ord, sizeof...(Dims)>(a_dims...)));
    }

  HOSTDEVICE void clear() noexcept
    {
      m_data = nullptr;
      m_stride = make_RVec<R-Sz::value, IIx>(UniformIdx<IIx>((IIx)0));
    }

  template <typename... Idxs>
  HOSTDEVICE void setAbsoluteLBIndices(const Idxs... a_absLB) noexcept
    { }

  template <typename... Idxs>
  HOSTDEVICE void setRelativeLBIndices(const Idxs... a_relLB) noexcept
    { }

  HOSTDEVICE IIx size(const unsigned a_dim) const noexcept
    {
      const IIx lower = (a_dim == 0) ? (IIx)1 : m_stride[a_dim - 1];
      return m_stride[a_dim]/lower;
    }

  HOSTDEVICE constexpr void checkRange(const unsigned a_dim,
                                       const IIx      a_idx) const noexcept
    { }

  HOSTDEVICE constexpr void checkRange0(const IIx a_idx) const noexcept
    { }

  T* m_data;
  RVec<R-Sz::value, IIx> m_stride;    ///< Stride for each higher dimension.
                                      ///< Unit stride 1 is not stored.
                                      ///< m_stride[0] is size of unit stride
                                      ///< dimension.  Whether or not full size
                                      ///< of array is stored in m_stride[R-1]
                                      ///< depends on Sz=full_size or
                                      ///< Sz=min_size
};

template <typename T, unsigned R, typename Ord, typename IIx, typename Sz>
struct BoundsImpl : public BareImpl<T, R, Ord, IIx, full_size>
{
private:
  using Bare = BareImpl<T, R, Ord, IIx, full_size>;
  // Note that Sz is replaced by full_size, this is the size requested by the
  // caller.  To construct a class with the same type as BoundsImpl, derived
  // classes need to know this.
  using BoundsDefSz = Sz;
  using RangeCheckDefSz = void;
  template <typename Impl>
  friend class ArrayBase;

public:
  using value_type = T;
  using index_type = IIx;
  using order_type = Ord;
  using stosize_type = full_size;
  static constexpr unsigned rank = R;

  HOSTDEVICE BoundsImpl() noexcept
    :
    Bare(),
    m_LB(make_RVec<R, IIx>(UniformIdx<IIx>((IIx) 0)))
#ifdef SHAPE_ARRAY_STOUB
    ,
    m_UB(make_RVec<R, IIx>(UniformIdx<IIx>((IIx)-1)))
#endif
    { }

  template <typename... Dims>
  HOSTDEVICE BoundsImpl(T *const a_data, const Dims... a_dims) noexcept
    :
    Bare(a_data, a_dims...),
    m_LB(make_RVec<R, IIx>(UniformIdx<IIx>((IIx) 0)))
#ifdef SHAPE_ARRAY_STOUB
    ,
    m_UB(make_RVec<R, IIx>(UBIdxDim<GetIdx<IIx, Ord, sizeof...(Dims)>>(
                             GetIdx<IIx, Ord, sizeof...(Dims)>(a_dims...))))
#endif
    { }

  BoundsImpl(const BoundsImpl&) = default;
  BoundsImpl(BoundsImpl&&) = default;
  BoundsImpl& operator=(const BoundsImpl&) = default;
  BoundsImpl& operator=(BoundsImpl&&) = default;
  ~BoundsImpl() = default;

  template <typename... Dims>
  HOSTDEVICE void define(T *const a_data, const Dims... a_dims) noexcept
    {
      Bare::define(a_data, a_dims...);
      m_LB = make_RVec<R, IIx>(UniformIdx<IIx>((IIx)0));
#ifdef SHAPE_ARRAY_STOUB
      m_UB = make_RVec<R, IIx>(
        UBIdxDim<GetIdx<IIx, Ord, sizeof...(Dims)>>(
          GetIdx<IIx, Ord, sizeof...(Dims)>(a_dims...)));
#endif
    }

  HOSTDEVICE void clear() noexcept
    {
      Bare::clear();
      m_LB = make_RVec<R, IIx>(UniformIdx<IIx>((IIx) 0));
#ifdef SHAPE_ARRAY_STOUB
      m_UB = make_RVec<R, IIx>(UniformIdx<IIx>((IIx)-1));
#endif
    }

  template <typename... Idxs>
  HOSTDEVICE void setAbsoluteLBIndices(const Idxs... a_absLB) noexcept
    {
      static_assert(sizeof...(Idxs) >= R, "Insufficient indices");
      SHAPE_ARRAY_DEFCHECK(Bare::m_stride[0] > 0);
      GetIdx<IIx, Ord, sizeof...(Idxs)> absLB(a_absLB...);
      forEachElement<R>([this, &absLB]
                        (const unsigned a_dim)
                          {
#ifdef SHAPE_ARRAY_STOUB
                            const IIx sizeDim = size(a_dim);
#endif
                            m_LB(a_dim) = absLB(a_dim);
#ifdef SHAPE_ARRAY_STOUB
                            m_UB(a_dim) = absLB(a_dim) + sizeDim - (IIx)1;
#endif
                          });
    }

  template <typename... Idxs>
  HOSTDEVICE void setRelativeLBIndices(const Idxs... a_relLB) noexcept
    {
      static_assert(sizeof...(Idxs) >= R, "Insufficient indices");
      SHAPE_ARRAY_DEFCHECK(Bare::m_stride[0] > 0);
      GetIdx<IIx, Ord, sizeof...(Idxs)> relLB(a_relLB...);
      forEachElement<R>([this, &relLB]
                        (const unsigned a_dim)
                          {
                            m_LB(a_dim) += relLB(a_dim);
#ifdef SHAPE_ARRAY_STOUB
                            m_UB(a_dim) += relLB(a_dim);
#endif
                          });
    }

#ifdef SHAPE_ARRAY_STOUB
  // Potentially cheaper size if UB is stored
  HOSTDEVICE IIx size(const unsigned a_dim) const noexcept
    {
      return m_UB(a_dim) - m_LB(a_dim) + (IIx)1;
    }
#endif

#ifndef SHAPE_ARRAY_STOUB
  // Used if UB is NOT stored
  HOSTDEVICE constexpr IIx m_UB(const unsigned a_dim) const noexcept
    {
      return m_LB(a_dim) + Bare::size(a_dim) - 1;
    }
#endif

  RVec<R, IIx> m_LB;        ///< Lower bounds of the array
#ifdef SHAPE_ARRAY_STOUB
  RVec<R, IIx> m_UB;        ///< Upper bounds of the array
#endif
};

template <typename T, unsigned R, typename Ord, typename IIx, typename Sz>
struct RangeCheckImpl : public BoundsImpl<T, R, Ord, IIx, full_size>
{
private:
  using Bounds = BoundsImpl<T, R, Ord, IIx, full_size>;
  using Bare   = BareImpl<T, R, Ord, IIx, full_size>;
  // Note that Sz is replaced by full_size, BoundsDefSz is the size requested by
  // the caller.  To construct a class with the same type as BoundsImpl or
  // RangeCheckImpl, derived classes need to know this.
  using BoundsDefSz = full_size;
  using RangeCheckDefSz = Sz;
  template <typename Impl>
  friend class ArrayBase;

public:
  using value_type = T;
  using index_type = IIx;
  using order_type = Ord;
  using stosize_type = full_size;
  static constexpr unsigned rank = R;

  HOSTDEVICE RangeCheckImpl() noexcept
    :
    Bounds()
    { }

  template <typename... Dims>
  HOSTDEVICE RangeCheckImpl(T *const a_data, const Dims... a_dims) noexcept
    :
    Bounds(a_data, a_dims...)
    { }

  RangeCheckImpl(const RangeCheckImpl&)            = default;
  RangeCheckImpl(RangeCheckImpl&&)                 = default;
  RangeCheckImpl& operator=(const RangeCheckImpl&) = default;
  RangeCheckImpl& operator=(RangeCheckImpl&&)      = default;
  ~RangeCheckImpl()                                = default;

  template <typename... Dims>
  HOSTDEVICE void define(T *const a_data, const Dims... a_dims) noexcept
    {
      Bounds::define(a_data, a_dims...);
    }

  HOSTDEVICE void clear() noexcept
    {
      Bounds::clear();
    }

  HOSTDEVICE constexpr void checkRange(const unsigned a_dim,
                                       const IIx      a_idx) const
#ifndef SHAPE_ARRAY_THROW
    noexcept
#endif
    {
      SHAPE_ARRAY_RNGCHECK(a_idx >= Bounds::m_LB(a_dim) &&
                           a_idx <= Bounds::m_UB(a_dim));
    }

  HOSTDEVICE constexpr void checkRange0(const IIx a_idx) const
#ifndef SHAPE_ARRAY_THROW
    noexcept
#endif
    {
      SHAPE_ARRAY_RNGCHECK(a_idx >= Bounds::m_LB(0) &&
                           a_idx <= Bounds::m_LB(0) + Bare::m_stride(0) - 1);
    }
};


/*******************************************************************************
 */
/// Base array has all operators except operator[]
/**
 *//*+*************************************************************************/

template <typename Impl>
class ArrayBase
{
private:
  using Impl_t = std::decay_t<Impl>;
  using T = typename Impl_t::value_type;
  using IIx = typename Impl_t::index_type;
  using Ord = typename Impl_t::order_type;
  using Sz = typename Impl_t::stosize_type;

public:
  using value_type = typename Impl_t::value_type;
  using index_type = typename Impl_t::index_type;
  using order_type = typename Impl_t::order_type;
  using stosize_type = typename Impl_t::stosize_type;
  static constexpr unsigned rank = Impl_t::rank;

  /// Default constructor (undefined array)
  HOSTDEVICE ArrayBase() noexcept
    :
    m_impl()
    { }

  /// Construct without dimensions
  /** Usually this implies use of CBareImpl
   *  \param[in]  a_data  Linear memory to reshaped
   */
  HOSTDEVICE ArrayBase(T *const a_data) noexcept
    :
    m_impl(a_data)
    {
    }

  /// Construct from integer sequence
  /** \tparam Dims...     Integer sequence of array dimensions
   *  \param[in]  a_data  Linear memory to reshaped
   *  \param[in]  a_dim0  Size of first dimension
   *  \param[in]  a_dims  Size of each remaining dimension
   */
  template <typename D0, typename... Dims,
            std::enable_if_t<std::is_integral<D0>::value, int> = 0>
  HOSTDEVICE ArrayBase(T *const a_data,
                       const D0 a_dim0, const Dims... a_dims) noexcept
    :
    m_impl(a_data, a_dim0, a_dims...)
    {
      static_assert(sizeof...(a_dims) + 1 >= rank, "Insufficient dimensions");
    }

  /// Construct from dimensions given in a vector-like object
  template <typename Vec,
            std::enable_if_t<std::is_pointer<Vec>::value ||
                             std::is_array<Vec>::value ||
                             std::is_class<Vec>::value, int> = 0>
  HOSTDEVICE ArrayBase(T *const a_data, const Vec& a_dims) noexcept
    :
    ArrayBase(a_data, a_dims, std::make_integer_sequence<unsigned, rank>{})
    {
    }

  /// Construct with a given reference to an implementation
  HOSTDEVICE ArrayBase(Impl&& a_impl) noexcept
    :
    m_impl(std::forward<Impl>(a_impl))
    { }

  /// Copy constructor (synthesized)
  ArrayBase(const ArrayBase&) = default;
  /// Move constructor (synthesized)
  ArrayBase(ArrayBase&&) = default;
  // Copy assignment (synthesized)
  ArrayBase& operator=(const ArrayBase&) = default;
  // Move assignment (synthesized)
  ArrayBase& operator=(ArrayBase&&) = default;
  // Destructor (sythesized)
  ~ArrayBase() = default;

  /// Define the array with pointer and integer sequence
  /** Old definition is discarded
   *  \tparam Dims...     Integer sequence of array dimensions
   *  \param[in]  a_data  Linear memory to reshaped
   *  \param[in]  a_dims  Size of each dimension
   */
  template <typename... Dims,
            std::enable_if_t<
              std::is_integral<
                typename std::tuple_element<0, std::tuple<Dims...>>::type
                >::value, int> = 0>
  HOSTDEVICE void define(T *const a_data, const Dims... a_dims) noexcept
    {
      static_assert(sizeof...(a_dims) >= rank, "Insufficient dimensions");
      m_impl.define(a_data, a_dims...);
    }

  /// Define the array with pointer and vector-like object
  /** Old definition is discarded
   *  \tparam Vec         A type with operator[]
   *  \param[in]  a_data  Linear memory to reshaped
   *  \param[in]  a_vec   Size of each dimension
   */
  template <typename Vec,
            std::enable_if_t<std::is_pointer<Vec>::value ||
                             std::is_array<Vec>::value ||
                             std::is_class<Vec>::value, int> = 0>
  HOSTDEVICE void define(T *const a_data, const Vec& a_dims) noexcept
    {
      define(a_data, a_dims, std::make_integer_sequence<unsigned, rank>{});
    }

  /// Clear the array (appears as undefined)
  HOSTDEVICE void clear() noexcept
    {
      m_impl.clear();
    }

//--Capabilities

  /// Are the sizes of all dimension known
  HOSTDEVICE constexpr bool is_full_size() const noexcept
    { return std::is_same<Sz, full_size>::value; }
  /// Is the input row-ordered?
  HOSTDEVICE constexpr bool is_row_ordered() const noexcept
    { return std::is_same<Ord, row_ordered>::value; }
  /// Is the input column-ordered?
  HOSTDEVICE constexpr bool is_column_ordered() const noexcept
    { return std::is_same<Ord, column_ordered>::value; }
  /// Are the bounds of the array stored?
  HOSTDEVICE constexpr bool has_bounds() const noexcept
    { 
      return std::is_base_of<
        BoundsImpl<T, rank, Ord, IIx, typename Impl_t::BoundsDefSz>, 
        Impl_t>::value;
    }
  /// Are indices range-checked against the bounds?
  HOSTDEVICE constexpr bool is_range_checked() const noexcept
    {
      return std::is_same<
        RangeCheckImpl<T, rank, Ord, IIx, typename Impl_t::RangeCheckDefSz>,
        Impl_t>::value;
    }

  /// Set a new absolute lower bound and reset the data pointer accordingly
  /** \tparam Idxs...     Integer sequence of lower bound indices
   *  \param[in]  a_absLB Absolute lower bound indices
   *  \note
   *  <ul>
   *    <li> This function is only available if the storage contains bounds.
   *  </ul>
   */
  template <typename... Idxs,
            typename _Impl_t = Impl_t,
            std::enable_if_t<std::is_base_of<BoundsImpl<T, rank, Ord, IIx, Sz>,
                                             _Impl_t>::value &&
                             std::is_integral<
                               typename std::tuple_element<
                                 0, std::tuple<Idxs...>>::type>::value,
                             int> = 0>
  HOSTDEVICE void setLB(const Idxs... a_absLB) noexcept
    {
      static_assert(sizeof...(Idxs) >= rank, "Insufficient indices");
      resetData(data()
                + linearIndexEachElement<rank>(m_impl.m_stride, m_impl.m_LB)
                - linearIndex(a_absLB...));
      m_impl.setAbsoluteLBIndices(a_absLB...);
    }

  /// Set a new absolute lower bound and reset the data pointer accordingly
  /** \tparam Vec         Vector-like object containing lower bound indices
   *  \param[in]  a_absLB Vector of absolute lower bound indices.  MUST have
   *                      size >= rank!
   *  \note
   *  <ul>
   *    <li> This function is only available if the storage contains bounds.
   *    <li> WARNING: There is no ability to check the size of the vector.
   *         Incorrect usage can lead to segmentation faults or worse.
   *  </ul>
   */
  template <typename Vec,
            typename _Impl_t = Impl_t,
            std::enable_if_t<std::is_base_of<BoundsImpl<T, rank, Ord, IIx, Sz>,
                                             _Impl_t>::value &&
                             (std::is_pointer<Vec>::value ||
                              std::is_array<Vec>::value ||
                              std::is_class<Vec>::value), int> = 0>
  HOSTDEVICE void setLB(const Vec& a_absLB) noexcept
    {
      setLB(a_absLB, std::make_integer_sequence<unsigned, rank>{});
    }

  /// Set a new relative lower bound and reset the data pointer accordingly
  /** \tparam Idxs...     Integer sequence of lower bound indices
   *  \param[in]  a_relLB Relative lower bound indices
   *  \note
   *  <ul>
   *    <li> This sets the lower bound relative to the current.  E.g, calling
   *         this with a_relLB = (0...) has no effect.
   *  </ul>
   */
  template <typename... Idxs,
            std::enable_if_t<
              std::is_integral<
                typename std::tuple_element<0, std::tuple<Idxs...>>::type
                >::value, int> = 0>
  HOSTDEVICE void setRelativeLB(const Idxs... a_relLB) noexcept
    {
      static_assert(sizeof...(Idxs) >= rank, "Insufficient indices");
      resetData(data() - linearIndex(a_relLB...));
      m_impl.setRelativeLBIndices(a_relLB...);
    }

  /// Set a new relative lower bound and reset the data pointer accordingly
  /** \tparam Vec         Vector-like object containing lower bound indices
   *  \param[in]  a_relLB Vector of relative lower bound indices.  MUST have
   *                      size >= rank!
   *  \note
   *  <ul>
   *    <li> This sets the lower bound relative to the current.  E.g, calling
   *         this with a_relLB = Vec(0) has no effect.
   *    <li> WARNING: There is no ability to check the size of the vector.
   *         Incorrect usage can lead to segmentation faults or worse.
   *  </ul>
   */
  template <typename Vec,
            std::enable_if_t<std::is_pointer<Vec>::value ||
                             std::is_array<Vec>::value ||
                             std::is_class<Vec>::value, int> = 0>
  HOSTDEVICE void setRelativeLB(const Vec& a_relLB) noexcept
    {
      setRelativeLB(a_relLB, std::make_integer_sequence<unsigned, rank>{});
    }

  /// Set new absolute lower bounds without modifying the data pointer
  /** \tparam Idxs...     Integer sequence of lower bound indices
   *  \param[in]  a_absLB Absolute lower bound indices
   *  \note
   *  <ul>
   *    <li> WARNING: This only sets the recorded lower bounds and does not
   *         adjust the data pointer.  Use setLB to also adjust the data
   *         pointer.  An example use-case is if the user manually resets the
   *         data pointer but wishes to use range-checking of indices.  The
   *         lower-bound indices then need to be properly set.
   *    <li> WARNING: This function is always available but only has effect
   *         if the storage contains bounds.
   *  </ul>
   */
  template <typename... Idxs,
            std::enable_if_t<
              std::is_integral<
                typename std::tuple_element<0, std::tuple<Idxs...>>::type
                >::value, int> = 0>
  HOSTDEVICE void setLBIndices(const Idxs... a_absLB) noexcept
    {
      static_assert(sizeof...(Idxs) >= rank, "Insufficient indices");
      m_impl.setAbsoluteLBIndices(a_absLB...);
    }

  /// Set new absolute lower bounds without modifying the data pointer
  /** \tparam Vec         Vector-like object containing lower bound indices
   *  \param[in]  a_absLB Vector of absolute lower bound indices.  MUST have
   *                      size >= rank!
   *  \note
   *  <ul>
   *    <li> WARNING: This only sets the recorded lower bounds and does not
   *         adjust the data pointer.  Use setLB to also adjust the data
   *         pointer.  An example use-case is if the user manually resets the
   *         data pointer but wishes to use range-checking of indices.  The
   *         lower-bound indices then need to be properly set.
   *    <li> WARNING: This function is always available but only has effect
   *         if the storage contains bounds.
   *    <li> WARNING: There is no ability to check the size of the vector.
   *         Incorrect usage can lead to segmentation faults or worse.
   *  </ul>
   */
  template <typename Vec,
            std::enable_if_t<std::is_pointer<Vec>::value ||
                             std::is_array<Vec>::value ||
                             std::is_class<Vec>::value, int> = 0>
  HOSTDEVICE void setLBIndices(const Vec& a_absLB) noexcept
    {
      setLBIndices(a_absLB, std::make_integer_sequence<unsigned, rank>{});
    }

  /// Set new relative lower bounds without modifying the data pointer
  /** \tparam Idxs...     Integer sequence of lower bound indices
   *  \param[in]  a_relLB Relative lower bound indices
   *  \note
   *  <ul>
   *    <li> WARNING: This only sets the recorded lower bounds and does not
   *         adjust the data pointer.  Use setRelativeLB to also adjust the data
   *         pointer.  An example use-case is if the user manually resets the
   *         data pointer but wishes to also use range-checking of indices.  The
   *         lower-bound indices then need to be properly set.
   *    <li> WARNING: This function is always available but only has effect
   *         if the storage contains bounds.
   *  </ul>
   */
  template <typename... Idxs,
            std::enable_if_t<
              std::is_integral<
                typename std::tuple_element<0, std::tuple<Idxs...>>::type
                >::value, int> = 0>
  HOSTDEVICE void setRelativeLBIndices(const Idxs... a_relLB) noexcept
    {
      static_assert(sizeof...(Idxs) >= rank, "Insufficient indices");
      m_impl.setRelativeLBIndices(a_relLB...);
    }

  /// Set new relative lower bounds without modifying the data pointer
  /** \tparam Vec         Vector-like object containing lower bound indices
   *  \param[in]  a_relLB Vector of relative lower bound indices.  MUST have
   *                      size >= rank!
   *  \note
   *  <ul>
   *    <li> WARNING: This only sets the recorded lower bounds and does not
   *         adjust the data pointer.  Use setRelativeLB to also adjust the data
   *         pointer.  An example use-case is if the user manually resets the
   *         data pointer but wishes to also use range-checking of indices.  The
   *         lower-bound indices then need to be properly set.
   *    <li> WARNING: This function is always available but only has effect
   *         if the storage contains bounds.
   *    <li> WARNING: There is no ability to check the size of the vector.
   *         Incorrect usage can lead to segmentation faults or worse.
   *  </ul>
   */
  template <typename Vec,
            std::enable_if_t<std::is_pointer<Vec>::value ||
                             std::is_array<Vec>::value ||
                             std::is_class<Vec>::value, int> = 0>
  HOSTDEVICE void setRelativeLBIndices(const Vec& a_relLB) noexcept
    {
      setRelativeLBIndices(a_relLB,
                           std::make_integer_sequence<unsigned, rank>{});
    }

  /// Reset only the data pointer
  /** This is often used if the bounds change but array shape remains the same
   *  \param[in]  a_data  New pointer to linear data
   *  \note
   *  <ul>
   *    <li> Use functions setLB (requires storage contains bounds) or 
   *         setRelativeLB to adjust the data pointer based on a set of
   *         lower-bound arguments if you do not want to calculate the
   *         offset manually.
   *    <li> WARNING: If storage contains bounds and/or range-checking, and the
   *         bounds of the array are not separately updated using setLBIndices
   *         or set RelativeLBIndices, expect problems.  It is safer to use
   *         setLB or setRelativeLB to change the lower bound.
   *  </ul>
   */
  HOSTDEVICE void resetData(T *const a_data) noexcept
    {
      m_impl.m_data = a_data;
    }

  /// Access the const data using operator()
  /** \tparam Idxs...     Integer sequence of indices
   *  \param[in]  a_idxs  Indices of element to access.
   *  \return             Const data at the element
   */
  template <typename... Idxs,
            std::enable_if_t<std::is_integral<
                               typename std::tuple_element<
                                 0, std::tuple<Idxs...>>::type>::value,
                             int> = 0>
  HOSTDEVICE const T& operator()(const Idxs... a_idxs) const noexcept
    {
      static_assert(sizeof...(a_idxs) >= rank, "Insufficient indices");
      // This looks excessive but saves 1 instruction on the GPU versus
      // computing the full linear index
      T* const data = m_impl.m_data + linearIndexFrom1RangeCheck(a_idxs...);
      const IIx i0 = GetIdx<IIx, Ord, sizeof...(Idxs)>(a_idxs...)(0);
      m_impl.checkRange0(i0);
      return data[i0];
    }

  /// Access the const data using operator() with vector-like object
  /** \tparam Vec         Vector-like type
   *  \param[in]  a_vec   Vector containing indices.  This must have a size
   *                      >= rank if column ordered and exactly = rank if
   *                      row-ordered.
   *  \return             Const data at the element
   *  \note
   *  <ul>
   *    <li> WARNING: There is no ability to check the size of the vector.
   *         Incorrect usage can lead to segmentation faults or worse.
   *  </ul>
   */
  template <typename Vec,
            std::enable_if_t<std::is_pointer<Vec>::value ||
                             std::is_array<Vec>::value ||
                             std::is_class<Vec>::value, int> = 0>
  HOSTDEVICE const T& operator()(const Vec& a_vec) const noexcept
    {
      T* const data = m_impl.m_data + linearIndexFrom1RangeCheckVec(a_vec);
      const IIx i0 = GetIdxVec<Vec, IIx, Ord, rank>(a_vec)(0);
      m_impl.checkRange0(i0);
      return data[i0];
    }

  /// Access the modifiable data using operator()
  /** \tparam Idxs...     Integer sequence of indices
   *  \param[in]  a_idxs  Indices of element to access.
   *  \return             Modifiable data at the element
   */
  template <typename... Idxs,
            std::enable_if_t<std::is_integral<
                               typename std::tuple_element<
                                 0, std::tuple<Idxs...>>::type>::value,
                             int> = 0>
  HOSTDEVICE T& operator()(const Idxs... a_idxs) noexcept
    {
      static_assert(sizeof...(a_idxs) >= rank, "Insufficient indices");
      // This looks excessive but saves 1 instruction on the GPU versus
      // computing the full linear index
      T* const data = m_impl.m_data + linearIndexFrom1RangeCheck(a_idxs...);
      const IIx i0 = GetIdx<IIx, Ord, sizeof...(Idxs)>(a_idxs...)(0);
      m_impl.checkRange0(i0);
      return data[i0];
    }

  /// Access the modifiable data using operator() with vector-like object
  /** \tparam Vec         Vector-like type
   *  \param[in]  a_vec   Vector containing indices.  This must have a size
   *                      >= rank if column ordered and exactly = rank if
   *                      row-ordered.
   *  \return             Modifiable data at the element
   *  \note
   *  <ul>
   *    <li> WARNING: There is no ability to check the size of the vector.
   *         Incorrect usage can lead to segmentation faults or worse.
   *  </ul>
   */
  template <typename Vec,
            std::enable_if_t<std::is_pointer<Vec>::value ||
                             std::is_array<Vec>::value ||
                             std::is_class<Vec>::value, int> = 0>
  HOSTDEVICE T& operator()(const Vec& a_vec) noexcept
    {
      T* const data = m_impl.m_data + linearIndexFrom1RangeCheckVec(a_vec);
      const IIx i0 = GetIdxVec<Vec, IIx, Ord, rank>(a_vec)(0);
      m_impl.checkRange0(i0);
      return data[i0];
    }

  /// Obtain an index to linear data (from a zero-based lower bound)
  template <typename... Idxs>
  HOSTDEVICE constexpr IIx linearIndex(const Idxs... a_idxs) const noexcept
    {
      static_assert(sizeof...(a_idxs) >= rank, "Insufficient indices");
      return linearIndexEachElement<rank>(
        m_impl.m_stride,
        GetIdx<IIx, Ord, sizeof...(Idxs)>(a_idxs...));
    }

  /// Obtain an index to linear data with a vector-like object (zero-based LB)
  template <typename Vec,
            std::enable_if_t<
              std::is_pointer<Vec>::value ||
              std::is_array<Vec>::value ||
              std::is_class<Vec>::value, int> = 0>
  HOSTDEVICE constexpr IIx linearIndex(const Vec& a_vec) const noexcept
    {
      return linearIndexEachElement<rank>(
        m_impl.m_stride,
        GetIdxVec<Vec, IIx, Ord, rank>(a_vec));
    }

  /// Return the total size of the array
  /** \return             Total number of elements in the array
   *  \note
   *  <ul>
   *    <li> This function is only available if the storage contains all sizes.
   *    <li> This function is fast as the total size is stored directly
   *  </ul>
   */
  template <typename _Sz = Sz,
            std::enable_if_t<std::is_same<_Sz, full_size>::value, int> = 0>
  // If you expect this function can be used, but the compiler says there is
  // no matching function and template argument deduction failed for this
  // candidate, it probably means you are not using storage containing sizes
  // for all dimensions and the enable_if above is disabling the function.
  HOSTDEVICE constexpr IIx size() const noexcept
    {
      return m_impl.m_stride[rank - 1];
    }

  /// Return the size of a dimension
  /** \param[in]  a_dim   Dimension to query size (0 = unit stride dimension)
   *  \return             Number of elements in dimension a_dim
   *  \note
   *  <ul>
   *    <li> This size for the dimension with largest stride is only available
   *         if the storage contains all sizes.
   *    <li> This function is slow as it requires division and should be
   *         avoided where expensive operators are avoided (only if
   *         SHAPE_ARRAY_STOUB is not defined.  However, be aware that
   *         defining SHAPE_ARRAY_STOUB increases the size of this object.
   *         See struct BoundsImpl for more information).
   *  </ul>
   */
  HOSTDEVICE constexpr IIx size(const unsigned a_dim) const
#ifndef SHAPE_ARRAY_THROW
    noexcept
#endif
    {
      SHAPE_ARRAY_DIMCHECK(a_dim < (rank-Sz::value));
      return m_impl.size(a_dim);
    }

  /// Return the stride of a dimension
  /** \param[in]  a_dim   Dimension to query size (0 = unit stride dimension)
   *  \return             Stride for dimension a_dim
   */
  HOSTDEVICE constexpr IIx stride(const unsigned a_dim) const
#ifndef SHAPE_ARRAY_THROW
    noexcept
#endif
    {
      SHAPE_ARRAY_DIMCHECK(a_dim < rank);
      return (a_dim == 0) ? (IIx)1 : m_impl.m_stride[a_dim - 1];
    }

  /// Query a lower bound
  /** \param[in]  a_dim   Dimension to query
   *  \return             Lower bound
   */
  template <typename _Impl_t = Impl_t,
            std::enable_if_t<
              std::is_base_of<
                BoundsImpl<T, rank, Ord, IIx, typename Impl_t::BoundsDefSz>, 
                _Impl_t>::value,
              int> = 0>
  // If you expect this function can be used, but the compiler says there is
  // no matching function and template argument deduction failed for this
  // candidate, it probably means you are not using storage containing bounds
  // and the enable_if above is disabling the function.
  HOSTDEVICE constexpr IIx lowerBound(const unsigned a_dim) const
#ifndef SHAPE_ARRAY_THROW
    noexcept
#endif
    {
      SHAPE_ARRAY_DIMCHECK(a_dim < rank);
      return m_impl.m_LB(a_dim);
    }

  /// Query an upper bound
  /** \param[in]  a_dim   Dimension to query
   *  \return             Upper bound
   *  \note
   *  <ul>
   *    <li> This function is slow as it requires division and should be
   *         avoided where expensive operators are avoided (only if
   *         SHAPE_ARRAY_STOUB is not defined.  However, be aware that
   *         defining SHAPE_ARRAY_STOUB increases the size of this object.
   *         See struct BoundsImpl for more information).
   *  </ul>
   */
  template <typename _Impl_t = Impl_t,
            std::enable_if_t<
              std::is_base_of<
                BoundsImpl<T, rank, Ord, IIx, typename Impl_t::BoundsDefSz>, 
                _Impl_t>::value,
              int> = 0>
  // If you expect this function can be used, but the compiler says there is
  // no matching function and template argument deduction failed for this
  // candidate, it probably means you are not using storage containing bounds
  // and the enable_if above is disabling the function.
  HOSTDEVICE constexpr IIx upperBound(const unsigned a_dim) const
#ifndef SHAPE_ARRAY_THROW
    noexcept
#endif
    {
      SHAPE_ARRAY_DIMCHECK(a_dim < rank);
      return m_impl.m_UB(a_dim);
    }

  /// Check that the range of an index is valid
  /** \param[in]  a_dim   Dimesion of the index
   *  \param[in]  a_idx   Index to range check
   *  There is no return value.  Possible side effects are assertion failure
   *  or throwing, depending on whether SHAPE_ARRAY_THROW is defined.  If
   *  storage does not contain range-checking or if SHAPE_ARRAY_NDEBUG is
   *  defined, there are no side-effects and this is a null function.
   */
  HOSTDEVICE constexpr void checkRange(const unsigned a_dim,
                                       const IIx      a_idx) const
#ifndef SHAPE_ARRAY_THROW
    noexcept
#endif
    {
      m_impl.checkRange(a_dim, a_idx);
    }

  /// Get pointer to the const linear data with offset for a LB
  /**
   *  \return             Pointer to the const data
   *  \note
   *  <ul>
   *    <li> WARNING: This will include an offset from the original pointer if
   *         a non-zero lower bound was set!
   *  </ul>
   */
  HOSTDEVICE const T* data() const noexcept
    {
      return m_impl.m_data;
    }


  /// Get pointer to the modifiable linear data with offset for a LB
  /** \return             Pointer to the modifiable data
   *  \note
   *  <ul>
   *    <li> WARNING: This will include an offset from the original pointer if
   *         a non-zero lower bound was set!
   *  </ul>
   */
  HOSTDEVICE T* data() noexcept
    {
      return m_impl.m_data;
    }

protected:

  /// Helper constructor that extracts dimensions from a vector-like object
  template <typename Vec, unsigned... Is>
  HOSTDEVICE ArrayBase(T *const   a_data,
                       const Vec& a_dims,
                       std::integer_sequence<unsigned, Is...>) noexcept
    :
    m_impl(a_data, a_dims[Is]...)
    {
    }

  /// Helper define that extracts dimensions from a vector-like object
  template <typename Vec, unsigned... Is>
  HOSTDEVICE void define(T *const   a_data,
                         const Vec& a_dims,
                         std::integer_sequence<unsigned, Is...>) noexcept
    {
      m_impl.define(a_data, a_dims[Is]...);
    }

  /// Helper that extracts absolute lower bound from a vector-like object
  template <typename Vec, unsigned... Is>
  HOSTDEVICE void setLB(const Vec& a_absLB,
                        std::integer_sequence<unsigned, Is...>) noexcept
    {
      setLB(a_absLB[Is]...);
    }

  /// Helper that extracts relative lower bound from a vector-like object
  template <typename Vec, unsigned... Is>
  HOSTDEVICE void setRelativeLB(const Vec& a_relLB,
                                std::integer_sequence<unsigned, Is...>) noexcept
    {
      setRelativeLB(a_relLB[Is]...);
    }

  /// Helper that extracts absolute lower bound from a vector-like object
  template <typename Vec, unsigned... Is>
  HOSTDEVICE void setLBIndices(const Vec& a_absLB,
                               std::integer_sequence<unsigned, Is...>) noexcept
    {
      m_impl.setAbsoluteLBIndices(a_absLB[Is]...);
    }

  /// Helper that extracts relative lower bound from a vector-like object
  template <typename Vec, unsigned... Is>
  HOSTDEVICE void setRelativeLBIndices(
    const Vec& a_relLB,
    std::integer_sequence<unsigned, Is...>) noexcept
    {
      m_impl.setRelativeLBIndices(a_relLB[Is]...);
    }

  /// Helper to find a linear index with range-checking exluding unit stride dim
  template <typename... Idxs>
  HOSTDEVICE constexpr IIx linearIndexFrom1RangeCheck(
    const Idxs... a_idxs) const
#ifndef SHAPE_ARRAY_THROW
    noexcept
#endif
    {
      static_assert(sizeof...(a_idxs) >= rank, "Insufficient indices");
      return linearIndexEachElementFrom1RC<rank, Impl_t>(
        m_impl,
        GetIdx<IIx, Ord, sizeof...(Idxs)>(a_idxs...));
    }

  /// Helper to find a linear index with range-checking exluding unit stride dim
  /** For vector-like argument
   */
  template <typename Vec>
  HOSTDEVICE constexpr IIx linearIndexFrom1RangeCheckVec(const Vec& a_vec) const
#ifndef SHAPE_ARRAY_THROW
    noexcept
#endif
    {
      return linearIndexEachElementFrom1RC<rank, Impl_t>(
        m_impl,
        GetIdxVec<Vec, IIx, Ord, rank>(a_vec));
    }

  Impl m_impl;                        ///< Implementation of data for the array
};


/*******************************************************************************
 */
/// Tally impementation
/**
 *//*+*************************************************************************/

template <typename Impl, unsigned D>
class ArrayTally : public ArrayBase<Impl>
{
private:
  using Impl_t = std::decay_t<Impl>;
  using T = typename Impl_t::value_type;
  using IIx = typename Impl_t::index_type;
  using Ord = typename Impl_t::order_type;

public:
  using value_type = typename Impl_t::value_type;
  using index_type = typename Impl_t::index_type;
  using order_type = typename Impl_t::order_type;
  static constexpr unsigned rank = Impl_t::rank;

  // Note: In a cascade, internal data such as a_impl is always an lvalue
  // reference
  HOSTDEVICE ArrayTally(Impl&& a_impl, T* const a_data) noexcept
    :
    ArrayBase<Impl>(std::forward<Impl>(a_impl)),
    m_data(a_data)
    { }

  ArrayTally(const ArrayTally&)            = default;
  ArrayTally(ArrayTally&&)                 = default;
  ArrayTally& operator=(const ArrayTally&) = default;
  ArrayTally& operator=(ArrayTally&&)      = default;
  ~ArrayTally()                            = default;

  /// Operator[]
  HOSTDEVICE const auto operator[](const IIx a_idx) const noexcept
    {
      ArrayBase<Impl>::checkRange(D-1, a_idx);
      return ArrayTally<const Impl&, D-1>(
        ArrayBase<Impl>::m_impl,
        m_data + a_idx*ArrayBase<Impl>::m_impl.m_stride[D-2]);
    }

  /// Operator[]
  HOSTDEVICE auto operator[](const IIx a_idx) noexcept
    {
      ArrayBase<Impl>::checkRange(D-1, a_idx);
      return ArrayTally<Impl&, D-1>(
        ArrayBase<Impl>::m_impl,
        m_data + a_idx*ArrayBase<Impl>::m_impl.m_stride[D-2]);
    }

  /// Linear Index (revised)
  template <typename... Idxs>
  HOSTDEVICE constexpr IIx linearIndex(const Idxs... a_idxs) const noexcept
    {
      static_assert(sizeof...(Idxs) >= D, "Insufficient indices");
      return (m_data - ArrayBase<Impl>::m_impl.m_data) +
        subLinearIndex(a_idxs...);
    }

  /// Subspace linear index
  template <typename... Idxs>
  HOSTDEVICE constexpr IIx subLinearIndex(const Idxs... a_idxs) const noexcept
    {
      static_assert(sizeof...(Idxs) >= D, "Insufficient indices");
      return linearIndexEachElement<D>(
        ArrayBase<Impl>::m_impl.m_stride,
        GetIdx<IIx, Ord, sizeof...(Idxs)>(a_idxs...));
    }

  /// Access the const sub-space data using operator()
  /** \tparam Idxs...     Integer sequence of indices
   *  \param[in]  a_idxs  Indices of element to access.
   *  \return             Const data at the element
   */
  template <typename... Idxs>
  HOSTDEVICE const T& operator()(const Idxs... a_idxs) const noexcept
    {
      static_assert(sizeof...(Idxs) >= D, "Insufficient indices");
      T* const data = m_data +
        linearIndexEachElementFrom1RC<D, Impl_t>(
          ArrayBase<Impl>::m_impl,
          GetIdx<IIx, Ord, sizeof...(Idxs)>(a_idxs...));
      const IIx i0 = GetIdx<IIx, Ord, sizeof...(Idxs)>(a_idxs...)(0);
      ArrayBase<Impl>::m_impl.checkRange0(i0);
      return data[i0];
    }

  /// Access the modifiable sub-space data using operator()
  /** \tparam Idxs...     Integer sequence of indices
   *  \param[in]  a_idxs  Indices of element to access.
   *  \return             Const data at the element
   */
  template <typename... Idxs>
  HOSTDEVICE T& operator()(const Idxs... a_idxs) noexcept
    {
      static_assert(sizeof...(Idxs) >= D, "Insufficient indices");
      T* const data = m_data +
        linearIndexEachElementFrom1RC<D, Impl_t>(
          ArrayBase<Impl>::m_impl,
          GetIdx<IIx, Ord, sizeof...(Idxs)>(a_idxs...));
      const IIx i0 = GetIdx<IIx, Ord, sizeof...(Idxs)>(a_idxs...)(0);
      ArrayBase<Impl>::m_impl.checkRange0(i0);
      return data[i0];
    }

private:
  T* const m_data;
};

// Specialization for D=1
template <typename Impl>
class ArrayTally<Impl, 1> : public ArrayBase<Impl>
{
private:
  using Impl_t = std::decay_t<Impl>;
  using T = typename Impl_t::value_type;
  using IIx = typename Impl_t::index_type;
  using Ord = typename Impl_t::order_type;
public:
  using value_type = typename Impl_t::value_type;
  using index_type = typename Impl_t::index_type;
  using order_type = typename Impl_t::order_type;
  static constexpr unsigned rank = Impl_t::rank;

public:
  HOSTDEVICE ArrayTally(Impl&& a_impl, T* const a_data) noexcept
    :
    ArrayBase<Impl>(std::forward<Impl>(a_impl)),
    m_data(a_data)
    { }

  ArrayTally(const ArrayTally&)            = default;
  ArrayTally(ArrayTally&&)                 = default;
  ArrayTally& operator=(const ArrayTally&) = default;
  ArrayTally& operator=(ArrayTally&&)      = default;
  ~ArrayTally()                            = default;

  // Operator[]
  HOSTDEVICE const T& operator[](const IIx a_idx) const noexcept
    {
      ArrayBase<Impl>::checkRange(0, a_idx);
      return m_data[a_idx];
    }

  // Operator[]
  HOSTDEVICE T& operator[](const IIx a_idx) noexcept
    {
      ArrayBase<Impl>::checkRange(0, a_idx);
      return m_data[a_idx];
    }

  /// Linear Index (revised)
  template <typename... Idxs>
  HOSTDEVICE constexpr IIx linearIndex(const Idxs... a_idxs) const noexcept
    {
      static_assert(sizeof...(Idxs) >= 1, "Insufficient indices");
      return (m_data - ArrayBase<Impl>::m_impl.m_data) +
        subLinearIndex(a_idxs...);
    }

  /// Subspace linear index
  template <typename... Idxs>
  HOSTDEVICE constexpr IIx subLinearIndex(const Idxs... a_idxs) const noexcept
    {
      static_assert(sizeof...(Idxs) >= 1, "Insufficient indices");
      return GetIdx<IIx, Ord, sizeof...(Idxs)>(a_idxs...)(0);
    }

  /// Access the const sub-space data using operator()
  /** \tparam Idxs...     Integer sequence of indices
   *  \param[in]  a_idxs  Indices of element to access.
   *  \return             Const data at the element
   */
  template <typename... Idxs>
  HOSTDEVICE const T& operator()(const Idxs... a_idxs) const noexcept
    {
      static_assert(sizeof...(Idxs) >= 1, "Insufficient indices");
      const IIx i0 = GetIdx<IIx, Ord, sizeof...(Idxs)>(a_idxs...)(0);
      ArrayBase<Impl>::m_impl.checkRange0(i0);
      return m_data[i0];
    }

  /// Access the modifiable sub-space data using operator()
  /** \tparam Idxs...     Integer sequence of indices
   *  \param[in]  a_idxs  Indices of element to access.
   *  \return             Const data at the element
   */
  template <typename... Idxs>
  HOSTDEVICE T& operator()(const Idxs... a_idxs) noexcept
    {
      static_assert(sizeof...(Idxs) >= 1, "Insufficient indices");
      const IIx i0 = GetIdx<IIx, Ord, sizeof...(Idxs)>(a_idxs...)(0);
      ArrayBase<Impl>::m_impl.checkRange0(i0);
      return m_data[i0];
    }

private:
  T* const m_data;
};


/*******************************************************************************
 */
/// Shapes linear memory according to some rank and dimensions and gives access
/**
 *  Template recursion is used to store the size of each dimension and at D==1,
 *  a pointer to the data is also stored.  Access is provided through
 *  multidimensional notation x[i2][i1][i0] (e.g, for rank = 3).  C (row)
 *  ordering is assumed, where i0 (D==1) is unit stride.  The size of the
 *  highest dimension (D==rank) is not required.
 *
 *  3-D (rank==3) example:
 *  We want an array of size i=4 x j=3 x k=2 where i is unit stride.
 *  Construct with a pointer to linear memory and the size of each dimension.
 *                                                                   \code{.cpp}
 *    double *data = new double[2*3*4];
 *    shape::array<double, 3> array3d(data, 2, 3, 4);
 *                                                                      \endcode
 *  Alternatively, use the helper function where template parameters are
 *  deduced.
 *                                                                   \code{.cpp}
 *    auto array3d = shape::make_array(data, 2, 3, 4);
 *                                                                      \endcode
 *  Access the array using operator[] for each dimension
 *                                                                   \code{.cpp}
 *    for (int k = 0; k != 2; ++k)
 *    for (int j = 0; j != 3; ++j)
 *    for (int i = 0; i != 4; ++i)
 *      {
 *        assert(&array3d[k][j][i] == &data[(k*3 + j)*4 + i]);
 *      }
 *                                                                      \endcode
 *
 *  Note that operator[] returns a helper 'Bracket' object of rank D-1.  When
 *  D=2, the data pointer is returned for standard 1-D access.  The cascade of
 *  Bracket objects retains a reference to the original shape::array and
 *  tabulates the indexing math.  The arrangement signficantly reduces the
 *  number of registers used in CUDA on GPUs, where the shape::array is most
 *  likely in shared memory and the 'Bracket' objects are in registers.  An
 *  alternative and simpler option is to have operator[] return a lower-
 *  dimension shape::array object.  While that works well on CPUs, it leads to
 *  excessive register use on GPUs (at least as of version 10.2 in CUDA).
 *
 *  Performance:
 *  Compilers are very effective at optimizing the access.  Constructing the
 *  shape::array is costly, so you would not want to do that for every access.
 *  However access is equivalent in performance to built-in Fortran arrays or
 *  in C, casting to pointers-to-VLAs to make use of built-in C-style
 *  multidimensional indexing.
 *
 *  Non-zero lower bounds:
 *  The indexing assumes a zero-based lower bound.  If otherwise, you must
 *  offset the array data.
 *
 *  \tparam T           Type of data in the array
 *  \tparam R           Rank of the array.
 *  \tparam D           Current dimension during recursion.
 *  \tparam Ord         Whether sequences of dimensions or indices are
 *                      row-ordered (unit stride on right) or column-orderd
 *                      (unit stride on left)
 *  \tparam IIx         Type of integer for indexing the array (default
 *                      std::ptrdiff_t)
 *  
 *//*+*************************************************************************/

template <typename Impl, unsigned D = std::decay_t<Impl>::rank>
class ArrayBracket : public ArrayBase<Impl>
{
private:
  using Impl_t = std::decay_t<Impl>;
  using T = typename Impl_t::value_type;
  using IIx = typename Impl_t::index_type;

public:
  using value_type = typename Impl_t::value_type;
  using index_type = typename Impl_t::index_type;
  using order_type = typename Impl_t::order_type;
  static constexpr unsigned rank = Impl_t::rank;

  HOSTDEVICE ArrayBracket() noexcept
    :
    ArrayBase<Impl>()
    { }

  template <typename... Dims>
  HOSTDEVICE ArrayBracket(T *const a_data, const Dims... a_dims) noexcept
    :
    ArrayBase<Impl>(a_data, a_dims...)
    { }

  ArrayBracket(const ArrayBracket&)            = default;
  ArrayBracket(ArrayBracket&&)                 = default;
  ArrayBracket& operator=(const ArrayBracket&) = default;
  ArrayBracket& operator=(ArrayBracket&&)      = default;
  ~ArrayBracket()                              = default;

  /// Access the constant data using square bracket indexing
  /** \param[in]  a_idx   Index for this dimension
   *  \return             A temporary array of one lower dimension
   */
  HOSTDEVICE const auto operator[](const IIx a_idx) const noexcept
    {
      ArrayBase<Impl>::checkRange(D-1, a_idx);
      return ArrayTally<const Impl&, D-1>(
        ArrayBase<Impl>::m_impl,
        ArrayBase<Impl>::m_impl.m_data +
        a_idx*ArrayBase<Impl>::m_impl.m_stride[D-2]);
    }

  /// Access the modifiable data using square bracket indexing
  /** \param[in]  a_idx   Index for this dimension
   *  \return             A temporary array of one lower dimension
   */
  HOSTDEVICE auto operator[](const IIx a_idx) noexcept
    {
      ArrayBase<Impl>::checkRange(D-1, a_idx);
      return ArrayTally<Impl&, D-1>(
        ArrayBase<Impl>::m_impl,
        ArrayBase<Impl>::m_impl.m_data +
        a_idx*ArrayBase<Impl>::m_impl.m_stride[D-2]);
    }
};


/*==============================================================================
 *
 * Class ArrayBracket, specialization for D=1
 *
 *============================================================================*/

template <typename Impl>
class ArrayBracket<Impl, 1> : public ArrayBase<Impl>
{
private:
  using Impl_t = std::decay_t<Impl>;
  using T = typename Impl_t::value_type;
  using IIx = typename Impl_t::index_type;

public:
  using value_type = typename Impl_t::value_type;
  using index_type = typename Impl_t::index_type;
  using order_type = typename Impl_t::order_type;
  static constexpr unsigned rank = Impl_t::rank;

  HOSTDEVICE ArrayBracket() noexcept
    :
    ArrayBase<Impl>()
    { }

  template <typename... Dims>
  HOSTDEVICE ArrayBracket(T *const a_data, const Dims... a_dims) noexcept
    :
    ArrayBase<Impl>(a_data, a_dims...)
    { }

  ArrayBracket(const ArrayBracket&)            = default;
  ArrayBracket(ArrayBracket&&)                 = default;
  ArrayBracket& operator=(const ArrayBracket&) = default;
  ArrayBracket& operator=(ArrayBracket&&)      = default;
  ~ArrayBracket()                              = default;

  HOSTDEVICE const T& operator[](const IIx a_idx) const noexcept
    {
      ArrayBase<Impl>::m_impl.checkRange0(a_idx);
      return ArrayBase<Impl>::data()[a_idx];
    }

  HOSTDEVICE T& operator[](const IIx a_idx) noexcept
    {
      ArrayBase<Impl>::m_impl.checkRange0(a_idx);
      return ArrayBase<Impl>::data()[a_idx];
    }
};

/*==============================================================================
 * Construction with separate implementation declaration
 *============================================================================*/

template <typename T,
          unsigned R,
          typename Ord = row_ordered,
          typename IIx = int,  // Consider std::ptrdiff_t
          typename Sz = min_size>
#ifndef SHAPE_ARRAY_DEBUG
using bare_impl = BareImpl<T, R, Ord, IIx, Sz>;
#else
using bare_impl = RangeCheckImpl<T, R, Ord, IIx, full_size>;
#endif
template <typename T,
          unsigned R,
          typename Ord = row_ordered,
          typename IIx = int,  // Consider std::ptrdiff_t
          typename Sz = full_size>
#ifndef SHAPE_ARRAY_DEBUG
using bounds_impl = BoundsImpl<T, R, Ord, IIx, full_size>;
#else
using bounds_impl = RangeCheckImpl<T, R, Ord, IIx, full_size>;
#endif
template <typename T,
          unsigned R,
          typename Ord = row_ordered,
          typename IIx = int,  // Consider std::ptrdiff_t
          typename Sz = full_size>
using range_check_impl = RangeCheckImpl<T, R, Ord, IIx, full_size>;

/// Make an array with an implementation
template<typename Impl, typename... Dims>
HOSTDEVICE inline auto
make_array_impl(typename Impl::value_type *a_data, Dims... a_dims) noexcept
{
  static_assert(sizeof...(Dims) >= Impl::rank, "Insufficient dimensions");
  return ArrayBracket<Impl>(a_data, a_dims...);
}

/*==============================================================================
 * Construction with implementation configuration specifier
 *============================================================================*/

/// Declare implementation type using an implementation configuration specifier
template<typename ImplSpec,
         typename T,
         unsigned R,
         typename Ord = row_ordered,
         typename IIx = int,  // Consider std::ptrdiff_t
         typename Sz = full_size>
struct array_config
{ };
template<typename T,
         unsigned R,
         typename Ord,
         typename IIx,
         typename Sz>
struct array_config<bare_config, T, R, Ord, IIx, Sz>
{
#ifndef SHAPE_ARRAY_DEBUG
  using type = BareImpl<T, R, Ord, IIx, Sz>;
#else
  using type = RangeCheckImpl<T, R, Ord, IIx, full_size>;
#endif
};
template<typename T,
         unsigned R,
         typename Ord,
         typename IIx,
         typename Sz>
struct array_config<bounds_config, T, R, Ord, IIx, Sz>
{
#ifndef SHAPE_ARRAY_DEBUG
  using type = BoundsImpl<T, R, Ord, IIx, full_size>;
#else
  using type = RangeCheckImpl<T, R, Ord, IIx, full_size>;
#endif
};
template<typename T,
         unsigned R,
         typename Ord,
         typename IIx,
         typename Sz>
struct array_config<range_check_config, T, R, Ord, IIx, Sz>
{
  using type = RangeCheckImpl<T, R, Ord, IIx, full_size>;
};

/// Make an array with an implementation configuration specifier
template<typename ImplSpec, typename T, typename... Dims>
HOSTDEVICE inline auto
make_array_config(T *a_data, Dims... a_dims) noexcept
{
  static_assert(std::is_same<bare_config, ImplSpec>::value ||
                std::is_same<bounds_config, ImplSpec>::value ||
                std::is_same<range_check_config, ImplSpec>::value,
                "Invalid array implementation configuration specifier");
  return ArrayBracket<
    typename array_config<ImplSpec, T, sizeof...(Dims)>::type>(
      a_data, a_dims...);
}

/// Make row-ordered array with an implementation configuration specifier
template<typename ImplSpec, typename T, typename... Dims>
HOSTDEVICE inline auto
make_row_array_config(T *a_data, Dims... a_dims) noexcept
{
  static_assert(std::is_same<bare_config, ImplSpec>::value ||
                std::is_same<bounds_config, ImplSpec>::value ||
                std::is_same<range_check_config, ImplSpec>::value,
                "Invalid array implementation configuration specifier");
  return ArrayBracket<
    typename array_config<ImplSpec, T, sizeof...(Dims), row_ordered>::type>(
      a_data, a_dims...);
}

/// Make column-ordered array with an implementation configuration specifier
template<typename ImplSpec, typename T, typename... Dims>
HOSTDEVICE inline auto
make_col_array_config(T *a_data, Dims... a_dims) noexcept
{
  static_assert(std::is_same<bare_config, ImplSpec>::value ||
                std::is_same<bounds_config, ImplSpec>::value ||
                std::is_same<range_check_config, ImplSpec>::value,
                "Invalid array implementation configuration specifier");
  return ArrayBracket<
    typename array_config<ImplSpec, T, sizeof...(Dims), column_ordered>::type>(
      a_data, a_dims...);
}

/*==============================================================================
 * Basic construction
 *============================================================================*/

/// Shape array class objects
template <typename T,
          unsigned R,
          typename Ord = row_ordered,
          typename IIx = int,  // Consider std::ptrdiff_t
          typename Sz = min_size>
#ifndef SHAPE_ARRAY_DEBUG
using array = ArrayBracket<BareImpl<T, R, Ord, IIx, Sz>>;
#else
using array = ArrayBracket<RangeCheckImpl<T, R, Ord, IIx, Sz>>;
#endif
template <typename T,
          unsigned R,
          typename Ord = row_ordered,
          typename IIx = int,  // Consider std::ptrdiff_t
          typename Sz = full_size>
#ifndef SHAPE_ARRAY_DEBUG
using array_bounds = ArrayBracket<BoundsImpl<T, R, Ord, IIx, Sz>>;
#else
using array_bounds = ArrayBracket<RangeCheckImpl<T, R, Ord, IIx, Sz>>;
#endif
template <typename T,
          unsigned R,
          typename Ord = row_ordered,
          typename IIx = int,  // Consider std::ptrdiff_t
          typename Sz = full_size>
using array_range_checked = ArrayBracket<RangeCheckImpl<T, R, Ord, IIx, Sz>>;

/// Helper function that builds a row-ordered shape::array using int-type index
/** All template parameters are deduced
 *  \param[in]  a_data  Linear memory to reshaped
 *  \param[in]  a_dims  Size of each dimension (unit stride on right)
 *  \return     shape::array with type given by a_data and rank defined by the
 *              number of dimension lengths.
 */
template <typename T, typename... Dims>
HOSTDEVICE inline auto
make_array(T* a_data, Dims... a_dims) noexcept
{
  return array<T, sizeof...(a_dims)>(
    a_data, a_dims...);
}

/// Helper function that builds a row-ordered shape::array using int-type index
/** All template parameters are deduced
 *  \param[in]  a_data  Linear memory to reshaped
 *  \param[in]  a_dims  Size of each dimension (unit stride on right)
 *  \return     shape::array with type given by a_data and rank defined by the
 *              number of dimension lengths.
 */
template <typename T, typename... Dims>
HOSTDEVICE inline auto
make_row_array(T* a_data, Dims... a_dims) noexcept
{
  return array<T, sizeof...(a_dims), row_ordered>(
    a_data, a_dims...);
}

/// Helper function that builds a column-ordered shape::array using int-type
/** All template parameters are deduced
 *  \param[in]  a_data  Linear memory to reshaped
 *  \param[in]  a_dims  Size of each dimension (unit stride on left)
 *  \return     shape::array with type given by a_data and rank defined by the
 *              number of dimension lengths.
 */
template <typename T, typename... Dims>
HOSTDEVICE inline auto
make_col_array(T* a_data, Dims... a_dims) noexcept
{
  return array<T, sizeof...(a_dims), column_ordered>(
    a_data, a_dims...);
}
/// Helper function that builds a row-ordered shape::array using int-type index
/** All template parameters are deduced
 *  \param[in]  a_data  Linear memory to reshaped
 *  \param[in]  a_dims  Size of each dimension (unit stride on right)
 *  \return     shape::array with type given by a_data and rank defined by the
 *              number of dimension lengths.
 */
template <typename T, typename... Dims>
HOSTDEVICE inline auto
make_array_bounds(T* a_data, Dims... a_dims) noexcept
{
  return array_bounds<T, sizeof...(a_dims)>(
    a_data, a_dims...);
}

/// Helper function that builds a row-ordered shape::array using int-type index
/** All template parameters are deduced
 *  \param[in]  a_data  Linear memory to reshaped
 *  \param[in]  a_dims  Size of each dimension (unit stride on right)
 *  \return     shape::array with type given by a_data and rank defined by the
 *              number of dimension lengths.
 */
template <typename T, typename... Dims>
HOSTDEVICE inline auto
make_row_array_bounds(T* a_data, Dims... a_dims) noexcept
{
  return array_bounds<T, sizeof...(a_dims), row_ordered>(
    a_data, a_dims...);
}

/// Helper function that builds a column-ordered shape::array using int-type
/** All template parameters are deduced
 *  \param[in]  a_data  Linear memory to reshaped
 *  \param[in]  a_dims  Size of each dimension (unit stride on left)
 *  \return     shape::array with type given by a_data and rank defined by the
 *              number of dimension lengths.
 */
template <typename T, typename... Dims>
HOSTDEVICE inline auto
make_col_array_bounds(T* a_data, Dims... a_dims) noexcept
{
  return array_bounds<T, sizeof...(a_dims), column_ordered>(
    a_data, a_dims...);
}
/// Helper function that builds a row-ordered shape::array using int-type index
/** All template parameters are deduced
 *  \param[in]  a_data  Linear memory to reshaped
 *  \param[in]  a_dims  Size of each dimension (unit stride on right)
 *  \return     shape::array with type given by a_data and rank defined by the
 *              number of dimension lengths.
 */
template <typename T, typename... Dims>
HOSTDEVICE inline auto
make_array_range_checked(T* a_data, Dims... a_dims) noexcept
{
  return array_range_checked<T, sizeof...(a_dims)>(
    a_data, a_dims...);
}

/// Helper function that builds a row-ordered shape::array using int-type index
/** All template parameters are deduced
 *  \param[in]  a_data  Linear memory to reshaped
 *  \param[in]  a_dims  Size of each dimension (unit stride on right)
 *  \return     shape::array with type given by a_data and rank defined by the
 *              number of dimension lengths.
 */
template <typename T, typename... Dims>
HOSTDEVICE inline auto
make_row_array_range_checked(T* a_data, Dims... a_dims) noexcept
{
  return array_range_checked<T, sizeof...(a_dims), row_ordered>(
    a_data, a_dims...);
}

/// Helper function that builds a column-ordered shape::array using int-type
/** All template parameters are deduced
 *  \param[in]  a_data  Linear memory to reshaped
 *  \param[in]  a_dims  Size of each dimension (unit stride on left)
 *  \return     shape::array with type given by a_data and rank defined by the
 *              number of dimension lengths.
 */
template <typename T, typename... Dims>
HOSTDEVICE inline auto
make_col_array_range_checked(T* a_data, Dims... a_dims) noexcept
{
  return array_range_checked<T, sizeof...(a_dims), column_ordered>(
    a_data, a_dims...);
}

}  // namespace shape

#include "BaseNamespaceFooter.H"

#endif  /* ! defined _SHAPEARRAY_H_ */
