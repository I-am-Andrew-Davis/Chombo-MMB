
#ifndef _DYNARRAY_H_
#define _DYNARRAY_H_

#include <tuple>

#include "ArrayAllocator.H"
#include "ShapeArray.H"
#include "StcVector.H"

#include "BaseNamespaceHeader.H"

namespace dyn
{

// Heap allocation (allocations and constructions)
using array_format = std::integral_constant<unsigned, 0>;
// Alias constructed data (no allocations or constructions)
using alias_format = std::integral_constant<unsigned, 1>;
// Alias raw memory (only contructions)
using raw_format   = std::integral_constant<unsigned, 2>;


/*******************************************************************************
 */
/// Dynamic array class based on shape::array
/**
 *  Supports array allocation on host and device and multiple types of indexing.
 *
 *  \tparam T           Array element type
 *  \tparam R           Rank (number of dimensions)
 *  \tparam Fmt         Format:
 *                      array_format: allocates on heap and constructs elements
 *                      alias_format: aliases constructed memory
 *                      raw_format:   use existing memory and construct elements
 *  \tparam IIx         Type for indexing
 *  \tparam Ord         Column (default) or row ordering
 *  \tparam Alloc       Allocator for memory:
 *                      DefaultArrayAlloc: (default if CH_GPU=FALSE).  This is
 *                        stateless and only allocates on the CPU.
 *                      CUDAArrayAlloc: (default if CH_GPU=TRUE).  Allocates on
 *                        both system and device memory and permits transfers
 *                        between the two.  Requires 16 bytes for state
 *                        (pointers to the system and device memory).
 *
 *  \note
 *  <ul>
 *    <li> Since the format is included in the type, a function parameter of,
 *         e.g., dyn::Array<Real, 3, dyn::array_format>& cannot be matched with
 *         a type having alias_format or raw_format.  Instead, make the function
 *         parameter have type dyn::Array<Real, 3, dyn::alias_format>.  Note
 *         this is not a reference.  Any array type can implicilty convert to
 *         dyn::Array<Real, 3, dyn::alias_format>.
 *  </ul>
 *//*+*************************************************************************/

template <typename T,
          unsigned R,
          typename IIx = int,
          typename Ord = shape::column_ordered,
          typename Cfg = shape::bare_config,
#ifdef CH_GPU
          template <typename, ArrayClassIndex> class Alloc = CUDAArrayAlloc
#else
          template <typename, ArrayClassIndex> class Alloc = DefaultArrayAlloc
#endif
          >
class Array : public Alloc<T, ArrayClassIndex::DynArray>,
              public shape::ArrayBracket<
                typename shape::array_config<
                  Cfg, T, R, Ord, IIx, shape::full_size>::type>
{


/*==============================================================================
 * Types
 *============================================================================*/

  using alloc_type = Alloc<T, ArrayClassIndex::DynArray>;
/*
  Note: to be able to call delete, the offset for a lower bound is stored in
  order to recover the original pointer.  A few bits therein also store the
  AllocOn and AllocBy.
 */
  using shape_array = shape::ArrayBracket<
    typename shape::array_config<Cfg, T, R, Ord, IIx, shape::full_size>::type>;

  using default_alloc_type = DefaultArrayAlloc<T, ArrayClassIndex::DynArray>;
#ifdef CH_GPU
  using cuda_alloc_type = CUDAArrayAlloc<T, ArrayClassIndex::DynArray>;
#endif

public:

  // Allocator_traits can provide defaults for missing features in the allocator
  // But its functions cannot be used on the GPU
  using AllocTr         = std::allocator_traits<alloc_type>;

  // Types
  using size_type       = typename AllocTr::size_type;
  using difference_type = typename AllocTr::difference_type;
  using pointer         = T*;
  using const_pointer   = const T*;
  using reference       = T&;
  using const_reference = const T&;
  using value_type      = T;

  // Non-types
  static constexpr unsigned rank = R;

  // Array access types
  using index_type      = IIx;

  static constexpr difference_type c_max_size =
    static_cast<difference_type>(~static_cast<size_type>(0) >> 5);
  static constexpr difference_type c_min_size = -c_max_size - 1;

/*--------------------------------------------------------------------*/
/// Default constructor
/**
 *//*-----------------------------------------------------------------*/

  HOSTDEVICE
  Array() noexcept
    :
    shape_array(),
    m_param(alloc_type::allocates_on)
    {
      staticChecks();
    }

#ifdef CH_GPU
/*--------------------------------------------------------------------*/
/// Perform default initialization
/** This should only be necessary for dynamic arrays created as
 *  shared variables on the GPU, e.g., a declaration prefixed with
 *  __shared__.  In that case, there is no initialization.  This must
 *  be used before one of the define methods.
 *//*-----------------------------------------------------------------*/

  DEVICE void
  INITIALIZEshared() noexcept
    {
      staticChecks();
      shape_array::clear();
      // Next also clears offset and allocBy
      m_param = static_cast<size_type>(AllocOn::gpu);
    }
#endif

/*--------------------------------------------------------------------*/
/// Copy constructor that swaps pointers to data on the GPU
/** This is a shallow copy that only copies and configures the class
 *  object.  The data must be moved separately.
 *//*-----------------------------------------------------------------*/


/*--------------------------------------------------------------------*/
/// Constructor with custom AllocOn (as AllocOn)
/**
 *//*-----------------------------------------------------------------*/

  HOSTDEVICE
  Array(const AllocOn a_allocOn) noexcept
    :
    m_param(static_cast<size_type>(a_allocOn)),
    shape_array()
    {
      staticChecks();
    }

/*--------------------------------------------------------------------*/
/// Constructor with custom AllocOn (as unsigned)
/**
 *//*-----------------------------------------------------------------*/

  HOSTDEVICE
  Array(const unsigned a_allocOn) noexcept
    :
    m_param(a_allocOn),
    shape_array()
    {
      staticChecks();
    }

/*--------------------------------------------------------------------*/
/// Constructor (array format)
/** \tparam Dim0        Type of first dimension (must be integral)
 *  \tparam DimArgs...  Remaining dimension types and optional
 *                      constructor types
 *  \param[in]  a_dim0  First dimension
 *  \param[in]  a_dimArgs
 *                      Remaining dimensions and optional arguments to
 *                      the constructor
 *  \note
 *  <ul>
 *    <li> Dim0 is separated to avoid conflict of a zero-sized
 *         parameter pack with the default constructor.  Also, testing
 *         the type of tuple_element 0 fails with no arguments (if
 *         just testing if this function is valid)
 *    <li> Constructor arguments are not supported on the GPU.
 *  </ul>
 *//*-----------------------------------------------------------------*/
  
  template <typename Dim0, typename... DimArgs,
            std::enable_if_t<std::is_integral<Dim0>::value, int> = 0>
  HOSTDEVICE
  Array(const Dim0& a_dim0, const DimArgs&... a_dimArgs) noexcept
    :
    shape_array(
      alloc_type::allocate(
        stc::product(stc::Vector<index_type, R>{ a_dim0, a_dimArgs... }),
        // Pass allocOn as "hint" for allocation.
        reinterpret_cast<const_pointer>(alloc_type::allocates_on)),
      stc::Vector<index_type, R>{ a_dim0, a_dimArgs... }),
    m_param(alloc_type::allocates_on | static_cast<unsigned>(AllocBy::array))
    {
      staticChecks();
      static_assert(1 + sizeof...(DimArgs) >= R, "Insufficient dimensions");
      constexpr int NArgs = 1 + sizeof...(DimArgs) - R;
      if (shape_array::data() != nullptr)
        {
          constructArray<R>(
            std::make_index_sequence<NArgs>{},
            std::integral_constant<bool, (NArgs > 0)>{},
            std::integral_constant<bool,
              (NArgs > 0) || std::is_same<T, Real>::value >{},
            a_dim0, a_dimArgs...);
        }
    }

/*--------------------------------------------------------------------*/
/// Constructor (array format) from lower and upper vectors of indices
/** The array indexing is a_lb <= i <= a_ub and the dimensions are
 *  dim = a_ub - a_lb + 1
 *  \tparam VecT        Type in the vector
 *  \tparam VecR        Rank of the vector
 *  \tparam Args...     Optional constructor types
 *  \param[in]  a_lb    Lower bound of the array
 *  \param[in]  a_ub    Upper bound of the array
 *  \param[in]  a_args  Optional arguments to the constructor
 *  \note
 *  <ul>
 *    <li> Constructor arguments are not supported on the GPU.
 *  </ul>
 *//*-----------------------------------------------------------------*/

  template <typename VecT, stc::array_size_type VecR, typename... Args,
            std::enable_if_t<std::is_integral<VecT>::value, int> = 0>
  HOSTDEVICE
  Array(const stc::Vector<VecT, VecR>& a_lb,
        const stc::Vector<VecT, VecR>& a_ub,
        const Args&...                 a_args) noexcept
    :
    shape_array(
      alloc_type::allocate(
        stc::product(a_ub - a_lb + 1),
        // Pass allocOn as "hint" for allocation.
        reinterpret_cast<const_pointer>(alloc_type::allocates_on)),
      stc::Vector<index_type, R>(a_ub - a_lb + 1)),
    m_param(alloc_type::allocates_on | static_cast<unsigned>(AllocBy::array))
    {
      staticChecks();
      static_assert(VecR >= R, "Insufficient dimensions");
      CH_assert((stc::Vector<index_type, R>(a_ub - a_lb + 1) > (index_type)0));
      constexpr int NArgs = sizeof...(Args);
      if (shape_array::data() != nullptr)
        {
          constructArray<0>(
            std::make_index_sequence<NArgs>{},
            std::integral_constant<bool, (NArgs > 0)>{},
            std::integral_constant<bool,
              (NArgs > 0) || std::is_same<T, Real>::value >{},
            a_args...);
        }
      setOffset(shape_array::linearIndex(a_lb));
      shape_array::resetData(shape_array::data() - offset());
      shape_array::setRelativeLBIndices(a_lb);
    }

/*--------------------------------------------------------------------*/
/// Constructor (alias format)
/** \tparam Dim0        Type for first dimension
 *  \tparam Dims...     Type for remaining dimensions
 *  \param[in]  a_addr  Address to alias.  T must already be
 *                      constructed on the address.
 *  \param[in]  a_dim0  Size of first dimension
 *  \param[in]  a_dims  Sizes of remaining dimensions
 *  \note
 *  <ul>
 *    <li> Pointers in the allocator are not defined.  If you need to
 *         set an alias to GPU memory, you have to configure that
 *         manually using setDeviceAddress.
 *  </ul>
 *//*-----------------------------------------------------------------*/

  template <typename Dim0, typename... Dims,
            std::enable_if_t<std::is_integral<Dim0>::value, int> = 0>
  HOSTDEVICE 
  Array(alias_format,
        pointer const  a_addr,
        const Dim0&    a_dim0, 
        const Dims&... a_dims) noexcept
    :
    shape_array(a_addr, a_dim0, a_dims...),
    m_param(alloc_type::allocates_on | static_cast<unsigned>(AllocBy::alias))
    {
      staticChecks();
      static_assert(1 + sizeof...(Dims) >= R, "Insufficient dimensions");
    }

/*--------------------------------------------------------------------*/
/// Constructor (alias format) from lower and upper vectors of indices
/** The array indexing is a_lb <= i <= a_ub and the dimensions are
 *  dim = a_ub - a_lb + 1
 *  \tparam VecT        Type in the vector
 *  \tparam VecR        Rank of the vector
 *  \param[in]  a_addr  Address to alias.  T must already be
 *                      constructed on the address.
 *  \param[in]  a_lb    Lower bound of the array
 *  \param[in]  a_ub    Upper bound of the array
 *  \note
 *  <ul>
 *    <li> Pointers in the allocator are not defined.  If you need to
 *         set an alias to GPU memory, you have to configure that
 *         manually using setDeviceAddress.
 *  </ul>
 *//*-----------------------------------------------------------------*/

  template <typename VecT, stc::array_size_type VecR,
            std::enable_if_t<std::is_integral<VecT>::value, int> = 0>
  HOSTDEVICE 
  Array(alias_format,
        pointer const                  a_addr,
        const stc::Vector<VecT, VecR>& a_lb,
        const stc::Vector<VecT, VecR>& a_ub) noexcept
    :
    shape_array(a_addr, stc::Vector<index_type, R>(a_ub - a_lb + 1)),
    m_param(alloc_type::allocates_on | static_cast<unsigned>(AllocBy::alias))
    {
      staticChecks();
      static_assert(VecR >= R, "Insufficient dimensions");
      CH_assert((stc::Vector<index_type, R>(a_ub - a_lb + 1) > (index_type)0));
      setOffset(shape_array::linearIndex(a_lb));
      shape_array::resetData(shape_array::data() - offset());
      shape_array::setRelativeLBIndices(a_lb);
    }

/*--------------------------------------------------------------------*/
/// Quasi copy constructor (to alias format)
/** The first argument forces explicit use from the perspective of the
 *  user, and acknowledges that this is not a deep copy.
 *  \param[in]  a_array Array to alias.
 *//*-----------------------------------------------------------------*/

  HOSTDEVICE
  Array(alias_format, const Array& a_array) noexcept
    :
    alloc_type(alloc_type::select_on_container_copy_construction(a_array)),
    shape_array(a_array),
    m_param(a_array.m_param)
    {
      setAllocBy(AllocBy::alias);
    }

/*--------------------------------------------------------------------*/
/// Constructor (raw format)
/** \tparam Dim0        Type for first dimension
 *  \tparam DimArgs...  Remaining dimension types and optional
 *                      constructor types
 *  \param[in]  a_addr  Address to alias.  T will be constructed on
 *                      the address.
 *  \param[in]  a_dim0  Size of first dimension
 *  \param[in]  a_dims  Sizes of remaining dimensions
 *  \note
 *  <ul>
 *  \note
 *  <ul>
 *    <li> Pointers in the allocator are not defined.  If you need to
 *         set an alias to GPU memory, you have to configure that
 *         manually using setDeviceAddress.
 *    <li> Constructor arguments are not supported on the GPU.
 *  </ul>
 *//*-----------------------------------------------------------------*/

  template <typename Dim0, typename... DimArgs,
            std::enable_if_t<std::is_integral<Dim0>::value, int> = 0>
  HOSTDEVICE 
  Array(raw_format,
        void *const       a_addr,
        const Dim0&       a_dim0, 
        const DimArgs&... a_dimArgs) noexcept
    :
    shape_array(static_cast<pointer>(a_addr),
                stc::Vector<index_type, R>{ a_dim0, a_dimArgs... }),
    m_param(alloc_type::allocates_on | static_cast<unsigned>(AllocBy::raw))
    {
      staticChecks();
      static_assert(1 + sizeof...(DimArgs) >= R, "Insufficient dimensions");
      constexpr int NArgs = 1 + sizeof...(DimArgs) - R;
      if (shape_array::data() != nullptr)
        {
          constructArray<R>(
            std::make_index_sequence<NArgs>{},
            std::integral_constant<bool, (NArgs > 0)>{},
            std::integral_constant<bool,
              (NArgs > 0) || std::is_same<T, Real>::value >{},
            a_dim0, a_dimArgs...);
        }
    }

/*--------------------------------------------------------------------*/
/// Constructor (raw format) from lower and upper vectors of indices
/** The array indexing is a_lb <= i <= a_ub and the dimensions are
 *  dim = a_ub - a_lb + 1
 *  \tparam VecT        Type in the vector
 *  \tparam VecR        Rank of the vector
 *  \tparam Args...     Optional constructor types
 *  \param[in]  a_addr  Address to alias.  T will be constructed on
 *                      the address.
 *  \param[in]  a_lb    Lower bound of the array
 *  \param[in]  a_ub    Upper bound of the array
 *  \param[in]  a_args  Optional arguments to the constructor
 *  \note
 *  <ul>
 *    <li> Pointers in the allocator are not defined.  If you need to
 *         set an alias to GPU memory, you have to configure that
 *         manually using setDeviceAddress.
 *    <li> Constructor arguments are not supported on the GPU.
 *  </ul>
 *//*-----------------------------------------------------------------*/

  template <typename VecT, stc::array_size_type VecR, typename... Args,
            std::enable_if_t<std::is_integral<VecT>::value, int> = 0>
  HOSTDEVICE 
  Array(raw_format,
        void *const                    a_addr,
        const stc::Vector<VecT, VecR>& a_lb,
        const stc::Vector<VecT, VecR>& a_ub,
        const Args&...                 a_args) noexcept
    :
    shape_array(static_cast<pointer>(a_addr),
                stc::Vector<index_type, R>(a_ub - a_lb + 1)),
    m_param(alloc_type::allocates_on | static_cast<unsigned>(AllocBy::raw))
    {
      staticChecks();
      static_assert(VecR >= R, "Insufficient dimensions");
      CH_assert((stc::Vector<index_type, R>(a_ub - a_lb + 1) > (index_type)0));
      constexpr int NArgs = sizeof...(Args);
      if (shape_array::data() != nullptr)
        {
          constructArray<0>(
            std::make_index_sequence<NArgs>{},
            std::integral_constant<bool, (NArgs > 0)>{},
            std::integral_constant<bool,
              (NArgs > 0) || std::is_same<T, Real>::value >{},
            a_args...);
        }
      setOffset(shape_array::linearIndex(a_lb));
      shape_array::resetData(shape_array::data() - offset());
      shape_array::setRelativeLBIndices(a_lb);
    }

/*--------------------------------------------------------------------*/
/// Move construction
/**
 *//*-----------------------------------------------------------------*/

  HOSTDEVICE
  Array(Array&& a_src) noexcept
    :
    alloc_type(std::move(a_src)),   // Allocator is moved
    shape_array(std::move(a_src)),  // shape_array is moved
    m_param(a_src.m_param)
    {
      a_src.clear();  // Make a_src degenerate
      a_src.m_param = 0;
    }

/*--------------------------------------------------------------------*/
/// Move assignment
/**
 *//*-----------------------------------------------------------------*/

  HOSTDEVICE Array&
  operator=(Array&& a_src) noexcept
    {
      constexpr bool moveStorage =
        AllocTr::propagate_on_container_move_assignment::value ||
        alloc_type::is_always_equal::value;  // Note: cannot use AllocTr until
                                             // c++17
      // Either the allocators are equal or the allocator is propagated in which
      // case the move is constant time.  If we actually have to copy each
      // element, well... we don't allow that.
      static_assert(moveStorage, "Dynamic arrays do not support copying "
                    "elements in a move operation to satisfy allocator "
                    "constraints");
      undefine();
      // Move the allocator (or not) using tag dispatch
      move_allocator(typename AllocTr::propagate_on_container_move_assignment{},
                     a_src.get_allocator());
      // Move the array
      *static_cast<shape_array*>(this) = std::move(a_src);
      m_param = a_src.m_param;
      // Make a_src degenerate
      a_src.clear();
      a_src.m_param = 0;
      return *this;
    }

private:

/*--------------------------------------------------------------------*/
/// Move allocator
/**
 *//*-----------------------------------------------------------------*/

  HOSTDEVICE void
  move_allocator(std::true_type, alloc_type& a_alloc)
    {
      get_allocator() = std::move(a_alloc);
    }

/*--------------------------------------------------------------------*/
/// Unnecessary to move allocator (keep using current)
/**
 *//*-----------------------------------------------------------------*/

  HOSTDEVICE void
  move_allocator(std::false_type, alloc_type& a_alloc)
    { }

public:

/*--------------------------------------------------------------------*/
/// Destructor
/**
 *//*-----------------------------------------------------------------*/

  HOSTDEVICE 
  ~Array()
    {
      undefine();
    }

#ifdef CH_GPU
/*--------------------------------------------------------------------*/
/// Copy constructor that swaps pointers to data on the GPU
/** This is a shallow copy that only copies and configures the class
 *  object.  The data must be moved separately.
 *//*-----------------------------------------------------------------*/

  Array(CH_Cuda::copy_to_device, const Array& a_array) noexcept
    :
    alloc_type(alloc_type::select_on_container_copy_construction(a_array)),
    shape_array(a_array),
    m_param(a_array.m_param)
    {
      // This becomes an alias because the source array on the host handles the
      // memory
      setAllocBy(AllocBy::alias);
      // This will fail at runtime if memory was not allocated on the GPU.
      CH_assert(a_array.allocOn_gpu());
      // Reset data to point to device
      shape_array::resetData(alloc_type::devicePtr() - offset());
    }
#endif


/*==============================================================================
 * Members functions
 *============================================================================*/

//--Defines

/*--------------------------------------------------------------------*/
/// Define (array format)
/** \tparam Dim0        Type of first dimension (must be integral)
 *  \tparam DimArgs...  Remaining dimension types and optional
 *                      constructor types
 *  \param[in]  a_dim0  First dimension
 *  \param[in]  a_dimArgs
 *                      Remaining dimensions and optional arguments
 *//*-----------------------------------------------------------------*/

  template <typename Dim0, typename... DimArgs,
            std::enable_if_t<std::is_integral<Dim0>::value, int> = 0>
  HOSTDEVICE void
  define(const Dim0& a_dim0, const DimArgs&... a_dimArgs) noexcept
    {
      static_assert(1 + sizeof...(DimArgs) >= R, "Insufficient dimensions");
      undefine();
      setAllocBy(AllocBy::array);
      const stc::Vector<index_type, R> dims{ a_dim0, a_dimArgs... };
      shape_array::define(
        alloc_type::allocate(stc::product(dims),
                             reinterpret_cast<const_pointer>(m_param)),
        dims);
      constexpr int NArgs = 1 + sizeof...(DimArgs) - R;
      if (shape_array::data() != nullptr)
        {
          constructArray<R>(
            std::make_index_sequence<NArgs>{},
            std::integral_constant<bool, (NArgs > 0)>{},
            std::integral_constant<bool,
              (NArgs > 0) || std::is_same<T, Real>::value >{},
            a_dim0, a_dimArgs...);
        }
    }

/*--------------------------------------------------------------------*/
/// Define (array format) from lower and upper vectors of indices
/** The array indexing is a_lb <= i <= a_ub and the dimensions are
 *  dim = a_ub - a_lb + 1
 *  \tparam VecT        Type in the vector
 *  \tparam VecR        Rank of the vector
 *  \tparam Args...     Optional constructor types
 *  \param[in]  a_lb    Lower bound of the array
 *  \param[in]  a_ub    Upper bound of the array
 *  \param[in]  a_args  Optional arguments to the constructor
 *//*-----------------------------------------------------------------*/

  template <typename VecT, stc::array_size_type VecR, typename... Args,
            std::enable_if_t<std::is_integral<VecT>::value, int> = 0>
  HOSTDEVICE void
  define(const stc::Vector<VecT, VecR>& a_lb,
         const stc::Vector<VecT, VecR>& a_ub,
         const Args&...                 a_args) noexcept
    {
      static_assert(VecR >= R, "Insufficient dimensions");
      undefine();
      setAllocBy(AllocBy::array);
      const stc::Vector<index_type, R> dims(a_ub - a_lb + 1);
      CH_assert(dims > (index_type)0);
      shape_array::define(
        alloc_type::allocate(stc::product(dims),
                             reinterpret_cast<const_pointer>(m_param)),
        dims);
      constexpr int NArgs = sizeof...(Args);
      if (shape_array::data() != nullptr)
        {
          constructArray<0>(
            std::make_index_sequence<NArgs>{},
            std::integral_constant<bool, (NArgs > 0)>{},
            std::integral_constant<bool,
              (NArgs > 0) || std::is_same<T, Real>::value >{},
            a_args...);
        }
      setOffset(shape_array::linearIndex(a_lb));
      shape_array::resetData(shape_array::data() - offset());
      shape_array::setRelativeLBIndices(a_lb);
    }

/*--------------------------------------------------------------------*/
/// Define (alias format).  Alias constructed T.
/** \tparam Dim0        Type for first dimension
 *  \tparam Dims...     Type for remaining dimensions
 *  \param[in]  a_addr  Address to alias.  T must already be
 *                      constructed on the address.
 *  \param[in]  a_dim0  Size of first dimension
 *  \param[in]  a_dims  Sizes of remaining dimensions
 *  \note
 *  <ul>
 *    <li> Pointers in the allocator are not defined.  If you need to
 *         set an alias to GPU memory, you have to configure that
 *         manually using setDeviceAddress.
 *  </ul>
 *//*-----------------------------------------------------------------*/

  template <typename Dim0, typename... Dims,
            std::enable_if_t<std::is_integral<Dim0>::value, int> = 0>
  HOSTDEVICE void
  defineAlias(pointer const  a_addr,
              const Dim0&    a_dim0,
              const Dims&... a_dims) noexcept
    {
      static_assert(1 + sizeof...(Dims) >= R, "Insufficient dimensions");
      undefine();
      setAllocBy(AllocBy::alias);
      shape_array::define(a_addr, a_dim0, a_dims...);
    }

/*--------------------------------------------------------------------*/
/// Define (alias format).  Alias constructed T.
/** The array indexing is a_lb <= i <= a_ub and the dimensions are
 *  dim = a_ub - a_lb + 1
 *  \tparam VecT        Type in the vector
 *  \tparam VecR        Rank of the vector
 *  \param[in]  a_addr  Address to alias.  T must already be
 *                      constructed on the address.
 *  \param[in]  a_lb    Lower bound of the array
 *  \param[in]  a_ub    Upper bound of the array
 *  \note
 *  <ul>
 *    <li> Pointers in the allocator are not defined.  If you need to
 *         set an alias to GPU memory, you have to configure that
 *         manually using setDeviceAddress.
 *  </ul>
 *//*-----------------------------------------------------------------*/

  template <typename VecT, stc::array_size_type VecR,
            std::enable_if_t<std::is_integral<VecT>::value, int> = 0>
  HOSTDEVICE void
  defineAlias(pointer const                  a_addr,
              const stc::Vector<VecT, VecR>& a_lb,
              const stc::Vector<VecT, VecR>& a_ub) noexcept
    {
      static_assert(VecR >= R, "Insufficient dimensions");
      undefine();
      setAllocBy(AllocBy::alias);
      const stc::Vector<index_type, R> dims(a_ub - a_lb + 1);
      CH_assert(dims > (index_type)0);
      shape_array::define(a_addr, dims);
      setOffset(shape_array::linearIndex(a_lb));
      shape_array::resetData(shape_array::data() - offset());
      shape_array::setRelativeLBIndices(a_lb);
    }

/*--------------------------------------------------------------------*/
/// Define (alias format).  Alias constructed T.
/** \param[in]  a_array Array to alias
 *  \note
 *  <ul>
 *    <li> Pointers in the allocator are not defined.  If you need to
 *         set an alias to GPU memory, you have to configure that
 *         manually using setDeviceAddress.
 *  </ul>
 *//*-----------------------------------------------------------------*/

  HOSTDEVICE void
  defineAlias(const Array& a_array) noexcept
    {
      undefine();
      setAllocBy(AllocBy::alias);
      shape_array::operator=(a_array);
      setOffset(a_array.offset());
    }

/*--------------------------------------------------------------------*/
/// Define (raw format).  Construct T on raw memory
/** \tparam Dim0        Type for first dimension
 *  \tparam DimsArgs... Remaining dimension types and optional
 *                      constructor types
 *  \param[in]  a_addr  Address to alias.  T must already be
 *                      constructed on the address.
 *  \param[in]  a_dim0  Size of first dimension
 *  \param[in]  a_dims  Sizes of remaining dimensions
 *  \note
 *  <ul>
 *    <li> Pointers in the allocator are not defined.  If you need to
 *         set an alias to GPU memory, you have to configure that
 *         manually using setDeviceAddress.
 *  </ul>
 *//*-----------------------------------------------------------------*/

  template <typename Dim0, typename... DimArgs,
            std::enable_if_t<std::is_integral<Dim0>::value, int> = 0>
  HOSTDEVICE void
  defineRaw(void *const       a_addr,
            const Dim0&       a_dim0,
            const DimArgs&... a_dimArgs) noexcept
    {
      static_assert(1 + sizeof...(DimArgs) >= R, "Insufficient dimensions");
      undefine();
      setAllocBy(AllocBy::raw);
      shape_array::define(static_cast<pointer>(a_addr),
                          stc::Vector<index_type, R>{ a_dim0, a_dimArgs... });
      constexpr int NArgs = 1 + sizeof...(DimArgs) - R;
      if (shape_array::data() != nullptr)
        {
          constructArray<R>(
            std::make_index_sequence<NArgs>{},
            std::integral_constant<bool, (NArgs > 0)>{},
            std::integral_constant<bool,
              (NArgs > 0) || std::is_same<T, Real>::value >{},
            a_dim0, a_dimArgs...);
        }
    }

/*--------------------------------------------------------------------*/
/// Define (raw format).  Alias constructed T.
/** The array indexing is a_lb <= i <= a_ub and the dimensions are
 *  dim = a_ub - a_lb + 1
 *  \tparam VecT        Type in the vector
 *  \tparam VecR        Rank of the vector
 *  \tparam Args...     Optional constructor types
 *  \param[in]  a_addr  Address to alias.  T must already be
 *                      constructed on the address.
 *  \param[in]  a_lb    Lower bound of the array
 *  \param[in]  a_ub    Upper bound of the array
 *  \param[in]  a_args  Optional arguments to the constructor
 *  \note
 *  <ul>
 *    <li> Pointers in the allocator are not defined.  If you need to
 *         set an alias to GPU memory, you have to configure that
 *         manually using setDeviceAddress.
 *  </ul>
 *//*-----------------------------------------------------------------*/

  template <typename VecT, stc::array_size_type VecR, typename... Args,
            std::enable_if_t<std::is_integral<VecT>::value, int> = 0>
  HOSTDEVICE void
  defineRaw(void *const                    a_addr,
            const stc::Vector<VecT, VecR>& a_lb,
            const stc::Vector<VecT, VecR>& a_ub,
            const Args&...                 a_args) noexcept
    {
      static_assert(VecR >= R, "Insufficient dimensions");
      undefine();
      setAllocBy(AllocBy::raw);
      const stc::Vector<index_type, R> dims(a_ub - a_lb + 1);
      CH_assert(dims > (index_type)0);
      shape_array::define(static_cast<pointer>(a_addr), dims);
      constexpr int NArgs = sizeof...(Args);
      if (shape_array::data() != nullptr)
        {
          constructArray<0>(
            std::make_index_sequence<NArgs>{},
            std::integral_constant<bool, (NArgs > 0)>{},
            std::integral_constant<bool,
              (NArgs > 0) || std::is_same<T, Real>::value >{},
            a_args...);
        }
      setOffset(shape_array::linearIndex(a_lb));
      shape_array::resetData(shape_array::data() - offset());
      shape_array::setRelativeLBIndices(a_lb);
    }

/*--------------------------------------------------------------------*/
/// Deallocate and destroy T depending on AllocBy
/**
 *//*-----------------------------------------------------------------*/

  HOSTDEVICE void
  undefine() noexcept
    {
      if (shape_array::data() != nullptr)
        {
          switch (allocBy())
            {
            case AllocBy::array:
              destroyArray(std::is_class<T>());
              alloc_type::deallocate(data(), shape_array::size());
              break;
            case AllocBy::raw:
              destroyArray(std::is_class<T>());
              break;
            case AllocBy::alias:
            case AllocBy::none:
              break;
            }
        }
      shape_array::clear();
      clearOffsetParam();
    }

/*--------------------------------------------------------------------*/
/// Shift the array indices
/** \tparam Idxs        Types for indices
 *  \param[in]  a_idxs  Values to shift by
 *//*-----------------------------------------------------------------*/

  template <typename... Idxs>
  HOSTDEVICE void
  shift(const Idxs&... a_idxs) noexcept
    {
      static_assert(
        std::is_integral<
          typename std::tuple_element<0, std::tuple<Idxs...>>::type
        >::value,
        "Shift argument must be integral");
      static_assert(sizeof...(Idxs) >= R, "Insufficient indices");
      const index_type relOffset = shape_array::linearIndex(a_idxs...);
      setOffset(offset() + relOffset);
      shape_array::resetData(shape_array::data() - relOffset);
      shape_array::setRelativeLBIndices(a_idxs...);
    }
 
/*--------------------------------------------------------------------*/
/// Shift the array indices
/** \tparam VecT        Type in the vector
 *  \tparam VecR        Rank of the vector
 *  \param[in]  a_shiftBy
 *                      Vector of indices to shift by
 *//*-----------------------------------------------------------------*/

  template <typename VecT, stc::array_size_type VecR>
  HOSTDEVICE void
  shift(const stc::Vector<VecT, VecR>& a_shiftBy) noexcept
    {
      static_assert(std::is_integral<VecT>::value,
                    "Shift argument must be integral");
      static_assert(VecR >= R, "Insufficient indices");
      const index_type relOffset = shape_array::linearIndex(a_shiftBy);
      setOffset(offset() + relOffset);
      shape_array::resetData(shape_array::data() - relOffset);
      shape_array::setRelativeLBIndices(a_shiftBy);
    }

/*--------------------------------------------------------------------*/
/// Initialize the array to a constant value
/** \param[in]  a_val   Value to initialize to
 *//*-----------------------------------------------------------------*/

  HOSTDEVICE void
  operator=(const T& a_val) noexcept
    {
      pointer p = data();
      for (size_type n = shape_array::size(); n--;)
        {
          *p++ = a_val;
        }
    }


/*==============================================================================
 * Accessors (all operator() are redefined because of new ones)
 *============================================================================*/

/*--------------------------------------------------------------------*/
/// Const access using indices (from shape_array)
/** 
 *//*-----------------------------------------------------------------*/

  template <typename... Idxs,
            std::enable_if_t<std::is_integral<
                               typename std::tuple_element<
                                 0, std::tuple<Idxs...>>::type>::value,
                             int> = 0>
  HOSTDEVICE const_reference
  operator()(const Idxs... a_idxs) const noexcept
    {
      return shape_array::operator()(a_idxs...);
    }

/*--------------------------------------------------------------------*/
/// Const access using a vector (from shape_array)
/** 
 *//*-----------------------------------------------------------------*/

  template <typename Vec,
            std::enable_if_t<std::is_pointer<Vec>::value ||
                             std::is_array<Vec>::value ||
                             std::is_class<Vec>::value, int> = 0>
  HOSTDEVICE const_reference
  operator()(const Vec& a_vec) const noexcept
    {
      return shape_array::operator()(a_vec);
    }

/*--------------------------------------------------------------------*/
/// Const access using a vector plus other indices
/** 
 *//*-----------------------------------------------------------------*/

  template <typename VecIIx, unsigned VecR, typename... Idxs,
            std::enable_if_t<std::is_integral<VecIIx>::value &&
                             (VecR < R) &&
                             (VecR + 1 + sizeof...(Idxs) >= R), int> = 0>
  HOSTDEVICE const_reference
  operator()(const stc::Vector<VecIIx, VecR>& a_vec,
             const IIx                        a_otherIdx0,
             const Idxs...                    a_otherIdxs) const noexcept
    {
      return shape_array::operator()(
        stc::Vector<IIx, R>(a_vec, a_otherIdx0, a_otherIdxs...));
    }

/*--------------------------------------------------------------------*/
/// Const access using a single leading index plus a vector
/** 
 *//*-----------------------------------------------------------------*/

  template <typename VecIIx, unsigned VecR,
            std::enable_if_t<std::is_integral<VecIIx>::value &&
                             (VecR < R) &&
                             (VecR + 1 >= R), int> = 0>
  HOSTDEVICE const_reference
  operator()(const IIx                        a_idx0,
             const stc::Vector<VecIIx, VecR>& a_vec) const noexcept
    {
      return shape_array::operator()(stc::Vector<IIx, R>(a_idx0, a_vec));
    }

/*--------------------------------------------------------------------*/
/// Const access using two leading indices plus a vector
/** 
 *//*-----------------------------------------------------------------*/

  template <typename VecIIx, unsigned VecR,
            std::enable_if_t<std::is_integral<VecIIx>::value &&
                             (VecR < R) &&
                             (VecR + 2 >= R), int> = 0>
  HOSTDEVICE const_reference
  operator()(const IIx                        a_idx0,
             const IIx                        a_idx1,
             const stc::Vector<VecIIx, VecR>& a_vec) const noexcept
    {
      return shape_array::operator()(
        stc::Vector<IIx, R>(a_idx0, a_idx1, a_vec));
    }

/*--------------------------------------------------------------------*/
/// Modifiable access using indices (from shape_array)
/** 
 *//*-----------------------------------------------------------------*/

  template <typename... Idxs,
            std::enable_if_t<std::is_integral<
                               typename std::tuple_element<
                                 0, std::tuple<Idxs...>>::type>::value,
                             int> = 0>
  HOSTDEVICE reference
  operator()(const Idxs... a_idxs) noexcept
    {
      return shape_array::operator()(a_idxs...);
    }

/*--------------------------------------------------------------------*/
/// Modifiable access using a vector (from shape_array)
/** 
 *//*-----------------------------------------------------------------*/

  template <typename Vec,
            std::enable_if_t<std::is_pointer<Vec>::value ||
                             std::is_array<Vec>::value ||
                             std::is_class<Vec>::value, int> = 0>
  HOSTDEVICE reference
  operator()(const Vec& a_vec) noexcept
    {
      return shape_array::operator()(a_vec);
    }

/*--------------------------------------------------------------------*/
/// Modifiable access using a vector plus other indices
/** 
 *//*-----------------------------------------------------------------*/

  template <typename VecIIx, unsigned VecR, typename... Idxs,
            std::enable_if_t<std::is_integral<VecIIx>::value &&
                             (VecR < R) &&
                             (VecR + 1 + sizeof...(Idxs) >= R), int> = 0>
  HOSTDEVICE reference
  operator()(const stc::Vector<VecIIx, VecR>& a_vec,
             const IIx                        a_otherIdx0,
             const Idxs...                    a_otherIdxs) noexcept
    {
      return shape_array::operator()(
        stc::Vector<IIx, R>(a_vec, a_otherIdx0, a_otherIdxs...));
    }

/*--------------------------------------------------------------------*/
/// Modifiable access using a single leading index plus a vector
/** 
 *//*-----------------------------------------------------------------*/

  template <typename VecIIx, unsigned VecR,
            std::enable_if_t<std::is_integral<VecIIx>::value &&
                             (VecR < R) &&
                             (VecR + 1 >= R), int> = 0>
  HOSTDEVICE reference
  operator()(const IIx                        a_idx0,
             const stc::Vector<VecIIx, VecR>& a_vec) noexcept
    {
      return shape_array::operator()(stc::Vector<IIx, R>(a_idx0, a_vec));
    }

/*--------------------------------------------------------------------*/
/// Modifiable access using two leading indices plus a vector
/** 
 *//*-----------------------------------------------------------------*/

  template <typename VecIIx, unsigned VecR,
            std::enable_if_t<std::is_integral<VecIIx>::value &&
                             (VecR < R) &&
                             (VecR + 2 >= R), int> = 0>
  HOSTDEVICE reference
  operator()(const IIx                        a_idx0,
             const IIx                        a_idx1,
             const stc::Vector<VecIIx, VecR>& a_vec) noexcept
    {
      return shape_array::operator()(
        stc::Vector<IIx, R>(a_idx0, a_idx1, a_vec));
    }

/*--------------------------------------------------------------------*/
/// Get pointer to the const linear data
/** \return             Pointer to the const data
 *//*-----------------------------------------------------------------*/

  HOSTDEVICE const T* data() const noexcept
    {
      return shape_array::data() + offset();
    }

/*--------------------------------------------------------------------*/
/// Get pointer to the modifiable linear data
/** \return             Pointer to the modifiable data
 *//*-----------------------------------------------------------------*/

  HOSTDEVICE T* data() noexcept
    {
      return shape_array::data() + offset();
    }

/*--------------------------------------------------------------------*/
/// Begin iterator
/**
 *//*-----------------------------------------------------------------*/

  HOSTDEVICE pointer
  begin() noexcept
    {
      return data();
    }

/*--------------------------------------------------------------------*/
/// Constant begin iterator
/** 
 *//*-----------------------------------------------------------------*/

  HOSTDEVICE const_pointer
  begin() const noexcept
    {
      return data();
    }

/*--------------------------------------------------------------------*/
/// Constant begin iterator
/**
 *//*-----------------------------------------------------------------*/

  HOSTDEVICE const_pointer
  cbegin() const noexcept
    {
      return data();
    }

/*--------------------------------------------------------------------*/
/// End iterator
/**
 *//*-----------------------------------------------------------------*/

  HOSTDEVICE pointer
  end() noexcept
    {
      return data() + shape_array::size();
    }

/*--------------------------------------------------------------------*/
/// Constant end iterator
/**
 *//*-----------------------------------------------------------------*/

  HOSTDEVICE const_pointer
  end() const noexcept
    {
      return data() + shape_array::size();
    }

/*--------------------------------------------------------------------*/
/// Constant end iterator
/**
 *//*-----------------------------------------------------------------*/

  HOSTDEVICE const_pointer
  cend() const noexcept
    {
      return data() + shape_array::size();
    }

#ifdef CH_GPU
/*--------------------------------------------------------------------*/
/// Synchronous copy from host to device
/**
 *//*-----------------------------------------------------------------*/

  template <typename _alloc_type = alloc_type,
            std::enable_if_t<
              std::is_same<_alloc_type, cuda_alloc_type>::value, int> = 0>
  void
  copyToDevice() const noexcept
    {
      alloc_type::copyToDevice(data(), shape_array::size());
    }

/*--------------------------------------------------------------------*/
/// Asynchronous copy from host to device
/** \param[in]  a_stream
 *                      Stream identifier
 *//*-----------------------------------------------------------------*/

  template <typename _alloc_type = alloc_type,
            std::enable_if_t<
              std::is_same<_alloc_type, cuda_alloc_type>::value, int> = 0>
  void
  copyToDeviceAsync(CUstream a_stream = CH_Cuda::c_defaultStream) const noexcept
    {
      alloc_type::copyToDeviceAsync(data(), shape_array::size(), a_stream);
    }

/*--------------------------------------------------------------------*/
/// Synchronous copy from device to host
/**
 *//*-----------------------------------------------------------------*/

  template <typename _alloc_type = alloc_type,
            std::enable_if_t<
              std::is_same<_alloc_type, cuda_alloc_type>::value, int> = 0>
  void
  copyToHost() noexcept
    {
      alloc_type::copyToHost(data(), shape_array::size());
    }

/*--------------------------------------------------------------------*/
/// Asynchronous copy from device to host
/** \param[in]  a_stream
 *                      Stream identifier
 *//*-----------------------------------------------------------------*/

  template <typename _alloc_type = alloc_type,
            std::enable_if_t<
              std::is_same<_alloc_type, cuda_alloc_type>::value, int> = 0>
  void
  copyToHostAsync(CUstream a_stream = CH_Cuda::c_defaultStream) noexcept
    {
      alloc_type::copyToHostAsync(data(), shape_array::size(), a_stream);
    }

/*--------------------------------------------------------------------*/
/// Configure for GPU
/** This should only be used to prepare for copying objects to the GPU
 *//*-----------------------------------------------------------------*/

  template <typename _alloc_type = alloc_type,
            std::enable_if_t<
              std::is_same<_alloc_type, cuda_alloc_type>::value, int> = 0>
  void
  configureForGPU() noexcept
    {
      CH_assert(allocOn_gpu());
      shape_array::resetData(alloc_type::devicePtr() - offset());
    }
#endif


/*==============================================================================
 * Internal policies
 *============================================================================*/

private:

/*--------------------------------------------------------------------*/
/// Static assertions which should be called for each building constructor
/** 
 *//*-----------------------------------------------------------------*/

  HOSTDEVICE constexpr void
  staticChecks() const noexcept
    {
      static_assert(std::is_unsigned<size_type>::value,
                    "size_type must be unsigned");
      static_assert(std::is_signed<difference_type>::value,
                    "difference_type must be signed");
      static_assert((static_cast<difference_type>(
                       static_cast<size_type>(c_min_size*16) &
                       (~static_cast<size_type>(0) << 4)) >> 4) ==
                    c_min_size,
                    "offset shifts failed for c_min_size");
      static_assert((static_cast<difference_type>(
                       static_cast<size_type>(
                         static_cast<difference_type>(-1)*16) &
                       (~static_cast<size_type>(0) << 4)) >> 4) ==
                    static_cast<difference_type>(-1),
                    "offset shifts failed for -1");
      static_assert((static_cast<difference_type>(
                       c_max_size*16 &
                       (~static_cast<size_type>(0) << 4)) >> 4) ==
                    c_max_size,
                    "offset shifts failed for c_max_size");
      // If one of the above failed, perhaps the architecture is strange?
      // Uncomment the following for more information:
      // static_assert(static_cast<size_type>(-1) ==
      //               ~static_cast<size_type>(0) ||
      //               static_cast<size_type>(-1) ==
      //               (~static_cast<size_type>(0) << 1),
      //               "integer format not 1 or 2s complement");
      // static_assert((static_cast<difference_type>(-2) >> 1) ==
      //               static_cast<difference_type>(-1),
      //               ">> (signed right shift) is not arithmetic");
    }

/*--------------------------------------------------------------------*/
/// Get allocator (unlike standard containers, this is modifiable)
/** 
 *//*-----------------------------------------------------------------*/

  HOSTDEVICE alloc_type&
  get_allocator() noexcept
    {
      return *static_cast<alloc_type*>(this);
    }

//--Accessing parameters

/*--------------------------------------------------------------------*/
/// Set allocBy
/** 
 *//*-----------------------------------------------------------------*/

  HOSTDEVICE void
  setAllocBy(const AllocBy a_allocBy) noexcept
    {
      m_param =
        (m_param & ~static_cast<size_type>(0x3)) |
        static_cast<size_type>(a_allocBy);
    }

/*--------------------------------------------------------------------*/
/// Set offset
/** 
 *//*-----------------------------------------------------------------*/

  HOSTDEVICE void
  setOffset(const difference_type a_offset) noexcept
    {
      CH_assert(a_offset >= c_min_size && a_offset <= c_max_size);
      m_param =
        (m_param & static_cast<size_type>(0xF)) |
        static_cast<size_type>(a_offset*16);
    }

/*--------------------------------------------------------------------*/
/// Get allocBy
/** 
 *//*-----------------------------------------------------------------*/

  HOSTDEVICE AllocBy
  allocBy() const noexcept
    {
      return static_cast<AllocBy>(m_param & static_cast<size_type>(0x3));
    }

/*--------------------------------------------------------------------*/
/// Get allocOn
/** 
 *//*-----------------------------------------------------------------*/

  HOSTDEVICE unsigned
  allocOn() const noexcept
    {
      return m_param & static_cast<size_type>(AllocOn::all);
    }

/*--------------------------------------------------------------------*/
/// Is there an allocation on the CPU
/** 
 *//*-----------------------------------------------------------------*/

  HOSTDEVICE bool
  allocOn_cpu() const noexcept
    {
      return m_param & static_cast<size_type>(AllocOn::cpu);
    }

/*--------------------------------------------------------------------*/
/// Is there an allocation on the GPU
/** 
 *//*-----------------------------------------------------------------*/

  HOSTDEVICE bool
  allocOn_gpu() const noexcept
    {
      return m_param & static_cast<size_type>(AllocOn::gpu);
    }

/*--------------------------------------------------------------------*/
/// Get offset
/** 
 *//*-----------------------------------------------------------------*/

  HOSTDEVICE difference_type
  offset() const noexcept
    {
      // Requires arithmetic right shift for negative values
      return static_cast<difference_type>(m_param) >> 4;
    }

/*--------------------------------------------------------------------*/
/// Clear the offset and allocBy (allocOn is untouched)
/** 
 *//*-----------------------------------------------------------------*/

  HOSTDEVICE void
  clearOffsetParam() noexcept
    {
      m_param &= static_cast<size_type>(AllocOn::all);
    }

//--Construction of elements

  /*
    Can we/should we do uninitialized if POD?  The answer, as I understand it,
    is no.  Placement new must be used to start the lifetime of array elements,
    even if POD.
  */

/*--------------------------------------------------------------------*/
/// Construct an array and initilalize (with constructor arguments)
/** 
 *//*-----------------------------------------------------------------*/

  template <unsigned NDim, typename... DimArgs, std::size_t... I>
  HOSTDEVICE void
  constructArray(std::index_sequence<I...>,
                 std::true_type,  // Have arguments (yes)
                 std::true_type,  // Always true if have arguments
                 const DimArgs&... a_args)
    {
      static_assert(std::is_class<T>::value || sizeof...(DimArgs) - NDim == 1,
                    "Too many arguments for a non-class type");
#ifdef __CUDACC__
      // The problem is using std::forward_as_tuple is not allowed
      static_assert(NDim == sizeof...(DimArgs),
                    "Constructor arguments not supported in device code");
      pointer p = data();
      for (size_type n = shape_array::size(); n--;)
        {
          ::new((void *)p++) T();
        }
#else
      const auto fullTuple = std::forward_as_tuple(a_args...);
      pointer p = data();
      auto& alloc = get_allocator();
      for (size_type n = shape_array::size(); n--;)
        {
          AllocTr::construct(alloc,
                             p++,
                             std::get<I+NDim>(fullTuple)...);
        }
#endif
    }

/*--------------------------------------------------------------------*/
  /// Construct an array and initilalize (no arguments)
/** 
 *//*-----------------------------------------------------------------*/

  template <unsigned NDim, typename... DimArgs, std::size_t... I>
  HOSTDEVICE void
  constructArray(std::index_sequence<I...>,
                 std::false_type,  // Have arguments (no)
                 std::true_type,   // T == Real (yes)
                 const DimArgs&... a_args)
    {
      pointer p = data();
      auto& alloc = get_allocator();
      for (size_type n = shape_array::size(); n--;)
        {
#ifdef __CUDACC__
          ::new((void *)p++) T(c_BaseFabRealSetVal);
#else
          AllocTr::construct(alloc, p++, c_BaseFabRealSetVal);
#endif
        }
    }

/*--------------------------------------------------------------------*/
  /// Construct an array and initilalize (no arguments)
/** 
 *//*-----------------------------------------------------------------*/

  template <unsigned NDim, typename... DimArgs, std::size_t... I>
  HOSTDEVICE void
  constructArray(std::index_sequence<I...>,
                 std::false_type,  // Have arguments (no)
                 std::false_type,  // T == Real (no)
                 const DimArgs&... a_args)
    {
      pointer p = data();
      auto& alloc = get_allocator();
      for (size_type n = shape_array::size(); n--;)
        {
#ifdef __CUDACC__
          ::new((void *)p++) T;
#else
          AllocTr::construct(alloc, p++);
#endif
        }
    }

//--Destruction of elements

/*--------------------------------------------------------------------*/
  /// Destroy an array invoking destructor on elements
/** 
 *//*-----------------------------------------------------------------*/

  HOSTDEVICE void
  destroyArray(std::true_type)
    {
      pointer p = data();
      auto& alloc = get_allocator();
      for (size_type n = shape_array::size(); n--;)
        {
#ifdef __CUDACC__
          p++->~T();
#else
          AllocTr::destroy(alloc, p++);
#endif
        }
    }

/*--------------------------------------------------------------------*/
  /// Destroy array (no need to call destructor, i.e., not class type)
/** 
 *//*-----------------------------------------------------------------*/

  HOSTDEVICE void
  destroyArray(std::false_type)
    { }


/*==============================================================================
 * Member data
 *============================================================================*/

  size_type m_param;                  ///< Pointer offset caused by non-zero
                                      ///< lower bound allocated in bits 4-end,
                                      ///< allocOn in bits 2 and 3, and
                                      ///< allocBy in bits 0 and 1
};

}  // namespace dyn

#include "BaseNamespaceFooter.H"

#ifdef CH_GPU
#include "CudaDriver.H"

#include "BaseNamespaceHeader.H"

namespace CH_Cuda
{

/*--------------------------------------------------------------------*/
/// Specialization of Converting for dyn::Array
/** This alters the data pointer to point to GPU memory when an array
 *  is passed as a parameter to a GPU kernel.
 *//*-----------------------------------------------------------------*/

template <typename T,
          unsigned R,
          typename IIx,
          typename Ord,
          typename Cfg>
struct Converting<dyn::Array<T, R, IIx, Ord, Cfg, CUDAArrayAlloc>&>
{
  using type = dyn::Array<T, R, IIx, Ord, Cfg, CUDAArrayAlloc>;
  static type builder(type& a_arg)
    {
      type local(copy_to_device{}, a_arg);
      // std::cout << "DID conversion: " << local.begin() << std::endl;
      return local;
    }
};

}  // namespace CH_Cuda

#include "BaseNamespaceFooter.H"

#endif  /* defined CH_GPU */

#endif  /* ! defined _DYNARRAY_H_ */
