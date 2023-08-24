#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _ARRAYALLOCATOR_H_
#define _ARRAYALLOCATOR_H_


/******************************************************************************/
/**
 *  \file
 *
 *  \brief Allocators and implementations for array containers such as BaseFab
 *
 *  Allocators are intended to be swappable.  The array implementation is
 *  intended to be a stable object that interfaces with the allocator to manage
 *  memory for the container.
 *
 *  \note
 *  <ul>
 *    <li> Lower case _cpu and _gpu are used in names since it is easier to
 *         distinguish between c and g than C and G
 *    <li> These definitions can be included in objects used on the GPU.
 *         However, very few of the functions are supported on the GPU, these
 *         are annotated with HOSTDEVICE.  Functions that are not annotated can
 *         only execute on the host.
 *  </ul>
 *
 *//***************************************************************************/

#include <cstring>
#include <memory>

#include "memtrack.H"
#include "memusage.H"

#ifdef CH_GPU
#include "CH_Cuda.H"
#endif

/*
  Notes:
    - tag __CUDACC__ means the file is processed with nvcc.  The class is
      modified whether on the cpu or gpu.
    - tag __CUDA_ARCH__ is defined for device code.  This only works in
      functions with the __device__ qualifier.  It cannot be used in the
      class declaration.
*/

#undef HOSTDEVICE
#ifdef __CUDACC__
#define HOSTDEVICE __host__ __device__
#else
#define HOSTDEVICE
#endif

#include "BaseNamespaceHeader.H"

// Defined by testing for verbose output
#ifdef DEBUGARRAYALLOC
  #define ARRAYALLOCDBG(x) x
#else
  #define ARRAYALLOCDBG(x) (void)0
#endif

// Define this to allocate memory with calloc instead of malloc
//#define USE_CALLOC

constexpr Real c_BaseFabRealSetVal = BASEFAB_REAL_SETVAL;

/// Record for tracking memory use
enum class ArrayClassIndex
{
  unknown        = 0,
  Array_impl     = 1,
  BaseFab        = 2,
  CHArray        = 3,
  HashIVS        = 4,
  HashGraphTable = 5,
  DynArray       = 6
};

/// Track source of memory for array (can only use bits 0 & 1)
enum class AllocBy : unsigned
{
  none  = 0,                        ///< Undefined
  array = 1,                        ///< Data allocated by new[]
  alias = 2,                        ///< Alias to constructed data
  raw   = 3                         ///< Alias to raw memory
};

/// Devices with an allocation (can only use bits 2 & 3)
enum class AllocOn : unsigned
{
  none      = 0,
  cpu       = 0x4,
  gpu       = 0x8,
  all       = 0xC
};

// Allow unions of AllocOn
constexpr unsigned operator |(const AllocOn a_x, const AllocOn a_y) noexcept
{
  return (unsigned)a_x | (unsigned)a_y;
}

constexpr unsigned operator |(const unsigned a_x, const AllocOn a_y) noexcept
{
  return a_x | (unsigned)a_y;
}

constexpr unsigned operator |(const AllocOn a_x, const unsigned a_y) noexcept
{
  return (unsigned)a_x | a_y;
}


/*============================================================================*/
/// Traits for array container classes
/*============================================================================*/

template <ArrayClassIndex I>
struct ArrayClassTr
{
  static const char* name()
    {
      return "Unknown";
    }
};

template <>
struct ArrayClassTr<ArrayClassIndex::Array_impl>
{
  static const char* name()
    {
      return "Array_impl";
    }
};

template <>
struct ArrayClassTr<ArrayClassIndex::BaseFab>
{
  static const char* name()
    {
      return "BaseFab";
    }
};

template <>
struct ArrayClassTr<ArrayClassIndex::CHArray>
{
  static const char* name()
    {
      return "CHArray";
    }
};

template <>
struct ArrayClassTr<ArrayClassIndex::HashIVS>
{
  static const char* name()
    {
      return "HashIVS";
    }
};

template <>
struct ArrayClassTr<ArrayClassIndex::HashGraphTable>
{
  static const char* name()
    {
      return "CH_Hash::StaticTable";
    }
};

template <>
struct ArrayClassTr<ArrayClassIndex::DynArray>
{
  static const char* name()
    {
      return "dyn::Array";
    }
};


/*******************************************************************************
 */
//  Default array allocator
/*
 ******************************************************************************/

template<typename T, ArrayClassIndex I = ArrayClassIndex::unknown>
class DefaultArrayAlloc;


/*******************************************************************************
 */
/// Default array allocator specialized for void
/*
 ******************************************************************************/

template <ArrayClassIndex I>
class DefaultArrayAlloc<void, I> : public std::allocator<void>
{
public:

  /// Propagate types from std::allocator
  using size_type       = typename std::allocator<void>::size_type;
  using difference_type = typename std::allocator<void>::difference_type;
  using pointer         = typename std::allocator<void>::pointer;
  using const_pointer   = typename std::allocator<void>::const_pointer;
  using value_type      = typename std::allocator<void>::value_type;

  template<typename U>
  struct rebind
  {
    using other = DefaultArrayAlloc<U, I>;
  };

  using propagate_on_container_copy_assignment = std::false_type;
  using propagate_on_container_move_assignment = std::false_type;
  using propagate_on_container_swap = std::false_type;
  using is_always_equal = std::true_type;

  /// Non standard flags
  using allocates_on_cpu = std::true_type;
  using allocates_on_gpu = std::false_type;
  static constexpr unsigned allocates_on = static_cast<unsigned>(AllocOn::cpu);
};


/*******************************************************************************
 */
/// Default array allocator
/**
 *  Default allocators used by array classes is based on std::allocator.
 *  This allocator is stateless and has true equivalence.  This means
 *  DefaultArrayAlloc A can deallocate the memory allocated by
 *  DefaultArrayAlloc B.  The implication of this is that the allocator does
 *  not need to be moved or copied with the underlying container.
 *
 *  We customize allocation in this class and also support registering memory
 *  with the memory tracker
 *
 *//*+*************************************************************************/

template <typename T, ArrayClassIndex I>
class DefaultArrayAlloc : public std::allocator<T>
{

/*==============================================================================
 * Types
 *============================================================================*/

public:

  // Propagate types from std::allocator
  using size_type       = typename std::allocator<T>::size_type;
  using difference_type = typename std::allocator<T>::difference_type;
  using pointer         = typename std::allocator<T>::pointer;
  using const_pointer   = typename std::allocator<T>::const_pointer;
  using reference       = typename std::allocator<T>::reference;
  using const_reference = typename std::allocator<T>::const_reference;
  using value_type      = typename std::allocator<T>::value_type;

  /// Bind this allocator to another type (ArrayClassIndex propagates)
  template<typename U>
  struct rebind
  {
    using other = DefaultArrayAlloc<U, I>;
  };

  // No propagation is necessary since the default allocator is stateless and
  // uses standard methods for allocation and deallocation.
  using propagate_on_container_copy_assignment = std::false_type;
  using propagate_on_container_move_assignment = std::false_type;
  using propagate_on_container_swap = std::false_type;
  using is_always_equal = std::true_type;

  /// Use this allocator on container copy
  static DefaultArrayAlloc
  select_on_container_copy_construction(const DefaultArrayAlloc& a_alloc)
    {
      return DefaultArrayAlloc{};
    }

  /// Void types (squirrely way of getting a void*)
  using void_allocator = typename rebind<void>::other;

  /// Non standard flags
  using allocates_on_cpu = std::true_type;
  using allocates_on_gpu = std::false_type;
  static constexpr unsigned allocates_on = static_cast<unsigned>(AllocOn::cpu);


/*==============================================================================
 * Constructors and destructors
 *============================================================================*/

public:

  /// Default constructor
  DefaultArrayAlloc()
    {
      ARRAYALLOCDBG(pout() << "DefaultArrayAllocator: constructor\n");
    }

  // Use synthesized copy, move, copy assignment, move assignment, and
  // destructor.


/*==============================================================================
 * Member functions
 *============================================================================*/

public:  // Would rather have these protected but that does not play well with
         // alloc_traits

  /// Allocate uninitialized memory (no construction)
  pointer allocate(size_type                              a_size,
                   typename void_allocator::const_pointer a_hint = nullptr);

  /// Deallocate memory (no destruction)
  void deallocate(pointer a_addr, const size_type a_size);

#ifdef CH_ALIGN
  /// Set alignment (default 2*sizeof(void*) or 16 bytes on 64-bit machines).
  static void setAlignment(const size_type a_alignment)
    {
      s_alignment = a_alignment;
    }
#endif


/*==============================================================================
 * Data members
 *============================================================================*/

protected:

#ifdef CH_ALIGN
  static size_type s_alignment;       ///< Alignment for allocations.  Default
                                      ///< is 2*sizeof(void*) or 16 bytes on
                                      ///< 64-bit machines
#endif

#ifdef CH_USE_MEMORY_TRACKING
  static MemtrackRegistry* s_memtrackRegistry;
                                      ///< A memory tracking registry for
                                      ///< counting the total number of bytes
                                      ///< allocated.  The actual registry
                                      ///< entry is unique to the T and I
                                      ///< template parameters
#endif
};


/*******************************************************************************
 */
///  CUDA array allocator
/**
 *   Allocates memory on both host and GPU for the array
 *
 *//*+*************************************************************************/

#ifdef CH_GPU

template<typename T, ArrayClassIndex I = ArrayClassIndex::unknown>
class CUDAArrayAlloc;


/*******************************************************************************
 */
/// CUDA array allocator specialized for void
/*
 ******************************************************************************/

template <ArrayClassIndex I>
class CUDAArrayAlloc<void, I> : public std::allocator<void>
{
public:

  /// Propagate types from std::allocator
  using size_type       = typename std::allocator<void>::size_type;
  using difference_type = typename std::allocator<void>::difference_type;
  using pointer         = typename std::allocator<void>::pointer;
  using const_pointer   = typename std::allocator<void>::const_pointer;
  using value_type      = typename std::allocator<void>::value_type;

  template<typename U>
  struct rebind
  {
    using other = CUDAArrayAlloc<U, I>;
  };

  using propagate_on_container_copy_assignment = std::true_type;
  using propagate_on_container_move_assignment = std::true_type;
  using propagate_on_container_swap = std::true_type;
  using is_always_equal = std::true_type;

  /// Non standard flags
  using allocates_on_cpu = std::true_type;
  using allocates_on_gpu = std::true_type;
#ifdef __CUDACC__
  static constexpr unsigned allocates_on = static_cast<unsigned>(AllocOn::gpu);
#else
  static constexpr unsigned allocates_on = static_cast<unsigned>(AllocOn::all);
#endif
};


/*******************************************************************************
 */
/// CUDA array allocator
/**
 *  Supports simultaneous allocations on both the host and device
 *
 *//*+*************************************************************************/

template <typename T, ArrayClassIndex I>
class CUDAArrayAlloc
#ifndef __CUDACC__
// nvcc will complain about allocator() being a host function
  : public std::allocator<T>
#endif
{

/*==============================================================================
 * Types
 *============================================================================*/

public:

  // Propagate types from std::allocator
  using size_type       = typename std::allocator<T>::size_type;
  using difference_type = typename std::allocator<T>::difference_type;
  using pointer         = typename std::allocator<T>::pointer;
  using const_pointer   = typename std::allocator<T>::const_pointer;
  using reference       = typename std::allocator<T>::reference;
  using const_reference = typename std::allocator<T>::const_reference;
  using value_type      = typename std::allocator<T>::value_type;

  /// Bind this allocator to another type (ArrayClassIndex propagates)
  template<typename U>
  struct rebind
  {
    using other = CUDAArrayAlloc<U, I>;
  };

  // Propagation is necessary since we store state in this class
  using propagate_on_container_copy_assignment = std::true_type;
  using propagate_on_container_move_assignment = std::true_type;
  using propagate_on_container_swap = std::true_type;
  // Allocators can allocate, deallocate each others memory
  using is_always_equal = std::true_type;

  /// Copy allocator on container copy
  HOSTDEVICE static CUDAArrayAlloc
  select_on_container_copy_construction(const CUDAArrayAlloc& a_alloc)
    {
      return a_alloc;
    }

  /// Void types (squirrely way of getting a void*)
  using void_allocator = typename rebind<void>::other;

  /// Non standard flags
  using allocates_on_cpu = std::true_type;
  using allocates_on_gpu = std::true_type;
#ifdef __CUDACC__
  static constexpr unsigned allocates_on = static_cast<unsigned>(AllocOn::gpu);
#else
  static constexpr unsigned allocates_on = static_cast<unsigned>(AllocOn::all);
#endif


/*==============================================================================
 * Constructors and destructors
 *============================================================================*/

public:

  /// Constructor
  HOSTDEVICE CUDAArrayAlloc()
    {
      ARRAYALLOCDBG(pout() << "CUDAArrayAllocator: constructor\n");
    }

  // Use synthesized copy, move, copy assignment, move assignment, and
  // destructor.


/*==============================================================================
 * Member functions
 *============================================================================*/

public:  // Would rather have these protected but that does not play well with
         // alloc_traits

  /// Allocate uninitialized memory (no construction)
  HOSTDEVICE pointer allocate(
    size_type                              a_size,
    typename void_allocator::const_pointer a_hint);

  /// Deallocate memory (no destruction)
  HOSTDEVICE void deallocate(pointer a_addr, const size_type a_size);

#ifdef CH_ALIGN
  /// Set alignment (default 2*sizeof(void*) or 16 bytes on 64-bit machines).
  static void setAlignment(const size_type a_alignment)
    {
      s_alignment = a_alignment;
    }
#endif

  /// Set whether to page-lock host memory if GPU is not used (default false)
  static void setAllocPageLocked(const bool a_allocPageLocked)
    {
      s_allocPageLocked = a_allocPageLocked;
    }

  /// Is the argument data configured for the GPU?
  constexpr bool isConfiguredFor_gpu(const pointer a_data) const noexcept;

  /// Synchronous copy from host to device
  void copyToDevice(const_pointer   a_hostPtr,
                    const size_type a_size) const noexcept;

  /// Asynchronous copy from host to device
  void copyToDeviceAsync(
    const_pointer   a_hostPtr,
    const size_type a_size,
    CUstream        a_stream = CH_Cuda::c_defaultStream) const noexcept;

  /// Synchronous copy from device to host
  void copyToHost(pointer         a_hostPtr,
                  const size_type a_size) noexcept;

  /// Asynchronous copy from device to host
  void copyToHostAsync(
    pointer         a_hostPtr,
    const size_type a_size,
    CUstream        a_stream = CH_Cuda::c_defaultStream) noexcept;

  /// Get the device pointer (as type pointer)
  HOSTDEVICE pointer devicePtr() const noexcept;

  /// Reset data pointer on GPU
  HOSTDEVICE void setDeviceAddress(
    const CH_Cuda::DevicePointer a_addr) noexcept;


/*==============================================================================
 * Data members
 *============================================================================*/

private:

  CH_Cuda::DevicePointer m_devicePtr; ///< Pointer to memory allocated on the
                                      ///< device

protected:

  static bool s_allocPageLocked;      ///< Considered for host memory when not
                                      ///< using GPU (i.e., m_allocGPU == false.
                                      ///< If, m_allocGPU == true, the host
                                      ///< memory is always page-locked.
                                      ///< Default is F.
                                      ///< T - allocate page-locked memory on
                                      ///<     the host, meaning the memory is
                                      ///<     locked in RAM.  This memory can
                                      ///<     still be efficiently move to the
                                      ///<     GPU
                                      ///< F - allocate pageable memory with
                                      ///<     malloc.  This memory still be
                                      ///<     moved to the GPU but with less
                                      ///<     performance.

#ifdef CH_ALIGN
  static size_type s_alignment;       ///< Alignment for allocations.  Default
                                      ///< is 2*sizeof(void*) or 16 bytes on
                                      ///< 64-bit machines
#endif

#ifdef CH_USE_MEMORY_TRACKING
  static MemtrackRegistry* s_memtrackRegistryHost;
                                      ///< A memory tracking registry for
                                      ///< counting the total number of bytes
                                      ///< allocated.  The actual registry
                                      ///< entry is unique to the T and I
                                      ///< template parameters
  static MemtrackRegistry* s_memtrackRegistryDevice;
                                      ///< Same as above but tracks memory on
                                      ///< device
#endif
};

#endif  /* defined CH_GPU */


/*******************************************************************************
 */
/// Array implementation
/**
 *  Interfaces with the allocator to provide a high-level management of
 *  memory used for arrays.  While allocators can be changed, the intention
 *  is that this class remains constant.
 *
 *  This class is low-level and only intended for use by other containers.  For
 *  a general-purpose dynamic array, use class dyn::Array instead.
 *
 *  \tparam     T       Type of element
 *  \tparam     Alloc   Allocator
 *
 *//*+*************************************************************************/

template <typename T, typename Alloc>
class Array_impl : public Alloc
{


/*==============================================================================
 * Types
 *============================================================================*/

public:

  // Allocator_traits can provide defaults for missing features in the allocator
  using AllocTr         = std::allocator_traits<Alloc>;

  // Types
  using size_type       = typename AllocTr::size_type;
  using difference_type = typename AllocTr::difference_type;
  using pointer         = T*;
  using const_pointer   = const T*;
  using reference       = T&;
  using const_reference = const T&;
  using value_type      = T;

  static constexpr size_type c_max_size =
    ~static_cast<size_type>(0) >> 4;  // size_type must be unsigned


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

public:

  /// Default constructor (default AllocOn varies)
  HOSTDEVICE Array_impl() noexcept
    :
#if defined(CH_GPU) && defined(__CUDACC__)
    Array_impl(AllocOn::gpu)         // Delegate to next constructor
#else
    Array_impl(Alloc::allocates_on)
#endif
    { }

#ifdef CH_GPU
  /// Perform default initialization
  DEVICE void INITIALIZEshared() noexcept;
#endif

  /// Constructor (takes single AllocOn)
  HOSTDEVICE Array_impl(const AllocOn a_allocOn) noexcept
    :
    Array_impl((unsigned)a_allocOn)  // Delegate to next constructor
    { }

  /// Constructor (takes union of AllocOn)
  HOSTDEVICE Array_impl(const unsigned a_allocOn) noexcept;

  // No copies allowed
  Array_impl(const Array_impl&) = delete;
  Array_impl& operator=(const Array_impl&) = delete;

  /// Move construction is allowed
  HOSTDEVICE Array_impl(Array_impl&& a_impl) noexcept;

  /// Move assignment is allowed
  HOSTDEVICE Array_impl& operator=(Array_impl&& a_impl) noexcept;

  /// Destructor
  HOSTDEVICE ~Array_impl()
    { undefine(); }

#ifdef CH_GPU
  /// Copy constructor that swaps pointers to data on the GPU
  Array_impl(CH_Cuda::copy_to_device, const Array_impl& a_impl) noexcept;

  /// Configure the array for use on the GPU
  void configureForGPU() noexcept;
#endif


/*==============================================================================
 * Members functions
 *============================================================================*/

//--Defines

  /// Allocate and construct T
  template <typename... Args>
  HOSTDEVICE void define(const size_type a_size, Args&&... a_args);

  /// Alias constructed T
  HOSTDEVICE void defineAlias(pointer a_addr, const size_type a_size);

  /// Construct T on raw memory
  template <typename... Args>
  HOSTDEVICE void defineRaw(void*           a_addr,
                            const size_type a_size,
                            Args&&...       a_args);

  /// Deallocate and destroy T
  HOSTDEVICE void undefine();

//--Miscellaneous

  /// Class is defined (but unusable if number of elements is zero)
  HOSTDEVICE bool isDefined() const noexcept;

  /// Class is defined an usable
  HOSTDEVICE bool isUsable() const noexcept;

//--Access

  /// Get allocator (unlike standard containers, this is modifiable)
  HOSTDEVICE Alloc& get_allocator() noexcept;

  /// Size (number of array elements)
  HOSTDEVICE size_type size() const noexcept;

  /// Allocation type
  HOSTDEVICE AllocBy allocBy() const noexcept;

  /// Get all processors with allocation
  HOSTDEVICE unsigned allocOn() const noexcept;

  /// Is there an allocation on the CPU
  HOSTDEVICE bool allocOn_cpu() const noexcept;

  /// Is there an allocation on the GPU
  HOSTDEVICE bool allocOn_gpu() const noexcept;

  /// Index
  HOSTDEVICE const_reference operator[](const size_type a_idx) const noexcept;

  /// Index
  HOSTDEVICE reference operator[](const size_type a_idx) noexcept;

  /// Access to the data pointer
  HOSTDEVICE const_pointer dataPtr() const noexcept;

  /// Access to the data pointer
  HOSTDEVICE pointer dataPtr() noexcept;

  /// Access to the data pointer
  HOSTDEVICE const_pointer data() const noexcept;

  /// Access to the data pointer
  HOSTDEVICE pointer data() noexcept;

  /// Begin iterator
  HOSTDEVICE pointer begin() noexcept;

  /// Constant begin iterator
  HOSTDEVICE const_pointer begin() const noexcept;

  /// Constant begin iterator
  HOSTDEVICE const_pointer cbegin() const noexcept;

  /// End iterator
  HOSTDEVICE pointer end() noexcept;

  /// Constant end iterator
  HOSTDEVICE const_pointer end() const noexcept;

  /// Constant end iterator
  HOSTDEVICE const_pointer cend() const noexcept;

#ifdef CH_GPU
  /// Synchronous copy from host to device
  void copyToDevice() const noexcept;

  /// Asynchronous copy from host to device
  void copyToDeviceAsync(
    CUstream a_stream = CH_Cuda::c_defaultStream) const noexcept;

  /// Synchronous copy from device to host
  void copyToHost() noexcept;

  /// Asynchronous copy from device to host
  void copyToHostAsync(CUstream a_stream = CH_Cuda::c_defaultStream) noexcept;
#endif



/*==============================================================================
 * Internal policies
 *============================================================================*/

private:

  friend int testArrayAllocator1();

//--Sizes and flags

  /// Set size, allocOn, and allocBy
  HOSTDEVICE void size(const size_type a_size,
                       const unsigned  a_allocOn,
                       const AllocBy   a_allocBy) noexcept;

  /// Set size and allocBy
  HOSTDEVICE void size(const size_type a_size,
                       const AllocBy   a_allocBy) noexcept;

  /// Set size
  HOSTDEVICE void size(const size_type a_size) noexcept;

  /// Set allocBy
  HOSTDEVICE void setAllocBy(const AllocBy a_allocBy) noexcept;

//--Movement of the allocator

  /// Move allocator
  HOSTDEVICE void move_allocator(std::true_type, Alloc& a_alloc) noexcept;

  /// Unnecessary to move allocator (keep using current)
  HOSTDEVICE void move_allocator(std::false_type, Alloc& a_alloc) noexcept;

//--Construction of elements

  /// Construct an array invoking constructor on elements
  template <typename... Args>
  HOSTDEVICE void constructArray(std::true_type, Args&&... a_args);

  /// Construct an array and initilalize (not class type and at least 1 arg)
  template <typename Arg0, typename... Args>
  HOSTDEVICE void constructArray(std::false_type,
                                 Arg0&&    a_arg0,
                                 Args&&... a_args);

  /// Construct array as uninitialized (i.e., not class type and no args)
  /*  This method has a specialization if T is Real in which case the array
   *  is initialized to c_BaseFabRealSetVal */
  HOSTDEVICE void constructArray(std::false_type);

//--Destruction of elements

  /// Destroy an array invoking destructor on elements
  HOSTDEVICE void destroyArray(std::true_type);

  /// Destroy array (no need to call destructor, i.e., not class type)
  HOSTDEVICE void destroyArray(std::false_type)
    { }

protected:

//==Data members of implementation.  Note that some allocators may add data.

  pointer m_data;                     ///< The data pointer.

private:

  // This must never be accessed directly since it holds more than just the size
  size_type m_param;                  ///< Total number of elements
                                      ///< allocated in bits 4-end,
                                      ///< allocOn in bits 2 and 3, and
                                      ///< allocBy in bits 0 and 1
};

#include "ArrayAllocatorI.H"
#include "BaseNamespaceFooter.H"

#ifdef CH_GPU
#include "CudaDriver.H"

#include "BaseNamespaceHeader.H"

namespace CH_Cuda
{

/// Specialization of Converting for Array_impl
/** This alters the data pointer to point to GPU memory when an array is passed
 *  as a parameter to a GPU kernel.
 */
template <typename T, typename Alloc>
struct Converting<Array_impl<T, Alloc>&>
{
  using type = Array_impl<T, Alloc>;
  static type builder(type& a_arg)
    {
      type local(copy_to_device{}, a_arg);
      // std::cout << "DID conversion: " << local.begin() << std::endl;
      return local;
    }
};

}  // namespace CH_Cuda

#include "BaseNamespaceFooter.H"

#endif


#endif  /* ! defined _ARRAYALLOCATOR_H_ */
