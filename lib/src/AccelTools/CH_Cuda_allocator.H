#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _CH_CUDA_ALLOCATOR_H_
#define _CH_CUDA_ALLOCATOR_H_


/******************************************************************************/
/**
 *  \file CH_Cuda_allocator.H
 *
 *  \brief A special allocator for creating objects on both the host and device.
 *
 *  This is only intended for use on the host but makes a dual allocation of
 *  memory on both the host and device.  You must separately keep track of the
 *  dual pointers for deallocation.  Think of this as an augmented new and
 *  delete.
 *
 *  \note
 *  <ul>
 *    <li> This allocator can be used without the GPU.  It will allocate on the
 *         host using malloc and otherwise behave like standard new.
 *    <li> If using the GPU, allocations are always made on both devices.
 *    <li> The SymbolPair type that is returned has methods for copying data
 *         between the host and device pointers.
 *  </ul>
 *
 *//***************************************************************************/

#include <cstring>
#include <memory>

#include "memtrack.H"
#include "memusage.H"
#include "CH_Cuda.H"

#include "BaseNamespaceHeader.H"

namespace CH_Cuda
{

template<typename T>
class Alloc;


/*******************************************************************************
 */
/// CUDA allocator specialized for void
/*
 ******************************************************************************/

template <>
class Alloc<void>
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
    using other = Alloc<U>;
  };

  using propagate_on_container_copy_assignment = std::false_type;
  using propagate_on_container_move_assignment = std::false_type;
  using propagate_on_container_swap = std::false_type;
  using is_always_equal = std::true_type;

  /// Non standard flags
  using allocates_on_cpu = std::true_type;
  using allocates_on_gpu = std::true_type;
  // See ArrayAllocator.H for definitions
  static constexpr unsigned allocates_on = static_cast<unsigned>(0xC);
};


/*******************************************************************************
 */
/// CUDA allocator
/**
 *  Supports simultaneous allocations on both the host and device
 *
 *//*+*************************************************************************/

template <typename T>
class Alloc
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
    using other = Alloc<U>;
  };

  using propagate_on_container_copy_assignment = std::false_type;
  using propagate_on_container_move_assignment = std::false_type;
  using propagate_on_container_swap = std::false_type;
  // Allocators can allocate, deallocate each others memory
  using is_always_equal = std::true_type;

  /// Use this allocator on container copy
  static Alloc
  select_on_container_copy_construction(const Alloc& a_alloc)
    {
      return Alloc{};
    }

  /// Void types (squirrely way of getting a void*)
  using void_allocator = typename rebind<void>::other;

  /// Non standard flags
  using allocates_on_cpu = std::true_type;
  using allocates_on_gpu = std::true_type;
  // See ArrayAllocator.H for definitions
  static constexpr unsigned allocates_on = static_cast<unsigned>(0xC);


/*==============================================================================
 * Constructors and destructors
 *============================================================================*/

  // Use synthesized default, copy, move, copy assignment, move assignment, and
  // destructor.


/*==============================================================================
 * Member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
/// Allocate uninitialized memory (no construction)
/**
 *//*-----------------------------------------------------------------*/

  static SymbolPair<T>
  allocate() noexcept
    {
      constexpr size_type numBytes = sizeof(T);
      SymbolPair<T> ptrs;
      void* hostPtr;
#ifdef CH_GPU
      // Allocate page-locked memory on host
      checkCudaErrors(cuMemHostAlloc(&hostPtr, numBytes, 0));
      // Allocate on device
      checkCudaErrors(cuMemAlloc(&ptrs.device, numBytes));
#else
      hostPtr = malloc(numBytes);
      if (hostPtr == nullptr)
        {
          print_memory_line("Out of memory");
          MayDay::Error("Out of memory in CH_Cuda::Alloc::allocate");
        }
#endif
      ptrs.host = static_cast<pointer>(hostPtr);
#ifdef CH_USE_MEMORY_TRACKING
      s_memtrackRegistryHost->addBytes(numBytes);
#ifdef CH_GPU
      s_memtrackRegistryDevice->addBytes(numBytes);
#endif
#endif
      return ptrs;
    }

/*--------------------------------------------------------------------*/
/// Allocate uninitialized memory only on host (no construction)
/**
 *//*-----------------------------------------------------------------*/

  static SymbolPair<T>
  allocateHost() noexcept
    {
      constexpr size_type numBytes = sizeof(T);
      SymbolPair<T> ptrs;
      void* hostPtr;
#ifdef CH_GPU
      // Allocate page-locked memory on host
      checkCudaErrors(cuMemHostAlloc(&hostPtr, numBytes, 0));
#else
      hostPtr = malloc(numBytes);
      if (hostPtr == nullptr)
        {
          print_memory_line("Out of memory");
          MayDay::Error("Out of memory in CH_Cuda::Alloc::allocate");
        }
#endif
      ptrs.host = static_cast<pointer>(hostPtr);
#ifdef CH_USE_MEMORY_TRACKING
      s_memtrackRegistryHost->addBytes(numBytes);
#endif
      return ptrs;
    }

/*--------------------------------------------------------------------*/
/// Deallocate raw memory (no destruction)
/** Used with a temporary argument
 *//*-----------------------------------------------------------------*/

  static void
  deallocate(SymbolPair<T>&& a_ptrs) noexcept
    {
      if (a_ptrs.host != nullptr)
        {
#ifdef CH_GPU
          checkCudaErrors(cuMemFreeHost(a_ptrs.host));
#else
          free(a_ptrs.host);
#endif
#ifdef CH_USE_MEMORY_TRACKING
          s_memtrackRegistryHost->subtractBytes(sizeof(T));
#endif
        }
#ifdef CH_GPU
      if (a_ptrs.device != 0)
        {
          checkCudaErrors(cuMemFree(a_ptrs.device));
#ifdef CH_USE_MEMORY_TRACKING
          s_memtrackRegistryDevice->subtractBytes(sizeof(T));
#endif
        }
#endif
    }

/*--------------------------------------------------------------------*/
/// Deallocate raw memory (no destruction)
/** Used with an l-value.  Sets a_ptrs to { nullptr, 0 }
 *//*-----------------------------------------------------------------*/

  /// Deallocate raw memory (no destruction)
  static void
  deallocate(SymbolPair<T>& a_ptrs) noexcept
    {
      if (a_ptrs.host != nullptr)
        {
#ifdef CH_GPU
          checkCudaErrors(cuMemFreeHost(a_ptrs.host));
#else
          free(a_ptrs.host);
#endif
          a_ptrs.host = nullptr;
#ifdef CH_USE_MEMORY_TRACKING
          s_memtrackRegistryHost->subtractBytes(sizeof(T));
#endif
        }
#ifdef CH_GPU
      if (a_ptrs.device != 0)
        {
          checkCudaErrors(cuMemFree(a_ptrs.device));
          a_ptrs.device = 0;
#ifdef CH_USE_MEMORY_TRACKING
          s_memtrackRegistryDevice->subtractBytes(sizeof(T));
#endif
        }
#endif
    }

/*--------------------------------------------------------------------*/
/// Construct a single object in the host memory
/** Note: you must construct even POD type to begin their lifetime
 *//*-----------------------------------------------------------------*/

  template <typename... Args>
  static void
  construct(SymbolPair<T>& a_ptrs, Args&&... a_args)
    {
      CH_assert(a_ptrs.host);
      ::new(static_cast<void*>(a_ptrs.host)) T(std::forward<Args>(a_args)...);
    }

/*--------------------------------------------------------------------*/
/// Destroy a single object in the host memory
/**
 *//*-----------------------------------------------------------------*/

  static void
  destroy(SymbolPair<T> a_ptrs)
    {
      CH_assert(a_ptrs.host);
      a_ptrs.host->~T();
    }

/*--------------------------------------------------------------------*/
/// Allocates and constructs an object (similar to new)
/** Note: the object is not constructed in GPU memory and must be
 *  copied there.
 *  \tparam Args        Argument types for the constructor
 *  \param[in]  a_args  Arguments for the constructor
 *  \return             Pointers to both host and device memory
 *//*-----------------------------------------------------------------*/

  template <typename... Args>
  static SymbolPair<T>
  newPair(Args&&... a_args)
    {
      SymbolPair<T> ptrs = allocate();
      construct(ptrs, std::forward<Args>(a_args)...);
      return ptrs;
    }

/*--------------------------------------------------------------------*/
/// Allocates and constructs an object only on host (similar to new)
/** \tparam Args        Argument types for the constructor
 *  \param[in]  a_args  Arguments for the constructor
 *  \return             Pointer to host memory
 *//*-----------------------------------------------------------------*/

  template <typename... Args>
  static SymbolPair<T>
  newHost(Args&&... a_args)
    {
      SymbolPair<T> ptrs = allocateHost();
      construct(ptrs, std::forward<Args>(a_args)...);
      return ptrs;
    }

/*--------------------------------------------------------------------*/
/// Allocates and constructs an object (similar to new)
/** Note: the object is not constructed in GPU memory and must be
 *  copied there.
 *  \tparam Args        Argument types for the constructor
 *  \param[in]  a_args  Arguments for the constructor
 *  \return             Pointers to both host and device memory
 *//*-----------------------------------------------------------------*/

  template <typename... Args>
  static SymbolPair<T>
  newSelectHost(std::false_type, Args&&... a_args)
    {
      return newPair(std::forward<Args>(a_args)...);
    }

/*--------------------------------------------------------------------*/
/// Allocates and constructs an object only on host (similar to new)
/** \tparam Args        Argument types for the constructor
 *  \param[in]  a_args  Arguments for the constructor
 *  \return             Pointer to host memory
 *//*-----------------------------------------------------------------*/

  template <typename... Args>
  static SymbolPair<T>
  newSelectHost(std::true_type, Args&&... a_args)
    {
      return newHost(std::forward<Args>(a_args)...);
    }

/*--------------------------------------------------------------------*/
/// Allocates and constructs an object (similar to new)
/** Note: the object is not constructed in GPU memory and must be
 *  copied there.
 *  \tparam Args        Argument types for the constructor
 *  \param[in]  a_args  Arguments for the constructor
 *  \return             Pointers to both host and device memory
 *//*-----------------------------------------------------------------*/

  template <typename... Args>
  static SymbolPair<T>
  newSelectPair(std::true_type, Args&&... a_args)
    {
      return newPair(std::forward<Args>(a_args)...);
    }

/*--------------------------------------------------------------------*/
/// Allocates and constructs an object only on host (similar to new)
/** \tparam Args        Argument types for the constructor
 *  \param[in]  a_args  Arguments for the constructor
 *  \return             Pointer to host memory
 *//*-----------------------------------------------------------------*/

  template <typename... Args>
  static SymbolPair<T>
  newSelectPair(std::false_type, Args&&... a_args)
    {
      return newHost(std::forward<Args>(a_args)...);
    }

/*--------------------------------------------------------------------*/
/// Destroys and deallocates an object (similar to delete)
/** Used with a temporary argument
 *  \param[in]   a_ptrs Pointers to both host and device memory
 *//*-----------------------------------------------------------------*/

  static void
  deletePair(SymbolPair<T>&& a_ptrs)
    {
      destroy(a_ptrs);
      deallocate(a_ptrs);
    }

/*--------------------------------------------------------------------*/
/// Destroys and deallocates an object (similar to delete)
/** Used with an l-value.  Sets a_ptrs to { nullptr, 0 }
 *  \param[in]  a_ptrs  Pointers to both host and device memory
 *  \param[out] a_ptrs  Set to { nullptr, 0 }
 *//*-----------------------------------------------------------------*/

  static void
  deletePair(SymbolPair<T>& a_ptrs)
    {
      destroy(a_ptrs);
      deallocate(a_ptrs);
    }


/*==============================================================================
 * Data members
 *============================================================================*/

#ifdef CH_USE_MEMORY_TRACKING
  static MemtrackRegistry* s_memtrackRegistryHost;
                                      ///< A memory tracking registry for
                                      ///< counting the total number of bytes
                                      ///< allocated.  The actual registry
                                      ///< entry is unique to the T template
                                      ///< parameter
#ifdef CH_GPU
  static MemtrackRegistry* s_memtrackRegistryDevice;
                                      ///< Same as above but tracks memory on
                                      ///< device
#endif
#endif
};


/*******************************************************************************
 *
 * Class Alloc: static member definitions
 *
 ******************************************************************************/

//--Initialization of memtrackRegistry (this memory is deallocated by the
//--dumpmemoryatexit() routine in memtrack.cpp after reporting)

#ifdef CH_USE_MEMORY_TRACKING
// For host
template <typename T>
MemtrackRegistry* Alloc<T>::s_memtrackRegistryHost =
  new MemtrackRegistry(
    // E.g. for int, results in CH_Cuda::Alloc<i> using gcc
    std::string("CH_Cuda::Alloc<")
    + std::string((typeid(T)).name())
    + ">");
#ifdef CH_GPU
// For device
template <typename T>
MemtrackRegistry* Alloc<T>::s_memtrackRegistryDevice =
  new MemtrackRegistry(
    // E.g. for int, results in CH_Cuda::Alloc_GPU<i> using gcc
    std::string("CH_Cuda::Alloc_GPU<")
    + std::string((typeid(T)).name())
    + ">",
    MemtrackRegistry::Device::GPU);
#endif
#endif

}  // namespace CH_Cuda

#include "BaseNamespaceFooter.H"

#endif  /* ! defined _CH_CUDA_ALLOCATOR_H_ */
