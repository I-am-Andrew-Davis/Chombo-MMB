
/*******************************************************************************
 *
 * Class DefaultArrayAlloc: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Allocate uninitialized memory (no construction)
/** Standard method is alloc.  Also has support for calloc and aligned
 *  memory based on preprocessor defines.
 *  \param[in]  a_size  Number of elements of type T to allocate
 *  \param[in]  a_hint  Hint as to where to place memory.  Default is
 *                      nullptr (no hint).
 *//*-----------------------------------------------------------------*/

template <typename T, ArrayClassIndex I>
inline auto
DefaultArrayAlloc<T, I>::allocate(size_type                              a_size,
                                  typename void_allocator::const_pointer a_hint)
  -> pointer
{
#if defined USE_CALLOC
  void* addr = calloc(1, a_size*sizeof(T));
  if (addr == nullptr)
    {
      print_memory_line("Out of memory");
      MayDay::Error("Out of memory in DefaultArrayAlloc::allocate");
    }
#elif defined CH_ALIGN
  void* addr;
  int stat = CH_System::memalign(&addr, s_alignment, a_size*sizeof(T));
  if (stat != 0)
    {
      switch (stat)
        {
        case EINVAL:
          print_memory_line("Invalid alignment");
          MayDay::Error("Invalid alignment in DefaultArrayAlloc::allocate");
          break;
        case ENOMEM:
          print_memory_line("Out of memory");
          MayDay::Error("Out of memory in DefaultArrayAlloc::allocate");
          break;
        default:
          print_memory_line("Error");
          MayDay::Error("Error in DefaultArrayAlloc::allocate");
          break;
        }
    }
#else
  void* addr = malloc(a_size*sizeof(T));
  if (addr == nullptr)
    {
      print_memory_line("Out of memory");
      MayDay::Error("Out of memory in DefaultArrayAlloc::allocate");
    }
#endif
#ifdef CH_USE_MEMORY_TRACKING
  s_memtrackRegistry->addBytes(a_size*sizeof(T));
#endif
  return static_cast<pointer>(addr);
}

/*--------------------------------------------------------------------*/
//  Deallocate memory (no destruction)
/** \param[in]  a_addr  Address of memory
 *  \param[in]  a_size  Number of elements of type T to deallocate.
 *//*-----------------------------------------------------------------*/

template <typename T, ArrayClassIndex I>
inline void
DefaultArrayAlloc<T, I>::deallocate(pointer a_addr, size_type a_size)
{
  free(a_addr);
#ifdef CH_USE_MEMORY_TRACKING
  s_memtrackRegistry->subtractBytes(a_size*sizeof(T));
#endif
}


/*******************************************************************************
 *
 * Class DefaultArrayAlloc: static member definitions
 *
 ******************************************************************************/

//--Default alignment is 2*sizeof(void*) or 16 bytes on 64-bit machines

#ifdef CH_ALIGN
template <typename T, ArrayClassIndex I>
typename DefaultArrayAlloc<T, I>::size_type
DefaultArrayAlloc<T, I>::s_alignment = 2*sizeof(void*);
#endif

//--Initialization of memtrackRegistry (this memory is deallocated by the
//--dumpmemoryatexit() routine in memtrack.cpp after reporting)

#ifdef CH_USE_MEMORY_TRACKING
template <typename T, ArrayClassIndex I>
MemtrackRegistry* DefaultArrayAlloc<T, I>::s_memtrackRegistry =
  new MemtrackRegistry(
    // E.g. for BaseFab<int>, results in BaseFab<i> using gcc
    std::string(ArrayClassTr<I>::name())
    + "<"
    + std::string((typeid(T)).name())
    + ">");
#endif


/*******************************************************************************
 *
 * Class DefaultArrayAlloc: external related functions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
/// True equivalence
/** Means two allocators of same type can deallocate each others
 *  allocations
 *//*-----------------------------------------------------------------*/

template <typename T1, typename T2, ArrayClassIndex I>
inline constexpr bool
operator==(const DefaultArrayAlloc<T1, I>&,
           const DefaultArrayAlloc<T2, I>&) noexcept
{
  return true;
}

/*--------------------------------------------------------------------*/
/// False non-equivalence
/** Same as true equivalence
 *//*-----------------------------------------------------------------*/

template <typename T1, typename T2, ArrayClassIndex I>
inline constexpr bool
operator!=(const DefaultArrayAlloc<T1, I>&,
           const DefaultArrayAlloc<T2, I>&) noexcept
{
  return false;
}


/*******************************************************************************
 *
 * Class CUDAArrayAlloc: inline member definitions
 *
 ******************************************************************************/

#ifdef CH_GPU

/*--------------------------------------------------------------------*/
//  Allocate uninitialized memory (no construction)
/** Standard method is alloc.  Also has support for calloc and aligned
 *  memory based on preprocessor defines.
 *  \param[in]  a_size  Number of elements of type T to allocate
 *  \param[in]  a_hint  Hint as to where to place memory.  Default is
 *                      nullptr (no hint).
 *//*-----------------------------------------------------------------*/

template <typename T, ArrayClassIndex I>
HOSTDEVICE inline auto
CUDAArrayAlloc<T, I>::allocate(size_type                              a_size,
                               typename void_allocator::const_pointer a_hint)
  -> pointer
{
  const size_type numBytes = a_size*sizeof(T);
#ifdef __CUDACC__
  setDeviceAddress(reinterpret_cast<CH_Cuda::DevicePointer>(
                     malloc(numBytes)));
  return devicePtr();
#else
  const bool allocOn_cpu =
    reinterpret_cast<uintptr_t>(a_hint) & static_cast<uintptr_t>(AllocOn::cpu);
  const bool allocOn_gpu =
    reinterpret_cast<uintptr_t>(a_hint) & static_cast<uintptr_t>(AllocOn::gpu);
  m_devicePtr = 0;
  void *hostPtr = nullptr;
  if (allocOn_cpu)
    {
      if (allocOn_gpu || s_allocPageLocked)
        {
          // Allocate page-locked memory on host
          checkCudaErrors(cuMemHostAlloc(&hostPtr, numBytes, 0));
        }
      else
        {
          // Allocate regular memory on host
#if defined USE_CALLOC
          void* hostPtr = calloc(1, a_size*sizeof(T));
          if (hostPtr == nullptr)
            {
              print_memory_line("Out of memory");
              MayDay::Error("Out of memory in CUDAArrayAlloc::allocate");
            }
#elif defined CH_ALIGN
          int stat =
            CH_System::memalign(&hostPtr, s_alignment, a_size*sizeof(T));
          if (stat != 0)
            {
              switch (stat)
                {
                case EINVAL:
                  print_memory_line("Invalid alignment");
                  MayDay::Error("Invalid alignment in "
                                "CUDAArrayAlloc::allocate");
                  break;
                case ENOMEM:
                  print_memory_line("Out of memory");
                  MayDay::Error("Out of memory in CUDAArrayAlloc::allocate");
                  break;
                default:
                  print_memory_line("Error");
                  MayDay::Error("Error in CUDAArrayAlloc::allocate");
                  break;
                }
            }
#else
          hostPtr = malloc(a_size*sizeof(T));
          if (hostPtr == nullptr)
            {
              print_memory_line("Out of memory");
              MayDay::Error("Out of memory in DefaultArrayAlloc::allocate");
            }
#endif
        }
#ifdef CH_USE_MEMORY_TRACKING
      s_memtrackRegistryHost->addBytes(numBytes);
#endif
    }
  if (allocOn_gpu)
    {
      // Allocate on device
      checkCudaErrors(cuMemAlloc(&(m_devicePtr), numBytes));
#ifdef CH_USE_MEMORY_TRACKING
      s_memtrackRegistryDevice->addBytes(numBytes);
#endif
    }
  return static_cast<pointer>(hostPtr);
#endif
}

/*--------------------------------------------------------------------*/
//  Deallocate memory (no destruction)
/** \param[in]  a_addr  Address of processor memory.  If this
 *                      allocator was used on the host, it should
 *                      point to system memory (if memory was
 *                      allocated for the CPU, otherwise set to
 *                      nullptr).  If this allocated was used on the
 *                      device, it should point to device memory.
 *  \param[in]  a_size  Number of elements of type T to deallocate.
 *//*-----------------------------------------------------------------*/

template <typename T, ArrayClassIndex I>
HOSTDEVICE inline void
CUDAArrayAlloc<T, I>::deallocate(pointer a_addr, size_type a_size)
{
#ifdef __CUDACC__
  CH_assert(a_addr == reinterpret_cast<pointer>(m_devicePtr));
  free(a_addr);
  m_devicePtr = 0;
#else
  const bool allocOn_gpu = m_devicePtr != 0;
  if (a_addr != nullptr)  // Indicates allocation made on CPU
    {
      if (allocOn_gpu || s_allocPageLocked)
        {
          // Deallocate page-locked memory on host
          checkCudaErrors(cuMemFreeHost(a_addr));
        }
      else
        {
          // Dealocate regular memory on host
          free(a_addr);
        }
#ifdef CH_USE_MEMORY_TRACKING
      s_memtrackRegistryHost->subtractBytes(a_size*sizeof(T));
#endif
    }
  if (allocOn_gpu)
    {
      // Deallocate on device
      checkCudaErrors(cuMemFree(m_devicePtr));
      m_devicePtr = 0;
#ifdef CH_USE_MEMORY_TRACKING
      s_memtrackRegistryDevice->subtractBytes(a_size*sizeof(T));
#endif
    }
#endif
}

/*--------------------------------------------------------------------*/
//  Is the argument data configured for the GPU?
/** \param[in]  a_data  The pointer used by the base object making use
 *                      of this allocator
 *//*-----------------------------------------------------------------*/

template <typename T, ArrayClassIndex I>
inline constexpr bool
CUDAArrayAlloc<T, I>::isConfiguredFor_gpu(const pointer a_data) const noexcept
{
  return (a_data != nullptr &&
          m_devicePtr == reinterpret_cast<CH_Cuda::DevicePointer>(a_data));
}

/*--------------------------------------------------------------------*/
//  Synchronous copy from host to device
/** \param[in]  a_hostPtr
 *                      Pointer to host memory
 *  \param[in]  a_size  Number of elements of type T to copy.
 *//*-----------------------------------------------------------------*/

template <typename T, ArrayClassIndex I>
inline void
CUDAArrayAlloc<T, I>::copyToDevice(const_pointer   a_hostPtr,
                                   const size_type a_size) const noexcept
{
  CH_assert(a_hostPtr != nullptr && m_devicePtr != 0);
  checkCudaErrors(cuMemcpyHtoD(m_devicePtr, a_hostPtr, a_size*sizeof(T)));
}

/*--------------------------------------------------------------------*/
//  Asynchronous copy from host to device
/** \param[in]  a_hostPtr
 *                      Pointer to host memory
 *  \param[in]  a_size  Number of elements of type T to copy.
 *  \param[in]  a_stream
 *                      Stream identifier
 *//*-----------------------------------------------------------------*/

template <typename T, ArrayClassIndex I>
inline void
CUDAArrayAlloc<T, I>::copyToDeviceAsync(const_pointer   a_hostPtr,
                                        const size_type a_size,
                                        CUstream        a_stream) const noexcept
{
  CH_assert(a_hostPtr != nullptr && m_devicePtr != 0);
  checkCudaErrors(cuMemcpyHtoDAsync(m_devicePtr,
                                    a_hostPtr,
                                    a_size*sizeof(T),
                                    a_stream));
}

/*--------------------------------------------------------------------*/
//  Synchronous copy from device to host
/** \param[in]  a_hostPtr
 *                      Pointer to host memory
 *  \param[in]  a_size  Number of elements of type T to copy.
 *//*-----------------------------------------------------------------*/

template <typename T, ArrayClassIndex I>
inline void
CUDAArrayAlloc<T, I>::copyToHost(pointer         a_hostPtr,
                                 const size_type a_size) noexcept
{
  CH_assert(a_hostPtr != nullptr && m_devicePtr != 0);
  checkCudaErrors(cuMemcpyDtoH(a_hostPtr,
                               m_devicePtr,
                               a_size*sizeof(T)));
}

/*--------------------------------------------------------------------*/
//  Asynchronous copy from device to host
/** \param[in]  a_hostPtr
 *                      Pointer to host memory
 *  \param[in]  a_size  Number of elements of type T to copy.
 *  \param[in]  a_stream
 *                      Stream identifier
 *//*-----------------------------------------------------------------*/

template <typename T, ArrayClassIndex I>
inline void
CUDAArrayAlloc<T, I>::copyToHostAsync(pointer         a_hostPtr,
                                      const size_type a_size,
                                      CUstream        a_stream) noexcept
{
  CH_assert(a_hostPtr != nullptr && m_devicePtr != 0);
  checkCudaErrors(cuMemcpyDtoHAsync(a_hostPtr,
                                    m_devicePtr,
                                    a_size*sizeof(T),
                                    a_stream));
}

/*--------------------------------------------------------------------*/
//  Get the device pointer (as type pointer)
/** \return             Address of data on the GPU
 *//*-----------------------------------------------------------------*/

template <typename T, ArrayClassIndex I>
HOSTDEVICE inline auto
CUDAArrayAlloc<T, I>::devicePtr() const noexcept -> pointer
{
  return reinterpret_cast<pointer>(m_devicePtr);
}

/*--------------------------------------------------------------------*/
//  Reset data pointer on GPU (will not alter allocOn markers)
/** \param[in]  a_addr  New address of data on the GPU
 *//*-----------------------------------------------------------------*/

template <typename T, ArrayClassIndex I>
HOSTDEVICE inline void
CUDAArrayAlloc<T, I>::setDeviceAddress(
  const CH_Cuda::DevicePointer a_addr) noexcept
{
  m_devicePtr = a_addr;
}


/*******************************************************************************
 *
 * Class CUDAArrayAlloc: static member definitions
 *
 ******************************************************************************/

//--Default is to not allocate page-locked memory if not using GPU for these
//--arrays.

template <typename T, ArrayClassIndex I>
bool CUDAArrayAlloc<T, I>::s_allocPageLocked = false;

//--Initialization of memtrackRegistry (this memory is deallocated by the
//--dumpmemoryatexit() routine in memtrack.cpp after reporting)

#ifdef CH_USE_MEMORY_TRACKING
// For host
template <typename T, ArrayClassIndex I>
MemtrackRegistry* CUDAArrayAlloc<T, I>::s_memtrackRegistryHost =
  new MemtrackRegistry(
    // E.g. for BaseFab<int>, results in BaseFab<i> using gcc
    std::string(ArrayClassTr<I>::name())
    + "<"
    + std::string((typeid(T)).name())
    + ">");
// For device
template <typename T, ArrayClassIndex I>
MemtrackRegistry* CUDAArrayAlloc<T, I>::s_memtrackRegistryDevice =
  new MemtrackRegistry(
    // E.g. for BaseFab<int>, results in BaseFab_GPU<i> using gcc
    std::string(ArrayClassTr<I>::name())
    + "_GPU<"
    + std::string((typeid(T)).name())
    + ">",
    MemtrackRegistry::Device::GPU);
#endif


/*******************************************************************************
 *
 * Class CUDAArrayAlloc: external related functions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
/// True equivalence
/** Means two allocators of same type can deallocate each others
 *  allocations
 *//*-----------------------------------------------------------------*/

template <typename T1, typename T2, ArrayClassIndex I>
inline constexpr bool
operator==(const CUDAArrayAlloc<T1, I>&,
           const CUDAArrayAlloc<T2, I>&) noexcept
{
  return true;
}

/*--------------------------------------------------------------------*/
/// False non-equivalence
/** Same as true equivalence
 *//*-----------------------------------------------------------------*/

template <typename T1, typename T2, ArrayClassIndex I>
inline constexpr bool
operator!=(const CUDAArrayAlloc<T1, I>&,
           const CUDAArrayAlloc<T2, I>&) noexcept
{
  return false;
}

#endif  /* defined CH_GPU */


/*******************************************************************************
 *
 * Class Array_impl: inline member definitions
 *
 ******************************************************************************/

#ifdef CH_GPU
/*--------------------------------------------------------------------*/
/// Perform default initialization
/** This should only be necessary for arrays created as shared
 *  variables on the GPU, e.g., a declaration prefixed with
 *  __shared__.  In that case, there is no initialization.  This must
 *  be used before one of the define methods.
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
DEVICE inline void
Array_impl<T, Alloc>::INITIALIZEshared() noexcept
{
  m_data = nullptr;
  // Next also clears size and allocBy
  m_param = static_cast<size_type>(AllocOn::gpu);
}
#endif

/*--------------------------------------------------------------------*/
//  Constructor (takes union of AllocOn)
/** All parameters have defaults to ensure the constructor can be
 *  called without any parameters
 *  \param[in]  a_allocOn
 *                      Where to allocate memory.  For all allocators,
 *                      the default is to only allocate on the CPU.
 *                      Options:
 *                      AllocOn::cpu - on cpu
 *                      AllocOn::gpu - on gpu
 *                      Flags can be "or"ed together.  Alternatively,
 *                      AllocOn::all = AllocOn::cpu | AllocOn::gpu.
 *                      Note that the allocator itself must support
 *                      GPUs to allocate memory on GPUs.  The
 *                      DefaulArrayAlloc allocator will disregard this
 *                      parameter and always allocate on the CPU.
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
HOSTDEVICE inline
Array_impl<T, Alloc>::Array_impl(const unsigned a_allocOn) noexcept
  :
  m_data(nullptr),
  m_param(static_cast<size_type>(a_allocOn) &
          static_cast<size_type>(AllocOn::all))
{
  static_assert(std::is_unsigned<size_type>::value,
                "size_type must be unsigned");
  // Make sure requested allocations are supported (note: cpu is
  // always supported so only test gpu)
  if (allocOn_gpu())
    {
      CH_assert(Alloc::allocates_on_gpu::value);
    }
}

/*--------------------------------------------------------------------*/
//  Move construction is allowed
/** Move is not allowed for aliased memory since the new
 *  implementation could have improper scope.
 *  \param[in]  a_impl  Array to move
 *  \param[out] a_impl  Defunct state
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
HOSTDEVICE inline
Array_impl<T, Alloc>::Array_impl(Array_impl&& a_impl) noexcept
  :
  Alloc(std::move(a_impl)),          // Allocator is moved
  m_data(a_impl.m_data),             // and everything else is copied
  m_param(a_impl.m_param)
{
  static_assert(std::is_unsigned<size_type>::value,
                "size_type must be unsigned");
  a_impl.size(0, (unsigned)AllocOn::none, AllocBy::none);
  a_impl.m_data = nullptr;           // Just for safety to avoid use
}

/*--------------------------------------------------------------------*/
//  Move assignment is allowed
/** But the allocator must not insist on copying elements (whether or
 *  not it is actually moved).  Otherwise a compile-time error is
 *  forced.  Use an explicit copy operator if you need a workaround.
 *  Move is not allowed for aliased memory since the new
 *  implementation could have improper scope.
 *  \param[in]  a_impl  Array to move
 *  \param[out] a_impl  Defunct state
 *  \return             Assigned implementation
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
HOSTDEVICE inline auto
Array_impl<T, Alloc>::operator=(Array_impl&& a_impl) noexcept -> Array_impl&
{
  constexpr bool moveStorage =
    AllocTr::propagate_on_container_move_assignment::value ||
    Alloc::is_always_equal::value;  // Note: cannot use AllocTr until c++17
  // Either the allocators are equal or the allocator is propagated in which
  // case the move is constant time.  If we actually have to copy each element
  // element, well... we don't allow that.
  static_assert(moveStorage, "Chombo arrays do not support copying elements "
                "in in a move operation to satisfy allocator constraints");
  undefine();  // Clear our resources
  // Move the allocator (or not) using tag dispatch
  move_allocator(typename AllocTr::propagate_on_container_move_assignment(),
                 a_impl.get_allocator());
  // and everything else is copy assigned
  m_data = a_impl.m_data;
  m_param = a_impl.m_param;
  a_impl.size(0, (unsigned)AllocOn::none, AllocBy::none);
  a_impl.m_data = nullptr;           // Just for safety to avoid use
  return *this;
}

#ifdef CH_GPU
/*--------------------------------------------------------------------*/
//  Copy constructor that swaps pointers to data on the GPU
/** This is a shallow copy that only copies and configures the class
 *  object.  The data must be moved separately.  The alloc by is
 *  converted to an alias.
 *  \param[in]  a_impl  Array to configure for the GPU
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
inline
Array_impl<T, Alloc>::Array_impl(CH_Cuda::copy_to_device,
                                 const Array_impl& a_impl) noexcept
  :
  Alloc(Alloc::select_on_container_copy_construction(a_impl)),
  m_data(Alloc::devicePtr()),
  m_param(a_impl.m_param)
{
  // This will fail at runtime if memory was not allocated on the GPU.
  CH_assert(allocOn_gpu());
  setAllocBy(AllocBy::alias);
}

/*--------------------------------------------------------------------*/
//  Configure the array for use on the GPU
/** The data must be moved separately.  The pointer to the host is
 *  lost so this cannot be reverted for later use on the CPU.
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
inline void
Array_impl<T, Alloc>::configureForGPU() noexcept
{
  // This will fail at runtime if memory was not allocated on the GPU.
  CH_assert(allocOn_gpu());
  m_data = Alloc::devicePtr();
  setAllocBy(AllocBy::alias);
}
#endif

/*--------------------------------------------------------------------*/
//  Allocate and construct T
/** Extra arguments are forwarded to the constructor of T.  Tag
 *  dispatch is used to optimize for uninitialized memory if T is POD.
 *  Optimizations are also used of T is POD and there is a single
 *  argument to which each element should be initialized.
 *  \param[in]  a_size  Number of elements in the array
 *  \param[in]  a_args  Arguments, if any, for the constructor of T
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
template <typename... Args>
HOSTDEVICE inline void
Array_impl<T, Alloc>::define(const size_type a_size, Args&&... a_args)
{
  CH_assert(a_size <= c_max_size);
  undefine();
  size(a_size, AllocBy::array);
#ifdef CH_GPU
  // Pass allocOn as "hint" for allocation
  m_data = Alloc::allocate(a_size, reinterpret_cast<const_pointer>(m_param));
#else
  m_data = Alloc::allocate(a_size);
#endif
  // Invoke construction policy using tag dispatch
  if (m_data != nullptr)  // Can be nullptr for allocation only on GPU
    {
      constructArray(std::is_class<T>(), std::forward<Args>(a_args)...);
    }
}

/*--------------------------------------------------------------------*/
//  Alias constructed T
/** T already exists on a_addr, so do nothing
 *  \param[in]  a_addr  Address of memory to alias
 *  \param[in]  a_size  Number of elements in the array
 *  \note
 *  <ul>
 *    <li> Pointers stored in the allocator (if stateful) are not
 *         defined and must be set manually if needed.
 *  </ul>
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
HOSTDEVICE inline void
Array_impl<T, Alloc>::defineAlias(pointer a_addr, const size_type a_size)
{
  CH_assert(a_size <= c_max_size);
  undefine();
  size(a_size, AllocBy::alias);
  m_data = a_addr;
}

/*--------------------------------------------------------------------*/
//  Construct T on raw memory
/** Memory exists, and only construction is necessary.  Extra
 *  arguments are forwarded to the constructor of T.  Tag dispatch is
 *  used to optimize for uninitialized memory if T is POD.
 *  Optimizations are also used of T is POD and there is a single
 *  argument to which each element should be initialized.
 *  \param[in]  a_addr  Address of raw memory
 *  \param[in]  a_size  Number of elements in the array
 *  \param[in]  a_args  Arguments, if any, for the constructor of T
 *  \note
 *  <ul>
 *    <li> Pointers stored in the allocator (if stateful) are not
 *         defined and must be set manually if needed.
 *  </ul>
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
template <typename... Args>
HOSTDEVICE inline void
Array_impl<T, Alloc>::defineRaw(
  void*           a_addr,
  const size_type a_size,
  Args&&...       a_args)
{
  CH_assert(a_size <= c_max_size);
  undefine();
  size(a_size, AllocBy::raw);
  m_data = static_cast<pointer>(a_addr);
  // Invoke construction policy using tag dispatch
  constructArray(std::is_class<T>(), std::forward<Args>(a_args)...);
}

/*--------------------------------------------------------------------*/
//  Deallocate and destroy T
/** For any define
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
HOSTDEVICE inline void
Array_impl<T, Alloc>::undefine()
{
  switch (allocBy())
    {
    case AllocBy::array:
      // Invoke destruction policy using tag dispatch
      destroyArray(std::is_class<T>());
      Alloc::deallocate(m_data, size());
      break;
    case AllocBy::raw:
      // Invoke destruction policy using tag dispatch
      destroyArray(std::is_class<T>());
      break;
    case AllocBy::alias:
    case AllocBy::none:
      break;
    }
  m_data = nullptr;
  size(0, AllocBy::none);
}

/*--------------------------------------------------------------------*/
//  Class is defined (but unusable if number of elements is zero)
/** It is legal to define the class with 0 elements making it defined
 *  but unusable.
 *  \return             T - define was called
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
HOSTDEVICE inline bool
Array_impl<T, Alloc>::isDefined() const noexcept
{
  return allocBy() != AllocBy::none;
}

/*--------------------------------------------------------------------*/
//  Class is defined an usable
/** \return             T - define was called and number of elements
 *                          is > 0.
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
HOSTDEVICE inline bool
Array_impl<T, Alloc>::isUsable() const noexcept
{
  return isDefined() && size() > 0;
}

/*--------------------------------------------------------------------*/
//  Get allocator (unlike standard containers, this is modifiable)
/** STL containers return by value.  A modifiable reference is
 *  returned because, although it is dangerous, there are cases where
 *  it is useful to modify the state of the allocator (see
 *  ArrayOfMatrixAlloc for CHArray).  Also, this is the internal
 *  method used to get a handle to the allocator.
 *  \return             Allocator for the implementation
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
HOSTDEVICE inline Alloc&
Array_impl<T, Alloc>::get_allocator() noexcept
{
  return *static_cast<Alloc*>(this);
}

/*--------------------------------------------------------------------*/
//  Size (number of array elements)
/** \return             Number of elements
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
inline auto
HOSTDEVICE Array_impl<T, Alloc>::size() const noexcept
  -> size_type
{
  return m_param >> 4;
}

/*--------------------------------------------------------------------*/
//  Allocation type
/** \return             Source of memory for array
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
HOSTDEVICE inline auto
Array_impl<T, Alloc>::allocBy() const noexcept
  -> AllocBy
{
  return static_cast<AllocBy>(m_param & 0x3);
}

/*--------------------------------------------------------------------*/
//  Get all processors with allocation
/** \return             Bit flags of all processors with an allocation
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
HOSTDEVICE inline unsigned
Array_impl<T, Alloc>::allocOn() const noexcept
{
  return m_param & static_cast<size_type>(AllocOn::all);
}

/*--------------------------------------------------------------------*/
//  Is there an allocation on the CPU
/** \return             T - Have allocation on CPU
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
HOSTDEVICE inline bool
Array_impl<T, Alloc>::allocOn_cpu() const noexcept
{
  return m_param & static_cast<size_type>(AllocOn::cpu);
}

/*--------------------------------------------------------------------*/
//  Is there an allocation on the GPU
/** \return             T - Have allocation on GPU
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
HOSTDEVICE inline bool
Array_impl<T, Alloc>::allocOn_gpu() const noexcept
{
  return m_param & static_cast<size_type>(AllocOn::gpu);
}

/*--------------------------------------------------------------------*/
//  Index
/** \return             Const reference to element
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
HOSTDEVICE inline auto
Array_impl<T, Alloc>::operator[](const size_type a_idx) const noexcept
  -> const_reference
{
  return m_data[a_idx];
}

/*--------------------------------------------------------------------*/
//  Index
/** \return             Reference to element
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
HOSTDEVICE inline auto
Array_impl<T, Alloc>::operator[](const size_type a_idx) noexcept
 -> reference
{
  return m_data[a_idx];
}

/*--------------------------------------------------------------------*/
//  Access to the data pointer
/** \return             Const pointer to beginning of array
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
HOSTDEVICE inline auto
Array_impl<T, Alloc>::dataPtr() const noexcept
  -> const_pointer
{
  return m_data;
}

/*--------------------------------------------------------------------*/
//  Access to the data pointer
/** \return             Pointer to beginning of array
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
HOSTDEVICE inline auto
Array_impl<T, Alloc>::dataPtr() noexcept
  -> pointer
{
  return m_data;
}

/*--------------------------------------------------------------------*/
//  Access to the data pointer
/** \return             Const pointer to beginning of array
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
HOSTDEVICE inline auto
Array_impl<T, Alloc>::data() const noexcept
  -> const_pointer
{
  return m_data;
}

/*--------------------------------------------------------------------*/
//  Access to the data pointer
/** \return             Pointer to beginning of array
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
HOSTDEVICE inline auto
Array_impl<T, Alloc>::data() noexcept
  -> pointer
{
  return m_data;
}

/*--------------------------------------------------------------------*/
//  Begin iterator
/** \return             Pointer to beginning of array
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
HOSTDEVICE inline auto
Array_impl<T, Alloc>::begin() noexcept
  -> pointer
{
  return m_data;
}

/*--------------------------------------------------------------------*/
//  Constant begin iterator
/** \return             Const pointer to beginning of array
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
HOSTDEVICE inline auto
Array_impl<T, Alloc>::begin() const noexcept
  -> const_pointer
{
  return m_data;
}

/*--------------------------------------------------------------------*/
//  Constant begin iterator
/** \return             Const pointer to beginning of array
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
HOSTDEVICE inline auto
Array_impl<T, Alloc>::cbegin() const noexcept
  -> const_pointer
{
  return m_data;
}

/*--------------------------------------------------------------------*/
//  End iterator
/** \return             Pointer to 1 past end of array
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
HOSTDEVICE inline auto
Array_impl<T, Alloc>::end() noexcept
  -> pointer
{
  return m_data + size();
}

/*--------------------------------------------------------------------*/
//  Constant end iterator
/** \return             Const pointer to 1 past end of array
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
HOSTDEVICE inline auto
Array_impl<T, Alloc>::end() const noexcept
  -> const_pointer
{
  return m_data + size();
}

/*--------------------------------------------------------------------*/
//  Constant end iterator
/** \return             Const pointer to 1 past end of array
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
HOSTDEVICE inline auto
Array_impl<T, Alloc>::cend() const noexcept
  -> const_pointer
{
  return m_data + size();
}

#ifdef CH_GPU
/*--------------------------------------------------------------------*/
//  Synchronous copy from host to device
/** 
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
inline void
Array_impl<T, Alloc>::copyToDevice() const noexcept
{
  Alloc::copyToDevice(m_data, size());
}

/*--------------------------------------------------------------------*/
//  Asynchronous copy from host to device
/** \param[in]  a_stream
 *                      Stream identifier
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
inline void
Array_impl<T, Alloc>::copyToDeviceAsync(CUstream a_stream) const noexcept
{
  Alloc::copyToDeviceAsync(m_data, size(), a_stream);
}

/*--------------------------------------------------------------------*/
//  Synchronous copy from device to host
/** 
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
inline void
Array_impl<T, Alloc>::copyToHost() noexcept
{
  Alloc::copyToHost(m_data, size());
}

/*--------------------------------------------------------------------*/
//  Asynchronous copy from device to host
/** \param[in]  a_stream
 *                      Stream identifier
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
inline void
Array_impl<T, Alloc>::copyToHostAsync(CUstream a_stream) noexcept
{
  Alloc::copyToHostAsync(m_data, size(), a_stream);
}
#endif

/*--------------------------------------------------------------------*/
//  Set size, allocOn, and allocBy
/** \param[in]  a_size  Array size
 *  \param[in]  a_allocOn
 *                      Processors to allocate on
 *  \param[in]  a_allocBy
 *                      How the array was allocated
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
HOSTDEVICE inline void
Array_impl<T, Alloc>::size(const size_type a_size,
                           const unsigned  a_allocOn,
                           const AllocBy   a_allocBy) noexcept
{
  m_param = (a_size << 4) |
    static_cast<size_type>(a_allocOn) |
    static_cast<size_type>(a_allocBy);
}

/*--------------------------------------------------------------------*/
//  Set size and allocBy
/** \param[in]  a_size  Array size
 *  \param[in]  a_allocBy
 *                      How the array was allocated
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
HOSTDEVICE inline void
Array_impl<T, Alloc>::size(const size_type a_size,
                           const AllocBy   a_allocBy) noexcept
{
  m_param = (a_size << 4) |
    (m_param & static_cast<size_type>(AllocOn::all)) |
    static_cast<size_type>(a_allocBy);
}

/*--------------------------------------------------------------------*/
//  Set size
/** \param[in]  a_size  Array size
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
HOSTDEVICE inline void
Array_impl<T, Alloc>::size(const size_type a_size) noexcept
{
  m_param = (a_size << 4) | (m_param & 0xF);
}

/*--------------------------------------------------------------------*/
//  Set allocBy
/** \param[in]  a_allocBy
 *                      How the array was allocated
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
HOSTDEVICE inline void
Array_impl<T, Alloc>::setAllocBy(const AllocBy a_allocBy) noexcept
{
  m_param = static_cast<size_type>(a_allocBy) |
    (m_param & ~static_cast<size_type>(0x3));
}

/*--------------------------------------------------------------------*/
//  Move allocator
/** Method for moving allocator
 *  \param[in]  a_alloc Source allocator
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
HOSTDEVICE inline void
Array_impl<T, Alloc>::move_allocator(std::true_type, Alloc& a_alloc) noexcept
{
  get_allocator() = std::move(a_alloc);
}

/*--------------------------------------------------------------------*/
//  Unnecessary to move allocator (keep using current)
/** If it is not necessary to move the allocator, just keep current
 *  \param[in]  a_alloc Source allocator
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
HOSTDEVICE inline void
Array_impl<T, Alloc>::move_allocator(std::false_type, Alloc& a_alloc) noexcept
{ }

/*--------------------------------------------------------------------*/
//  Construct an array invoking constructor on elements
/** This is the general approach.  Any number of constructor
 *  arguments, including none, are taken and forwarded.
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
template <typename... Args>
HOSTDEVICE inline void
Array_impl<T, Alloc>::constructArray(std::true_type, Args&&... a_args)
{
  pointer p = m_data;
  Alloc& alloc = get_allocator();
  for (size_type n = size(); n--;)
    {
#ifdef __CUDACC__
      ::new((void *)p++) T(std::forward<Args>(a_args)...);
#else
      AllocTr::construct(alloc, p++, std::forward<Args>(a_args)...);
#endif
    }
}

/*--------------------------------------------------------------------*/
//  Construct an array and initilalize (not class type and at least 1
//  argument).
/** \param[in]  a_arg0  First argument to the constructor.
 *  \param[in]  a_args  Additional arguments to the constructor.
 *//*-----------------------------------------------------------------*/

/*
  This is not different from the above but is kept for legacy reasons if
  one wanted to do something different for a non-class type
 */

template <typename T, typename Alloc>
template <typename Arg0, typename... Args>
HOSTDEVICE inline void
Array_impl<T, Alloc>::constructArray(std::false_type,
                                     Arg0&&    a_arg0,
                                     Args&&... a_args)
{
  pointer p = m_data;
  Alloc& alloc = get_allocator();
  for (size_type n = size(); n--;)
    {
#ifdef __CUDACC__
      ::new((void *)p++) T(std::forward<Arg0>(a_arg0),
                           std::forward<Args>(a_args)...);
#else
      AllocTr::construct(alloc, p++,
                         std::forward<Arg0>(a_arg0),
                         std::forward<Args>(a_args)...);
#endif
    }
}

/*--------------------------------------------------------------------*/
//  Construct array as uninitialized (i.e., not class type and no
//  args).
/** This method has a specialization if T is Real in which case the
 *  array is initialized to c_BaseFabRealSetVal
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
HOSTDEVICE inline void
Array_impl<T, Alloc>::constructArray(std::false_type)
{
  pointer p = m_data;
  Alloc& alloc = get_allocator();
  for (size_type n = size(); n--;)
    {
#ifdef __CUDACC__
      ::new((void *)p++) T;
#else
      AllocTr::construct(alloc, p++);
#endif
    }
}

/*--------------------------------------------------------------------*/
//  Destroy an array invoking destructor on elements
/**            
 *//*-----------------------------------------------------------------*/

template <typename T, typename Alloc>
HOSTDEVICE inline void
Array_impl<T, Alloc>::destroyArray(std::true_type)
{
  pointer p = m_data;
  Alloc& alloc = get_allocator();
  for (size_type n = size(); n--;)
    {
#ifdef __CUDACC__
      p++->~T();
#else
      AllocTr::destroy(alloc, p++);
#endif
    }
}


/*******************************************************************************
 *
 * Class Array_impl: member specializations
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
//  Specialization for BaseFab of Real with no arguments using
//  DefaultArrayAlloc
/** Initialize to c_BaseFabRealSetVal
 *//*-----------------------------------------------------------------*/

template<>
inline void
Array_impl<Real, DefaultArrayAlloc<Real, ArrayClassIndex::BaseFab> >::
constructArray(std::false_type)
{
  Real* p = m_data;
  DefaultArrayAlloc<Real, ArrayClassIndex::BaseFab>& alloc = get_allocator();
  for (size_type n = size(); n--;)
    {
#ifdef __CUDACC__
      // To quiet the compiler, DefaultArrayAlloc cannot be used on GPU
      ::new((void *)p++) Real(c_BaseFabRealSetVal);
#else
      AllocTr::construct(alloc, p++, c_BaseFabRealSetVal);
#endif
    }
}

/*--------------------------------------------------------------------*/
//  Specialization for BaseFab of Real with no arguments using
//  CUDAArrayAlloc
/** Initialize to c_BaseFabRealSetVal
 *//*-----------------------------------------------------------------*/

#ifdef CH_GPU

template<>
inline void
Array_impl<Real, CUDAArrayAlloc<Real, ArrayClassIndex::BaseFab> >::
constructArray(std::false_type)
{
  Real* p = m_data;
  CUDAArrayAlloc<Real, ArrayClassIndex::BaseFab>& alloc = get_allocator();
  for (size_type n = size(); n--;)
    {
#ifdef __CUDACC__
      ::new((void *)p++) Real(c_BaseFabRealSetVal);
#else
      AllocTr::construct(alloc, p++, c_BaseFabRealSetVal);
#endif
    }
}

#endif  /* defined CH_GPU */
