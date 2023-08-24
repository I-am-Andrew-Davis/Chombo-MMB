#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _CH_CUDA_SHARED_PTR_H_
#define _CH_CUDA_SHARED_PTR_H_


/******************************************************************************/
/**
 * \file CH_Cuda_shared_ptr.H
 *
 * \brief Shared pointer supporting CUDA
 *
 * This class is designed for creating shared pointers on the CPU (perhaps as
 * member data of classes) and then using the pointers (or classes) on the GPU.
 * It is primarily designed for use with types T that have the copy constructor
 * T(copy_to_device, const T&)
 *
 *  \note
 *  <ul>
 *    <li> use_count() is not available in device code.  When the shared_ptr
 *         is copied to the GPU, it behaves more like a weak_ptr.
 *    <li> the pointer should only be accessed on the GPU.  It is not expected
 *         or supported to create new shared pointers there.
 *    <li> if not using the GPU, this class reverts to the behavior and code
 *         of std::shared_ptr
 *  </ul>
 *
 *//*+*************************************************************************/

#include <iostream>
#include <cstdio>
#include <cstdint>
#include <atomic>

#include "CH_Cuda.H"
#include "CH_Cuda_allocator.H"
#include "CH_assert.H"

/*
  Notes:
    - tag __CUDACC__ means the file is processed with nvcc.  The class can
      be modified whether on the cpu or gpu.
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

namespace CH_Cuda
{

/*******************************************************************************
 */
/// Allows for shared pointers that can be used on the GPU.
/**
 *  \tparam T           Element type given by the pointer
 *
 *  \note
 *  <ul>
 *    <li> Any details regarding thread safety are the same as for
 *         std::shared_ptr
 *    <li> A custom destructor is not supported.  CH_Cuda::Alloc<T>::deallocate
 *         will be used for deallocation of GPU memory and
 *         CH_Cuda::Alloc<T>::newPair should be used for allocation.  It will
 *         still work if you do otherwise, but the memory register will be
 *         inaccurate (tracking is enabled in CH_Cuda::Alloc).
 *    <li> use_count() is not available in device code.  When the shared_ptr
 *         is copied to the GPU, it behaves more like a weak_ptr.
 *    <li> the pointer should only be accessed on the GPU.  It is not expected
 *         or supported to create new shared pointers there.
 *    <li> if not using the GPU, this class reverts to the code of
 *         std::shared_ptr
 *  </ul>
 *//*+*************************************************************************/

template <typename T>
class shared_ptr : public std::shared_ptr<T>
{

public:

  using element_type    = T;
  using pointer         = T*;

/*--------------------------------------------------------------------*/
/// Default constructor
/** Leaves pointer in a null state
 *//*-----------------------------------------------------------------*/

  shared_ptr() noexcept
    :
    std::shared_ptr<T>()
#ifdef CH_GPU
    ,
    m_devicePtr(0),
    m_count(nullptr)
#endif
    { }

/*--------------------------------------------------------------------*/
/// Create a shared pointer to memory in a_ptrs
/**
 *//*-----------------------------------------------------------------*/

  shared_ptr(SymbolPair<T> a_ptrs)
    :
    std::shared_ptr<T>(a_ptrs.host,
                       []
                       (auto p)
                         {
                           CH_Cuda::Alloc<T>::deletePair(SymbolPair<T>(p));
                         })
#ifdef CH_GPU
    ,
    m_devicePtr(a_ptrs.device),
    m_count(nullptr)
#endif
    {
#ifdef CH_GPU
      if (a_ptrs.device)
        {
          m_count = new std::atomic<int>(1);
        }
#endif
    }

/*--------------------------------------------------------------------*/
/// Create a shared pointer to memory in a_hostPtr and a_devicePtr
/** The device pointer can be null
 *//*-----------------------------------------------------------------*/

  shared_ptr(pointer a_hostPtr, CH_Cuda::DevicePointer a_devicePtr = 0)
    :
    shared_ptr(SymbolPair<T>(a_hostPtr, a_devicePtr))
    { }

/*--------------------------------------------------------------------*/
//--Use synthesized move and move assignment
/*--------------------------------------------------------------------*/

  shared_ptr(shared_ptr&&) = default;
  shared_ptr& operator=(shared_ptr&&) = default;

/*--------------------------------------------------------------------*/
/// Create a copy of a shared pointer
/** Increments use_count by 1
 *//*-----------------------------------------------------------------*/

  shared_ptr(const shared_ptr& a_src) noexcept
    :
    std::shared_ptr<T>(a_src)
#ifdef CH_GPU
    ,
    m_devicePtr(a_src.m_devicePtr),
    m_count(a_src.m_count)
#endif
    {
#ifdef CH_GPU
      if (m_count)
        {
          m_count->fetch_add(1);
        }
#endif
    }

/*--------------------------------------------------------------------*/
/// Assignment of a shared pointer
/** Releases current pointer (if not nullptr).  Increments use_count
 *  of new pointer by 1.
 *//*-----------------------------------------------------------------*/

  shared_ptr& operator=(const shared_ptr& a_src) noexcept
    {
      if (&a_src != this)
        {
          std::shared_ptr<T>::operator=(a_src);
#ifdef CH_GPU
          localReset(0);
          m_devicePtr = a_src.m_devicePtr;
          m_count     = a_src.m_count;
          if (m_count)
            {
              m_count->fetch_add(1);
            }
#endif
        }
      return *this;
    }

/*--------------------------------------------------------------------*/
/// Destructor
/** Decrements use_count by 1 and deletes if last holder
 *//*-----------------------------------------------------------------*/

  ~shared_ptr() noexcept
    {
#ifdef CH_GPU
      localReset(0);
#endif
    }

#ifdef CH_GPU
/*--------------------------------------------------------------------*/
/// Copy constructor that swaps pointers to data on the GPU
/**
 *//*-----------------------------------------------------------------*/

  shared_ptr(copy_to_device, const shared_ptr& a_src) noexcept
    :
    std::shared_ptr<T>(),
    m_devicePtr(a_src.m_devicePtr),
    m_count(nullptr)  // Disables deallocations
    { }
#endif

/*--------------------------------------------------------------------*/
/// Reset to nullptr
/**
 *//*-----------------------------------------------------------------*/

  void
  reset() noexcept
    {
      std::shared_ptr<T>::reset();
#ifdef CH_GPU
      localReset(0);
#endif
    }

/*--------------------------------------------------------------------*/
/// Reset to memory in a_ptrs
/**
 *//*-----------------------------------------------------------------*/

  void
  reset(SymbolPair<T> a_ptrs) noexcept
    {
      std::shared_ptr<T>::reset(
        a_ptrs.host,
        []
        (auto p)
          {
            CH_Cuda::Alloc<T>::deletePair(SymbolPair<T>(p));
          });
#ifdef CH_GPU
      if (!a_ptrs.host)
        {
          CH_assert(a_ptrs.device == 0);
        }
      localReset(a_ptrs.device);
#endif
    }

/*--------------------------------------------------------------------*/
/// Reset to memory in a_hostPtr and a_devicePtr
/** The device pointer can be null
 *//*-----------------------------------------------------------------*/

  void
  reset(pointer a_hostPtr, CH_Cuda::DevicePointer a_devicePtr = 0)
    {
      std::shared_ptr<T>::reset(
        a_hostPtr,
        []
        (auto p)
          {
            CH_Cuda::Alloc<T>::deletePair(SymbolPair<T>(p));
          });
#ifdef CH_GPU
      if (!a_hostPtr)
        {
          CH_assert(a_devicePtr == 0);
        }
      localReset(a_devicePtr);
#endif
    }

#ifdef CH_GPU
/*--------------------------------------------------------------------*/
/// Get the pointer
/** The memory to access is selected here based on the compiler.  This
 *  code is designed for use with the driver API
 *//*-----------------------------------------------------------------*/

  HOSTDEVICE pointer
  get() const noexcept
    {
#ifdef __CUDACC__
      return reinterpret_cast<pointer>(m_devicePtr);
#else
      return std::shared_ptr<T>::get();
#endif      
    }
// #else use std::shared_ptr<T>::get()

/*--------------------------------------------------------------------*/
/// Get the device pointer (as a pointer)
/** 
 *//*-----------------------------------------------------------------*/

  HOSTDEVICE pointer
  devicePtr() const noexcept
    {
      return reinterpret_cast<pointer>(m_devicePtr);
    }

/*--------------------------------------------------------------------*/
/// Get the symbol pair
/**
 *//*-----------------------------------------------------------------*/

  SymbolPair<T>
  getPair() const noexcept
    {
      return SymbolPair<T>(std::shared_ptr<T>::get(), m_devicePtr);
    }

/*--------------------------------------------------------------------*/
/// Dereference the pointer
/**
 *//*-----------------------------------------------------------------*/

  HOSTDEVICE T&
  operator*() const noexcept
    {
      return *get();
    }
// #else use std::shared_ptr<T>::operator*()

/*--------------------------------------------------------------------*/
/// Get the pointer
/**
 *//*-----------------------------------------------------------------*/

  HOSTDEVICE pointer
  operator->() const noexcept
    {
      return get();
    }
// #else use std::shared_ptr<T>::operator->()

/*--------------------------------------------------------------------*/
/// Check if non-null
/**
 *//*-----------------------------------------------------------------*/

  HOSTDEVICE explicit
  operator bool() const noexcept
    {
      return get() != nullptr;
    }
// #else use std::shared_ptr<T>::operator bool()

/*--------------------------------------------------------------------*/
/// Synchronous copy of the shared object from host to device
/** This makes a copy of the object, converting it for using on the
 *  GPU, and then copying the bytes.  The object must support the
 *  copy constructor T(copy_to_device, const T&).
 *
 *  \note
 *  <ul>
 *    <li> This does not otherwise prepare the object for use on the
 *         GPU.  It may be necessary to additionally call
 *         get()->copyToDevice() or get()->copyToDeviceAsync()
 *         depending on the nature of T.
 *         The reason these are kept separate is because copying back
 *         to the host is likely more nuanced for most T.
 *    <li> The copy is necessarily synchronous so no asynchronous
 *         version is provided.
 *  </ul>
 *//*-----------------------------------------------------------------*/

  void
  copyToDevice() const noexcept
    {
      // Make a temporary object configured for the GPU
      element_type obj(CH_Cuda::copy_to_device{}, *get());
      SymbolPair<T> ptrs(&obj, m_devicePtr);
      ptrs.copyHtoD(1);
    }

/*--------------------------------------------------------------------*/
/// Synchronous copy of the shared object from device to host
/** WARNING: unless T is elementary (stores data itself) this routine
 *  is not necessary.  E.g., if T is a DynArray, there is no sense in
 *  copying the meta-data (the class structure) back to the host since
 *  it is unchanged from what is already on the host.  The data
 *  pointed to by DynArray would need to be copied back via
 *  get()->copyToHost().
 *  \note
 *  <ul>
 *    <li> This does not otherwise prepare the object for use on the
 *         host.  It may be necessary to additionally call
 *         get()->copyToHost() or get()->copyToHostAsync() depending
 *         on the nature of T.
 *    <li> The copy is not necessarily synchronous but to be
 *         consistent with copyToDevice, no ansynchronous version is
 *         provided.
 *  </ul>
 *//*-----------------------------------------------------------------*/

  void
  copyToHost() noexcept
    {
      getPair().copyDtoH(1);
    }

private:

/*--------------------------------------------------------------------*/

  void
  localReset(const CH_Cuda::DevicePointer a_devicePtr) noexcept
    {
      if (m_count)
        {
          if (m_count->fetch_sub(1) == 1)
            {
              delete m_count;
              CH_Cuda::Alloc<T>::deallocate(
                SymbolPair<T>(nullptr, m_devicePtr));
            }
        }
      m_devicePtr = a_devicePtr;
      m_count = nullptr;
      if (a_devicePtr)
        {
          m_count = new std::atomic<int>(1);
        }
    }
#endif  /* CH_GPU */

//--Data members

private:

#ifdef CH_GPU
  CH_Cuda::DevicePointer m_devicePtr; ///< Pointer on the device
  std::atomic<int>* m_count;          ///< Unfortunately needed if we want to
                                      ///< guarantee deallocation of
                                      ///< m_devicePtr in a multithreaded
                                      ///< environment.
                                      ///< std::shared_ptr::use_count() can
                                      ///< only be approximately used.
#endif
};

/*--------------------------------------------------------------------*/
/// Specialization of Converting for shared_ptr
/** This alters the data pointer to point to GPU memory when a pointer
 *  is passed as a parameter to a GPU kernel.
 *//*-----------------------------------------------------------------*/

#ifdef CH_GPU
template <typename T>
struct Converting<shared_ptr<T>&>
{
  using type = shared_ptr<T>;
  static type builder(type& a_arg)
    {
      type local(copy_to_device{}, a_arg);
      // std::cout << "DID conversion: " << local.begin() << std::endl;
      return local;
    }
};
#endif

}  // End of namespace CH_Cuda

#endif  /* ! defined _CH_CUDA_SHARED_PTR_H_ */
