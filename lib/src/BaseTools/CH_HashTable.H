#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _CH_HASHTABLE_H_
#define _CH_HASHTABLE_H_


/******************************************************************************/
/**
 * \file CH_HashTable.H
 *
 * \brief Hash tables constructed as a hash graph
 *
 * The HashGraphTable arranges all elements into buckets and stores these in
 * contiguous memory.  An indirection array, here denoted as the offset array,
 * is used to find the extents of a bucket.  Given a key, it is hashed to
 * find location idxBckt in the offset array.  The location in the element
 * array is given by m_offset[idxBckt] and the number of elements in the bucket
 * is given by m_offset[idxBckt+1] - m_offset[idxBckt].
 *
 * In general, consider to first use std::unordered_map to build a table.  Next
 * convert the unordered_map into a HashGraphTable for use on GPUs.  To move the
 * table to a GPU, two arrays, m_offset and m_elem, have to be copied.  Once
 * the HashGraphTable is constructued, no further insertions or deletions are
 * allowed.
 *
 * \note
 * <ul>
 *   <li> std::unordered_map uses open hashing (separate chaining), basically
 *        the opposite of open addressing.  The default max_load_factor is 1.0
 *        and good performance is obeserved in this range.  Resizing is
 *        automatic when load_factor becomes > max_load_factor.  For building
 *        a hash table of unknown size, this is the preferred choice.
 *   <li> In the HashGraphTable, using a load factor of 1 will lead to excessive
 *        empty elements in the offset array.  The default load factor is 2
 *        which leads to an offset array that is half the size of the number of
 *        elements.  This recommended value has a good balance between memory
 *        use and performance (on a GPU, there is essentially no impact to
 *        performance.
 *   <li> The hash tables in this file use the default hashing type of
 *        CH_Hash::XXHash<Key> instead of std::hash<Key>.   You can
 *        choose any hashing algorithm including the defaults when defining
 *        the table type.  CH_Hash::google_CityHash<Key> is also supported but
 *        XXHash appears to be slightly faster with similar quality.
 *   <li> The hash tables in this file use the default equal_to type of
 *        CH_Hash::equal_to<Key> instead of std::equal_to<Key>.  You can
 *        choose any hashing algorithm including the defaults when defining
 *        the table type.  For most purposes, these are identical and depend on
 *        operator==.
 *   <li> The tables extract the size of the hash from Hash::result_type.  This
 *        is depreciated in std::hash but we keep it to support 32-bit
 *        operations on the GPU.
 *   <li> Not all of the std::unordered_map interface is supported.  For lookup
 *        only find and insert are provided.  What is provided is identical
 *        in behavior to std::unordered_map (or fails an assertion).
 * </ul>
 *
 *//*+*************************************************************************/

#include <cmath>
#include <algorithm>
#ifndef __CUDACC__
#include <unordered_map>
#endif

#include "CH_Hash.H"
#include "ArrayAllocator.H"

/*
  Notes:
    - tag __CUDACC__ means the file is processed with nvcc.  The class is
      modified whether on the cpu or gpu.
    - tag __CUDA_ARCH__ is defined for device code.  This only works in
      functions with the __device__ qualifier.  It cannot be used in the
      class declaration.
*/

#undef DEVICE
#undef HOSTDEVICE
#ifdef __CUDACC__
#define DEVICE __device__
#define HOSTDEVICE __host__ __device__
#else
#define DEVICE
#define HOSTDEVICE
#endif

#include "BaseNamespaceHeader.H"

namespace CH_Hash
{

// Alias constructed data (no allocations or constructions)
using alias_format = std::integral_constant<unsigned, 1>;

//--Basic reimplementations for use on the GPU


/******************************************************************************/
/// pair
/** Provides basic pair
 *  \note
 *  <ul>
 *    <li> Nothing fancy including no piecewise_construct so try to keep the
 *         types of smaller size.
 *  </ul>
 *
 *//*+*************************************************************************/

template<typename T1, typename T2>
struct pair
{
  using first_type = T1;
  using second_type = T2;
  T1 first;
  T2 second;  
};

template<typename T1, typename T2>
HOSTDEVICE constexpr pair<std::decay_t<T1>, std::decay_t<T2>> make_pair(
  T1&& a_x, T2&& a_y)
{
  return pair<std::decay_t<T1>, std::decay_t<T2>>
    { std::forward<T1>(a_x), std::forward<T2>(a_y) };
}


/******************************************************************************/
/// equal_to
/** Provides basic equal_to
 *  \note
 *  <ul>
 *    <li> Does not derive from binary_function
 *  </ul>
 *
 *//*+*************************************************************************/

template<typename T>
struct equal_to
{
  HOSTDEVICE constexpr bool operator()(const T& a_x, const T& a_y) const
    { return a_x == a_y; }
};


/******************************************************************************/
/// Hash graph table using a load factor of 2
/** No insertions, deletions, or resizing is allowed.  A recommended strategy is
 *  to first build with a std::unordered_map and then use that to construct
 *  this class.
 *
 *//*+*************************************************************************/

template <typename Key,
          typename Mapped,
          typename Hash = CH_Hash::XXHash32<Key>,
          typename KeyEqual = CH_Hash::equal_to<Key>,
#ifdef CH_GPU
          typename Alloc = CUDAArrayAlloc<CH_Hash::pair<Key, Mapped>,
                                          ArrayClassIndex::HashGraphTable>
#else
          typename Alloc = DefaultArrayAlloc<CH_Hash::pair<Key, Mapped>,
                                             ArrayClassIndex::HashGraphTable>
#endif
          >
class HashGraphTable
{
  static constexpr float c_defaultMaxLoadFactor = 2.0f;
public:
  using key_type = Key;
  using mapped_type = Mapped;
  using value_type = CH_Hash::pair<key_type, mapped_type>;
  using size_type = unsigned;
  using index_type = typename std::make_signed<size_type>::type;
  using hasher = Hash;
  using key_equal = KeyEqual;
  using iterator = value_type*;
  using const_iterator = const value_type*;

  /// Default constructor
  HashGraphTable(
    const hasher&    a_hasher = hasher{},
    const key_equal& a_eq = key_equal{},
#ifdef CH_GPU
    const unsigned   a_allocOn =
    (unsigned)AllocOn::cpu | (unsigned)AllocOn::gpu
#else
    const unsigned   a_allocOn = (unsigned)AllocOn::cpu
#endif
    )
    :
    m_hasher(a_hasher),
    m_eq(a_eq),
    m_offset(a_allocOn),
    m_elem(a_allocOn)
    { }

  /// Construct from another container of value_type (probably another table)
  template <typename Container,
            std::enable_if_t<
              !std::is_same<std::decay_t<Container>,
                            HashGraphTable>::value, int> = 0>
  HashGraphTable(
    const Container& a_container,
    const key_type&  a_emptyKey,
    const float      a_maxLoadFactor = c_defaultMaxLoadFactor,
    const hasher&    a_hasher = hasher{},
    const key_equal& a_eq = key_equal{},
#ifdef CH_GPU
    const unsigned   a_allocOn =
    (unsigned)AllocOn::cpu | (unsigned)AllocOn::gpu
#else
    const unsigned   a_allocOn = (unsigned)AllocOn::cpu
#endif
    )
    :
    m_hasher(a_hasher),
    m_eq(a_eq),
    m_offset(a_allocOn),
    m_elem(a_allocOn)
    {
      define(a_container, a_emptyKey, a_maxLoadFactor);
    }

  /// Move construct from another container of value_types
  /** The elements of a_container should be considered degenerate as they
   *  have been moved to this object.  But the structure of a_container itself
   *  is unmodified and it should be cleared afterwards to make sure the
   *  the elements are not used.
   */
  template <typename Container,
            std::enable_if_t<
              !std::is_same<std::decay_t<Container>,
                            HashGraphTable>::value, int> = 0>
  HashGraphTable(
    Container&&      a_container,
    const key_type&  a_emptyKey,
    const float      a_maxLoadFactor = c_defaultMaxLoadFactor,
    const hasher&    a_hasher = hasher{},
    const key_equal& a_eq = key_equal{},
#ifdef CH_GPU
    const unsigned   a_allocOn =
    (unsigned)AllocOn::cpu | (unsigned)AllocOn::gpu
#else
    const unsigned   a_allocOn = (unsigned)AllocOn::cpu
#endif
    )
    :
    m_hasher(a_hasher),
    m_eq(a_eq),
    m_offset(a_allocOn),
    m_elem(a_allocOn)
    {
      define(std::move(a_container), a_emptyKey, a_maxLoadFactor);
    }

  /// Construct from another container of value_type (probably another table)
  /** This constructor does not require an empty key.  However, the cost is
   *  that another array of counts must be created during construction, so this
   *  will use more memory during construction.
   */
  template <typename Container,
            std::enable_if_t<
              !std::is_same<std::decay_t<Container>,
                            HashGraphTable>::value, int> = 0>
  HashGraphTable(
    const Container& a_container,
    const float      a_maxLoadFactor = c_defaultMaxLoadFactor,
    const hasher&    a_hasher = hasher{},
    const key_equal& a_eq = key_equal{},
#ifdef CH_GPU
    const unsigned   a_allocOn =
    (unsigned)AllocOn::cpu | (unsigned)AllocOn::gpu
#else
    const unsigned   a_allocOn = (unsigned)AllocOn::cpu
#endif
    )
    :
    m_hasher(a_hasher),
    m_eq(a_eq),
    m_offset(a_allocOn),
    m_elem(a_allocOn)
    {
      size_type numBucket;
      const size_type numElem =
        defineOffset(a_container, a_maxLoadFactor, numBucket);
      // Define a counter array for the number of elements in each bucket
      std::vector<size_type> cntElem(numBucket, (size_type)0);
      // Define and load the elements
      m_elem.define(numElem + 1, value_type{ key_type{}, mapped_type{} });
      for (const auto& value : a_container)
        {
          const unsigned idxBckt = bucket(value.first);
          const unsigned idxElem = m_offset[idxBckt] + cntElem[idxBckt];
          ++cntElem[idxBckt];
          m_elem[idxElem] = { value.first, value.second };
        }
    }

  /// Alias an existing HashGraphTable
  HashGraphTable(alias_format,
                 HashGraphTable& a_table)
    :
    m_hasher(a_table.m_hasher),
    m_eq(a_table.m_eq),
    m_offset(a_table.m_offset.allocOn()),
    m_elem(a_table.m_elem.allocOn())
    {
      m_offset.defineAlias(a_table.m_offset.data(), a_table.m_offset.size());
      m_elem.defineAlias(a_table.m_elem.data(), a_table.m_elem.size());
    }

  /// Weak construction from an existing container
  template <typename Container>
  void define(
    const Container& a_container,
    const key_type&  a_emptyKey,
    const float      a_maxLoadFactor = c_defaultMaxLoadFactor) noexcept
    {
      // Define the offset array
      size_type numBucket;
      const size_type numElem =
        defineOffset(a_container, a_maxLoadFactor, numBucket);
      // Define and load the elements
      m_elem.define(numElem + 1, value_type{ a_emptyKey, mapped_type{} });
      for (const auto& value : a_container)
        {
          const unsigned idxBckt = bucket(value.first);
          unsigned idxElem = m_offset[idxBckt];
          while (!m_eq(m_elem[idxElem].first, a_emptyKey)) ++idxElem;
          m_elem[idxElem].first  = value.first;
          m_elem[idxElem].second = value.second;
        }
    }

  /// Weak move construction from an existing container
  /** The elements of a_container should be considered degenerate as they
   *  have been moved to this object.  But the structure of a_container itself
   *  is unmodified and it should be cleared afterwards to make sure the
   *  the elements are not used.
   */
  template <typename Container>
  void define(
    Container&&     a_container,
    const key_type& a_emptyKey,
    const float     a_maxLoadFactor = c_defaultMaxLoadFactor) noexcept
    {
      // Define the offset array
      size_type numBucket;
      const size_type numElem =
        defineOffset(a_container, a_maxLoadFactor, numBucket);
      // Define and load the elements
      m_elem.define(numElem + 1, value_type{ a_emptyKey, mapped_type{} });
      for (const auto& value : a_container)
        {
          const unsigned idxBckt = bucket(value.first);
          unsigned idxElem = m_offset[idxBckt];
          while (!m_eq(m_elem[idxElem].first, a_emptyKey)) ++idxElem;
          m_elem[idxElem].first  = std::move(value.first);
          m_elem[idxElem].second = std::move(value.second);
        }
    }

private:

  /// Define the offset array
  /** hashed keys are used to index this array so its size (-1) is the number of
   *  buckets.
   */
  template <typename Container>
  size_type
  defineOffset(const Container& a_container,
               const float      a_maxLoadFactor,
               size_type&       a_numBucket) noexcept
    {
      const size_type numElem = a_container.size();
      const size_type numBucket =
        std::ceil((double)numElem/(double)a_maxLoadFactor);
      m_offset.define(numBucket + 1u, (size_type)0);
      // Count the number of elements in each bucket
      for (const auto& value : a_container)
        {
          ++m_offset[bucket(value.first)];
        }
      // Convert to a prefix sum array.  This is the offset array
      unsigned c = 0;
      for (unsigned i = 0; i != numBucket; ++i)
        {
          c += m_offset[i];
          m_offset[i] = c - m_offset[i];
        }
      CH_assert(c == numElem);
      m_offset[numBucket] = numElem;
      a_numBucket = numBucket;
      return numElem;
    }

public:

#ifdef CH_GPU
  /// Copy constructor that swaps pointers to data on the GPU
  /** This is a shallow copy that only copies the class object to the device
   *  (either through parameters or to constant memory).  The data (in m_elem)
   *  must be moved separately.
   */
  HashGraphTable(CH_Cuda::copy_to_device,
                 const HashGraphTable& a_table) noexcept
    :
    m_hasher(a_table.m_hasher),
    m_eq(a_table.m_eq),
    m_offset(AllocOn::gpu),  // This alias is on the GPU only
    m_elem(AllocOn::gpu)     // This alias is on the GPU only
    {
      // This will fail at runtime if memory was not allocated on the GPU.
      CH_assert(a_table.m_offset.allocOn_gpu());
      CH_assert(a_table.m_elem.allocOn_gpu());
      // The data pointer for the arrays point to device memory allocated for
      // them
      m_offset.defineAlias(a_table.m_offset.devicePtr(), a_table.m_offset.size());
      m_elem.defineAlias(a_table.m_elem.devicePtr(), a_table.m_elem.size());
    }

  /// Synchronous copy from host to device
  void copyToDevice() const noexcept
    {
      m_offset.copyToDevice();
      m_elem.copyToDevice();
    }

  /// Asynchronous copy from host to device
  /** \param[in]  a_stream
   *                      Stream identifier
   */
  void copyToDeviceAsync(
    CUstream a_stream = CH_Cuda::c_defaultStream) const noexcept
    {
      m_offset.copyToDeviceAsync(a_stream);
      m_elem.copyToDeviceAsync(a_stream);
    }

  /// Synchronous copy from device to host
  void copyToHost() noexcept
    {
      m_offset.copyToHost();
      m_elem.copyToHost();
    }

  /// Asynchronous copy from device to host
  /** \param[in]  a_stream
   *                      Stream identifier
   */
  void copyToHostAsync(
    CUstream a_stream = CH_Cuda::c_defaultStream) noexcept
    {
      m_offset.copyToHostAsync(a_stream);
      m_elem.copyToHostAsync(a_stream);
    }
#endif

//--Iterators

  /// Begin
  HOSTDEVICE iterator begin() noexcept
    {
      return m_elem.data();
    }

  /// Constant begin
  HOSTDEVICE const_iterator begin() const noexcept
    {
      return m_elem.data();
    }

  /// Constant begin
  HOSTDEVICE const_iterator cbegin() const noexcept
    {
      return m_elem.data();
    }

  /// End
  HOSTDEVICE iterator end() noexcept
    {
      return m_elem.data() + size();
    }

  /// Constant end
  HOSTDEVICE const_iterator end() const noexcept
    {
      return m_elem.data() + size();
    }

  /// Constant end
  HOSTDEVICE const_iterator cend() const noexcept
    {
      return m_elem.data() + size();
    }

//--Capacity

  /// Is the container empty
  HOSTDEVICE bool empty() const noexcept
    {
      return size() == (size_type)0;
    }

  /// Number of used elements
  HOSTDEVICE size_type size() const noexcept
    {
      return m_elem.size() - (size_type)1;
    }

//--Lookup

  /// Const find
  HOSTDEVICE const_iterator find(const key_type& a_key) const noexcept
    {
      const size_type idxBckt = bucket(a_key);
      for (int idxElem = m_offset[idxBckt], idxElem_end = m_offset[idxBckt + 1];
           idxElem != idxElem_end; ++idxElem)
        {
          if (m_eq(m_elem[idxElem].first, a_key)) return cbegin() + idxElem;
        }
      return cend();
    }

  /// Find
  HOSTDEVICE iterator find(const key_type& a_key) noexcept
    {
      const size_type idxBckt = bucket(a_key);
      for (int idxElem = m_offset[idxBckt], idxElem_end = m_offset[idxBckt + 1];
           idxElem != idxElem_end; ++idxElem)
        {
          if (m_eq(m_elem[idxElem].first, a_key)) return begin() + idxElem;
        }
      return end();
    }

#ifdef __CUDACC__

  /// Parallel const find
  template <unsigned NumThr, typename T>
  DEVICE const_iterator find(const key_type& a_key,
                             const unsigned  a_thrIdxQuery,
                             T*              a_index) const noexcept
    {
      static_assert(std::is_integral<T>::value, "T must be of integral type");
/* Note: For .target sm_6x or below, all threads in mask must execute the same __syncwarp() in convergence, and the union of all values in mask must be equal to the active mask. Otherwise, the behavior is undefined. */
      // Query the elements
      if (a_thrIdxQuery < NumThr)
        {
          // Prep the collaborative memory
          if (a_thrIdxQuery == 0)
            {
              *a_index = size();
            }
          const unsigned thrIdx =
            (threadIdx.z*blockDim.y + threadIdx.y)*blockDim.x + threadIdx.x;
          const unsigned shift = (thrIdx/NumThr)*NumThr;
          const unsigned mask =  ((1u << NumThr) - 1u) << shift;
          // printf("mask %u %u\n", mask, shift);
          __syncwarp(mask);
          size_type idxBckt = bucket(a_key);
          for (int idxElem = m_offset[idxBckt] + a_thrIdxQuery,
                 idxElem_end = m_offset[idxBckt + 1];
               idxElem < idxElem_end; idxElem += NumThr)
            {
              if (m_eq(m_elem[idxElem].first, a_key))
                {
                  *a_index = (T)idxElem;
                }
            }
        }
      __syncwarp(__activemask());
      return cbegin() + *a_index;
    }

  /// Parallel find
  template <unsigned NumThr, typename T>
  DEVICE iterator find(const key_type& a_key,
                       const unsigned  a_thrIdxQuery,
                       T*              a_index) noexcept
    {
      static_assert(std::is_integral<T>::value, "T must be of integral type");
/* Note: For .target sm_6x or below, all threads in mask must execute the same __syncwarp() in convergence, and the union of all values in mask must be equal to the active mask. Otherwise, the behavior is undefined. */
      // Query the elements
      if (a_thrIdxQuery < NumThr)
        {
          // Prep the collaborative memory
          if (a_thrIdxQuery == 0)
            {
              *a_index = size();
            }
          const unsigned thrIdx =
            (threadIdx.z*blockDim.y + threadIdx.y)*blockDim.x + threadIdx.x;
          const unsigned shift = (thrIdx/NumThr)*NumThr;
          const unsigned mask =  ((1u << NumThr) - 1u) << shift;
          // printf("mask %u %u\n", mask, shift);
          __syncwarp(mask);
          size_type idxBckt = bucket(a_key);
          for (int idxElem = m_offset[idxBckt] + a_thrIdxQuery,
                 idxElem_end = m_offset[idxBckt + 1];
               idxElem < idxElem_end; idxElem += NumThr)
            {
              if (m_eq(m_elem[idxElem].first, a_key))
                {
                  *a_index = (T)idxElem;
                }
            }
        }
      __syncwarp(__activemask());
      return begin() + *a_index;
    }

#endif

//--Bucket interface

  /// Number of buckets
  HOSTDEVICE size_type bucket_count() const noexcept
    {
      return m_offset.size() - 1;
    }

  /// Size of a bucket (all have the same size)
  HOSTDEVICE constexpr size_type bucket_size(
    const size_type a_idx) const noexcept
    {
      return m_offset[a_idx + 1] - m_offset[a_idx];
    }

  /// Hash to the offset array
  /**
   */
  HOSTDEVICE size_type bucket(const key_type& a_key) const noexcept
    {
      return m_hasher(a_key) % bucket_count();
    }

//--Observers

  /// Returns function used to hash the keys
  HOSTDEVICE hasher hash_function() const noexcept
    {
      return m_hasher;
    }

  /// Returns the function used to compare keys for equality
  HOSTDEVICE key_equal key_eq() const noexcept
    {
      return m_eq;
    }

//--Others

  void print_stats(std::ostream& a_out) const
    {
     size_type num0 = 0;
     size_type numV = 0;
     size_type maxC = 0;
     size_type num1 = 0;
     size_type num4p = 0;
     size_type num8p = 0;
     size_type totalE = 0;
     for (size_type i = 0, i_end = bucket_count(); i != i_end; ++i)
       {
         size_type c = bucket_size(i);
         if (c == 0) ++num0;
         if (c != 0) ++numV;
         maxC = std::max(maxC, c);
         if (c == 1) ++num1;
         if (c > 4) ++num4p;
         if (c > 8) ++num8p;
       }
     a_out << "Total elements           : " << size() << std::endl;
     a_out << "Total buckets            : " << bucket_count() << std::endl;
     a_out << "Populated buckets        : " << numV << std::endl;
     a_out << "Buckets with 0 elements  : " << num0 << std::endl;
     a_out << "Buckets with 1 element   : " << num1 << std::endl;
     a_out << "Buckets with >4 elements : " << num4p << std::endl;
     a_out << "Buckets with >8 elements : " << num8p << std::endl;
     a_out << "Maximum bucket size      : " << maxC << std::endl;
     a_out << "Average Pop. bucket size : " << (float)size()/(float)numV
           << std::endl;
    }

//--Data members

//protected:
public:

  hasher m_hasher;                    ///< Hashing algorithm
  key_equal m_eq;                     ///< Key comparison algorithm
  Array_impl<unsigned,
             typename Alloc::template rebind<unsigned>::other> m_offset;
  Array_impl<value_type, Alloc> m_elem;
                                      ///< Vector of elements.  There is extra
                                      ///< space for an end marker and to
                                      ///< support parallel queries across a
                                      ///< complete bucket.
};

//--Definitions (concerns about ODR use)

// c_defaultMaxLoadFactor
template <typename Key,
          typename Mapped,
          typename Hash,
          typename KeyEqual,
          typename Alloc>
constexpr float
HashGraphTable<Key, Mapped, Hash, KeyEqual, Alloc>::c_defaultMaxLoadFactor;

//--Preferred tables

#ifndef __CUDACC__
template <typename Key, typename Mapped>
using DynamicTable = std::unordered_map<Key,
                                        Mapped,
                                        CH_Hash::XXHash<Key>,
                                        CH_Hash::equal_to<Key>>;
#else
template <typename Key, typename Mapped>
using DynamicTable = void;
#endif

template <typename Key, typename Mapped>
using StaticTable = CH_Hash::HashGraphTable<Key,
                                            Mapped,
                                            CH_Hash::XXHash32<Key>,
                                            CH_Hash::equal_to<Key>>;

}  // namespace CH_Hash

#include "BaseNamespaceFooter.H"

#ifdef CH_GPU
#include "CudaDriver.H"

#include "BaseNamespaceHeader.H"

namespace CH_Cuda
{

/// Specialization of Converting for HashGraphTable
/** This alters the data pointer to point to GPU memory when a table is passed
 *  as a parameter to a GPU kernel.
 */
template <typename Key,
          typename Mapped,
          typename Hash,
          typename KeyEqual>
struct Converting<CH_Hash::HashGraphTable<
                    Key, Mapped, Hash, KeyEqual,
                    CUDAArrayAlloc<CH_Hash::pair<Key, Mapped>,
                                   ArrayClassIndex::HashGraphTable>>&>
{
  using type = CH_Hash::HashGraphTable<
    Key, Mapped, Hash, KeyEqual,
    CUDAArrayAlloc<CH_Hash::pair<Key, Mapped>, ArrayClassIndex::HashGraphTable>>;
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

#endif  /* _CH_HASHTABLE_H_ */
