#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _HASHIVSDATA_H_
#define _HASHIVSDATA_H_


/******************************************************************************/
/**
 * \file HashIVSData.H
 *
 * \brief Data for a HashIVS layout
 *
 *//*+*************************************************************************/

#include "Pool.H"
#include "ShapeArray.H"
#include "CH_HashTable.H"
#include "HashIVSLayout.H"

/*
  Notes:
    - tag __CUDACC__ means the file is processed with nvcc.  The class is
      modified whether on the cpu or gpu.
    - tag __CUDA_ARCH__ is defined for device code.  This only works in
      functions with the __device__ qualifier.  It cannot be used in the
      class declaration.
*/

#undef CUHOST
#undef HOSTDEVICE
#undef DEVICE
#ifdef __CUDACC__
#define CUHOST __host__
#define HOSTDEVICE __host__ __device__
#define DEVICE __device__
#else
#define CUHOST
#define HOSTDEVICE
#define DEVICE
#endif

#include "NamespaceHeader.H"

namespace HashIVS
{

// Be static if constructor takes layout (always compact)
// Be dynamic otherwise with option for compaction and moving to GPU

using static = std::integral_constant<unsigned, 0u>;
                                      ///< The layout must be separately
                                      ///< defined, synchronized, and cannot be
                                      ///< altered.  This option supports
                                      ///< parallel distribution of the layout
                                      ///< and data.  If the layout is packed,
                                      ///< this can be used on GPUs.
// using open = std::integral_constant<unsigned, 1u>;
                                      ///< Layout and data can be defined at the
                                      ///< same time.  However the hash-tables
                                      ///< are fixed in size and cannot be
                                      ///< resized.  Supported on CPUs and GPUs
                                      ///< but only in serial.  This
                                      ///< configuration is best for temporary
                                      ///< or intermediate data of a known size.
using dynamic = std::integral_constant<unsigned, 2u>;
                                      ///< Layout and data can be defined at the
                                      ///< same time.  The hash-tables are
                                      ///< dynamic and can be resized.  Only
                                      ///< supported on CPUs in serial.  This
                                      ///< configuration is best for temporary
                                      ///< or intermediate data of an unknown
                                      ///< size.
// using bare_static = std::integral_constant<unsigned, 1u + 4u>;
                                      ///< Only stores data.  The layout must be
                                      ///< externally defined if needed. The
                                      ///< hash-tables are dynamic and can be
                                      ///< resized.  Only supported on CPUs in
                                      ///< serial.  This configuration is best
                                      ///< for temporary or intermediate data of
                                      ///< a known size.
using bare_dynamic = std::integral_constant<unsigned, 2u + 4u>;
                                      ///< Only stores data.  The layout must be
                                      ///< externally defined if needed. The
                                      ///< hash-tables are dynamic and can be
                                      ///< resized.  Only supported on CPUs in
                                      ///< serial.  This configuration is best
                                      ///< for temporary or intermediate data of
                                      ///< an unknown size.

template <typename Config>
struct DataTraits;

template <>
struct DataTraits<closed>
{
  template <typename Key, typename T, typename Hash>
  using hashmap_type = CH_Hash::CompactHashTable<Key, T, Hash>;
};

template <>
struct DataTraits<open>
{
  template <typename Key, typename T, typename Hash>
  using hashmap_type = CH_Hash::OpenAddressHashTable<Key, T, Hash>;
};

template <>
struct DataTraits<dynamic>
{
  template <typename Key, typename T, typename Hash>
  using hashmap_type = std::unordered_map<Key, T, Hash>;
};

template <>
struct DataTraits<bare_open>
{
  template <typename Key, typename T, typename Hash>
  using hashmap_type = CH_Hash::OpenAddressHashTable<Key, T, Hash>;
};

template <>
struct DataTraits<bare_dynamic>
{
  template <typename Key, typename T, typename Hash>
  using hashmap_type = std::unordered_map<Key, T, Hash>;
};


template <typename MDDataSet>
struct MDDataSetInfo
{
  MDDataSet m_dataSet;
};


/*******************************************************************************
 */
///  Hash layout incorporating meta and data
/**
 *
 ******************************************************************************/

/* __META__
 * Is distributed with synchronization function
 * Bitset of 
 * 
 *
 * __DATA__
 * Data pointers must be separate
 * 
 */

// MD_HASHLOOP(HashIVSLayout, ivkey)
// {
//   MD_STAGE(A, ivkey, s_a);
//   MD_STAGE(B, ivkey, s_b);
//   MD_WORKLOOP(HashIVSLayout, ivkey, ivw)
// }

/*
  We do want expression template for this... with indexing.  Indexing always
  references the base (there will always be an offset applied to each index).
  Specify ranges when applying (having a range + and offset would be weird for
  components).
*/

// S=source, D=destination
// S receives okay to write to D
// S creates epoch using shared lock and writes data to D
// S also update contribution marker in same epoch
// S unlocks
// D can check contribution marker and, when up to date, gain *exclusive* lock
// on data to read it.
// When done, D sends new write authorization to S.

template <typename Config, typename T, typename UIx, UIx c_numBytes>
struct DataPool
{
  static void* get();
  static void remit(void *const ptr);
};

template <typename T, typename UIx, UIx c_numBytes>
struct DataPool<dynamic, T, UIx, c_numBytes>
{
  static void* get()
    { return s_pool.getPtr(); }
  static void remit(void *const a_ptr)
    { s_pool.returnPtr(a_ptr); }
  static Pool s_pool;                 ///< Pool for the data
};

template <typename T, typename UIx, UIx c_numBytes>
Pool DataPool<dynamic, T, UIx, c_numBytes>::s_pool(
  c_numBytes,
  "HashIVSDataPool",
  32,                         // Number of DataSets per chunk
  std::min(32ul, 4*sizeof(T)),  // Desired alignement (bytes)
  true);                      // Allow unaligned allocation

template <typename T, typename UIx, UIx c_numBytes>
struct DataPool<bare_dynamic, T, UIx, c_numBytes>
{
  static void* get()
    { return s_pool.getPtr(); }
  static void remit(void *const a_ptr)
    { s_pool.returnPtr(a_ptr); }
  static Pool s_pool;                 ///< Pool for the data
};

template <typename T, typename UIx, UIx c_numBytes>
Pool DataPool<bare_dynamic, T, UIx, c_numBytes>::s_pool(
  c_numBytes,
  "HashIVSDataPool",
  32,                         // Number of DataSets per chunk
  std::min(32ul, 4*sizeof(T)),  // Desired alignement (bytes)
  true);                      // Allow unaligned allocation

/*
  - Allow post dims for components, orientations, etc...
  - Need an iv (bracket) ability?
  - Linear, reduce remaining IV and mult by strides
  - Only has data dims ... workbox dims should be in assembly?

  ... you could hash more than in the IVS
  The hash is fast
  hash -> array -> array is nuts!
*/

template <typename T,
          typename Config,
          typename HashIVSLayout,
          typename HashIVSLayout::uix_type... Dims>  // DataSet dims
class Data
{
public:
  using value_type = T;
  using iix_type = typename HashIVSLayout::iix_type;
  using uix_type = typename HashIVSLayout::uix_type;
  static constexpr uix_type rank = sizeof...(Dims);
  using ivsvec_type = typename HashIVSLayout::iixvec_type;
  using iixvec_type = stc::Vector<iix_type, rank>;
  using uixvec_type = stc::Vector<uix_type, rank>;

  using layout_type = HashIVSLayout;

  static constexpr uixvec_type c_dims{ Dims... };
  static constexpr iixvec_type c_iix_dims = c_dims;
  static constexpr uixvec_type c_log2dims = stc::clog2(c_dims);
  static constexpr uix_type c_dataSz = c_dims.product();
  static constexpr uix_type c_dataSzB = c_dataSz*
    sizeof(typename std::aligned_storage<sizeof(T), alignof(T)>::type);

  using MDDataSet_type =
    shape::ArrayBracket<shape::CBareImpl<T,
                                         rank,
                                         shape::column_ordered,
                                         iix_type,
                                         Dims...>>;

  using MDDataSetInfo_type = MDDataSetInfo<MDDataSet_type>;

  using hashmap_type =
    typename DataTraits<Config>::template hashmap_type<
      iixvec_type,
      MDDataSetInfo_type,
      CH_Hash::google_CityHash<iixvec_type>>;

  using datapool_type = DataPool<Config, T, uix_type, c_dataSzB>;

  Data()
    :
    Data(layout_type(std::integral_constant<bool, !isBareLayout()>{}))
    {
      static_assert(isModifiableTable(),
                    "Default construction with a closed layout is not "
                    "permitted");
    }

  // template <typename... Dims>  // Sizes > hashed for range checking
  Data(const layout_type& a_layout)
    :
    m_layout(a_layout),
    m_cachedKey(std::numeric_limits<iix_type>::max()),
    m_cachedMapped(nullptr),
    m_tableOpen(isModifiableTable())
    {
      if (Config::value & 1u)  // Open addressing hash table
        {
          max_load_factor(0.7);
        }
    }

  // template <typename C = Config,
  //           std::enable_if_t<(C::value & 4u) == 0u, int> = 0>
  const layout_type& layout() const noexcept
    {
      static_assert(!isBareLayout(),
                    "Accessing a bare layout is not permitted");
      return m_layout;
    }

  // template <typename C = Config,
  //           std::enable_if_t<(C::value & 4u) == 0u, int> = 0>
  layout_type& layout() noexcept
    { 
      static_assert(!isBareLayout(),
                    "Accessing a bare layout is not permitted");
      return m_layout;
    }

  /// Whether empty or not
  bool hashmap_empty() const
    {
      return m_dataSets.empty();
    }

  /// Number of MDBitSets
  typename hashmap_type::size_type hashmap_size() const
    {
      return m_dataSets.size();
    }

  /// Number of buckets in the hash map
  typename hashmap_type::size_type bucket_count() const
    {
      return m_dataSets.bucket_count();
    }

  float load_factor() const
    {
      return m_dataSets.load_factor();
    }

  float max_load_factor() const
    {
      return m_dataSets.max_load_factor();
    }

  void max_load_factor(const float a_loadFactor)
    {
      m_dataSets.max_load_factor(a_loadFactor);
    }

  void reserve(const uix_type a_count)
    {
      m_dataSets.reserve(a_count);
    }

  template <typename... Idxs,
            std::enable_if_t<std::is_integral<
                               typename std::tuple_element<
                                 0, std::tuple<Idxs...>>::type>::value &&
                             (rank > layout_type::rank),
                             int> = 0>
  T& operator()(const ivsvec_type& a_ivs, Idxs... a_otherIdxs) noexcept
    {
      if (!isBareLayout() && isModifiableTable() && isTableOpen())
        {
          m_layout |= a_ivs;
        }
      iixvec_type iv(a_ivs, a_otherIdxs...);
      const iixvec_type key = stc::ccoarsen_log2(iv, c_log2dims);
      MDDataSetInfo_type& dataSetInfo = insertAndReturnDataSet(key);
      iv -= key*c_iix_dims;
      return dataSetInfo.m_dataSet(iv);
    }

  HOSTDEVICE MDDataSetInfo_type& insertAndReturnDataSet(
    const iixvec_type& a_key)
    {
      if (a_key != m_cachedKey)
        {
          m_cachedKey = a_key;
          if (isModifiableTable() && isTableOpen())
            {
              const auto ins =
                m_dataSets.insert({ a_key,
                                    MDDataSetInfo_type{ MDDataSet_type{} } });
              m_cachedMapped = &(ins.first->second);
              if (ins.second) // Insertion succeeded
                {
                  m_cachedMapped->m_dataSet.resetData(getNewDataSet());
                }
              // else            Already in place
            }
          else
            {
              const auto iter = m_dataSets.find(a_key);
              CH_assert(iter != m_dataSets.end());
              m_cachedMapped = &(iter->second);
            }
        }
      return *m_cachedMapped;
    }

  /// Erase
  /*
    if modifiable and open, remove cell
      if last iv in dataset, remove array
   */

  static constexpr bool isBareLayout()
    { return (Config::value & 4u); }

  static constexpr bool isModifiableTable()
    { return (Config::value & 3u); }
  
  /// Check if the hash table can be modified
  bool isTableOpen() const noexcept
    { return m_tableOpen; }

  /// Check if the hash table can not be modified
  bool isTableClosed() const noexcept
    { return !m_tableOpen; }

  void closeTable() noexcept
    { m_tableOpen = false; }

  void openTable() noexcept
    {
      static_assert(!std::is_same<Config, closed>::value,
                    "It is not permitted to open the hash table for a "
                    "HashIVS::Data type configured as closed");
      m_tableOpen = true;
    }


private:

  T* getNewDataSet() noexcept
    {
      return new(datapool_type::get()) T[c_dataSz];
    }

  void remitDataSet(T *const a_dataSet) noexcept
    {
      datapool_type::remit(a_dataSet);
    }

  layout_type m_layout;
  typename hashmap_type::key_type m_cachedKey;
  typename hashmap_type::mapped_type* m_cachedMapped;
  hashmap_type m_dataSets;
  bool m_tableOpen;
};

//--Definitions (concerns about ODR use)

// c_dims
template <typename T,
          typename Config,
          typename HashIVSLayout,
          typename HashIVSLayout::uix_type... Dims>
constexpr typename Data<T, Config, HashIVSLayout, Dims...>::uixvec_type
Data<T, Config, HashIVSLayout, Dims...>::c_dims;
// c_iix_dims
template <typename T,
          typename Config,
          typename HashIVSLayout,
          typename HashIVSLayout::uix_type... Dims>
constexpr typename Data<T, Config, HashIVSLayout, Dims...>::iixvec_type
Data<T, Config, HashIVSLayout, Dims...>::c_iix_dims;
// c_log2dims
template <typename T,
          typename Config,
          typename HashIVSLayout,
          typename HashIVSLayout::uix_type... Dims>
constexpr typename Data<T, Config, HashIVSLayout, Dims...>::uixvec_type
Data<T, Config, HashIVSLayout, Dims...>::c_log2dims;

}  // namespace HashIVS

#include "NamespaceFooter.H"

#endif  /* ! defined _HASHIVSDATA_H_ */
