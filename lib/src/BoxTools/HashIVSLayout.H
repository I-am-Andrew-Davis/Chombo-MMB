#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _HASHIVSLAYOUT_H_
#define _HASHIVSLAYOUT_H_


/******************************************************************************/
/**
 * \file HashIVSLayout.H
 *
 * \brief Meta data for a hashed IVS layout
 *
 *//*+*************************************************************************/

/*
 Priorities:
 numPts(), construction with box, construction with IVS/copy construction, makeEmptyBits(), intersection operator (&=), union operator (|=), contains(IntVect), contains(Box). And I also use IVSIterator quite often.
*/

#include <cstdint>
#include <unordered_map>
#include <memory>

#include "SPMD.H"
#include "SPACE.H"
#include "StcVector.H"
#include "MDBitSet.H"
#include "CH_Hash.H"
#include "CH_HashTable.H"
#include "Box.H"
#include "ShapeArray.H"
#include "ArrayAllocator.H"
#include "CH_Cuda_shared_ptr.H"
#include "CH_Cuda_allocator.H"

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

#include "NamespaceHeader.H"

/*
  To use CH_Hash, you have to indicate that the type you are hashing is hashable
  (i.e., you have considered the impact of padding) and provide the number of
  bytes to be hashed.  This needs to be done for the type IIx and rank used to
  build stc::Vector<IIx, rank>.  Note that IntVect = stc::Vector<int, SpaceDim>,
  but this class permits constructions at any rank.  The most common ranks
  (1--8) for 32-bit and 64-bit signed integers are declared hashable below.  You
  can still use any other stc::Vector<IIx, rank> but you will have to declare it
  as hashable separately.

  IntVect also enables at SpaceDim and SpaceDim+1 which is reason for if macros
  based on CH_SPACEDIM
*/
namespace CH_Hash
{
  // int32_t
#if CH_SPACEDIM>1
  template<>
  struct isHashable<stc::Vector<int32_t, 1>> : public std::true_type
  {
    static constexpr int c_hashSize = sizeof(stc::Vector<int32_t, 1>);
  };
#endif
#if CH_SPACEDIM>2
  template<>
  struct isHashable<stc::Vector<int32_t, 2>> : public std::true_type
  {
    static constexpr int c_hashSize = sizeof(stc::Vector<int32_t, 2>);
  };
#endif
#if CH_SPACEDIM<2 || CH_SPACEDIM>3
  template<>
  struct isHashable<stc::Vector<int32_t, 3>> : public std::true_type
  {
    static constexpr int c_hashSize = sizeof(stc::Vector<int32_t, 3>);
  };
#endif
#if CH_SPACEDIM<3 || CH_SPACEDIM>4
  template<>
  struct isHashable<stc::Vector<int32_t, 4>> : public std::true_type
  {
    static constexpr int c_hashSize = sizeof(stc::Vector<int32_t, 4>);
  };
#endif
#if CH_SPACEDIM<4 || CH_SPACEDIM>5
  template<>
  struct isHashable<stc::Vector<int32_t, 5>> : public std::true_type
  {
    static constexpr int c_hashSize = sizeof(stc::Vector<int32_t, 5>);
  };
#endif
#if CH_SPACEDIM<5 || CH_SPACEDIM>6
  template<>
  struct isHashable<stc::Vector<int32_t, 6>> : public std::true_type
  {
    static constexpr int c_hashSize = sizeof(stc::Vector<int32_t, 6>);
  };
#endif
#if CH_SPACEDIM!=6
  template<>
  struct isHashable<stc::Vector<int32_t, 7>> : public std::true_type
  {
    static constexpr int c_hashSize = sizeof(stc::Vector<int32_t, 7>);
  };
#endif
  template<>
  struct isHashable<stc::Vector<int32_t, 8>> : public std::true_type
  {
    static constexpr int c_hashSize = sizeof(stc::Vector<int32_t, 8>);
  };
  // int64_t
  template<>
  struct isHashable<stc::Vector<int64_t, 1>> : public std::true_type
  {
    static constexpr int c_hashSize = sizeof(stc::Vector<int64_t, 1>);
  };
  template<>
  struct isHashable<stc::Vector<int64_t, 2>> : public std::true_type
  {
    static constexpr int c_hashSize = sizeof(stc::Vector<int64_t, 2>);
  };
  template<>
  struct isHashable<stc::Vector<int64_t, 3>> : public std::true_type
  {
    static constexpr int c_hashSize = sizeof(stc::Vector<int64_t, 3>);
  };
  template<>
  struct isHashable<stc::Vector<int64_t, 4>> : public std::true_type
  {
    static constexpr int c_hashSize = sizeof(stc::Vector<int64_t, 4>);
  };
  template<>
  struct isHashable<stc::Vector<int64_t, 5>> : public std::true_type
  {
    static constexpr int c_hashSize = sizeof(stc::Vector<int64_t, 5>);
  };
  template<>
  struct isHashable<stc::Vector<int64_t, 6>> : public std::true_type
  {
    static constexpr int c_hashSize = sizeof(stc::Vector<int64_t, 6>);
  };
  template<>
  struct isHashable<stc::Vector<int64_t, 7>> : public std::true_type
  {
    static constexpr int c_hashSize = sizeof(stc::Vector<int64_t, 7>);
  };
  template<>
  struct isHashable<stc::Vector<int64_t, 8>> : public std::true_type
  {
    static constexpr int c_hashSize = sizeof(stc::Vector<int64_t, 8>);
  };
}  // namespace CH_Hash

/// Trimming for neighbor iterators
enum
{
  TrimCenter = (1<<0),
  TrimFace   = (1<<1),
  TrimEdge   = (1<<2),
  TrimCorner = (1<<3)
};

namespace HashIVS
{

using active_layout = std::true_type; ///< Constructs the layout as active.
                                      ///< This is the default state
using inactive_layout = std::false_type;
                                      ///< Constructs an inactive layout.  The
                                      ///< hash table is not allocated and no
                                      ///< members should be used

// #ifdef CH_GPU
// template <typename T>
// using StaticVector = Array_impl<T, CUDAArrayAlloc<T, ArrayClassIndex::HashIVS>>;
// #else
// template <typename T>
// using StaticVector =
//   Array_impl<T, DefaultArrayAlloc<T, ArrayClassIndex::HashIVS>>;
// #endif
// template <typename T>
// using DynamicVector = std::vector<T>;

#ifdef CH_MPI
template <typename T>
struct MPI_Type;

template <>
struct MPI_Type<int>
MPI_Type
{
  static constexpr MPI_Datatype iix_type = MPI_INT;
  static constexpr MPI_Datatype uix_type = MPI_UNSIGNED;
};
template <>
struct MPI_Type<unsigned>
MPI_Type
{
  static constexpr MPI_Datatype iix_type = MPI_INT;
  static constexpr MPI_Datatype uix_type = MPI_UNSIGNED;
};
template <>
struct MPI_Type<long long>
MPI_Type
{
  static constexpr MPI_Datatype iix_type = MPI_LONG_LONG;
  static constexpr MPI_Datatype uix_type = MPI_UNSIGNED_LONG_LONG;
};
template <>
struct MPI_Type<unsigned long long>
MPI_Type
{
  static constexpr MPI_Datatype iix_type = MPI_LONG_LONG;
  static constexpr MPI_Datatype uix_type = MPI_UNSIGNED_LONG_LONG;
};
template <>
struct MPI_Type<int16_t>
MPI_Type
{
  static constexpr MPI_Datatype iix_type = MPI_INT16_T;
  static constexpr MPI_Datatype uix_type = MPI_UINT16_T;
};
template <>
struct MPI_Type<uint16_t>
MPI_Type
{
  static constexpr MPI_Datatype iix_type = MPI_INT16_T;
  static constexpr MPI_Datatype uix_type = MPI_UINT16_T;
};
template <>
struct MPI_Type<int32_t>
MPI_Type
{
  static constexpr MPI_Datatype iix_type = MPI_INT32_T;
  static constexpr MPI_Datatype uix_type = MPI_UINT32_T;
};
template <>
struct MPI_Type<uint32_t>
MPI_Type
{
  static constexpr MPI_Datatype iix_type = MPI_INT32_T;
  static constexpr MPI_Datatype uix_type = MPI_UINT32_T;
};
template <>
struct MPI_Type<int64_t>
MPI_Type
{
  static constexpr MPI_Datatype iix_type = MPI_INT64_T;
  static constexpr MPI_Datatype uix_type = MPI_UINT64_T;
};
template <>
struct MPI_Type<uint64_t>
MPI_Type
{
  static constexpr MPI_Datatype iix_type = MPI_INT64_T;
  static constexpr MPI_Datatype uix_type = MPI_UINT64_T;
};
#endif

template<typename IIx, typename UIx>
struct AuxInfoDefault
{
  using iix_type = IIx;
  using uix_type = UIx;
  constexpr bool compareAuxEqual(const AuxInfoDefault& a_other) const noexcept
    { return true; }
  constexpr bool compareAuxLess(const AuxInfoDefault& a_other) const noexcept
    { return false; }
};

template<typename IIx, typename UIx>
struct AuxInfoBlock
{
  using iix_type = IIx;
  using uix_type = UIx;
  constexpr bool compareAuxEqual(const AuxInfoBlock& a_other) const noexcept
    { return m_blockIdx == a_other.m_blockIdx; }
  constexpr bool compareAuxLess(const AuxInfoBlock& a_other) const noexcept
    { return m_blockIdx < a_other.m_blockIdx; }
  iix_type m_blockIdx;
};

template <typename MDBitSet>
struct MDBitSetInfo
{
  using MDBitSet_type = MDBitSet;
  using iix_type = typename MDBitSet_type::iix_type;

  MDBitSetInfo()
    :
    m_idxBitSet(0),
    m_copyLevel(0)
#ifdef CH_MPI
    ,
    m_procID(-1),
    m_tag(-1)
#endif
    { }
  // template <typename... AuxArgs>
  // MDBitSetInfo(MDBitSet_type *const                   a_bitSet,
  //              const typename MDBitSet_type::iix_type a_procID,
  //              AuxArgs&&...                           a_auxArgs)
  //   :
  //   AuxInfo(std::forward<AuxArgs>(a_auxArgs)...),
  //   m_bitSet(a_bitSet),
  //   m_procID(a_procID)
  //   { }

  //**This should be a size_t giving a difference from a base address
  //**For new, refer to 0x0, for arrays refer to base of array depending on
  //**copyLevel
  std::uintptr_t m_idxBitSet;
  //  MDBitSet_type* m_bitSet;
//**FIXME how to do this without needing a hierarchy?
  int m_copyLevel;                    ///< Level of the current copy (used to
                                      ///< implement copy on write for copies
                                      ///< of the HashIVSLayout)
#ifdef CH_MPI
  std::forward_list<int>::iterator m_motionPuts;
                                      ///< Puts required of this bitset for
                                      ///< communicating data with neighbors.
                                      ///< This is only to be used on the host
                                      ///< but the data (the iterator itself)
                                      ///< can reside on the GPU
  int m_numMotionsPuts;               ///< Number of puts
  int m_procID;                       ///< Process to which this bitset is
                                      ///< assigned
  int m_tag;                          ///< The tag is used for MPI and is
                                      ///< unique on the process that owns the
                                      ///< bitset.  It is the index of the key
                                      ///< in m_globalInfo relative to the first
                                      ///< key on the process (given by m_global
                                      ///< Displ[procID]).
#endif
};

template <typename AuxInfo, typename AuxInfo::uix_type R>
struct MDBitSetGlobalInfo : public AuxInfo
{
  using iix_type = typename AuxInfo::iix_type;
  using uix_type = typename AuxInfo::uix_type;
  static constexpr uix_type rank = R;
  using iixvec_type = typename stc::Vector<iix_type, rank>;
  MDBitSetGlobalInfo(const iixvec_type& a_key, const float a_weight = 1.f)
    :
    m_key(a_key),
    m_weight(a_weight)
    { }
  static constexpr uix_type sizeBytes() noexcept
    {
      return sizeof(MDBitSetGlobalInfo);
    }
  constexpr bool compareMorton(const iixvec_type&        a_offset,
                               const MDBitSetGlobalInfo& a_other) const noexcept
    {
      return stc::morton(m_key + a_offset, a_other.m_key + a_offset);
    }
  iixvec_type m_key;
  float m_weight;
};

template <typename T>
struct MortonBalance
{
  using MDBitSetGlobalInfo_type = T;
  using uix_type = typename MDBitSetGlobalInfo_type::uix_type;
  using iixvec_type = typename MDBitSetGlobalInfo_type::iixvec_type;
  /**
   *  \param[in]  a_offset
   *                      The minimum negative key (or zero).  Because
   *                      of twos complement, if this offset is not
   *                      applied, there would be a discontinuity at 0.
   */
  MortonBalance(const MDBitSetGlobalInfo_type* a_globalInfo,
                const iixvec_type&             a_offset)
    :
    m_globalInfo(a_globalInfo),
    m_offset(a_offset)
    { }
  bool operator()(const uix_type& a_x, const uix_type& a_y) const noexcept
    {
      const MDBitSetGlobalInfo_type& x = m_globalInfo[a_x];
      const MDBitSetGlobalInfo_type& y = m_globalInfo[a_y];
      if (x.compareAuxEqual(y))
        {
          return x.compareMorton(m_offset, y);
        }
      else
        {
          return x.compareAuxLess(y);
        }
    }
  const MDBitSetGlobalInfo_type* m_globalInfo;
                                      ///< A pointer array of MDBitSetGlobalInfo
                                      ///< to be sorted using Morton ordering
  const iixvec_type m_offset;         ///< The minimum negative key (or zero).
                                      ///< Because of twos complement, if this
                                      ///< offset is not applied, there would be
                                      ///< a discontinuity at 0.
};


/*******************************************************************************
 */
///  Hashed IntVectSet layout
/**
 *   Sync - find neighbors
 *   Grow?
 *
 ******************************************************************************/

/*
  Template on these boxes?
  DataBox <= WorkBox <= IVSBox
  E.g 4^3    16^3       16^3

  __META__
  Is distributed with synchronization function
  Bitset of

  __DATA__
  Data pointers must be separate

  Sync
    Share all meta in vector
 */

// template <typename IIxVec>
// struct LocalProcMap
// {
//   int operator()(const IIxVec& a_key) const noexcept
//     {
//       return procID();
//     }
// };

template <typename Table>
struct HashTableTraits;

template <typename Key, typename Mapped>
struct HashTableTraits<CH_Hash::StaticTable<Key, Mapped>>
{
  using is_static = std::true_type;
  using is_dynamic = std::false_type;
  // This next type is an array if the table is static
#ifdef CH_GPU
  using GlobalInfoCont_type =
    Array_impl<Mapped, CUDAArrayAlloc<Mapped, ArrayClassIndex::HashIVS>>;
#else
  using GlobalInfoCont_type =
    Array_impl<Mapped, DefaultArrayAlloc<Mapped, ArrayClassIndex::HashIVS>>;
#endif
};

template <typename Key, typename Mapped>
struct HashTableTraits<CH_Hash::DynamicTable<Key, Mapped>>
{
  using is_static = std::false_type;
  using is_dynamic = std::true_type;
  // This next type is also an unordred_map if the table is dynamic
  using GlobalInfoCont_type = CH_Hash::DynamicTable<Key, Mapped>;
};


template <typename W, typename IIx, typename UIx,
          template <typename, typename> class HashTable,
          template <typename, typename> class AuxInfo,
          UIx D0b, UIx... Dimsb>
class Layout
{
public:
  // For external use
  using word_type = W;
  // A signed integer type with good performance on all devices (32-bit
  // recommended)
  using iix_type = IIx;
  // An unsigned integer type with good performance on all devices (32-bit
  // recommended)
  using uix_type = UIx;
  // An unsigned integer type typically representing storage sizes (set as
  // std::size_t and 64-bit is recommended)
  using size_type = std::size_t;
  // The MPI type for diplacements in the MPI implementation.  Usually this is
  // int.  Do not change unless you are sure your MPI implementation supports
  // this.  The maximum value dictates the total number of global bitsets.
  using displ_type = int;
  
  static constexpr uix_type rank = 1 + sizeof...(Dimsb);

  // Mostly for internal use but made public
  using MDBitSet_type = MDBitSet<word_type,
                                 iix_type,
                                 uix_type,
                                 sizeof(word_type),
                                 D0b, Dimsb...>;
  using AuxInfo_type = AuxInfo<iix_type, uix_type>;
  using MDBitSetInfo_type = MDBitSetInfo<MDBitSet_type>;
  using MDBitSetGlobalInfo_type = MDBitSetGlobalInfo<AuxInfo_type, rank>;

  using uixvec_type = typename stc::Vector<uix_type, rank>;
  using iixvec_type = typename stc::Vector<iix_type, rank>;
  static constexpr uixvec_type c_dimsb{ D0b, Dimsb... };
  static constexpr uixvec_type c_log2dimsb = stc::clog2(c_dimsb);
  static constexpr iixvec_type c_iix_keyStep =
    stc::ccoarsen_log2(MDBitSet_type::c_iix_dimsb, c_log2dimsb);
  static constexpr iixvec_type c_emptyKey =
    iixvec_type(std::numeric_limits<iix_type>::max());

  using hashmap_type = HashTable<iixvec_type, MDBitSetInfo_type>;
  // If the hash table is dynamic, this type is also a dynamic hash table.  If
  // the hash table is static, this is an array of MDBitSetGlobalInfo_type
  using GlobalInfoCont_type =
    typename HashTableTraits<HashTable<iixvec_type, MDBitSetGlobalInfo_type>>
    ::GlobalInfoCont_type;

#ifdef CH_GPU
  template <typename T>
  using alloc_type = CUDAArrayAlloc<T, ArrayClassIndex::HashIVS>;
#else
  template <typename T>
  using alloc_type = DefaultArrayAlloc<T, ArrayClassIndex::HashIVS>;
#endif

  Layout()
    {
      reset();
    }

  Layout(active_layout)
    {
    }

  Layout(inactive_layout)
    :
    m_cachedKey(std::numeric_limits<iix_type>::max()),
    m_cachedMapped(nullptr),
    m_numLocalBS(0),
    m_numGlobalBS(0)
    {
    }

  /// Copy constructor
  Layout(const Layout& a_layout)
    :
    m_cachedKey(std::numeric_limits<iix_type>::max()),
    m_cachedMapped(nullptr),
    m_numLocalBS(0),
    m_numGlobalBS(0)
    {
    }

  // /// Number of MDBitSets
  // typename HashMap::size_type size() const
  //   {
  //     return m_bitSets->size();
  //   }

  /// Whether empty or not
  bool hashmap_empty() const
    {
      return m_bitSets->empty();
    }

  /// Number of MDBitSets
  typename hashmap_type::size_type hashmap_size() const
    {
      return m_bitSets->size();
    }

  /// Number of buckets in the hash map
  typename hashmap_type::size_type bucket_count() const
    {
      return m_bitSets->bucket_count();
    }

  // float load_factor() const
  //   {
  //     return m_bitSets->load_factor();
  //   }

  // float max_load_factor() const
  //   {
  //     return m_bitSets->max_load_factor();
  //   }

  // void max_load_factor(const float a_loadFactor)
  //   {
  //     m_bitSets->max_load_factor(a_loadFactor);
  //   }

  void reserve(const UIx a_count)
    {
      if (HashTableTraits<hashmap_type>::is_dynamic::value)
        {
          m_bitSets->reserve(a_count);
          m_globalInfo->reserve(a_count);
        }
    }

  void reserve(const Box&  a_domain,
               const int   a_numProc = numProc(),
               const float a_mult = 1.0)
    {
      const iixvec_type ivLo = a_domain.smallEnd();
      const iixvec_type ivHi = a_domain.bigEnd();
      reserve(ivLo, ivHi, a_numProc, a_mult);
    }

  void reserve(const iixvec_type& a_ivLo,
               const iixvec_type& a_ivHi,
               const int     a_numProc = numProc(),
               const float   a_mult = 1.0)
    {
      CH_assert(a_ivHi >= a_ivLo);
      CH_assert(a_numProc >= 0);
      CH_assert(a_mult > 0.f);
      const iixvec_type keyLo = stc::ccoarsen_log2(a_ivLo, c_log2dimsb);
      const iixvec_type keyHi = stc::ccoarsen_log2(a_ivHi, c_log2dimsb);
      const uix_type count =
        (std::min(1, static_cast<iix_type>(
                    stc::product(keyHi - keyLo + 1)*a_mult)) +
         a_numProc - 1)/a_numProc;
      reserve(count);
    }

/*--------------------------------------------------------------------*/
/**  \name Union operators
 *//*-----------------------------------------------------------------*/
//@{

  /// Unions a single IntVect into this IntVectSet
#ifndef __CUDACC__
  Layout& operator|=(iixvec_type a_iv)
    {
      const iixvec_type key = stc::ccoarsen_log2(a_iv, c_log2dimsb);
      MDBitSetInfo_type& bsInfo = insertAndReturnBitSet(key);
      a_iv -= key*MDBitSet_type::c_iix_dimsb;  // Zero-based in BS
      getBitSet(bsInfo).setb(a_iv);
      return *this;
    }

  /// Unions the IntVects in the Box into this IntVectSet
  Layout& operator|=(const Box& a_box)
    {
      forEachWordInBox(a_box,
                       []
                       (const iixvec_type&                      a_ivw,
                        const typename MDBitSet_type::word_type a_mask,
                        typename MDBitSet_type::word_type&      a_word)
                         {
                           a_word |= a_mask;
                         });
      return *this;
    }
#endif

  /// Unions the IntVects in the Box \b into this IntVectSet
  // HashIntVectSet& operator|=(const Box& a_box);
//@}

/*--------------------------------------------------------------------*/
/**  \name Local queries
 *//*-----------------------------------------------------------------*/
//@{
  /// Return the number of local IntVects in the layout
  HOSTDEVICE size_type
  numPts() const noexcept
    {
      size_type cnt = 0;
      for (const auto& value : m_bitSets)
        {
          const MDBitSetInfo_type& bsInfo = value.second;
          if (bsInfo.m_procID == procID())
            {
              cnt += getBitSet(bsInfo).count();
            }
        }
      return cnt;
    }
//@}

  /// Get the global info struct for a bitset.
  /** The argument a_iv must belong to an already existing bitset but it does
   *  not need to be set.
   */
  template <typename H = hashmap_type,
            std::enable_if_t<HashTableTraits<H>::is_dynamic::value, int> = 0>
  MDBitSetGlobalInfo_type&
  globalInfo(const iixvec_type& a_iv) noexcept
    {
      const iixvec_type key = stc::ccoarsen_log2(a_iv, c_log2dimsb);
      const auto iter = m_globalInfo->find(key);
      CH_assert(iter != m_globalInfo->end());
      return iter->second;
    }

  // f(iv, mask, word)
  /** If a_box has less dimensions than the rank, extra are set to zero
   *  If a_box has more dimensions than rank, only rank dimensions are used
   */
  template <typename F>
  void
  forEachWordInBox(const Box& a_box, F&& f)
    {
      const iixvec_type ivLo = a_box.smallEnd();
      const iixvec_type ivHi = a_box.bigEnd();
      const iixvec_type keyLo = stc::ccoarsen_log2(ivLo, c_log2dimsb);
      const iixvec_type keyHi = stc::ccoarsen_log2(ivHi, c_log2dimsb);
      stc::nestedLoop(
        keyLo,
        keyHi,
        c_iix_keyStep,
        [this, &ivLo, &ivHi, &f]
        (const iixvec_type& a_key)
          {
            const iixvec_type iv = a_key*MDBitSet_type::c_iix_dimsb;
            const iixvec_type loRng = stc::max(iv, ivLo);
            const iixvec_type hiRng = stc::min(
              iv + MDBitSet_type::c_iix_dimsb - iixvec_type(1), ivHi);
            MDBitSetInfo_type& bsInfo = insertAndReturnBitSet(a_key);
            getBitSet(bsInfo).forEachWord(
              iv,
              loRng - iv,
              hiRng - iv,
              std::forward<F>(f));
          });
    }

  HOSTDEVICE MDBitSetInfo_type&
  insertAndReturnBitSet(const iixvec_type& a_key)
    {
      if (a_key != m_cachedKey)
        {
          m_cachedKey = a_key;
          const auto ins =
            m_bitSets->insert({ a_key, MDBitSetInfo_type{} });
          m_cachedMapped = &(ins.first->second);
          if (ins.second) // Insertion succeeded
            {
              m_cachedMapped->m_idxBitSet = getNewBitSet();
#ifdef CH_MPI
              m_cachedMapped->m_procID = procID();
#endif
              const auto ins2 =
                m_globalInfo->insert({ a_key, MDBitSetGlobalInfo_type(a_key) });
              CH_assert(ins2.second);  // Must succeed if insertion succeeded
                                       // for m_bitSets.
            }
          // else            Already in place
        }
      return *m_cachedMapped;
    }

  HOSTDEVICE const MDBitSetInfo_type& insertAndReturnBitSet(
    const iixvec_type& a_key) const
    {
      if (a_key != m_cachedKey)
        {
          m_cachedKey = a_key;
          const auto iter = m_bitSets->find(a_key);
          CH_assert(iter != m_bitSets->end());  // Must exist, insertion not
                                                // allowed
          m_cachedMapped = &(iter->second);
        }
      return *m_cachedMapped;
    }

  // f(iv)
  template <typename F>
  void forEachIV(F&& f)
    {
      for (const auto& value : *m_bitSets)
        {
          const MDBitSetInfo_type& bsInfo = value.second;
#ifdef CH_MPI
          if (bsInfo.m_procID == procID())
#endif
            {
              getBitSet(bsInfo).forEachTrueBit(
                value.first*MDBitSet_type::c_iix_dimsb, std::forward<F>(f));
            }
        }
    }

  /// Synchronizing constructor
  /*  Could be a resync

      1) Number of BS
      2) Global info
      3) Rebalance BS

      2a) Communicate all local BS
      3) Find candidate neighbors and store in hashmap
      2b) Finish communication
      4a) Optional (rebalance)
      5) Loop over global BS and match with neighbors
      4b) Finish communication
      6) Send neighbor information (move for rebalance?)
  */

#if 0
  template <template <typename, typename> class SrcHashTable>
  void sync(
    Layout<
      W, IIx, UIx, SrcHashTable, AuxInfo, D0b, Dimsb...
    >&             a_layout,
    const iix_type a_numGhost,
    const bool     a_balance,
    const unsigned a_trim = 0)
    {
      static_assert(HashTableTraits<hashmap_type>::is_static::value,
                    "A synced layout must be static.");

//--Initialization

//**FIXME what about replacing existing layout?  Make sure 'this' is fully
//**      cleared
      m_cachedKey = c_emptyKey;
      m_cachedMapped = nullptr;
      m_bitSets.reset(CH_Cuda::Alloc<hashmap_type>::newPair());
      m_bitSetArray.reset(
        CH_Cuda::Alloc<Array_impl<MDBitSet_type, alloc_type<MDBitSet_type>>>
        ::newPair());
      // m_bitSetArrayCOW not constructed
      m_bitSetArrayCOW.reset();
      m_bitSetData[0] = reinterpret_cast<unsigned char*>(m_bitSetArray->data());
      std::uintptr_t idxDstBitSet = 0;
      m_bitSetData[1] = reinterpret_cast<unsigned char*>(0);
      m_globalInfo.reset(CH_Cuda::Alloc<GlobalInfoCont_type>::newPair());
      // m_globalInfo is defined at comment 'Define m_globalInfo' below
//**FIXME m_globalDispl may need a definition if a_layout is dynamic (see
//**      comments in reset())
//**      It does need a define below.
      m_globalDispl = a_layout.m_globalDispl;
      m_numLocalBS = 0;
      m_numGlobalBS = 0;

//#ifdef CH_MPI
      if (numProc() > 1)
        {

//**Why would we have anything other than local BS?
//**Can't use container if we need the block index

          // uix_type numLocalBS_old = a_layout.size();

          // for (const auto& value : *m_bitSets)
          //   {
          //     if (value.second.m_procID[0] == procID()) ++numLocalBS;
          //   }
          // Start gather of number of BS on each proc

//**Do we really are about this?  Annoying to set globalDispl if a_layout is
//**dynamic
          CH_assert(a_layout.m_numLocalBS == ((*m_globalDispl)[procID() + 1] -
                                              (*m_globalDispl)[procID()]));

//--1a) Communicate number of local BS

          /*
            Note: a separate vector is created here because both sizes and
            displacements are needed for a future MPI call
          */

          std::vector<displ_type> numGlobalBSvec(numProc());
          MPI_Request request;
          MPI_Iallgather(&a_layout.m_numLocalBS, 1,
                         MPI_TYPE<displ_type>::iix_type,
                         numGlobalBSvec.data(),
                         MPI_TYPE<displ_type>::iix_type,
                         Chombo_MPI::comm, &request);

//--2a) Pack local set of global info in preparation for distributing globally

          /*
            Note: m_globalInfo may have been a hash table.  Also, repacking
            avoids overlapping sendbuf with recvbuf.
          */

          std::vector<MDBitSetGlobalInfo_type> localInfo;
          localInfo.reserve(a_layout.m_numLocalBS);
          if (HashTableTraits<SrcHashTable<iixvec_type, MDBitSetInfo_type>>
              ::is_dynamic::value)
            {

              /*
                Note: if the hash table is dynamic, there is an expectation that
                all bitsets are local
              */

              // Copy from an unordered_map
              CH_assert(a_layout.m_bitSets->size() == a_layout.m_numLocalBS);
              for (const auto& value : a_layout->m_bitSets)
                {
                  const auto iter = m_globalInfo->find(value.first);
                  CH_assert(iter != m_globalInfo->end());
                  localInfo.push_back(iter->second);
                }
            }
          else
            {
              // Copy from a vector
              CH_assert(a_layout.m_numLocalBS ==
                        ((*m_globalDispl)[procID() + 1] -
                         (*m_globalDispl)[procID()]));
              for (displ_type i = (*m_globalDispl)[procID()],
                     i_end = (*m_globalDispl)[procID() + 1]; i != i_end; ++i)
                {
                  localInfo.push_back((*a_layout.m_globalInfo)[i]);
                }
            }

//--1b) Finish messages from 1a

          MPI_Wait(&request, MPI_STATUS_IGNORE);

//--1c) We can now configure the displacements for each proc into global info

//**displ = m_globalDispl->data();

          // Count the number of global BS and make sure it is within the limits
          // of an int (otherwise need an MPI implementation supporting 64-bit
          // displacements)
          size_type numGlobalBS = 0;
          for (int idxProc = 0, idxProc_end = numProc(); idxProc != idxProc_end;
               ++idxProc)
            {
              (*m_globalDispl)[idxProc] = (displ_type)numGlobalBS;
              numGlobalBS += (size_type)numGlobalBSvec[idxProc];
            }
          // So you can always find the number of BS on proc i by
          // (*m_globalDispl)[i+1] - (*m_globalDispl)[i]
          (*m_globalDispl)[numProc()] = (displ_type)numGlobalBS;
          // If this fails, the only option is to use an MPI implementation
          // that can support 64-bit displacements
          CH_assert(numGlobalBS <=
                    (size_type)std::numeric_limits<displ_type>::max());
          m_numGlobalBS = numGlobalBS;

//--2a) Communicate global info from local BS info

          // Find displacements and count number of global BS
          constexpr int globalInfoSizeByte =
            MDBitSetGlobalInfo_type::sizeBytes();

//--m_globalInfo is defined here

          // Define m_globalInfo
          m_globalInfo->define(m_numGlobalBS);

          // Gather global info to all
          MPI_Datatype userMPI_GLOBALINFO;
          MPI_Type_contiguous(MDBitSetGlobalInfo_type::sizeBytes(),
                              MPI_BYTE, &userMPI_GLOBALINFO);
          MPI_Type_commit(&userMPI_GLOBALINFO);
          // Note: this is blocking
          MPI_allgatherv(localInfo.data(),
                         a_layout.m_numLocalBS,
                         userMPI_GLOBALINFO,
                         m_globalInfo->data(),
                         numGlobalBSvec.data(),
                         m_globalDispl->data(),
                         userMPI_GLOBALINFO,
                         Chombo_MPI::comm);
          MPI_Type_free(&userMPI_GLOBALINFO);

//--3) Optional (rebalance)

          if (a_balance)
            {
              // Save the original displacements (m_globalDispl becomes dst)

              /*
                Note: the source displacements were originally placed in
                m_globalDispl in case a rebalance was not called
              */

              std::vector<displ_type> srcDispl(numProc() + 1);
              std::memcpy(srcDispl.data(), m_globalDispl->data(),
                          (numProc() + 1)*sizeof(displ_type));

//--3a) Sort indices to the array using Morton ordering.  Have old layout (src)
//--    and new (dst).

              // The indirection array 'balanced' is sorted
              std::vector<displ_type> balanced(m_numGlobalBS);
              for (int i = 0; i != m_numGlobalBS; ++i)
                {
                  balanced[i] = i;
                }
              // Find the average weight and smallest negative key.  All keys
              // need to be positive for Morton ordering otherwise there will be
              // a discontinuity at 0 (for twos complement)
              iixvec_type minNegKey(0);
              double avgWeight = 0.;
              for (const auto& bsInfo : *m_globalInfo)
                {
                  minNegKey.min(bsInfo.m_key);
                  avgWeight += static_cast<double>(bsInfo.m_weight);
                }
              avgWeight /= numProc();
              // Note: if this takes too long, recommend using std::partition
              // and then parallelize partitions.
              std::sort(balanced.begin(), balanced.end(),
                        MortonBalance(m_globalInfo->data(), minNegKey));
              // New processor assignments (intervals of
              // m_globalInfo[balanced[:]])
              const double halfAvgWeight = 0.5*avgWeight;
              // This is adjusted back and forth depending on how much the
              // allocations to the previous processes differed from the average
              // weight
              double adjAvgWeight = avgWeight;
              double totalWeight = 0.;
              displ_type idxBS = 0;
              (*m_globalDispl)[0] = idxBS;
              for (int idxProc = 0, idxProc_end = numProc();
                   idxProc != idxProc_end; ++idxProc)
                {
                  double weight = 0.;
                  while (weight < adjAvgWeight && idxBS < m_numGlobalBS)
                    {
                      weight += static_cast<double>(
                        (*m_globalInfo)[balanced[idxBS++]].weight);
                    }
                  (*m_globalDispl)[idxProc+1] = idxBS;
                  // If this proc has a high weight, reduce the weight for the
                  // next proc to try to get back to average.
                  totalWeight += weight;
                  // If distributions to previous procs has higher than average
                  // weight, reduce the the weight for the next proc (and vice-
                  // versa).
                  adjAvgWeight = std::max(halfAvgWeight,
                                          (idxProc+2)*avgWeight - totalWeight);
//**FIXME print load balancing diagnostics (which is the reason we are keeping
//**      per proc weight and adjAvgWeight versus totals).
                }
//**CHECK how close is idxBS to m_numGlobalBS?
              (*m_globalDispl)[numProc()] = m_numGlobalBS;

              const displ_type srcRangeBeg = srcDispl[procID()];
              const displ_type srcRangeEnd = srcDispl[procID() + 1];
              const displ_type dstRangeBeg = (*m_globalDispl)[procID()];
              const displ_type dstRangeEnd = (*m_globalDispl)[procID()+1];
//**FIXME use m_numLocalBS instead of dstNumLocalBS?
              const displ_type dstNumLocalBS = dstRangeEnd - dstRangeBeg;

              // Reverse lookup of src index and processor for each dst index on
              // this processor
              std::vector<displ_type> invBalanced(dstNumLocalBS);
              std::vector<int> srcProc(dstNumLocalBS);
              for (int idxProc = 0, idxProc_end = numProc();
                   idxProc != idxProc_end; ++idxProc)
                {
                  for (displ_type iSrc = srcDispl[idxProc],
                         iSrc_end = srcDispl[idxProc+1]; iSrc != iSrc_end;
                       ++iSrc)
                    {
                      const displ_type iDst = balanced[iSrc];
                      if (iDst < dstRangeEnd && iDst >= dstRangeBeg)
                        {
                          srcProc[iDst - dstRangeBeg] = idxProc;
                          invBalanced[iDst - dstRangeBeg] = iSrc;
                        }
                    }
                }

//-- 3b) Redistribute the bitsets.  Post sends and receives.

              /*
                Recv: know destination index, source index, and source process
                Send: know source index and destination index (a binary search
                      is used to find the destination process)
              */

              // Post receives
              std::vector<MPI_Request> recvRequests(dstNumLocalBS);
              int idxRecvReq = 0;
              for (displ_type iDst = dstRangeBeg; iDst != dstRangeEnd; ++iDst)
                {
                  const displ_type iSrc = invBalanced[iDst - dstRangeBeg];
                  auto ins = m_bitSets->insert(
                    { m_globalInfo[iSrc].m_key, MDBitSetInfo_type{} });
                  CH_assert(ins.second);  // Insertion must succeed
                  auto iterDst = ins.first;
                  iterDst->second.m_idxBitSet = idxDstBitSet++;
                  iterDst->second.m_procID = procID();
                  iterDst->second.m_tag = iDst - dstRangeBeg;
                  auto iterSrc =
                    a_layout.m_bitSets->find(m_globalInfo[iSrc].m_key);
                  if (iterSrc != a_layout.m_bitSets->end())
                    // This is a local copy
                    {
                      getBitSet(iterDst->second) =
                        a_layout.getBitSet(iterSrc->second);
                    }
                  else
                    // Post receive
                    {
                      MPI_Irecv(getBitSet(iterDst->second).data(),
                                MDBitSet_type::NW,
                                MPI_Type<word_type>::uix_type,
                                srcProc[iDst - dstRangeBeg],
                                iterDst->second.m_tag,
                                Chombo_MPI::comm,
                                &recvRequests[idxRecvReq++]);
                    }
                }

              // Post sends
              std::vector<MPI_Request> sendRequests(a_layout.m_numLocalBS);
              int idxSendReq = 0;
              for (displ_type iSrc = srcRangeBeg; iSrc != srcRangeEnd; ++iSrc)
                {
                  auto iterSrc =
                    a_layout.m_bitSets->find(m_globalInfo[iSrc].m_key);
                  if (iterSrc == a_layout.m_bitSets->end())
                    // Not local copy, must post send
                    {
                      const displ_type iDst = balanced[iSrc];
                      const auto upper = std::upper_bound(m_globalDispl.begin(),
                                                          m_globalDispl.end(),
                                                          iDst);
                      // The processor index is given by the distance of the
                      // iterator
                      const int dstProc = std::distance(m_globalDispl.begin(),
                                                        upper) - 1;
                      const int tag = iDst - dstRangeBeg;
                      MPI_Isend(a_layout.getBitSet(iterSrc->second).data(),
                                MDBitSet_type::NW,
                                MPI_Type<W>::uix_type,
                                dstProc,
                                tag,
                                Chombo_MPI::comm,
                                &sendRequests[idxSendReq++]);
                    }
                }

//--3c) While communicating, reorganize globalInfo

              // Note, this will reset the balance array
              for (displ_type i = 0, i_end = m_numGlobalBS; i != i_end; ++i)
                {
                  if (i != balanced[i])
                    {
                      MDBitSetGlobalInfo_type tmp = m_globalInfo[i];
                      displ_type jp = i;
                      for (displ_type j = balanced[i]; j != i; j = balanced[j])
                        {
                          balanced[jp] = jp;
                          m_globalInfo[jp] = m_globalInfo[j];
                          jp = j;
                        }
                      balanced[jp] = jp;
                      m_globalInfo[jp] = tmp;
                    }
                }

//--3d) Finish all receives

              int mpierr = MPI_Waitall(idxRecvReq, recvRequests,
                                       MPI_STATUSES_IGNORE);
              CH_assert(!mpierr);
HERE
//**FIXME Finish all sends
              

//--3) Populate neighbors of bitsets with meta-data (the bitset itself is still
//--   a null pointer).

          hashmap_type neighborBitSets;
          for (const auto& value : *m_bitSets)
            {
              const iixvec_type baseKey = value.first;
              MDBS::nestedLoop(
                iixvec_type(-1), iixvec_type(2),
                [&baseKey]
                (const iixvec_type& a_ivDir)
                  {
                    const iixvec_type nbrKey = baseKey + a_ivDir;
                    if (m_bitSets->find(nbrKey) == m_bitSets->end())
                      // Not on this processor so this has to be marked as a
                      // candidate neighbor
                      {
                        neighborBitSets.insert({ nbrKey, MDBitSetInfo_type{} });
                      }
                  });
            }
/*
  DELETEME You can't do this kind of detailed matching unless both sides know
  the details of the neighbor.  Otherwise one side might accept the neighbor
  and the other reject it leading to mismatched send/receives.
 *
          hashmap_type neighborBitSets;
          for (const auto& value : *m_bitSets)
            {
              // Using minBox is a bit too crude: two IVs in opposing corners
              // will require ghost from the complete surroundings.  Explore
              // sections based on the neighbor and check to see if range is
              // zero (if so, neglect neighbor).  Check range using mask!
              // const MDBitSetInfo_type& bsInfo = value.second;
              const iixvec_type baseKey = value.first;
              const MDBitSet_type& baseBS = *value.second.m_bitSet;
              MDBS::nestedLoop(
                iixvec_type(-1), iixvec_type(2),
                [&baseKey, &baseBS]
                (const iixvec_type& a_ivDir)
                  {
                    const iixvec_type nbrKey = baseKey + a_ivDir;
                    if (m_bitSets->find(nbrKey) == m_bitSets->end())
                      // Not on this processor
                      {
                        // 0 unless neighbor at upper side, then
                        // dimsb - numGhost
                        const iixvec_type lb = stc::eq(a_ivDir, (iix_type)1)*
                          (MDBitSet_type::c_iix_dimsb - a_numGhost);
                        // dimsb-1 unless neighbor at lower side, then
                        // numGhost-1
                        const iixvec_type ub = MDBitSet_type::c_iix_dimsb - 1 -
                          stc::eq(a_ivDir, (iix_type)-1)*
                          (MDBitSet_type::c_iix_dimsb - a_numGhost);
                        bool haveNeighbor = false;
                        baseBS.forEachCWord(
                          iixvec_type(0), lb, ub,
                          [&haveNeighbor]
                          (const iixvec_type& a_ivw,
                           const word_type&   a_mask,
                           const word_type&   a_word)
                            {
                              if (a_word & a_mask)
                                {
                                  haveNeighbor = true;
                                }
                            });
                        if (haveNeighbor)
                          // Might have neighbor in a_iv direction
                          {
                            const iixvec_type key = value.first + a_ivDir;
                            neighborBitSets.insert(
                              { key, MDBitSetInfo_type{} });
                          }
                      }
                  });
            }
*/


          /*
            Use Morton ordering and move the BS.  Beware in step 5 that source
            locations may change.
           */

//--5) Loop over global BS and match with neighbors

          iix_type numRecv = 0;
          for (iix_type idxProc = 0, idxProc_end = numProc();
               idxProc != idxProc_end; ++idxProc)
            {
              if (idxProc != procID())
                {
                  for (iix_type*
                         keyData =
                         &m_globalKeys[m_globalKeysDispl[idxProc]],
                         keyData_end =
                         &m_globalKeys[m_globalKeysDispl[idxProc + 1]];
                       keyData < keyData_end; keyData += rank)
                    {
                      const iixvec_type key(keyData);
                      auto nbrIter = neighborBitSets.find(key);
                      if (nbrIter != neighborBitSets.end())
                        {
                          ++numRecv;
                          nbrIter->second.m_bitSet = getNewBitSet();
                          nbrIter->second.m_procID[0] = idxProc;
                          nbrIter->second.m_tag = (keyData - m_globalKeys)/rank;
                        }
                    }
                }
              else
                {
                  for (iix_type*
                         keyData =
                         &m_globalKeys[m_globalKeysDispl[idxProc]],
                         keyData_end =
                         &m_globalKeys[m_globalKeysDispl[idxProc + 1]];
                       keyData < keyData_end; keyData += rank)
                    {
                      const iixvec_type key(keyData);
                      auto lclIter = m_bitSets.find(key);
                      CH_assert(lclIter != m_bitSets.end());  // Must exist
                      lclIter->second.m_tag = (keyData - m_globalKeys)/rank;
                    }
                }
              
              // Divide out rank from m_globalKeysDispl to get displ into
              // keys instead of integer elements
              m_globalKeysDispl[idxProc] /= rank;
            }
          m_globalKeysDispl[numProc()] /= rank;

//--4b) Finish messages from 4a

//--6) Send neighbor information (move for rebalance?)

          // To any neighbor, post send and receives
          //**FIXME For now, this does not combine data sent to each processor
          //**FIXME and each bitet is sent/received separately

          // Requests for all messages
          std::vector<MPI_Request> requests(2*numRecv);
          int idxReq = 0;
          // Handles to received bitsets that match requests
          std::vector<hashmap_type::iterator> handles(numRecv);
          // Post receives
          for (hashmap_type::iterator iter = neighborBitSets.begin();
               iter != neighborBitSets.end(); ++iter)
            {
              if (iter->second.m_procID != -1)
                // Expect a receive
                {
                  handles[idxReq] = iter;
                  MPI_Irecv(iter->second.m_bitSet->data(),
                            MDBitSet_type::NW,
                            MPI_Type<W>::uix_type,
                            rmtProcID,
                            iter->second.m_tag,
                            Chombo_MPI::comm,
                            &requests[idxReq++]);
                }
            }
          CH_assert(idxReq == numRecv);

          // Post sends
          iix_type numSend = 0;
          std::set<iixvec_type> rmtProcs;  //**FIXME Reservation would be nice
          for (const auto& value : *m_bitSets)
            {
              rmtProcs.clear();
              const iixvec_type baseKey = value.first;
              MDBS::nestedLoop(
                iixvec_type(-1), iixvec_type(2),
                [this, &&baseKey]
                (const iixvec_type& a_ivDir)
                  {
                    const iixvec_type nbrKey = baseKey + a_ivDir;
                    auto nbrIter = neighborBitSets.find(nbrKey);
                    if (nbrIter != neighborBitSets.end())
                      {
                        const iix_type rmtProcID = nbrIter->second.m_procID;
                        if (rmtProcID != -1 &&
                            rmtProcs.find(rmtProcID) == rmtProcs.end())
                          // Send to this process
                          {
                            rmtProcs.insert(rmtProcID);
                            MPI_Isend(value.second.m_bitSet->data(),
                                      MDBitSet_type::NW,
                                      MPI_Type<W>::uix_type,
                                      rmtProcID,
                                      value.second.m_tag,
                                      Chombo_MPI::comm,
                                      &requests[idxReq++]);
                          }
                      }
                  });
            }

          // Wait for completion
//**FIXME all or process as the receives complete.
          // int mpierr = MPI_Waitany(idxReq,
          //                          requests.data(),
          //                          MPI_STATUSES_IGNORE);
          for (int iReq = 0; iReq != idxReq; ++iReq)
            {
              int ridx;  // Request index
              int mpierr = MPI_Waitany(idxReq,
                                       requests.data(),
                                       &ridx,
                                       MPI_STATUS_IGNORE);
              CH_assert(!mpierr);
              if (ridx < numRecv)
                {
                  hashmap_type::iterator iter = handles[ridx];
                  // Check if we want this neighbor (it has data adjacent within
                  // ghosts to local data).  If we don't want it, delete the
                  // bitset and set m_procID to -1.
                }
            }
          

/*
  How do we match sends and receives?  It is not paired.  Can simply pull if all data is in a window... Or send list of receives.  Or pair afterwards.

  Okay, now recheck neighbors
 */
        }
//#endif
    }
#endif

private:

  /// Fully decouples to an unsynced empty state
  void
  reset()
    {
//**FIXME maybe only allocate if dynamic -- otherwise all sync.  If this is
//**      only dynamic, change newSelectHost to newHost
      m_cachedKey = c_emptyKey;
      m_cachedMapped = nullptr;
      m_bitSets.reset(
        CH_Cuda::Alloc<hashmap_type>::newSelectHost(
          typename HashTableTraits<hashmap_type>::is_dynamic{}));
      // m_bitSetArray not constructed
      m_bitSetArray.reset();
      // m_bitSetArrayCOW not constructed
      m_bitSetArrayCOW.reset();
      m_bitSetData[0] = reinterpret_cast<unsigned char*>(0);
      m_bitSetData[1] = reinterpret_cast<unsigned char*>(0);
      m_globalInfo.reset(
        CH_Cuda::Alloc<GlobalInfoCont_type>::newSelectHost(
          typename HashTableTraits<hashmap_type>::is_dynamic{}));
      m_globalDispl.reset(
        CH_Cuda::Alloc<Array_impl<displ_type, alloc_type<displ_type>>>
        ::newPair());
//**FIXME Why is this necessary for dynamic?  Isn't all useful information in
//**      m_numLocalBS?  It's actually completely useless since m_globalInfo is
//**      likely a hash table.
      m_globalDispl->define(numProc() + 1);
      (*m_globalDispl)[procID()    ] = 0;
      (*m_globalDispl)[procID() + 1] = 0;
      m_numLocalBS = 0;
      m_numGlobalBS = 0;
    }

  /// Get an existing MDBitSet
  HOSTDEVICE MDBitSet_type&
  getBitSet(const MDBitSetInfo_type& a_bsInfo) const noexcept
    {
      return *reinterpret_cast<MDBitSet_type*>(
        m_bitSetData[a_bsInfo.m_copyLevel] + a_bsInfo.m_idxBitSet);
    }

  /// Get an MDBitSet
  //  Consider using pool after performance testing
  std::uintptr_t
  getNewBitSet()
    {
      static_assert(HashTableTraits<hashmap_type>::is_dynamic::value,
                    "New bitsets can only be added to dynamic layouts");
      MDBitSet_type* addr = new MDBitSet_type;
      return reinterpret_cast<std::uintptr_t>(addr);
    }
  
  //**FIXME deleting the bitsets?

  /*
    Single hash map
    Local vector of iterators <- fast compute
    Neighbor vector of iterators

   */

  typename hashmap_type::key_type m_cachedKey;
  typename hashmap_type::mapped_type* m_cachedMapped;

/*
  If this is a pointer, copy on write is possible but you have to swap every
  pointer on the GPU.  If not a pointer, copy on write is not possible.
  Currently implemented as an offset from m_bitSetData
*/
//**FIXME If not a pointer, copy on write is not possible for new schedules
//**FIXME This needs to be convertible
  CH_Cuda::shared_ptr<hashmap_type> m_bitSets;

//--The following are only used for synced hash tables

  // The primary array (for a compacted and synced layout)
  CH_Cuda::shared_ptr<
    Array_impl<MDBitSet_type, alloc_type<MDBitSet_type>>> m_bitSetArray;
  // The secondary copy-on-write array (for a compacted and synced layout).  
  CH_Cuda::shared_ptr<
    Array_impl<MDBitSet_type, alloc_type<MDBitSet_type>>> m_bitSetArrayCOW;
  // Pointers to the beginnings of the above arrays
  unsigned char* m_bitSetData[2];

//--Used for both unsynced and synced...

//**FIXME This should be shared and never change, except that a sync should
//**      fully decouple from parents
  CH_Cuda::shared_ptr<GlobalInfoCont_type> m_globalInfo;
  // CH_Cuda::shared_ptr<Vector<MDBitSetGlobalInfo_type>> m_globalInfo;
                                      ///< On startup, this vector contains
                                      ///< additional information about local BS
                                      ///< including auxiliary user data.  After
                                      ///< sync(), this vector contains all
                                      ///< global information and is used for
                                      ///< load balancing
  /*
    The following are used in MPI where type int is required
  */
//**FIXME This should be shared and never change, except that a sync should
//**      fully decouple from parents
  CH_Cuda::shared_ptr<Array_impl<displ_type, alloc_type<displ_type>>>
    m_globalDispl;
                                      ///< Displacements for each process into
                                      ///< m_globalInfo.  This is always sized
                                      ///< to numProc()+1.  The number of BS
                                      ///< assigned to process i can be
                                      ///< determined from m_globalDispl[i+1] -
                                      ///< m_globalDispl[i].  This must always
                                      ///< be valid, even on startup
  int m_numLocalBS;                   ///< Number of BS local to this processor.
                                      ///< Only set after sync.
  size_type m_numGlobalBS;            ///< Total number of BS across all
                                      ///< processes
  // int status;                         ///< 0 - unsynced primary
  //                                     ///< 1 - synced primary
  //                                     ///< 2 - uncompacted secondary
  //                                     ///< 3 - compacted secondary
};

//--Definitions (concerns about ODR use)

// c_dimsb
template <typename W, typename IIx, typename UIx,
          template <typename, typename> class HashTable,
          template <typename, typename> class AuxInfo, 
          UIx D0b, UIx... Dimsb>
constexpr
typename Layout<W, IIx, UIx, HashTable, AuxInfo, D0b, Dimsb...>::uixvec_type
Layout<W, IIx, UIx, HashTable, AuxInfo, D0b, Dimsb...>::c_dimsb;
// c_log2dimsb
template <typename W, typename IIx, typename UIx,
          template <typename, typename> class HashTable,
          template <typename, typename> class AuxInfo, 
          UIx D0b, UIx... Dimsb>
constexpr
typename Layout<W, IIx, UIx, HashTable, AuxInfo, D0b, Dimsb...>::uixvec_type
Layout<W, IIx, UIx, HashTable, AuxInfo, D0b, Dimsb...>::c_log2dimsb;
// c_iix_keyStep
template <typename W, typename IIx, typename UIx,
          template <typename, typename> class HashTable,
          template <typename, typename> class AuxInfo, 
          UIx D0b, UIx... Dimsb>
constexpr
typename Layout<W, IIx, UIx, HashTable, AuxInfo, D0b, Dimsb...>::iixvec_type
Layout<W, IIx, UIx, HashTable, AuxInfo, D0b, Dimsb...>::c_iix_keyStep;
// c_emptyKey
template <typename W, typename IIx, typename UIx,
          template <typename, typename> class HashTable,
          template <typename, typename> class AuxInfo, 
          UIx D0b, UIx... Dimsb>
constexpr
typename Layout<W, IIx, UIx, HashTable, AuxInfo, D0b, Dimsb...>::iixvec_type
Layout<W, IIx, UIx, HashTable, AuxInfo, D0b, Dimsb...>::c_emptyKey;


}  // namespace HashIVS

#include "NamespaceFooter.H"

#endif  /* ! defined _HASHIVSLAYOUT_H_ */
