#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _CH_BOXKEYS_H_
#define _CH_BOXKEYS_H_

#include "DataIndex.H"


/*******************************************************************************
 *
 * Combinations of BoxTools objects as keys for hashing
 *
 ******************************************************************************/


/*==============================================================================
 * Define the number of bytes that are hashed for each type while also
 * indicating that these types can be byte-hashed.  This is done first so
 * the sizes can be runtime-asserted in the class definitions
 *============================================================================*/

#include "NamespaceHeader.H"
// Forward declarations
namespace CH_BoxKeys
{
  struct BoxInt;
  struct DataIndex;
  struct DataIndexInt;
}
#include "NamespaceFooter.H"

#include "BaseNamespaceHeader.H"
#include "UsingNamespace.H"
namespace CH_Hash
{

// Forward declaration
template <typename T>
struct isHashable;

/// Remark that BoxInt can be byte-hashed and give bytes to hash
template <>
struct isHashable<CH_BoxKeys::BoxInt> : public std::true_type
{
  static_assert(sizeof(IndexType) == sizeof(unsigned char),
                "Cannot guarantee that IndexType is not padded");
  static constexpr int c_hashSize =
    2*isHashable<IntVect>::c_hashSize + sizeof(int) + sizeof(IndexType);
};

/// Remark that CH_BoxKeys::DataIndex can be byte-hashed and give bytes to hash
template <>
struct isHashable<CH_BoxKeys::DataIndex> : public std::true_type
{
  static constexpr int c_hashSize = sizeof(int*) + sizeof(int);
};

/// Remark that DataIndexInt can be byte-hashed and give bytes to hash
template <>
struct isHashable<CH_BoxKeys::DataIndexInt> : public std::true_type
{
  static constexpr int c_hashSize = sizeof(int*) + 2*sizeof(int);
};

}  // namespace CH_Hash
#include "BaseNamespaceFooter.H"

#include "NamespaceHeader.H"

/*==============================================================================
 * The key objects.  Each must assert that there is no padding (rearrange data
 * as required).
 *============================================================================*/

namespace CH_BoxKeys
{

/// A Box and an integer
struct BoxInt
{
  BoxInt(const Box& a_box, int a_int)
    :
    m_lo(a_box.smallEnd()),
    m_hi(a_box.bigEnd()),
    m_int(a_int),
    m_ixType(a_box.ixType())
    {
      // Ensure no padding in BoxInt
      CH_assert(reinterpret_cast<char*>(this) +
                CH_Hash::isHashable<CH_BoxKeys::BoxInt>::c_hashSize ==
                reinterpret_cast<char*>(&m_ixType) + sizeof(IndexType));
    }
  // All other copy, move, assignment, destructor are synthesized
  BoxInt(const BoxInt&) = default;
  BoxInt(BoxInt&&) = default;
  BoxInt& operator=(const BoxInt&) = default;
  BoxInt& operator=(BoxInt&&) = default;
  ~BoxInt() = default;
  // Data
  IntVect   m_lo;
  IntVect   m_hi;
  int       m_int;
  IndexType m_ixType;
};

/// A DataIndex
struct DataIndex
{
  DataIndex(const ::DataIndex& a_didx)
    :
    m_intPtr(a_didx.layoutIntPtr()),
    m_intCode(a_didx.intCode())
    {
      // Ensure no padding in DataIndex
      CH_assert(reinterpret_cast<char*>(this) +
                CH_Hash::isHashable<CH_BoxKeys::DataIndex>::c_hashSize ==
                reinterpret_cast<char*>(&m_intCode) + sizeof(int));
    }
  // All other copy, move, assignment, destructor are synthesized
  DataIndex(const DataIndex&) = default;
  DataIndex(DataIndex&&) = default;
  DataIndex& operator=(const DataIndex&) = default;
  DataIndex& operator=(DataIndex&&) = default;
  ~DataIndex() = default;
  // Data
  const int* m_intPtr;
  int m_intCode;
};

/// A DataIndex and an integer
struct DataIndexInt
{
  DataIndexInt(const ::DataIndex& a_didx, int a_int)
    :
    m_intPtr(a_didx.layoutIntPtr()),
    m_intCode(a_didx.intCode()),
    m_int(a_int)
    { 
      // Ensure no padding in DataIndexInt
      CH_assert(reinterpret_cast<char*>(this) +
                CH_Hash::isHashable<CH_BoxKeys::DataIndexInt>::c_hashSize ==
                reinterpret_cast<char*>(&m_int) + sizeof(int));
    }
  // All other copy, move, assignment, destructor are synthesized
  DataIndexInt(const DataIndexInt&) = default;
  DataIndexInt(DataIndexInt&&) = default;
  DataIndexInt& operator=(const DataIndexInt&) = default;
  DataIndexInt& operator=(DataIndexInt&&) = default;
  ~DataIndexInt() = default;
  // Data
  const int* m_intPtr;
  int m_intCode;
  int m_int;
};


/*==============================================================================
 * Binary operator == for each key
 *============================================================================*/

/// Equals for BoxInt
inline bool
operator==(const CH_BoxKeys::BoxInt& a_x, const CH_BoxKeys::BoxInt& a_y)
{
  return (a_x.m_lo == a_y.m_lo &&
          a_x.m_hi == a_y.m_hi &&
          a_x.m_int == a_y.m_int &&
          a_x.m_ixType == a_y.m_ixType);
}

/// Equals for DataIndex
inline bool
operator==(const CH_BoxKeys::DataIndex& a_x,
           const CH_BoxKeys::DataIndex& a_y)
{
  return (a_x.m_intPtr  == a_y.m_intPtr &&
          a_x.m_intCode == a_y.m_intCode);
}

/// Equals for DataIndexInt
inline bool
operator==(const CH_BoxKeys::DataIndexInt& a_x,
           const CH_BoxKeys::DataIndexInt& a_y)
{
  return (a_x.m_intPtr  == a_y.m_intPtr &&
          a_x.m_intCode == a_y.m_intCode &&
          a_x.m_int     == a_y.m_int);
}

}  // namespace CH_BoxKeys

#include "NamespaceFooter.H"

#endif  /* _CH_BOXKEYS_H_ */
