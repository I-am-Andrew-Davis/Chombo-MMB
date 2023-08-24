#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _CH_HASH_H_
#define _CH_HASH_H_

#include <cstdint>
#include <limits>
#include <string>
/* Only for debugging */
// #include <iostream>
// #include <iomanip>

#include "google_city.H"
#include "xxhash.H"

#undef CUHOST
#undef HOSTDEVICE
#ifdef __CUDACC__
#define CUHOST __host__
#define HOSTDEVICE __host__ __device__
#else
#define CUHOST
#define HOSTDEVICE
#endif

#include "BaseNamespaceHeader.H"

//--Forward declarations

namespace stc
{
  using array_size_type = unsigned;
  template <typename T,
            array_size_type N,
            template<typename, array_size_type> class Impl>
  class VectorImpl;
  template <array_size_type N,
            template <typename, array_size_type> class Impl>
  class String;
}

namespace CH_Hash
{


/*******************************************************************************
 *
 * IMPORTANT: These hashes work on the consecutive bytes in a type.  If there
 * is any padding, it can make the hash non-deterministic.  Other types,
 * such as DataIndex, may have extraneous, partially filled data such as the
 * data index (versus the layout index).  By default the hash functions fail for
 * all types and must be explicitly enabled by specializing isHashable for the
 * type.  This is to force you to verify that the data layout is safe for
 * hashing.  Alternatively, one can also specialize struct google_CityHash or
 * XXHash for the type.  The latter is required for dynamically sized types such
 * as std::string and stc::String.
 *
 ******************************************************************************/

/// By default, all types cannot be hashed
template <typename T>
struct isHashable : public std::false_type
{
  static constexpr int c_hashSize = 0;
};

// All POD types can be hashed.  Derive from true_type and set c_hashSize.
template <>
struct isHashable<int> : public std::true_type
{
  static constexpr int c_hashSize = sizeof(int);
};

template <>
struct isHashable<unsigned> : public std::true_type
{
  static constexpr int c_hashSize = sizeof(unsigned);
};

template <>
struct isHashable<long long> : public std::true_type
{
  static constexpr int c_hashSize = sizeof(long long);
};

template <>
struct isHashable<long long unsigned> : public std::true_type
{
  static constexpr int c_hashSize = sizeof(long long unsigned);
};


/*******************************************************************************
 *
 * Google city hash collection
 *
 ******************************************************************************/


/*==============================================================================
 * Function objects
 *============================================================================*/


/*==============================================================================
 */
/// Hash object recommended for use with std::unordered_map
/**
 *  Return value is std::size_t and a 32-bit hashing algorithm is always used.
 *  If std::size_t is 64-bit, the hash value is folded into a 64-bit value.
 *  Testing suggests similar hash quality but much faster speed using 32-bit
 *  hash for most Box-based applications.
 *
 *//*=========================================================================*/

template <typename T>
struct google_CityHash
{
  using result_type = std::size_t;

private:

//--Hashes with dispatch based on hash size

  HOSTDEVICE static std::size_t hash(
    const T& a_key,
    std::integral_constant<int, sizeof(uint64_t)>) noexcept
    {
      std::size_t hashval = static_cast<std::size_t>(
        google::CityHash32(reinterpret_cast<const char*>(&a_key),
                           isHashable<T>::c_hashSize));
      /*
        Keep this debugging code which is especially useful for finding
        problems with the key
      */
      // std::cout << "\n*HASH* Size: " << isHashable<T>::c_hashSize
      //           << ", Key: ";
      // const uint8_t* p = reinterpret_cast<const uint8_t*>(&a_key);
      // for (int i = 0; i != isHashable<T>::c_hashSize; ++i)
      //   {
      //     std::cout << ' ' << std::setw(2) << std::setfill('0') << std::hex
      //               << (unsigned)(*p++);
      //   }
      // std::cout << std::dec << ", Hash: " << ((hashval << 32) | hashval)
      //           << std::endl;
      return (hashval << 32) | hashval;
    }
  HOSTDEVICE static std::size_t hash(
    const T&     a_key,
    std::integral_constant<int, sizeof(uint32_t)>) noexcept
    {
      return google::CityHash32(reinterpret_cast<const char*>(&a_key),
                                isHashable<T>::c_hashSize);
    }

public:
  /// Hashing function of the object
  HOSTDEVICE std::size_t operator()(const T& a_key) const noexcept
    {
      static_assert(isHashable<T>::value,
                    "Type not declared hashable with google_CityHash");
      static_assert(isHashable<T>::c_hashSize > 0, "Invalid hash size");
      return hash(a_key,
                  std::integral_constant<int, sizeof(std::size_t)>{});
    }
};

/*--------------------------------------------------------------------*/
/// Specialization for std::string
/*--------------------------------------------------------------------*/

template <>
struct google_CityHash<std::string>
{
  using result_type = std::size_t;

private:

//--Hashes with dispatch based on hash size

  HOSTDEVICE static std::size_t hash(
    const std::string& a_key,
    std::integral_constant<int, sizeof(uint64_t)>) noexcept
    {
      std::size_t hashval = static_cast<std::size_t>(
        google::CityHash32(a_key.c_str(), a_key.size()));
      return (hashval << 32) | hashval;
    }
  HOSTDEVICE static std::size_t hash(
    const std::string& a_key,
    std::integral_constant<int, sizeof(uint32_t)>) noexcept
    {
      return google::CityHash32(a_key.c_str(), a_key.size());
    }

public:

  /// Hashing function of the object
  HOSTDEVICE std::size_t operator()(const std::string& a_key) const noexcept
    {
      return hash(a_key, std::integral_constant<int, sizeof(std::size_t)>{});
    }
};

/*--------------------------------------------------------------------*/
/// Specialization for stc::String
/*--------------------------------------------------------------------*/

template <stc::array_size_type N,
          template <typename, stc::array_size_type> class Impl>
struct google_CityHash<stc::String<N, Impl>>
{
  using result_type = std::size_t;

private:

//--Hashes with dispatch based on hash size

  HOSTDEVICE static std::size_t hash(
    const stc::String<N, Impl>& a_key,
    std::integral_constant<int, sizeof(uint64_t)>) noexcept
    {
      std::size_t hashval = static_cast<std::size_t>(
        google::CityHash32(a_key.c_str(), a_key.length()));
      return (hashval << 32) | hashval;
    }
  HOSTDEVICE static std::size_t hash(
    const stc::String<N, Impl>& a_key,
    std::integral_constant<int, sizeof(uint32_t)>) noexcept
    {
      return google::CityHash32(a_key.c_str(), a_key.length());
    }

public:

  /// Hashing function of the object
  HOSTDEVICE std::size_t operator()(
    const stc::String<N, Impl>& a_key) const noexcept
    {
      return hash(a_key, std::integral_constant<int, sizeof(std::size_t)>{});
    }
};


/*==============================================================================
 */
/// City-Hash object using 64-bit hash and returning 64-bit value
/**
 *//*=========================================================================*/

template <typename T>
struct google_CityHash64
{
  using result_type = uint64_t;
  CUHOST uint64_t operator()(const T& a_key) const noexcept
    {
      static_assert(isHashable<T>::value,
                    "Type not declared hashable with google_CityHash");
      static_assert(isHashable<T>::c_hashSize > 0, "Invalid hash size");
      return google::CityHash64(reinterpret_cast<const char*>(&a_key),
                                isHashable<T>::c_hashSize);
    }
};


/*==============================================================================
 */
/// City-Hash object using 32-bit hash and returning 32-bit value
/**
 *//*=========================================================================*/

template <typename T>
struct google_CityHash32
{
  using result_type = uint32_t;
  HOSTDEVICE uint32_t operator()(const T& a_key) const noexcept
    {
      static_assert(isHashable<T>::value,
                    "Type not declared hashable with google_CityHash");
      static_assert(isHashable<T>::c_hashSize > 0, "Invalid hash size");
      return google::CityHash32(reinterpret_cast<const char*>(&a_key),
                                isHashable<T>::c_hashSize);
    }
};


/*==============================================================================
 */
/// City-Hash object using 32-bit hash and returning 64-bit value
/**
 *  The 32-bit hash is folded into a 64-bit value
 *
 *//*=========================================================================*/

template <typename T>
struct google_CityHash32_r64
{
  using result_type = uint64_t;
  HOSTDEVICE uint64_t operator()(const T& a_key) const noexcept
    {
      static_assert(isHashable<T>::value,
                    "Type not declared hashable with google_CityHash");
      static_assert(isHashable<T>::c_hashSize > 0, "Invalid hash size");
      uint64_t hashval = static_cast<uint64_t>(
        google::CityHash32(reinterpret_cast<const char*>(&a_key),
                           isHashable<T>::c_hashSize));
      return (hashval << 32) | hashval;
    }
};


/*==============================================================================
 * Helper functions
 *============================================================================*/

/// City-Hash function returning hash value into std::size_t
/** Note that if std::size_t is 64-bit, the 32-bit hash is folded into a 64-bit
 *  value
 */
template <typename T>
HOSTDEVICE std::size_t
hash_google_CityHash(const T& a_key) noexcept
{
  return (google_CityHash<T>{}).operator()(a_key);
}

/// Genuine 64-bit City-Hash of some key
template <typename T>
CUHOST uint64_t
hash_google_CityHash64(const T& a_key) noexcept
{
  return (google_CityHash64<T>{}).operator()(a_key);
}

/// 32-bit City-Hash of some key
template <typename T>
HOSTDEVICE uint32_t
hash_google_CityHash32(const T& a_key) noexcept
{
  return (google_CityHash32<T>{}).operator()(a_key);
}


/*******************************************************************************
 *
 * XXHash collection
 *
 ******************************************************************************/


/*==============================================================================
 * Function objects
 *============================================================================*/

/*==============================================================================
 */
/// Hash object recommended for use with std::unordered_map
/**
 *  Return value is std::size_t and a 32-bit hashing algorithm is always used.
 *  If std::size_t is 64-bit, the hash value is folded into a 64-bit value.
 *
 *//*=========================================================================*/

template <typename T>
struct XXHash
{
  using result_type = std::size_t;

private:

//--Hashes with dispatch based on hash size

  HOSTDEVICE static std::size_t hash(
    const T&       a_key,
    const uint32_t a_seed,
    std::integral_constant<int, sizeof(uint64_t)>) noexcept
    {
      std::size_t hashval = static_cast<std::size_t>(
        xxhash::XXH32(reinterpret_cast<const char*>(&a_key),
                      isHashable<T>::c_hashSize,
                      a_seed));
      /*
        Keep this debugging code which is especially useful for finding
        problems with the key
      */
      // std::cout << "\n*HASH* Size: " << isHashable<T>::c_hashSize
      //           << ", Key: ";
      // const uint8_t* p = reinterpret_cast<const uint8_t*>(&a_key);
      // for (int i = 0; i != isHashable<T>::c_hashSize; ++i)
      //   {
      //     std::cout << ' ' << std::setw(2) << std::setfill('0') << std::hex
      //               << (unsigned)(*p++);
      //   }
      // std::cout << std::dec << ", Hash: " << ((hashval << 32) | hashval)
      //           << std::endl;
      return (hashval << 32) | hashval;
    }
  HOSTDEVICE static std::size_t hash(
    const T&       a_key,
    const uint32_t a_seed,
    std::integral_constant<int, sizeof(uint32_t)>) noexcept
    {
      return xxhash::XXH32(reinterpret_cast<const char*>(&a_key),
                           isHashable<T>::c_hashSize,
                           a_seed);
    }

public:
  /// Hashing function of the object
  HOSTDEVICE std::size_t operator()(const T&       a_key,
                                    const uint32_t a_seed = 0u) const noexcept
    {
      static_assert(isHashable<T>::value,
                    "Type not declared hashable with google_CityHash");
      static_assert(isHashable<T>::c_hashSize > 0, "Invalid hash size");
      return hash(a_key,
                  a_seed,
                  std::integral_constant<int, sizeof(std::size_t)>{});
    }
};

/*--------------------------------------------------------------------*/
/// Specialization for std::string
/*--------------------------------------------------------------------*/

template <>
struct XXHash<std::string>
{
  using result_type = std::size_t;

private:

//--Hashes with dispatch based on hash size

  HOSTDEVICE static std::size_t hash(
    const std::string& a_key,
    const uint32_t     a_seed,
    std::integral_constant<int, sizeof(uint64_t)>) noexcept
    {
      std::size_t hashval = static_cast<std::size_t>(
        xxhash::XXH32(a_key.c_str(), a_key.size(), a_seed));
      return (hashval << 32) | hashval;
    }
  HOSTDEVICE static std::size_t hash(
    const std::string& a_key,
    const uint32_t     a_seed,
    std::integral_constant<int, sizeof(uint32_t)>) noexcept
    {
      return xxhash::XXH32(a_key.c_str(), a_key.size(), a_seed);
    }

public:

  /// Hashing function of the object
  HOSTDEVICE std::size_t operator()(
    const std::string& a_key,
    const uint32_t     a_seed = 0u) const noexcept
    {
      return hash(a_key,
                  a_seed,
                  std::integral_constant<int, sizeof(std::size_t)>{});
    }
};

/*--------------------------------------------------------------------*/
/// Specialization for stc::String
/*--------------------------------------------------------------------*/

template <stc::array_size_type N,
          template <typename, stc::array_size_type> class Impl>
struct XXHash<stc::String<N, Impl>>
{
  using result_type = std::size_t;

private:

//--Hashes with dispatch based on hash size

  HOSTDEVICE static std::size_t hash(
    const stc::String<N, Impl>& a_key,
    const uint32_t              a_seed,
    std::integral_constant<int, sizeof(uint64_t)>) noexcept
    {
      std::size_t hashval = static_cast<std::size_t>(
        xxhash::XXH32(a_key.c_str(), a_key.length(), a_seed));
      return (hashval << 32) | hashval;
    }
  HOSTDEVICE static std::size_t hash(
    const stc::String<N, Impl>& a_key,
    const uint32_t              a_seed,
    std::integral_constant<int, sizeof(uint32_t)>) noexcept
    {
      return xxhash::XXH32(a_key.c_str(), a_key.length(), a_seed);
    }

public:

  /// Hashing function of the object
  HOSTDEVICE std::size_t operator()(
    const stc::String<N, Impl>& a_key,
    const uint32_t              a_seed = 0u) const noexcept
    {
      return hash(a_key,
                  a_seed,
                  std::integral_constant<int, sizeof(std::size_t)>{});
    }
};


/*==============================================================================
 */
/// xxhash object using 32-bit hash and returning 32-bit value
/**
 *//*=========================================================================*/

template <typename T>
struct XXHash32
{
  using result_type = uint32_t;
  HOSTDEVICE uint32_t operator()(const T&       a_key,
                                 const uint32_t a_seed = 0u) const noexcept
    {
      static_assert(isHashable<T>::value,
                    "Type not declared hashable with xxhash");
      static_assert(isHashable<T>::c_hashSize > 0, "Invalid hash size");
      return xxhash::XXH32(reinterpret_cast<const char*>(&a_key),
                           isHashable<T>::c_hashSize,
                           a_seed);
    }
};

/*--------------------------------------------------------------------*/
/// Specialization for std::string
/*--------------------------------------------------------------------*/

template <>
struct XXHash32<std::string>
{
  using result_type = uint32_t;
  HOSTDEVICE uint32_t operator()(
    const std::string& a_key,
    const uint32_t     a_seed = 0u) const noexcept
    {
      return xxhash::XXH32(a_key.c_str(), a_key.size(), a_seed);
    }
};

/*--------------------------------------------------------------------*/
/// Specialization for stc::String
/*--------------------------------------------------------------------*/

template <stc::array_size_type N,
          template <typename, stc::array_size_type> class Impl>
struct XXHash32<stc::String<N, Impl>>
{
  using result_type = uint32_t;
  HOSTDEVICE uint32_t operator()(
    const stc::String<N, Impl>& a_key,
    const uint32_t              a_seed = 0u) const noexcept
    {
      return xxhash::XXH32(a_key.c_str(), a_key.length(), a_seed);
    }
};


/*==============================================================================
 */
/// Google City-Hash object using 32-bit hash and returning 64-bit value
/** 
 *  The 32-bit hash is folded into a 64-bit value
 *
 *//*=========================================================================*/

template <typename T>
struct XXHash32_r64
{
  using result_type = uint64_t;
  HOSTDEVICE uint64_t operator()(const T&       a_key,
                                 const uint32_t a_seed = 0u) const noexcept
    {
      static_assert(isHashable<T>::value,
                    "Type not declared hashable with xxhash");
      static_assert(isHashable<T>::c_hashSize > 0, "Invalid hash size");
      uint64_t hashval = static_cast<uint64_t>(
        xxhash::XXH32(reinterpret_cast<const char*>(&a_key),
                      isHashable<T>::c_hashSize,
                      a_seed));
      return (hashval << 32) | hashval;
    }
};

/*--------------------------------------------------------------------*/
/// Specialization for std::string
/*--------------------------------------------------------------------*/

template <>
struct XXHash32_r64<std::string>
{
  using result_type = uint64_t;
  HOSTDEVICE uint64_t operator()(
    const std::string& a_key,
    const uint32_t     a_seed = 0u) const noexcept
    {
      std::size_t hashval = static_cast<std::size_t>(
        xxhash::XXH32(a_key.c_str(), a_key.size(), a_seed));
      return (hashval << 32) | hashval;
    }
};

/*--------------------------------------------------------------------*/
/// Specialization for stc::String
/*--------------------------------------------------------------------*/

template <stc::array_size_type N,
          template <typename, stc::array_size_type> class Impl>
struct XXHash32_r64<stc::String<N, Impl>>
{
  using result_type = uint64_t;
  HOSTDEVICE uint64_t operator()(
    const stc::String<N, Impl>& a_key,
    const uint32_t              a_seed = 0u) const noexcept
    {
      std::size_t hashval = static_cast<std::size_t>(
        xxhash::XXH32(a_key.c_str(), a_key.length(), a_seed));
      return (hashval << 32) | hashval;
    }
};

/*==============================================================================
 * Helper functions
 *============================================================================*/

/// xxhash function returning hash value into std::size_t
/** Note that if std::size_t is 64-bit, the 32-bit hash is folded into a 64-bit
 *  value
 */
template <typename T>
HOSTDEVICE std::size_t
hash_XXHash(const T& a_key, const uint32_t a_seed = 0u) noexcept
{
  return (XXHash<T>{}).operator()(a_key);
}

/// Genuine 64-bit City-Hash of some key
template <typename T>
CUHOST uint64_t
hash_XXHash64(const T& a_key, const uint32_t a_seed = 0u) noexcept
{
  return (XXHash32_r64<T>{}).operator()(a_key);
}

/// 32-bit City-Hash of some key
template <typename T>
HOSTDEVICE uint32_t
hash_XXHash32(const T& a_key, const uint32_t a_seed = 0u) noexcept
{
  return (XXHash32<T>{}).operator()(a_key);
}

}  // namespace CH_Hash

#include "BaseNamespaceFooter.H"

#endif  /* _CH_HASH_H_ */
