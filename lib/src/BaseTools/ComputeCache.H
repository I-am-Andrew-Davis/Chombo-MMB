#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _CH_COMPUTECACHE_H_
#define _CH_COMPUTECACHE_H_


/******************************************************************************/
/**
 * \file ComputeCache.H
 *
 * \brief Implements caching of a computed value based on key
 *
 *//*+*************************************************************************/

#include <unordered_map>

#include "CH_Hash.H"
#include "BaseNamespaceHeader.H"


/*******************************************************************************
 */
/// Implements caching of value based on key using an unordered map
/**
 *  \tparam Key         Key for lookup in cache
 *  \tparam T           Value returned from lookup in cache
 *  \tparam MakeT       MakeT::operator() used to construct T
 *  \tparam Hash        Method for hashing keys
 *  \tparam KeyEqual    Method for testing if keys are equal
 *
 *  This class provides 'T' for a given 'Key'.  It is expected that T is
 *  expensive to produce and is used repeatedly.  The main intent of this class
 *  is to provide a mechanism for caching T and retrieving it.  Most often, 'T'
 *  is somehow constructed from 'Key'.  However, this is not required and 'T'
 *  can be constructed independent from 'Key' as MakeT::operator()(a_makeTArgs).
 *  The onus is on the user for maintaining some relation between 'Key' and
 *  'a_makeTArgs'.  You are allowed to modify the retrieved cached values at
 *  your own risk.
 *
 *  MakeT is expected to be a lambda constructued as a temporary.  The object
 *  is passed in as a universal reference meaning it's type should be deduced.
 *  It is recommended to use the helper make_ComputeCacheUMap routines.  See
 *  testCache.cpp for an example of using the constructor directly.  If the
 *  MakeT object is an rvalue, it is stored herein; if it is an lvalue, only a
 *  reference is stored and the object must exist for as long as this class.
 *
 *//*+*************************************************************************/

template <typename Key,
          typename T,
          typename MakeT,
          typename Hash = CH_Hash::google_CityHash<Key>,
          typename KeyEqual = std::equal_to<Key>>
class ComputeCacheUMap
{

//--Types

public:

  using key_type   = Key;
  using value_type = T;
  using maker_type = MakeT;
  using cache_type = std::unordered_map<Key, T, Hash, KeyEqual>;

//--Constructors

  /// Constructor takes MakeT object (usually lambda or std::function)
  //  Note: it is recommended to use make_ComputeCacheUMap routines instead.
  ComputeCacheUMap(MakeT&& a_makeT)
    :
    m_makeT(std::forward<MakeT>(a_makeT))
    { }

  // Use synthesized copy, move, and destructor
  // Copy and move assignment are deleted because they are probably not
  // supported if make_type is a lambda
  ComputeCacheUMap(const ComputeCacheUMap&)            = default;
  ComputeCacheUMap(ComputeCacheUMap&&)                 = default;
  ComputeCacheUMap& operator=(const ComputeCacheUMap&) = delete;
  ComputeCacheUMap& operator=(ComputeCacheUMap&&)      = delete;
  ~ComputeCacheUMap()                                  = default;

//--Member functions

  /// Clear the cache
  void clear()
    { m_cache.clear(); }

  /// Request a value for a key
  template <typename... Args>
  auto operator()(const Key& a_key, Args&&... a_makeTargs)
    -> value_type&;

  /// Request a value for a key (if 'a_makeTargs' is same as 'a_key')
  auto operator()(const Key& a_key)
    -> value_type&;

  /// Access the maker
  maker_type& get_makeT()
    { return m_makeT; }

  /// Acces the cache
  cache_type& get_cache()
    { return m_cache; }

//--Data members

protected:

  maker_type m_makeT;                 ///< Operator() makes value_type from
                                      ///< either key_type or specific arguments
  cache_type m_cache;                 ///< Cache of value_type from key_type
};


/*******************************************************************************
 *
 * Class ComputeCacheUMap: construction helpers
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
/// Make a ComputeCacheUMap object with default hash and equal
/** Do not specify the last template so that it can be deduced.
 *  \param[in]  a_makeT a_makeT::operator() used to generate value
 *  \return             Constructed cache
 *//*-----------------------------------------------------------------*/

template <typename Key, typename T, typename MakeT>
auto make_ComputeCacheUMap(MakeT&& a_makeT)
{
  return ComputeCacheUMap<Key, T, MakeT>(std::forward<MakeT>(a_makeT));
}

/*--------------------------------------------------------------------*/
/// Make a ComputeCacheUMap object specifying hash and equal classes
/** Do not specify the last template so that it can be deduced.
 *  \param[in]  a_makeT a_makeT::operator() used to generate value
 *  \return             Constructed cache
 *//*-----------------------------------------------------------------*/

template <typename Key,
          typename T,
          typename Hash,
          typename KeyEqual,
          typename MakeT>
auto make_ComputeCacheUMap(MakeT&& a_makeT)
{
  return ComputeCacheUMap<Key, T, MakeT, Hash, KeyEqual>(
    std::forward<MakeT>(a_makeT));
}

/*******************************************************************************
 *
 * Class ComputeCacheUMap: inline member definitions
 *
 ******************************************************************************/

/*--------------------------------------------------------------------*/
/// Request a value for a key
/** \param[in]  a_key   Key to lookup
 *  \param[in]  a_makeTargs
 *                      Arguments that are passed to function
 *                      m_makeT::operator() if the value must be made.
 *                      This may be the same as 'a_key', in which case
 *                      the sister routine can be used.  See
 *                      testCache.cpp for an example of using a
 *                      std::function pointing to a member fuction
 *                      where both the class object and key must be
 *                      passed as arguments.
 *  \return             The value, either retrieved from cache or made
 *                      from m_makeT::operator()(a_makeTargs) and now
 *                      stored in cache.
 *
 *  T (aka value_type) does not need to be default-constructible.
 *//*-----------------------------------------------------------------*/

template <typename Key,
          typename T,
          typename Hash,
          typename KeyEqual,
          typename MakeT>
template <typename... Args>
inline auto
ComputeCacheUMap<Key, T, Hash, KeyEqual, MakeT>::operator()(
  const Key& a_key,
  Args&&...  a_makeTargs)
  -> value_type&
{
  typename cache_type::iterator iter = m_cache.find(a_key);
  if (iter == m_cache.end())  // Key not found, must create T
    {
      auto ins = m_cache.insert(
        { a_key, m_makeT(std::forward<Args>(a_makeTargs)...) });
      iter = ins.first;
      CH_assert(ins.second);  // Insert must have succeeded.
    }
  return iter->second;
}

/*--------------------------------------------------------------------*/
/// Request a value for a key
/** \param[in]  a_key   Key to lookup and sole argument to
 *                      m_makeT::operator()
 *  \return             The value, either retrieved from cache or made
 *                      from m_makeT::operator()(a_key) and now
 *                      stored in cache.
 *
 *  T (aka value_type) does not need to be default-constructible.
 *//*-----------------------------------------------------------------*/

template <typename Key,
          typename T,
          typename Hash,
          typename KeyEqual,
          typename MakeT>
inline auto
ComputeCacheUMap<Key, T, Hash, KeyEqual, MakeT>::operator()(
  const Key& a_key)
  -> value_type&
{
  return this->operator()(a_key, a_key);
}

#include "BaseNamespaceFooter.H"

#endif  /* _CH_CACHE_H_ */
