#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _CH_BOXHASH_H_
#define _CH_BOXHASH_H_

#include <cstdint>
#include <limits>

#include "CH_Hash.H"
#include "Box.H"

#include "NamespaceHeader.H"

namespace CH_Hash
{


/*******************************************************************************
 *
 * Google City Hash collection, specializations for box objects
 *
 ******************************************************************************/


/*==============================================================================
 * Function objects
 *============================================================================*/

// Specialize for Box
template <>
struct google_CityHash<Box>
{
private:
  static std::size_t hash(const stc::IVec<2*SpaceDim>& a_key,
                          std::integral_constant<int, sizeof(uint64_t)>)
    {
      std::size_t hashval = static_cast<std::size_t>(
        google::CityHash32(reinterpret_cast<const char*>(a_key.dataPtr()),
                           a_key.size()*sizeof(int)));
      return (hashval << 32) | hashval;
    }
  static std::size_t hash(const stc::IVec<2*SpaceDim>& a_key,
                          std::integral_constant<int, sizeof(uint32_t)>)
    {
      return google::CityHash32(reinterpret_cast<const char*>(a_key.dataPtr()),
                                a_key.size()*sizeof(int));
    }
public:
  std::size_t operator()(const Box& a_box) const
    {
      stc::IVec<2*SpaceDim> key;
      for (const int dir : EachDir)
        {
          key[dir]          = a_box.smallEnd(dir);
          key[dir+SpaceDim] = a_box.bigEnd(dir);
        }
      return hash(key, std::integral_constant<int, sizeof(std::size_t)>{});
    }
};

}  // namespace CH_Hash

#include "NamespaceFooter.H"

#endif  /* _CH_BOXHASH_H_ */
