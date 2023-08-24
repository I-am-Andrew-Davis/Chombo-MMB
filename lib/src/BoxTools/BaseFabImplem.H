#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BASEFABIMPLEM_H_
#define _BASEFABIMPLEM_H_

#include <utility>

#include "BaseFabMacros.H"
#include "CH_Timer.H"
#include "SliceSpec.H"
#include "MayDay.H"
#include "parstream.H"
#include "NamespaceHeader.H"

// template <typename Alloc>  ??FIXME
// int BaseFab<int, Alloc>::test();

template <typename T, typename Alloc>
int
BaseFab<T, Alloc>::testBoxAndComp()
{
  MayDay::Warning("pretty minimal test");
  Box b(IntVect::Zero, IntVect::Unit);
  BaseFab<T> blerg;
  blerg.define(b, 1);

  if (blerg.nComp() != 1)
    {
      pout() << "ncomp busted" << endl;
      return -2;
    }
  if (blerg.box() != b)
    {
      pout() << "box return busted" << endl;
      return -2;
    }
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (blerg.size()[idir] != 2)
        {
          pout() << "size busted" <<endl;
          return -3;
        }
    }

  if (blerg.smallEnd() != IntVect::Zero)
    {
      pout() << "smallend busted" <<endl;
      return -4;
    }
  if (blerg.bigEnd() != IntVect::Unit)
    {
      pout() << "bigend busted" <<endl;
      return -5;
    }
  return 0;
}///

template <typename T, typename Alloc>
int
BaseFab<T, Alloc>::test()
{
  int retval = testBoxAndComp();
  return retval;
}

//
// Implementation.=====================================
//

template <typename T, typename Alloc>
inline
BaseFab<T, Alloc>::BaseFab(const unsigned a_allocOn)
  :
  m_domain(Box()),
  m_nvar(0),
  m_numpts(0),
  m_baseFabImpl(a_allocOn)
{
}

template <typename T, typename Alloc>
template <typename... Args>
inline
BaseFab<T, Alloc>::BaseFab(const unsigned a_allocOn,
                           const Box&     a_bx,
                           int            a_n,
                           T*             a_alias,
                           Args&&...      a_args)
  :
  m_domain(a_bx),
  m_nvar(a_n),
  m_numpts(a_bx.numPts()),
  m_baseFabImpl(a_allocOn)
{
  if (a_alias == nullptr)
    {
      m_baseFabImpl.define(m_numpts*m_nvar, std::forward<Args>(a_args)...);
    }
  else
    {
      m_baseFabImpl.defineAlias(a_alias, m_numpts*m_nvar);
    }
  T* mddata = (m_baseFabImpl.isUsable()) ?
    m_baseFabImpl.dataPtr() + m_domain.index(IntVect_zero) : nullptr;
  const IntVect boxDims = m_domain.size();
  m_mdarray.define(mddata, MD_EXPANDIV(boxDims), m_nvar);
  if (m_baseFabImpl.isUsable())
    {
      m_mdarray.setLBIndices(MD_EXPANDIV(m_domain.smallEnd()), 0);
    }
}

// alias
template <typename T, typename Alloc>
inline
BaseFab<T, Alloc>::BaseFab(const Interval& a_comps,
                           BaseFab<T>&     a_original)
  :
  m_domain(a_original.m_domain),
  m_nvar(a_comps.size()),
  m_numpts(a_original.m_numpts),
  m_baseFabImpl()
{
  m_baseFabImpl.defineAlias(a_original.dataPtr(a_comps.begin()),
                            m_numpts*m_nvar);
  T* mddata = (m_baseFabImpl.isUsable()) ?
    m_baseFabImpl.dataPtr() + m_domain.index(IntVect_zero) : nullptr;
  const IntVect boxDims = m_domain.size();
  m_mdarray.define(mddata, MD_EXPANDIV(boxDims), m_nvar);
  if (m_baseFabImpl.isUsable())
    {
      m_mdarray.setLBIndices(MD_EXPANDIV(m_domain.smallEnd()), 0);
    }
}

//--Move constructors

template <typename T, typename Alloc>
inline
BaseFab<T, Alloc>::BaseFab(BaseFab&& a_fab) noexcept
  :
  m_domain(std::move(a_fab.m_domain)),
  m_nvar(a_fab.m_nvar),
  m_numpts(a_fab.m_numpts),
  m_baseFabImpl(std::move(a_fab.m_baseFabImpl)),
  m_mdarray(a_fab.m_mdarray)
{ }

template <typename T, typename Alloc>
inline BaseFab<T, Alloc>&
BaseFab<T, Alloc>::operator=(BaseFab&& a_fab) noexcept
{
  if (&a_fab != this)
    {
      m_domain = std::move(a_fab.m_domain);
      m_nvar = a_fab.m_nvar;
      m_numpts = a_fab.m_numpts;
      m_baseFabImpl = std::move(a_fab.m_baseFabImpl);
      m_mdarray = a_fab.m_mdarray;
    }
  return *this;
}

// Destructor
template <typename T, typename Alloc>
inline
BaseFab<T, Alloc>::~BaseFab()
{
  // Destruction of m_baseFabImpl will release memory if allocated
}

template <typename T, typename Alloc>
template <typename... Args>
inline void
BaseFab<T, Alloc>::resize(const Box& a_b,
                          int        a_n,
                          T*         a_alias,
                          Args&&...  a_args)
{
  m_domain = a_b;
  m_nvar   = a_n;
  m_numpts = m_domain.numPts();
  if (a_alias == nullptr)
  {
    m_baseFabImpl.define(m_numpts*m_nvar, std::forward<Args>(a_args)...);
  }
  else
  {
    m_baseFabImpl.defineAlias(a_alias, m_numpts*m_nvar);
  }
  T* mddata = (m_baseFabImpl.isUsable()) ?
    m_baseFabImpl.dataPtr() + m_domain.index(IntVect_zero) : nullptr;
  const IntVect boxDims = m_domain.size();
  m_mdarray.define(mddata, MD_EXPANDIV(boxDims), m_nvar);
  if (m_baseFabImpl.isUsable())
    {
      m_mdarray.setLBIndices(MD_EXPANDIV(m_domain.smallEnd()), 0);
    }
}

template <typename T, typename Alloc>
inline void
BaseFab<T, Alloc>::define(const Interval& a_comps,
                          BaseFab<T>&     a_original)
{
  m_domain   = a_original.m_domain;
  m_nvar     = a_comps.size();
  m_numpts   = a_original.m_numpts;
  m_baseFabImpl.defineAlias(a_original.dataPtr(a_comps.begin()),
                            m_numpts*m_nvar);
  T* mddata = (m_baseFabImpl.isUsable()) ?
    m_baseFabImpl.dataPtr() + m_domain.index(IntVect_zero) : nullptr;
  const IntVect boxDims = m_domain.size();
  m_mdarray.define(mddata, MD_EXPANDIV(boxDims), m_nvar);
  if (m_baseFabImpl.isUsable())
    {
      m_mdarray.setLBIndices(MD_EXPANDIV(m_domain.smallEnd()), 0);
    }
}

#ifdef USE_PROTO
template < class T>
template<unsigned int C>
void BaseFab<T>::define(Proto::BoxData<T,C,1,1>&   a_original)
{
    undefine();
    m_domain = getBox(a_original.box());
    m_numpts = a_original.box().size();
    m_truesize = a_original.size();
    m_nvar = C;
    m_dptr = &a_original(a_original.box().low());
    m_aliased = true;
}
#endif

template <typename T, typename Alloc>
inline void
BaseFab<T, Alloc>::clear()
{
  m_domain = Box();
  m_nvar   = 0;
  m_numpts = 0;
  m_baseFabImpl.undefine();
  m_mdarray.clear();
}

template <typename T, typename Alloc>
inline int
BaseFab<T, Alloc>::nComp() const
{
  return m_nvar;
}

template <typename T, typename Alloc>
inline const Box&
BaseFab<T, Alloc>::box() const
{
  return m_domain;
}

template <typename T, typename Alloc>
inline IntVect
BaseFab<T, Alloc>::size() const
{
  return m_domain.size();
}

template <typename T, typename Alloc>
inline const IntVect&
BaseFab<T, Alloc>::smallEnd() const
{
  return m_domain.smallEnd();
}

template <typename T, typename Alloc>
inline const IntVect&
BaseFab<T, Alloc>::bigEnd() const
{
  return m_domain.bigEnd();
}

template <typename T, typename Alloc>
inline T&
BaseFab<T, Alloc>::operator()(const IntVect& a_p,
                              int            a_n)
{
  CH_assert(a_n >= 0);
  CH_assert(a_n < m_nvar);
  CH_assert(m_baseFabImpl.isUsable());
  CH_assert(m_domain.contains(a_p));

  const size_type idx =
    (size_type)m_domain.index(a_p) + m_numpts * (size_type)a_n;

  return m_baseFabImpl[idx];
}

template <typename T, typename Alloc>
inline T&
BaseFab<T, Alloc>::operator()(const IntVect& a_p)
{
  CH_assert(m_baseFabImpl.isUsable());
  CH_assert(m_domain.contains(a_p));

  return m_baseFabImpl[m_domain.index(a_p)];
}

template <typename T, typename Alloc>
inline const T&
BaseFab<T, Alloc>::operator()(const IntVect& a_p,
                              int            a_n) const
{
  CH_assert(a_n >= 0);
  CH_assert(a_n < m_nvar);
  CH_assert(m_baseFabImpl.isUsable());
  CH_assert(m_domain.contains(a_p));

  const size_type idx =
    (size_type)m_domain.index(a_p) + m_numpts * (size_type)a_n;

  return m_baseFabImpl[idx];
}

template <typename T, typename Alloc>
inline const T&
BaseFab<T, Alloc>::operator()(const IntVect& a_p) const
{
  CH_assert(m_baseFabImpl.isUsable());
  CH_assert(m_domain.contains(a_p));

  return m_baseFabImpl[m_domain.index(a_p)];
}

template <typename T, typename Alloc>
inline void
BaseFab<T, Alloc>::getVal(T*             a_data,
                          const IntVect& a_pos,
                          int            a_n,
                          int            a_numcomp) const
{
  CH_assert(m_baseFabImpl.isUsable());
  CH_assert(a_n >= 0 && a_n + a_numcomp <= m_nvar);

  const size_type loc  = m_domain.index(a_pos);
  const size_type size = m_domain.numPts();

  for (int k = 0; k < a_numcomp; k++)
  {
    a_data[k] = m_baseFabImpl[loc+(a_n+k)*size];
  }
}

template <typename T, typename Alloc>
inline void
BaseFab<T, Alloc>::getVal(T*             a_data,
                          const IntVect& a_pos) const
{
  getVal(a_data,a_pos,0,m_nvar);
}

template <typename T, typename Alloc>
inline const int*
BaseFab<T, Alloc>::loVect() const
{
  return m_domain.loVect();
}

template <typename T, typename Alloc>
inline const int*
BaseFab<T, Alloc>::hiVect() const
{
  return m_domain.hiVect();
}

template <typename T, typename Alloc>
inline const int*
BaseFab<T, Alloc>::nCompPtr() const
{
  CH_assert(m_baseFabImpl.isUsable());

  return &m_nvar;
}

template <typename T, typename Alloc>
inline T*
BaseFab<T, Alloc>::dataPtr(int a_n)
{
  CH_assert(m_baseFabImpl.isUsable());
  CH_assert((a_n >= 0) && (a_n < m_nvar));

  return &m_baseFabImpl[a_n * m_numpts];
}

template <typename T, typename Alloc>
inline const T*
BaseFab<T, Alloc>::dataPtr(int a_n) const
{
  CH_assert(m_baseFabImpl.isUsable());
  CH_assert((a_n >= 0) && (a_n < m_nvar));

  return &m_baseFabImpl[a_n * m_numpts];
}

template <typename T, typename Alloc>
inline bool
BaseFab<T, Alloc>::contains(const BaseFab<T>& a_fab) const
{
  return box().contains(a_fab.box()) && m_nvar <= a_fab.m_nvar;
}

template <typename T, typename Alloc>
inline bool
BaseFab<T, Alloc>::contains (const Box& a_bx) const
{
  return box().contains(a_bx);
}

template <typename T, typename Alloc>
inline void
BaseFab<T, Alloc>::setVal(T          a_x,
                          const Box& a_bx,
                          int        a_nstart,
                          int        a_numcomp)
{
  performSetVal(a_x,a_bx,a_nstart,a_numcomp);
}

template <typename T, typename Alloc>
inline void
BaseFab<T, Alloc>::setVal(T          a_x,
                          const Box& a_bx,
                          int        a_n)
{
  performSetVal(a_x,a_bx,a_n,1);
}

template <typename T, typename Alloc>
inline void
BaseFab<T, Alloc>::setVal(T   a_x,
                          int a_n)
{
  performSetVal(a_x,m_domain,a_n,1);
}

template <typename T, typename Alloc>
inline void
BaseFab<T, Alloc>::setVal(T a_x)
{
  setVal(typename std::is_pod<T>(), a_x);
}

template <typename T, typename Alloc>
inline void
BaseFab<T, Alloc>::setVal(std::true_type, T a_x)
{
  if ((T)0 == a_x)
    {
      std::memset(dataPtr(), 0, m_baseFabImpl.size()*sizeof(T));
    }
  else
    {
      performSetVal(a_x,box(),0,m_nvar);
    }
}

template <typename T, typename Alloc>
inline void
BaseFab<T, Alloc>::setVal(std::false_type, T a_x)
{
  performSetVal(a_x,box(),0,m_nvar);
}

template <typename T, typename Alloc>
inline BaseFab<T>&
BaseFab<T, Alloc>::copy(const BaseFab<T>& a_src,
                        const Box&        a_srcbox,
                        int               a_srccomp,
                        const Box&        a_destbox,
                        int               a_destcomp,
                        int               a_numcomp)
{
  CH_assert(a_srcbox.sameSize(a_destbox));
  CH_assert(a_src.box().contains(a_srcbox));
  CH_assert(m_domain.contains(a_destbox));
  CH_assert(a_srccomp >= 0 && a_srccomp+a_numcomp <= a_src.nComp());
  CH_assert(a_destcomp >= 0 && a_destcomp+a_numcomp <= m_nvar);

  performCopy(a_src,a_srcbox,a_srccomp,a_destbox,a_destcomp,a_numcomp);

  return *this;
}

template <typename T, typename Alloc>
inline BaseFab<T>&
BaseFab<T, Alloc>::copy(const BaseFab<T>& a_src,
                        int               a_srccomp,
                        int               a_destcomp,
                        int               a_numcomp)
{
  CH_assert(a_srccomp  >= 0 && a_srccomp  + a_numcomp <= a_src.m_nvar);
  CH_assert(a_destcomp >= 0 && a_destcomp + a_numcomp <= m_nvar);

  Box overlap(m_domain);
  overlap &= a_src.m_domain;

  if (!overlap.isEmpty())
  {
    performCopy(a_src,overlap,a_srccomp,overlap,a_destcomp,a_numcomp);
  }

  return *this;
}

template <typename T, typename Alloc>
inline BaseFab<T>&
BaseFab<T, Alloc>::copy(const BaseFab<T>& a_src,
                        const Box&        a_destbox)
{
  CH_assert(m_nvar <= a_src.m_nvar);
  CH_assert(m_domain.contains(a_destbox));

  Box overlap(a_destbox);
  overlap &= a_src.m_domain;

  if (!overlap.isEmpty())
  {
    performCopy(a_src,overlap,0,overlap,0,m_nvar);
  }

  return *this;
}

template <typename T, typename Alloc>
inline BaseFab<T>&
BaseFab<T, Alloc>::copy(const BaseFab<T>& a_src)
{
  CH_assert(m_nvar <= a_src.m_nvar);
  CH_assert(m_domain.sameType(a_src.m_domain));

  Box overlap(m_domain);
  overlap &= a_src.m_domain;

  if (!overlap.isEmpty())
  {
    performCopy(a_src,overlap,0,overlap,0,m_nvar);
  }

  return *this;
}

template <typename T, typename Alloc>
inline
void BaseFab<T, Alloc>::copy(const Box&        a_RegionFrom,
                             const Interval&   a_Cdest,
                             const Box&        a_RegionTo,
                             const BaseFab<T>& a_src,
                             const Interval&   a_Csrc)
{
  if ((this == &a_src) && (a_RegionFrom == a_RegionTo) && (a_Cdest == a_Csrc) )
  {
    return;
  }

  CH_assert(a_Cdest.size() == a_Csrc.size());

  copy(a_src, a_RegionFrom, a_Csrc.begin(), a_RegionTo,
       a_Cdest.begin(), a_Cdest.size());
}

template <typename T, typename Alloc>
inline BaseFab<T>&
BaseFab<T, Alloc>::shift(const IntVect& a_v)
{
  CH_assert(m_baseFabImpl.isUsable());
  m_domain += a_v;
  m_mdarray.resetData(m_baseFabImpl.dataPtr() + m_domain.index(IntVect_zero));
  m_mdarray.setLBIndices(MD_EXPANDIV(m_domain.smallEnd()), 0);
  return *this;
}

template <typename T, typename Alloc>
inline BaseFab<T>&
BaseFab<T, Alloc>::shift(int a_idir,
                         int a_ncells)
{
  CH_assert(m_baseFabImpl.isUsable());
  m_domain.shift(a_idir,a_ncells);
  m_mdarray.resetData(m_baseFabImpl.dataPtr() + m_domain.index(IntVect_zero));
  m_mdarray.setLBIndices(MD_EXPANDIV(m_domain.smallEnd()), 0);
  return *this;
}

template <typename T, typename Alloc>
inline BaseFab<T>&
BaseFab<T, Alloc>::shiftHalf(int a_idir,
                             int a_numHalfs)
{
  CH_assert(m_baseFabImpl.isUsable());
  m_domain.shiftHalf(a_idir,a_numHalfs);
  m_mdarray.resetData(m_baseFabImpl.dataPtr() + m_domain.index(IntVect_zero));
  m_mdarray.setLBIndices(MD_EXPANDIV(m_domain.smallEnd()), 0);
  return *this;
}

template <typename T, typename Alloc>
inline BaseFab<T>&
BaseFab<T, Alloc>::shiftHalf(const IntVect& a_v)
{
  CH_assert(m_baseFabImpl.isUsable());
  m_domain.shiftHalf(a_v);
  m_mdarray.resetData(m_baseFabImpl.dataPtr() + m_domain.index(IntVect_zero));
  m_mdarray.setLBIndices(MD_EXPANDIV(m_domain.smallEnd()), 0);
  return *this;
}

template <typename T, typename Alloc>
inline auto
BaseFab<T, Alloc>::size(const Box&      a_box,
                        const Interval& a_comps) const
  -> size_type
{
  return a_box.numPts() * (sizeof(T) * a_comps.size());
}

template <typename T, typename Alloc>
inline void
BaseFab<T, Alloc>::linearOut(void*           a_buf,
                             const Box&      a_R,
                             const Interval& a_comps) const
{
  linearOut2(a_buf, a_R, a_comps);
}

template <typename T, typename Alloc>
inline void*
BaseFab<T, Alloc>::linearOut2(void*           a_buf,
                              const Box&      a_R,
                              const Interval& a_comps) const
{
  CH_assert(Misc::isAligned<T>(a_buf));
  T* buffer = static_cast<T*>(a_buf);
  for (int iComp = a_comps.begin(), iComp_end = a_comps.end() + 1;
       iComp != iComp_end; ++iComp)
    {
      MD_BOXLOOP(a_R, i)
        {
          *buffer++ = m_mdarray[MD_IX(i, iComp)];
        }
    }
  return static_cast<void*>(buffer);
}

template <typename T, typename Alloc>
inline void
BaseFab<T, Alloc>::linearIn(void*           a_buf,
                            const Box&      a_R,
                            const Interval& a_comps)
{
  linearIn2(a_buf, a_R, a_comps);
}

template <typename T, typename Alloc>
inline void*
BaseFab<T, Alloc>::linearIn2(void*           a_buf,
                             const Box&      a_R,
                             const Interval& a_comps)
{
  CH_assert(Misc::isAligned<T>(a_buf));
  T* buffer = static_cast<T*>(a_buf);
  for (int iComp = a_comps.begin(), iComp_end = a_comps.end() + 1;
       iComp != iComp_end; ++iComp)
    {
      MD_BOXLOOP(a_R, i)
        {
          m_mdarray[MD_IX(i, iComp)] = *buffer++;
        }
    }
  return static_cast<void*>(buffer);
}

//--Linear in/out for broadcast and gather; the Box is in the msg

template <typename T, typename Alloc>
inline void
BaseFab<T, Alloc>::linearOut(void* a_buf) const
{
  char* bytebuf = static_cast<char*>(a_buf);

  CH_XD::linearOut(bytebuf, m_domain);
  bytebuf += CH_XD::linearSize(m_domain);

  CH_XD::linearOut(bytebuf, m_nvar);
  bytebuf += sizeof(int) ;

  linearOut(bytebuf, m_domain, interval());
}

template <typename T, typename Alloc>
inline void
BaseFab<T, Alloc>::linearIn(const void* const a_buf)
{
  // first get the box, then the number of components, resize, then the data
  const char* bytebuf = static_cast<const char*>(a_buf);

  Box domain;
  CH_XD::linearIn(domain, bytebuf);
  bytebuf += CH_XD::linearSize(domain);

  int ncomps;
  CH_XD::linearIn(ncomps, bytebuf);
  bytebuf += sizeof(ncomps);

  resize(domain, ncomps);

  CH_assert(m_domain == domain && m_nvar == ncomps);
  linearIn(bytebuf, domain, ncomps);
}

template <typename T, typename Alloc>
inline int
BaseFab<T, Alloc>::linearSize() const
{
  // the linearization contains the Box, #components, and data
  return CH_XD::linearSize(m_domain) + sizeof(int) + size(m_domain, interval());
}

template <typename T, typename Alloc>
inline std::string
BaseFab<T, Alloc>::name()
{
  std::string rtn = (typeid(T)).name();

  return rtn;
}

template <typename T, typename Alloc>
inline void
BaseFab<T, Alloc>::performCopy(const BaseFab<T>& a_src,
                               const Box&        a_srcbox,
                               int               a_srccomp,
                               const Box&        a_destbox,
                               int               a_destcomp,
                               int               a_numcomp)
{
  CH_assert(a_src.box().contains(a_srcbox));
  CH_assert(box().contains(a_destbox));
  CH_assert(a_destbox.sameSize(a_srcbox));
  CH_assert(a_srccomp  >= 0 && a_srccomp  + a_numcomp <= a_src.nComp());
  CH_assert(a_destcomp >= 0 && a_destcomp + a_numcomp <= nComp());
  // CH_TIME("BaseFab::performCopy")
  for (int iCdst = a_destcomp, iCdst_end = a_destcomp + a_numcomp,
         iCsrc = a_srccomp; iCdst != iCdst_end; ++iCdst, ++iCsrc)
    {
      MD_2BOXLOOP(a_destbox, idst, a_srcbox, isrc)
        {
          m_mdarray[MD_IX(idst, iCdst)] = a_src[MD_IX(isrc, iCsrc)];
        }
    }
}

template <typename T, typename Alloc>
inline void
BaseFab<T, Alloc>::performSetVal(T          a_x,
                                 const Box& a_bx,
                                 int        a_nstart,
                                 int        a_numcomp)
{
  CH_assert(m_domain.contains(a_bx));
  CH_assert(a_nstart >= 0);
  CH_assert(a_nstart + a_numcomp <= m_nvar);

  if (a_bx == m_domain)
    {
      T* data = dataPtr(a_nstart);
      for (size_type n = a_numcomp*m_numpts; n--;)
        {
          *data++ = a_x;
        }
    }
  else
    {
      for (int iC = a_nstart, iC_end = a_nstart + a_numcomp; iC != iC_end; ++iC)
        {
          MD_BOXLOOP(a_bx, i)
            {
              m_mdarray[MD_IX(i, iC)] = a_x;
            }
        }
    }
}

template <typename T, typename Alloc>
inline bool
BaseFab<T, Alloc>::isAliased() const
{
  return m_baseFabImpl.allocBy() == Array_impl<T, Alloc>::AllocBy::alias;
}

template <typename T, typename Alloc>
inline void
BaseFab<T, Alloc>::degenerate(BaseFab<T>&      a_slice,
                              const SliceSpec& a_sliceSpec ) const
{
  bool outofbounds;
  Box degenerateBox;
  this->box().degenerate( degenerateBox, a_sliceSpec, &outofbounds );
  if ( outofbounds )
    {
      MayDay::Error( "Tried to slice out-of-bounds." );
    }
  a_slice.define( degenerateBox, this->nComp() );
  a_slice.copy( *this, degenerateBox );
}

#include "NamespaceFooter.H"

#endif
