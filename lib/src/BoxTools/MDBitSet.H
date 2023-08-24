#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MDBITSET_H_
#define _MDBITSET_H_


/******************************************************************************/
/**
 * \file BitSetFxSz.H
 *
 * \brief
 *
 *//*+*************************************************************************/

#include <bitset>
#include <iomanip>

#include "StcVector.H"

#include "NamespaceHeader.H"

namespace MDBS
{

/// Useful for debugging
const char* writeBitsInInt(const unsigned i)
{
  const int numBit = 8*sizeof(unsigned);
  static char buffer[numBit+1];
  static const char charBit[2] = { '0', '1' };
  for (int n = numBit; n--;)
    {
      buffer[numBit-n-1] = charBit[((i & (1 << n)) != 0)];
    }
  buffer[numBit] = '\0';
  return buffer;
}

/// A nested loop over some integer vector of rank N
template <typename IIx, unsigned I, unsigned rank, typename F>
struct NestedLoop
{
  static constexpr void loop(const stc::Vector<IIx, rank>& lo,
                             const stc::Vector<IIx, rank>& hi,
                             const stc::Vector<IIx, rank>& step,
                             stc::Vector<IIx, rank>& iv,
                             F&& f)
    {
      IIx& i = iv[I-1];
      for (i = lo[I-1]; i < hi[I-1]; i += step[I-1])
        {
          NestedLoop<IIx, I-1, rank, F>::loop(lo, hi, step, iv,
                                              std::forward<F>(f));
        }
    }
};

// Terminate when I = 1 (execute the body of the loop).
template <typename IIx, unsigned rank, typename F>
struct NestedLoop<IIx, 1, rank, F>
{
  static constexpr void loop(const stc::Vector<IIx, rank>& lo,
                             const stc::Vector<IIx, rank>& hi,
                             const stc::Vector<IIx, rank>& step,
                             stc::Vector<IIx, rank>& iv,
                             F&& f)
    {
      IIx& i = iv[0];
      for (i = lo[0]; i < hi[0]; i += step[0])
        {
          f(iv);
        }
    }
};

/// Execute function object 'f(iv)' within nested loop over bounds
/** Note: the bounds are inclusive lo and exclusive hi, lo <= iv < hi
 */
template <typename IIx, unsigned rank, typename F>
constexpr void
nestedLoop(const stc::Vector<IIx, rank>& lo,
           const stc::Vector<IIx, rank>& hi,
           const stc::Vector<IIx, rank>& step,
           F&& f)
{
  stc::Vector<IIx, rank> iv;
  NestedLoop<IIx, rank, rank, F>::loop(lo, hi, step, iv, std::forward<F>(f));
}

/// Fold bits into least significant positions
template <typename T, int FEnd = 4*sizeof(T), int Factor = 1>
struct FoldBits
{
  static constexpr T eval(const T a_word)
    {
      return FoldBits<T, FEnd, 2*Factor>::eval(a_word | (a_word >> Factor));
    }
};

// Terminate when Factor = FEnd
template <typename T, int FEnd>
struct FoldBits<T, FEnd, FEnd>
{
  static constexpr T eval(const T a_word)
    {
      return a_word | (a_word >> FEnd);
    }
};

/// Fold bits into least significant positions
/** E.g.
 *  01000010 Start
 *  01100011 w |= (w >> 1)
 *  01111011 w |= (w >> 2)
 *  01111111 w |= (w >> 4)
 */
template <typename T>
constexpr T
foldBits(const T a_word)
{
  return FoldBits<T>::eval(a_word);
}

}  // namespace MDBS


/*******************************************************************************
 */
/// Fixed-size multi-dimensional bitset
/** \tparam W           Type for word
 *  \tparam IIx         Type for signed integer index
 *  \tparam UIx         Type for unsigned integer index
 *  \tparam D0inWb      Number of bits in dimension zero in a word.  This is
 *                      mostly for internal use and users can always set
 *                      D0inWb == D0b.  One should set D0inWb <= Wszb where
 *                      Wszb is the number of bits in type W.  D0inWb must be
 *                      <= D0b and if <, D0b must be evenly divisible by D0inWb.
 *                      Finally, this value may be neglected, e.g., if
 *                      D0b % Wszb < D0inWb.  One should normally set D0inWb as
 *                      a power of 2
 *  \tparam D0b         The number of bits in the direction having unit stride.
 *                      This should be a power of 2.  If possible, set D0b ==
 *                      Wszb where Wszb is the number of bits in type W.  But do
 *                      not necessarily choose a small type for W just to allow
 *                      this.  E.g., a common 3-D setting is 16x16x16 using
 *                      W=unsigned with 32 bits.
 *  \tparam Dimsb       All remaining dimensions > 0 with the highest having
 *                      the largest stride.
 *
 * Strategy: the multidimensional loops are over words.  For sanity,
 * only the bits in dimension 0 (D0) are returned in a word, if D0b is < Wszb.
 * This is the meaning of c_D0inWszb which is set from parameter D0inWb but can
 * be adjusted due to other considerations.  For example, consider a 2-D bitset
 * with size 16x2 using W=unsigned (32 bits).  A loop returning words would
 * return bit ranges (0:15, 0), (16:31, 0), (0:15, 1), (16:31, 1) in that order.
 * Each range is accessible in the word in bits 0:15.
 *
 * Strategy: for operations the involve neighboring bitsets, a bitset with
 * triple size is created and initialized with all info.  After modifications
 * focused on the central set, the data on the triple is typically or'd back
 * into the original bitsets.
 *
 ******************************************************************************/

template <typename W, typename IIx, typename UIx,  // types
          UIx D0inWb,                              // bits in D0
          UIx D0b, UIx... Dimsb>                   // dimensions
class MDBitSet
{
public:
  // For external use
  using word_type = W;
  using iix_type = IIx;
  using uix_type = UIx;
  static constexpr UIx rank = 1 + sizeof...(Dimsb);

  // Mostly for internal use but made public
  using UIxVec = typename stc::Vector<UIx, rank>;
  using IIxVec = typename stc::Vector<IIx, rank>;
  static constexpr UIxVec c_dimsb{ D0b, Dimsb... };
  static constexpr IIxVec c_iix_dimsb = c_dimsb;
  static constexpr UIxVec c_strideb = stc::cstride(c_dimsb);
  static constexpr IIxVec c_iix_strideb = c_strideb;
  static constexpr UIxVec c_tallyb = c_strideb*c_dimsb;
  static constexpr UIx Nb = c_dimsb.product();
  static constexpr UIx WszB = sizeof(W);
  static constexpr UIx Wszb = 8*WszB;
  static constexpr UIx NW = (Nb + Wszb - 1)/Wszb;
  // Number of bits in last word if not completely filled
  static constexpr UIx Nbrem = Nb % Wszb;
  // Number of surplus (unused) bits
  static constexpr UIx Nbsrp = NW*Wszb - Nb;
  // Assuming at least 1 word covering dimension 0, the number of used bits in
  // that word.
  static constexpr UIx c_D0inWszb =
    std::min(c_dimsb[0] % Wszb, D0inWb) +
    (c_dimsb[0] % Wszb == 0)*std::min(D0inWb, Wszb);
  static constexpr UIxVec c_stepD0inWb =
    stc::make_vector<UIx, rank>::unit() +
    (c_D0inWszb - (UIx)1)*stc::make_vector<UIx, rank>::basis(0);
  static constexpr IIxVec c_iix_stepD0inWb = c_stepD0inWb;
  // Mask for the bits in dim0
  // static constexpr W c_maskD0inW = ((W)1 << c_D0inWszb) - (W)1;
  // Because various compilers complain about shifting more than 8*sizeof(W) in
  // above
  static constexpr W c_maskD0inW =
    (c_D0inWszb < Wszb)*((W)1 << std::min(c_D0inWszb, Wszb - (UIx)1)) - (W)1;
  // The number of words in each direction.  0 if a word doesn't fit into this
  // direction.  This is computed as follows:
  //   if (Wszb > c_tallyb[i])
  //     c_dimsW[i] = 0;
  //   else if (c_strideb[i] < Wszb)
  //     c_dimsW[i] = c_tallyb[i]/Wszb;  // Must fit exactly
  //   else
  //     c_dimsW[i] = c_dimsb[i];
  // which can be expressed as
  //   (Wszb <= c_tallyb)*((c_strideb < Wszb)*(c_tallyb/Wszb) +
  //                       (Wszb <= c_strideb)*c_dimsb);
  static constexpr UIxVec c_dimsW = stc::le(Wszb, c_tallyb)*(
    stc::lt(c_strideb, Wszb)*(c_tallyb/Wszb) +
    stc::le(Wszb, c_strideb)*c_dimsb);
  // The number of words in each direction with 1 or more words in direction 0.
  // During access the bits in direction 0 are copied to a larger word if
  // necessary.  This way, one can assume at least 1 word in direction 0.
  static constexpr UIxVec c_dimsWD0inW = c_dimsb*(
    stc::make_vector<UIx, rank>::unit() -
    stc::make_vector<UIx, rank>::basis(0)) +
    std::max((UIx)1, c_dimsb[0]/Wszb)*
    stc::make_vector<UIx, rank>::basis(0);
#if 0
  // Multiply by index in bits then divide by Wszb to get cW
  static constexpr UIxVec c_strideWallb =
    stc::cstride(stc::max((UIx)1, c_dimsW));
  // Multiply by index in words
  static constexpr UIxVec c_strideW = c_strideWallb*stc::min((UIx)1, c_dimsW);
  // Multiply by index in words
  static constexpr UIxVec c_strideWin0 = stc::cstride(c_dimsWin0);
#endif

  static constexpr UIx size() noexcept
    {
      return Nb;
    }

  static constexpr UIxVec dimensions() noexcept
    {
      return c_dimsb;
    }

  static constexpr UIx dimension(const UIx a_idx) noexcept
    {
      return c_dimsb[a_idx];
    }

  /// Mask a single bit
  template <typename T = W>
  static constexpr T mask1bit(const UIx a_cb) noexcept
    {
      return (T)1 << a_cb;
    }

  /// Mask all bits >= a_cb
  template <typename T = W>
  static constexpr T maskUpper(const UIx a_cb) noexcept
    {
      return ~(mask1bit(a_cb) - (T)1);
    }

  /// Mask all bits <= a_cb
  template <typename T = W>
  static constexpr T maskLower(const UIx a_cb) noexcept
    {
      return mask1bit(a_cb + (T)1) - (T)1;
    }

  /// Returns a word with only the least significant bit from the argument
  template <typename T>
  static constexpr T leastSigBit(const T a_word) noexcept
    {
      return (a_word & -a_word);
    }

  /// Position of the least significant bit
  /** Gives 8*sizeof(T) if a_word = 0
   */
  template <typename T>
  static IIx posLeastSigBit(const T a_word)
    {
      return Misc::popcnt(leastSigBit<T>(a_word) - (T)1);
    }

  /// Pop least significant bit in word and return its position
  /** Gives 8*sizeof(T) if a_word = 0
   */
  template <typename T>
  static IIx posLeastSigBitAndPop(T& a_word) noexcept
    {
      const T i = leastSigBit<T>(a_word);
      a_word -= i;
      return Misc::popcnt(i - (T)1);
    }

  /// Position of the most significant bit (aka floor of log2 with 0 -> 0)
  /** Gives 0 if a_word = 0
   */
  template <typename T>
  static IIx posMostSigBit(const T a_word) noexcept
    {
      const T i = MDBS::foldBits(a_word);
      return Misc::popcnt(i - (T)1);
    }

  MDBitSet()
    {
      clearAll();
    }

  // Allow assignment
  MDBitSet& operator=(const MDBitSet&) = default;

  auto makeThreeWay() const noexcept
    {
      MDBitSet<W, IIx, UIx, c_D0inWszb, 3*D0b, Dimsb...> bs;
      bs.clearAll();
      return bs;
    }

  bool getb(const IIxVec& a_iv) const noexcept
    {
      CH_assert(a_iv >= (stc::make_vector<IIx, rank>::zero()) &&
                a_iv <  c_iix_dimsb);
      const UIx cidx = stc::sum(a_iv*c_iix_strideb);
      const UIx cW = cidx/Wszb;
      W word = m_words[cW];
      const UIx cb = cidx % Wszb;
      return word & mask1bit(cb);
    }

  void setb(const IIxVec& a_iv) noexcept
    {
      CH_assert(a_iv >= (stc::make_vector<IIx, rank>::zero()) &&
                a_iv <  c_iix_dimsb);
      const UIx cidx = stc::sum(a_iv*c_iix_strideb);
      const UIx cW = cidx/Wszb;
      W& word = m_words[cW];
      const UIx cb = cidx % Wszb;
      word |= mask1bit(cb);
    }

  UIx count() const noexcept
    {
      UIx cnt = 0;
      for (IIx cW = 0; cW != (IIx)NW; ++cW)
        {
          W word = m_words[cW];
          cnt += Misc::popcnt(word);
        }
      return cnt;
    }

  void setAll() noexcept
    {
      std::memset(m_words, ~(unsigned char)0, NW*WszB);
      if (Nb % Wszb)
        {
          m_words[NW-1] &= maskLower(Nb % Wszb);
        }
    }

  void clearAll() noexcept
    {
      std::memset(m_words, (unsigned char)0, NW*WszB);
    }

  /// Get a word (mostly used for testing)
  /** Most usage will use loadD0inW instead
   */
  W getW(const UIx a_cW) noexcept
    {
      CH_assert(a_cW < NW);
      return m_words[a_cW];
    }

  /// Load a word with bits max(c_dimsb[0], Wszb)
  /** \tparam W2          Word to store into.  Must have number of bits >=
   *                      min(Wszb, c_D0inWszb).  Recommended to use W.
   *  \param[in]  a_cb    Linear index to a bit in the word
   *  \return             The word with c_D0inWszb bits in range 0:c_D0inWszb-1.
   *                      Other bits are zeroed.
   */
  template <typename W2 = W>
  W2 loadD0inW(const UIx a_cb) noexcept
    {
      const UIx cW = a_cb/Wszb;
      W word = m_words[cW];
      if (c_D0inWszb < Wszb)
        {
          const UIx shift = ((a_cb % Wszb)/c_D0inWszb)*c_D0inWszb;
          word = (word >> shift) & c_maskD0inW;
        }
      return static_cast<W2>(word);
    }

  /// Store (copy) a word with bits max(c_dimsb[0], Wszb)
  /** \param[in] a_cb     Linear index to a bit in the word
   */
  template <typename W2 = W>
  void storeD0inW(const W2 a_word, const UIx a_cb) noexcept
    {
      const UIx cW = a_cb/Wszb;
      W word = static_cast<W>(a_word);
      if (c_D0inWszb < Wszb)
        {
          const UIx shift = ((a_cb % Wszb)/c_D0inWszb)*c_D0inWszb;
          W wordOld = m_words[cW] & ~(c_maskD0inW << shift);  // Clear existing
          word = wordOld | ((word & c_maskD0inW) << shift);   // Copy new
        }
      m_words[cW] = word;
    }

  /// Store (copy) a word with bits max(c_dimsb[0], Wszb) and additional mask
  /** \param[in] a_cb     Linear index to a bit in the word
   */
  template <typename W2 = W>
  void storeD0inW(const W2 a_word, const UIx a_cb, W a_mask) noexcept
    {
      const UIx cW = a_cb/Wszb;
      W word = static_cast<W>(a_word);
      a_mask &= c_maskD0inW;
      const UIx shift = ((a_cb % Wszb)/c_D0inWszb)*c_D0inWszb;
      W wordOld = m_words[cW] & ~(a_mask << shift);  // Clear existing
      word = wordOld | ((word & a_mask) << shift);   // Copy new
      m_words[cW] = word;
    }

  /// Store (or) a word with bits max(c_dimsb[0], Wszb)
  /** \param[in] a_cb     Linear index to a bit in the word
   */
  template <typename W2 = W>
  void orD0inW(const W2 a_word, const UIx a_cb) noexcept
    {
      const UIx cW = a_cb/Wszb;
      W word = static_cast<W>(a_word);
      if (c_D0inWszb < Wszb)
        {
          const UIx shift = ((a_cb % Wszb)/c_D0inWszb)*c_D0inWszb;
          word = (word & c_maskD0inW) << shift;
        }
      m_words[cW] |= word;
    }

  /// Returns 'true' if the entire bitset is zero
  bool isEmpty() const
    {
      for (IIx cW = 0; cW != (IIx)NW; ++cW)
        if (m_words[cW]) return false;
      return true;
    }

  /// Returns 'true' if entire bitset is 1
  bool isFull() const
    {
      constexpr W g = ~(W)0;
      constexpr IIx cW_end = NW - (Nbrem > 0);
      for (IIx cW = 0; cW != cW_end; ++cW)
        if (m_words[cW] != g) return false;
      if (Nbrem)
        {
          if (m_words[NW-1] ^ maskLower(Nbrem)) return false;
        }
      return true;
    }

  /// Union with another bitset
  void operator|=(const MDBitSet& a_bs)
    {
      for (IIx cW = 0; cW != (IIx)NW; ++cW)
        m_words[cW] |= a_bs.m_words[cW];
    }

  /// Complement of argument bitset within this bitset
  void operator-=(const MDBitSet& a_bs)
    {
      for (IIx cW = 0; cW != (IIx)NW; ++cW)
        m_words[cW] &= (~(a_bs.m_words[cW]));
    }

  /// Intersection of argument bitset with this bitset
  void operator&=(const MDBitSet& a_bs)
    {
      for (IIx cW = 0; cW != (IIx)NW; ++cW)
        m_words[cW] &= a_bs.m_words[cW];
    }

  template <typename F>
  void forEachBitInD0inWordLoad(const W a_word, const IIxVec& a_ivw, F&& f)
    {
      IIxVec iv(a_ivw);
      for (IIx cb = 0; cb != (IIx)c_D0inWszb; ++cb)
        {
          const bool getb = a_word & mask1bit(cb);
          iv[0] = a_ivw[0] + cb;
          f(getb, iv);
        }
    }

  template <typename F>
  void forEachBitInD0inWordStore(W& a_word, const IIxVec& a_ivw, F&& f)
    {
      W stoWord = W(0);
      IIxVec iv(a_ivw);
      for (IIx cb = 0; cb != (IIx)c_D0inWszb; ++cb)
        {
          const W mask = mask1bit(cb);
          const bool getb = a_word & mask;
          iv[0] = a_ivw[0] + cb;
          bool setb = f(getb, iv);
          stoWord |= ((W)setb)*mask;
        }
      a_word = stoWord;
    }

  template <typename F>
  void forEachTrueBitInD0inWord(W a_word, const IIxVec& a_ivw, F&& f)
    {
      IIxVec iv(a_ivw);
      while (a_word)
        {
          iv[0] = a_ivw[0] + posLeastSigBitAndPop<W>(a_word);
          f(iv);
        }
    }

  // template <typename F>
  // void forEachTrueBitInWordRange(W word, const stc::Vector<IIx, rank>& a_iv,
  //                                IIx i,
  //                                F&& f)
  //   {
  //     UIx ivw0 = a_iv[0]/c_dimsb[0];
  //     UIx ivwLo0 = (!ivw0)
  //     const UIx cidx = stc::sum(a_iv*c_iix_strideb);
  //     const UIx cW = cidx/Wszb;
  //     W word = getD0InW<W>(a_cW);
  //     stc::Vector<IIx, rank> iv(a_ivw);
  //     word &= (c_maskUpper(iv[0]) & c_maskLower(iv[0]))
  //     while (word)
  //       {
  //         iv[0] = a_ivw[0] + posLeastSigBitAndPop<W>(word);
  //         f(iv);
  //       }
  //   }

  // f(iv)
  template <typename F>
  void forEachTrueBit(const IIxVec& a_baseIV, F&& f)
    {
      // Lambda g finds and iterates on the word
      auto g =
        [this, &a_baseIV, &f]
        (const IIxVec& a_ivw)
          {
            const UIx cb = stc::sum(a_ivw*c_iix_strideb);
            W word = loadD0inW<W>(cb);
            forEachTrueBitInD0inWord<F>(word,
                                        IIxVec(a_baseIV + a_ivw),
                                        std::forward<F>(f));
          };
      MDBS::nestedLoop<IIx, rank>(
        stc::make_vector<IIx, rank>::zero(),
        c_iix_dimsb,
        c_iix_stepD0inWb,
        g);
    }

  template <typename F>
  void forEachTrueBit(const IIxVec& a_baseIV,
                      const IIxVec& a_ivLo, const IIxVec& a_ivHi, F&& f)
    {
      CH_assert(a_ivLo >= (stc::make_vector<IIx, rank>::zero()));
      CH_assert(a_ivHi < c_iix_dimsb);
      CH_assert(a_ivLo <= a_ivHi);
      
      // Lambda g finds and iterates on the word
      auto g =
        [this, &a_baseIV, &f, i0Lo = a_ivLo[0], i0Hi = a_ivHi[0]]
        (IIxVec a_ivb)
          {
            const UIx cb = stc::sum(a_ivb*c_iix_strideb);
            W word = loadD0inW<W>(cb);
            // Find ivw (same as ivb in dims > 0)
            const IIx ivw0 = (a_ivb[0]/c_D0inWszb)*c_D0inWszb;
            // Lop the word
            CH_assert(i0Lo - ivw0 < (IIx)c_D0inWszb);  // By loop bounds
            const W maskU = maskUpper(std::max((IIx)0,
                                               i0Lo - ivw0));
            CH_assert(i0Hi - ivw0 >= (IIx)0);          // By loop bounds
            const W maskL = maskLower(std::min((IIx)c_D0inWszb - (IIx)1,
                                               i0Hi - ivw0));
            word &= (maskU & maskL);
            a_ivb[0] = ivw0;
            // Now iterate on the word
            forEachTrueBitInD0inWord<F>(word,
                                        IIxVec(a_baseIV + a_ivb),
                                        std::forward<F>(f));
          };
      MDBS::nestedLoop<IIx, rank>(
        a_ivLo,
        IIxVec(a_ivHi + stc::make_vector<IIx, rank>::unit()),
        c_iix_stepD0inWb,
        g);
    }

  /// Iterate over all the words (modifiable access)
  /** F
   *  \param[in]  f        Function object with signature
   *                       f(ivw, mask, word)
   *                       where ivw is the IIxVec defining the
   *                       start of the word.  Mask covers all bits for
   *                       dimension 0 in the word.  Word is modifiable.
   *                       Unused bits in the word are zero.  They can
   *                       be modified but this will have no effect.
   */
  template <typename F>
  void forEachWord(const IIxVec& a_baseIV, F&& f)
    {
      // Lambda g finds and iterates on the word
      auto g =
        [this, &a_baseIV, &f]
        (const IIxVec& a_ivw)
          {
            const UIx cb = stc::sum(a_ivw*c_iix_strideb);
            W word = loadD0inW<W>(cb);
            f(IIxVec(a_baseIV + a_ivw), c_maskD0inW, word);
            storeD0inW(word, cb);
          };
      MDBS::nestedLoop<IIx, rank>(
        stc::make_vector<IIx, rank>::zero(),
        c_iix_dimsb,
        c_iix_stepD0inWb,
        g);
    }

  /// Iterate over all the words (constant access)
  /** F
   *  \param[in]  f        Function object with signature
   *                       f(ivw, mask, word)
   *                       where ivw is the IIxVec defining the
   *                       start of the word.  Mask covers all bits for
   *                       dimension 0 in the word.  Word is not
   *                       modifiable.  Unused bits in the word are
   *                       zero.  Bits cannot be modified.
   */
  template <typename F>
  void forEachCWord(const IIxVec& a_baseIV, F&& f)
    {
      // Lambda g finds and iterates on the word
      auto g =
        [this, &a_baseIV, &f]
        (const IIxVec& a_ivw)
          {
            const UIx cb = stc::sum(a_ivw*c_iix_strideb);
            const W word = loadD0inW<W>(cb);
            f(IIxVec(a_baseIV + a_ivw), c_maskD0inW, word);
          };
      MDBS::nestedLoop<IIx, rank>(
        stc::make_vector<IIx, rank>::zero(),
        c_iix_dimsb,
        c_iix_stepD0inWb,
        g);
    }

  /*--------------------------------------------------------------------*/
  /// Iterate over all the words intersecting range (modifiable access)
  /** \param[in]  a_baseIV
   *                      Base IV of the bitset (set to zero for zero-
   *                      based
   *  \param[in]  a_ivLo  Lower bound of range
   *  \param[in]  a_ivHi  Upper bound of range
   *  \param[in]  a_f     Function object with signature
   *                      f(ivw, mask, word) where ivw is the IIxVec
   *                      defining the start of the word.  The mask
   *                      covers only valid bits that are in the range.
   *                      Invalid bits in the word are zeroed.  Valid
   *                      bits that are outside of the range are not
   *                      zeroed.  Bits outside of the range can be
   *                      modified but this will not have any effect.
   */
  template <typename F>
  void forEachWord(const IIxVec& a_baseIV,
                   const IIxVec& a_ivLo, const IIxVec& a_ivHi, F&& f)
    {
      CH_assert(a_ivLo >= (stc::make_vector<IIx, rank>::zero()));
      CH_assert(a_ivHi < c_iix_dimsb);
      CH_assert(a_ivLo <= a_ivHi);
      
      // Lambda g finds and iterates on the word
      auto g =
        [this, &a_baseIV, &f, i0Lo = a_ivLo[0], i0Hi = a_ivHi[0]]
        (IIxVec a_ivb)
          {
            const UIx cb = stc::sum(a_ivb*c_iix_strideb);
            W word = loadD0inW<W>(cb);
            // Find ivw (same as ivb in dims > 0)
            const IIx ivw0 = (a_ivb[0]/c_D0inWszb)*c_D0inWszb;
            // Lop the word
            CH_assert(i0Lo - ivw0 < (IIx)c_D0inWszb);  // By loop bounds
            const W maskU = maskUpper(std::max((IIx)0,
                                               i0Lo - ivw0));
            CH_assert(i0Hi - ivw0 >= (IIx)0);          // By loop bounds
            const W maskL = maskLower(std::min((IIx)c_D0inWszb - (IIx)1,
                                               i0Hi - ivw0));
            const W mask = maskU & maskL;
            a_ivb[0] = ivw0;
            f(IIxVec(a_baseIV + a_ivb), mask, word);
            storeD0inW(word, cb, mask);
          };
      MDBS::nestedLoop<IIx, rank>(
        a_ivLo,
        IIxVec(a_ivHi + stc::make_vector<IIx, rank>::unit()),
        c_iix_stepD0inWb,
        g);
    }

  /*--------------------------------------------------------------------*/
  /// Iterate over all the words intersecting range (constant access)
  /** \param[in]  a_baseIV
   *                      Base IV of the bitset (set to zero for zero-
   *                      based
   *  \param[in]  a_ivLo  Lower bound of range
   *  \param[in]  a_ivHi  Upper bound of range
   *  \param[in]  a_f     Function object with signature
   *                      f(ivw, mask, word) where ivw is the IIxVec
   *                      defining the start of the word.  The mask
   *                      covers only valid bits that are in the range.
   *                      Invalid bits in the word are zeroed.  Valid
   *                      bits that are outside of the range are not
   *                      zeroed.  Bits cannot be modified.
   */
  template <typename F>
  void forEachCWord(const IIxVec& a_baseIV,
                    const IIxVec& a_ivLo, const IIxVec& a_ivHi, F&& f)
    {
      CH_assert(a_ivLo >= (stc::make_vector<IIx, rank>::zero()));
      CH_assert(a_ivHi < c_iix_dimsb);
      CH_assert(a_ivLo <= a_ivHi);

      // Lambda g finds and iterates on the word
      auto g =
        [this, &a_baseIV, &f, i0Lo = a_ivLo[0], i0Hi = a_ivHi[0]]
        (IIxVec a_ivb)
          {
            const UIx cb = stc::sum(a_ivb*c_iix_strideb);
            const W word = loadD0inW<W>(cb);
            // Find ivw (same as ivb in dims > 0)
            const IIx ivw0 = (a_ivb[0]/c_D0inWszb)*c_D0inWszb;
            // Lop the word
            CH_assert(i0Lo - ivw0 < (IIx)c_D0inWszb);  // By loop bounds
            const W maskU = maskUpper(std::max((IIx)0,
                                               i0Lo - ivw0));
            CH_assert(i0Hi - ivw0 >= (IIx)0);          // By loop bounds
            const W maskL = maskLower(std::min((IIx)c_D0inWszb - (IIx)1,
                                               i0Hi - ivw0));
            const W mask = maskU & maskL;
            a_ivb[0] = ivw0;
            f(IIxVec(a_baseIV + a_ivb), mask, word);
          };
      MDBS::nestedLoop<IIx, rank>(
        a_ivLo,
        IIxVec(a_ivHi + stc::make_vector<IIx, rank>::unit()),
        c_iix_stepD0inWb,
        g);
    }

  /// Zero-based lower and upper bounds of minimum box containing all IV
  /** \param[out] a_min   Zero-based lower bound of minbox
   *  \param[out] a_max   Zero-based upper bound of minbox
   *  If no bits are set, a_max = -1 and a_min = dimensions of bitset
   */
  void minBox(IIxVec& a_min, IIxVec& a_max)
    {
      a_min = c_iix_dimsb;
      a_max = (iix_type)(-1);
      forEachCWord(IIxVec(0),
                   [this, &a_min, & a_max]
                   (IIxVec          a_ivw,
                    const word_type a_mask,
                    const word_type a_word)
                     {
                       if (a_word)
                         {
                           const iix_type lo0 =
                             a_ivw[0] + posLeastSigBit(a_word);
                           const iix_type hi0 =
                             a_ivw[0] + posMostSigBit(a_word);
                           a_ivw[0] = lo0;
                           a_min = stc::min(a_min, a_ivw);
                           a_ivw[0] = hi0;
                           a_max = stc::max(a_max, a_ivw);
                         }
                     });
    }

  /// Direct access to the underlying data
  W* data() noexcept
    { return m_words; }

  /// Direct access to the underlying constant data
  const W* data() const noexcept
    { return m_words; }
  
//--Data members

  W m_words[NW];
};

//--Definitions (concerns about ODR use)

// c_dimsb
template <typename W, typename IIx, typename UIx,
          UIx D0inWb,
          UIx D0b, UIx... Dimsb>
constexpr typename MDBitSet<W, IIx, UIx, D0inWb, D0b, Dimsb...>::UIxVec
MDBitSet<W, IIx, UIx, D0inWb, D0b, Dimsb...>::c_dimsb;
// c_iix_dimsb
template <typename W, typename IIx, typename UIx,
          UIx D0inWb,
          UIx D0b, UIx... Dimsb>
constexpr typename MDBitSet<W, IIx, UIx, D0inWb, D0b, Dimsb...>::IIxVec
MDBitSet<W, IIx, UIx, D0inWb, D0b, Dimsb...>::c_iix_dimsb;
// c_strideb
template <typename W, typename IIx, typename UIx,
          UIx D0inWb,
          UIx D0b, UIx... Dimsb>
constexpr typename MDBitSet<W, IIx, UIx, D0inWb, D0b, Dimsb...>::UIxVec
MDBitSet<W, IIx, UIx, D0inWb, D0b, Dimsb...>::c_strideb;
// c_iix_strideb
template <typename W, typename IIx, typename UIx,
          UIx D0inWb,
          UIx D0b, UIx... Dimsb>
constexpr typename MDBitSet<W, IIx, UIx, D0inWb, D0b, Dimsb...>::IIxVec
MDBitSet<W, IIx, UIx, D0inWb, D0b, Dimsb...>::c_iix_strideb;
// c_tallyb
template <typename W, typename IIx, typename UIx,
          UIx D0inWb,
          UIx D0b, UIx... Dimsb>
constexpr typename MDBitSet<W, IIx, UIx, D0inWb, D0b, Dimsb...>::UIxVec
MDBitSet<W, IIx, UIx, D0inWb, D0b, Dimsb...>::c_tallyb;
// c_stepD0inWb
template <typename W, typename IIx, typename UIx,
          UIx D0inWb,
          UIx D0b, UIx... Dimsb>
constexpr typename MDBitSet<W, IIx, UIx, D0inWb, D0b, Dimsb...>::UIxVec
MDBitSet<W, IIx, UIx, D0inWb, D0b, Dimsb...>::c_stepD0inWb;
// c_iix_stepD0inWb
template <typename W, typename IIx, typename UIx,
          UIx D0inWb,
          UIx D0b, UIx... Dimsb>
constexpr typename MDBitSet<W, IIx, UIx, D0inWb, D0b, Dimsb...>::IIxVec
MDBitSet<W, IIx, UIx, D0inWb, D0b, Dimsb...>::c_iix_stepD0inWb;
// c_dimsW
template <typename W, typename IIx, typename UIx,
          UIx D0inWb,
          UIx D0b, UIx... Dimsb>
constexpr typename MDBitSet<W, IIx, UIx, D0inWb, D0b, Dimsb...>::UIxVec
MDBitSet<W, IIx, UIx, D0inWb, D0b, Dimsb...>::c_dimsW;
// c_dimsWD0inW
template <typename W, typename IIx, typename UIx,
          UIx D0inWb,
          UIx D0b, UIx... Dimsb>
constexpr typename MDBitSet<W, IIx, UIx, D0inWb, D0b, Dimsb...>::UIxVec
MDBitSet<W, IIx, UIx, D0inWb, D0b, Dimsb...>::c_dimsWD0inW;

#define MD_WORDLOOP(_bs, _iw, _cw)                                      \
  _bs::uix_type _cw = 0;                                                \
  D_INVTERM6(                                                           \
    for (_bs:iix_type _iw##0 = 0; _iw##0 <  (_bs).c_dimsb[0); _iw##0 += Win0szb, ++ _cw), \
    for (_bs:iix_type _iw##1 = 0; _iw##1 != (_bs).c_dimsb[1); ++ _iw##1), \
    for (_bs:iix_type _iw##2 = 0; _iw##2 != (_bs).c_dimsb[2); ++ _iw##2), \
    for (_bs:iix_type _iw##3 = 0; _iw##3 != (_bs).c_dimsb[3); ++ _iw##3), \
    for (_bs:iix_type _iw##4 = 0; _iw##4 != (_bs).c_dimsb[4); ++ _iw##4), \
    for (_bs:iix_type _iw##5 = 0; _iw##5 != (_bs).c_dimsb[5); ++ _iw##5))

#define MD_GETWORD(word_t, _bs, _iw, _cw)        \
  (_bs).getInW<word_t>(_cw)

#include "NamespaceFooter.H"

#endif  /* ! defined _BITSETFXSZ_H_ */
