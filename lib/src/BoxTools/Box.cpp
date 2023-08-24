#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <climits>

#include "MayDay.H"
#include "Box.H"
#include "parstream.H"
#include "SliceSpec.H"
#include "NamespaceHeader.H"



using std::endl;
using std::ws;
using std::ostream;
using std::istream;

// const Box Box::Empty = Box(IntVect::Unit,
//                            IntVect::Zero);

CH_XDIR::IntVect
convertOldToNew(const IntVect& a_ivOld,
                const IntVect& a_permutation,
                const IntVect& a_sign,
                const IntVect& a_translation)
{
  IntVect ivNew;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      ivNew[idir] =
        a_sign[idir]*a_ivOld[a_permutation[idir]] + a_translation[idir];
    }
  return ivNew;
}


CH_XDIR::IntVect
convertNewToOld(const IntVect& a_ivNew,
                const IntVect& a_permutation,
                const IntVect& a_sign,
                const IntVect& a_translation)
{
  IntVect ivOld;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      ivOld[a_permutation[idir]] = a_sign[idir] *
        (a_ivNew[idir] - a_translation[idir]);
    }
  return ivOld;
}


///multiblock stuff.
void
Box::
convertOldToNew(const IntVect& a_permutation,
                const IntVect& a_sign,
                const IntVect& a_translation)
{
  IntVect ivNewLo = CH_XDIR::convertOldToNew(smallEnd(), a_permutation, a_sign, a_translation);
  IntVect ivNewHi = CH_XDIR::convertOldToNew(bigEnd()  , a_permutation, a_sign, a_translation);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      int iLo = Min(ivNewLo[idir], ivNewHi[idir]);
      int iHi = Max(ivNewLo[idir], ivNewHi[idir]);
      ivNewLo[idir] = iLo;
      ivNewHi[idir] = iHi;
    }
  //  Box bxNewNodes(ivNewLo, ivNewHi, IndexType::TheNodeType());
  //  Box bxNewCells = enclosedCells(bxNewNodes);
  *this = Box(ivNewLo, ivNewHi);
}


///multiblock stuff
void
Box::
convertNewToOld(const IntVect& a_permutation,
                const IntVect& a_sign,
                const IntVect& a_translation)
{
  IntVect ivOldLo = CH_XDIR::convertNewToOld(smallEnd(), a_permutation, a_sign, a_translation);
  IntVect ivOldHi = CH_XDIR::convertNewToOld(bigEnd()  , a_permutation, a_sign, a_translation);
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      int iLo = Min(ivOldLo[idir], ivOldHi[idir]);
      int iHi = Max(ivOldLo[idir], ivOldHi[idir]);
      ivOldLo[idir] = iLo;
      ivOldHi[idir] = iHi;
    }
  //  Box bxOldNodes(ivOldLo, ivOldHi, IndexType::TheNodeType());
  //  Box bxOldCells = enclosedCells(bxOldNodes);
  *this =  Box(ivOldLo, ivOldHi);
}

IndexType
IndexType::TheCellType ()
{
    static const IndexType Cell(D_DECL6(IndexType::CELL,
                                        IndexType::CELL,
                                        IndexType::CELL,
                                        IndexType::CELL,
                                        IndexType::CELL,
                                        IndexType::CELL));
    return Cell;
}

IndexType
IndexType::TheNodeType ()
{
    static const IndexType Node(D_DECL6(IndexType::NODE,
                                        IndexType::NODE,
                                        IndexType::NODE,
                                        IndexType::NODE,
                                        IndexType::NODE,
                                        IndexType::NODE));

    return Node;
}

ostream&
operator<< (ostream&         os,
            const IndexType& it)
{
    os << '('
       << D_TERM6( (it.test(0)?'N':'C'),
                   << ',' << (it.test(1)?'N':'C'),
                   << ',' << (it.test(2)?'N':'C'),
                   << ',' << (it.test(3)?'N':'C'),
                   << ',' << (it.test(4)?'N':'C'),
                   << ',' << (it.test(5)?'N':'C')) << ')' ;
    if (os.fail())
        MayDay::Error("operator<<(ostream&,IndexType&) failed");

    return os;
}

//
// Copied from <Utility.H>
//
#define CH_IGNORE_MAX 100000

istream&
operator>> (istream&   is,
            IndexType& it)
{
  char D_DECL6(t0,t1,t2,
               t3,t4,t5);

  D_EXPR6( is.ignore(CH_IGNORE_MAX, '(') >> t0,
           is.ignore(CH_IGNORE_MAX, ',') >> t1,
           is.ignore(CH_IGNORE_MAX, ',') >> t2,
           is.ignore(CH_IGNORE_MAX, ',') >> t3,
           is.ignore(CH_IGNORE_MAX, ',') >> t4,
           is.ignore(CH_IGNORE_MAX, ',') >> t5);
  is.ignore(CH_IGNORE_MAX, ')');
  D_TERM6(
          CH_assert(t0 == 'C' || t0 == 'N'); t0=='N'?it.set(0):it.unset(0); ,
          CH_assert(t1 == 'C' || t1 == 'N'); t1=='N'?it.set(1):it.unset(1); ,
          CH_assert(t2 == 'C' || t2 == 'N'); t2=='N'?it.set(2):it.unset(2); ,
          CH_assert(t3 == 'C' || t3 == 'N'); t3=='N'?it.set(3):it.unset(3); ,
          CH_assert(t4 == 'C' || t4 == 'N'); t4=='N'?it.set(4):it.unset(4); ,
          CH_assert(t5 == 'C' || t5 == 'N'); t5=='N'?it.set(5):it.unset(5));

  if (is.fail())
        MayDay::Error("operator>>(ostream&,IndexType&) failed");

    return is;
}

// const Box&
// Box::TheUnitBox ()
// {
//     static const Box Unit(IntVect::Zero, IntVect::Unit);
//     return Unit;
// }

//
// Administrative functions.
//
Box::Box ()
    : smallend(IntVect::Unit),
      bigend(IntVect::Zero),
      //      len(IntVect::Zero),
      btype()
{
}

Box::Box (const IntVect& small,
          const int*     vec_len)
    : smallend(small),
      bigend(D_DECL6(small[0]+vec_len[0]-1,
                     small[1]+vec_len[1]-1,
                     small[2]+vec_len[2]-1,
                     small[3]+vec_len[3]-1,
                     small[4]+vec_len[4]-1,
                     small[5]+vec_len[5]-1))
{
  CH_assert(D_TERM6(vec_len[0] >= 0, && vec_len[1] >= 0, && vec_len[2] >= 0,
                    && vec_len[3] >= 0, && vec_len[4] >= 0, && vec_len[5] >= 0));

  //  D_EXPR(len[0] = vec_len[0], len[1] = vec_len[1], len[2] = vec_len[2]);
}

Box::Box (const IntVect&   small,
          const IntVect&   big,
          const IndexType& t)
    : smallend(small),
      bigend(big),
      btype(t)
{
  CH_assert (small <= big);

  computeBoxLen();
}

void Box::define(const IntVect&   small,
                 const IntVect&   big,
                 const IndexType& t)
{
  CH_assert (small <= big);
  smallend = small;
  bigend   = big;
  btype    = t;

  computeBoxLen();
}

void Box::define (const IntVect& small,
                  const IntVect& big)
{
  smallend = small;
  bigend=big;
  CH_assert (small <= big);
  computeBoxLen();
}

Box::Box (const IntVect& small,
          const IntVect& big,
          const IntVect& typ)
    : smallend(small),
      bigend(big),
      btype(typ)
{
  CH_assert (small <= big);
  CH_assert(typ >= IntVect::Zero && typ <= IntVect::Unit);

  computeBoxLen();
}

void Box::define(const IntVect& small,
                 const IntVect& big,
                 const IntVect& typ)
{
  CH_assert (small <= big);
  CH_assert(typ >= IntVect::Zero && typ <= IntVect::Unit);
  smallend = small;
  bigend   = big;
  btype    = typ;
  computeBoxLen();
}

// note that if b is undefined, then this will be undefined as well
// (but is not an error)
void Box::define(const Box& b)
{
  smallend = b.smallend;
  bigend = b.bigend;
  //  len = b.len;
  btype = b.btype;
}

/*
#if CH_SPACEDIM >= 1
    long M = (bigend[0]-smallend[0]+1);
    N = M;
#endif
#if CH_SPACEDIM >= 2
    M = (bigend[1]-smallend[1]+1);
    if (M <= LONG_MAX/N)
    {
      N *= M;
    }
    else
    {
      N = 0;
      return false;
    }
#endif
#if CH_SPACEDIM >= 3
    M = (bigend[2]-smallend[2]+1);
    if (M <= LONG_MAX/N)
    {
      N *= M;
    }
    else
    {
      N = 0;
      return false;
    }
#endif
#if CH_SPACEDIM >= 4
    M = (bigend[3]-smallend[3]+1);
    if (M <= LONG_MAX/N)
    {
      N *= M;
    }
    else
    {
      N = 0;
      return false;
    }
#endif
#if CH_SPACEDIM >= 5
    M = (bigend[4]-smallend[4]+1);
    if (M <= LONG_MAX/N)
    {
      N *= M;
    }
    else
    {
      N = 0;
      return false;
    }
#endif
#if CH_SPACEDIM >= 6
    M = (bigend[5]-smallend[5]+1);
    if (M <= LONG_MAX/N)
    {
      N *= M;
    }
    else
    {
      N = 0;
      return false;
    }
#endif
#if CH_SPACEDIM < 1 || CH_SPACEDIM > 6
    MayDay::Error("Box::numPtsOK() not implmented");
#endif

    return true;
}

long
Box::numPts () const
{
  if (bigend[0] < smallend[0]) return 0;

  long result;
  if (!numPtsOK(result))
    {
      cerr << "Bad Box::numPts:  box = " << *this << endl;
      MayDay::Error("Arithmetic overflow in Box::numPts()");
    }
  return result;
}

bool
Box::volumeOK (long& N) const
{
  return(numPtsOK(N));
}

long
Box::volume () const
{
  return(numPts());
}
*/

//
// Next:step through the rectangle with unit increment.
//

void
Box::next (IntVect& p) const
{
    CH_assert(contains(p));

    p.shift(0,1);
#if CH_SPACEDIM == 1
    // nothing more to do
#elif CH_SPACEDIM==2
    if (!(p <= bigend))
    {
        p.setVal(0,smallend[0]);
        p.shift(1,1);
    }
#elif CH_SPACEDIM==3
    if (!(p <= bigend))
    {
        p.setVal(0,smallend[0]);
        p.shift(1,1);
        if (!(p <= bigend))
        {
            p.setVal(1,smallend[1]);
            p.shift(2,1);
        }
    }
#elif CH_SPACEDIM == 4
    if (!(p <= bigend))
    {
        p.setVal(0,smallend[0]);
        p.shift(1,1);
        if (!(p <= bigend))
        {
            p.setVal(1,smallend[1]);
            p.shift(2,1);
            if (!(p <= bigend))
              {
                p.setVal(2,smallend[2]);
                p.shift(3,1);
              }
        }
    }
#elif CH_SPACEDIM == 5
    if (!(p <= bigend))
    {
        p.setVal(0,smallend[0]);
        p.shift(1,1);
        if (!(p <= bigend))
        {
            p.setVal(1,smallend[1]);
            p.shift(2,1);
            if (!(p <= bigend))
              {
                p.setVal(2,smallend[2]);
                p.shift(3,1);
                if (!(p <= bigend))
                  {
                    p.setVal(3,smallend[3]);
                    p.shift(4,1);
                  }
              }
        }
    }
#elif CH_SPACEDIM == 6
    if (!(p <= bigend))
    {
        p.setVal(0,smallend[0]);
        p.shift(1,1);
        if (!(p <= bigend))
        {
            p.setVal(1,smallend[1]);
            p.shift(2,1);
            if (!(p <= bigend))
              {
                p.setVal(2,smallend[2]);
                p.shift(3,1);
                if (!(p <= bigend))
                  {
                    p.setVal(3,smallend[3]);
                    p.shift(4,1);
                    if (!(p <= bigend))
                      {
                        p.setVal(4,smallend[4]);
                        p.shift(5,1);
                      }
                  }
              }
        }
    }
#else
    bad spacedim!
#endif
}

//
// Scan point over region of object Box with a vector incrment
// point incrments by 0 direction portion of vector.  When end of
// Box is reached, an increment is made with the 1 direction portion
// of the vector, and the 0 direction scan is resumed.
// effectively, we are scanning a Box, whose length vector is the argument
// vector over the object Box.
// when scan is over, the argument point is over edge of object Box
// this is the signal that we can go no further.
//

void
Box::next (IntVect&   p,
           const int* shv) const
{
    CH_assert(contains(p));

#if   CH_SPACEDIM==1
    p.shift(0,shv[0]);
#elif CH_SPACEDIM==2
    p.shift(0,shv[0]);
    if (!(p <= bigend))
    {
        //
        // Reset 1 coord is on edge, and 2 coord is incremented.
        //
        p.setVal(0,smallend[0]);
        p.shift(1,shv[1]);
    }
#elif CH_SPACEDIM==3
    p.shift(0,shv[0]);
    if (!(p <= bigend))
    {
        //
        // Reset 1 coord is on edge, and 2 coord is incremented.
        //
        p.setVal(0,smallend[0]);
        p.shift(1,shv[1]);
        if (!(p <= bigend))
        {
            p.setVal(1,smallend[1]);
            p.shift(2,shv[2]);
        }
    }
#elif CH_SPACEDIM==4
    p.shift(0,shv[0]);
    if (!(p <= bigend))
    {
        //
        // Reset 1 coord is on edge, and 2 coord is incremented.
        //
        p.setVal(0,smallend[0]);
        p.shift(1,shv[1]);
        if (!(p <= bigend))
        {
            p.setVal(1,smallend[1]);
            p.shift(2,shv[2]);
            if (!(p <= bigend))
              {
                p.setVal(2,smallend[2]);
                p.shift(3,shv[3]);
              }
        }
    }
#elif CH_SPACEDIM==5
    p.shift(0,shv[0]);
    if (!(p <= bigend))
    {
        //
        // Reset 1 coord is on edge, and 2 coord is incremented.
        //
        p.setVal(0,smallend[0]);
        p.shift(1,shv[1]);
        if (!(p <= bigend))
        {
            p.setVal(1,smallend[1]);
            p.shift(2,shv[2]);
            if (!(p <= bigend))
              {
                p.setVal(2,smallend[2]);
                p.shift(3,shv[3]);
                if (!(p <= bigend))
                  {
                    p.setVal(3,smallend[3]);
                    p.shift(4,shv[4]);
                  }
              }
        }
    }
#elif CH_SPACEDIM==6
    p.shift(0,shv[0]);
    if (!(p <= bigend))
    {
        //
        // Reset 1 coord is on edge, and 2 coord is incremented.
        //
        p.setVal(0,smallend[0]);
        p.shift(1,shv[1]);
        if (!(p <= bigend))
        {
            p.setVal(1,smallend[1]);
            p.shift(2,shv[2]);
            if (!(p <= bigend))
              {
                p.setVal(2,smallend[2]);
                p.shift(3,shv[3]);
                if (!(p <= bigend))
                  {
                    p.setVal(3,smallend[3]);
                    p.shift(4,shv[4]);
                    if (!(p <= bigend))
                      {
                        p.setVal(4,smallend[4]);
                        p.shift(5,shv[5]);
                      }
                  }
              }
        }
    }
#else
    bad_spacedim!
#endif
}

//
// Modified Box is low end, returned Box is high end.
// If CELL: chop_pnt included in high end.
// If NODE: chop_pnt included in both Boxes.
//

Box
Box::chop (int dir,
           int chop_pnt)
{
  CH_assert(!isEmpty());
  //
  // Define new high end Box including chop_pnt.
  //
  IntVect sm(smallend);
  IntVect bg(bigend);
  sm.setVal(dir,chop_pnt);
  if (btype[dir])
  {
    //
    // NODE centered Box.
    //
    CH_assert(chop_pnt > smallend[dir] && chop_pnt < bigend[dir]);
    //
    // Shrink original Box to just contain chop_pnt.
    //
    bigend.setVal(dir,chop_pnt);
  }
  else
  {
    //
    // CELL centered Box.
    //
    CH_assert(chop_pnt > smallend[dir] && chop_pnt <= bigend[dir]);
    //
    // Shrink origional Box to one below chop_pnt.
    //
    bigend.setVal(dir,chop_pnt-1);
  }
  computeBoxLen();
  return Box(sm,bg,btype);
}


void
Box::degenerate( Box& a_to, const SliceSpec& a_sliceSpec,
                 bool* a_outofbounds /*=0*/ ) const
{
  CH_assert(!isEmpty());
  CH_assert( (a_sliceSpec.direction >= 0) && (a_sliceSpec.direction < CH_SPACEDIM) );
  IntVect smTo, bgTo;
  for ( int i=0;i<CH_SPACEDIM;++i )
    {
      smTo[i] = this->smallend[i];
      bgTo[i] = this->bigend[i];
    }
  smTo[a_sliceSpec.direction] = a_sliceSpec.position;
  bgTo[a_sliceSpec.direction] = a_sliceSpec.position;
  Box result( smTo, bgTo );

  if ( a_outofbounds )
    {
      if ( (a_sliceSpec.position < this->smallend[a_sliceSpec.direction])
      ||  (a_sliceSpec.position > this->bigend[a_sliceSpec.direction]) )
        {
          *a_outofbounds = true;
        }
      else
        {
          *a_outofbounds = false;
        }
    }
  a_to = result;
}

//
// I/O functions.
//

ostream&
operator<< (ostream&   os,
            const Box& b)
{
  if ( Box::s_tempestOutputFormat )
  {
    os << '['
       << b.smallend << ','
       << b.bigend
       << ']';
  } else
  {
    os << '('
       << b.smallend << ' '
       << b.bigend   << ' '
       << b.btype.ixType()
       << ')';
  }
  if (os.fail())
      MayDay::Error("operator<<(ostream&,Box&) failed");
  return os;
}

void Box::p() const
{
  pout() << *this << "\n";
}

//
// Moved out of Utility.H
//
#define CH_IGNORE_MAX 100000

istream&
operator>> (istream& is,
            Box&     b)
{
    is >> ws;
    char c;
    is >> c;
    is.putback(c);
    if (c == '(')
    {
        is.ignore(CH_IGNORE_MAX, '(');
        is >> b.smallend ;
        is >> b.bigend ;
        IntVect v;
        is >> v;
        b.btype = IndexType(v);
        CH_assert(b.btype.ok());
        is.ignore(CH_IGNORE_MAX,')');
        b.computeBoxLen();
    }
    else if (c == '<')
    {
        is >> b.smallend;
        is >> b.bigend;
        IntVect v;
        is >> v;
        b.btype = IndexType(v);
        CH_assert(b.btype.ok());
        b.computeBoxLen();
    }
    else
        MayDay::Error("operator>>(istream&,Box&): expected \'<\'");

    if (is.fail())
        MayDay::Error("operator>>(istream&,Box&) failed");

    return is;
}

void
Box::dumpOn (ostream& strm) const
{
    strm << "Box "
         << smallend
         << " to "
         << bigend
         << " type ["
         << btype.ixType()
         << "]"
         << '\n';

    if (strm.fail())
        MayDay::Error("Box::dumpOn(ostream&) failed");
}

/*static*/ void Box::setTempestOutputFormat( bool b )
{
  Box::s_tempestOutputFormat = b;
}
bool Box::s_tempestOutputFormat = false;

#include "NamespaceFooter.H"
#include "UsingNamespace.H"
#include "BaseNamespaceHeader.H"

template < >
void linearIn<Box>(Box& a_outputT, const void* const a_inBuf)
{
  const unsigned char* buf = static_cast<const unsigned char*>(a_inBuf);
  IntVect lo, hi;
  CH_XD::linearIn(lo, buf);
  buf += CH_XD::linearSize(lo);
  CH_XD::linearIn(hi, buf);
  if (lo <= hi)
    {
      a_outputT = Box(lo,hi);
    }
  else
    {
      a_outputT = Box();
    }
}

template < >
void linearOut<Box>(void* const a_outBuf, const Box& a_inputT)
{
  unsigned char* buf = static_cast<unsigned char*>(a_outBuf);
  CH_XD::linearOut(buf, a_inputT.smallEnd());
  buf += CH_XD::linearSize(a_inputT.smallEnd());
  CH_XD::linearOut(buf, a_inputT.bigEnd());
}

template < >
int linearSize<Box>(const Box& a_input)
{
  //box is stored as 2*spaceDim integers
  return(2*SpaceDim*sizeof(int));
}

//Vector<Box>  specialization
template < > int linearSize(const Vector<Box>& a_input)
{
  return linearListSize(a_input);
}
template < > void linearIn(Vector<Box>& a_outputT, const void* const inBuf)
{
  linearListIn(a_outputT, inBuf);
}
template < > void linearOut(void* const a_outBuf, const Vector<Box>& a_inputT)
{
  linearListOut(a_outBuf, a_inputT);
}

//Vector<Vector<Box> >  specialization
template < > int linearSize(const Vector<Vector<Box> >& a_input)
{
  return linearListSize(a_input);
}
template < > void linearIn(Vector<Vector<Box> >& a_outputT, const void* const inBuf)
{
  linearListIn(a_outputT, inBuf);
}
template < > void linearOut(void* const a_outBuf, const Vector<Vector<Box> >& a_inputT)
{
  linearListOut(a_outBuf, a_inputT);
}



#include "BaseNamespaceFooter.H"
