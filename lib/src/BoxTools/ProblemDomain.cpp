#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "ProblemDomain.H"
#include "MayDay.H"
#include <list>
#include "NamespaceHeader.H"

using std::istream;
using std::ws;

ShiftIterator::ShiftIterator(const bool* a_isPeriodic)
{
  m_index = 0;
  computeShifts(a_isPeriodic);
}

ShiftIterator::ShiftIterator(const bool* a_isPeriodic, const IntVect& a_numWraps)
{
  m_index = 0;
  computeShifts(a_isPeriodic, a_numWraps);
}

ShiftIterator::~ShiftIterator()
{
  m_shift_vectors.clear();
}

void
ShiftIterator::computeShifts(const bool* a_isPeriodic)
{
  IntVect numWraps = IntVect::Unit;
  computeShifts(a_isPeriodic, numWraps);
}


void
ShiftIterator::computeShifts(const bool* a_isPeriodic, const IntVect& a_numWraps)
{  
  m_shift_vectors.resize(0);

  // handle no directions are periodic
  if (D_TERM6(!a_isPeriodic[0],
              && !a_isPeriodic[1],
              && !a_isPeriodic[2],
              && !a_isPeriodic[3],
              && !a_isPeriodic[4],
              && !a_isPeriodic[5]))
    {
      // return with no periodic directions defined
    }
  // now handle if all directions are periodic
  else if (D_TERM6(a_isPeriodic[0],
                   && a_isPeriodic[1],
                   && a_isPeriodic[2],
                   && a_isPeriodic[3],
                   && a_isPeriodic[4],
                   && a_isPeriodic[5]))
    {
      // loop over all possible shifts in this case (including corners)
      D_TERM6(
              for (int xVect=-a_numWraps[0]; xVect<=a_numWraps[0]; xVect++),
              for (int yVect=-a_numWraps[1]; yVect<=a_numWraps[1]; yVect++),
              for (int zVect=-a_numWraps[2]; zVect<=a_numWraps[2]; zVect++),
              for (int uVect=-a_numWraps[3]; uVect<=a_numWraps[3]; uVect++),
              for (int vVect=-a_numWraps[4]; vVect<=a_numWraps[4]; vVect++),
              for (int wVect=-a_numWraps[5]; wVect<=a_numWraps[5]; wVect++) )
        {
          IntVect shiftVect(D_DECL6(xVect,yVect,zVect,
                                    uVect,vVect,wVect));
          // no point in adding the zero vector
          if (shiftVect != IntVect::Zero)
            {
              m_shift_vectors.push_back(shiftVect);
            }
        } // end loop over possible shift vectors
    }
  // case where some but not all directions are periodic
  // this could probably also swallow the completely periodic case
  else
    {
      // i think this should work
      // this is designed to filter out non-periodic directions
      D_TERM6(int xMult = (a_isPeriodic[0]?a_numWraps[0]:0);,
              int yMult = (a_isPeriodic[1]?a_numWraps[1]:0);,
              int zMult = (a_isPeriodic[2]?a_numWraps[2]:0);,
              int uMult = (a_isPeriodic[3]?a_numWraps[3]:0);,
              int vMult = (a_isPeriodic[4]?a_numWraps[4]:0);,
              int wMult = (a_isPeriodic[5]?a_numWraps[5]:0););

      // in non-periodic dirctions, shift vector component should
      // always be 0
      D_TERM6(
              for (int xVect=-1*xMult; xVect<1+xMult; xVect++),
              for (int yVect=-1*yMult; yVect<1+yMult; yVect++),
              for (int zVect=-1*zMult; zVect<1+zMult; zVect++),
              for (int uVect=-1*uMult; uVect<1+uMult; uVect++),
              for (int vVect=-1*vMult; vVect<1+vMult; vVect++),
              for (int wVect=-1*wMult; wVect<1+wMult; wVect++)  )
        {
          IntVect shiftVect(D_DECL6(xVect, yVect, zVect,
                                    uVect, vVect, wVect));
          // no point in adding the zero vector
          if (shiftVect != IntVect::Zero)
            {
              m_shift_vectors.push_back(shiftVect);
            }
        } // end loop over possible shift vectors
    } // end case where some directions are periodic
  // set index to 0...
  begin();
}

ProblemDomain::ProblemDomain(const Box& a_domBox)
  : m_isConnected(false)
{
  define(a_domBox);
}

ProblemDomain::ProblemDomain(const Box& a_domBox,
                             const bool* a_isPeriodic)
  : m_isConnected(false)
{
  define(a_domBox, a_isPeriodic);
}

ProblemDomain::ProblemDomain(const Box&           a_domBox,
                             const ProblemDomain& a_src)
{
  define(a_domBox, a_src.m_isPeriodic);
  m_isConnected = a_src.m_isConnected;
}

ProblemDomain::ProblemDomain(const IntVect& small, const IntVect& big)
  : m_isConnected(false)
{
  define(small, big);
}

ProblemDomain::ProblemDomain(const IntVect& small, const IntVect& big,
                             const bool* a_isPeriodic)
  : m_isConnected(false)
{
  define(small, big, a_isPeriodic);
}

ProblemDomain::ProblemDomain(const IntVect& small, const int* vec_len)
  : m_isConnected(false)
{
  define(small, vec_len);
}

ProblemDomain::ProblemDomain(const IntVect& small, const int* vec_len,
                             const bool* a_isPeriodic)
  : m_isConnected(false)
{
  define(small, vec_len, a_isPeriodic);
}

void
ProblemDomain::define(const Box& a_domBox)
{
  m_isConnected = false;
  m_domainBox = a_domBox;
  // default is a non-periodic domain
  for (int dir=0; dir<SpaceDim; dir++)
  {
    m_isPeriodic[dir] = false;
  }
  m_shiftIt.computeShifts(m_isPeriodic);
}

void
ProblemDomain::define(const Box& a_domBox,
                      const bool* a_isPeriodic)
{
  m_isConnected = false;
  m_domainBox = a_domBox;
  for (int dir=0; dir<SpaceDim; dir++)
  {
    m_isPeriodic[dir] = a_isPeriodic[dir];
  }
  m_shiftIt.computeShifts(m_isPeriodic);
}

void
ProblemDomain::define(const IntVect& small, const IntVect& big)
{
  m_isConnected = false;
  m_domainBox.define(small, big);
  // default is a non-periodic domain
  for (int dir=0; dir<SpaceDim; dir++)
  {
    m_isPeriodic[dir] = false;
  }
  m_shiftIt.computeShifts(m_isPeriodic);
}

void
ProblemDomain::define(const IntVect& small, const IntVect& big,
                      const bool* a_isPeriodic)
{
  m_isConnected = false;
  m_domainBox.define(small, big);
  // default is a non-periodic domain
  for (int dir=0; dir<SpaceDim; dir++)
  {
    m_isPeriodic[dir] = a_isPeriodic[dir];
  }
  m_shiftIt.computeShifts(m_isPeriodic);
}

void
ProblemDomain::define(const IntVect& small, const int* vec_len)
{
  m_isConnected = false;
  m_domainBox = Box(small, vec_len);
  // default is a non-periodic domain
  for (int dir=0; dir<SpaceDim; dir++)
  {
    m_isPeriodic[dir] = false;
  }
  m_shiftIt.computeShifts(m_isPeriodic);
}

void
ProblemDomain::define(const IntVect& small, const int* vec_len,
                      const bool* a_isPeriodic)
{
  m_isConnected = false;
  m_domainBox= Box(small, vec_len);
  // default is a non-periodic domain
  for (int dir=0; dir<SpaceDim; dir++)
  {
    m_isPeriodic[dir] = a_isPeriodic[dir];
  }
  m_shiftIt.computeShifts(m_isPeriodic);
}

// returns true if box intersects problem domain
bool
ProblemDomain::intersects(const Box& a_box) const
{
  if (isEmpty() || a_box.isEmpty()) return false;

  CH_assert (m_domainBox.sameType(a_box));
  IntVect lo(m_domainBox.smallEnd());
  IntVect hi(m_domainBox.bigEnd());
  lo.max(a_box.smallEnd());
  hi.min(a_box.bigEnd());

  return stc::andEachElement<SpaceDim>(
    [&]
    (const stc::array_size_type a_dir)
    {
      if (isConnected(a_dir, Side::Lo)) lo[a_dir] = a_box.smallEnd(a_dir);
      if (isConnected(a_dir, Side::Hi)) hi[a_dir] = a_box.bigEnd  (a_dir);
      return m_isPeriodic[a_dir] || (lo[a_dir] <= hi[a_dir]);
    });
}

bool
ProblemDomain::intersectsNotEmpty(const Box& a_box) const
{
  CH_assert (m_domainBox.sameType(a_box));
  IntVect lo(m_domainBox.smallEnd());
  IntVect hi(m_domainBox.bigEnd());
  lo.max(a_box.smallEnd());
  hi.min(a_box.bigEnd());

  return stc::andEachElement<SpaceDim>(
    [&]
    (const stc::array_size_type a_dir)
    {
      if (isConnected(a_dir, Side::Lo)) lo[a_dir] = a_box.smallEnd(a_dir);
      if (isConnected(a_dir, Side::Hi)) hi[a_dir] = a_box.bigEnd  (a_dir);
      return m_isPeriodic[a_dir] || (lo[a_dir] <= hi[a_dir]);
    });
}

bool ProblemDomain::periodicAdjacent(const Box& a_box) const
{
  for (int i=0; i<CH_SPACEDIM; i++)
    {
      if (m_isPeriodic[i] && (a_box.smallEnd(i)== m_domainBox.smallEnd(i) ||
                             a_box.bigEnd(i)  == m_domainBox.bigEnd(i)))
        {
          return true;
        }
    }
  return false;
}

void ProblemDomain::insertImages(std::list<Box>& a_list, const Box& a_box) const
{
  CH_assert(m_domainBox.contains(a_box));
  if (!isPeriodic()) return;

  IntVect shiftMult(m_domainBox.size());
  ShiftIterator shiftIt = m_shiftIt;
  for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
    {
      Box shiftBox(a_box);
      IntVect shiftVect = shiftMult*shiftIt();
      shiftBox.shift(shiftVect);
      a_list.push_back(shiftBox);
    }

}

bool
ProblemDomain::intersects(const Box& box1, const Box& box2) const
{

  bool boxesIntersect = false;
  // first do the obvious and easy check
  boxesIntersect = box1.intersects(box2);

  // now do the harder check whether periodic images intersect
  if (!boxesIntersect && isPeriodic())
    {
      if (m_domainBox.contains(box1) && m_domainBox.contains(box2))
        {
          // no need to do any periodic checking
        }
      else
        {
          // loop over periodic shift directions, shift box1, and
          // check for intersections
          IntVect shiftMult(m_domainBox.size());
          ShiftIterator shiftIt = m_shiftIt;
          for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
            {
              Box shiftBox(box1);
              IntVect shiftVect = shiftMult*shiftIt();
              shiftBox.shift(shiftVect);
              if (shiftBox.intersects(box2)) boxesIntersect = true;
            }
        }
    } // end if we need to do periodic check
  return boxesIntersect;
}

void
ProblemDomain::setPeriodic(int a_dir, bool a_isPeriodic)
{
  CH_assert (a_dir < SpaceDim);
  m_isPeriodic[a_dir] = a_isPeriodic;
  m_shiftIt.computeShifts(m_isPeriodic);
}

Box
bdryLo(const ProblemDomain& a_pd, int a_dir, int a_len)
{
  Box bndryBox;
  if (!a_pd.isPeriodic(a_dir) && !a_pd.isConnected(a_dir, Side::Lo))
    {
      bndryBox = bdryLo(a_pd.domainBox(),a_dir,a_len);
    }
  return bndryBox;
}

Box
bdryHi(const ProblemDomain& a_pd, int a_dir, int a_len)
{
  Box bndryBox;
  if (!a_pd.isPeriodic(a_dir) && !a_pd.isConnected(a_dir, Side::Hi))
    {
      bndryBox = bdryHi(a_pd.domainBox(),a_dir,a_len);
    }
  return bndryBox;
}

Box
adjCellLo(const ProblemDomain& a_pd, int a_dir, int a_len)
{
  Box adjBox;
  if (!a_pd.isPeriodic(a_dir) && !a_pd.isConnected(a_dir, Side::Lo))
    {
      adjBox = adjCellLo(a_pd.domainBox(),a_dir, a_len);
    }
  return adjBox;

}

Box
adjCellHi(const ProblemDomain& a_pd, int a_dir, int a_len)
{
  Box adjBox;
  if (!a_pd.isPeriodic(a_dir) && !a_pd.isConnected(a_dir, Side::Hi))
    {
      adjBox = adjCellHi(a_pd.domainBox(),a_dir, a_len);
    }
  return adjBox;
}

Box
ProblemDomain::operator&(const Box& b) const
{
  Box intersectBox;
  // check for empty -- if either box is empty, return empty intersectBox
  if (!b.isEmpty() && !isEmpty())
    {
      intersectBox = m_domainBox;

      if (b.type() != intersectBox.type())
      {
        // Check b's centering and adjust intersectBox as needed
        for (int dir = 0; dir < SpaceDim; dir++)
        {
          if (b.type(dir) != intersectBox.type(dir))
          {
            if (b.type(dir) == IndexType::NODE)
            {
              intersectBox.surroundingNodes(dir);
            }
            else
            {
              intersectBox.enclosedCells(dir);
            }
          }
        }
      }

      // do this by creating an intersection box which will
      // take the periodic case into account
      if (isPeriodic())
        {
          for (int dir=0; dir<SpaceDim; dir++)
            {
              if (m_isPeriodic[dir])
                {
                  // the idea behind this is that intersect does nothing in
                  // periodic directions
                  intersectBox.setRange(dir, b.smallEnd(dir), b.size(dir));
                }
            }
        }

      // consider connections which are defined per side.  the min and max
      // handle situations where 'b' is completely outside the domain.
      if (isConnected())
        {
          for (int dir = 0; dir != SpaceDim; ++dir)
            {
              if (isConnected(dir, Side::Lo))
                {
                  intersectBox.setSmall(
                    dir, std::min(intersectBox.bigEnd(dir), b.smallEnd(dir)));
                }
              if (isConnected(dir, Side::Hi))
                {
                  intersectBox.setBig(
                    dir, std::max(intersectBox.smallEnd(dir), b.bigEnd(dir)));
                }
            }
        }

      intersectBox &= b;
    } // end if neither box is empty

  return intersectBox;
}

ProblemDomain&
ProblemDomain::refine(int a_refinement_ratio)
{
  m_domainBox.refine(a_refinement_ratio);
  return *this;
}

ProblemDomain
refine(const ProblemDomain& a_probdomain, int a_refinement_ratio)
{

  if (a_probdomain.isEmpty()) return a_probdomain;
  Box newdomain = a_probdomain.domainBox();
  newdomain.refine(a_refinement_ratio);

  return ProblemDomain(newdomain, a_probdomain);
}

ProblemDomain&
ProblemDomain::refine(const IntVect& a_refinement_ratio)
{
  m_domainBox.refine(a_refinement_ratio);
  return *this;
}

ProblemDomain
refine(const ProblemDomain& a_probdomain, const IntVect& a_refinement_ratio)
{
  Box newdomain(a_probdomain.domainBox());
  newdomain.refine(a_refinement_ratio);

  return ProblemDomain(newdomain, a_probdomain);
}

ProblemDomain&
ProblemDomain::coarsen(int a_refinement_ratio)
{
  m_domainBox.coarsen(a_refinement_ratio);
  return *this;
}

ProblemDomain
coarsen(const ProblemDomain& a_probdomain, int a_refinement_ratio)
{
  Box newdomain(a_probdomain.domainBox());

  newdomain.coarsen(a_refinement_ratio);
  return ProblemDomain(newdomain, a_probdomain);
}

ProblemDomain&
ProblemDomain::coarsen(const IntVect& a_refinement_ratio)
{
  m_domainBox.coarsen(a_refinement_ratio);
  return *this;
}

ProblemDomain
coarsen(const ProblemDomain& a_probdomain,
        const IntVect& a_refinement_ratio)
{

  Box newdomain(a_probdomain.domainBox());
  newdomain.coarsen(a_refinement_ratio);

  return ProblemDomain(newdomain, a_probdomain);
}

ostream& operator<< (ostream& os, const ProblemDomain& a_probdomain)
{
  os << a_probdomain.domainBox()
     << ", " << D_TERM6( a_probdomain.m_isPeriodic[0], << ' ' <<
                         a_probdomain.m_isPeriodic[1], << ' ' <<
                         a_probdomain.m_isPeriodic[2], << ' ' <<
                         a_probdomain.m_isPeriodic[3], << ' ' <<
                         a_probdomain.m_isPeriodic[4], << ' ' <<
                         a_probdomain.m_isPeriodic[5])
     << ", " << D_TERM6( a_probdomain.isConnected(0, Side::Lo) << ' ' <<
                         a_probdomain.isConnected(0, Side::Hi), << ' ' <<
                         a_probdomain.isConnected(1, Side::Lo) << ' ' <<
                         a_probdomain.isConnected(1, Side::Hi), << ' ' <<
                         a_probdomain.isConnected(2, Side::Lo) << ' ' <<
                         a_probdomain.isConnected(2, Side::Hi), << ' ' <<
                         a_probdomain.isConnected(3, Side::Lo) << ' ' <<
                         a_probdomain.isConnected(3, Side::Hi), << ' ' <<
                         a_probdomain.isConnected(4, Side::Lo) << ' ' <<
                         a_probdomain.isConnected(4, Side::Hi), << ' ' <<
                         a_probdomain.isConnected(5, Side::Lo) << ' ' <<
                         a_probdomain.isConnected(5, Side::Hi));

  if (os.fail())
    MayDay::Error("operator<<(ostream&,ProblemDomain&) failed");

  return os;
}

//
// Copied from <Utility.H>
//
#define CH_IGNORE_MAX 100000

istream&
operator>> (istream& is, ProblemDomain& a_probdomain)
{

  is >> a_probdomain.m_domainBox;
  is >> ws;
  char c;
  is >> c;
  is.putback(c);
  if (c == ',')
    {
      is.ignore(CH_IGNORE_MAX, ',');
      for (int dir=0; dir<SpaceDim; dir++)
        {
          is >> a_probdomain.m_isPeriodic[dir];
        }
    }
  is >> c;
  is.putback(c);
  if (c == ',')
    {
      is.ignore(CH_IGNORE_MAX, ',');
      for (int dir=0; dir<SpaceDim; dir++)
        {
          is >> a_probdomain.m_isConnected[2*dir + 0]
             >> a_probdomain.m_isConnected[2*dir + 1];
        }
    }

  else
    MayDay::Error("operator>>(istream&,ProblemDomain&): expected \',\'");

  if (is.fail())
    MayDay::Error("operator>>(istream&,ProblemDomain&) failed");

  return is;
}

void
ProblemDomain::dumpOn(ostream& strm) const
{

  strm << "domainBox: ";
  m_domainBox.dumpOn(strm);
  strm << " isPeriodic = "
       << D_TERM6(m_isPeriodic[0], << " " <<
                  m_isPeriodic[1], << " " <<
                  m_isPeriodic[2], << " " <<
                  m_isPeriodic[3], << " " <<
                  m_isPeriodic[4], << " " <<
                  m_isPeriodic[5]);
  strm << " isConnected = "
       << D_TERM6(isConnected(0, Side::Lo) << ':' <<
                  isConnected(0, Side::Hi), << ' ' <<
                  isConnected(1, Side::Lo) << ':' <<
                  isConnected(1, Side::Hi), << ' ' <<
                  isConnected(2, Side::Lo) << ':' <<
                  isConnected(2, Side::Hi), << ' ' <<
                  isConnected(3, Side::Lo) << ':' <<
                  isConnected(3, Side::Hi), << ' ' <<
                  isConnected(4, Side::Lo) << ':' <<
                  isConnected(4, Side::Hi), << ' ' <<
                  isConnected(5, Side::Lo) << ':' <<
                  isConnected(5, Side::Hi))
       << '\n';

  if (strm.fail())
    MayDay::Error("ProblemDomain::dumpOn(ostream&) failed");
}

void operator &= (Box& a_box, const ProblemDomain& a_probdomain)
{
  a_box = a_probdomain & a_box;
}

Box operator & (const Box& a_box, const ProblemDomain& a_probdomain)
{
  Box returnBox = a_probdomain & a_box;
  return returnBox;
}

void ImageIterator::define(const ProblemDomain& a_domain)
{
  m_domain = a_domain;
  IntVect size = a_domain.domainBox().size();
  ShiftIterator it(a_domain.shiftIterator());
  int index=0;
  for (it.begin(); it.ok(); ++index, ++it)
    {
      m_shifter[index] = -(it()*size);
      m_quadrant[index] = a_domain.domainBox();
      m_quadrant[index].shift(it()*size);
    }
  m_shifter[index]=IntVect::Zero;
}

void ImageIterator::operator++()
{
  ++m_counter;
  if (ok())
  {
    m_current = m_box & m_quadrant[m_counter];
    m_current.shift(m_shifter[m_counter]);
  }
}
#include "NamespaceFooter.H"
