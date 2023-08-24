#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _REMOTECOPIER_H_
#define _REMOTECOPIER_H_

#include "Copier.H"
#include <map>
#include "NamespaceHeader.H"

/*============================================================================*/
/// Record of remote motions
/** The data indices must be filled before this can be used.  Both source and
 *  destination indices are required in all cases.
 */
/*============================================================================*/

struct RemoteMotion
{
  using receive_only = std::true_type;
  using send_only = std::false_type;

  /// Will produce a send and receive message if all boxes are populated
  RemoteMotion(const Box& a_lclSto, const Box& a_lclSnd, const Box& a_lclRcv,
               const Box& a_rmtSto, const Box& a_rmtSnd, const Box& a_rmtRcv,
               const int a_rmtProc)
    :
    m_lclSto(a_lclSto),
    m_lclSnd(a_lclSnd),
    m_lclRcv(a_lclRcv),
    m_rmtSto(a_rmtSto),
    m_rmtSnd(a_rmtSnd),
    m_rmtRcv(a_rmtRcv),
    m_rmtProc(a_rmtProc),
    m_lclStoDidx(),
    m_rmtStoDidx()
    {
#ifndef NDEBUG
      if (!a_lclSto.isEmpty())
        {
          CH_assert(a_lclSto.contains(a_lclSnd));
          CH_assert(a_lclSto.contains(a_lclRcv));
        }
      if (!a_rmtSto.isEmpty())
        {
          CH_assert(a_rmtSto.contains(a_rmtSnd));
          CH_assert(a_rmtSto.contains(a_rmtRcv));
        }
      CH_assert(a_lclSnd.numPts() == a_rmtRcv.numPts());
      CH_assert(a_lclRcv.numPts() == a_rmtSnd.numPts());
#endif
    }
  /// Will produce a receive message
  /** Call with first argument as RemoteMotion::receive_only{} to distinguish
   *  from following constructor
   */
  RemoteMotion(receive_only, const Box& a_lclSto, const Box& a_lclRcv,
               const Box& a_rmtSto, const Box& a_rmtSnd,
               const int a_rmtProc)
    :
    m_lclSto(a_lclSto),
    m_lclSnd(),
    m_lclRcv(a_lclRcv),
    m_rmtSto(a_rmtSto),
    m_rmtSnd(a_rmtSnd),
    m_rmtRcv(),
    m_rmtProc(a_rmtProc),
    m_lclStoDidx(),
    m_rmtStoDidx()
    {
#ifndef NDEBUG
      if (!a_lclSto.isEmpty())
        {
          CH_assert(a_lclSto.contains(a_lclRcv));
        }
      if (!a_rmtSto.isEmpty())
        {
          CH_assert(a_rmtSto.contains(a_rmtSnd));
        }
      CH_assert(a_lclRcv.numPts() == a_rmtSnd.numPts());
#endif
    }
#ifdef CH_MPI
  /// Will produce a send message (interface only available with MPI)
  /** Call with first argument as RemoteMotion::send_only{} to distinguish
   *  from previous constructor
   */
  RemoteMotion(send_only, const Box& a_lclSto, const Box& a_lclSnd,
               const Box& a_rmtSto, const Box& a_rmtRcv,
               const int a_rmtProc)
    :
    m_lclSto(a_lclSto),
    m_lclSnd(a_lclSnd),
    m_lclRcv(),
    m_rmtSto(a_rmtSto),
    m_rmtSnd(),
    m_rmtRcv(a_rmtRcv),
    m_rmtProc(a_rmtProc),
    m_lclStoDidx(),
    m_rmtStoDidx()
    {
#ifndef NDEBUG
      if (!a_lclSto.isEmpty())
        {
          CH_assert(a_lclSto.contains(a_lclSnd));
        }
      if (!a_rmtSto.isEmpty())
        {
          CH_assert(a_rmtSto.contains(a_rmtRcv));
        }
      CH_assert(a_lclSnd.numPts() == a_rmtRcv.numPts());
#endif
    }
#endif
  // Note: m_lclSto and m_lclRcv are not used to construct the motion items in
  // the copier.  However, they are often needed to find the data indices.  If
  // superfluous, just set as the minbox containing send and receive boxes.
  Box m_lclSto;                       ///< Storage defined for m_lclStoDidx
  Box m_lclSnd;                       ///< Region of lclSto sent
  Box m_lclRcv;                       ///< Region of lclSto received
  Box m_rmtSto;                       ///< Storage defined for m_rmtStoDidx
  Box m_rmtSnd;                       ///< Region of rmtSto sent
  Box m_rmtRcv;                       ///< Region of rmtSto received
  int m_rmtProc;                      ///< ID of remote process
  // Retrieve both of these from a layout index
  DataIndex m_lclStoDidx;             ///< DataIndex of local storage
  DataIndex m_rmtStoDidx;             ///< DataIndex of remote storage
};

/*============================================================================*/

class RemoteCopier: public Copier
{
  public:

  RemoteCopier();
  RemoteCopier(const DisjointBoxLayout& a_level, const BoxLayout& a_dest, const IntVect& a_ghostVect = IntVect::Zero, bool a_exchange = false);

  virtual ~RemoteCopier();

  void define(const DisjointBoxLayout& a_level,
              const BoxLayout& a_dest,
              const IntVect& a_destGhost,
              bool  a_exchange = false);

  void define(const DisjointBoxLayout& a_level,
              const BoxLayout& a_dest,
              bool a_exchange = false)
  {
    define(a_level, a_dest, IntVect::Zero, a_exchange);
  }

  //! This allows one to define a transfer of information from boxes that are contained by
  //! those in a given disjoint box layout.
  void define(const DisjointBoxLayout& a_from,
              const DisjointBoxLayout& a_to,
              const DisjointBoxLayout& a_sourceBoxes,
              const BoxLayout& a_destBoxes,
              const std::map<Box, Box>& sourceToDestMapping,
              const IntVect& a_destGhost,
              bool a_exchange = false);

  //! This allows one to define a transfer of information from boxes that are contained by
  //! those in a given disjoint box layout.
  void define(const DisjointBoxLayout& a_from,
              const DisjointBoxLayout& a_to,
              const DisjointBoxLayout& a_sourceBoxes,
              const BoxLayout& a_destBoxes,
              const std::map<Box, Box>& sourceToDestMapping,
              bool a_exchange = false)
  {
    define(a_from, a_to, a_sourceBoxes, a_destBoxes, sourceToDestMapping, IntVect::Zero, a_exchange);
  }

  //! Same as others but all motions are already given
  void define(const Vector<RemoteMotion>& a_rMotions);

  private:

  void define(const DisjointBoxLayout& a_level,
              const BoxLayout& a_dest,
              const ProblemDomain& a_domain,
              bool  a_exchange = false)
  {
  }

  void define(const BoxLayout& a_level,
              const BoxLayout& a_dest,
              const ProblemDomain& a_domain,
              const IntVect& a_destGhost,
              bool  a_exchange = false)
  {
  }

  void exchangeDefine(const DisjointBoxLayout& a_grids,
                      const IntVect& a_ghost)
  {
  }

  void ghostDefine(const DisjointBoxLayout& a_src,
                   const DisjointBoxLayout& a_dest,
                   const ProblemDomain& a_domain,
                   const IntVect& a_srcGhost)
  {
  }

  void reverse()
  {
  }

};

#include "NamespaceFooter.H"
#endif
