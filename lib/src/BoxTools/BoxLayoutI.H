#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _BOXLAYOUTI_H_
#define _BOXLAYOUTI_H_

#include "NamespaceHeader.H"


/*--------------------------------------------------------------------*/
//  Create a new layout with references from an older layout
/** Most commonly, this is used to create a more compact data
 *  structure.  From the compact layout, one can find the origin info
 *  using the returned a_origin.  From the orignal layout, one can
 *  find an index in the compact layout using the returned a_didxMap.
 *
 *  \tparam     T       Information about identifying origin data.
 *                      Type key_type must be defined in T, which is
 *                      used as the key for a_didxMap.  Member
 *                      T::key_type key() must be defined to convert
 *                      T into T::key_type.  key_type can be simply T.
 *  \param[in]  a_boxes The boxes for the new layout on this
 *                      processor.  The map gives for each box an
 *                      identifier of source data (e.g., commonly
 *                      including at least a DataIndex from an
 *                      original layout).
 *  \param[out] a_origin
 *                      The origin identity for each location in the
 *                      new layout.  This is given as a LayoutData.
 *                      Set to nullptr if not required as output.
 *  \param[out] a_didxMap
 *                      Gives the DataIndex in the new compact layout
 *                      based on origin T.  This is returned as a map
 *                      since it is expected that the compact layout
 *                      is significantly smaller than the original
 *                      layout.  a_didxMap->reserve() is called.  If
 *                      desired, set a_didxMap->max_load_factor(x) in
 *                      caller.  Set to nullptr if not required as
 *                      output.
 *//*-----------------------------------------------------------------*/

template <typename T>
void
BoxLayout::define(const UMap_Box_T<T>&                          a_boxes,
                  LayoutData<T> *const                          a_origin,
                  UMap_T_DataIndex<typename T::key_type> *const a_didxMap)
{
  allGatherDefine(a_boxes,
                  []
                  (const typename UMap_Box_T<T>::value_type& a_val)
                  -> const Box&
                  {
                    return a_val.first;
                  });
  const int numLocalBox = a_boxes.size();
  // Relate the indices of this new layout to something (T) that identifies
  // the origin data, this is achieved by comparing the boxes which must be
  // unique.
  if (a_origin || a_didxMap)
    {
      if (a_origin)
        {
          a_origin->define(*this);
        }
      if (a_didxMap)
        {
          a_didxMap->clear();
          a_didxMap->reserve(numLocalBox);
        }
      // Only a forward declaration of DataIterator is available so we cannot
      // use it.  If using MPI, the DataIndex are in m_dataIndex.  Otherwise,
      // this list would be rebuilt in construction of LayoutIterator.  We can
      // just construct based on the list of Entry instead.
#ifdef CH_MPI
      for (const DataIndex& didx : *m_dataIndex)
#else
      // Curious: buildDataIndex() uses unsigned for the index, but the
      // LayoutIndex constructor takes ints!  We'll just use ints here
      for (int index = 0, index_end = m_boxes->size(); index != index_end;
           ++index)
#endif
        {
#ifndef CH_MPI
          const DataIndex didx(index, index, m_layout);
#endif
          const Box& box = this->operator[](didx);
          const auto iter = a_boxes.find(box);
          CH_assert(iter != a_boxes.end());  // The box must be found
          if (a_origin)
            {
              (*a_origin)[didx] = iter->second;
            }
          if (a_didxMap)
            {
              auto ins = a_didxMap->insert({ iter->second.key(), didx });
              // Insert must have succeeded for unique keys
              CH_assert(ins.second);
            }
        }
    }
}

/*--------------------------------------------------------------------*/
//  Create a new layout with boxes just for the local processor
/** Generates a required lambda and calls generic allGatherDefine
 *  \tparam     C       A container of boxes where
 *                      C::value_type == Box
 *  \param[in]  a_boxes The boxes local to this processor
 *//*-----------------------------------------------------------------*/

template <typename C>
void
BoxLayout::allGatherDefine(const C& a_boxes)
{
  allGatherDefine(a_boxes,
                  []
                  (const Box& a_val) -> const Box&
                  {
                    return a_val;
                  });
}

/*--------------------------------------------------------------------*/
//
/** Create a new layout with boxes just for the local processor
 *
 *  \tparam     C       A container of boxes
 *  \tparam     F       Deduced from a_vt2box which is typically a
 *                      lambda.
 *  \param[in]  a_boxes The boxes local to this processor
 *  \param[in]  a_vt2box
 *                      Lambda or function object that extracts a box
 *                      from C::value_type.  E.g., if
 *                      C = std::map<Box, int>, then a_vtbox might be      \code
  []
  (const std::map<Box, int>::value_type& a_val) -> const Box&
  {
    return a_val.first;
  });                                                                   \endcode
 *                      If C::value_type == Box, do not specify a_vt2box
 *                      and an overload (above) will generate the
 *                      lambda.
 *//*-----------------------------------------------------------------*/

template <typename C, typename F>
void
BoxLayout::allGatherDefine(const C& a_boxes, F&& a_vt2box)
{
  const int numLocalBox = a_boxes.size();
#ifdef CH_MPI
  const int sizeBytesBox = linearSize(Box{});
  Vector<int> globalBoxCount(numProc());
  MPI_Allgather(&numLocalBox, 1, MPI_INT, globalBoxCount.data(), 1, MPI_INT,
                Chombo_MPI::comm);
  // Find displacements in bytes, count number of global boxes, and recast
  // globalBoxCount into bytes
  Vector<int> displacements(numProc());
  int numGlobalBox = 0;
  for (int idxProc = 0, idxProc_end = numProc(); idxProc != idxProc_end;
       ++idxProc)
    {
      displacements[idxProc] = numGlobalBox*sizeBytesBox;
      numGlobalBox += globalBoxCount[idxProc];
      globalBoxCount[idxProc] *= sizeBytesBox;
    }
  CH_assert(static_cast<long>(numGlobalBox)*sizeBytesBox <
            std::numeric_limits<int>::max());
  Vector<char> localBoxBuffer(numLocalBox*sizeBytesBox);
  int idxBox = 0;
  for (const auto& vtBox : a_boxes)
    {
      linearOut(&localBoxBuffer[idxBox*sizeBytesBox], a_vt2box(vtBox));
      ++idxBox;
    }
  Vector<char> globalBoxBuffer(numGlobalBox*sizeBytesBox);
  // Gather boxes to all
  MPI_Allgatherv(localBoxBuffer.data(), numLocalBox*sizeBytesBox,
                 MPI_BYTE, globalBoxBuffer.data(), globalBoxCount.data(),
                 displacements.data(), MPI_BYTE, Chombo_MPI::comm);
  // Now we can finally populate m_boxes.  Boxes retain processor assignments
  // from the argument layout
  m_boxes->resize(numGlobalBox);
  int idxGlobalBox = 0;
  for (int idxProc = 0, idxProc_end = numProc(); idxProc != idxProc_end;
       ++idxProc)
    {
      const int numProcBox = globalBoxCount[idxProc]/sizeBytesBox;
      for (int idxBox = 0; idxBox != numProcBox; ++idxBox)
        {
          Entry& entry = m_boxes->operator[](idxGlobalBox);
          linearIn(entry.box,
                   &globalBoxBuffer[idxGlobalBox*sizeBytesBox]);
          entry.m_procID = idxProc;
          ++idxGlobalBox;
        }
    }
  CH_assert(idxGlobalBox == numGlobalBox);
  close();
#else
  CH_assert(numProc() == 1);
  // Store the set into a vector
  Vector<Box> boxes;
  boxes.reserve(numLocalBox);
  for (const auto& vtBox : a_boxes)
    {
      boxes.push_back(a_vt2box(vtBox));
    }
  Vector<int> procIDs;
  define(boxes, procIDs);
#endif
}

#include "NamespaceFooter.H"
#endif  /* _BOXLAYOUTI_H_ */
