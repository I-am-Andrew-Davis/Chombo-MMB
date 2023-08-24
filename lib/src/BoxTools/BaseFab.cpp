#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cstring>

#include "BaseFab.H"
#include "BoxIterator.H"
#include "NamespaceHeader.H"

Real BaseFabRealSetVal = BASEFAB_REAL_SETVAL;

//?? Why not in testing folder?
template <>
int
BaseFab<int, DefaultArrayAlloc<int, ArrayClassIndex::BaseFab> >::test()
{
  int retbox = testBoxAndComp();
  if (retbox != 0)
    {
      pout() << "testboxandcomp failed" << endl;
      return retbox;
    }

  Box b(IntVect::Zero, IntVect::Unit);
  BaseFab<int> blerg;
  blerg.define(b, 1);
  int val = 4;
  blerg.setVal(val);
  for (BoxIterator bit(b); bit.ok();++bit)
    {
      if (blerg(bit(), 0) != val)
        {
          pout() << "setval or index busted" << endl;
          return -1;
        }
    }

  Box bcop(IntVect::Unit, 2*IntVect::Unit);
  BaseFab<int> bfcopy(bcop, 1);
  bfcopy.setVal(2*val);
  bfcopy.copy(blerg, 0, 0, 1);
  Box binter =bcop;
  binter &= b;
  for (BoxIterator bit(binter); bit.ok();++bit)
    {
      if (bfcopy(bit(), 0) != val)
        {
          pout() << "copy busted" << endl;
          return -2;
        }
    }

  return 0;
}

#include "NamespaceFooter.H"
