#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _LOHICENTERAPPLYOP_H_
#define _LOHICENTERAPPLYOP_H_

#include "Box.H"
#include "LoHiCenter.H"
#include "BaseFab.H"
#include "BaseFabMacros.H"
#include "ProblemDomain.H"
#include "NamespaceHeader.H"

template<typename T,
         typename PrePointOp,
         typename CenterStencilOp,
         typename OffsetStencilOp,
         typename PostPointOp>
inline void
loHiCenter3ApplyOp(const Box&             a_box,
                   const ProblemDomain&   a_domain,
                   PrePointOp&      a_prePointOp,
                   CenterStencilOp& a_centerStencilOp,
                   OffsetStencilOp& a_offsetStencilOp,
                   PostPointOp&     a_postPointOp)
{
  CH_assert(a_domain.contains(a_box));
  
  const Box box1 = grow(a_box, 1) & a_domain;  // Touched cells in a_oprFab
  const Box cenBox = grow(box1, -1);           // Cells in result updated by
                                               // a centered stencil

  if (cenBox.contains(a_box))
    // Can use an optimized operation since all cells are centered
    {
      MD_BOXLOOP(a_box, i)
        {
          a_prePointOp(MD_EXPANDIX(i));
          T stencilTerm = MD_DIRSUM(
            [=, &centerStencilOp = a_centerStencilOp]
            (const int            a_dir,
             MD_DECLIX(const int, a_ii))
            {
              return centerStencilOp(MD_EXPANDIX(i),
                                     a_dir,
                                     MD_EXPANDIX(a_ii));
            });
          a_postPointOp(MD_EXPANDIX(i), stencilTerm);
        }
    }
  else
    // Have one-sided operations to perform
    {
      BaseFab<T> stencilTerm(AllocOn::cpu, a_box, 1, nullptr, 0.);
      // Do stencil ops in all cells direction by direction
      for (int dir = 0; dir != SpaceDim; ++dir)
        {
          Box loBox, hiBox, centerBox, entireBox;
          int hasLo, hasHi;
          // Generate the domain boundary boxes, loBox and hiBox, if there
          // are domain boundaries here
          Box box1dir(a_box);
          box1dir.grow(dir, 1);
          loHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox, box1dir,
                     a_domain, dir);
          const int MD_ID(ii, dir);
          if (hasLo)
            {
              MD_BOXLOOP(loBox, i)
                {
                  a_prePointOp(MD_EXPANDIX(i));
                  stencilTerm[MD_IX(i, 0)] +=
                    a_offsetStencilOp(-1, MD_EXPANDIX(i), dir, MD_EXPANDIX(ii));
                }
            }
          MD_BOXLOOP(centerBox, i)
            {
              a_prePointOp(MD_EXPANDIX(i));
              stencilTerm[MD_IX(i, 0)] +=
                a_centerStencilOp(MD_EXPANDIX(i), dir, MD_EXPANDIX(ii));
            }
          if (hasHi)
            {
              MD_BOXLOOP(hiBox, i)
                {
                  a_prePointOp(MD_EXPANDIX(i));
                  stencilTerm[MD_IX(i, 0)] +=
                    a_offsetStencilOp(1, MD_EXPANDIX(i), dir, MD_EXPANDIX(ii));
                }
            }
        }
      // Do point ops in all cells
      MD_BOXLOOP(a_box, i)
        {
          a_prePointOp(MD_EXPANDIX(i));
          a_postPointOp(MD_EXPANDIX(i), stencilTerm[MD_IX(i, 0)]);
        }
    }
}

#include "NamespaceFooter.H"
#endif  /* _LOHICENTERAPPLYOP_H_ */
