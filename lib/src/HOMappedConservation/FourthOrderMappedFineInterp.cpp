#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


/******************************************************************************/
/**
 * \file
 *
 * \brief Non-inline definitions for classes in FourthOrderMappedFineInterp.H
 *
 *//*+*************************************************************************/

#include "FourthOrderMappedFineInterp.H"
#include "FourthOrderMappedFineInterpF_F.H"
#include "FourthOrderMappedFineInterpSup.H"
#include "RealVect.H"
#include "BaseFabMacros.H"
#include "InsertionSort.H"
#include "BoxIterator.H"
#include "ProblemDomain.H"
#include "FourthOrderUtil.H"
#include "LevelGridMetrics.H"
#include "CHMatrixOps.H"

#include "NamespaceHeader.H"

// // Blas and Lapack routines needed for Least Squares operations
// #ifdef CH_USE_DOUBLE
//   #define BLAS1_DOT FORTRAN_NAME(DDOT, ddot)
//   #define BLAS2_GEMV FORTRAN_NAME(DGEMV, dgemv)
//   #define BLAS3_GEMM FORTRAN_NAME(DGEMM, dgemm)
//   #define BLAS3_SYRK FORTRAN_NAME(DSYRK, dsyrk)
//   #define BLAS3_SYMM FORTRAN_NAME(DSYMM, dsymm)
// #else
//   #define BLAS1_DOT FORTRAN_NAME(SDOT, sdot)
//   #define BLAS1_NRM2 FORTRAN_NAME(SNRM2, snrm2)
//   #define BLAS2_GEMV FORTRAN_NAME(SGEMV, sgemv)
//   #define BLAS3_GEMM FORTRAN_NAME(SGEMM, sgemm)
//   #define BLAS3_SYRK FORTRAN_NAME(SSYRK, ssyrk)
//   #define BLAS3_SYMM FORTRAN_NAME(SSYMM, ssymm)
// #endif

// extern "C" Real BLAS1_DOT(int*, Real*, int*, Real*, int*);
// extern "C" void BLAS2_GEMV(char*, int*, int*, Real*, Real*, int*, Real*, int*,
//                            Real*, Real*, int*);
// extern "C" void BLAS3_GEMM(char*, char*, int*, int*, int*, Real*, Real*, int*,
//                            Real*, int*, Real*, Real*, int*);
// extern "C" void BLAS3_SYRK(char*, char*, int*, int*, Real*, Real*, int*, Real*,
//                            Real*, int*);
// extern "C" void BLAS3_SYMM(char*, char*, int*, int*, Real*, Real*, int*, Real*,
//                            int*, Real*, Real*, int*);

/*******************************************************************************
 *
 * Class FourthOrderMappedFineInterp: member definitions
 *
 ******************************************************************************/


/*==============================================================================
 * Constructors and destructors
 *============================================================================*/

/*--------------------------------------------------------------------*/
//   Non-defining constructor
/**  \param[in]  a_degree
 *                      Degree of the polynomial
 *   \param[in]  a_nComp
 *                      Number of components in \f$<U>\f$ and
 *                      \f$<JU>\f$
 *//*-----------------------------------------------------------------*/

FourthOrderMappedFineInterp::FourthOrderMappedFineInterp(const int a_degree,
                                                         const int a_nComp)
  :
  m_degree(a_degree),
  // This is hard-coded because the code for marking the cells in the stencil
  // is quite dependent on a value of 2.  An extent of 2 means that the
  // interior stencil uses a maximum of 2 extra cells in a direction.
  m_stencilExtent(2),
  m_nComp(a_nComp),
  m_willFillGhostsWithCrFnJU(false),
  m_CrFnNumGhost(-1*IntVect::Unit),
  m_CrFnStatus(CrFnStatusInvalid),
  m_isMultiblock(false),
  m_useClipping(false),
  m_useHOClipping(false),
  m_usePostSmoothing(false)
{
  m_defined[0] = false;
  m_defined[1] = false;
  m_defined[2] = false;
  m_defined[3] = false;
}

/*--------------------------------------------------------------------*/
//   Constructor that initializes the stencils
/**  \param[in]  a_degree
 *                      Degree of the polynomial
 *   \param[in]  a_nComp
 *                      Number of components in \f$<U>\f$ and
 *                      \f$<JU>\f$
 *   \param[in]  a_refRatio
 *                      Refinement ratio in each direction
 *   \param[in]  a_h    Vector of mesh spacing in each direction on
 *                      the coarse mesh in computation space
 *//*-----------------------------------------------------------------*/

FourthOrderMappedFineInterp::FourthOrderMappedFineInterp(
  const int       a_degree,
  const int       a_nComp,
  const IntVect&  a_refRatio,
  const RealVect& a_h)
  :
  m_degree(a_degree),
  // This is hard-coded because the code for marking the cells in the stencil
  // is quite dependent on a value of 2.  An extent of 2 means that the
  // interior stencil uses a maximum of 2 extra cells in a direction.
  m_stencilExtent(2),
  m_nComp(a_nComp),
  m_willFillGhostsWithCrFnJU(false),
  m_CrFnNumGhost(-1*IntVect::Unit),
  m_CrFnStatus(CrFnStatusInvalid),
  m_isMultiblock(false),
  m_useClipping(false),
  m_useHOClipping(false),
  m_usePostSmoothing(false)
{
  m_defined[0] = false;
  m_defined[1] = false;
  m_defined[2] = false;
  m_defined[3] = false;
  defineStencils(a_refRatio, a_h);
}

/*--------------------------------------------------------------------*/
//   Destructor
/*--------------------------------------------------------------------*/

FourthOrderMappedFineInterp::~FourthOrderMappedFineInterp()
{
  if (isStencilDefined())
    {
      delete[] m_fnXiLSx;
      delete[] m_fnGradXiLSx;
    }
}


/*==============================================================================
 * Member functions
 *============================================================================*/

/*--------------------------------------------------------------------*/
//   Define that initializes or redefines the stencils
/**  Stencils are selected.  The number is minimized by considering
 *   symmetries.  Then, A^T is computed for each stencil.  Also, the
 *   powers of displacments and gradients are computed for all the
 *   fine cells in a coarse cell.
 *   \param[in]  a_refRatio
 *                      Refinement ratio in each direction
 *   \param[in]  a_h    Vector of mesh spacing in each direction on
 *                      the coarse mesh in computation space.
 *//*-----------------------------------------------------------------*/

void
FourthOrderMappedFineInterp::defineStencils(
  const IntVect&  a_refRatio,
  const RealVect& a_h)
{
  // Ensure memory is cleared if redefining the stencils.
  if (isStencilDefined())
    {
      delete[] m_fnXiLSx;
      delete[] m_fnGradXiLSx;
    }

  m_refRatio = a_refRatio;
  m_mapStencilIndex.clear();

  // Number of powers.  This is number of columns in 'A'.
  const int lenN = binomial(m_degree + SpaceDim, SpaceDim);
  // Maximum number of cells in a stencil
  const int lenMMax = ((int)std::pow((float)3, SpaceDim)) + 2*SpaceDim;
  int seqIndex = 0;

//--Loop through all the stencils.

  // Mark which indices share the same dx
  m_uniqueDx[0] = 0;
  m_numUniqueDx = 1;
  for (int i = 1; i != SpaceDim; ++i)
    {
      int j = 0;
      while (j != i && a_h[i] != a_h[j]) ++j;
      m_uniqueDx[i] = (j == i) ? m_numUniqueDx++ : m_uniqueDx[j];
    }

  // Predict the number of stencils
  m_numStencil = 1;
  for (int i = 0; i != m_numUniqueDx; ++i)
    {
      int c = 0;
      for (int j = 0; j != SpaceDim; ++j)
        {
          if (m_uniqueDx[j] == i) ++c;
        }
      m_numStencil *= binomial(c + m_stencilExtent, m_stencilExtent);
    }
  
  // Set the maximum matrix size in the allocator
  m_At.getAllocator().define(lenN, lenMMax);
  // Now define the array of matrices
  m_At.define(m_numStencil);

//--Compute powers of displacements for all fine cells in a coarse cell

  RealVect fineDx;
  for (int i = 0; i != SpaceDim; ++i)
    {
      fineDx[i] = 1./m_refRatio[i];
    }
  // The actual box where we need <xi>^p and the gradients
  Box fineCellBox(IntVect::Zero, m_refRatio - IntVect::Unit);
  m_avgXipFine.define(lenN, fineCellBox);
  m_avgXiKpFine.define(lenN, fineCellBox);
  m_avgGradXipFine.define(lenN*SpaceDim, fineCellBox);
  BoxIterator fineCellBIt1(fineCellBox);
  for (fineCellBIt1.begin(); fineCellBIt1.ok(); ++fineCellBIt1)
    {
      const IntVect& fineCellIV1 = fineCellBIt1();
      RealVect rdispl;
      for (int i = 0; i != SpaceDim; ++i)
        {
          // The displacement from the center of the coarse cell
          rdispl[i] = (fineCellIV1[i] - 0.5*(m_refRatio[i] - 1))*fineDx[i];
        }
      loadAvgXipA(&m_avgXipFine(0, fineCellIV1), rdispl, fineDx);
      loadAvgXiKpA(&m_avgXiKpFine(0, fineCellIV1), rdispl, fineDx);
      loadAvgGradXip(&m_avgGradXipFine(0, fineCellIV1), rdispl, fineDx);
    }

//--Prepare matrices for each stencil

  // 'At' with constant K is temporary.  This is the transpose of matrix A
  // for solving \f$A\mathbf{x} = \mathbf{b}\f$.  I.e., the powers of the
  // displacements \f$\delta\xi\f$ to the neighbour cells.  This matrix is
  // created from XiK^p on the coarse cells.
  CHMatrix AKt(lenN, lenMMax);

  m_fnXiLSx = new CHArray<Real, SpaceDim+1, ArRangeCol>[m_numStencil];
  for (int i = 0; i != m_numStencil; ++i)
    {
      m_fnXiLSx[i].define(lenMMax, fineCellBox);
    }
  m_fnGradXiLSx = new CHArray<Real, SpaceDim+1, ArRangeCol>[m_numStencil];
  for (int i = 0; i != m_numStencil; ++i)
    {
      m_fnGradXiLSx[i].define(lenMMax*SpaceDim, fineCellBox);
    }

  // The maximum stencil extent is 1 beyond the standard interior
  // stencil (see Fig 3c for 2D in "Fourth-order conservation law example")
  const int maxStencilExtent = m_stencilExtent + 1;
  const Box maxStencilBox(-maxStencilExtent*IntVect::Unit,
                          maxStencilExtent*IntVect::Unit);
  m_ivsStencilMask.define(m_numStencil);

  // soi[] (Stencil Offset from Interior) is the displacement of the center
  // of the stencil from the standard interior stencil
  int soi[SpaceDim];

  // This describes the dependence between the dimensions.  Dimensions with
  // different dx are independent and loop to the full stencil extent.
  // Dimensions with the same dx can be interchanged so that order is not
  // important --- they are "combinations" so only consider
  // i0 >= i1 >= i2 >= i3 >= i4 >= i5
  // for all dimensions with the same dx
  const int *dimDependence[SpaceDim];
  dimDependence[0] = &m_stencilExtent;
  for (int i = 1; i != SpaceDim; ++i)
    {
      int j = i-1;
      while (j >= 0 && a_h[i] != a_h[j]) --j;
      dimDependence[i] = (j < 0) ? &m_stencilExtent : &soi[j];
    }
  D_TERM6(for (soi[0] = 0; soi[0] <= (*dimDependence[0]); ++soi[0]),
          for (soi[1] = 0; soi[1] <= (*dimDependence[1]); ++soi[1]),
          for (soi[2] = 0; soi[2] <= (*dimDependence[2]); ++soi[2]),
          for (soi[3] = 0; soi[3] <= (*dimDependence[3]); ++soi[3]),
          for (soi[4] = 0; soi[4] <= (*dimDependence[4]); ++soi[4]),
          for (soi[5] = 0; soi[5] <= (*dimDependence[5]); ++soi[5]))
    {
      // Store the index to this matrix
      m_mapStencilIndex.insert(std::make_pair(offset2Key(soi), seqIndex));

//--Mark cells in this stencil.

      // Walls appear in the positive direction and push the stencil the
      // other way.

      // Use an IntVectSet to save the stencil
      m_ivsStencilMask(seqIndex).define(maxStencilBox);
      m_ivsStencilMask(seqIndex).makeEmptyBits();
      // Inner set
      IntVect innerCenter = IntVect::Zero;
      for (int iDir = 0; iDir != SpaceDim; ++iDir)
        {
          if (soi[iDir] == 2)
            {
              innerCenter.shift(iDir, -1);
            }
        }
      m_ivsStencilMask(seqIndex) |= Box(diagShift(innerCenter, -1),
                                        diagShift(innerCenter,  1));
      // Outer set
      for (int iDir = 0; iDir != SpaceDim; ++iDir)
        {
          const int innerShiftDir = innerCenter[iDir];
          // Always add to the -ve direction
          m_ivsStencilMask(seqIndex) |=
            IntVect(IntVect::Zero).shift(iDir, innerShiftDir - 2);
          // Only to +ve direction if no wall
          if (soi[iDir] == 0)
            {
              m_ivsStencilMask(seqIndex) |=
                IntVect(IntVect::Zero).shift(iDir, innerShiftDir + 2);
            }
        }

//--Compute matrix A^T for this stencil

      int idxCenCrCell = -1;  // Seq. index in M of center coarse cell
      // Load A^T (row = powers, col = cells, and column-ordered storage)
      IVSIterator displ(m_ivsStencilMask(seqIndex));
      int iCol = 0;
      for (displ.begin(); displ.ok(); ++displ)
        {
          if (displ() == IntVect::Zero)
            {
              idxCenCrCell = iCol;
            }
          const RealVect rdispl(displ());
          const RealVect rdx(IntVect::Unit);
          int iColIdx = iCol++;
          loadAvgXipA(&m_At(seqIndex)(0, iColIdx), rdispl, rdx);
          loadAvgXiKpA(&AKt(0, iColIdx), rdispl, rdx);
        }
      // Note lenM = iCol

//--Compute Xi (A^TA) AT (constrained least squares intepolation of U).  The
//--constraint is on U, not JU for filling ghost cells.

      int ier;
      FORT_MAPLSATAINVAT(
        CHF_RCHARRAY(RANK_SPACEDIM_PLUS_1, m_fnXiLSx[seqIndex]),
        CHF_RCHARRAY(RANK_SPACEDIM_PLUS_1, m_fnGradXiLSx[seqIndex]),
        CHF_MATRIX(AKt),
        CHF_BOX(fineCellBox),
        CHF_CONST_RCHARRAY(RANK_SPACEDIM_PLUS_1, m_avgXiKpFine),
        CHF_CONST_RCHARRAY(RANK_SPACEDIM_PLUS_1, m_avgGradXipFine),
        CHF_CONST_INT(idxCenCrCell),
        CHF_CONST_INT(iCol),
        CHF_INT(ier));
      if (ier != 0)
        {
          pout() << "Lapack error code: " << ier << std::endl;
          MayDay::Error("Failure in matrix inversion in FORT_MAPINVERSEATA");
        }

      ++seqIndex;
    }  // Loop over stencils
  // Sanity check
  CH_assert(seqIndex == m_numStencil);

  m_defined[0] = true;
}

/*--------------------------------------------------------------------*/
//   Define that allocates data for coarse representations of the fine
//   grid
/**  Note that this data hangs around forever.
 *   \param[in]  a_FnGrid
 *                      Layout of the fine boxes.
 *   \param[in]  a_CrGrid
 *                      Layout of the coarse boxes.
 *   \param[in]  a_CrMetricsPtr
 *                      Used to get problem domain on the coarse grid
 *                      Only needed for single-block and can be null
 *                      otherwise
 *   \param[in]  a_FnInterpRadVec
 *                      Number of fine cells to fill in each direction
 *                      (default 0)
 *   \param[in]  a_willFillGhostsWithCrFnJU
 *                      T - \<JU\> will be provided on the CrFn mesh
 *                          for filling ghost cells.  This implies
 *                          that \<U\> needs to be computed and
 *                          therefore one extra ghost is needed on
 *                          the CrFn layout (default true).  You
 *                          don't _have_ to fill with \<JU\> if set
 *                          true but you definitely cannot if set to
 *                          false.
 *                      F - Memory for \<JU\> on the CrFn mesh is not
 *                          even allocated
 *   \param[in]  a_mbgeoPtr
 *                      Information about the multiblock geometry.  If
 *                      provided, this triggers a hack where the CrFn
 *                      mesh for \<JU\> is actually the coarse mesh.
 *                      This allows us to compute \<U\> before
 *                      invoking a special multiblock copyTo.  Much of
 *                      this hack, relating to time interpolation
 *                      instead of space interpolation, is within this
 *                      class simply because the relevant data
 *                      structures are here.  And since it's probably
 *                      just a temporary hack.
 *   \param[in]  a_ghostVectU
 *                      Only required for multiblock.  This is the
 *                      number of ghosts for <U>.  m_CrLevU must be
 *                      defined with this many ghosts because a
 *                      multi-block exchange will be invoked on the
 *                      level.
 *//*-----------------------------------------------------------------*/

void
FourthOrderMappedFineInterp::defineGrids(
  const DisjointBoxLayout&         a_FnGrid,
  const DisjointBoxLayout&         a_CrGrid,
  LevelGridMetrics *const          a_CrMetricsPtr,
  const IntVect&                   a_FnInterpRadVec,
  const bool                       a_willFillGhostsWithCrFnJU,
  const MultiBlockLevelGeom *const a_mbgeoPtr,
  const IntVect *const             a_ghostVectU)
{
  CH_TIME("FourthOrderMappedFineInterp::defineGrids");
  m_isMultiblock = (a_mbgeoPtr != 0);
  // Ghost radius
  m_FnInterpRadVec = a_FnInterpRadVec;
  D_TERM6(m_CrInterpRadVec[0] = (a_FnInterpRadVec[0] + m_refRatio[0] - 1)/m_refRatio[0];,
          m_CrInterpRadVec[1] = (a_FnInterpRadVec[1] + m_refRatio[1] - 1)/m_refRatio[1];,
          m_CrInterpRadVec[2] = (a_FnInterpRadVec[2] + m_refRatio[2] - 1)/m_refRatio[2];,
          m_CrInterpRadVec[3] = (a_FnInterpRadVec[3] + m_refRatio[3] - 1)/m_refRatio[3];,
          m_CrInterpRadVec[4] = (a_FnInterpRadVec[4] + m_refRatio[4] - 1)/m_refRatio[4];,
          m_CrInterpRadVec[5] = (a_FnInterpRadVec[5] + m_refRatio[5] - 1)/m_refRatio[5];)
  m_willFillGhostsWithCrFnJU = a_willFillGhostsWithCrFnJU;
  //** The following requires an int but we generally assume that refRatio is
  //** a vector
  m_CrFnGrid = DisjointBoxLayout();
  coarsen(m_CrFnGrid, a_FnGrid, m_refRatio[0]);
  int numExtraGhosts = m_stencilExtent;
  if (a_willFillGhostsWithCrFnJU && !isMultiblock())
    {
      // The CrFn JU needs to be 1 ghost cell larger than U because when we are
      // filling ghosts we need to first compute U.  Note however that
      // m_CrFnLevU is allocated with the same number of ghosts so that we can
      // use the same copier.
      // However, if multiblock, we *always* copy U from the coarse grid to the
      // coarsened fine mesh.  Hence, we don't need this extra ghost.  THIS MAY
      // CHANGE ONCE THE WORKAROUND FOR TIME INTERPOLATION IN MULTIBLOCK IS
      // RESOLVED!!!
      numExtraGhosts += 1;
    }
  m_CrFnNumGhost = (m_CrInterpRadVec + numExtraGhosts*IntVect::Unit);
  m_CrFnLevU.define(m_CrFnGrid, m_nComp, m_CrFnNumGhost);
  if (a_willFillGhostsWithCrFnJU)
    {
      if (a_mbgeoPtr)
        {
          // Special hack for multiblock.  Since this covers the entire coarse
          // grid, only 1 ghost cell is need to allow for computing <U>
          m_CrFnLevJU.define(a_CrGrid, m_nComp, IntVect::Unit);
          // But <U> will need enough ghost cells for multi-block exchange
          CH_assert(a_ghostVectU != nullptr);
          m_CrLevU.define(a_CrGrid, m_nComp, *a_ghostVectU);
        }
      else
        {
          // Assuming single block
          m_CrFnLevJU.define(m_CrFnGrid, m_nComp, m_CrFnNumGhost);
        }
    }
  if (isMultiblock())
    {
      // Make the disjoint source include all extra-block ghosts.  The boxes
      // are still disjoint since the blocks are spaced throughout the domain
      // to provide unique indexing for the ghost cells.
      DisjointBoxLayout modCrGrid;
      // Make a deep copy of the original coarse layout
      modCrGrid.deepCopy(a_CrGrid);
      Vector<Entry>& layoutEntries = *modCrGrid.rawPtr();
      MultiBlockCoordSys& coordSys = *a_mbgeoPtr->coordSysPtr();
      
      for (int idx = 0, idx_end = modCrGrid.size(); idx != idx_end; ++idx)
        {
          Box box = layoutEntries[idx].box;
          const Box& domainBox = coordSys.problemDomain(box).domainBox();
          if (!domainBox.contains(grow(box, m_CrFnNumGhost)))
            {
              Box modBox(box);
              // Grow along faces adjacent to the block boundary
              for (const auto side : EachSide)
                {
                  const IntVect ivSideBox    = box.sideEnd(side);
                  const IntVect ivSideDomain = domainBox.sideEnd(side);
                  for (const int dir : EachDir)
                    {
                      if (ivSideBox[dir] == ivSideDomain[dir])
                        {
                          modBox.growDir(dir, side, m_CrFnNumGhost[dir]);
                        }
                    }
                }
              layoutEntries[idx].box = modBox;
            }
        }
      
      // Close the layout using same neighbor information as previous.  This
      // will re-sort the boxes.
      modCrGrid.closeN(a_CrGrid);

      // Now define the copier
      m_copier.define(modCrGrid, m_CrFnGrid, m_CrFnNumGhost);
    }
  else
    {
      const ProblemDomain& crProblemDomain =
        a_CrMetricsPtr->getCoordSys().problemDomain(0);
      m_copier.define(a_CrGrid, m_CrFnGrid, crProblemDomain, m_CrFnNumGhost);
    }
  m_defined[1] = true;
  m_defined[2] = false;              // Invalidates CFInterface.
  m_defined[3] = false;              // Invalidates m_CrFnLevJ.,
  m_CrFnStatus = CrFnStatusInvalid;  // Invalidate CrFn Data
}

/*--------------------------------------------------------------------*/
//   Define that saves the coarse <J> on the coarsened-fine layout
/**  This is useful if we are given \<JU\> on the coarsened-fine
 *   and need to convert it to </U\> before interpolating to fill any
 *   ghost cells.  The grid must be defined so we know the coarsened-
 *   fine layout.  It is necessary to call this routine before using
 *   'fillGhostsFromCrJU'
 *   \param[in]  a_CrJ  \<J\> from the coarse layout
 *   \note
 *   <ul>
 *     <li> You must have indicated that you will fill ghosts with
 *          CrFnJU when defineGrids was called.  This ensures there
 *          are sufficient layers of ghosts to calculate \<U\> on all
 *          cells in a stencil.
 *   </ul>
 *//*-----------------------------------------------------------------*/

void
FourthOrderMappedFineInterp::defineCrFnJ(const LevelData<FArrayBox>& a_CrJ)
{
  CH_TIME("FourthOrderMappedFineInterp::defineCrFnJ");
  CH_assert(isGridDefined());
  CH_assert(m_willFillGhostsWithCrFnJU);
  if (isMultiblock())
    {
      // Special hack for multiblock.
      LevelData<FArrayBox> *const fullCrJPtr =
        const_cast<LevelData<FArrayBox>*>(&a_CrJ);
      aliasLevelData(m_CrFnLevJ, fullCrJPtr, a_CrJ.interval());
    }
  else
    {
      m_CrFnLevJ.define(m_CrFnLevJU.getBoxes(), 1, m_CrFnLevJU.ghostVect());
      a_CrJ.copyTo(m_CrFnLevJ, m_copier);
    }
  m_defined[3] = true;
}

/*--------------------------------------------------------------------*/
//   Define that initializes the coarse-fine interface for filling
//   ghosts
/**  \param[in]  a_CrMetricsPtr
 *                      Used to get problem domain on the coarse mesh
 *//*-----------------------------------------------------------------*/

void
FourthOrderMappedFineInterp::defineCFInterface(
  LevelGridMetrics *const a_CrMetricsPtr)
{
  CH_TIME("FourthOrderMappedFineInterp::defineCFInterface");
  // Need some stencil data
  CH_assert(isStencilDefined());
  // The grid has to be defined to use this functionality
  CH_assert(isGridDefined());

//--Initializations

  DataIterator dit = m_CrFnGrid.dataIterator();

  bool isSorted = m_CrFnGrid.isSorted();

  // Re-definitions
  m_numCrFnFillGhost = 0;
  m_CrFnFillInterp.define(m_CrFnGrid);
  m_ArIdxFillInterp.define(m_CrFnGrid);
  for (dit.begin(); dit.ok(); ++dit)
    {
      m_ArIdxFillInterp[dit].define(grow(m_CrFnGrid[dit], m_CrInterpRadVec), 1);
    }

  // A box which will determine whether a given box adjoins a periodic boundary
  Box crPeriodicTestBox;
  bool isPeriodicDomain = false;
  if (!isMultiblock())
    {
      CH_assert(!a_CrMetricsPtr->getCoordSys().isMultiBlock());
      const ProblemDomain& crProblemDomain =
        a_CrMetricsPtr->getCoordSys().problemDomain(0);
      crPeriodicTestBox = crProblemDomain.domainBox();
      if (crProblemDomain.isPeriodic())
        {
          isPeriodicDomain = true;
          for (int iDir = 0; iDir < SpaceDim; ++iDir)
            {
              if (crProblemDomain.isPeriodic(iDir))
                {
                  crPeriodicTestBox.grow(iDir, -1);
                }
            }
        }
    }

//--Mark all coarse cells on the coarsened fine layout that will be used to fill
//--ghost cells on the fine mesh

  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box crFnBox = m_CrFnGrid[dit];
      int idxBlk = 0;
      if (isMultiblock())
        {
          idxBlk = m_CrFnGrid.blockIndex(dit);
          CH_assert(idxBlk != -1);
        }
      const ProblemDomain crBlockDomain =
        a_CrMetricsPtr->getCoordSys().problemDomain(idxBlk);

      const Box crFnGhostBox =
        grow(crFnBox, m_CrInterpRadVec) & crBlockDomain;
      IntVectSet& crFnFillGhostsIVS = m_CrFnFillInterp[dit];

      // Initialize the IVS as the whole crFnGhostBox
      crFnFillGhostsIVS.define(crFnGhostBox);

      // Check if any of the ghost cells extend across a periodic domain
      bool notAcrossPeriodic = true;
      if (isPeriodicDomain)
        {
          notAcrossPeriodic =
            crBlockDomain.domainBox().contains(crFnGhostBox);
        }

      // Now we iterate over all the boxes in m_CrFnGrid and subtract off all
      // the real cells leaving only the invalid ghost cells
      LayoutIterator lit = m_CrFnGrid.layoutIterator();
      for (lit.begin(); lit.ok(); ++lit)
        {
          const Box crFnRealBox = m_CrFnGrid[lit];

          if (notAcrossPeriodic && isSorted &&
              (crFnRealBox.bigEnd(0) < crFnGhostBox.smallEnd(0)))
            {
              // Nothing from this box
              continue;
            }

          if (notAcrossPeriodic && isSorted &&
              (crFnRealBox.smallEnd(0) > crFnGhostBox.bigEnd(0)))
            {
              // Can break out of loop, since we know that the smallEnd of all
              // the remaining boxes are lexigraphically beyond the ghost box.
              break;
            }

          crFnFillGhostsIVS -= crFnRealBox;

          // We also need to remove any periodic images of real cells
          if (isPeriodicDomain
              && !crPeriodicTestBox.contains(crFnRealBox)
              && !crPeriodicTestBox.contains(crFnGhostBox))
            {
              ShiftIterator shiftIt = crBlockDomain.shiftIterator();
              IntVect shiftMult(crBlockDomain.domainBox().size());
              Box shiftedBox(crFnRealBox);
              for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
                {
                  IntVect shiftVect = shiftMult*shiftIt();
                  shiftedBox.shift(shiftVect);
                  crFnFillGhostsIVS -= shiftedBox;
                  shiftedBox.shift(-shiftVect);
                }
            }
        }  // Second loop over all boxes

      /* Consider block boundaries:
         - Currently we interpolate to any fine box adjacent to an AMR interface
           coincident with a block boundary.  If there are fine cells on the
           other side, MBEx will overwrite interpolated values.  This may be a
           good strategy.
         - First, if the disjoint crFn box is not
           adjacent to the block boundary, no ghosts are cropped (all fine cells
           are interpolated).  This mean all invalid ghost cells are filled by
           level interpolation if there is any adjacency to a block boundary.
           However, if there is a valid cell on the other side, MBEx might
           overwrite this.
         - Otherwise, crop the ghosts if all fine cells pair with fine cells
           across the boundary.  MultiBlockRegions stores this information.
      */
//       if (isMultiblock())
//         {
//           const testBlock = grow(crBlockDomain.box(), -1);
//           if (!testBlock.contains(crFnBox))
//             {
//               stc::Vector<Vector<Box>, 2*SpaceDim>& adjFn =
//                 a_fnMBRegions.getLevelLocation(dit());
//               using FaceTag = MultiBlockRegions::FaceTag;
//               for (const int dir : EachDir)
//                 {
//                   for (const auto side : EachSide)
//                     {
//                       const Vector<Box>& adjFn = a_fnMBRegions.getLevelLocations
//                     }
//                 }
//             }
//         }

      // Set a sequential index for all coarse ghost cells in the IVS:
      BaseFab<int>& arIdxFillInterp = m_ArIdxFillInterp[dit];
      arIdxFillInterp.setVal(-1);
      IVSIterator crFnFillGhostsIt(crFnFillGhostsIVS);
      for (crFnFillGhostsIt.begin(); crFnFillGhostsIt.ok(); ++crFnFillGhostsIt)
        {
          arIdxFillInterp(crFnFillGhostsIt()) = m_numCrFnFillGhost++;
        }
    }  // First loop over boxes on the processor

//--Set a box of fine ghost cells to be filled for each coarse cell

  m_FnBoxFillInterp.define(m_numCrFnFillGhost);
  for (dit.begin(); dit.ok(); ++dit)
    {
      Box fnGhostBox =
        grow(refine(m_CrFnGrid[dit], m_refRatio), m_FnInterpRadVec);
      BaseFab<int>& arIdxFillInterp = m_ArIdxFillInterp[dit];
      IVSIterator crFnFillGhostsIt(m_CrFnFillInterp[dit]);
      for (crFnFillGhostsIt.begin(); crFnFillGhostsIt.ok(); ++crFnFillGhostsIt)
        {
          const IntVect crFnFillGhostIV = crFnFillGhostsIt();
          int arIdx = arIdxFillInterp(crFnFillGhostIV);
          m_FnBoxFillInterp(arIdx) = fnGhostBox &
            refine(Box(crFnFillGhostIV, crFnFillGhostIV), m_refRatio);
        }
    }

  m_defined[2] = true;
  m_CrFnStatus = CrFnStatusInvalid;  // Invalidate CrFn Data
}

/*--------------------------------------------------------------------*/
// Set options for interpolation clipping
/**
 * \param[in] m_useClipping
 *                      Clip solutions to prevent solution overshoot
 *                      maintains conservation
 * \param[in] a_useHOClipping
 *                      Clipping not applied to smooth data
 * \param[in] a_usePostSmoothing
 *                      Smooth the interpolation after clipping
 */
/*--------------------------------------------------------------------*/
void
FourthOrderMappedFineInterp::defineClipping(const bool a_useClipping,
                                            const bool a_useHOClipping,
                                            const bool a_usePostSmoothing)
{
  m_useClipping = a_useClipping;
  m_useHOClipping = a_useHOClipping;
  m_usePostSmoothing = a_usePostSmoothing;
}

/*--------------------------------------------------------------------*/
//  Interpolate from a coarse to fine level on a mapped grid
/** The copier will be used to get the coarse data onto a coarse
 *  representation of the fine mesh.  If the grids have been defined,
 *  the coarse data will be placed there.  Otherwise it will be
 *  placed in a temporary and deallocated.
 *  \param[out] a_FnLevJU
 *                      \<JU\> interpolated on the fine grid
 *  \param[in]  a_FnLevJ
 *                      The metrics \<J\> on the fine grid
 *  \param[in]  a_CrLevU
 *                      \<U\> on the coarse grid.  If multi-block, the
 *                      extra-block ghost must have been exchanged.
 *  \param[in]  a_CrLevJU
 *                      \<JU\> on the coarse grid
 *  \param[in]  a_CrMetricsPtr
 *                      Used to get problem domain on the coarse grid
 *  \param[in]  a_vectorIntv
 *                      An interval in components of U consisting of
 *                      vector data
 *  \note
 *   <ul>
 *     <li> Stencils must have been defined.
 *     <li> If the grids have been defined, the coarse data will be
 *          copied there.
 *     <li> If the grids have not been defined, the coarse data will
 *          be copied to a temporary and deleted upon exit.
 *   </ul>
 *//*-----------------------------------------------------------------*/

void
FourthOrderMappedFineInterp::interpToFine(
  LevelData<FArrayBox>&       a_FnLevJU,
  const LevelData<FArrayBox>& a_FnLevJ,
  const LevelData<FArrayBox>& a_CrLevU,
  const LevelData<FArrayBox>& a_CrLevJU,
  LevelGridMetrics *const     a_CrMetricsPtr,
  const Interval&             a_vectorIntv,
  const Interval&             a_interpInterval)
{
  CH_TIME("FourthOrderMappedFineInterp::interpToFine");
  CH_assert(isStencilDefined());
  //**FIXME m_nComp and m_CrFnNumGhost are defined, maybe use those?
  const int nComp   = a_FnLevJU.nComp();
  const IntVect numGhost(m_stencilExtent*IntVect::Unit);

  LevelData<FArrayBox> *CrFnLevU;
  LevelData<FArrayBox> *CrFnLevJU;
  Copier *copier;
  // Modifiable
  DisjointBoxLayout *CrFnGridPtr = 0;
  // Unmodifiable if we only reference
  const DisjointBoxLayout *CrFnGrid = 0;
  if (!isGridDefined())
    {
      // Proceeding without a grid defined is only possible for single-block
      CH_assert(!isMultiblock());
      CH_assert(!a_CrMetricsPtr->getCoordSys().isMultiBlock());
      CrFnGridPtr = new DisjointBoxLayout;
      CrFnLevU = new LevelData<FArrayBox>;
      //** The following requires an int but we generally assume that refRatio
      //** is a vector
      coarsen(*CrFnGridPtr, a_FnLevJU.getBoxes(), m_refRatio[0]);
      CrFnGrid = CrFnGridPtr;
      CrFnLevU->define(*CrFnGrid, nComp, numGhost);
      CrFnLevJU = new LevelData<FArrayBox>;
      CrFnLevJU->define(*CrFnGrid, nComp, numGhost);
      const ProblemDomain& crProblemDomain =
        a_CrMetricsPtr->getCoordSys().problemDomain(0);
      copier = new Copier(a_CrLevU.getBoxes(),
                          *CrFnGrid,
                          crProblemDomain,
                          numGhost);
    }
  else
    {
      CrFnGrid = &m_CrFnGrid;
      CrFnLevU = &m_CrFnLevU;
      CrFnLevJU = &m_CrFnLevJU;
      copier = &m_copier;
    }
  if (isMultiblock())
    {
      // Special hack for multiblock.  CrFnLevJU is redefined because the stored
      // m_CrFnLevJU covers the complete coarse grid and uses the coarse-grid
      // layout.
      CrFnLevJU = new LevelData<FArrayBox>;
      CrFnLevJU->define(*CrFnGrid, m_nComp, IntVect::Zero);
    }

  // Check for proper nesting FIXME :: doesn't work for MMB
  CH_assert(properNesting(a_CrLevJU.getBoxes(),
                          *CrFnGrid,
                          a_CrMetricsPtr));

  // Copy data to the coarsened fine grid
  // For multiblock, <JU> is only required in the valid cells.  We cannot use
  // the copier since it expects that extra-block ghosts are part of the
  // disjoint layout.  To copy <U>, we use the presetCrFnLevU routine
  // which copies to m_CrFnLevU and requires that extra-block ghosts have been
  // filled
  if (isMultiblock())
    {
      presetCrFnLevU(a_CrLevU, a_vectorIntv);
      a_CrLevJU.copyTo(*CrFnLevJU);
    }
  else
    {
      a_CrLevU.copyTo(*CrFnLevU, *copier);
      a_CrLevJU.copyTo(*CrFnLevJU, *copier);
    }

  // Now continue with the least-squares interpolation
  interpCrFnToFine(
    a_FnLevJU,
    a_FnLevJ,
    *CrFnGrid,
    *CrFnLevU,
    *CrFnLevJU,
    a_CrMetricsPtr,
    a_interpInterval);

  // Clean up memory if we allocated it
  if (!isGridDefined())
    {
      delete CrFnGridPtr;
      delete CrFnLevU;
      delete CrFnLevJU;
      delete copier;
    }
  else if (isMultiblock())
    {
      // Special hack for multiblock.
      delete CrFnLevJU;
    }
}

/*--------------------------------------------------------------------*/
//   Fill invalid ghost cells from preset coarsened-fine data on a
//   mapped grid
/**  To call this, one of the "::presetCrFnLev*" routines must have
 *   been used to get the coarse data on a coarse representation of
 *   the fine mesh.  Also, the coarse-fine interface must have been
 *   identified.
 *   \param[out] a_FnLevU
 *                      \<U\> interpolated on the fine grid
 *   \param[out] a_FnLevJU
 *                      \<JU\> interpolated on the fine grid.  Can
 *                      have less ghosts than a_FnLevU -- only updated
 *                      for cells that exist
 *   \param[in]  a_FnLevJ
 *                      The metrics \<J\> on the fine grid
 *   \param[in]  a_CrMetricsPtr
 *                      Used to get problem domain on the coarse grid
 *
 *  \note
 *   <ul>
 *     <li> Stencils must have been defined.
 *     <li> The internal coarsened fine representation of \<U\> must
 *          have been set either by 'fillGhosts()' or using
 *          'presetCrFnLevU()' or 'presetCrFnLevJU()' to set it
 *          externally.
 *   </ul>
 *//*-----------------------------------------------------------------*/

void
FourthOrderMappedFineInterp::fillGhosts(
  LevelData<FArrayBox>&       a_FnLevU,
  LevelData<FArrayBox>&       a_FnLevJU,
  const LevelData<FArrayBox>& a_FnLevJ,
  LevelGridMetrics *const     a_CrMetricsPtr)
{
  CH_TIME("FourthOrderMappedFineInterp::fillGhosts");
  CH_assert(isStencilDefined());
  CH_assert(isGridDefined());
  CH_assert(isCFInterfaceDefined());
  if (m_CrFnStatus == CrFnStatusHaveJU)
    {
      // Deferred calculations of <U> is only possible for single-block
      CH_assert(!isMultiblock());
      CH_assert(!a_CrMetricsPtr->getCoordSys().isMultiBlock());
      const ProblemDomain& crProblemDomain =
        a_CrMetricsPtr->getCoordSys().problemDomain(0);
      // This means we have <JU>, but not <U> so compute <U>
      for (DataIterator dit = m_CrFnLevU.dataIterator(); dit.ok(); ++dit)
        {
          // Compute CrFn <U> in the domain
          FArrayBox& crFnFabU = m_CrFnLevU[dit];
          const FArrayBox& crFnFabJU = m_CrFnLevJU[dit];
          const Box UBox = grow(crFnFabJU.box(), -1) & crProblemDomain;
          cellFGToCellF(crFnFabU,
                        crFnFabJU,
                        m_CrFnLevJ[dit],
                        UBox,
                        crProblemDomain,
                        true);
        }
      m_CrFnStatus |= CrFnStatusHaveU;
    }
  CH_assert(m_CrFnStatus & CrFnStatusHaveU);

  // Number of ghosts required for full interpolation stencil
  // int numGhostStencil = m_CrInterpRadVec[0] + m_stencilExtent;

  // Loop over the boxes and call fillGhosts
  for (DataIterator dit = m_CrFnGrid.dataIterator(); dit.ok(); ++dit)
    {
      const Box crFnBox = m_CrFnGrid[dit];
      const ProblemDomain& crBlockDomain =
        a_CrMetricsPtr->getCoordSys().problemDomain(crFnBox);
      fillGhosts(a_FnLevU[dit],
                 a_FnLevJU[dit],
                 a_FnLevJ[dit],
                 m_CrFnFillInterp[dit],
                 m_ArIdxFillInterp[dit],
                 m_CrFnLevU[dit],
                 crBlockDomain);
    }
}

/*--------------------------------------------------------------------*/
//   Interpolate from a coarsened-fine to fine level on a mapped grid
/**  To call this, the copier must have already been used to get the
 *   coarse data on a coarse representation of the fine mesh.
 *   \param[out] a_FnLevJU
 *                      \<JU\> interpolated on the fine grid
 *   \param[in]  a_FnLevJ
 *                      The metrics \<J\> on the fine grid
 *   \param[in]  a_CrFnGrid
 *                      A coarse representation of the fine grid.
 *                      This is required instead of m_CrFnGrid
 *                      directly because the grid may have been
 *                      temporarily defined in routine interpToFine
 *   \param[in]  a_CrFnLevU
 *                      \<U\> on a coarse representation of the fine
 *                      grid
 *   \param[in]  a_CrFnLevJU
 *                      \<JU\> on a coarse representation of the fine
 *                      grid
 *   \param[in]  a_CrMetricsPtr
 *                      Used to get problem domain on the coarse grid
 *
 *   \note
 *   <ul>
 *     <li> Stencils must have been defined.
 *     <li> The internal coarsened fine grids are unused and do not
 *          need to be defined.
 *   </ul>
 *//*-----------------------------------------------------------------*/

void
FourthOrderMappedFineInterp::interpCrFnToFine(
  LevelData<FArrayBox>&       a_FnLevJU,
  const LevelData<FArrayBox>& a_FnLevJ,
  const DisjointBoxLayout&    a_CrFnGrid,
  const LevelData<FArrayBox>& a_CrFnLevU,
  const LevelData<FArrayBox>& a_CrFnLevJU,
  LevelGridMetrics *const     a_CrMetricsPtr,
  const Interval&             a_interpInterval) const
{
  CH_TIME("FourthOrderMappedFineInterp::interpCrFnToFine");
  CH_assert(isStencilDefined());

  // Simply loop over the boxes and call interpolate
  for (DataIterator dit = a_CrFnGrid.dataIterator(); dit.ok(); ++dit)
    {
      const int numInterpComps = a_interpInterval.size();
      const int numJUComps = a_FnLevJU.nComp();
      int ier;
      if (numInterpComps < numJUComps)
        {
          // Hack for SGS KE interpolation -- interp only some components
          const int numInterpComps = a_interpInterval.size();
          const int srcComp = a_interpInterval.begin();
          FArrayBox FnLevJU((a_FnLevJU[dit]).box(), numInterpComps);
          FArrayBox CrFnLevU((a_CrFnLevU[dit]).box(), numInterpComps);
          FArrayBox CrFnLevJU((a_CrFnLevJU[dit]).box(), numInterpComps);
          FnLevJU.copy(a_FnLevJU[dit], srcComp, 0, numInterpComps);
          CrFnLevU.copy(a_CrFnLevU[dit], srcComp, 0, numInterpComps);
          CrFnLevJU.copy(a_CrFnLevJU[dit], srcComp, 0, numInterpComps);

          const Box& crFnBox = a_CrFnGrid[dit];
          const ProblemDomain& crBlockDomain =
            a_CrMetricsPtr->getCoordSys().problemDomain(crFnBox);
          ier = interpToFine(FnLevJU,
                             a_FnLevJ[dit],
                             crFnBox,
                             CrFnLevU,
                             CrFnLevJU,
                             crBlockDomain);
         a_FnLevJU[dit].copy(FnLevJU, 0, srcComp, numInterpComps);
        }
      else
        {
          const Box& crFnBox = a_CrFnGrid[dit];
          const ProblemDomain& crBlockDomain =
            a_CrMetricsPtr->getCoordSys().problemDomain(crFnBox);
          ier = interpToFine(a_FnLevJU[dit],
                             a_FnLevJ[dit],
                             crFnBox,
                             a_CrFnLevU[dit],
                             a_CrFnLevJU[dit],
                             crBlockDomain);
        }
      if (ier)
        {
          MayDay::Error("Constrained least-squares interpolation failed");
        }
    }
}

/*--------------------------------------------------------------------*/
//   Interpolate from a coarse to fine box on a mapped grid
/**  An appropriate stencil is selected depending on the proximity of
 *   a coarse cell to a boundary.  Data required for the least-squares
 *   problem is packed for the stencil before handing off to a Fortran
 *   routine.
 *   \param[out] a_FnFabJU
 *                      Interpolation of \<JU\> on the fine cells
 *   \param[in]  a_FnFabJ
 *                      Metrics on the fine mesh
 *   \param[in]  a_CrBox
 *                      Coarse box to interpolate from, a coarsened
 *                      representation of the fine box
 *   \param[in]  a_CrFabU
 *                      Coarse \<U\> on the box.  Mapped to LSbl for
 *                      each LS problem.
 *   \param[in]  a_CrFabJU
 *                      Coarse \<JU\> on the box.  This is the
 *                      constraint 'LSd' for each coarse cell.
 *   \param[in]  a_CrProblemDomain
 *                      The problem domain as specified on the coarse
 *                      level
 *   \return            Error status.  >0 means a Lapack routine
 *                      returned an error code.
 *                       0 - success
 *                      >0 - failure
 *//*-----------------------------------------------------------------*/

int
FourthOrderMappedFineInterp::interpToFine(
  FArrayBox&           a_FnFabJU,
  const FArrayBox&     a_FnFabJ,
  const Box&           a_CrBox,
  const FArrayBox&     a_CrFabU,
  const FArrayBox&     a_CrFabJU,
  const ProblemDomain& a_CrProblemDomain) const
{
  CH_TIME("FourthOrderMappedFineInterp::interpToFine");
  // Some definitions
  const int lenMMax = m_At(0).size(1);
  const int nComp   = a_FnFabJU.nComp();

  // Some temporaries
  CHArray<Real, 2> LSbl(lenMMax, nComp);
  Vector<Real> LSd(nComp);
  int iStencil;             // Stencil index
  CoordTransform ctStoP;    // Transformation from stencil coordinate sytem
                            // to physical coordinate system
  CoordTransform ctPtoS;    // Transformation from physical coordinate sytem
                            // to stencil coordinate system
  IntVect offset;           // Offset between (0,0) in S and P coordinates
  IntVect fineSmlEnd;       // What to add to a fine cell location transformed
                            // from S to P to get the actual location in Z

  // A box of valid fine cells where fine-grid data can be accessed
  ProblemDomain FnProblemDomain(a_CrProblemDomain);
  FnProblemDomain.refine(m_refRatio);
  Box validBox(a_FnFabJU.box());
  validBox &= FnProblemDomain;
  
  int ierBox = 0;  // Accumulated over the box
  // For each coarse cell
  BoxIterator bIt(a_CrBox);
  for (bIt.begin(); bIt.ok(); ++bIt)
    {
      const IntVect center = bIt();

      // Select a stencil and get indexes and transformation info
      selectStencil(center, a_CrProblemDomain, iStencil, ctStoP, ctPtoS, offset,
                    fineSmlEnd);

      // Get the refinement ratio in the stencil coordinate system
      const IntVect refRatioS = absolute(ctPtoS.transform(m_refRatio));

      // Box of fine cells to update (in stencil space)
      Box fineBoxInCrS(IntVect::Zero, refRatioS - IntVect::Unit);

      // Pack b = <U>_{j^l} according to the stencil
      int iN = 0;
      IVSIterator displ(m_ivsStencilMask(iStencil));
      for (displ.begin(); displ.ok(); ++displ)
        {
          IntVect tdispl = ctStoP.transform(displ());
          tdispl += center;
          for (int iComp = 0; iComp != nComp; ++iComp)
            {
              LSbl(iN, iComp) = a_CrFabU(tdispl, iComp);
            }
          ++iN;
        }

      // Pack d = <JU>_{i^l}
      for (int iComp = 0; iComp != nComp; ++iComp)
        {
          LSd[iComp] = a_CrFabJU(center, iComp);
        }

      // Stencil size
      const int lenM = m_ivsStencilMask(iStencil).numPts();

      int ier;
      FORT_MAPLSINTERP(
        CHF_FRA(a_FnFabJU),
        CHF_CONST_FRA(a_FnFabJ),
        CHF_BOX(validBox),
        CHF_CONST_I1D(ctStoP.xctm(), 2*SpaceDim),
        CHF_BOX(fineBoxInCrS),
        CHF_CONST_INTVECT(fineSmlEnd),
        CHF_CONST_RCHARRAY(2, LSbl),
        CHF_CONST_R1D(&LSd[0], nComp),
        CHF_CONST_MATRIX(m_At(iStencil)),
        CHF_CONST_RCHARRAY(RANK_SPACEDIM_PLUS_1, m_avgXipFine),
        CHF_CONST_RCHARRAY(RANK_SPACEDIM_PLUS_1, m_avgGradXipFine),
        CHF_CONST_INT(lenM),
        CHF_INT(ier));
      // mapLSinterp(a_FnFabJU,
      //             a_FnFabJ,
      //             validBox,
      //             ctStoP,
      //             fineBoxInCrS,
      //             fineSmlEnd,
      //             LSbl,
      //             LSd,
      //             m_At(iStencil),
      //             m_avgXipFine,
      //             m_avgGradXipFine,
      //             lenM,
      //             ier);
      if (ier != 0)
        {
          if (ierBox == 0)
            {
              ierBox = -1;
              MayDay::Warning("Possible failure in constrained least-squares "
                              "computation in FORT_MAPLSINTERP");
              pout() << "Error code: " << ier << std::endl;
            }
          ierBox = std::max(ierBox, ier);
          // But continue anyways
        }
    }  // Loop over coarse cells

  // Clip JU
  if (m_useClipping)
    {
      applyMappedClipping(a_FnFabJU, a_FnFabJ, a_CrBox, a_CrFabU, a_CrFabJU, a_CrProblemDomain);
    }

  return ierBox;
}

/*--------------------------------------------------------------------*/
//   Fill any ghost cells on a finer mesh
/**  An appropriate stencil is selected depending on the proximity of
 *   a coarse cell to a boundary.  Data required for the least-squares
 *   problem is packed for the stencil before handing off to a Fortran
 *   routine.
 *   \param[out] a_FnFabU
 *                      Interpolation of \<U\> into the invalid fine
 *                      ghost cells
 *   \param[out] a_FnFabJU
 *                      Interpolation of \<JU\> into the invalid fine
 *                      ghost cells.  Can have less ghosts than
 *                      a_FnLevU -- only updated for cells that exist
 *   \param[in]  a_FnFabJ
 *                      Metrics on the fine mesh
 *   \param[in]  a_CrCells
 *                      Coarse cells underlying the fine ghost cells
 *                      in this box
 *   \param[in]  a_CrCellsArIdx
 *                      A basefab giving an array index for all the
 *                      cells in a_CrCells.  This is the location of
 *                      any auxiliary data associated with the coarse
 *                      cell.
 *   \param[in]  a_CrFabU
 *                      Coarse \<U\> on the box.  Mapped to LSbl for
 *                      each LS problem.
 *   \param[in]  a_CrProblemDomain
 *                      The problem domain as specified on the coarse
 *                      level
 *   \note
 *   <ul>
 *     <li> The coarsened fabs and IVS must exist on a coarsened
 *          representation of the fine mesh
 *   </ul>
 *//*-----------------------------------------------------------------*/

void
FourthOrderMappedFineInterp::fillGhosts(
  FArrayBox&           a_FnFabU,
  FArrayBox&           a_FnFabJU,
  const FArrayBox&     a_FnFabJ,
  const IntVectSet&    a_CrCells,
  const BaseFab<int>&  a_CrCellsArIdx,
  const FArrayBox&     a_CrFabU,
  const ProblemDomain& a_CrProblemDomain) const
{
  CH_TIME("FourthOrderMappedFineInterp::fillGhosts");
  // There is a requirement here that J and JU have the same number of ghost
  // cells.  Both boxes should have equal sizes.
  CH_assert(a_FnFabJU.box() == a_FnFabJ.box());

  // A box of valid fine cells where fine-grid data can be accessed
  ProblemDomain FnProblemDomain(a_CrProblemDomain);
  FnProblemDomain.refine(m_refRatio);
  Box validBox(a_FnFabJU.box());
  validBox &= FnProblemDomain;

  // Some definitions
  const int lenMMax = m_At(0).size(1);
  const int nComp   = a_FnFabJU.nComp();

  // Some temporaries
  CHArray<Real, 2> LSbl(lenMMax, nComp);
  int iStencil;             // Stencil index
  CoordTransform ctStoP;    // Transformation from stencil coordinate sytem
                            // to physical coordinate system
  CoordTransform ctPtoS;    // Transformation from physical coordinate sytem
                            // to stencil coordinate system
  IntVect offset;           // Offset between (0,0) in S and P coordinates
  IntVect fineSmlEnd;       // What to add to a fine cell location transformed
                            // from S to P to get the actual location in Z

  int ierBox = 0;  // Accumulated over the box
  // For each coarse cell
  IVSIterator crCellIt(a_CrCells);
  for (crCellIt.begin(); crCellIt.ok(); ++crCellIt)
    {
      IntVect center = crCellIt();

      // Select a stencil and get indexes and transformation info
      selectStencil(center, a_CrProblemDomain, iStencil, ctStoP, ctPtoS, offset,
                    fineSmlEnd);

      // Pack b = <U>_{j^l} according to the stencil
      int iN = 0;
      auto& ivsStencil = m_ivsStencilMask(iStencil);
      IVSIterator displ(ivsStencil);
      for (displ.begin(); displ.ok(); ++displ)
        {
          IntVect tdispl = ctStoP.transform(displ());
          tdispl += center;
          // make sure the stencil is domain aware
          CH_assert(a_CrProblemDomain.contains(tdispl));
          for (int iComp = 0; iComp != nComp; ++iComp)
            {
              LSbl(iN, iComp) = a_CrFabU(tdispl, iComp);
            }
          ++iN;
        }
      CH_assert(iN <= lenMMax); // make sure matrix is not over filled

      // Get the box of fine cells to update (in stencil space)
      // Start with the box in physical space
      const int iArIdx = a_CrCellsArIdx(center);
      Box fineBoxInCrP = m_FnBoxFillInterp(iArIdx);
      // Get indicies w.r.t. small corner of coarse cell
      const IntVect loP(fineBoxInCrP.smallEnd() - center*m_refRatio);
      const IntVect hiP(fineBoxInCrP.bigEnd() - center*m_refRatio);
      // Transform to stencil space
      const IntVect loPtoS(ctPtoS.transform(loP - offset));
      const IntVect hiPtoS(ctPtoS.transform(hiP - offset));
      // Find the small and big corners in stencil space
#ifdef USE_STCVECTOR
      const IntVect loS(stc::min(loPtoS, hiPtoS));
      const IntVect hiS(stc::max(loPtoS, hiPtoS));
#else
      const IntVect loS(min(loPtoS, hiPtoS));
      const IntVect hiS(max(loPtoS, hiPtoS));
#endif
      Box fineBoxInCrS(loS, hiS);

      // Stencil size
      const int lenM = m_ivsStencilMask(iStencil).numPts();

      // LS Fill U ghosts
      // The old fortran  method
      // int ier = 0;
      // FORT_MAPLSFILLGHOSTS(
      //   CHF_FRA(a_FnFabU),
      //   CHF_FRA(a_FnFabJU),
      //   CHF_CONST_FRA(a_FnFabJ),
      //   CHF_BOX(validBox),
      //   CHF_CONST_I1D(ctStoP.xctm(), 2*SpaceDim),
      //   CHF_BOX(fineBoxInCrS),
      //   CHF_CONST_INTVECT(fineSmlEnd),
      //   CHF_CONST_RCHARRAY(2, LSbl),
      //   CHF_CONST_RCHARRAY(RANK_SPACEDIM_PLUS_1, m_fnXiLSx[iStencil]),
      //   CHF_CONST_RCHARRAY(RANK_SPACEDIM_PLUS_1, m_fnGradXiLSx[iStencil]),
      //   CHF_CONST_INT(lenM),
      //   CHF_INT(ier));
      // if (ier != 0)
      //   {
      //     if (ierBox == 0)
      //       {
      //         ierBox = -1;
      //         MayDay::Warning("Possible failure in least-squares computation "
      //                         "in FORT_MAPLSFILLGHOST");
      //         pout() << "Error code: " << ier << std::endl;
      //       }
      //     ierBox = std::max(ierBox, ier);
      //     // But continue anyways
      //   }

      // c++ version, should work identically
      // mapLSfillGhosts(a_FnFabU,
      //                 a_FnFabJU,
      //                 a_FnFabJ,
      //                 validBox,
      //                 ctStoP,
      //                 fineBoxInCrS,
      //                 fineSmlEnd,
      //                 LSbl,
      //                 m_fnXiLSx[iStencil],
      //                 m_fnGradXiLSx[iStencil],
      //                 lenM);
      // Only fill U for now, JU update after clipping
      mapLSfillGhostsU(a_FnFabU,
                       validBox,
                       ctStoP,
                       fineBoxInCrS,
                       fineSmlEnd,
                       LSbl,
                       m_fnXiLSx[iStencil],
                       m_fnGradXiLSx[iStencil],
                       lenM);
    }  // Loop over coarse cells

  if (m_useClipping)
    {
      // Clip U
      applyClipping(a_FnFabU, a_CrCells, a_CrCellsArIdx, a_CrFabU, a_CrProblemDomain);
    }
  // Solve for JU from U, using gradients of U from coarse interpolation stencil
  for (crCellIt.begin(); crCellIt.ok(); ++crCellIt)
    {
      IntVect center = crCellIt();

      // Select a stencil and get indexes and transformation info
      selectStencil(center, a_CrProblemDomain, iStencil, ctStoP, ctPtoS, offset,
                    fineSmlEnd);

      // Get the box of fine cells to update (in stencil space)
      // Start with the box in physical space
      const int iArIdx = a_CrCellsArIdx(center);
      Box fineBoxInCrP = m_FnBoxFillInterp(iArIdx);
      // Get indicies w.r.t. small corner of coarse cell
      const IntVect loP(fineBoxInCrP.smallEnd() - center*m_refRatio);
      const IntVect hiP(fineBoxInCrP.bigEnd() - center*m_refRatio);
      // Transform to stencil space
      const IntVect loPtoS(ctPtoS.transform(loP - offset));
      const IntVect hiPtoS(ctPtoS.transform(hiP - offset));
      // Find the small and big corners in stencil space
#ifdef USE_STCVECTOR
      const IntVect loS(stc::min(loPtoS, hiPtoS));
      const IntVect hiS(stc::max(loPtoS, hiPtoS));
#else
      const IntVect loS(min(loPtoS, hiPtoS));
      const IntVect hiS(max(loPtoS, hiPtoS));
#endif
      Box fineBoxInCrS(loS, hiS);

      // Stencil size
      const int lenM = m_ivsStencilMask(iStencil).numPts();
      mapLSfillGhostsJU(a_FnFabJU,
                        a_FnFabU,
                        a_FnFabJ,
                        validBox,
                        ctStoP,
                        fineBoxInCrS,
                        fineSmlEnd,
                        LSbl,
                        m_fnXiLSx[iStencil],
                        m_fnGradXiLSx[iStencil],
                        lenM);
    } // Loop over coarse cells
}

/// Clip and redistribute an interpolation to prevent new extrema
void
FourthOrderMappedFineInterp::applyClipping(
  FArrayBox&           a_FnFabU,
  const IntVectSet&    a_CrCells,
  const BaseFab<int>&  a_CrCellsArIdx,
  const FArrayBox&     a_CrFabU,
  const ProblemDomain& a_CrProblemDomain) const
{
  CH_TIME("FourthOrderMappedFineInterp::applyClipping");
  int iStencil;             // Stencil index
  CoordTransform ctStoP;    // Transformation from stencil coordinate sytem
                            // to physical coordinate system
  CoordTransform ctPtoS;    // Transformation from physical coordinate sytem
                            // to stencil coordinate system
  IntVect offset;           // Offset between (0,0) in S and P coordinates
  IntVect fineSmlEnd;       // What to add to a fine cell location transformed
                            // from S to P to get the actual location in Z

  // Coarsened cells plus a radius of 1
  // For smoothing, the clipped values need to exist this far out
  IntVectSet CrCellsP1(a_CrCells);
  CrCellsP1.grow(1);
  CrCellsP1 &= a_CrProblemDomain;
  IVSIterator crCellP1It(CrCellsP1);

  // fine domain
  Box CrBox = a_CrFabU.box();
  Box FnBox = a_FnFabU.box();
  Box CrFnBox = refine(CrBox, m_refRatio);
  ProblemDomain FnDomain = refine(a_CrProblemDomain, m_refRatio);

  // Make clipping bounds
  FArrayBox minBnd(CrFnBox, m_nComp);
  FArrayBox maxBnd(CrFnBox, m_nComp);
  // populate clipping bounds using the least squares stencils over the coarse level
  CH_assert(CrFnBox.contains(FnBox));
  IVSIterator crCellIt(a_CrCells);
  for (crCellIt.begin(); crCellIt.ok(); ++crCellIt)
    {
      const IntVect center = crCellIt();
      // Select a stencil and get indexes and transformation info
      selectStencil(center, a_CrProblemDomain, iStencil, ctStoP, ctPtoS, offset,
                    fineSmlEnd);

      // find the min and max of each stencil
      auto& ivsStencil = m_ivsStencilMask(iStencil);
      IVSIterator displ(ivsStencil);
      for (int comp = 0; comp != m_nComp; comp++)
        {
          Real minVal = std::numeric_limits<Real>::max();
          Real maxVal = std::numeric_limits<Real>::lowest();
          for (displ.begin(); displ.ok(); ++displ)
            {
              const IntVect tdispl = ctStoP.transform(displ()) + center;
              // make sure the stencil is domain aware
              CH_assert(a_CrProblemDomain.contains(tdispl));
              // find the min and max
              minVal = std::min(minVal, a_CrFabU[MD_IV(tdispl, comp)]);
              maxVal = std::max(maxVal, a_CrFabU[MD_IV(tdispl, comp)]);
            }

          // Get the box of fine cells to update (in stencil space)
          // Start with the box in physical space
          const int iArIdx = a_CrCellsArIdx(center);
          const Box fineBoxInCrP = m_FnBoxFillInterp(iArIdx);
          MD_BOXLOOP(fineBoxInCrP, i)
            {
              minBnd[MD_IX(i, comp)] = minVal;
              maxBnd[MD_IX(i, comp)] = maxVal;
            }
        }
    }

  // Smooth bounds
  // use refinement ratio number of iterations, since that requires 1 extra
  // coarse ghost cell
  // applySmoothing(minBnd, FnBox, FnDomain, m_refRatio[0]);
  // applySmoothing(maxBnd, FnBox, FnDomain, m_refRatio[0]);

  // Compute second order derivatives for smoothness check
  // Leave undivided because only relative values matter
  FArrayBox dUdXi2;
  if (m_useHOClipping)
    {
      dUdXi2.define(CrBox, m_nComp*SpaceDim);
      //MD_BOXLOOP(CrBoxP1, j)
      for (crCellIt.begin(); crCellP1It.ok(); ++crCellP1It)
        {
          //IntVect iv = MD_GETIV(j);
          const IntVect jv = crCellP1It();
          IntVect iv = jv;
          for (const auto dir : EachDir)
            {
              const IntVect p1 = IntVect_basis(dir);
              if (!a_CrProblemDomain.domainBox().contains(iv+p1))
                {
                  iv -= p1;
                }
              if (!a_CrProblemDomain.domainBox().contains(iv-p1))
                {
                  iv += p1;
                }
              // 2nd order if iv did not shift, 1st otherwise
              for (int comp = 0; comp != m_nComp; comp++)
                {
                  dUdXi2[MD_IV(jv, comp + m_nComp*dir)] = (a_CrFabU[MD_IV(iv-p1, comp)]
                                                           - 2*a_CrFabU[MD_IV(iv, comp)]
                                                           + a_CrFabU[MD_IV(iv+p1, comp)]);
                }
            }
        }
    }

  // Do the clipping and redistributing
  for (int comp = 0; comp != m_nComp; comp++)
    {
      for (crCellIt.begin(); crCellIt.ok(); ++crCellIt)
        {
          const IntVect center = crCellIt();
          const int iArIdx = a_CrCellsArIdx(center);
          const Box fineBoxInCrP = m_FnBoxFillInterp(iArIdx);
          const bool isFullFineBox = fineBoxInCrP.size() == IntVect::Unit*m_refRatio;
          // A smoothness check
          if (m_useHOClipping)
            {
              // Check if all second derivatives have the same sign
              bool smooth = true;
              for (const auto dir : EachDir)
                {
                  const IntVect p1 = IntVect_basis(dir);
                  for (const auto deriv : EachDir)
                    {
                      const int dComp = comp + deriv*m_nComp;
                      if (a_CrProblemDomain.domainBox().contains(center+p1))
                        {
                          smooth &= (std::signbit(dUdXi2[MD_IV(center, dComp)]) ==
                                     std::signbit(dUdXi2[MD_IV(center + p1, dComp)]));
                        }
                      if (a_CrProblemDomain.domainBox().contains(center-p1))
                        {
                          smooth &= (std::signbit(dUdXi2[MD_IV(center, dComp)]) ==
                                     std::signbit(dUdXi2[MD_IV(center - p1, dComp)]));
                        }
                      if (a_CrProblemDomain.domainBox().contains(center+2*p1))
                        {
                          smooth &= (std::signbit(dUdXi2[MD_IV(center, dComp)]) ==
                                     std::signbit(dUdXi2[MD_IV(center + 2*p1, dComp)]));
                        }
                      if (a_CrProblemDomain.domainBox().contains(center-2*p1))
                        {
                          smooth &= (std::signbit(dUdXi2[MD_IV(center, dComp)]) ==
                                     std::signbit(dUdXi2[MD_IV(center - 2*p1, dComp)]));
                        }
                    }
                }
              // skip clipping for smooth regions
              if (smooth) continue;
            }
          
          // Clip the overshoots and store the excess
          Real extraMass = 0;
          MD_BOXLOOP(fineBoxInCrP, i)
            {
              const Real minBndDif = (a_FnFabU[MD_IX(i, comp)] -
                                      minBnd[MD_IX(i, comp)]);
              const Real maxBndDif = (a_FnFabU[MD_IX(i, comp)] -
                                      maxBnd[MD_IX(i, comp)]);
              if (minBndDif < 0)
                {
                  extraMass += minBndDif;
                  a_FnFabU[MD_IX(i, comp)] = minBnd[MD_IX(i, comp)];
                }
              else if (maxBndDif > 0)
                {
                  extraMass += maxBndDif;
                  a_FnFabU[MD_IX(i, comp)] = maxBnd[MD_IX(i, comp)];
                }
            }
          if (extraMass == 0) continue; // there is no clipping to be done
          // Solve redistribution allowances
          Real sumAllowance = 0;
          FArrayBox allowance(fineBoxInCrP, 1);
          if (extraMass > 0)
            {
              MD_BOXLOOP(fineBoxInCrP, i)
                {
                  allowance[MD_IX(i, 0)] = std::max(maxBnd[MD_IX(i, comp)]
                                                    - a_FnFabU[MD_IX(i, comp)],
                                                    0.);
                  sumAllowance += allowance[MD_IX(i, 0)];
                }
            }
          else if (extraMass < 0)
            {
              MD_BOXLOOP(fineBoxInCrP, i)
                {
                  allowance[MD_IX(i, 0)] = std::min(minBnd[MD_IX(i, comp)]
                                                    - a_FnFabU[MD_IX(i, comp)],
                                                    0.);
                  sumAllowance += allowance[MD_IX(i, 0)];
                }
            }
          // Make sure redistribution is possible
          if (sumAllowance == 0) continue;
          
          // Partial boxes can redistribute outside the existing solution
          if (!isFullFineBox)
            {
              // pout() << "Pushing " << sumAllowance - extraMass
              //        << " outside of solved ghost cells" << std::endl;
              extraMass = sumAllowance;
            }
          // Check if there is grossly more mass to distribute than possible
          // The bounds are loose, so it only asserts for catastrophic failure
      // CH_assert((std::abs(extraMass) - std::abs(sumAllowance)) <
      //           std::abs(sumAllowance)*std::sqrt(std::numeric_limits<Real>::epsilon()));

          // Do the redistribution
          MD_BOXLOOP(fineBoxInCrP, i)
            {
              a_FnFabU[MD_IX(i, comp)] += extraMass * allowance[MD_IX(i, 0)] / sumAllowance;
            }
        }
    }
}

/// Clip and redistribute an interpolation to prevent new extrema
void
FourthOrderMappedFineInterp::applyMappedClipping(
  FArrayBox&           a_FnFabJU,
  const FArrayBox&     a_FnFabJ,
  const Box&           a_CrBox,
  const FArrayBox&     a_CrFabU,
  const FArrayBox&     a_CrFabJU,
  const ProblemDomain& a_CrProblemDomain) const
{
  CH_TIME("FourthOrderMappedFineInterp::applyMappedClipping");
  int iStencil;             // Stencil index
  CoordTransform ctStoP;    // Transformation from stencil coordinate sytem
                            // to physical coordinate system
  CoordTransform ctPtoS;    // Transformation from physical coordinate sytem
                            // to stencil coordinate system
  IntVect offset;           // Offset between (0,0) in S and P coordinates
  IntVect fineSmlEnd;       // What to add to a fine cell location transformed
                            // from S to P to get the actual location in Z

  // coarse grown by 1 cell
  Box CrBoxP1(a_CrBox);
  CrBoxP1.grow(1);
  CrBoxP1 &= a_CrProblemDomain;
  Box CrBoxP2(a_CrBox);
  CrBoxP2.grow(2);
  CrBoxP2 &= a_CrProblemDomain;

  // fine domain
  Box FnBox = refine(a_CrBox, m_refRatio);
  Box FnBoxP1 = refine(CrBoxP1, m_refRatio);
  ProblemDomain FnDomain = refine(a_CrProblemDomain, m_refRatio);

  // Make clipping bounds
  FArrayBox minBnd(FnBoxP1, m_nComp);
  FArrayBox maxBnd(FnBoxP1, m_nComp);
  // populate clipping bounds using the least squares stencils over the coarse level
  CH_assert(a_CrFabU.box().contains(CrBoxP1));
  for (int comp = 0; comp != m_nComp; comp++)
    {
      MD_BOXLOOP(CrBoxP1, j)
        {
          const IntVect center = MD_GETIV(j);
          // Select a stencil and get indexes and transformation info
          selectStencil(center, a_CrProblemDomain, iStencil, ctStoP, ctPtoS, offset,
                        fineSmlEnd);

          // Get the refinement ratio in the stencil coordinate system
          const IntVect refRatioS = absolute(ctPtoS.transform(m_refRatio));

          // find the min and max of each stencil
          auto& ivsStencil = m_ivsStencilMask(iStencil);
          IVSIterator displ(ivsStencil);

          Real minVal = std::numeric_limits<Real>::max();
          Real maxVal = std::numeric_limits<Real>::lowest();
          for (displ.begin(); displ.ok(); ++displ)
            {
              const IntVect tdispl = ctStoP.transform(displ()) + center;
              // make sure the stencil is domain aware
              CH_assert(a_CrProblemDomain.contains(tdispl));
              // find the min and max
              minVal = std::min(minVal, a_CrFabU[MD_IV(tdispl, comp)]);
              maxVal = std::max(maxVal, a_CrFabU[MD_IV(tdispl, comp)]);
            }

          // Get the box of fine cells to update
          const Box fineBoxInCrc(center*m_refRatio,
                                 (center+1)*m_refRatio - IntVect::Unit);
          MD_BOXLOOP(fineBoxInCrc, i)
            {
              minBnd[MD_IX(i, comp)] = minVal;
              maxBnd[MD_IX(i, comp)] = maxVal;
            }
        }
    }

  // Smooth bounds
  // use refinement ratio number of iterations, since that requires 1 extra
  // coarse ghost cell
  applySmoothing(minBnd, FnBox, FnDomain, m_refRatio[0]);
  applySmoothing(maxBnd, FnBox, FnDomain, m_refRatio[0]);

  // Compute second order derivatives for smoothness check
  // Leave undivided because only relative values matter
  FArrayBox dUdXi2;
  if (m_useHOClipping)
    {
      dUdXi2.define(CrBoxP2, m_nComp*SpaceDim);
      MD_BOXLOOP(CrBoxP1, j)
        {
          IntVect iv = MD_GETIV(j);
          for (const auto dir : EachDir)
            {
              const IntVect p1 = IntVect_basis(dir);
              if (!a_CrProblemDomain.domainBox().contains(iv+p1))
                {
                  iv -= p1;
                }
              if (!a_CrProblemDomain.domainBox().contains(iv-p1))
                {
                  iv += p1;
                }
              // 2nd order if iv did not shift, 1st otherwise
              for (int comp = 0; comp != m_nComp; comp++)
                {
                  dUdXi2[MD_IX(j, comp + m_nComp*dir)] = (a_CrFabU[MD_IV(iv-p1, comp)]
                                                          - 2*a_CrFabU[MD_IV(iv, comp)]
                                                          + a_CrFabU[MD_IV(iv+p1, comp)]);
                }
            }
        }
    }

  // Get second order U
  FArrayBox U2nd(a_FnFabJU.box(), m_nComp);
  for (int comp = 0; comp!=m_nComp; comp++)
    {
      MD_BOXLOOP(U2nd.box(), j)
        {
          U2nd[MD_IX(j, comp)] = a_FnFabJU[MD_IX(j, comp)]/a_FnFabJ[MD_IX(j, 0)];
        }
    }
  
  // Do the clipping and redistributing
  for (int comp = 0; comp != m_nComp; comp++)
    {
      MD_BOXLOOP(a_CrBox, j)
        {
          // get the box of fine cells to redistribute over
          const IntVect center = MD_GETIV(j);
          const Box fineBoxInCrc(center*m_refRatio,
                                 (center+1)*m_refRatio - IntVect::Unit);
          // A smoothness check
          if (m_useHOClipping)
            {
              // Check if all second derivatives have the same sign
              bool smooth = true;
              for (const auto dir : EachDir)
                {
                  const IntVect p1 = IntVect_basis(dir);
                  for (const auto deriv : EachDir)
                    {
                      const int dComp = comp + deriv*m_nComp;
                      if (a_CrProblemDomain.domainBox().contains(center+p1))
                        {
                          smooth &= (std::signbit(dUdXi2[MD_IX(j, dComp)]) ==
                                     std::signbit(dUdXi2[MD_OFFSETIV(j, +, p1, dComp)]));
                        }
                      if (a_CrProblemDomain.domainBox().contains(center-p1))
                        {
                          smooth &= (std::signbit(dUdXi2[MD_IX(j, dComp)]) ==
                                     std::signbit(dUdXi2[MD_OFFSETIV(j, -, p1, dComp)]));
                        }
                      if (a_CrProblemDomain.domainBox().contains(center+2*p1))
                        {
                          smooth &= (std::signbit(dUdXi2[MD_IX(j, dComp)]) ==
                                     std::signbit(dUdXi2[MD_OFFSETIV(j, +, 2*p1, dComp)]));
                        }
                      if (a_CrProblemDomain.domainBox().contains(center-2*p1))
                        {
                          smooth &= (std::signbit(dUdXi2[MD_IX(j, dComp)]) ==
                                     std::signbit(dUdXi2[MD_OFFSETIV(j, -, 2*p1, dComp)]));
                        }
                    }
                }
              // skip clipping for smooth regions
              if (smooth) continue;
            }

          // Clip the overshoots and store the excess
          Real extraMass = 0;
          MD_BOXLOOP(fineBoxInCrc, i)
            {
              const Real minBndDif = (U2nd[MD_IX(i, comp)] -
                                      minBnd[MD_IX(i, comp)]);
              const Real maxBndDif = (U2nd[MD_IX(i, comp)] -
                                      maxBnd[MD_IX(i, comp)]);

              if (minBndDif < 0)
                {
                  extraMass += minBndDif * a_FnFabJ[MD_IX(i, 0)];
                  U2nd[MD_IX(i, comp)] = minBnd[MD_IX(i, comp)];
                  a_FnFabJU[MD_IX(i, comp)] = (minBnd[MD_IX(i, comp)]
                                               * a_FnFabJ[MD_IX(i, 0)]);
                }
              else if (maxBndDif > 0)
                {
                  extraMass += maxBndDif * a_FnFabJ[MD_IX(i, 0)];
                  U2nd[MD_IX(i, comp)] = maxBnd[MD_IX(i, comp)];
                  a_FnFabJU[MD_IX(i, comp)] = (maxBnd[MD_IX(i, comp)]
                                               * a_FnFabJ[MD_IX(i, 0)]);
                }
            }
          if (extraMass == 0) continue; // there is no clipping to be done
          // Solve redistribution allowances
          Real sumAllowance = 0;
          FArrayBox allowance(fineBoxInCrc, 1);
          if (extraMass > 0)
            {
              MD_BOXLOOP(fineBoxInCrc, i)
                {
                  allowance[MD_IX(i, 0)] = std::max(maxBnd[MD_IX(i, comp)]
                                                    - U2nd[MD_IX(i, comp)], 0.);
                  sumAllowance += allowance[MD_IX(i, 0)];
                }
            }
          else if (extraMass < 0)
            {
              MD_BOXLOOP(fineBoxInCrc, i)
                {
                  allowance[MD_IX(i, 0)] = std::min(minBnd[MD_IX(i, comp)]
                                                    - U2nd[MD_IX(i, comp)], 0.);
                  sumAllowance += allowance[MD_IX(i, 0)];
                }
            }

          // Make sure redistribution is possible
          if (sumAllowance == 0) continue;
          // Do the redistribution
          MD_BOXLOOP(fineBoxInCrc, i)
            {
              a_FnFabJU[MD_IX(i, comp)] += extraMass * allowance[MD_IX(i, 0)] / sumAllowance;
            }
        }
    }

  // post smoothing
  if (m_usePostSmoothing)
    {
      applySmoothing(a_FnFabJU, FnBox, FnDomain, m_refRatio[0]);
    }
}


void
FourthOrderMappedFineInterp::applySmoothing(FArrayBox&           a_U,
                                            const Box&           a_validBox,
                                            const ProblemDomain& a_problemDomain,
                                            const int            a_iterations) const
{
  CH_TIME("FourthOrderMappedFineInterp::applySmoothing");
  // interval
  const Interval intv(0, a_U.nComp());
  // Temp working space
  FArrayBox intermediate(a_U.box(), a_U.nComp());
  if (a_iterations > 1)
    {
      intermediate.copy(a_U);
    }
  // smoothing data to operate on
  // these are aliases
  FArrayBox raw;
  FArrayBox smoothed;
  
  // Each smoothing pass
  for (int iter = 0; iter != a_iterations; iter++)
    {
      // Swap smoothing in/out locations, no data movement
      if (iter%2 == 0)
        {
          raw.define(intv, a_U);
          smoothed.define(intv, intermediate);
        }
      else
        {
          raw.define(intv, intermediate);
          smoothed.define(intv, a_U);
        }
      // do a smoothing pass
      for (int comp = 0; comp != a_U.nComp(); comp++)
        {
          MD_BOXLOOP(a_validBox, i)
            {
              IntVect iv = MD_GETIV(i);
              smoothed[MD_IX(i, comp)] = 0;
              for (const auto dir : EachDir)
                {
                  // center contribution
                  smoothed[MD_IX(i, comp)] += 2*raw[MD_IX(i, comp)];
                  // contribution from neighbors
                  const int MD_ID(ii, dir);
                  const IntVect p1 = IntVect_basis(dir);
                  bool hasLo = a_problemDomain.contains(iv-p1);
                  bool hasHi = a_problemDomain.contains(iv+p1);
                  if (hasLo && hasHi)
                    {
                      smoothed[MD_IX(i, comp)] += 1*raw[MD_OFFSETIX(i,+,ii,comp)];
                      smoothed[MD_IX(i, comp)] += 1*raw[MD_OFFSETIX(i,-,ii,comp)];
                    }
                  else if (hasHi)
                    {
                      smoothed[MD_IX(i, comp)] += 2*raw[MD_OFFSETIX(i,+,ii,comp)];
                    }
                  else if (hasLo)
                    {
                      smoothed[MD_IX(i, comp)] += 2*raw[MD_OFFSETIX(i,-,ii,comp)];
                    }
                }
              // Divide my number of additions to get average
              smoothed[MD_IX(i, comp)] /= 4*SpaceDim;
            }
        }
    }
  // swap the output if needed, does move data
  if (a_iterations%2 == 1)
    {
      a_U.copy(intermediate);
    }
}

// ----------------------------------------------------------------------
//  Interpolate the solution in the fine cells within a single coarse
//  cell
//  JUFine   <=  solution on the fine cells
//  JFine     => metric terms on the fine mesh
//  xctm      => coordinate transformation matrix
//  fineBoxInCr
//            => box of fine cells in each coarse cell
//               (Zero:refRatio-Unit)
//  fineSmlEnd
//            => where the lower corner of fineBoxInCr is located for
//               the fine FArrays
//  LSbl      => 'b' for the minimization of Ax=b (0:lenN-1,0:nComp-1)
//  LSd       => 'd' for the constraint Bx=d (0:nComp-1)
//  At        => transpose of 'A' for the minimization of Ax=b.  Size
//               (0:lenN-1,0:lenMMax-1) but filled to
//               (0:lenN-1,0:lenM-1)
//  Xip       => <xi_i^p> for the fine cells
//                 (p[0:lenN-1], cell[fineBox])
//  gradXip   => gradients (note that dir and p are together in 1
//               dimensions to support 6 space dimensions)
//                 (dir[0:SpaceDim-1] (,/*) p[0:lenN-1],
//                  cell[fineBox])
//  lenM      => number of cells in stencil (At = lenN x lenM)
//  stat     <=  return status
//               -1 - mismatch between Lapack computation and internal
//                    computation
//                0 - success
//               >0 - Lapack error
// ----------------------------------------------------------------------

void
FourthOrderMappedFineInterp::mapLSinterp(FArrayBox& a_JUFine,
                                         const FArrayBox& a_JFine,
                                         const Box& a_validBox,
                                         const CoordTransform& a_crStoP,
                                         const Box& a_fineBoxInCr,
                                         const IntVect a_fineSmlEnd,
                                         const CHArray<Real, 2>& a_LSbl,
                                         const Vector<Real>& a_LSd,
                                         const CHMatrix& a_At,
                                         const CHArrayRD& a_Xip,
                                         const CHArrayRD& a_gradXip,
                                         const int a_lenM,
                                         int a_stat) const
{
//---- Declarations
  // No errors yet
  a_stat = 0;

  // Error
  int ier;

  // All arrays are defined based on matrix At which has size
  // lenN x lenMMax.  Components are taken from LSd

  // Matrix sizes
  int lenN    = a_At.size(0);
//  int lenMMax = a_At.size(1);
  int lenNmP  = lenN - 1;
  
// //    Matrix At is lenN x lenM but stored as lenN x lenMMax
//       int lenN, lenMMax, lenNmP, nComp
// //    Number of fine cells in a coarse cell
//       int numFn

  // Gradients of J in a fine cell (dir, fine cell)  The gradients of J are
  // stored if we are filling ghost cells, in which case <U> is interpolated
  // and <JU> is computed via a product average
  FArrayBox gradJ(a_fineBoxInCr, SpaceDim);

  // Product of <J Xip> (power, fine cell)
  FArrayBox JXip(a_fineBoxInCr, lenN);

  // B = Sum <J Xip> (power)
  CHVector LSBu(lenN);
  // Q
  CHMatrix Q(lenN, lenN);
  // Q2tAt
  CHMatrix Q2tAt(lenNmP, lenN);
  // Q2tAtA
  CHMatrix Q2tAtA(lenNmP, lenN);
  // Q2tAtAQ2
  CHMatrix M1(lenNmP, lenNmP);
  // Q2tAtAQ1
  CHVector vQ(lenNmP);
  // vb = Q2tAtb - Q2tAtAQ1y
  CHVector vb(lenNmP);
  // Qtx
  CHVector vQtx(lenN);
  // x
  CHVector vLSx(lenN);

#ifdef LSE_USE_LAPACK
  int lwork;
  CHMatrix A_x(N, M);
  CHVector LSBu_x(M);
  CHVector LSbl_x(N);
  CHVector LSd_x(10); //?
  CHMatrix LSx_x(lenN, m_nComp);
//  Real, allocatable, dimension(:) :: work;
//    Floating point comparison
  //bool FPCompare1;
#endif

//---- Constants
  Real factor = 1./12.;

// DBG
//     write(*,*) 'Vals, ', lenN, lenM, lenMMax, iAtlo0, iAthi0,
//    &     iAtlo1, iAthi1
//     write(*,*) 'Working with iFineSmall = ', fineSmlEnd(0),
//    &     fineSmlEnd(1)
// DBG-end

  //  numFn = CHF_DTERM[
  // &            (ifineBoxInCrhi0 - ifineBoxInCrlo0 + 1);
  // &           *(ifineBoxInCrhi1 - ifineBoxInCrlo1 + 1);
  // &           *(ifineBoxInCrhi2 - ifineBoxInCrlo2 + 1);
  // &           *(ifineBoxInCrhi3 - ifineBoxInCrlo3 + 1);
  // &           *(ifineBoxInCrhi4 - ifineBoxInCrlo4 + 1);
  // &           *(ifineBoxInCrhi5 - ifineBoxInCrlo5 + 1)]

// //    Bounds of JFine in an array, for one sided gradients
//       CHF_DTERM[
//          JFnLoa(0) = CHF_LBOUND[validBox;0];
//          JFnLoa(1) = CHF_LBOUND[validBox;1];
//          JFnLoa(2) = CHF_LBOUND[validBox;2];
//          JFnLoa(3) = CHF_LBOUND[validBox;3];
//          JFnLoa(4) = CHF_LBOUND[validBox;4];
//          JFnLoa(5) = CHF_LBOUND[validBox;5]]
//       CHF_DTERM[
//          JFnHia(0) = CHF_UBOUND[validBox;0];
//          JFnHia(1) = CHF_UBOUND[validBox;1];
//          JFnHia(2) = CHF_UBOUND[validBox;2];
//          JFnHia(3) = CHF_UBOUND[validBox;3];
//          JFnHia(4) = CHF_UBOUND[validBox;4];
//          JFnHia(5) = CHF_UBOUND[validBox;5]]

// ---- Compute <JXi^p> for the fine cells
  for (int p = 0; p != lenN; p++)
    {
      LSBu(p) = 0;
    }

  // Loop over the fine cells in a coarse cell
  MD_BOXLOOP(a_fineBoxInCr, ifb)
    {
      const IntVect ifbS = MD_GETIV(ifb);
      // Transform and add offset.  'ifn' then indexes fine FArrayBox
      const IntVect ifbP = a_crStoP.transform(ifbS);
      const IntVect ifn = a_fineSmlEnd + ifbP;

      // Loop over the space dimensions to get the gradients.  We
      // loop over the stencil directions and convert this to
      // physical directions.
      for (int d = 0; d !=SpaceDim; d++)
        {
          const IntVect i1S = IntVect_basis(d); //CHF_ID(0,d);
          const IntVect i1P = a_crStoP.transform(i1S);

          // Compute the gradients of J in gradJ(d)
          // Compute the gradients of J for this direction using first-order
          // gradient orthogonal to coarse-fine interface (if J doesn't exist
          // at +2 ghost outside fine grid)
          gradJ[MD_IX(ifb, d)] = 0.;
          int gradOffsetCount = 0;
          const int dP = a_crStoP.indexAt(d);
          const int dSign = a_crStoP.signAt(d);
          if ((ifn(dP) + dSign - a_validBox.smallEnd(dP))*
              (ifn(dP) + dSign - a_validBox.bigEnd(dP)) < 0)
            {
              // The selection of JFine using i1P is always on the high side
              // according to the gradient in the stencil direction.  The
              // other side must be subtracted which is why gradOffsetCount is
              // decremented
              gradOffsetCount -= 1;
              gradJ[MD_IX(ifb, d)] += a_JFine[MD_IV(ifn + i1P, 0)];
            }
          if ((ifn(dP) - dSign - a_validBox.smallEnd(dP))*
              (ifn(dP) - dSign - a_validBox.bigEnd(dP)) < 0)
            {
              // The selection of JFine using i1P is always on the low side
              // according to the gradient in the stencil direction.  The
              // other side must be added which is why gradOffsetCount is
              // incremented
              gradOffsetCount += 1;
              gradJ[MD_IX(ifb, d)] += a_JFine[MD_IV(ifn - i1P, 0)];
            }
          if (gradOffsetCount == 0)
            {
              // If zero, this is a centered difference
              gradJ[MD_IX(ifb, d)] /= 2.;
            }
          else
            {
              // Otherwise this is a one-sided difference
              gradJ[MD_IX(ifb, d)] += gradOffsetCount*a_JFine[MD_IV(ifn, 0)];
            }
        }

      // For each power (this is completely in stencil directions)
      for (int p = 0; p != lenN; p++)
        {
          // Compute the products <J Xip>_{i^l+1} in JXip(p, ifb)
          Real JXipGradSum = 0.;
          for (int d = 0; d != SpaceDim; d++)
            {
              JXipGradSum += gradJ[MD_IX(ifb, d)]*a_gradXip(p + d*lenN, ifbS);
            }
          JXip[MD_IX(ifb, p)] = a_JFine[MD_IV(ifn, 0)]*a_Xip(p, ifbS) +
            factor*JXipGradSum;

          // Add to the vector B
          LSBu(p) += JXip[MD_IX(ifb, p)];
        }
    } // Loop over fine cells
  // Weight the averages by volume.
  for (int p = 0; p != lenN; p++)
    {
      LSBu(p) /= stc::product(m_refRatio);
    }

//DBG
#ifdef LSE_USE_LAPACK

//---- See what LAPACK can do

  // Get work size
  Real rtmp;
  LAPACK_GGLSE(a_lenM, lenN, 1, A_x, lenMMax, LSBu_x, 1,
               LSbl_x, LSd_x, LSx_x, rtmp, -1, ier);
  if (ier != 0.)
    {
      a_stat = ier;
      return;
    }
  lwork = rtmp1;
  allocate(work(lwork));

  for (int iComp = 0; iComp != m_nComp; iComp++)
    {
      for (int ir = 0; ir != lenN; ir++)
        {
          LSBu_x(ir) = LSBu(ir);
        }
      for (int ir = 0; ir != a_lenM; ir++)
        {
          LSbl_x(ir) = LSbl(ir, iComp);
          for (int ic = 0; ic != lenN; ic++)
            {
              A_x(ir, ic) = At(ic, ir);
            }
        }
      LSd_x(0) = LSd(iComp);

//DBG
//        write(*,*) 'For component', iComp
//        write(*,*) 'This is A'
//        do ir = 0, lenM-1
//           write(*, "(13(ES9.2, 1X))") (A_x(ir, ic), ic= 0, lenN-1)
//        }
//        write(*,*) 'This is B'
//        write(*, "(10(ES9.2, 1X))") (LSBu_x(ic), ic= 0, lenN-1)
//        write(*,*) 'This is b'
//        write(*, "(13(ES9.2, 1X))") (LSbl_x(ic), ic= 0, lenM-1)
//        write(*,*) 'This is d'
//        write(*, "(ES9.1)") LSd_x(0)
//DBG-end

//        DGELS('N', lenM, lenN, 1, A_x, lenMMax, LSbl_x, lenM,
//               work, lwork, ier)
//        do ic= 0, lenN-1
//           LSx_x(ic, iComp) = LSbl_x(ic)
//        }
      LAPACK_GGLSE(lenM, lenN, 1, A_x, lenMMax, LSBu_x, 1,
                   LSbl_x, LSd_x(0), LSx_x(0, iComp), work, lwork, ier);
      if (ier != 0.)
        {
          a_stat = ier;
          return;
        }
//DBG
//        write(*, "(/A, /10(ES13.6, 1X))") 'LAPACK X:',
//            (LSx_x(ic, iComp), ic= 0, lenN-1)
//DBG-end

//---- If we are generating a solution with GGLSE, use the coefficients x to
//---- interpolate the solution everywhere.

#ifndef LSE_USE_INTERNAL
//        write(*, '(/''For component '', I3)') iComp

      // Loop over the fine cells in a coarse cell
      MD_BOXLOOP(a_fineBoxInCr, ifb)
        {
          // Transform and add offset. 'ifn' then indexes fine FArrayBox.
          const IntVect ifbP = a_crStoP.transform(ifbS);
          const IntVect ifn = a_fineSmlEnd + ifbP;

          a_JUFine[MD_IV(ifn, iComp)] =
            CHdot(lenN,
                  const_cast<Real*>(&LSx_x(0, iComp)), 1,
                  const_cast<Real*>(&JXip(0, CHF_AUTOIX[ifb])), 1);
//DBG
//           if (iComp == 0) then
//           write(*, '(''Interpolated '', ES22.15, '' at cell '', '//
//               '6(I4, :, '',''))')
//               JUFine(CHF_AUTOIX[ifn], iComp), CHF_AUTOIX[ifn]
//           }
//DBG-end
        }
#endif
    }
  deallocate(work);

// End if using LAPACK
#endif

//DBG-end

#ifdef LSE_USE_INTERNAL

//---- Solve the constrained least-square probles for the coefficients
//---- 'x'

//DBG
//     write(*,*) 'This is A'
//     do ir = 0, lenM-1
//        write(*, "(13(ES9.2, 1X))") (At(ic, ir), ic= 0, lenN-1)
//     }
//     write(*,*) 'This is B'
//     write(*, "(10(ES9.2, 1X))") (LSBu(ic), ic= 0, lenN-1)
//DBG-end

  // For the QR factorization of B,
  // LSR is the L_2 norm of B
  Real LSR = BLAS1_NRM2(lenN, LSBu, 1);

  // Q is directly specified from LSBu -- split into Q1_{(p,0)} and Q2_{(p,1:n)}
  Real alpha = LSBu(0) - LSR;
  Real beta = 2*alpha*(alpha - LSBu(0));

  // ir = row of Q and ic = col of Q
  // row 0 and col 0
  Q(0,0) = 1. - 2*std::pow(alpha, 2)/beta;
  for (int ir = 1, ir != lenN-1, ir++)
    {
      Real tmpQ = -2*alpha*LSBu(ir)/beta;
      Q(ir, 0) = tmpQ;
      Q(0, ir) = tmpQ;
    }

  // remaining rows and columns
  for (int ir = 1; ir != lenN; ir++)
    {
      Q(ir, ir) = 1. - 2*LSBu(ir)**2/beta
        for (int ic = ir+1; ic != lenN; ic++)
          {
            Real tmpQ = -2*LSBu(ir)*LSBu(ic)/beta;
            Q(ir, ic) = tmpQ;
            Q(ic, ir) = tmpQ;
          }
    }
//DBG
//     write(*,*) 'This is our Q'
//     do ir = 0, lenN-1
//        write(*, "(10(ES9.2, 1X))") (Q(ir, ic), ic = 0, lenN-1)
//     }

//     ! Check that Q is orthogonal
//     BLAS3_GEMM('T', 'N', lenN, lenN, lenN, 1., Q, lenN, Q,
//                    lenN, 0., Axx, lenN)
//     write(*,*) 'Is Q orthogonal?'
//     do ir = 0, lenN-1
//        write(*, "(10(ES9.2, 1X))") (Axx(ir, ic), ic = 0, lenN-1)
//     }
//DBG-end

  // Q2tAt
  // BLAS3_GEMM('T', 'N', lenNmP, a_lenM, lenN, 1., Q(0,1), lenN,
  //            At, lenN, 0., Q2tAt, lenNmP);
  CHgemm(Q, a_At, Q2tAt, 'T', 'N');
  // Q2tAtA
  // BLAS3_GEMM('N', 'T', lenNmP, lenN, a_lenM, 1., Q2tAt, lenNmP,
  //            At, lenN, 0., Q2tAtA, lenNmP);
  CHgemm(Q2tAt, a_At, Q2tAtA 'N', 'T');
  // M1 = Q2tAtAQ2
  // BLAS3_GEMM('N', 'N', lenNmP, lenNmP, lenN, 1., Q2tAtA,
  //            lenNmP, Q(0,1), lenN, 0., M1, lenNmP);
  CHgemm(Q2tAtA, Q, M1);
  // Compute M1^{-1} using a Cholesky factorization
  CHinverseCholesky(M1, 'A');
  // LAPACK_POTRF('L', lenNmP, M1, lenNmP, ier);
  // if (ier != 0.)
  //   {
  //     a_stat = ier;
  //     return;
  //   }
  // LAPACK_POTRI('L', lenNmP, M1, lenNmP, ier);
  // if (ier != 0.)
  //   {
  //     a_stat = ier;
  //     return;
  //   }
  // // Symmetric BLAS routines seem quite slow -- populate the upper section of
  // // M1 so we can use the general routines.
  // sy2ge(lenNmP, M1, lenNmP);

  // Compute vQ = Q2tAtAQ1
  // BLAS2_GEMV('N', lenNmP, lenN, 1., Q2tAtA, lenNmP, Q, 1,
  //            0., vQ, 1);
  CHgemv(Q2tAtA, Q, vO);

  // Loop over the components
  for (int iComp = 0; iComp != m_nComp; iComp++)
    {
//DBG
//        write(*,*) 'For component', iComp
//        write(*,*) 'This is b'
//        write(*, "(13(ES9.2, 1X))") (LSbl(ic, iComp), ic= 0, lenM-1)
//        write(*,*) 'This is d'
//        write(*, "(ES9.1)") LSd(iComp)
//DBG-end

      Real LSy = LSd(iComp)/LSR;

      // Compute vb = Q2tAtb - Q2tAtAQ1y
      //            = Q2tAtb - vQy
      for (int ir = 0; ir != lenNmP; ir++)
        {
          vb(ir) = vQ(ir);
        }
      // BLAS2_GEMV('N', lenNmP, a_lenM, 1., Q2tAt, lenNmP,
      //            a_LSbl(0, iComp), 1, -LSy, vb, 1);
      CHgemv(Q2tAt, a_LSbl, vb);

      // Qtx (composed of y and least squares solution, z)
      vQtx(0) = LSy;
      // BLAS2_GEMV('N', lenNmP, lenNmP, 1., M1, lenNmP, vb, 1,
      //            0., vQtx(1), 1);
      CHgemv(M1, vb, vQtx(1));

      // Now finally, vLSx = Q.Qtx
      // BLAS2_GEMV('N', lenN, lenN, 1., Q, lenN, vQtx, 1, 0.,
      //            vLSx, 1);
      CHgemv(Q, vQtx, vLsx);
//DBG
//        write(*, "(/A, /10(ES13.6, 1X))") 'Internal X:',
//            (vLSx(ir), ir = 0, lenN-1)
//DBG-end

//DBG
#ifdef LSE_USE_LAPACK

//---- Compare LAPACK with internal if solutions were computed with both
      bool cmptest = false;
      for (int ir = 0; ir != lenN; ir++)
        {
          if (FPCompare1(vLSx(ir), LSx_x(ir, iComp), LSE_COMPARE_TOL))
            {
              cmptest = true;
            }
        }
      if (cmptest)
        {
          a_stat = -1;
          pout() << "LSE solutions by internal and lapack do not agree for component "
                 << iComp << " at lower fine "
                 << fineSmlEnd(d) << std::endl;
          pout() << "A LAPACK " << LSx_x << std::endl;
          pout() << "A internal " << vLSx << std::endl;
        }
#endif
//        write(*, '(/''For component '', I3)') iComp
//DBG-endx

//---- Now use the coefficients x to interpolate the solution, <JU>, everywhere.

      // Loop over the fine cells in a coarse cell
      MD_BOXLOOP(a_fineBoxInCr, ifb)
        {
          // Transform and add offset.  'ifn' then indexes fine FArrayBox.
          const IntVect ifbS = MD_GETIV(ifb);
          const IntVect ifbP = a_crStoP.transform(ifb);
          const IntVect ifn = a_fineSmlEnd + ifbP;
          ctTransformStoP(ifbPa, ifbSa, xctm);
          
          // a_JUFine[MD_IV(ifn, iComp)] =
          //   BLAS1_DOT(lenN, vLSx, 1, JXip[MD_IV(ifb, 0)], 1);
          a_JUFine[MD_IV(ifn, iComp)] =
            CHdot(lenN,
                  const_cast<Real*>(&vLSx), 1,
                  const_cast<Real*>(&JXip[MD_IV(ifb, 0)]), 1);
//DBG
//           write(*, '(''Interpolated '', ES22.15, '' at cell '', '//
//               '6(I4, :, '',''))')
//               JUFine(CHF_AUTOIX[ifn], iComp), CHF_AUTOIX[ifn]
//DBG-end
        }
    }  // Loop over components

// End if using internal
#endif

  return;
}

// ----------------------------------------------------------------------
//  Compute F * (A^T * A)^{-1} * A^T where F is the displacements or
//  gradients of the displacements in the fine cells.  The answers, when
//  multiplied by LSbl, yield the interpolated solution in the fine cells
//  The inverse of (A^T * A) is determined using a Cholesky decomposition
//  All RCHARRAY arrays store information continguously per fine cell
//  fnXiLSx     <=  If the final answer in a fine cell is Xi * LSx * LSb,
//                  this is (Xi * LSx)
//  fnGradXiLSx <=  If the gradients in a fine cell is grad(Xi) * LSx *
//                  LSb, this is (grad(Xi) * LSx).  The 1D component
//                  index has storage (lenN, SpaceDim) following
//                  column ordering
//  At           => N x M matrix At from coarse cells - essentially
//                  the displacements of the coarse cells from the
//                  center coarse cell.
//              <=  Modified (so do not reuse!)
//  fineBoxInCr  => A box of all the fine cells in a coarse cell
//  fnXip        => Average displacements of the fine cell (with constant
//                  modification K)
//  fnGradXip    => Average of gradients of displacements for each fine
//                  cell (gradient removes constant K).
//  idxCenCrCell => The index of the center coarse cell in M (counting
//                  from zero).  The order of M is whatever an IVS
//                  of all the coarse cells produces using an iterator.
//  lenM         => Actual number of coarse cells used to construct
//                  the polynomial
//  stat        <=  return status
//                   0 - success
//                  !0 - Lapack error
// ----------------------------------------------------------------------
void
FourthOrderMappedFineInterp::mapLSAtAinvAt(CHArrayRD& a_fnXiLSx,
                                           CHArrayRD& a_fnGradXiLSx,
                                           CHMatrix& a_At,
                                           const Box& a_fineBoxInCr,
                                           const CHArrayRD& a_fnXip,
                                           const CHArrayRD& a_fnGradXip,
                                           const int a_idxCenCrCell,
                                           const int a_lenM,
                                           int a_stat) const
{
  // // Arguments needed for Blas/Lapack (must be passed as pointers)
  // char transpose = 'T';
  // int stride = 1;
  // Real one = 1;
  // Real zero = 0;
  // All arrays are defined based on matrix At which has size lenN x lenMMax
  // Matrix sizes
  int lenN    = a_At.size(0);
  int lenMMax = a_At.size(1);
  int lenNmP  = lenN - 1;

  // Matrix AtAinv is lenNmP x lenNmP
  CHMatrix AtAinv(lenNmP, lenNmP);
  // Matrix AtAinvAt is lenN x lenM
  CHMatrix AtAinvAt(lenN, lenN);

  // Get (A^T * A)
  // Row of A^T is removed.  This is a result of splitting A via A*Q but here
  // Q is diagonal.  This removes the row of 1s and now size is lenNmP
  // REVIEW :: make sure this is the same
  // BLAS3_SYRK('L', 'N', lenNmP, a_lenM, 1, a_At(1, 0), lenN, 0, AtAinv, lenNmP);
  CHsyrk(a_At, AtAinv);

  // Compute (A^T * A)^{-1} using a Cholesky factorization
  CHinverseCholesky(AtAinv, 'L');

//DBG
//      if (debug == 1) then
//         write(*, '(''This is (A^T A)^{-1}'')')
//         do ir = 0, lenNmP-1
//            write(*, '(10('' '', ES17.10))')
//              (AtAInv(ir, ic), ic = 0, lenNmP-1)
//         end do
//      end if
//DBG-End

  // This routine also populates the upper section of AtAinvAt during the
  // multiplication.  Note that we leave row 0 blank in AtAinvAt.
  // REVIEW :: make sure this is the same
  // BLAS3_SYMM('L', 'L', lenNmP, a_lenM, 1., AtAinv, lenNmP,
  //            a_At(1, 0), lenN, 0., AtAinvAt(1, 0), lenN);
  CHsymm(AtAinv, a_At, AtAinvAt);

  // Sum the rows and subtract from the column associated with the principal
  // cell.  This is really the - A_1 y part of the operation (b - A_1 y) where
  // A1 is just a column of 1s and y is the solution in the principal cell.
  // Note:  y is really x/R^T but Q is diagonal and R is 1 so these
  // transformations are unnecessary
  for (int ir = 1; ir != lenN; ir++)
    {
      Real sumrow = 0;
      for (int ic = 0; ic != a_lenM; ic++)
        {
          sumrow += AtAinvAt(ir, ic);
        }
      AtAinvAt(ir, a_idxCenCrCell) -= sumrow;
    }

  // Add a row to represent the known solution 'y'.  This is equivalent to
  // reconstructing x from Q*[y, z]^T
  // Add in row of 0s
  for (int ic = 0; ic != a_lenM; ic++)
    {
      AtAinvAt(0, ic) = 0.;
    }
  // Set the principal cell location to 1
  AtAinvAt(0, a_idxCenCrCell) = 1.;
//DBG
//      if (debug == 1) then
//         write(*, '(''This is (A^T A)^{-1} A^T'')')
//         do ir = 0, lenN-1
//            write(*, '(13('' '', ES12.5))')
//              (AtAInvAt(ir, ic), ic = 0, a_lenM-1)
//         end do
//      end if
//DBG-End

//DBG
//     AtAInvAt is an N x M matrix.  All rows should sum to zero except row
//     0 which must sum to one.
//      do ir = 0, lenN-1
//         sumrow = 0.
//         do ic = 0, a_lenM-1
//            sumrow = sumrow + AtAinvAt(ir, ic)
//         }
//         write(*, '(''Sum for row '', I2, 1X, ES12.5)'), ir, sumrow
//      }
//DBG-End

  // Pre-multiply with fine cell displacements to set a_fnXiLSx
  MD_BOXLOOP(a_fineBoxInCr, ifb)
    {
      const IntVect ifbV = MD_GETIV(ifb);
      // REVIEW :: make sure this is the same
      // BLAS2_GEMV(&transpose, &lenN, const_cast<int*>&a_lenM, 1,
      //            &AtAinvAt, &lenN,
      //            &a_fnXip(0, ifbV), 1, 0,
      //            &a_fnXiLSx(0, ifbV), 1);

      CHgemv(AtAinvAt.begin(),
             const_cast<Real*>(&a_fnXip(0, ifbV)),
             const_cast<Real*>(&a_fnXiLSx(0, ifbV)),
             'T', lenN, a_lenM);
    }
// DBG
//      if (debug == 1) then
//         write(*, '(''This is a_fnXiLSx'')')
//         CHF_AUTOMULTIDO[fineBoxInCr; ifb]
//            write(*, '(13('' '', ES12.5))')
//              (a_fnXiLSx(ic, CHF_AUTOIX[ifb]), ic = 0, a_lenM-1)
//         CHF_}
//      end if
//DBG-End

  // Pre-multiply with gradients of fine cell displacements to set fnGradXiLSx
  for (int iDir = 0; iDir != SpaceDim; iDir++)
    {
      MD_BOXLOOP(a_fineBoxInCr, ifb)
        {
          const IntVect ifbV = MD_GETIV(ifb);
          // REVIEW :: make sure this is the same
          // BLAS2_GEMV(&transpose, &lenN, &a_lenM, 1,
          //            &AtAinvAt, &lenN,
          //            &a_fnGradXip(lenN*iDir, ifbV), 1, 0,
          //            &a_fnGradXiLSx(lenMMax*iDir, ifbV), 1);

          CHgemv(AtAinvAt.begin(),
                 const_cast<Real*>(&a_fnGradXip(lenN*iDir, ifbV)),
                 const_cast<Real*>(&a_fnGradXiLSx(lenMMax*iDir, ifbV)),
                 'T', lenN, a_lenM);
        }
    }
  return;
}

// ----------------------------------------------------------------------
//  Interpolate the solution in the fine cells within a single coarse
//  cell
//  UFine    <=  solution on the fine cells
//  JUFine   <=  metrics  solution product on the fine cells.  The
//               box defining JUFine can be smaller than UFine and not
//               contain fineBoxInCr (when transformed and offset in
//               physical space).  JUFine is not computed for any cells
//               that are not contained.
//  JFine     => metric terms on the fine mesh
//  validBox  => box of valid cells for JFine and JUFine
//  xctm      => x coordinate transformation matrix
//  fineBoxInCr
//            => box of fine cells in each coarse cell
//               (Zero:refRatio-Unit)
//  fineSmlEnd
//            => where the lower corner of fineBoxInCr is located for
//               the fine FArrays
//  LSbl      => 'b' for the minimization of Ax=b (0:lenN-1,0:nComp-1)
//  At        => transpose of 'A' for the minimization of Ax=b.  Size
//               (0:lenN-1,0:lenMMax-1) but filled to
//               (0:lenN-1,0:lenM-1)
//  Xip       => <xi_i^p> for the fine cells
//                 (p[0:lenN-1], cell[fineBox])
//  gradXip   => gradients (note that dir and p are together in 1
//               dimensions to support 6 space dimensions)
//                 (dir[0:SpaceDim-1] (,/*) p[0:lenN-1],
//                  cell[fineBox])
//  lenM      => number of cells in stencil (At = lenN x lenM)
//  stat     <=  return status
//               -1 - mismatch between Lapack computation and internal
//                    computation
//                0 - success
//               >0 - Lapack error
// ----------------------------------------------------------------------
void
FourthOrderMappedFineInterp::mapLSfillGhosts(FArrayBox& a_UFine,
                                             FArrayBox& a_JUFine,
                                             const FArrayBox& a_JFine,
                                             const Box& a_validBox,
                                             const CoordTransform& a_crStoP,
                                             const Box& a_fineBoxInCr,
                                             const IntVect& a_fineSmlEnd,
                                             const CHArray<Real, 2>& a_LSbl,
                                             const CHArrayRD& a_fnXiLSx,
                                             const CHArrayRD& a_fnGradXiLSx,
                                             const int a_lenM) const
{
  // Matrix sizes
  // Matrix LSbl is a_lenM x nComp but stored as lenMMax x nComp
  int lenMMax = a_LSbl.size(0);
//  int nComp   = a_LSbl.size(1);

  // Constants
  Real factor = 1./12.;

  // Gradients of J in a fine cell (dir)
  std::array<Real, SpaceDim> gradJ;

  // ---- Now use the coefficients x to interpolate <U>, d<U>dxi and compute <JU>
  // ---- for filling fine ghost cells

  // Loop over the fine cells in a coarse cell
  MD_BOXLOOP(a_fineBoxInCr, ifb)
    {
      const IntVect ifbS = MD_GETIV(ifb);
      // Transform and add offset. 'ifn' then indexes fine FArrayBox.
      const IntVect ifbP = a_crStoP.transform(ifbS);
      const IntVect ifn = a_fineSmlEnd + ifbP;

      // Check if this fine cell is in JUFine (if not, only return U)
      bool fnCellInJU = a_JUFine.box().contains(ifn);
      if (fnCellInJU)
        {
          // Get the gradients of J
          for (const auto d : EachDir)
            {
              // Convert stencil directions to physical directions for computing
              // the gradients
              const IntVect i1S = IntVect_basis(d); //CHF_ID(0,d);
              const IntVect i1P = a_crStoP.transform(i1S);
              // Compute the gradients of J for this direction using first-order
              // gradient orthogonal to coarse-fine interface (if J doesn't exist
              // at +2 ghost outside fine grid)
              gradJ[d] = 0.;
              // A count of which sides can be used for a gradient
              int gradOffsetCount = 0;
              // Get the coordinate direction in physical space and the sign of
              // the transformation
              const int dP = a_crStoP.indexAt(d);
              const int dSign = a_crStoP.signAt(d);
              if ((ifn(dP) + dSign - a_validBox.smallEnd(dP))*
                  (ifn(dP) + dSign - a_validBox.bigEnd(dP)) < 0)
                {
                  // The selection of JFine using i1P is always on the high side
                  // according to the gradient in the stencil direction.  The
                  // other side must be subtracted which is why gradOffsetCount is
                  // decremented
                  gradOffsetCount -= 1;
                  gradJ[d] += a_JFine[MD_IV(ifn + i1P, 0)];
                }
              if ((ifn(dP) - dSign - a_validBox.smallEnd(dP))*
                  (ifn(dP) - dSign - a_validBox.bigEnd(dP)) < 0)
                {
                  // The selection of a_JFine using i1P is always on the low side
                  // according to the gradient in the stencil direction.  The
                  // other side must be added which is why gradOffsetCount is
                  // incremented
                  gradOffsetCount += 1;
                  gradJ[d] -= a_JFine[MD_IV(ifn - i1P, 0)];
                }
              if (gradOffsetCount == 0)
                {
                  // If zero, this is a centered difference
                  gradJ[d] /= 2;
                }
              else
                {
                  // Otherwise this is a one-sided difference
                  gradJ[d] += gradOffsetCount*a_JFine[MD_IV(ifn, 0)];
                }
            }
        }

      // Compute <U> for each component and the product of <U> and <J>
      // if 'ifn' is contained in a_JUFine
      for (int iComp = 0; iComp != m_nComp; iComp++)
        {
          // Evaluate <U>
          const Real avgU = CHdot(a_lenM,
                                  &a_fnXiLSx(0, ifbS), 1,
                                  &a_LSbl(0, iComp), 1);
          a_UFine[MD_IV(ifn, iComp)] = avgU;
          // Evaluate <JU>
          if (fnCellInJU)
            {
              // Accumulate products of gradients J and U
              Real JUGradSum = 0.;
              // Add to product of gradU and gradJ
              for (const auto d : EachDir)
                {
                  // solve <dJ_d>*<dU_d>
                  JUGradSum += gradJ[d]*
                    CHdot(a_lenM,
                          &a_fnGradXiLSx(d*lenMMax, ifbS), 1,
                          &a_LSbl(0, iComp), 1);
                }
              a_JUFine[MD_IV(ifn, iComp)] =
                a_JFine[MD_IV(ifn, 0)]*avgU + factor*JUGradSum;
            }
        }
    }
  return;
}

void
FourthOrderMappedFineInterp::mapLSfillGhostsU(FArrayBox& a_UFine,
                                              const Box& a_validBox,
                                              const CoordTransform& a_crStoP,
                                              const Box& a_fineBoxInCr,
                                              const IntVect& a_fineSmlEnd,
                                              const CHArray<Real, 2>& a_LSbl,
                                              const CHArrayRD& a_fnXiLSx,
                                              const CHArrayRD& a_fnGradXiLSx,
                                              const int a_lenM) const
{
  // ---- Now use the coefficients x to interpolate <U>, d<U>dxi and compute <JU>
  // ---- for filling fine ghost cells

  // Loop over the fine cells in a coarse cell
  MD_BOXLOOP(a_fineBoxInCr, ifb)
    {
      const IntVect ifbS = MD_GETIV(ifb);
      // Transform and add offset. 'ifn' then indexes fine FArrayBox.
      const IntVect ifbP = a_crStoP.transform(ifbS);
      const IntVect ifn = a_fineSmlEnd + ifbP;

      // Compute <U> for each component
      for (int iComp = 0; iComp != m_nComp; iComp++)
        {
          // Evaluate <U>
          a_UFine[MD_IV(ifn, iComp)] = CHdot(a_lenM,
                                             &a_fnXiLSx(0, ifbS), 1,
                                             &a_LSbl(0, iComp), 1);
        }
    }
  return;
}

void
FourthOrderMappedFineInterp::mapLSfillGhostsJU(FArrayBox& a_JUFine,
                                               const FArrayBox& a_UFine,
                                               const FArrayBox& a_JFine,
                                               const Box& a_validBox,
                                               const CoordTransform& a_crStoP,
                                               const Box& a_fineBoxInCr,
                                               const IntVect& a_fineSmlEnd,
                                               const CHArray<Real, 2>& a_LSbl,
                                               const CHArrayRD& a_fnXiLSx,
                                               const CHArrayRD& a_fnGradXiLSx,
                                               const int a_lenM) const
{
  // Matrix sizes
  // Matrix LSbl is a_lenM x nComp but stored as lenMMax x nComp
  int lenMMax = a_LSbl.size(0);
//  int nComp   = a_LSbl.size(1);

  // Constants
  Real factor = 1./12.;

  // Gradients of J in a fine cell (dir)
  std::array<Real, SpaceDim> gradJ;

  // ---- Now use the coefficients x to interpolate <U>, d<U>dxi and compute <JU>
  // ---- for filling fine ghost cells

  // Loop over the fine cells in a coarse cell
  MD_BOXLOOP(a_fineBoxInCr, ifb)
    {
      const IntVect ifbS = MD_GETIV(ifb);
      // Transform and add offset. 'ifn' then indexes fine FArrayBox.
      const IntVect ifbP = a_crStoP.transform(ifbS);
      const IntVect ifn = a_fineSmlEnd + ifbP;

      // Check if this fine cell is in JUFine (if not, only return U)
      bool fnCellInJU = a_JUFine.box().contains(ifn);
      if (fnCellInJU)
        {
          // Get the gradients of J
          for (const auto d : EachDir)
            {
              // Convert stencil directions to physical directions for computing
              // the gradients
              const IntVect i1S = IntVect_basis(d); //CHF_ID(0,d);
              const IntVect i1P = a_crStoP.transform(i1S);
              // Compute the gradients of J for this direction using first-order
              // gradient orthogonal to coarse-fine interface (if J doesn't exist
              // at +2 ghost outside fine grid)
              gradJ[d] = 0.;
              // A count of which sides can be used for a gradient
              int gradOffsetCount = 0;
              // Get the coordinate direction in physical space and the sign of
              // the transformation
              const int dP = a_crStoP.indexAt(d);
              const int dSign = a_crStoP.signAt(d);
              if ((ifn(dP) + dSign - a_validBox.smallEnd(dP))*
                  (ifn(dP) + dSign - a_validBox.bigEnd(dP)) < 0)
                {
                  // The selection of JFine using i1P is always on the high side
                  // according to the gradient in the stencil direction.  The
                  // other side must be subtracted which is why gradOffsetCount is
                  // decremented
                  gradOffsetCount -= 1;
                  gradJ[d] += a_JFine[MD_IV(ifn + i1P, 0)];
                }
              if ((ifn(dP) - dSign - a_validBox.smallEnd(dP))*
                  (ifn(dP) - dSign - a_validBox.bigEnd(dP)) < 0)
                {
                  // The selection of a_JFine using i1P is always on the low side
                  // according to the gradient in the stencil direction.  The
                  // other side must be added which is why gradOffsetCount is
                  // incremented
                  gradOffsetCount += 1;
                  gradJ[d] -= a_JFine[MD_IV(ifn - i1P, 0)];
                }
              if (gradOffsetCount == 0)
                {
                  // If zero, this is a centered difference
                  gradJ[d] /= 2;
                }
              else
                {
                  // Otherwise this is a one-sided difference
                  gradJ[d] += gradOffsetCount*a_JFine[MD_IV(ifn, 0)];
                }
            }
        }

      // Compute <U> for each component and the product of <U> and <J>
      // if 'ifn' is contained in a_JUFine
      for (int iComp = 0; iComp != m_nComp; iComp++)
        {
          // Evaluate <JU>
          if (fnCellInJU)
            {
              // Accumulate products of gradients J and U
              Real JUGradSum = 0.;
              // Add to product of gradU and gradJ
              for (const auto d : EachDir)
                {
                  // solve <dJ_d>*<dU_d>
                  JUGradSum += gradJ[d]*
                    CHdot(a_lenM,
                          &a_fnGradXiLSx(d*lenMMax, ifbS), 1,
                          &a_LSbl(0, iComp), 1);
                }
              a_JUFine[MD_IV(ifn, iComp)] =
                a_JFine[MD_IV(ifn, 0)]*a_UFine[MD_IV(ifn, iComp)] + factor*JUGradSum;
            }
        }
    }
  return;
}
/*--------------------------------------------------------------------*/
//   Preset the coarsened-fine \<U\> from coarse level data
/**  \param[in]  a_CrU  \<U\> on the coarse level
 *   \param[in]  a_vectorIntv
 *                      An interval in components of U consisting of
 *                      vector data
 *//*-----------------------------------------------------------------*/

void
FourthOrderMappedFineInterp::presetCrFnLevU(
  const LevelData<FArrayBox>& a_CrLevU,
  const Interval&             a_vectorIntv)
{
  CH_TIME("FourthOrderMappedFineInterp::presetCrFnLevU");
  CH_assert(isGridDefined());
  m_CrFnStatus |= CrFnStatusHaveU;
  // Normal copy.  If multiBlock, it is expected that extra-block ghosts in
  // a_CrU have been filled by means of a multi-block exchange.
  a_CrLevU.copyTo(m_CrFnLevU, m_copier);
}

// Should not be used
// /*--------------------------------------------------------------------*/
// //   Preset the coarsened-fine \<JU\> from coarse level data
// /**  \param[in]  a_CrJU \<JU\> on the coarse level
//  *//*-----------------------------------------------------------------*/

// void
// FourthOrderMappedFineInterp::presetCrFnLevJU(
//   const LevelData<FArrayBox>& a_CrLevJU)
// {
//   CH_assert(isGridDefined());
//   CH_assert(!m_mbcp.isDefined());  // Must not be used with multi-block
//   m_CrFnStatus |= CrFnStatusHaveJU;
//   a_CrLevJU.copyTo(m_CrFnLevJU, m_copier);
// }

/*--------------------------------------------------------------------*/
//   Select a stencil
/**  Selects a stencil and returns a number of useful parameters
 *   describing the stencil and associated transformations
 *   \param[in]  a_center
 *                      Location of coarse cell
 *   \param[in]  a_CrProblemDomain
 *                      The problem domain as specified on the coarse
 *                      level
 *   \param[out] a_iStencil
 *                      Index of the stencil to use
 *   \param[out] a_ctStoP
 *                      Transforms from stencil (S) to physical (P)
 *                      space
 *   \param[out] a_ctPtoS
 *                      Transforms from physical (P) to stencil (S)
 *                      space
 *   \param[out] a_offset
 *                      The offset defines the variation in location
 *                      (0,0) between S and P space.  We need the
 *                      offset since the location of (0,0) for a fine
 *                      box of cells in a coarse cells in stencil
 *                      space is a corner and is not the same in S and
 *                      P spaces (unlike the stencil box for which the
 *                      location of (0,0) is at the center of the box
 *                      and does not change in S and P spaces).
 *   \param[out] a_fineSmlEnd
 *                      When we transform an index from S to P, we
 *                      will need to add this to get the actual
 *                      IntVect in Z.
 *
 *   \note
 *   <ul>
 *     <li> For transforming the location of a fine cell, with (0,0)
 *          normalized to the location of the small end in the coarse
 *          cell, the operations are essentially
 *            P = Ts->p(S) + offset
 *          where
 *            offset[i] = 0            if Tp->s[i] is positive
 *                      = n_ref[i] - 1 if Tp->s[i] is negative
 *          and
 *            S = Tp->s(P - offset)
 *   </ul>
 *//*-----------------------------------------------------------------*/

void
FourthOrderMappedFineInterp::selectStencil(
  const IntVect        a_center,
  const ProblemDomain& a_CrProblemDomain,
  int&                 a_iStencil,
  CoordTransform&      a_ctStoP,
  CoordTransform&      a_ctPtoS,
  IntVect&             a_offset,
  IntVect&             a_fineSmlEnd) const
{
  // Most cells will be at the interior so use that as default
  int soi[SpaceDim] =
  {
    D_DECL6(0, 0, 0, 0, 0, 0)
  };
  int transform[SpaceDim] =
  {
    D_DECL6(1, 2, 3, 4, 5, 6)
  };
  Box offsetBox(a_center - m_stencilExtent*IntVect::Unit,
                a_center + m_stencilExtent*IntVect::Unit);
  if (!a_CrProblemDomain.domainBox().contains(offsetBox))
    {
      for (const auto iDir : EachDir)
        {
          if (!a_CrProblemDomain.isPeriodic(iDir))
            {
              // We assume the box is big enough such that we can only
              // be at one boundary in a given direction.
              // negative offset is inside the domain
              int offsetDirHi =
                offsetBox.bigEnd(iDir) -
                a_CrProblemDomain.domainBox().bigEnd(iDir);
              int offsetDirLo =
                a_CrProblemDomain.domainBox().smallEnd(iDir) -
                offsetBox.smallEnd(iDir);
              // connected boundaries are valid regions for stencils
              if (a_CrProblemDomain.isConnected(iDir, Side::Lo))
                {
                  offsetDirLo = 0;
                }
              if (a_CrProblemDomain.isConnected(iDir, Side::Hi))
                {
                  offsetDirHi = 0;
                }
              // update using the offset
              const int maxOffset = std::max(offsetDirHi, offsetDirLo);
              if (maxOffset > 0)
                {
                  soi[iDir] = maxOffset;
                  if (offsetDirLo > offsetDirHi)
                    {
                      transform[iDir] = -transform[iDir];
                    }
                }
            }
        }
      // This gives us both the stencil index and the coordinate
      // transformation
      sortEachSetDx(soi, transform);
    }
  // Look up the stencil in
  auto stenIdx = m_mapStencilIndex.find(offset2Key(soi));
  CH_assert(stenIdx != m_mapStencilIndex.end());
  // we want the stencil index
  a_iStencil = stenIdx->second;

/*** Now we are in a transformed space.  Any offset o(ii, jj, kk) from      ***
 *** cell u(i, j, k) as prescribed by a stencil is now given by             ***
 *** ct.transform(o)                                                        ***/

  a_ctStoP.expand(transform);           // Transformation from stencil
                                        // coordinate sytem to physical
                                        // coordinate system
  a_ctPtoS = a_ctStoP;                  // Transformation from physical
  a_ctPtoS.reverse();                   // coordinate sytem to stencil
                                        // coordinate system

  // Define the offset
  a_offset = IntVect::Zero;
  for (const auto iDir : EachDir)
    {
      if (a_ctPtoS.signAt(iDir) < 0)
        {
          a_offset[iDir] += m_refRatio[iDir] - 1;
        }
    }

  // When we transform an index from S to P, we will need to add this to
  // get the actual IntVect in Z.
  a_fineSmlEnd = a_center*m_refRatio + a_offset;
}

/*--------------------------------------------------------------------*/
//   Check for proper nesting of the fine level in the coarse level
/**  This is really expensive but does perform a thorough check.  If
 *   ghost cells are to be filled, then defineCFInterface should have
 *   already been called.  Otherwise, only sufficient nesting to
 *   interpolate all interior cells on the fine mesh is guaranteed.
 *   \param[in]  a_CrBoxes
 *                      Layout of the fine boxes.
 *   \param[in]  a_CrFnBoxes
 *                      Layout of the coarsened fine boxes.
 *   \param[in]  a_CrMetricsPtr
 *                      Used to get problem domain on the coarse mesh
 *//*-----------------------------------------------------------------*/

bool
FourthOrderMappedFineInterp::properNesting(
  const DisjointBoxLayout& a_CrBoxes,
  const DisjointBoxLayout& a_CrFnBoxes,
  LevelGridMetrics *const  a_CrMetricsPtr) const
{
  CH_TIME("FourthOrderMappedFineInterp::properNesting");
  IntVect totalRad = m_stencilExtent*IntVect::Unit;
  if (isCFInterfaceDefined())
    {
      totalRad += m_CrInterpRadVec;
    }
  if (m_willFillGhostsWithCrFnJU && !isMultiblock())
    {
      totalRad += IntVect::Unit;
    }

  // A box which will determine whether a given box adjoins a periodic boundary
  Box crPeriodicTestBox;
  bool isPeriodicDomain = false;
  if (!isMultiblock())
    {
      CH_assert(!a_CrMetricsPtr->getCoordSys().isMultiBlock());
      const ProblemDomain& crProblemDomain =
        a_CrMetricsPtr->getCoordSys().problemDomain(0);
      crPeriodicTestBox = crProblemDomain.domainBox();
      if (crProblemDomain.isPeriodic())
        {
          isPeriodicDomain = true;
          for (int iDir = 0; iDir < SpaceDim; ++iDir)
            {
              if (crProblemDomain.isPeriodic(iDir))
                {
                  crPeriodicTestBox.grow(iDir, -1);
                }
            }
        }
    }

  for (DataIterator ditCrFn = a_CrFnBoxes.dataIterator(); ditCrFn.ok();
       ++ditCrFn)
    {
      const Box crFnBox = a_CrFnBoxes[ditCrFn];
      const ProblemDomain& crBlockDomain =
        a_CrMetricsPtr->getCoordSys().problemDomain(crFnBox);
      //**FIXME currently, there is no testing for if there is sufficient
      //**      nesting across multiblock boundaries
      // This is the number of cells in coarsened fine boxes that need to be
      // filled from real cells in the coarse grid
      Box crFnGhostBox = grow(crFnBox, totalRad);
      //**FIXME this really isn't sufficient.  You either need to be adjacent
      //**      to the boundary or m_CrInterpRadVec cells away.  In meshRefine,
      //**      the block factor generally ensures this.
      crFnGhostBox &= crBlockDomain;
      // FIXME a quick fix for MMB by just treating it like a boundary
      if (a_CrMetricsPtr->getCoordSys().isMultiBlock())
        {
          crFnGhostBox &= crBlockDomain.domainBox();
        }

      // Initialize the IVS as the whole crFnGhostBox
      IntVectSet crFnGhostIVS(crFnGhostBox);

      // Now we iterate over all the boxes in the coarse grid and subtract off
      // all the real cells.  All cells in the coarsened-fine box should be
      // removed
      LayoutIterator litCr = a_CrBoxes.layoutIterator();
      for (litCr.begin(); litCr.ok(); ++litCr)
        {
          const Box crRealBox = a_CrBoxes[litCr];

          crFnGhostIVS -= crRealBox;

          // We also need to remove any periodic images of real cells
          if (isPeriodicDomain
              && !crPeriodicTestBox.contains(crFnGhostBox)
              && !crPeriodicTestBox.contains(crRealBox))
            {
              ShiftIterator shiftIt = crBlockDomain.shiftIterator();
              IntVect shiftMult(crBlockDomain.domainBox().size());
              Box shiftedBox(crRealBox);
              for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
                {
                  IntVect shiftVect = shiftMult*shiftIt();
                  shiftedBox.shift(shiftVect);
                  crFnGhostIVS -= shiftedBox;
                  shiftedBox.shift(-shiftVect);
                }
            }
        }  // Second loop over all boxes

      // crFnGhostIVS should now be empty
      if (!crFnGhostIVS.isEmpty())
        {
          MayDay::Warning("Insufficient nesting detected.");
          return false;
        }
    }  // First loop over boxes on the processor
  return true;
}

/*--------------------------------------------------------------------*/
//   Sort sets of dimensions with the same mesh spacing
/**  Allows us to find a unique stencil index and also gives the
 *   transformation required to orient the physical coordinates in
 *   the coordinate system of this stencil
 *   \param[in]  a_a    The array to sort
 *   \param[out] a_a    Sorted from greatest to lowest for each dx
 *   \param[in]  a_b    An array to rearrange along with the sorting of a_a
 *   \param[out] a_b    Rearranged the same as a_a
 *//*-----------------------------------------------------------------*/

void
FourthOrderMappedFineInterp::sortEachSetDx(int *const a_a,
                                           int *const a_b) const
{
  int index[SpaceDim] =
  {
    D_DECL6(0, 1, 2, 3, 4, 5)
  };
  // Sort each unique dx separately
  bool mask[SpaceDim];
  for (int i = 0; i != m_numUniqueDx; ++i)
    {
      for (int j = 0; j != SpaceDim; ++j)
        {
          mask[j] = (m_uniqueDx[j] == i);
        }
      Sort::insertion(SpaceDim,
                      index,
                      Sort::CmpGreaterIndex<int>(a_a),
                      mask);
    }
  // Now rearrange both a and b to the index
  Sort::Move2Array<int, int> mover(a_a, a_b);
  Sort::rearrangeToIndex(SpaceDim, index, mover);
}

/*--------------------------------------------------------------------*/
//   Get the stencil key from a list of offsets
/**  The offsets must be sorted according to 'sortEachSetDx'
 *   \param[in]  a_offsets
 *                      Sorted offsets of the stencil from the
 *                      interior
 *   \return            A unique stencil index
 *//*-----------------------------------------------------------------*/

int
FourthOrderMappedFineInterp::offset2Key(const int *const a_offsets) const
{
  int mult = 1;
  int key = 0;
  for (int i = SpaceDim; i--;)
    {
      CH_assert(a_offsets[i] >= 0);
      CH_assert(a_offsets[i] <= m_stencilExtent);
      key += mult*a_offsets[i];
      mult *= m_stencilExtent + 1;
    }
  return key;
}

/*--------------------------------------------------------------------*/
//   Load <\delta\xi^\mathbf{p}> for all p
/**  We evaluate
 *   \f[
 *     \langle \delta\xi^{\mathbf{p}} \rangle
 *     =
 *     \prod_d \left(\frac
 *       {\displaystyle \int_{\delta\xi_d - \frac{h_d}{2}}
 *                          ^{\delta\xi_d + \frac{h_d}{2}}
 *         (\xi_d)^{p_d} \mathrm{d}\xi_d}
 *       {\displaystyle h_d}\right)
 *     =
 *     \prod_d \left(\frac
 *       {\displaystyle \left[
 *         \left(\delta\xi_d + \frac{h_d}{2}\right)^{p_d + 1} -
 *         \left(\delta\xi_d - \frac{h_d}{2}\right)^{p_d + 1}\right]}
 *       {\displaystyle h_d\left(p_d + 1\right)}\right)
 *   \f]
 *   for all values of \f$\|\mathbf{p}\|_1 \leq \mbox{m\_degree}\f$.
 *   \param[out] a_A    Vector of results
 *   \param[in]  a_dxi  Displacement of \f$\xi\f$ from center cell
 *   \param[in]  a_h    Mesh spacing
 *//*-----------------------------------------------------------------*/

void
FourthOrderMappedFineInterp::loadAvgXipA(
  Real*           a_A,
  const RealVect& a_dxi,
  const RealVect  a_h) const
{
  RealVect mhi;
  RealVect mlo;
  for (int iDir = 0; iDir != SpaceDim; ++iDir)
    {
      mhi[iDir] = a_dxi[iDir] + 0.5*a_h[iDir];
      mlo[iDir] = a_dxi[iDir] - 0.5*a_h[iDir];
    }
  // For direction '#'
  int D_DECL6(p0, p1, p2, p3, p4, p5);
                                      // p# is the power
  RealVect hi;                        // hi# is the integrand evaluated at the
                                      // upper bound (= mhi#^p#)
  RealVect lo;                        // lo# is the integrand evaluated at the
                                      // lower bound (= mlo#^p#)
  RealVect num;                       // num# is the factor of the numerator for
                                      // this and all lower dimensions
                                      // (= Prod[(hi# - lo#)(...)(hi0 - lo0)])
  RealVect den;                       // den# is the factor of the denominator
                                      // for this and all lower dimension
                                      // (= Prod[(p# + 1)(...)(p0 + 1)])
  Real hVol = a_h.product();          // cell volume
  // And the result is in num[SpaceDim-1]/den[SpaceDim-1]
  D_TERM6(for (p0 = 0, hi[0] = mhi[0], lo[0] = mlo[0], num[0] = (hi[0] - lo[0])       , den[0] = hVol  ; p0                          <= m_degree; ++p0, hi[0] *= mhi[0], lo[0] *= mlo[0], num[0] = (hi[0] - lo[0])       , den[0] = (p0 + 1)*hVol  ),
          for (p1 = 0, hi[1] = mhi[1], lo[1] = mlo[1], num[1] = (hi[1] - lo[1])*num[0], den[1] = den[0]; p0 + p1                     <= m_degree; ++p1, hi[1] *= mhi[1], lo[1] *= mlo[1], num[1] = (hi[1] - lo[1])*num[0], den[1] = (p1 + 1)*den[0]),
          for (p2 = 0, hi[2] = mhi[2], lo[2] = mlo[2], num[2] = (hi[2] - lo[2])*num[1], den[2] = den[1]; p0 + p1 + p2                <= m_degree; ++p2, hi[2] *= mhi[2], lo[2] *= mlo[2], num[2] = (hi[2] - lo[2])*num[1], den[2] = (p2 + 1)*den[1]),
          for (p3 = 0, hi[3] = mhi[3], lo[3] = mlo[3], num[3] = (hi[3] - lo[3])*num[2], den[3] = den[2]; p0 + p1 + p2 + p3           <= m_degree; ++p3, hi[3] *= mhi[3], lo[3] *= mlo[3], num[3] = (hi[3] - lo[3])*num[2], den[3] = (p3 + 1)*den[2]),
          for (p4 = 0, hi[4] = mhi[4], lo[4] = mlo[4], num[4] = (hi[4] - lo[4])*num[3], den[4] = den[3]; p0 + p1 + p2 + p3 + p4      <= m_degree; ++p4, hi[4] *= mhi[4], lo[4] *= mlo[4], num[4] = (hi[4] - lo[4])*num[3], den[4] = (p4 + 1)*den[3]),
          for (p5 = 0, hi[5] = mhi[5], lo[5] = mlo[5], num[5] = (hi[5] - lo[5])*num[4], den[5] = den[4]; p0 + p1 + p2 + p3 + p4 + p5 <= m_degree; ++p5, hi[5] *= mhi[5], lo[5] *= mlo[5], num[5] = (hi[5] - lo[5])*num[4], den[5] = (p5 + 1)*den[4]))
    {
      *a_A++ = num[SpaceDim-1]/den[SpaceDim-1];
    }
}

/*--------------------------------------------------------------------*/
//   Load <\delta\xi^\mathbf{p}> for all p
/**  We evaluate
 *   \f[
 *     \langle \delta\xi^{\mathbf{p}} \rangle
 *     =
 *     \prod_d \left(\frac
 *       {\displaystyle \int_{\delta\xi_d - \frac{h_d}{2}}
 *                          ^{\delta\xi_d + \frac{h_d}{2}}
 *         (\xi_d)^{p_d} \mathrm{d}\xi_d}
 *       {\displaystyle h_d}\right)
 *     =
 *     \prod_d \left(\frac
 *       {\displaystyle \left[
 *         \left(\delta\xi_d + \frac{h_d}{2}\right)^{p_d + 1} -
 *         \left(\delta\xi_d - \frac{h_d}{2}\right)^{p_d + 1}\right]}
 *       {\displaystyle h_d\left(p_d + 1\right)}\right)
 *   \f]
 *   for all values of \f$\|\mathbf{p}\|_1 \leq \mbox{m\_degree}\f$.
 *   \param[out] a_A    Vector of results
 *   \param[in]  a_dxi  Displacement of \f$\xi\f$ from center cell
 *   \param[in]  a_h    Mesh spacing
 *//*-----------------------------------------------------------------*/

void
FourthOrderMappedFineInterp::loadAvgXiKpA(
  Real*           a_A,
  const RealVect& a_dxi,
  const RealVect  a_h) const
{
  RealVect mhi;
  RealVect mlo;
  
  for (int iDir = 0; iDir != SpaceDim; ++iDir)
    {
      mhi[iDir] = a_dxi[iDir] + 0.5*a_h[iDir];
      mlo[iDir] = a_dxi[iDir] - 0.5*a_h[iDir];
    }
  // For direction '#'
  int D_DECL6(p0, p1, p2, p3, p4, p5);
                                      // p# is the power
  RealVect hi;                        // hi# is the integrand evaluated at the
                                      // upper bound (= mhi#^p#)
  RealVect lo;                        // lo# is the integrand evaluated at the
                                      // lower bound (= mlo#^p#)
  RealVect num;                       // num# is the factor of the numerator for
                                      // this and all lower dimensions
                                      // (= Prod[(hi# - lo#)(...)(hi0 - lo0)])
  RealVect den;                       // den# is the factor of the denominator
                                      // for this and all lower dimension
                                      // (= Prod[(p# + 1)(...)(p0 + 1)])
  Real hVol = a_h.product();          // cell volume
  // And the result is in num[SpaceDim-1]/den[SpaceDim-1]
  D_TERM6(for (p0 = 0, hi[0] = mhi[0], lo[0] = mlo[0], num[0] = (hi[0] - lo[0])       , den[0] = hVol  ; p0                          <= m_degree; ++p0, hi[0] *= mhi[0], lo[0] *= mlo[0], num[0] = (hi[0] - lo[0] - (p0 && p0%2 == 0)*std::pow(0.5, p0)*(mhi[0] - mlo[0]))       , den[0] = (p0 + 1)*hVol  ),
          for (p1 = 0, hi[1] = mhi[1], lo[1] = mlo[1], num[1] = (hi[1] - lo[1])*num[0], den[1] = den[0]; p0 + p1                     <= m_degree; ++p1, hi[1] *= mhi[1], lo[1] *= mlo[1], num[1] = (hi[1] - lo[1] - (p1 && p1%2 == 0)*std::pow(0.5, p1)*(mhi[1] - mlo[1]))*num[0], den[1] = (p1 + 1)*den[0]),
          for (p2 = 0, hi[2] = mhi[2], lo[2] = mlo[2], num[2] = (hi[2] - lo[2])*num[1], den[2] = den[1]; p0 + p1 + p2                <= m_degree; ++p2, hi[2] *= mhi[2], lo[2] *= mlo[2], num[2] = (hi[2] - lo[2] - (p2 && p2%2 == 0)*std::pow(0.5, p2)*(mhi[2] - mlo[2]))*num[1], den[2] = (p2 + 1)*den[1]),
          for (p3 = 0, hi[3] = mhi[3], lo[3] = mlo[3], num[3] = (hi[3] - lo[3])*num[2], den[3] = den[2]; p0 + p1 + p2 + p3           <= m_degree; ++p3, hi[3] *= mhi[3], lo[3] *= mlo[3], num[3] = (hi[3] - lo[3] - (p3 && p3%2 == 0)*std::pow(0.5, p3)*(mhi[3] - mlo[3]))*num[2], den[3] = (p3 + 1)*den[2]),
          for (p4 = 0, hi[4] = mhi[4], lo[4] = mlo[4], num[4] = (hi[4] - lo[4])*num[3], den[4] = den[3]; p0 + p1 + p2 + p3 + p4      <= m_degree; ++p4, hi[4] *= mhi[4], lo[4] *= mlo[4], num[4] = (hi[4] - lo[4] - (p4 && p4%2 == 0)*std::pow(0.5, p4)*(mhi[4] - mlo[4]))*num[3], den[4] = (p4 + 1)*den[3]),
          for (p5 = 0, hi[5] = mhi[5], lo[5] = mlo[5], num[5] = (hi[5] - lo[5])*num[4], den[5] = den[4]; p0 + p1 + p2 + p3 + p4 + p5 <= m_degree; ++p5, hi[5] *= mhi[5], lo[5] *= mlo[5], num[5] = (hi[5] - lo[5] - (p5 && p5%2 == 0)*std::pow(0.5, p5)*(mhi[5] - mlo[5]))*num[4], den[5] = (p5 + 1)*den[4]))
    {
      *a_A++ = num[SpaceDim-1]/den[SpaceDim-1];
    }
}

/*--------------------------------------------------------------------*/
//   Load gradients of <\delta\xi^\mathbf{p}> for all p
/**  These gradients are the same for the displacements determined
 *   from either loadAvgXip or loadAvgXiKp
 *   We evaluate
 *   \f[
 *     \langle \delta\frac{\mathrm{d}\xi^{\mathbf{p}}}
 *                        {\mathrm{d}\xi^{d}} \rangle
 *   \f]
 *   for all values of \f$\|\mathbf{p}\|_1 \leq \mbox{m\_degree}\f$.
 *   This is the same as loadAvgXipA except that a new coefficient is
 *   added in front (from the deriviate), p is reduced for
 *   direction d, and the constant K is removed.
 *   \param[out] a_A    Vector of results
 *   \param[in]  a_dxi  Displacement of \f$\xi\f$ from center cell
 *   \param[in]  a_h    Mesh spacing
 *//*-----------------------------------------------------------------*/

void
FourthOrderMappedFineInterp::loadAvgGradXip(
  Real*           a_gradXip,
  const RealVect& a_dxi,
  const RealVect  a_h) const
{
  RealVect mhi;
  RealVect mlo;
  
  for (int iDir = 0; iDir != SpaceDim; ++iDir)
    {
      mhi[iDir] = a_dxi[iDir] + 0.5*a_h[iDir];
      mlo[iDir] = a_dxi[iDir] - 0.5*a_h[iDir];
    }
  // For direction '#'
  IntVect p;
  D_TERM6(int& p0 = p[0];,
          int& p1 = p[1];,
          int& p2 = p[2];,
          int& p3 = p[3];,
          int& p4 = p[4];,
          int& p5 = p[5];)            // p# is the power
  RealVect hi;                        // hi# is the integrand evaluated at the
                                      // upper bound (= mhi#^p#)
  RealVect lo;                        // lo# is the integrand evaluated at the
                                      // lower bound (= mlo#^p#)
  RealVect num;                       // num# is the factor of the numerator for
                                      // this and all lower dimensions
                                      // (= Prod[(hi# - lo#)(...)(hi0 - lo0)])
  RealVect den;                       // den# is the factor of the denominator
                                      // for this and all lower dimension
                                      // (= Prod[(p# + 1)(...)(p0 + 1)])
  Real hVol = a_h.product();          // cell volume
  for (int dir = 0; dir != SpaceDim; ++dir)
    {

      // And the result is in p[dir]*num[SpaceDim-1]/den[SpaceDim-1]
      D_TERM6(for (p0 = 0, hi[0] = (dir != 0)*mhi[0] + (dir == 0), lo[0] = (dir != 0)*mlo[0] + (dir == 0), num[0] = (hi[0] - lo[0])       , den[0] = hVol  ; p0                          <= m_degree; ++p0, hi[0] *= mhi[0], lo[0] *= mlo[0], num[0] = (hi[0] - lo[0])       , den[0] = (p0 + (dir != 0))*hVol  ),
              for (p1 = 0, hi[1] = (dir != 1)*mhi[1] + (dir == 1), lo[1] = (dir != 1)*mlo[1] + (dir == 1), num[1] = (hi[1] - lo[1])*num[0], den[1] = den[0]; p0 + p1                     <= m_degree; ++p1, hi[1] *= mhi[1], lo[1] *= mlo[1], num[1] = (hi[1] - lo[1])*num[0], den[1] = (p1 + (dir != 1))*den[0]),
              for (p2 = 0, hi[2] = (dir != 2)*mhi[2] + (dir == 2), lo[2] = (dir != 2)*mlo[2] + (dir == 2), num[2] = (hi[2] - lo[2])*num[1], den[2] = den[1]; p0 + p1 + p2                <= m_degree; ++p2, hi[2] *= mhi[2], lo[2] *= mlo[2], num[2] = (hi[2] - lo[2])*num[1], den[2] = (p2 + (dir != 2))*den[1]),
              for (p3 = 0, hi[3] = (dir != 3)*mhi[3] + (dir == 3), lo[3] = (dir != 3)*mlo[3] + (dir == 3), num[3] = (hi[3] - lo[3])*num[2], den[3] = den[2]; p0 + p1 + p2 + p3           <= m_degree; ++p3, hi[3] *= mhi[3], lo[3] *= mlo[3], num[3] = (hi[3] - lo[3])*num[2], den[3] = (p3 + (dir != 3))*den[2]),
              for (p4 = 0, hi[4] = (dir != 4)*mhi[4] + (dir == 4), lo[4] = (dir != 4)*mlo[4] + (dir == 4), num[4] = (hi[4] - lo[4])*num[3], den[4] = den[3]; p0 + p1 + p2 + p3 + p4      <= m_degree; ++p4, hi[4] *= mhi[4], lo[4] *= mlo[4], num[4] = (hi[4] - lo[4])*num[3], den[4] = (p4 + (dir != 4))*den[3]),
              for (p5 = 0, hi[5] = (dir != 5)*mhi[5] + (dir == 5), lo[5] = (dir != 5)*mlo[5] + (dir == 5), num[5] = (hi[5] - lo[5])*num[4], den[5] = den[4]; p0 + p1 + p2 + p3 + p4 + p5 <= m_degree; ++p5, hi[5] *= mhi[5], lo[5] *= mlo[5], num[5] = (hi[5] - lo[5])*num[4], den[5] = (p5 + (dir != 5))*den[4]))
        {
          // p is from derivative, a_h is to make the derivative "undivided"
          *a_gradXip++ = p[dir]*a_h[dir]*num[SpaceDim-1]/den[SpaceDim-1];
        }
    }
}

#include "NamespaceFooter.H"
