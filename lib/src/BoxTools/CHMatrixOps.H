#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef _CHMATRIXOPS_H_
#define _CHMATRIXOPS_H_


/******************************************************************************/
/**
 * \file CHMatrixOps.H
 *
 * \brief Blas/Lapack matrix operations for the CHMatrix, CHVector types
 *
 * Simpler interfaces to Blas/Lapack if the array types are CHMatrix and/or
 * CHVector. This implies real ArColZero type arrays with the default allocator.
 * For other types of arrays, you have to call Blas/Lapack directly.
 *
 *//*+*************************************************************************/

//--Forward declarations
#include "CHArray.H"

#ifdef CH_SPACEDIM // In Chombo land
#include "NamespaceHeader.H"
#endif

/// vector-vector multiplication
Real CHdot(const CHVector& a_X,
           const CHVector& a_Y);

/// vector-vector multiplication
Real CHdot(int size,
           const Real* a_X,
           int a_xStride,
           const Real* a_Y,
           int a_yStride);

/// Matrix-vector multiplication
void CHgemv(const CHMatrix& a_A,
            const CHVector& a_X,
            CHVector&       a_Y,
            char            a_transposeA = 'N',
            int             a_M = 0,
            int             a_N = 0,
            Real            a_alpha = 1.,
            Real            a_beta = 0.);

/// Matrix-vector multiplication
void CHgemv(const Real* a_A,
            const Real* a_X,
            Real* a_Y,
            char            a_transposeA,
            int             a_M,
            int             a_N,
            Real            a_alpha = 1.,
            Real            a_beta = 0.);

/// Matrix-matrix multiplication
void CHgemm(const CHMatrix& a_A,
            const CHMatrix& a_B,
            CHMatrix&       a_C,
            char            a_transposeA = 'N',
            char            a_transposeB = 'N',
            int             a_M = 0,
            int             a_N = 0,
            int             a_K = 0,
            Real            a_alpha = 1.,
            Real            a_beta = 0.);

/// Symmetric matrix creation by multiplication with self
void CHsyrk(const CHMatrix& a_A,
            CHMatrix&       a_C,
            char            a_transposeA = 'N',
            char            a_uploC = 'L',
            int             a_N = 0,
            int             a_K = 0,
            Real            a_alpha = 1.,
            Real            a_beta = 0.);

/// Symmetric matrix-matrix multiplication
void CHsymm(const CHMatrix& a_A,
            const CHMatrix& a_B,
            CHMatrix&       a_C,
            char            a_sideA = 'L',
            char            a_uploC = 'L',
            int             a_M = 0,
            int             a_N = 0,
            Real            a_alpha = 1.,
            Real            a_beta = 0.);

/// Invert using LU factorization
void CHinverse(CHMatrix& a_A,
               CHArray<int, 1>& a_ipiv,
               CHVector& a_work);

/// Invert using a Cholesky factorization
/// (requires a_A to be symmetric positive definite)
void CHinverseCholesky(CHMatrix& a_A,
                       char      a_uploA = 'A');

/// Cholesky factorization
void CHcholesky(CHMatrix& a_A,
                char      a_uploA = 'A');

/// QR factorization (Q & R matrices returned)
void CHgeqrf(const CHMatrix& a_A,
             CHMatrix&       a_Q,
             CHMatrix&       a_R,
             CHVector&       a_work,
             int             a_M = 0,
             int             a_N = 0);

/// QR factorization (LAPACK result)
void CHgeqrf(CHMatrix& a_A,
             CHVector& a_tau,
             CHVector& a_work,
             int       a_M,
             int       a_N);

/// Find Cholesky factorization
void CHpbtrf(CHMatrix& a_A,
             char      a_uploA = 'L');

/// Compute the solution to a real system of linear equation A*x = b
void CHgesv(CHMatrix& a_A,
            CHVector& a_B,
            CHArray<int, 1>& a_ipiv);

/// Compute the solution for multiple systems of linear equations A*x = b
void CHgesv(CHMatrix& a_A,
            CHMatrix& a_B,
            CHArray<int, 1>& a_ipiv);

/// Compute the solution to a real system of linear equation A*x = b
void CHsgesv(CHMatrix& a_A,
             CHVector& a_B,
             CHVector& a_X,
             CHVector& a_work,
             CHVector& a_swork,
             CHArray<int, 1>& a_ipiv,
             int a_N);

/// Compute the solution to a banded system of linear equation A*x = b
void CHgbsv(CHMatrix& a_A,
            int a_kl,
            int a_ku,
            CHVector& a_b);		

/// Compute the eigenvalues of A
void CHgeev(CHMatrix& a_A,
            CHVector& a_WR,
            CHVector& a_WI);

/// Compute the eigenvalues and eigenvectors of A
void CHgeev(CHMatrix& a_A,
            CHVector& a_WR,
            CHVector& a_WI,
            CHMatrix& a_VR,
            char      a_jobR = 'V');

/// Compute the eigenvalues and eigenvectors of A
void CHgeev(CHMatrix& a_A,
            CHVector& a_WR,
            CHVector& a_WI,
            CHMatrix& a_VR,
            CHMatrix& a_VL,
            char      a_jobR = 'V',
            char      a_jobL = 'V');

/// Compute the singular value decomposition
void CHgesvd(CHMatrix& a_A,
             CHMatrix& a_U,
             CHMatrix& a_S,
             char      a_jobU = 'A',
             char      a_jobVT = 'N');

/// Compute the solution to the Least Squares problem
void CHgels(CHMatrix& a_A,
            CHVector& a_b,
            char      a_trans = 'N');

///  Compute the determinant
/**  Warning this does not call a LAPACK function. This because LAPACK
 *   does not have methods for determinants due to accuracy concerns.
 *   This exists here but should be avoided for all by small matrices.
 *   \param[in]  a_A   Matrix to solve determinate of
 *   \param[in]  N     Size of square matrix to solve determinate of
 */
inline Real determinate(const CHMatrix& a_A,
                        const int N)
{
  CH_assert(a_A.size(0) >= N);
  CH_assert(a_A.size(1) >= N);
  Real det = 0.;
  
  if(N == 1)
    {
      det = a_A(0,0);
    }
  else if(N == 2)
    {
      det = a_A(0,0)*a_A(1,1) - a_A(1,0)*a_A(0,1);
      return det;
    }
  else
    {
      CHMatrix C(N,N);
      for(int p = 0; p != N; ++p)
        {
          int h = 0;
          int k = 0;
          for(int i = 1; i != N; ++i)
            {
              for(int j = 0; j != N; ++j)
                {
                  if(j == p)
                    {
                      continue;
                    }
                  C(h,k) = a_A(i,j);
                  k++;
                  if(k == N-1)
                    {
                      h++;
                      k = 0;
                    }
                }
            }
          det = det + a_A(0,p)*pow(-1.,p)*determinate(C, N-1);
        }
    }
  return det;
}

/// Fill the transpose of a matrix
/** \param[out] a_At    Transpose of matrix a_A
 *  \param[in]  a_A     Matrix input
 */
inline void transpose(CHMatrix& a_At,
                      const CHMatrix& a_A)
{
  CH_assert(a_At.size(0) >= a_A.size(1));
  CH_assert(a_At.size(1) >= a_A.size(0));
  for(int row = 0; row != a_A.size(0); ++row)
    {
      for(int col = 0; col != a_A.size(1); ++col)
        {
          a_At(col,row) = a_A(row,col);
        }
    }
}

/// Transpose a matrix in place. must be square
/** \param[out] a_A   Transpose of matrix a_A
 *  \param[in]  a_A   Matrix input
 */
inline void transpose(CHMatrix& a_A)
{
  CH_assert(a_A.size(0) == a_A.size(1));
  for(int row = 0; row != a_A.size(0); ++row)
    {
      for(int col = row; col != a_A.size(1); ++col)
        {
          Real tmp = a_A(col,row);
          a_A(col,row) = a_A(row,col);
          a_A(row,col) = tmp;
        }
    }
}

/// Convert a symmetric matrix stored only on LD to full LDU
/** \param[in]  a_N     Order of matrix
 *  \param[in]  a_A     Symmetric matrix stored on LD
 *  \param[out] a_A     Symmetric matric stored on LDU
 */
inline void sy2ge(const int a_n, CHMatrix &a_A)
{
  CH_assert(a_A.size(0) >= a_n);
  CH_assert(a_A.size(1) >= a_n);
  for (int j = 1; j != a_n; ++j)
    {
      const int iEnd = j;
      for (int i = 0; i != iEnd; ++i)
        {
          a_A(i, j) = a_A(j, i);
        }
    }
}

#ifdef CH_SPACEDIM // In Chombo land
#include "NamespaceFooter.H"
#endif

#endif
