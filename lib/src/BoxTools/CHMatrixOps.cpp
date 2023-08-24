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
 * \file CHMatrixOps.cpp
 *
 * \brief Blas/Lapack matrix operations for the CHMatrix, CHVector types
 *
 *//*+*************************************************************************/

#include <cmath>
#include <cstring>
#include "CHMatrixOps.H"
#include "CHArray.H"

#ifdef CH_SPACEDIM // In Chombo land
#include "NamespaceHeader.H"
#endif

/*==============================================================================
 * Define how is a fortran function named
 *============================================================================*/
#ifdef CH_FORT_UPPERCASE
  #ifdef CH_FORT_UNDERSCORE
    #define FORTRAN_NAME( NAME ,name ) NAME ## _
  #else
    #define FORTRAN_NAME( NAME ,name ) NAME
  #endif
#else
  #ifdef CH_FORT_UNDERSCORE
    #define FORTRAN_NAME( NAME ,name ) name ## _
  #else
    #define FORTRAN_NAME( NAME ,name ) name
  #endif
#endif

/*==============================================================================
 * BLAS & LAPACK routines
 *============================================================================*/
#ifdef CH_USE_DOUBLE
  #define BLAS1_DOT FORTRAN_NAME(DDOT, ddot)
  #define BLAS1_NRM2 FORTRAN_NAME(DNRM2, dnrm2)
  #define BLAS2_GEMV FORTRAN_NAME(DGEMV, dgemv)
  #define BLAS3_GEMM FORTRAN_NAME(DGEMM, dgemm)
  #define BLAS3_SYRK FORTRAN_NAME(DSYRK, dsyrk)
  #define BLAS3_SYMM FORTRAN_NAME(DSYMM, dsymm)
  #define LAPACK_GETRF FORTRAN_NAME(DGETRF, dgetrf)
  #define LAPACK_GETRI FORTRAN_NAME(DGETRI, dgetri)
  #define LAPACK_POTRF FORTRAN_NAME(DPOTRF, dpotrf)
  #define LAPACK_POTRI FORTRAN_NAME(DPOTRI, dpotri)
  #define LAPACK_GEQRF FORTRAN_NAME(DGEQRF, dgeqrf)
  #define LAPACK_GESV FORTRAN_NAME(DGESV, dgesv)
  #define LAPACK_SGESV FORTRAN_NAME(DSGESV, dsgesv)
  #define LAPACK_PBTRF FORTRAN_NAME(DPBTRF, dpbtrf)
  #define LAPACK_GEEV FORTRAN_NAME(DGEEV, dgeev)
  #define LAPACK_GESVD FORTRAN_NAME(DGEEV, dgesvd)
  #define LAPACK_GBSV FORTRAN_NAME(DGBSV, dgbsv)
  #define LAPACK_GELS FORTRAN_NAME(DGELS, dgels)
#else
  #define BLAS1_DOT FORTRAN_NAME(SDOT, sdot)
  #define BLAS1_NRM2 FORTRAN_NAME(SNRM2, snrm2)
  #define BLAS2_GEMV FORTRAN_NAME(SGEMV, sgemv)
  #define BLAS3_GEMM FORTRAN_NAME(SGEMM, sgemm)
  #define BLAS3_SYRK FORTRAN_NAME(SSYRK, ssyrk)
  #define BLAS3_SYMM FORTRAN_NAME(SSYMM, ssymm)
  #define LAPACK_GETRF FORTRAN_NAME(SGETRF, sgetrf)
  #define LAPACK_GETRI FORTRAN_NAME(SGETRI, sgetri)
  #define LAPACK_POTRF FORTRAN_NAME(SPOTRF, spotrf)
  #define LAPACK_POTRI FORTRAN_NAME(SPOTRI, spotri)
  #define LAPACK_GEQRF FORTRAN_NAME(SGEQRF, sgeqrf)
  #define LAPACK_GESV FORTRAN_NAME(SGESV, sgesv)
  #define LAPACK_SGESV FORTRAN_NAME(SSGESV, ssgesv)
  #define LAPACK_PBTRF FORTRAN_NAME(SPBTRF, spbtrf)
  #define LAPACK_GEEV FORTRAN_NAME(SGEEV, sgeev)
  #define LAPACK_GESVD FORTRAN_NAME(SGEEV, sgesvd)
  #define LAPACK_GBSV FORTRAN_NAME(SGBSV, sgbsv)
  #define LAPACK_GELS FORTRAN_NAME(SGELS, sgels)
#endif

extern "C" Real BLAS1_NRM2(int*, Real*, int*);
extern "C" Real BLAS1_DOT(int*, Real*, int*, Real*, int*);
extern "C" void BLAS2_GEMV(char*, int*, int*, Real*, Real*, int*, Real*, int*,
                           Real*, Real*, int*);
extern "C" void BLAS3_GEMM(char*, char*, int*, int*, int*, Real*, Real*, int*,
                           Real*, int*, Real*, Real*, int*);
extern "C" void BLAS3_SYRK(char*, char*, int*, int*, Real*, Real*, int*, Real*,
                           Real*, int*);
extern "C" void BLAS3_SYMM(char*, char*, int*, int*, Real*, Real*, int*, Real*,
                           int*, Real*, Real*, int*);
extern "C" void LAPACK_GETRF(int*, int*, Real*, int*, int*, int*);
extern "C" void LAPACK_GETRI(int*, Real*, int*, int*, Real*, int*, int*);
extern "C" void LAPACK_POTRF(char*, int*, Real*, int*, int*);
extern "C" void LAPACK_POTRI(char*, int*, Real*, int*, int*);
extern "C" void LAPACK_GEQRF(int*, int*, Real*, int*, Real*, Real*, int*, int*);
extern "C" void LAPACK_GESV(int*, int*, Real*, int*, int*, Real*, int*, int*);
extern "C" void LAPACK_SGESV(int*, int*, Real*, int*, int*, Real*, int*, Real*,
                             int*, Real*, Real*, int*, int*);
extern "C" void LAPACK_PBTRF(char*, int*, int*, Real*, int*, int*);
extern "C" void LAPACK_GEEV(char*, char*, int*, Real*, int*, Real*, Real*,
                            Real*, int*, Real*, int*, Real*, int*, int*);
extern "C" void LAPACK_GESVD(char*, char*, int*, int*, Real*, int*, Real*,
                             Real*, int*, Real*, int*, Real*, int*, int*);
extern "C" void LAPACK_GBSV(int* , int*, int*, int*, Real*, int*, int*,
                            Real*, int*, int*);
extern "C" void LAPACK_GELS(char*, int*, int*, int*, Real*, int*,
                            Real*, int*, Real*, int*, int*);

/*==============================================================================
 * Helper routines
 *============================================================================*/

//--Parse transpose characters

inline bool isTranspose(const char a_ch)
{
  bool transpose = false;
  switch (a_ch)
    {
    case 'n':
    case 'N':
      break;
    case 't':
    case 'T':
    case 'c':
    case 'C':
      transpose = true;
      break;
    default:
      MayDay::Error("Invalid transpose character");
      break;
    }
  return transpose;
}

//--Parse up-lo characters

inline bool isFull(char& a_ch)
{
  bool full = false;
  switch (a_ch)
    {
    case 'l':
    case 'L':
    case 'u':
    case 'U':
      break;
    case 'a':
    case 'A':
      full = true;
      a_ch = 'L';
      break;
    default:
      MayDay::Error("Invalid up-lo character");
      break;
    }
  return full;
}

//--Parse usage characters

inline bool isJobz(char& a_ch)
{
  bool jobz = false;
  switch (a_ch)
    {
    case 'n':
    case 'N':
      break;
    case 'v':
    case 'V':
      jobz = true;
      break;
    default:
      MayDay::Error("Invalid jobz character");
      break;
    }
  return jobz;
}

/*--------------------------------------------------------------------*/
///  Vector-Vector multiplication
/** Solve for
 *  \f[
 *    y = x\mbox{op} y
 *  \f]
 *
 *   \param[in]  a_X    Input Vector X
 *   \param[in]  a_Y    Input Vector Y
 *   \return            Scalar product
 *//*-----------------------------------------------------------------*/

Real CHdot(const CHVector& a_X,
           const CHVector& a_Y)
{
  CH_assert(a_X.size() == a_Y.size());
  int size = a_X.size(0);
  int stride = 1;

  return BLAS1_DOT(&size,
                   const_cast<Real*>(a_X.begin()), &stride,
                   const_cast<Real*>(a_Y.begin()), &stride);
}

Real CHdot(int size,
           const Real* a_X,
           int a_xStride,
           const Real* a_Y,
           int a_yStride)
{
  return BLAS1_DOT(&size,
                   const_cast<Real*>(a_X), &a_xStride,
                   const_cast<Real*>(a_Y), &a_yStride);
}

/*--------------------------------------------------------------------*/
///  Matrix-vector multiplication
/** Solve for
 *  \f[
 *    y = \alpha\mbox{op}(A)x + \beta y
 *  \f]
 *  where \f$\mbox{op}(A) = A\f$ if a_transposeA is 'n' or 'N' and
 *  \f$\mbox{op}(A) = A^T\f$ if a_transposeA is 't', 'T', 'c', or 'C'.
 *
 *   \param[in]  a_A    Matrix A
 *   \param[in]  a_X    Vector X
 *   \param[in]  a_Y    Input Vector Y
 *   \param[out] a_Y    Resulting vector Y
 *   \param[in]  a_transposeA
 *                      Flag indicating if A is transposed or not,
 *                      default not transposed
 *   \param[in]  a_M    Rows of A, default current size of A
 *   \param[in]  a_N    Columns of A, default current size of A
 *   \param[in]  a_alpha
 *                      Scalar multiplied by X, default 1
 *   \param[in]  a_beta
 *                      Scalar multiplied by Y, default 0
 *//*-----------------------------------------------------------------*/

void CHgemv(const CHMatrix& a_A,
            const CHVector& a_X,
            CHVector&       a_Y,
            char            a_transposeA,
            int             a_M,
            int             a_N,
            Real            a_alpha,
            Real            a_beta)
{
  const int transposeA = isTranspose(a_transposeA);
  if (a_M == 0) a_M = a_A.size(0);
  if (a_N == 0) a_N = a_A.size(1);
  CH_assert(a_A.size(0) >= a_M);
  CH_assert(a_A.size(1) >= a_N);
  CH_assert(a_X.size() >= a_A.size(!transposeA));
  CH_assert(a_Y.size() >= a_A.size(transposeA));
  int lfi_LDA = a_A.size(0);
  int lfi_one = 1;

  BLAS2_GEMV(&a_transposeA,
             &a_M, &a_N,
             &a_alpha,
             const_cast<Real*>(a_A.begin()), &lfi_LDA,
             const_cast<Real*>(a_X.begin()), &lfi_one,
             &a_beta,
             a_Y.begin(), &lfi_one);
}


void CHgemv(const Real* a_A,
            const Real* a_X,
            Real*       a_Y,
            char            a_transposeA,
            int             a_M,
            int             a_N,
            Real            a_alpha,
            Real            a_beta)
{
  int lfi_one = 1;

  BLAS2_GEMV(&a_transposeA,
             &a_M, &a_N,
             &a_alpha,
             const_cast<Real*>(a_A), &a_M,
             const_cast<Real*>(a_X), &lfi_one,
             &a_beta,
             a_Y, &lfi_one);
}

/*--------------------------------------------------------------------*/
///  Matrix-matrix multiplication
/** Solves of
 *  \f[
 *    C = \alpha\mbox{op}(A)\mbox{op}(B) + \beta C
 *  \f]
 *  where \f$\mbox{op}(A) = A\f$ if a_transposeA is 'n' or 'N' and
 *  \f$\mbox{op}(A) = A^T\f$ if a_transposeA is 't', 'T', 'c', or 'C'.
 *  \f$\mbox{op}(B) = B^T\f$ if a_transposeA is 't', 'T', 'c', or 'C'.
 *
 *   \param[in]  a_A    Matrix A
 *   \param[in]  a_B    Matrix B
 *   \param[in]  a_C    Input matrix C
 *   \param[out] a_C    Resulting matrix C
 *   \param[in]  a_transposeA
 *                      Flag indicating if A is transposed or not,
 *                      default not transposed
 *   \param[in]  a_transposeB
 *                      Flag indicating if B is transposed or not,
 *                      default not transposed
 *   \param[in]  a_M    Rows of A, default current size of A
 *   \param[in]  a_N    Columns of B, default current size of B
 *   \param[in]  a_K    Columns of A, default current size of A
 *   \param[in]  a_alpha
 *                      Scalar multiplied by B, default 1
 *   \param[in]  a_beta
 *                      Scalar multiplied by C, default 0
 *//*-----------------------------------------------------------------*/

void CHgemm(const CHMatrix& a_A,
            const CHMatrix& a_B,
            CHMatrix&       a_C,
            char            a_transposeA,
            char            a_transposeB,
            int             a_M,
            int             a_N,
            int             a_K,
            Real            a_alpha,
            Real            a_beta)
{
  const int transposeA = isTranspose(a_transposeA);
  const int transposeB = isTranspose(a_transposeB);
  if (a_M == 0) a_M = a_A.size(transposeA);
  if (a_N == 0) a_N = a_B.size(!transposeB);
  if (a_K == 0) a_K = a_A.size(!transposeA);
  CH_assert(a_A.size(transposeA)  >= a_M);
  CH_assert(a_A.size(!transposeA) >= a_K);
  CH_assert(a_B.size(!transposeB) >= a_N);
  CH_assert(a_B.size(transposeB)  >= a_K);
  int lfi_LDA = a_A.size(0);
  int lfi_LDB = a_B.size(0);
  int lfi_LDC = a_C.size(0);

  BLAS3_GEMM(&a_transposeA, &a_transposeB,
             &a_M, &a_N, &a_K,
             &a_alpha,
             const_cast<Real*>(a_A.begin()), &lfi_LDA,
             const_cast<Real*>(a_B.begin()), &lfi_LDB,
             &a_beta,
             a_C.begin(), &lfi_LDC);
}

/*--------------------------------------------------------------------*/
/// Symmetric matrix creation by multiplication with self
/**
 *  Performs one of the symmetric rank k operations
 *    C := alpha*A*A**T + beta*C,
 *  or
 *    C := alpha*A**T*A + beta*C,
 * where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
 * and  A  is an  n by k  matrix in the first case and a  k by n  matrix
 * in the second case.
 *
 *  \param[in]  a_A     General matrix
 *  \param[in]  a_C     Symmetric n by n matrix
 *  \param[in]  a_transposeA
 *                      Flag indicating if A is transposed or not,
 *                      default not transposed
 *  \param[in]  a_uploC Flag indicating to use the 'U' upper or 'L' lower
 *                      triangular portion of A
 *                      default is 'L'
 *  \param[in]  a_N     Rows of A, default current size of A
 *  \param[in]  a_K     Columns of A, default current size of A
 *  \param[in]  a_alpha Scalar multiplied by A*A**T, default 1
 *  \param[in]  a_beta  Scalar multiplied by C, default 0
 *  \param[out] a_C     Resulting matix
 *//*-----------------------------------------------------------------*/
void CHsyrk(const CHMatrix& a_A,
            CHMatrix&       a_C,
            char            a_transposeA,
            char            a_uploC,
            int             a_N,
            int             a_K,
            Real            a_alpha,
            Real            a_beta)
{
  const int transposeA = isTranspose(a_transposeA);
  if (a_N == 0) a_N = a_A.size(transposeA);
  if (a_K == 0) a_K = a_A.size(!transposeA);
  CH_assert(a_A.size(transposeA)  >= a_N);
  CH_assert(a_A.size(!transposeA) >= a_K);
  CH_assert(a_C.size(0) >= a_N);
  CH_assert(a_C.size(0) >= a_N);
  int lfi_LDA = a_A.size(0);
  int lfi_LDC = a_C.size(0);

  BLAS3_SYRK(&a_uploC,
             &a_transposeA,
             &a_N, &a_K,
             &a_alpha,
             const_cast<Real*>(a_A.begin()), &lfi_LDA,            
             &a_beta,
             a_C.begin(), &lfi_LDC);
}

/*--------------------------------------------------------------------*/
///  Symmetric matrix-matrix multiplication
/**
 *  Performs one of the matrix-matrix operations
 *    C := alpha*A*B + beta*C,
 *  or
 *    C := alpha*B*A + beta*C,
 *  where alpha and beta are scalars, A is a symmetric matrix and B and
 *  C are m by n matrices.
 *
 *  \param[in]  a_A     Symmetric matrix
 *  \param[in]  a_B     General m by n matrix
 *  \param[in]  a_C     General m by n matrix
 *  \param[in]  a_sideA Flag indicating to solve A*B or B*A
 *                      'L' is A*B, and 'R' is B*A
 *                      default is 'L'
 *  \param[in]  a_uploC Flag indicating to use the 'U' upper or 'L' lower
 *                      triangular portion of A
 *                      default is 'L'
 *  \param[in]  a_M     Rows of C, default current size of C
 *  \param[in]  a_N     Columns of C, default current size of C
 *  \param[in]  a_alpha Scalar multiplied by A*B, default 1
 *  \param[in]  a_beta  Scalar multiplied by C, default 0
 *  \param[out] a_C     Resulting matix
 *//*-----------------------------------------------------------------*/
void CHsymm(const CHMatrix& a_A,
            const CHMatrix& a_B,
            CHMatrix&       a_C,
            char            a_sideA,
            char            a_uploC,
            int             a_M,
            int             a_N,
            Real            a_alpha,
            Real            a_beta)
{
  if (a_N == 0) a_N = a_C.size(0);
  if (a_M == 0) a_M = a_C.size(1);
  CH_assert(a_A.size(0) >= a_M);
  CH_assert(a_A.size(1) >= a_N);
  CH_assert(a_B.size(0) >= a_M);
  CH_assert(a_B.size(1) >= a_N);
  CH_assert(a_C.size(0) >= a_M);
  CH_assert(a_C.size(0) >= a_N);
  int lfi_LDA = a_A.size(0);
  int lfi_LDB = a_B.size(0);
  int lfi_LDC = a_C.size(0);

  BLAS3_SYMM(&a_sideA,
             &a_uploC,
             &a_M, &a_N,
             &a_alpha,
             const_cast<Real*>(a_A.begin()), &lfi_LDA,
             const_cast<Real*>(a_B.begin()), &lfi_LDA,
             &a_beta,
             a_C.begin(), &lfi_LDC);
}

/*--------------------------------------------------------------------*/
///  Invert using LU factorization
/**
 *  \param[in]  a_A     General matrix to invert
 *  \param[out] a_A     Inverse matrix
 *  \param[out] a_ipiv  An integer array used as working space
 *  \param[out] a_work  A vector used as working space
 *//*-----------------------------------------------------------------*/

void CHinverse(CHMatrix& a_A,
               CHArray<int, 1>& a_ipiv,
               CHVector& a_work)
{
  int ier;
  int lfi_N = a_A.size(0);
  CH_assert(a_A.size(1) >= lfi_N);
  int lfi_LDA = lfi_N;

  // Setup pivots if required
  if ((!a_ipiv.isAllocated()) || a_ipiv.size() < lfi_N)
    {
      a_ipiv.define(lfi_N);
    }

  // Setup work if required
  int lwork = a_work.size();
  if ((!a_work.isAllocated()) || a_work.size() < lfi_N)
    {
      lwork = -1;
      Real query;
      // Get the optimal work array size
      LAPACK_GETRI(&lfi_N, a_A.begin(), &lfi_LDA, a_ipiv.begin(),
                   &query, &lwork, &ier);
      if (ier != 0)
        {
          pout() << "GETRI failure status: " << ier << std::endl;
          MayDay::Error("LAPACK routine @GETRI failed");
        }
      lwork = (int)query;
      a_work.define(lwork);
    }

  LAPACK_GETRF(&lfi_N, &lfi_N, a_A.begin(), &lfi_LDA, a_ipiv.begin(), &ier);
  if (ier != 0)
    {
      pout() << "GETRF failure status: " << ier << std::endl;
      MayDay::Error("LAPACK routine @GETRF failed");
    }

  LAPACK_GETRI(&lfi_N, a_A.begin(), &lfi_LDA, a_ipiv.begin(), a_work.begin(),
               &lwork, &ier);
  if (ier != 0)
    {
      pout() << "GETRI failure status: " << ier << std::endl;
      MayDay::Error("LAPACK routine @GETRI failed");
    }
}

/*--------------------------------------------------------------------*/
/// Invert using a Cholesky factorization
/** Matrix 'a_A' must be symmetric positive definite
 *
 *  \param[in]  a_A     Symmetric positive definite matrix.  The order
 *                      of the matrix is taken from dimension 0.
 *                      Dimension 1 must be >= dimension 0.
 *  \param[out] a_A     Inverse matrix
 *  \param[in]  a_uploA 'l' or 'L': I/O of a_A is on lower triangular
 *                      'u' or 'U': I/O of a_A is on upper triangular
 *                      'a' or 'A': (default) Input of a_A from lower
 *                                  triangluar part.  On output, a_A
 *                                  has correct values everywhere.
 *//*-----------------------------------------------------------------*/

void CHinverseCholesky(CHMatrix& a_A,
                       char      a_uploA)
{
  const int fullA = isFull(a_uploA);
  int lfi_N = a_A.size(0);
  CH_assert(a_A.size(1) == lfi_N);
  int lfi_LDA = lfi_N;

  int ier = 0;
  LAPACK_POTRF(&a_uploA, &lfi_N, a_A.begin(), &lfi_LDA, &ier);
  if (ier != 0)
    {
      pout() << "POTRF failure status: " << ier << std::endl;
      MayDay::Error("LAPACK routine @POTRF failed");
    }
  LAPACK_POTRI(&a_uploA, &lfi_N, a_A.begin(), &lfi_LDA, &ier);
  if (ier != 0)
    {
      pout() << "POTRI failure status: " << ier << std::endl;
      MayDay::Error("LAPACK routine @POTRI failed");
    }

  if (fullA)  // Populate upper triangular part of symmetric matrix
    {
      sy2ge(lfi_N, a_A);
    } 
}

/*--------------------------------------------------------------------*/
/// Find Cholesky factorization
/** Matrix 'a_A' must be symmetric positive definite
 *
 *  \param[in]  a_A     Symmetric positive definite matrix.  The order
 *                      of the matrix is taken from dimension 0.
 *                      Dimension 1 must be >= dimension 0.
 *  \param[out] a_A     factored matrix
 *  \param[in]  a_uploA 'a' or 'A': (default) Input of a_A from lower
 *                                  triangluar part.  On output, a_A
 *                                  has correct values everywhere.
 *//*-----------------------------------------------------------------*/

void CHcholesky(CHMatrix& a_A,
                char      a_uploA)
{
  const int fullA = isFull(a_uploA);
  int lfi_N = a_A.size(0);
  CH_assert(a_A.size(1) == lfi_N);
  int lfi_LDA = lfi_N;

  int ier = 0;
  LAPACK_POTRF(&a_uploA, &lfi_N, a_A.begin(), &lfi_LDA, &ier);
  if (ier != 0)
    {
      pout() << "POTRF failure status: " << ier << std::endl;
      MayDay::Error("LAPACK routine @POTRF failed");
    }
    if (fullA)  // Populate upper triangular part of symmetric matrix
    {
      sy2ge(lfi_N, a_A);
    } 
  // zero lower triangle
  for(int i = 0; i < a_A.size(0); i++)
    {
      for(int j = 0; j < i; j++)
        {
          a_A(i, j) = 0.0;
        }
    }
}

/*--------------------------------------------------------------------*/
/// QR factorization (Q & R matrices returned)
/** This routine unpacks the LAPACK results (returned in A and tau)
 *  to populate Q & R.  The work array will be redefined if it is not
 *  allocated or of insufficient size.
 *
 *  \param[in]  a_A    Input matrix A
 *  \param[out] a_Q    Matrix Q
 *  \param[out] a_R    Matrix R
 *  \param[out] a_work Vector used as working space
 *  \param[in]  a_M    Rows of A
 *  \param[in]  a_N    Columns of A
 *//*-----------------------------------------------------------------*/

void CHgeqrf(const CHMatrix& a_A,
             CHMatrix&       a_Q,
             CHMatrix&       a_R,
             CHVector&       a_work,
             int             a_M,
             int             a_N)
{
  int ier = 0;
  if (a_M == 0) a_M = a_A.size(0);
  if (a_N == 0) a_N = a_A.size(1);
  int minMN = std::min(a_M, a_N);
  CH_assert(a_Q.size(0) >= a_M);
  CH_assert(a_Q.size(1) >= a_M);
  CH_assert(a_R.size(0) >= minMN);
  CH_assert(a_R.size(1) >= a_N);

  // Copy a_A into A so it can be modified
  CHMatrix A(a_M, a_N);
  {
    const int MByte = a_M*sizeof(Real);
    for (int j = 0; j != a_N; ++j)
    {
      std::memcpy(&A(0, j), &a_A(0, j), MByte);
    }
  }
  CHVector tau(minMN);
  // Setup work if required
  int lwork = a_work.size();
  // if ((!a_work.isAllocated()) || a_work.size() < a_N)
  //   {
  //     lwork = -1;
  //     Real query;
  //     // Get the optimal work array size
  //     LAPACK_GEQRF(&a_M, &a_N, A.begin(), &a_M, tau.begin(),
  //                  &query, &lwork, &ier);
  //     if (ier != 0)
  //       {
  //         pout() << "GEQRF failure status: " << ier << std::endl;
  //         MayDay::Error("LAPACK routine @GEQRF failed");
  //       }
  //     lwork = (int)query;
  //     a_work.define(lwork);
  //   }

  // Do the factorization
  LAPACK_GEQRF(&a_M, &a_N, A.begin(), &a_M, tau.begin(), a_work.begin(), &lwork,
               &ier);
  if (ier != 0)
    {
      pout() << "GEQRF failure status: " << ier << std::endl;
      MayDay::Error("LAPACK routine @GEQRF failed");
    }

  // Populate R
  a_R = 0.;
  for (int j = 0; j != a_N; ++j)
    {
      const int imax = std::min(j, minMN);
      for (int i = 0; i <= imax; ++i)
        {
          a_R(i, j) = A(i, j);
        }
    }

  // Populate Q

  // Initialize Q to identity
  CHMatrix I(a_M, a_M);
  I   = 0.;
  a_Q = 0.;
  for (int i = 0; i != a_M; ++i)
    {
      I(i, i)   = 1.;
      a_Q(i, i) = 1.;
    }
  CHMatrix Qtmp(a_M, a_M);
  // Note that these pointers are swapped when entering the loop below
  CHMatrix *Qn = &a_Q;   // Next update of Q
  CHMatrix *Qc = &Qtmp;  // Current Q

  CHMatrix H(a_M, a_M);
  for (int k = 0; k != minMN; ++k)
    {
      const Real ctau = tau(k);
      CHMatrix *Qt = Qc;
      Qc = Qn;
      Qn = Qt;

      // Perform H(i) = I - tau * v * v**T
      std::memcpy(H.begin(), I.begin(), H.size()*sizeof(Real));
      H(k, k) -= ctau;
      if (k+1 != a_M)
        {
          CHArray<Real, 1, ArRangeCol> v;
          v.define(&A(k+1, k), CHRange(k+1, a_M-1));
          // Remainder of (k,:)
          for (int i = k+1; i != a_M; ++i)
            {
              H(k, i) -= ctau*v(i);
            }
          // Remainder of (:,k)
          for (int i = k+1; i != a_M; ++i)
            {
              H(i, k) -= ctau*v(i);
            }
          // Everything else
          for (int j = k+1; j != a_M; ++j)
            {
              const Real vjctau = ctau*v(j);
              for (int i = k+1; i != a_M; ++i)
                {
                  H(i, j) -= vjctau*v(i);
                }
            }
        }

      // Qn = Qc.H
      CHgemm((*Qc), H, (*Qn), 'N', 'N', a_M, a_M, a_M);
    }

  // But Qn might point to Qtmp.  If so, copy Qtmp to Q
  if (Qn != &a_Q)
    {
      const int MByte = a_M*sizeof(Real);
      for (int j = 0; j != a_M; ++j)
        {
          std::memcpy(&a_Q(0, j), &Qtmp(0, j), MByte);
        }
    }
}

/*--------------------------------------------------------------------*/
///  QR factorization (LAPACK result)
/** This simple but efficient routine returns the results as LAPACK
 *  itself does (in a_A and a_tau).  All arguments must be fully
 *  and correctly specified
 *
 *  \param[in]  a_A    Input matrix A
 *  \param[out] a_A    Output matrix X
 *  \param[out] a_tau  Vector
 *  \param[out] a_work Vector used as working space
 *  \param[in]  a_M    Rows of A
 *  \param[in]  a_N    Columns of A
 *//*-----------------------------------------------------------------*/

void CHgeqrf(CHMatrix& a_A,
             CHVector& a_tau,
             CHVector& a_work,
             int       a_M,
             int       a_N)
{
  int ier = 0;
  CH_assert(a_M <= a_A.size(0) && a_N <= a_A.size(1));
  CH_assert(std::min(a_M, a_N) < a_tau.size());
  CH_assert(a_N <= a_work.size());
  int lfi_LDA   = a_A.size(0);
  int lfi_lwork = a_work.size();
  // Do the factorization
  LAPACK_GEQRF(&a_M, &a_N, a_A.begin(), &lfi_LDA, a_tau.begin(), a_work.begin(),
               &lfi_lwork, &ier);
  if (ier != 0)
    {
      pout() << "GEQRF failure status: " << ier << std::endl;
      MayDay::Error("LAPACK routine @GEQRF failed");
    } 
}

/*--------------------------------------------------------------------*/
///  QR factorization (LAPACK result)
/**  Computes the Cholesky factorization of a real symmetric
 *   positive definite band matrix A.
 *   The factorization has the form
 *   A = U**T * U,  if a_uploA = 'U', or
 *   A = L  * L**T,  if a_uploA = 'L',
 *   where U is an upper triangular matrix and L is lower triangular.
 *
 *  \param[in]  a_A    Input matrix A
 *  \param[out] a_A    Output matrix A, fully filled
 *  \param[in]  a_uploA
 *                     Flag for if using upper of lower triangle
 *//*-----------------------------------------------------------------*/

void CHpbtrf(CHMatrix& a_A,
             char      a_uploA)
{
  const int fullA = isFull(a_uploA);
  int ier = 0;
  int lfi_N = a_A.size(0);
  // int KD = lfi_N - 1;
  int lfi_LDA = lfi_N;

  //LAPACK_PBTRF(&a_uploA, &lfi_N, &KD, a_A.begin(), &lfi_LDA, &ier);
  LAPACK_POTRF(&a_uploA, &lfi_N, a_A.begin(), &lfi_LDA, &ier);
  if (ier != 0)
    {
      pout() << "PBTRF failure status: " << ier << std::endl;
      MayDay::Error("LAPACK routine @PBTRF failed");
    }
  // Populate upper triangular part of symmetric matrix
  if (fullA)
    {
      sy2ge(lfi_N, a_A);
    } 
}

/*--------------------------------------------------------------------*/
/// Solve a system of linear equations
/** computes the solution to a real system of linear equations
 *  A * X = B,
 *  where A is an N-by-N matrix and X and B are N-by-1 vectorss.
 *  
 *  The LU decomposition with partial pivoting and row interchanges is
 *  used to factor A as
 *      A = P * L * U,
 *  where P is a permutation matrix, L is unit lower triangular, and U is
 *  upper triangular.  The factored form of A is then used to solve the
 *  system of equations A * X = B.
 *
 *  \param[in]  a_A    Input square matrix A
 *  \param[out] a_A    Output matrix with L and U from factorization
 *                     A = P*L*U; The unit diagonal of L is not stored.
 *  \param[in]  a_B    Input vector B
 *  \param[out] a_B    Output vector X
 *  \param[in]  a_ipiv Pivot array
 *//*-----------------------------------------------------------------*/

void CHgesv(CHMatrix& a_A,
            CHVector& a_B,
            CHArray<int, 1>& a_ipiv)
{
  int a_N = a_A.size(0);
  int a_NRHS = 1;
  int lfi_N = a_A.size(0);
  CH_assert(a_A.size(1) >= lfi_N);
  CH_assert(a_A.size(0) == a_A.size(1));
  CH_assert(a_A.size(0) == a_B.size(0));
  int lfi_LDA = a_A.size(0);
  int lfi_LDB = a_A.size(0);
  int ier = 0;
    
  // Setup pivots if required
  if ((!a_ipiv.isAllocated()) || a_ipiv.size() < lfi_N)
    {
      a_ipiv.define(lfi_N);
    }

  LAPACK_GESV(&a_N, &a_NRHS, a_A.begin(), &lfi_LDA, a_ipiv.begin(), a_B.begin(), &lfi_LDB, &ier);
  if (ier != 0)
    {
      pout() << "GESV failure status: " << ier << std::endl;
      MayDay::Error("LAPACK routine @GESV failed");
    }
}

/// Same as above except now B and X are matrices
/*--------------------------------------------------------------------*/
/// Solve multiple systems of linear equations
/** computes the solution to a real system of linear equations
 *  A * X = B,
 *  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
 *  
 *  The LU decomposition with partial pivoting and row interchanges is
 *  used to factor A as
 *      A = P * L * U,
 *  where P is a permutation matrix, L is unit lower triangular, and U is
 *  upper triangular.  The factored form of A is then used to solve the
 *  system of equations A * X = B.
 *
 *  \param[in]  a_A    Input square matrix A
 *  \param[out] a_A    Output matrix with L and U from factorization
 *                     A = P*L*U; The unit diagonal of L is not stored.
 *  \param[in]  a_B    Input matrix B
 *  \param[out] a_B    Output matrix X
 *  \param[in]  a_ipiv Pivot array
 *//*-----------------------------------------------------------------*/

void CHgesv(CHMatrix& a_A,
            CHMatrix& a_B,
            CHArray<int, 1>& a_ipiv)
{
  int lfi_N = a_A.size(1);
  int lfi_NB = a_B.size(0);
  CH_assert(lfi_N == lfi_NB);
  CH_assert(a_A.size(0) == a_A.size(1));
  CH_assert(a_B.size(0) == a_B.size(1));
  CH_assert(a_A.size() == a_B.size());
  int lfi_NRHS = a_B.size(1);
  int lfi_LDA = a_A.size(0);
  int lfi_LDB = a_B.size(0);
  int ier = 0;
  
  // Setup pivots if required
  if ((!a_ipiv.isAllocated()) || a_ipiv.size() < lfi_N)
    {
      a_ipiv.define(lfi_N);
    }

  LAPACK_GESV(&lfi_N, &lfi_NRHS, a_A.begin(), &lfi_LDA, a_ipiv.begin(), a_B.begin(), &lfi_LDB, &ier);
  if (ier != 0)
    {
      pout() << "GESV failure status: " << ier << std::endl;
      MayDay::Error("LAPACK routine @GESV failed");
    }
}

/*--------------------------------------------------------------------*/
///  Solve a system of linear equations
/** computes the solution to a real system of linear equations
 *  A * X = B,
 *  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
 *  
 *  The LU decomposition with partial pivoting and row interchanges is
 *  used to factor A as
 *      A = P * L * U,
 *  where P is a permutation matrix, L is unit lower triangular, and U is
 *  upper triangular.  The factored form of A is then used to solve the
 *  system of equations A * X = B.
 *
 *  \param[in]  a_A    Input matrix A
 *  \param[in]  a_B    Input vector B
 *  \param[out] a_B    Output vector X
 *  \param[in]  a_work Matrix of for working space
 *  \param[in]  a_swork
 *                     Matrix of for working space
 *  \param[in]  a_ipiv Pivot array
 *  \param[in]  a_N    Size of matrix A
 *//*-----------------------------------------------------------------*/

void CHsgesv(CHMatrix& a_A,
             CHVector& a_B,
             CHVector& a_X,
             CHVector& a_work,
             CHVector& a_swork,
             CHArray<int, 1>& a_ipiv,
             int a_N)
{
  int ier = 0;
  int iter = 30;
  int a_NRHS = 1;
  int lfi_N = a_A.size(0);
  CH_assert(a_A.size(1) >= lfi_N);
  int lfi_LDA = lfi_N;
  int lfi_LDB = a_B.size(0);
  int lfi_LDX = a_X.size(0);
  
  
  // Setup pivots if required
  if ((!a_ipiv.isAllocated()) || a_ipiv.size() < lfi_N)
    {
      a_ipiv.define(lfi_N);
    }

  if ((!a_work.isAllocated()) || a_work.size() < lfi_N)
    {
      a_work.define(lfi_N);
    }

  if ((!a_swork.isAllocated()) || a_swork.size() < lfi_N*(a_N + a_NRHS))
    {
      a_swork.define(lfi_N*(a_N + a_NRHS));
    }

  LAPACK_SGESV(&a_N, &a_NRHS, a_A.begin(), &lfi_LDA, a_ipiv.begin(), a_B.begin(), &lfi_LDB, a_X.begin(), &lfi_LDX, a_work.begin(), a_swork.begin(), &iter, &ier);
  if (ier != 0)
    {
      pout() << "SGESV failure status: " << ier << std::endl;
      MayDay::Error("LAPACK routine @SGESV failed");
    }

  if (iter < 0)
    {
      pout() << "SGESV failure status at iteration: " << iter << std::endl;
      MayDay::Error("LAPACK routine @SGESV failed");
    }
}

/*--------------------------------------------------------------------*/
///  Compute eigenvalues of square matrix
/**
 *  The right eigenvector v(j) of A satisfies
 *  A * v(j) = lambda(j) * v(j)
 *  where lambda(j) is its eigenvalue.
 *  The left eigenvector u(j) of A satisfies
 *  u(j)**H * A = lambda(j) * u(j)**H
 *  where u(j)**H denotes the conjugate-transpose of u(j).
 *  
 *  The computed eigenvectors are normalized to have Euclidean norm
 *  equal to 1 and largest component real.
 *
 *  \param[in]  a_A    Input square matrix A
 *  \param[out] a_WR   Real part of the eigenvalues
 *  \param[out] a_WI   Imaginary part of the eigenvalues
 *//*-----------------------------------------------------------------*/

void CHgeev(CHMatrix& a_A,
            CHVector& a_WR,
            CHVector& a_WI)
{
  char a_jobR = 'N';
  char a_jobL = 'N';
  CHMatrix a_VR(1,1);
  CHMatrix a_VL(1,1);
  CHgeev(a_A, a_WR, a_WI, a_VR, a_VL, a_jobR, a_jobL);
}

/*--------------------------------------------------------------------*/
///  Compute eigenvalues and eigenvectors of square matrix
/**
 *  The right eigenvector v(j) of A satisfies
 *  A * v(j) = lambda(j) * v(j)
 *  where lambda(j) is its eigenvalue.
 *  The left eigenvector u(j) of A satisfies
 *  u(j)**H * A = lambda(j) * u(j)**H
 *  where u(j)**H denotes the conjugate-transpose of u(j).
 *  
 *  The computed eigenvectors are normalized to have Euclidean norm
 *  equal to 1 and largest component real.
 *
 *  \param[in]  a_A    Input square matrix A
 *  \param[out] a_WR   Real part of the eigenvalues
 *  \param[out] a_WI   Imaginary part of the eigenvalues
 *  \param[out] a_VR   Right eigenvectors, stored in columns
 *  \param[in]  a_jobR Flag is right eigenvectors should be solved
 *                     'V' to solve, and 'N' to ignore
 *//*-----------------------------------------------------------------*/

void CHgeev(CHMatrix& a_A,
            CHVector& a_WR,
            CHVector& a_WI,
            CHMatrix& a_VR,
            char      a_jobR)
{
  char a_jobL = 'N';
  CHMatrix a_VL(1,1);
  CHgeev(a_A, a_WR, a_WI, a_VR, a_VL, a_jobR, a_jobL);
}

/*--------------------------------------------------------------------*/
///  Compute eigenvalues and eigenvectors of square matrix
/**
 *  The right eigenvector v(j) of A satisfies
 *  A * v(j) = lambda(j) * v(j)
 *  where lambda(j) is its eigenvalue.
 *  The left eigenvector u(j) of A satisfies
 *  u(j)**H * A = lambda(j) * u(j)**H
 *  where u(j)**H denotes the conjugate-transpose of u(j).
 *  
 *  The computed eigenvectors are normalized to have Euclidean norm
 *  equal to 1 and largest component real.
 *
 *  \param[in]  a_A    Input square matrix A
 *  \param[out] a_A    Overwritten as workspace
 *  \param[out] a_WR   Real part of the eigenvalues
 *  \param[out] a_WI   Imaginary part of the eigenvalues
 *  \param[out] a_VR   Right eigenvectors, stored in columns
 *  \param[out] a_VL   Left eigenvectors, stored in columns
 *  \param[in]  a_jobR Flag is right eigenvectors should be solved
 *  \param[in]  a_jobL Flag is left eigenvectors should be solved
 *                     'V' to solve, and 'N' to ignore
 *//*-----------------------------------------------------------------*/

void CHgeev(CHMatrix& a_A,
            CHVector& a_WR,
            CHVector& a_WI,
            CHMatrix& a_VR,
            CHMatrix& a_VL,
            char      a_jobR,
            char      a_jobL)
{
  const int isJobR = isJobz(a_jobR);
  const int isJobL = isJobz(a_jobL);
  int ier = 0;
  int lfi_LDA = a_A.size(0);
  int lfi_N = a_A.size(1);
  int lfi_LDVR = isJobR ? a_VR.size(0) : 1;
  int lfi_LDVL = isJobL ? a_VL.size(0) : 1;
  CH_assert(a_A.size(0) == a_A.size(1));
  CH_assert(a_WR.size(0) >= lfi_N);
  CH_assert(a_WI.size(0) >= lfi_N);
  CH_assert(a_VR.size(0) >= lfi_N*isJobR);
  CH_assert(a_VL.size(0) >= lfi_N*isJobL);
  
  int lfi_LWORK = 5*lfi_N;
  CHVector work(lfi_LWORK);
  
  LAPACK_GEEV(&a_jobL, &a_jobR, &lfi_N, a_A.begin(), &lfi_LDA, a_WR.begin(), a_WI.begin(), a_VL.begin(), &lfi_LDVL, a_VR.begin(), &lfi_LDVR, work.begin(), &lfi_LWORK, &ier);
  
  if (ier != 0)
    {
      pout() << "GEEV failure status: " << ier << std::endl;
      MayDay::Error("LAPACK routine @SYEV failed");
    }
}

/*--------------------------------------------------------------------*/
///  Computes the singular value decomposition of a real M-by-N matrix
/**
 *  Compute SVD, optionally computing the left and/or right singular
 *  vectors. The SVD is written
 *     A = U * SIGMA * transpose(V)
 *  where SIGMA is an M-by-N matrix which is zero except for its
 *  min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
 *  V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
 *  are the singular values of A; they are real and non-negative, and
 *  are returned in descending order.  The first min(m,n) columns of
 *  U and V are the left and right singular vectors of A.
 *  
 *  Note that the routine returns V**T, not V.
 *
 *  \param[in]  a_A    Input matrix A
 *  \param[out] a_A    if JOBU = 'O', A is overwritten with the first
 *                      min(m,n) columns of U (the left singular vectors);
 *                     if JOBVT = 'O', A is overwritten with the first
 *                      min(m,n) rows of V**T (the right singular vectors);
 *                     if JOBU != 'O' and JOBVT != 'O', the contents of A
                          are destroyed.
 *  \param[out] a_U    If JOBU = 'A', U contains M*M orthogonal matrix U;
 *                     if JOBU = 'S', U contains the first min(m,n)
 *                      columns of U (the left singular vectors);
 *                     if JOBU = 'N' or 'O', U is not referenced.
 *  \param[out] a_S    Singular values of A, sorted so S(i) >= S(i+1)
 *  \param[in]  a_jobU Flag for how to fill U
 *  \param[in]  a_jobVT
 *                     Flag for how to fill V^t
 *//*-----------------------------------------------------------------*/

void CHgesvd(CHMatrix& a_A,
             CHMatrix& a_U,
             CHMatrix& a_S,
             char      a_jobU,
             char      a_jobVT)
{
  int ier = 0;
  int m = a_A.size(0);
  int n = a_A.size(1);
  int min = 0;
  if (m <= n) min=m;
  else min=n;
  CHMatrix a_VT(n, n);
  int LVDT = n;
  CHVector S(min);
  int LWORK = 5*(m+n);
  CHVector WORK(LWORK);
 
  LAPACK_GESVD(&a_jobU, &a_jobVT, &m, &n, a_A.begin(), &m, S.begin(), a_U.begin(), &m, a_VT.begin(), &LVDT, WORK.begin(), &LWORK, &ier);
  //put S in matrix form
  a_S = 0.0;
  for(int i = 0; i!=min; i++){
    a_S(i, i) = S(i);
  }

  if (ier != 0)
    {
      pout() << "GESVD failure status: " << ier << std::endl;
      MayDay::Error("LAPACK routine @GESVD failed");
    }
}

/*--------------------------------------------------------------------*/
/// Compute the solution to a Banded system of linear equation A*x=b
/** matix looks as follows, with kl lower diagonals, and ku upper diagonals
 * As an example with kl = ku = 1
 *    a11 a12                  *   *   *   *   *
 *    a21 a22 a23              *   a12 a23 a34 a45
 *        a32 a33 a34      ->  a11 a22 a33 a44 a55
 *            a43 a44 a45      a21 a32 a43 a54 *
 *                a54 a55
 *                             must be size (2*kl + ku + 1)
 *                             Stars are placeholders needed, more are okay
 *
 *  \param[in]  a_A    Input banded matrix A, with above format
 *  \param[out] a_A    Output factorized matrix
 *  \param[in]  a_kl   Width of lower band
 *  \param[in]  a_ku   Width of upper band
 *  \param[in]  a_b    Input right hand side vector b
 *  \param[out] a_b    Output vector x
 *//*-----------------------------------------------------------------*/

void CHgbsv(CHMatrix& a_A,
            int a_kl,
            int a_ku,
            CHVector& a_b)
{
  int ier = 0;
  int N = a_A.size(1);
  int NRHS = 1;
  int LDAB = a_A.size(0);
  CHArray<int, 1> ipiv(N);
  int LDB = a_b.size(0);
  CH_assert(a_kl >= 0);
  CH_assert(a_ku >= 0);
  CH_assert(a_b.size(0) == a_A.size(1));
  CH_assert(LDAB >= 2*a_kl+a_ku+1);
  
  LAPACK_GBSV(&N, &a_kl, &a_ku, &NRHS, a_A.begin(), &LDAB, ipiv.begin(), a_b.begin(), &LDB, &ier);
  
  if (ier != 0)
    {
      pout() << "GBSV failure status: " << ier << std::endl;
      MayDay::Error("LAPACK routine @GBSV failed");
    }
}

/*----------------------------------------------------------------------------*/
/// Compute the least squares solution to the possibly deficient problem A*x = b
/** With A an mxn matrix, b the nx1 rhs vector, and x the nx1 solution
 *  vector. In this problem m >= n.
 *
 *  \param[in]  a_A      Input matrix A
 *  \param[in]  a_b      The RHS of the matrix equation Ax = b
 *  \param[out] a_b      The solution x of LLS problem Ax = b
 *  \param[in]  a_trans  Is the matrix A transpose. 'N' for not (default),
 *                       'T' for A involves the transpose.
 *//*-------------------------------------------------------------------------*/
void CHgels(CHMatrix& a_A,
            CHVector& a_x,
            CHVector& a_b,
            char      a_trans)
{
  CH_assert(a_x.size(0) == a_A.size(1));
  CH_assert(a_b.size(0) == a_A.size(0));

  int M    = a_A.size(0);          // number of rows of A
  int N    = a_A.size(1);          // number of columns of A
  int NRHS = a_b.size(1);          // number of columns of b
  int LDA  = a_A.size(0);          // leading dimension of A
  int LDB  = a_b.size(0);          // leading dimension of b

  int MN = std::min(M,N);
  int LWORK = MN + std::max(MN,NRHS);
  CHVector WORK (LWORK);
  WORK = 0.0;
  int ier = 0;

  LAPACK_GELS(&a_trans, &M, &N, &NRHS, a_A.begin(), &LDA,
              a_b.begin(), &LDB, WORK.begin(), &LWORK, &ier);

  // Copy the result of the computation a_b back to a_x
  Real* xptr = a_x.begin();
  Real* bptr = a_b.begin();
  for(; xptr != a_x.end(); ++xptr, ++bptr)
    {
      *xptr = *bptr;
    }

  if(ier != 0)
    {
      pout() << "GELS failure status: " << ier << std::endl;
      MayDay::Error("LAPACK routine @GELS failed");
    }
}

#ifdef CH_SPACEDIM // In Chombo land
#include "NamespaceFooter.H"
#endif
