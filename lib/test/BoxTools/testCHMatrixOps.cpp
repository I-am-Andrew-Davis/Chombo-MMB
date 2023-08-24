#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

/*------------------------------------------------------------------------------
 * This file runs test to make sure calls to Lapack/Blas work as expected
 *----------------------------------------------------------------------------*/

#include <iostream>
#include <cstring>
using std::endl;
#include "parstream.H"
#include "CHMatrixOps.H"
#ifdef CH_MPI
#include <mpi.h>
#endif
#include "UsingNamespace.H"

/// Prototypes:

void
parseTestOptions( int argc ,char* argv[]);

bool
matrixEqual(const CHMatrix& a_A, const CHMatrix& a_B);

bool
vectorEqual(const CHVector& a_A, const CHVector& a_B);

int
testCHMatrixOps();

/// Global variables for handling output:
static const char *pgmname = "testCHMatrixOps";
static const char *indent = "   ", *indent2 = "      " ;
static bool verbose = true ;

/// Code:

int
main(int argc, char* argv[])
{
#ifdef CH_MPI
  MPI_Init(&argc, &argv);
#endif
  parseTestOptions( argc ,argv ) ;

  if (verbose)
    pout() << indent2 << "Beginning " << pgmname << " ..." << endl ;

  ///
  // Run the tests
  ///
  int ret = testCHMatrixOps() ;

  if (ret == 0)
    {
      pout() << indent << pgmname << " passed all tests" << endl ;
    }
  else
    {
      pout() << indent << pgmname << " failed " << ret << " test(s)" << endl ;
    }

#ifdef CH_MPI
  MPI_Finalize();
#endif
  return ret;
}


/*******************************************************************************
 *
 * Routine testCHMatrixOps
 *
 ******************************************************************************/

int
testCHMatrixOps()
{
  int errors = 0;

/*
 * Parameters to use an test inputs
 * Default matrix values are zero
 */
  CHMatrix A3(3,3);
  CHMatrix B3(3,3);
  CHMatrix I3(3,3); // an identity matrix
  CHMatrix Z3(3,3); // a zero matrix
  CHVector z3(3);   // a zero vector
  CHVector one3(3); // vector of ones
  one3 = 1;
  CHVector x3(3);
  for (int r = 0; r!= A3.size(0); ++r)
    {
      for (int c = 0; c!= A3.size(1); ++c)
        {
          A3(r, c) = (1 + c) + r*3;
          B3(r, c) = (1 + r) + c*3;
        }
      I3(r, r) = 1;
      x3(r) = (A3.size(0) - r)*(A3.size(0) - r);
    }

  CHMatrix A3inv(3,3);
  CHMatrix R3(3,3);
  CHVector r3(3);
  // working vars
  CHArray<int, 1> ipiv;
  CHVector work;
/*
 * Test Matrix transpose
 */
  {
    int err = 0;
    CHMatrix transp(3,3);
    
    transpose(transp, A3);
    if (!matrixEqual(B3, transp)) err++;
    
    transpose(transp);
    if (!matrixEqual(A3, transp)) err++;
    
    if (verbose)
      {
        pout() << indent2 << "Test transpose ";
        if (err == 0)
          {
            pout() << "passed" << endl;
          }
        else
          {
            pout() << "failed, with " << err << " errors" << endl;
            ++errors;
          }
      }
  }

/*
 * Test Matrix determinate
 */
  {
    int err = 0;
    const int sz = 3;
    Real dI = determinate(I3, sz);
    if (dI != 1) err++;

    Real dA = determinate(A3, sz);
    Real dB = determinate(B3, sz);
    if (dA != dB) err++;
    if (dA != 0) err++;
    
    if (verbose)
      {
        pout() << indent2 << "Test determinate ";
        if (err == 0)
          {
            pout() << "passed" << endl;
          }
        else
          {
            pout() << "failed, with " << err << " errors" << endl;
            ++errors;
          }
      }
  }
  
/*
 * Test Matrix-Vector multiplication (CHgemv)
 */
  {
    int err = 0;

    // multiply by zero matrix
    CHgemv(Z3, x3, r3);
    if (!vectorEqual(r3, z3)) err++;

    // multiply by zero vector
    CHgemv(A3, z3, r3);
    if (!vectorEqual(r3, z3)) err++;
    
    // multiply by identity matrix
    CHgemv(I3, x3, r3);
    if (!vectorEqual(x3, r3)) err++;

    // multiply an arbitrary matrix and vector
    CHgemv(A3, x3, r3);
    CHVector AX(3); // solution vector
    AX(0) = 9*1 + 4*2 + 1*3;
    AX(1) = 9*4 + 4*5 + 1*6;
    AX(2) = 9*7 + 4*8 + 1*9;
    if (!vectorEqual(AX, r3)) err++;
    
    if (verbose)
      {
        pout() << indent2 << "Test CHgemv ";
        if (err == 0)
          {
            pout() << "passed" << endl;
          }
        else
          {
            pout() << "failed, with " << err << " errors" << endl;
            ++errors;
          }
        // show multiplication
        pout() << "A*x=b\n" << A3 << " * \n" << x3 << "\n = \n" << r3 << endl;
      }
  }

/*
 * Test Matrix-Matrix multiplication (CHgemm)
 */
  {
    int err = 0;

    // multiply by zero matrix
    CHgemm(Z3, A3, R3);
    if (!matrixEqual(R3, Z3)) err++;

    CHgemm(A3, Z3, R3);
    if (!matrixEqual(R3, Z3)) err++;
    
    // multiply by identity matrix
    CHgemm(I3, A3, R3);
    if (!matrixEqual(A3, R3)) err++;

    CHgemm(A3, I3, R3);
    if (!matrixEqual(A3, R3)) err++;
    
    // multiply an arbitrary matrix and vector
    CHgemm(A3, B3, R3);
    CHMatrix AX(3, 3); // solution
    AX(0, 0) = 1*1 + 2*2 + 3*3;
    AX(1, 0) = 4*1 + 5*2 + 6*3;
    AX(2, 0) = 7*1 + 8*2 + 9*3;
    AX(0, 1) = 1*4 + 2*5 + 3*6;
    AX(1, 1) = 4*4 + 5*5 + 6*6;
    AX(2, 1) = 7*4 + 8*5 + 9*6;
    AX(0, 2) = 1*7 + 2*8 + 3*9;
    AX(1, 2) = 4*7 + 5*8 + 6*9;
    AX(2, 2) = 7*7 + 8*8 + 9*9;
    if (!matrixEqual(AX, R3)) err++;
    
    if (verbose)
      {
        pout() << indent2 << "Test CHgemm ";
        if (err == 0)
          {
            pout() << "passed" << endl;
          }
        else
          {
            pout() << "failed, with " << err << " errors" << endl;
            ++errors;
          }
        // show multipication
        pout() << "A*B=C\n" << A3 << " * \n" << B3 << " = \n" << R3 << endl;
      }
  }

/*
 * Test Matrix Inverse
 */
  {
    int err = 0;
    A3(2,2) = 0; // alter A3 to be invertable
    // copy A3 to A3inv
    for (int r = 0; r!= A3.size(0); ++r)
      {
        for (int c = 0; c!= A3.size(1); ++c)
          {
            A3inv(r, c) = A3(r, c);
          }
      }
    // solve the inverse
    CHinverse(A3inv, ipiv, work);
    // check A * Ainv = I
    CHgemm(A3, A3inv, R3);
    if (!matrixEqual(R3, I3)) err++;
     
    if (verbose)
      {
        pout() << indent2 << "Test inverse ";
        if (err == 0)
          {
            pout() << "passed" << endl;
          }
        else
          {
            pout() << "failed, with " << err << " errors" << endl;
            ++errors;
          }
        // show inverse
        pout() << "A*A^-1=I\n" << A3 << " * \n" << A3inv << endl;
      }
  }
  
/*
 * Test Ax = b solution (CHgesv)
 */
  {
    int err = 0;

    // copy x3 into r3, and A3 in tmpA3
    CHMatrix tmpA3(3,3);
    for (int r = 0; r!= A3.size(0); ++r)
      {
        for (int c = 0; c!= A3.size(1); ++c)
          {
            tmpA3(r, c) = A3(r, c);
          }
        r3(r) = x3(r);
      }
    // solve for r in A*r = x
    CHgesv(tmpA3, r3, ipiv);
    
    // check, r = Ainv*x
    CHVector tmpr3(3);
    CHgemv(A3inv, x3, tmpr3);
    if (!vectorEqual(r3, tmpr3)) err++;
     
    if (verbose)
      {
        pout() << indent2 << "Test CHgesv ";
        if (err == 0)
          {
            pout() << "passed" << endl;
          }
        else
          {
            pout() << "failed, with " << err << " errors" << endl;
            ++errors;
          }
        // show solution
        pout() << "A*x=b\n" << A3 << " * \n" << r3 << "\n = \n" << x3 << endl;
      }
  }
  
/*
 * Test Ax = b solution for banded matrix (CHgbsv)
 */
  {
    int err = 0;

    // Create a banded matrix to solve, use the classic B(1,-2,1)
    // test both the full, and banded structure
    int sz = 5, kl = 1, ku = 1, l = 2*kl+ku+1;
    CHMatrix fullB(sz, sz), B(l, sz), tmpB(sz, sz);
    CHVector bb(sz), xb(sz), xFull(sz);
    bb(0) = 1;
    bb(sz-1) = 2;
    xb(0) = 1;
    xb(sz-1) = 2;
    xFull(0) = 1;
    xFull(sz-1) = 2;
    for (int r = 0; r!= sz; ++r)
      {
        // fill the full matrix
        fullB(r, r) = -2;
        tmpB(r, r) = -2;
        if (r+1 < sz)
          {
            fullB(r, r+1) = 1;
            fullB(r+1, r) = 1;
            tmpB(r, r+1) = 1;
            tmpB(r+1, r) = 1;
          }
        // fill banded
        B(l-1-kl, r) = -2;
        B(l-1-kl-ku, r) = 1;
        B(l-1, r) = 1;
      }
    // solve for x in A*x = b, using banded solver
    CHgbsv(B, kl, ku, xb);
    // solve for x in A*x = b, full solver
    CHgesv(tmpB, xFull, ipiv);
    //check solution match
    if (!vectorEqual(xb, xFull)) err++;
     
    if (verbose)
      {
        pout() << indent2 << "Test CHgbsv ";
        if (err == 0)
          {
            pout() << "passed" << endl;
          }
        else
          {
            pout() << "failed, with " << err << " errors" << endl;
            ++errors;
          }
        // show solution
        pout() << "A*x=b\n" << fullB << " * \n" << xb << "\n = \n" << bb << endl;
      }
  }

/*
 * Test eigenvalues and vectors of matrix (CHgeev)
 */
  {
    int err = 0;
    
    CHVector wr3(3), wi3(3);
    CHMatrix Vr3(3,3), Vl3(3,3);
    CHMatrix tmpAV3(3,3), tmpwV3(3,3);
    // copy matrices
    CHMatrix tmpA3(3,3), tmpI3(3,3);
    tmpI3=0;
    for (int r = 0; r!= A3.size(0); ++r)
      {
        for (int c = 0; c!= A3.size(1); ++c)
          {
            tmpA3(r, c) = A3(r, c);
          }
        tmpI3(r, r) = 1;
      }
    
    // test an identity matrix
    CHgeev(tmpI3, wr3, wi3, Vr3, Vl3);
    if (!vectorEqual(wr3, one3)) err++;
    if (!vectorEqual(wi3, z3)) err++;
    if (!matrixEqual(Vr3, I3)) err++;
    if (!matrixEqual(Vl3, I3)) err++;

    // test the A matrix
    CHgeev(tmpA3, wr3, wi3, Vr3, Vl3);
    if (!vectorEqual(wi3, z3)) err++; // eigenvalues are real
    // diagonal matrix of eigenvalues
    CHMatrix Wr3(3,3);
    for (int r = 0; r!= A3.size(0); ++r)
      {
        Wr3(r,r) = wr3(r);
      }
    // check that A*Vr = Vr*W
    CHgemm(A3, Vr3, tmpAV3);
    CHgemm(Vr3, Wr3, tmpwV3);
    if (!matrixEqual(tmpAV3, tmpwV3)) err++;
    // check that Vl^T*A = Vl*W
    transpose(Vl3);
    CHgemm(Vl3, A3, tmpAV3);
    CHgemm(Wr3, Vl3, tmpwV3);
    if (!matrixEqual(tmpAV3, tmpwV3)) err++;

    // test the A matrix, for just the eigenvalues
    CHgeev(tmpA3, wr3, wi3);
    if (!vectorEqual(wi3, z3)) err++; // eigenvalues are real
    
    if (verbose)
      {
        pout() << indent2 << "Test CHgeev ";
        if (err == 0)
          {
            pout() << "passed" << endl;
          }
        else
          {
            pout() << "failed, with " << err << " errors" << endl;
            ++errors;
          }
      }
  }

  /*
   * Test Least Squares solve (CHgels)
   */
  {
    int err = 0;

    CHVector ones3(3), zeros3(3), V3(3);
    CHVector ones2(2), zeros2(2), V2(2);

  }


  return errors;
}


/*--------------------------------------------------------------------*
 * Check if two matrices are the same
 *--------------------------------------------------------------------*/
bool
matrixEqual(const CHMatrix& a_A, const CHMatrix& a_B)
{
  const Real eps = std::numeric_limits<Real>::epsilon();
  const Real tol = eps*std::numeric_limits<Real>::digits*std::numeric_limits<Real>::digits;
  if (a_A.size() != a_B.size())
    {
      return false;
    }
  for (int r = 0; r != a_A.size(0); ++r)
    {
      for (int c = 0; c != a_A.size(0); ++c)
        {
          // this should be a better check for floats
          // this only works for scales greater than 1
          if ((std::abs(a_A(r,c) - a_B(r,c))) >=
              (tol*std::max(1.0, std::abs(a_A(r,c) + a_B(r,c)))))
            {
              return false;
            }
        }
    }
  return true;
}

/*--------------------------------------------------------------------*
 * Check if two vectors are the same
 *--------------------------------------------------------------------*/
bool
vectorEqual(const CHVector& a_A, const CHVector& a_B)
{
  const Real eps = std::numeric_limits<Real>::epsilon();
  const Real tol = eps*std::numeric_limits<Real>::digits*std::numeric_limits<Real>::digits;
  if (a_A.size() != a_B.size())
    {
      return false;
    }
  for (int r = 0; r != a_A.size(0); ++r)
    {
      // this should be a better check for floats
      // this only works for scales greater than 1
      if ((std::abs(a_A(r) - a_B(r))) >=
          (tol*std::max(1.0, std::abs(a_A(r) + a_B(r)))))
        {
          return false;
        }
    }
  return true;
}

/*--------------------------------------------------------------------*
 * Parse the standard test options (-v -q) out of the command line.
 * Stop parsing when a non-option argument is found.
 *--------------------------------------------------------------------*/
void
parseTestOptions( int argc ,char* argv[] )
{
  for ( int i = 1 ; i < argc ; ++i )
    {
      if (argv[i][0] == '-') //if it is an option
        {
          // compare 3 chars to differentiate -x from -xx
          if (strncmp( argv[i] ,"-v" ,3 ) == 0)
            {
              verbose = true ;
              // argv[i] = "" ;
            }
          else if (strncmp( argv[i] ,"-q" ,3 ) == 0)
            {
              verbose = false ;
              // argv[i] = "" ;
            }
          else
            {
              break ;
            }
        }
    }
  return ;
}
