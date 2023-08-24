C      _______              __
C     / ___/ /  ___  __ _  / /  ___
C    / /__/ _ \/ _ \/  V \/ _ \/ _ \
C    \___/_//_/\___/_/_/_/_.__/\___/
C    Please refer to Copyright.txt, in Chombo's root directory.

#include "REAL.H"
#include "CONSTANTS.H"

      subroutine solvenormalequations(neqns, nvars, at, findc, cond)
!
!     Given matrix A(1:neqns, 1:nvars), this subroutine finds
!     M(1:nvars, 1:neqns) satisfying the normal equations
!     (A^T * A) * M = A^T.
!
!     Why these equations:  Our goal is to find a least-squares
!     solution, x(1:nvars), to the overdetermined system
!     A * x = b
!     given A(1:neqns, 1:nvars) and b(1:neqns).
!     That is, find x minimizing ||A * x - b||.
!     We do this by solving for x in:
!     A^T * A * x = A^T * b.
!     Solution is x = inv(A^T * A) * A^T * b.
!     This is x = M * b if we set
!     M = inv(A^T * A) * A^T.
!     Another way to write this is:
!     (A^T * A) * M = A^T.
!     We solve for M(1:nvars, 1:neqns) because we will be
!     using different values of b.
!
!     Also, we start with A^T(1:nvars, 1:neqns), because in filling
!     in the matrix, we change variables (coefficients) faster
!     than equations (neighbors).
      
      implicit none
!     --- arguments ---
      integer neqns, nvars
      integer findc
      real_t at(nvars, neqns)
      real_t cond
!     --- locals ---
      real_t ata(nvars, nvars)
      integer info
      character*1 prec
!     --- locals used iff findc .eq. 1
      real_t af(nvars, nvars), scales(nvars), soln(nvars, neqns), rcond
      real_t ferr(neqns), berr(neqns), work(3*nvars)
      integer iwork(nvars)
      integer ivar, ieqn
      character*1 equed

#ifdef CH_USE_DOUBLE
      prec = 'D'
#else
      prec = 'S'
#endif

!     Multiply A^T by A.

!     First matrix is NOT transposed (A^T), but second matrix (A) is.
!     Number of rows of A^T is nvars, same as number of rows of A^T*A.
!     Number of columns of A is nvars, same as number of columns of A^T*A.
!     Number of columns of A^T is neqns, same as number of rows of A.
#ifdef CH_USE_DOUBLE
      call DGEMM
#else
      call SGEMM
#endif
     &     ('N', 'T', nvars, nvars, neqns, one,
     &     at, nvars,
     &     at, nvars,
     &     zero, ata, nvars)

!     Solve for M in (A^T*A)*M = A^T.
!     We overwrite the array "at" with M.

!     Use upper triangular part of ata.
!     Matrix ata is nvars by nvars.
!     Right-hand side A^T has neqns columns.
      if (findc .ne. 1) then
#ifdef CH_USE_DOUBLE
         call DPOSV
#else
         call SPOSV
#endif
     &        ('U', nvars, neqns, ata, nvars, at, nvars, info)

         if (info .lt. 0) then
            print *, 'Error:  LAPACK ', prec, 'POSV returned ', info
            print *, 'Illegal value in argument number ', -info
            call MAYDAY_ERROR()
         elseif (info .gt. 0) then
            print *, 'Error:  LAPACK ', prec, 'POSV returned ', info
            print *, 'Singular matrix, leading minor of order ', info
            call MAYDAY_ERROR()
         endif
      else
!     findc .eq. 1

!         print *, 'N = ', nvars
!         print *, 'NRHS = ', neqns
!         print *, 'LDA = ', nvars
!         print *, 'LDAF = ', nvars
!         print *, 'LDB = ', nvars
!         print *, 'LDX = ', nvars

#ifdef CH_USE_DOUBLE
         call DPOSVX
#else
         call SPOSVX
#endif
     &        ('N', 'U', nvars, neqns, ata, nvars, 
     &        af, nvars, equed, scales, 
     &        at, nvars, 
     &        soln, nvars, rcond, ferr, berr, work, iwork, info)

         if (info .lt. 0) then
            print *, 'Error:  LAPACK ', prec, 'POSVX returned ', info
            print *, 'Illegal value in argument number ', -info
            call MAYDAY_ERROR()
         elseif (info .eq. (nvars + 1)) then
            print *, 'Error:  LAPACK ', prec, 'POSVX returned ', info
            print *, 'Sub-machine-precision RCOND = ', rcond
            call MAYDAY_ERROR()
         elseif (info .gt. 0) then
            print *, 'Error:  LAPACK ', prec, 'POSVX returned ', info
            print *, 'Singular matrix, leading minor of order ', info
            call MAYDAY_ERROR()
         endif
!     Copy "soln" to "at".
         do ieqn = 1, neqns
            do ivar = 1, nvars
               at(ivar, ieqn) = soln(ivar, ieqn)
            enddo
         enddo
         cond = one / rcond 
      endif

      return
      end
