C      _______              __
C     / ___/ /  ___  __ _  / /  ___
C    / /__/ _ \/ _ \/  V \/ _ \/ _ \
C    \___/_//_/\___/_/_/_/_.__/\___/
C    Please refer to Copyright.txt, in Chombo's root directory.

#include "REAL.H"
#include "CONSTANTS.H"

      subroutine solveqr(neqns, nvars, A, Xt, findc, cond)
!
!     Given matrix A(1:neqns, 1:nvars), this subroutine finds
!     matrices Q(1:neqns, 1:neqns) and R(1:neqns, 1:nvars)
!     such that A = Q * R
!     and then solves for X(1:nvars, 1:neqns) in R * X = Q^T.
!     Then A * X = Q * R * X = Q * Q^T = I
!     where I has length neqns.
!      
      implicit none
!     --- arguments ---
      integer neqns, nvars
      real_t A(neqns, nvars)
      real_t Xt(neqns, nvars)
      integer findc
      real_t cond
!     --- locals ---
      integer info, i, j
      character*1 prec
      integer lwork
      real_t tau(nvars), work(nvars*neqns)
!     --- locals used if and only if findc .eq. 1
      real_t condwk(3*nvars), rcond
      integer iwork(nvars)

#ifdef CH_USE_DOUBLE
      prec = 'D'
#else
      prec = 'S'
#endif

!     DGEQRF finds the QR decomposition of A,
!     and stores the result in a rearranged A and tau.
!     For optimum performance LWORK >= N*NB, where NB is the optimal blocksize.
!     I just set NB to neqns so we'll have some slack.
      lwork = nvars * neqns
#ifdef CH_USE_DOUBLE
      call DGEQRF
#else
      call SGEQRF
#endif
     &     (neqns, nvars, A, neqns, tau, work, lwork, info)
      if (info .lt. 0) then
         print *, 'Error:  LAPACK ', prec, 'GEQRF returned ', info
         print *, 'Illegal value in argument number ', -info
         call MAYDAY_ERROR()
      endif

      if (findc .eq. 1) then
!     DTRCON finds the reciprocal of the condition number of R (hence of A).
!     '1' for 1-norm;
!     'U' for upper triangular matrix;
!     'N' for not assuming matrix to be unit triangular.
#ifdef CH_USE_DOUBLE
         call DTRCON
#else
         call STRCON
#endif
     &        ('1', 'U', 'N', nvars,
     &        A, neqns, rcond, condwk, iwork, info)
         if (info .lt. 0) then
            print *, 'Error:  LAPACK ', prec, 'TRCON returned ', info
            print *, 'Illegal value in argument number ', -info
            call MAYDAY_ERROR()
         endif
         cond = one / rcond
      endif

!     Set Xt = A.
      do i = 1, neqns
         do j = 1, nvars
            Xt(i, j) = A(i, j)
         enddo
      enddo

!     DORGQR generates the real matrix Q(1:neqns, nvars),
!     and stores it in Xt.
#ifdef CH_USE_DOUBLE
      call DORGQR
#else
      call SORGQR
#endif
     &     (neqns, nvars, nvars, Xt, neqns, tau, work, lwork, info)
      if (info .lt. 0) then
         print *, 'Error:  LAPACK ', prec, 'ORGQR returned ', info
         print *, 'Illegal value in argument number ', -info
         call MAYDAY_ERROR()
      endif

!     DTRSM solves for X^T in X^T * R^T = one * Q.
!     (This is how we solve for X in R * X = Q^T.)
!     'R' for multiplication on the right;
!     'U' for upper triangular matrix;
!     'T' for multiplication by matrix transpose;
!     'N' for not assuming matrix to be unit triangular.
#ifdef CH_USE_DOUBLE
      call DTRSM
#else
      call STRSM
#endif
     &     ('R', 'U', 'T', 'N', neqns, nvars,
     &     one, A, neqns, Xt, neqns)

      return
      end
