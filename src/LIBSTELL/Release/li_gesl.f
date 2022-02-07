/* Copyright (C) 1991-2021 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <https://www.gnu.org/licenses/>.  */


/* This header is separate from features.h so that the compiler can
   include it implicitly at the start of every compilation.  It must
   not itself include <features.h> or any other header that includes
   <features.h> because the implicit include comes before any feature
   test macros that may be defined in a source file before it first
   explicitly includes a system header.  GCC knows the name of this
   header in order to preinclude it.  */

/* glibc's intent is to support the IEC 559 math functionality, real
   and complex.  If the GCC (4.9 and later) predefined macros
   specifying compiler intent are available, use them to determine
   whether the overall intent is to support these features; otherwise,
   presume an older compiler has intent to support these features and
   define these macros by default.  */



/* wchar_t uses Unicode 10.0.0.  Version 10.0 of the Unicode Standard is
   synchronized with ISO/IEC 10646:2017, fifth edition, plus
   the following additions from Amendment 1 to the fifth edition:
   - 56 emoji characters
   - 285 hentaigana
   - 3 additional Zanabazar Square characters */

      SUBROUTINE sgesl1 (A, LDA, N, IPVT, B, JOB)
      USE LIPREC, ONLY: WP => SP
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER, INTENT(IN) :: LDA, N, JOB
      INTEGER, DIMENSION(N), INTENT(IN) :: IPVT
      REAL(WP), DIMENSION(*), INTENT(IN) :: A
      REAL(WP), DIMENSION(N), INTENT(INOUT) :: B

      CALL sgesl(A, LDA, N, IPVT, B, JOB)

      END SUBROUTINE sgesl1

      SUBROUTINE dgesl1 (A, LDA, N, IPVT, B, JOB)
      USE LIPREC, ONLY: WP => DP
      IMPLICIT NONE
      
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER, INTENT(IN) :: LDA, N, JOB
      INTEGER, DIMENSION(N), INTENT(IN) :: IPVT
      REAL(WP), DIMENSION(*), INTENT(IN) :: A
      REAL(WP), DIMENSION(N), INTENT(INOUT) :: B

      CALL dgesl(A, LDA, N, IPVT, B, JOB)

      END SUBROUTINE dgesl1

      SUBROUTINE sgesl (A, LDA, N, IPVT, B, JOB)
      USE LIPREC, ONLY: WP => SP
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER, INTENT(IN) :: LDA, N, JOB
      INTEGER, DIMENSION(N), INTENT(IN) :: IPVT
      REAL(WP), DIMENSION(LDA,N), INTENT(IN) :: A
      REAL(WP), DIMENSION(N), INTENT(INOUT) :: B
C-----------------------------------------------
C   L O C A L   V A R I A B L E S
C-----------------------------------------------
      INTEGER :: K, L, NM1
      REAL(WP) :: ELEMENT
C-----------------------------------------------
c
c     sgesl solves the real system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by sgeco or sgefa.
c
c     on ENTRY
c
c        a       REAL(lda, n)
c                the output from sgeco or sgefa.
c
c        lda     INTEGER
c                the leading dimension of the array  a .
c
c        n       INTEGER
c                the order of the matrix  a .
c
c        ipvt    INTEGER(n)
c                the pivot vector from sgeco or sgefa.
c
c        b       REAL(n)
c                the right hand side vector.
c
c        job     INTEGER
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b  where
c                            trans(a)  is the transpose.
c
c     on RETURN
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if sgeco has set rcond.gt.0.0
c        or sgefa has set info.eq.0 .
c
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c
      nm1 = n - 1
      IF (job .eq. 0) THEN
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         IF (nm1 .ge. 1) THEN
            DO k = 1, nm1
               l = ipvt(k)
               element = b(l)
               IF (l .ne. k) THEN
                  b(l) = b(k)
                  b(k) = element
               END IF
!              CALL bla_axpy(n-k,element,a(k+1,k),1,b(k+1),1)
               b(k+1:n) = b(k+1:n) + element*a(k+1:n,k)
            END DO
         END IF
c
c        now solve  u*x = y
c
         DO k = n, 1, -1
            b(k) = b(k)/a(k,k)
            element = -b(k)
!           CALL bla_axpy(k-1,element,a(1,k),1,b,1)
            b(1:k-1) = b(1:k-1) + element*a(1:k-1,k)
         END DO
      ELSE
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         DO k = 1, n
            element = SUM(a(:k-1,k)*b(:k-1))
            b(k) = (b(k)-element)/a(k,k)
         END DO
c
c        now solve trans(l)*x = y
c
         IF (nm1 .ge. 1) THEN
            DO k = nm1, 1, -1
               b(k) = b(k) + SUM(a(k+1:n,k)*b(k+1:n))
               l = ipvt(k)
               IF (l .ne. k) THEN
                  element = b(l)
                  b(l) = b(k)
                  b(k) = element
               END IF
            END DO
         END IF
      END IF

      END SUBROUTINE sgesl

      SUBROUTINE dgesl (A, LDA, N, IPVT, B, JOB)
      USE LIPREC, ONLY: WP => DP
      IMPLICIT NONE
C-----------------------------------------------
C   D U M M Y   A R G U M E N T S
C-----------------------------------------------
      INTEGER, INTENT(IN) :: LDA, N, JOB
      INTEGER, DIMENSION(N), INTENT(IN) :: IPVT
      REAL(WP), DIMENSION(LDA,N), INTENT(IN) :: A
      REAL(WP), DIMENSION(N), INTENT(INOUT) :: B
C-----------------------------------------------
C   L O C A L   V A R I A B L E S
C-----------------------------------------------
      INTEGER :: K, L, NM1
      REAL(WP) :: ELEMENT
C-----------------------------------------------
c
c     sgesl solves the REAL system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by sgeco or sgefa.
c
c     on ENTRY
c
c        a       REAL(lda, n)
c                the output from sgeco or sgefa.
c
c        lda     INTEGER
c                the leading DIMENSION of the array  a .
c
c        n       INTEGER
c                the order of the matrix  a .
c
c        ipvt    INTEGER(n)
c                the pivot vector from sgeco or sgefa.
c
c        b       REAL(n)
c                the right hand side vector.
c
c        job     INTEGER
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b  WHERE
c                            trans(a)  is the TRANSPOSE.
c
c     on RETURN
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if sgeco has set rcond.gt.0.0
c        or sgefa has set info.eq.0 .
c
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c
      nm1 = n - 1
      IF (job .eq. 0) THEN
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         IF (nm1 .ge. 1) THEN
            DO k = 1, nm1
               l = ipvt(k)
               element = b(l)
               IF (l .ne. k) THEN
                  b(l) = b(k)
                  b(k) = element
               END IF
!              CALL bla_axpy(n-k,element,a(k+1,k),1,b(k+1),1)
               b(k+1:n) = b(k+1:n) + element*a(k+1:n,k)
            END DO
         END IF
c
c        now solve  u*x = y
c
         DO k = n, 1, -1
            b(k) = b(k)/a(k,k)
            element = -b(k)
!           CALL bla_axpy(k-1,element,a(1,k),1,b,1)
            b(1:k-1) = b(1:k-1) + element*a(1:k-1,k)
         END DO
      ELSE
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         DO k = 1, n
            element = SUM(a(:k-1,k)*b(:k-1))
            b(k) = (b(k)-element)/a(k,k)
         END DO
c
c        now solve trans(l)*x = y
c
         IF (nm1 .ge. 1) THEN
            DO k = nm1, 1, -1
               b(k) = b(k) + SUM(a(k+1:n,k)*b(k+1:n))
               l = ipvt(k)
               IF (l .ne. k) THEN
                  element = b(l)
                  b(l) = b(k)
                  b(k) = element
               END IF
            END DO
         END IF
      END IF

      END SUBROUTINE dgesl
