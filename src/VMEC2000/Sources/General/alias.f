      SUBROUTINE alias(gcons, gcona, ztemp, gcs, gsc, gcc)
      USE vmec_main
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: p5 = 0.5_dp
      REAL(rprec), DIMENSION(ns*nzeta,ntheta3), INTENT(out) ::
     1   gcons, gcona
      REAL(rprec), DIMENSION(ns*nzeta,ntheta3), INTENT(in)  :: ztemp
      REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1) :: gcs, gsc, gcc
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: m, i, ir, jk, jka, n, k, js, l
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: work
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: gss
C-----------------------------------------------
      ALLOCATE (work(ns*nzeta,4), gss(ns,0:ntor,0:mpol1))

      gcons = 0
      gcona = 0

      gcs(:,:ntor,:mpol1) = 0;  gsc(:,:ntor,:mpol1) = 0
      gcc(:,:ntor,:mpol1) = 0;  gss(:,:ntor,:mpol1) = 0

      DO m = 1, mpol1 - 1
         work = 0
         DO i = 1, ntheta2
            DO jk = 1, ns*nzeta
               work(jk,1) = work(jk,1) + ztemp(jk,i)*cosmui(i,m)
               work(jk,2) = work(jk,2) + ztemp(jk,i)*sinmui(i,m)
            END DO
            IF (.not.lasym) CYCLE
            ir = ntheta1 + 2 - i
            IF (i .eq. 1) ir = 1
            DO jk = 1, ns*nzeta
               jka = ireflect(jk)
               work(jk,3) = work(jk,3) + ztemp(jka,ir)*cosmui(i,m)
               work(jk,4) = work(jk,4) + ztemp(jka,ir)*sinmui(i,m)
            END DO
         END DO

         DO n = 0, ntor
            DO k = 1, nzeta
               l = ns*(k-1)
               IF (.not.lasym) THEN
               DO js = 2,ns
                  gcs(js,n,m) = gcs(js,n,m) + tcon(js)*work(js+l,1)*
     1               sinnv(k,n)
                  gsc(js,n,m) = gsc(js,n,m) + tcon(js)*work(js+l,2)*
     1               cosnv(k,n)
               END DO
               ELSE
               DO js = 2,ns
                  gcs(js,n,m) = gcs(js,n,m) + p5*tcon(js)*sinnv(k,n)*
     1               (work(js+l,1)-work(js+l,3))
                  gsc(js,n,m) = gsc(js,n,m) + p5*tcon(js)*cosnv(k,n)*
     1               (work(js+l,2)-work(js+l,4))
                  gss(js,n,m) = gss(js,n,m) + p5*tcon(js)*sinnv(k,n)*
     1               (work(js+l,2)+work(js+l,4))
                  gcc(js,n,m) = gcc(js,n,m) + p5*tcon(js)*cosnv(k,n)*
     1               (work(js+l,1)+work(js+l,3))
               END DO
               END IF
            END DO
         END DO
!
!        INVERSE FOURIER TRANSFORM DE-ALIASED GCON
!
         work = 0

         DO n = 0, ntor
            DO k = 1, nzeta
               l = ns*(k-1)
               DO js = 2, ns
                  work(js+l,3) = work(js+l,3) + gcs(js,n,m)*sinnv(k,n)
                  work(js+l,4) = work(js+l,4) + gsc(js,n,m)*cosnv(k,n)
               END DO
               IF (.not.lasym) CYCLE
               DO js = 2, ns
                  work(js+l,1) = work(js+l,1) + gcc(js,n,m)*cosnv(k,n)
                  work(js+l,2) = work(js+l,2) + gss(js,n,m)*sinnv(k,n)
               END DO
            END DO
         END DO
         DO i = 1, ntheta2
            DO jk = 1, ns*nzeta
               gcons(jk,i) = gcons(jk,i) + (work(jk,3)*cosmu(i,m)
     1                     + work(jk,4)*sinmu(i,m))*faccon(m)
            END DO
            IF (.not.lasym) CYCLE
            DO jk = 1, ns*nzeta
               gcona(jk,i) = gcona(jk,i) + (work(jk,1)*cosmu(i,m)
     1                     + work(jk,2)*sinmu(i,m))*faccon(m)
            END DO
         END DO
      END DO

      IF (lasym) THEN
!
!     EXTEND GCON INTO THETA = PI,2*PI DOMAIN
!
         DO i = 1 + ntheta2, ntheta1
            ir = ntheta1 + 2 - i
            DO jk = 1, ns*nzeta
               jka = ireflect(jk)
               gcons(jk,i) = -gcons(jka,ir) + gcona(jka,ir)
            END DO
         END DO
!
!     ADD SYMMETRIC, ANTI-SYMMETRIC PIECES IN THETA = 0,PI DOMAIN
!
         gcons(:,:ntheta2) = gcons(:,:ntheta2) + gcona(:,:ntheta2)

      END IF

      DEALLOCATE (work, gss)

      END SUBROUTINE alias
