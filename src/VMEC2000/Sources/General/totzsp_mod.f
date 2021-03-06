      MODULE totzsp_mod

      CONTAINS

      SUBROUTINE totzsps(rzl_array, r11, ru1, rv1, z11, zu1, zv1,
     1                   lu1, lv1, rcn1, zcn1)
      USE vmec_main
      USE vmec_params, ONLY: jmin1, jlam, ntmax, rcc, rss, zsc, zcs
      USE precon2d, ONLY: ictrl_prec2d
      USE xstuff, ONLY: xc
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,3*ntmax),
     1   TARGET, INTENT(inout) :: rzl_array
      REAL(rprec), DIMENSION(ns*nzeta*ntheta3,0:1),
     1   INTENT(out) :: r11, ru1,
     1   rv1, z11, zu1,  zv1, lu1, lv1, rcn1, zcn1
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: m0 = 0, m1 = 1, n0 = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: n, m, mparity, k, i, j1, l, j1l, nsl
      INTEGER :: ioff, joff, mj, ni, nsz
      REAL(rprec), DIMENSION(:,:,:), POINTER ::
     1           rmncc, rmnss, zmncs, zmnsc, lmncs, lmnsc
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: work1
      REAL(rprec) :: cosmux, sinmux
C-----------------------------------------------
!
!     WORK1 Array of inverse transforms in toroidal angle (zeta), for all radial positions
!     NOTE: ORDERING OF LAST INDEX IS DIFFERENT HERE THAN IN PREVIOUS VMEC2000 VERSIONS
!
!     CONVERT FROM INTERNAL XC REPRESENTATION FOR m=1 MODES, R+(stored at rss) = .5(rss + zcs),
!     R-(stored at zcs) = .5(rss - zcs), TO EXTERNAL ("PHYSICAL") rss, zcs FORMS. NEED THIS EVEN 
!     WHEN COMPUTING HESSIAN FOR FREE BOUNDARY (rmnss, zmncs at JS=NS needed in vacuum call)
!
!     WHEN COMPUTING PRECONDITIONER, USE FASTER HESSIAN VERSION (totzsps_hess) INSTEAD.

      rmncc => rzl_array(:,:,:,rcc)               !!COS(mu) COS(nv)
      zmnsc => rzl_array(:,:,:,zsc+ntmax)         !!SIN(mu) COS(nv)
      lmnsc => rzl_array(:,:,:,zsc+2*ntmax)       !!SIN(mu) COS(nv)
      IF (lthreed) THEN
         rmnss => rzl_array(:,:,:,rss)               !!SIN(mu) SIN(nv)
         zmncs => rzl_array(:,:,:,zcs+ntmax)         !!COS(mu) SIN(nv)
         lmncs => rzl_array(:,:,:,zcs+2*ntmax)       !!COS(mu) SIN(nv)
         CALL convert_sym (rmnss, zmncs)
      END IF

!
!     ORIGIN EXTRAPOLATION (JS=1) FOR M=1 MODES
!     NOTE: PREVIOUS VERSIONS OF VMEC USED TWO-POINT EXTRAPOLATION 
!           FOR R,Z. HOWEVER,THIS CAN NOT BE USED TO COMPUTE THE 
!           TRI-DIAG 2D PRECONDITIONER
!
      rzl_array(1,:,m1,:)  = rzl_array(2,:,m1,:)

      ioff = LBOUND(rmncc,2)
      joff = LBOUND(rmncc,3)

!
!     ORIGIN EXTRAPOLATION OF M=0 MODES FOR LAMBDA 
!
      IF (lthreed .and. jlam(m0).gt.1) 
     1   lmncs(1,:,m0+joff) = lmncs(2,:,m0+joff)
 
      IF (ictrl_prec2d .eq. 3) RETURN
     
      nsz = ns*nzeta
      ALLOCATE (work1(nsz,12), stat=i)
      IF (i .ne. 0) STOP 'Allocation error in VMEC2000 totzsps'

      r11 = 0;  ru1 = 0;  rv1 = 0;  rcn1 = 0
      z11 = 0;  zu1 = 0;  zv1 = 0;  zcn1 = 0
      lu1 = 0;  lv1 = 0

!
!     COMPUTE R, Z, AND LAMBDA IN REAL SPACE
!     NOTE: LU = d(Lam)/du, LV = -d(Lam)/dv
!

      DO m = 0, mpol1
         mparity = MOD(m,2)
         mj = m+joff
         work1 = 0
         j1 = jmin1(m)
!
!        INVERSE TRANSFORM IN N-ZETA, FOR FIXED M
!
         DO n = 0, ntor
            ni = n+ioff
            DO k = 1, nzeta
               l = ns*(k-1)
               j1l = j1+l;  nsl = ns+l
               work1(j1l:nsl,1) = work1(j1l:nsl,1) 
     1                          + rmncc(j1:ns,ni,mj)*cosnv(k,n)
               work1(j1l:nsl,6) = work1(j1l:nsl,6)
     1                          + zmnsc(j1:ns,ni,mj)*cosnv(k,n)
               work1(j1l:nsl,10) = work1(j1l:nsl,10) 
     1                          + lmnsc(j1:ns,ni,mj)*cosnv(k,n)

               IF (.not.lthreed) CYCLE
               
               work1(j1l:nsl,4) = work1(j1l:nsl,4) 
     1                          + rmnss(j1:ns,ni,mj)*cosnvn(k,n)
               work1(j1l:nsl,7) = work1(j1l:nsl,7) 
     1                          + zmncs(j1:ns,ni,mj)*cosnvn(k,n)
               work1(j1l:nsl,11) = work1(j1l:nsl,11)
     1                          + lmncs(j1:ns,ni,mj)*cosnvn(k,n)

               work1(j1l:nsl,2) = work1(j1l:nsl,2) 
     1                          + rmnss(j1:ns,ni,mj)*sinnv(k,n)
               work1(j1l:nsl,5) = work1(j1l:nsl,5) 
     1                          + zmncs(j1:ns,ni,mj)*sinnv(k,n)
               work1(j1l:nsl,9) = work1(j1l:nsl,9) 
     1                          + lmncs(j1:ns,ni,mj)*sinnv(k,n)

               work1(j1l:nsl,3) = work1(j1l:nsl,3) 
     1                          + rmncc(j1:ns,ni,mj)*sinnvn(k,n)
               work1(j1l:nsl,8) = work1(j1l:nsl,8) 
     1                          + zmnsc(j1:ns,ni,mj)*sinnvn(k,n)
               work1(j1l:nsl,12) = work1(j1l:nsl,12)
     1                          + lmnsc(j1:ns,ni,mj)*sinnvn(k,n)
            END DO
         END DO
!
!        INVERSE TRANSFORM IN M-THETA, FOR ALL RADIAL, ZETA VALUES
!
         l = 0
         DO i = 1, ntheta2
            j1l = l+1;  nsl = nsz+l
            l = l + nsz
            cosmux = xmpq(m,1)*cosmu(i,m)
            sinmux = xmpq(m,1)*sinmu(i,m)
           
            r11(j1l:nsl,mparity)  = r11(j1l:nsl,mparity)  
     1                            + work1(1:nsz,1)*cosmu(i,m)
            ru1(j1l:nsl,mparity)  = ru1(j1l:nsl,mparity)  
     1                            + work1(1:nsz,1)*sinmum(i,m)
            rcn1(j1l:nsl,mparity) = rcn1(j1l:nsl,mparity) 
     1                            + work1(1:nsz,1)*cosmux
            z11(j1l:nsl,mparity)  = z11(j1l:nsl,mparity)  
     1                            + work1(1:nsz,6)*sinmu(i,m)

            zu1(j1l:nsl,mparity)  = zu1(j1l:nsl,mparity)  
     1                            + work1(1:nsz,6)*cosmum(i,m)
            zcn1(j1l:nsl,mparity) = zcn1(j1l:nsl,mparity) 
     1                            + work1(1:nsz,6)*sinmux
             
            lu1(j1l:nsl,mparity)  = lu1(j1l:nsl,mparity)  
     1                            + work1(1:nsz,10)*cosmum(i,m)

            IF (.not.lthreed) CYCLE

            r11(j1l:nsl,mparity)  = r11(j1l:nsl,mparity)  
     1                            + work1(1:nsz,2)*sinmu(i,m)
            ru1(j1l:nsl,mparity)  = ru1(j1l:nsl,mparity)  
     1                            + work1(1:nsz,2)*cosmum(i,m)
            rcn1(j1l:nsl,mparity) = rcn1(j1l:nsl,mparity) 
     1                            + work1(1:nsz,2)*sinmux
            rv1(j1l:nsl,mparity)  = rv1(j1l:nsl,mparity)  
     1                            + work1(1:nsz,3)*cosmu(i,m) 
     1                            + work1(1:nsz,4)*sinmu(i,m)
            z11(j1l:nsl,mparity)  = z11(j1l:nsl,mparity)  
     1                            + work1(1:nsz,5)*cosmu(i,m)

            zu1(j1l:nsl,mparity)  = zu1(j1l:nsl,mparity)  
     1                            + work1(1:nsz,5)*sinmum(i,m)
            zcn1(j1l:nsl,mparity) = zcn1(j1l:nsl,mparity) 
     1                            + work1(1:nsz,5)*cosmux
            zv1(j1l:nsl,mparity)  = zv1(j1l:nsl,mparity)  
     1                            + work1(1:nsz,7)*cosmu(i,m) 
     1                            + work1(1:nsz,8)*sinmu(i,m)

            lu1(j1l:nsl,mparity)  = lu1(j1l:nsl,mparity)  
     1                            + work1(1:nsz,9)*sinmum(i,m)
            lv1(j1l:nsl,mparity)  = lv1(j1l:nsl,mparity)  
     1                            - (work1(1:nsz,11)*cosmu(i,m) 
     1                            +  work1(1:nsz,12)*sinmu(i,m))
         END DO

      END DO

      DEALLOCATE (work1)

      z01(1:ns) = zmnsc(1:ns,n0+ioff,m1+joff)
      r01(1:ns) = rmncc(1:ns,n0+ioff,m1+joff)
      IF (r01(1) .eq. zero) STOP 'r01(0) = 0 in totzsps'
      dkappa = z01(1)/r01(1)

      END SUBROUTINE totzsps

 
      SUBROUTINE totzsps_hess(rzl_array, r11, ru1, rv1, z11, zu1, 
     1                        zv1, lu1, lv1, rcn1, zcn1)
      USE vmec_main
      USE vmec_params, ONLY: jmin1, jlam, ntmax, rcc, rss, zsc, zcs
      USE precon2d
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,3*ntmax),
     1   TARGET, INTENT(inout) :: rzl_array
      REAL(rprec), DIMENSION(ns*nzeta,ntheta3,0:1),
     1   INTENT(inout) :: r11, ru1,
     1   rv1, z11, zu1,  zv1, lu1, lv1, rcn1, zcn1
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: m0 = 0, m1 = 1, n0 = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: n, m, mparity, k, i, j1, l, j1l, nsl
      INTEGER :: ioff, joff, mj, ni
      REAL(rprec), DIMENSION(:,:,:), POINTER ::
     1           rmncc, rmnss, zmncs, zmnsc, lmncs, lmnsc
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: work1
      REAL(rprec) :: cosmux, sinmux
      LOGICAL :: logl, logr, logz
C-----------------------------------------------
!
!     SAME AS totzsps, BUT ONLY COMPUTES perturbation (rzl_array is perturbation) TO EXISTING R,Z FOR
!     A PARTICULAR (m,n,ntype) = (m_2d,n_2d,ntype_2d) VALUE
!
!     
!     STORE (FOR ictrl_prec2d=2) OR RESTORE (ictrl_prec2d=3) VARIABLES NOT RECOMPUTED 
!     IN TOTZSP DURING HESSIAN JOGS
!
      IF (ictrl_prec2d .eq. 2) THEN
         lus_save = lu1;  lvs_save = lv1;
         rus_save = ru1;  rvs_save = rv1; r1s_save = r11
         rcons_save = rcn1
         zus_save = zu1;  zvs_save = zv1; z1s_save = z11
         zcons_save = zcn1
         RETURN
      ELSE
         lu1 = lus_save;  lv1 = lvs_save
         ru1 = rus_save;  rv1 = rvs_save; r11 = r1s_save
         rcn1 = rcons_save
         zu1 = zus_save;  zv1 = zvs_save; z11 = z1s_save
         zcn1 = zcons_save
      END IF

      logr = ntype_2d .le. ntmax
      logz = (ntype_2d.gt.ntmax) .and. (ntype_2d.le.2*ntmax)
      logl = ntype_2d .gt. 2*ntmax

      rmncc => rzl_array(:,:,:,rcc)               !!COS(mu) COS(nv)
      zmnsc => rzl_array(:,:,:,zsc+ntmax)         !!SIN(mu) COS(nv)
      lmnsc => rzl_array(:,:,:,zsc+2*ntmax)       !!SIN(mu) COS(nv)
      IF (lthreed) THEN
         rmnss => rzl_array(:,:,:,rss)            !!SIN(mu) SIN(nv)
         zmncs => rzl_array(:,:,:,zcs+ntmax)      !!COS(mu) SIN(nv)
         lmncs => rzl_array(:,:,:,zcs+2*ntmax)    !!COS(mu) SIN(nv)
      END IF

      ioff = LBOUND(rmncc,2)
      joff = LBOUND(rmncc,3)

      ALLOCATE (work1(ns*nzeta,12), stat=i)
      IF (i .ne. 0) STOP 'Allocation error in VMEC2000 totzsps'

!
!     ENFORCE CONSTRAINT ON m=1 MODES FOR 3D (CONSISTENT WITH gcz(zcs) = 0 IN RESIDUE)
!     NOTE: Since r,z variations are coupled, must vary BOTH whenever ONE is varied
!

      IF (lthreed .and. m_2d.eq.1 .and. n_2d.ne.0 .and. 
     1   (ntype_2d.eq.rss .or. ntype_2d.eq.(zcs+ntmax))) THEN
!         zmncs(:,:,m1+joff) = 0
!         IF (logz) logz = .false.
!         IF (logr .and. lconm1) logz = .true.
         logr = .true.;  logz = .true.
         CALL convert_sym (rmnss, zmncs)
      END IF

!
!     EXTRAPOLATION AT JS=1 FOR M=1 MODES
!
      rzl_array(1,:,m1,:) = rzl_array(2,:,m1,:)   

!
!     EXTRAPOLATION AT JS=1 OF M=0 MODES FOR LAMBDA
!
      IF (lthreed .and. jlam(m0).gt.1) 
     1   lmncs(1,:,m0+joff) = lmncs(2,:,m0+joff)

!
!     COMPUTE R, Z, AND LAMBDA IN REAL SPACE
!     NOTE: LU = d(Lam)/du, LV = -d(Lam)/dv
!
         m = m_2d
         mparity = MOD(m,2)
         mj = m+joff
         work1 = 0
         j1 = jmin1(m)
!
!        INVERSE TRANSFORM IN N-ZETA, FOR FIXED M
!
         n = n_2d
            ni = n+ioff
            DO k = 1, nzeta
               l = ns*(k-1)
               j1l = j1+l;  nsl = ns+l
               IF (logr)
     1         work1(j1l:nsl,1) = work1(j1l:nsl,1) 
     2                          + rmncc(j1:ns,ni,mj)*cosnv(k,n)
               IF (logz)
     1         work1(j1l:nsl,6) = work1(j1l:nsl,6)
     2                          + zmnsc(j1:ns,ni,mj)*cosnv(k,n)
               IF (logl)
     1         work1(j1l:nsl,10) = work1(j1l:nsl,10) 
     2                          + lmnsc(j1:ns,ni,mj)*cosnv(k,n)

               IF (.not.lthreed) CYCLE
     
               IF (logr) THEN          
               work1(j1l:nsl,4) = work1(j1l:nsl,4) 
     1                          + rmnss(j1:ns,ni,mj)*cosnvn(k,n)
               work1(j1l:nsl,2) = work1(j1l:nsl,2) 
     1                          + rmnss(j1:ns,ni,mj)*sinnv(k,n)
               work1(j1l:nsl,3) = work1(j1l:nsl,3) 
     1                          + rmncc(j1:ns,ni,mj)*sinnvn(k,n)
               END IF

               IF (logz) THEN
               work1(j1l:nsl,7) = work1(j1l:nsl,7) 
     1                          + zmncs(j1:ns,ni,mj)*cosnvn(k,n)
               work1(j1l:nsl,5) = work1(j1l:nsl,5) 
     1                          + zmncs(j1:ns,ni,mj)*sinnv(k,n)
               work1(j1l:nsl,8) = work1(j1l:nsl,8) 
     1                          + zmnsc(j1:ns,ni,mj)*sinnvn(k,n)
               END IF

               IF (logl) THEN
               work1(j1l:nsl,11) = work1(j1l:nsl,11)
     1                          + lmncs(j1:ns,ni,mj)*cosnvn(k,n)
               work1(j1l:nsl,9) = work1(j1l:nsl,9) 
     1                          + lmncs(j1:ns,ni,mj)*sinnv(k,n)
               work1(j1l:nsl,12) = work1(j1l:nsl,12)
     1                          + lmnsc(j1:ns,ni,mj)*sinnvn(k,n)
               END IF

            END DO
!
!        INVERSE TRANSFORM IN M-THETA, FOR ALL RADIAL, ZETA VALUES
!
         DO i = 1, ntheta2
            cosmux = xmpq(m,1)*cosmu(i,m)
            sinmux = xmpq(m,1)*sinmu(i,m)
            IF (logr) THEN
            r11(:,i,mparity) = r11(:,i,mparity) + work1(:,1)*
     1            cosmu(i,m)
            ru1(:,i,mparity) = ru1(:,i,mparity) + work1(:,1)*
     1            sinmum(i,m)
            rcn1(:,i,mparity) = rcn1(:,i,mparity) + work1(:,1)*
     1            cosmux
            END IF

            IF (logz) THEN
            z11(:,i,mparity) = z11(:,i,mparity) + work1(:,6)*
     1            sinmu(i,m)
            zu1(:,i,mparity) = zu1(:,i,mparity) + work1(:,6)*
     1            cosmum(i,m)
            zcn1(:,i,mparity) = zcn1(:,i,mparity) + work1(:,6)*
     1            sinmux
            END IF
 
            IF (logl)
     1      lu1(:,i,mparity) = lu1(:,i,mparity) + work1(:,10)*
     2            cosmum(i,m)

            IF (.not.lthreed) CYCLE

            IF (logr) THEN
            r11(:,i,mparity) = r11(:,i,mparity) + work1(:,2)*
     1            sinmu(i,m)
            ru1(:,i,mparity) = ru1(:,i,mparity) + work1(:,2)*
     1            cosmum(i,m)
            rcn1(:,i,mparity) = rcn1(:,i,mparity) + work1(:,2)*
     1            sinmux
            rv1(:,i,mparity) = rv1(:,i,mparity) + work1(:,3)*
     1            cosmu(i,m) + work1(:,4)*sinmu(i,m)
            END IF

            IF (logz) THEN
            z11(:,i,mparity) = z11(:,i,mparity) + work1(:,5)*
     1            cosmu(i,m)
            zu1(:,i,mparity) = zu1(:,i,mparity) + work1(:,5)*
     1            sinmum(i,m)
            zcn1(:,i,mparity) = zcn1(:,i,mparity) + work1(:,5)*
     1            cosmux
            zv1(:,i,mparity) = zv1(:,i,mparity) + work1(:,7)*
     1            cosmu(i,m) + work1(:,8)*sinmu(i,m)
            END IF

            IF (logl) THEN
            lu1(:,i,mparity) = lu1(:,i,mparity) + work1(:,9)*
     1            sinmum(i,m)
            lv1(:,i,mparity) = lv1(:,i,mparity) - (work1(:,11)*
     1            cosmu(i,m) + work1(:,12)*sinmu(i,m))
            END IF

         END DO

      DEALLOCATE (work1)

      END SUBROUTINE totzsps_hess

      SUBROUTINE convert_sym(rmnss, zmncs)
      USE vmec_main, p5 => cp5
      USE precon2d, ONLY: ictrl_prec2d
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1), INTENT(inout) :: 
     1                                           rmnss, zmncs
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), DIMENSION(ns,0:ntor) :: temp
C-----------------------------------------------
!
!     CONVERT FROM INTERNAL REPRESENTATION TO "PHYSICAL" RMNSS, ZMNCS FOURIER FORM
!
      IF (.not.lconm1) RETURN

      temp(:,:) = rmnss(:,:,1)
      rmnss(:,:,1) = temp(:,:) + zmncs(:,:,1)
      zmncs(:,:,1) = temp(:,:) - zmncs(:,:,1)

      END SUBROUTINE convert_sym


      SUBROUTINE totzspa(rzl_array, r11, ru1, rv1, z11, zu1, zv1, lu1,
     1   lv1, rcn1, zcn1)
      USE vmec_main
      USE vmec_params, ONLY: jmin1, jlam, ntmax, rcs, rsc, zcc, zss
      USE precon2d, ONLY: ictrl_prec2d
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,3*ntmax),
     1   TARGET, INTENT(inout) :: rzl_array
      REAL(rprec), DIMENSION(ns*nzeta,ntheta3,0:1), INTENT(out) ::
     1   r11, ru1, rv1, z11, zu1, zv1, lu1, lv1, rcn1, zcn1
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: m0 = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: m, n, mparity, k, i, l, j1, j1l, nsl
      INTEGER :: ioff, joff, mj, ni
      REAL(rprec), DIMENSION(:,:,:), POINTER ::
     1           rmncs, rmnsc, zmncc, zmnss, lmncc, lmnss
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: work1
      REAL(rprec) :: cosmux, sinmux
C-----------------------------------------------
!     WHEN COMPUTING PRECONDITIONER, USE FASTER HESSIAN VERSION (totzsps_hess) INSTEAD.

      rmnsc => rzl_array(:,:,:,rsc)               !!SIN(mu) COS(nv)
      zmncc => rzl_array(:,:,:,zcc+ntmax)         !!COS(mu) COS(nv)
      lmncc => rzl_array(:,:,:,zcc+2*ntmax)       !!COS(mu) COS(nv)
      IF (lthreed) THEN
         rmncs => rzl_array(:,:,:,rcs)               !!COS(mu) SIN(nv)
         zmnss => rzl_array(:,:,:,zss+ntmax)         !!SIN(mu) SIN(nv)
         lmnss => rzl_array(:,:,:,zss+2*ntmax)       !!SIN(mu) SIN(nv)
      END IF

!
!     CONVERT FROM INTERNAL XC REPRESENTATION FOR m=1 MODES, R+(at rsc) = .5(rsc + zcc),
!     R-(at zcc) = .5(rsc - zcc), TO REQUIRED rsc, zcc FORMS
!
      CALL convert_asym (rmnsc, zmncc)

      IF (ictrl_prec2d .eq. 3) RETURN

      ioff = LBOUND(rmnsc,2)
      joff = LBOUND(rmnsc,3)

      z00b = zmncc(ns,ioff,joff)

      ALLOCATE (work1(ns*nzeta,12), stat=i)
      IF (i .ne. 0) STOP 'Allocation error in VMEC totzspa'

!
!     INITIALIZATION BLOCK
!
      r11 = 0;  ru1 = 0;  rv1 = 0;  z11 = 0;  zu1 = 0
      zv1 = 0;  lu1 = 0;  lv1 = 0;  rcn1 = 0; zcn1 = 0

      IF (jlam(m0) .gt. 1) lmncc(1,:,m0+joff) = lmncc(2,:,m0+joff)

      DO m = 0, mpol1
         mparity = MOD(m,2)
         mj = m+joff
         work1 = 0
         j1 = jmin1(m)
         DO n = 0, ntor
            ni = n+ioff
            DO k = 1, nzeta
               l = ns*(k-1)
               j1l = j1+l;  nsl = ns+l
               work1(j1l:nsl,1) = work1(j1l:nsl,1)
     1                          + rmnsc(j1:ns,ni,mj)*cosnv(k,n)
               work1(j1l:nsl,6) = work1(j1l:nsl,6) 
     1                          + zmncc(j1:ns,ni,mj)*cosnv(k,n)
               work1(j1l:nsl,10) = work1(j1l:nsl,10)
     1                          + lmncc(j1:ns,ni,mj)*cosnv(k,n)

               IF (.not.lthreed) CYCLE

               work1(j1l:nsl,2) = work1(j1l:nsl,2)
     1                          + rmncs(j1:ns,ni,mj)*sinnv(k,n)
               work1(j1l:nsl,3) = work1(j1l:nsl,3) 
     1                          + rmnsc(j1:ns,ni,mj)*sinnvn(k,n)
               work1(j1l:nsl,4) = work1(j1l:nsl,4) 
     1                          + rmncs(j1:ns,ni,mj)*cosnvn(k,n)
               work1(j1l:nsl,5) = work1(j1l:nsl,5) 
     1                          + zmnss(j1:ns,ni,mj)*sinnv(k,n)
               work1(j1l:nsl,7) = work1(j1l:nsl,7)
     1                          + zmnss(j1:ns,ni,mj)*cosnvn(k,n)
               work1(j1l:nsl,8) = work1(j1l:nsl,8) 
     1                          + zmncc(j1:ns,ni,mj)*sinnvn(k,n)
               work1(j1l:nsl,9) = work1(j1l:nsl,9) 
     1                          + lmnss(j1:ns,ni,mj)*sinnv(k,n)
               work1(j1l:nsl,11) = work1(j1l:nsl,11)
     1                          + lmnss(j1:ns,ni,mj)*cosnvn(k,n)
               work1(j1l:nsl,12) = work1(j1l:nsl,12) 
     1                          + lmncc(j1:ns,ni,mj)*sinnvn(k,n)
            END DO
         END DO

!
!        INVERSE TRANSFORM IN M-THETA
!
         DO i = 1, ntheta2
            cosmux = xmpq(m,1)*cosmu(i,m)
            sinmux = xmpq(m,1)*sinmu(i,m)
            r11(:,i,mparity) = r11(:,i,mparity) + work1(:,1)*
     1            sinmu(i,m)
            ru1(:,i,mparity) = ru1(:,i,mparity) + work1(:,1)*
     1            cosmum(i,m)
            z11(:,i,mparity) = z11(:,i,mparity) + work1(:,6)*
     1            cosmu(i,m)
            zu1(:,i,mparity) = zu1(:,i,mparity) + work1(:,6)*
     1            sinmum(i,m)
            lu1(:,i,mparity) = lu1(:,i,mparity) + work1(:,10)*
     1            sinmum(i,m)
            rcn1(:,i,mparity) = rcn1(:,i,mparity) + work1(:,1)*
     1            sinmux
            zcn1(:,i,mparity) = zcn1(:,i,mparity) + work1(:,6)*
     1            cosmux

            IF (.not.lthreed) CYCLE
               
            r11(:,i,mparity) = r11(:,i,mparity) + work1(:,2)*
     1               cosmu(i,m)
            ru1(:,i,mparity) = ru1(:,i,mparity) + work1(:,2)*
     1               sinmum(i,m)
            z11(:,i,mparity) = z11(:,i,mparity) + work1(:,5)*
     1               sinmu(i,m)
            zu1(:,i,mparity) = zu1(:,i,mparity) + work1(:,5)*
     1               cosmum(i,m)
            lu1(:,i,mparity) = lu1(:,i,mparity) + work1(:,9)*
     1               cosmum(i,m)
            rcn1(:,i,mparity) = rcn1(:,i,mparity) + work1(:,2)*
     1               cosmux
            zcn1(:,i,mparity) = zcn1(:,i,mparity) + work1(:,5)*
     1               sinmux
            rv1(:,i,mparity) = rv1(:,i,mparity) + work1(:,3)*
     1               sinmu(i,m) + work1(:,4)*cosmu(i,m)
            zv1(:,i,mparity) = zv1(:,i,mparity) + work1(:,7)*
     1               sinmu(i,m) + work1(:,8)*cosmu(i,m)
            lv1(:,i,mparity) = lv1(:,i,mparity) - (work1(:,11)*
     1               sinmu(i,m)+work1(:,12)*cosmu(i,m))
         END DO
      END DO

      DEALLOCATE (work1)

      END SUBROUTINE totzspa

      SUBROUTINE totzspa_hess(rzl_array, r11, ru1, rv1, z11, zu1, 
     1                        zv1, lu1, lv1, rcn1, zcn1)
      USE vmec_main
      USE vmec_params, ONLY: jmin1, jlam, ntmax, rcs, rsc, zcc, zss
      USE precon2d
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,3*ntmax),
     1   TARGET, INTENT(inout) :: rzl_array
      REAL(rprec), DIMENSION(ns*nzeta,ntheta3,0:1),
     1   INTENT(inout) :: r11, ru1,
     1   rv1, z11, zu1,  zv1, lu1, lv1, rcn1, zcn1
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      INTEGER, PARAMETER :: m0 = 0, m1 = 1, n0 = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: n, m, mparity, k, i, j1, l, j1l, nsl
      INTEGER :: ioff, joff, mj, ni
      REAL(rprec), DIMENSION(:,:,:), POINTER ::
     1           rmncs, rmnsc, zmncc, zmnss, lmncc, lmnss
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: work1
      REAL(rprec) :: cosmux, sinmux
      LOGICAL :: logl, logr, logz
C-----------------------------------------------
!
!     SAME AS totzspa, BUT ONLY COMPUTES perturbation (rzl_array is perturbation) TO EXISTING R,Z, L FOR
!     A PARTICULAR (m,n,ntype) = (m_2d,n_2d,ntype_2d) VALUE
!
!     
!     STORE (FOR ictrl_prec2d=2) OR RESTORE (ictrl_prec2d=3) VARIABLES NOT RECOMPUTED 
!     IN TOTZSP DURING HESSIAN JOGS
!
      IF (ictrl_prec2d .eq. 2) THEN
         lua_save = lu1;  lva_save = lv1;
         rua_save = ru1;  rva_save = rv1; r1a_save = r11 
         rcona_save = rcn1
         zua_save = zu1;  zva_save = zv1; z1a_save = z11 
         zcona_save = zcn1
         RETURN
      ELSE
         lu1 = lua_save;  lv1 = lva_save
         ru1 = rua_save;  rv1 = rva_save; r11 = r1a_save 
         rcn1 = rcona_save
         zu1 = zua_save;  zv1 = zva_save; z11 = z1a_save
         zcn1 = zcona_save
      END IF

      logr = ntype_2d .le. ntmax
      logz = (ntype_2d.gt.ntmax) .and. (ntype_2d.le.2*ntmax)
      logl = ntype_2d .gt. 2*ntmax

      rmnsc => rzl_array(:,:,:,rsc)               !!SIN(mu) COS(nv)
      zmncc => rzl_array(:,:,:,zcc+ntmax)         !!COS(mu) COS(nv)
      lmncc => rzl_array(:,:,:,zcc+2*ntmax)       !!COS(mu) COS(nv)
      IF (lthreed) THEN
         rmncs => rzl_array(:,:,:,rcs)               !!COS(mu) SIN(nv)
         zmnss => rzl_array(:,:,:,zss+ntmax)         !!SIN(mu) SIN(nv)
         lmnss => rzl_array(:,:,:,zss+2*ntmax)       !!SIN(mu) SIN(nv)
      END IF

      ioff = LBOUND(rmnsc,2)
      joff = LBOUND(rmnsc,3)

!
!     ENFORCE CONSTRAINT ON m=1 MODES FOR 3D (CONSISTENT WITH gcz(zcc) = 0 IN RESIDUE)
!     NOTE: Since r,z variations are coupled, must vary BOTH whenever ONE is varied
!
      IF (m_2d.eq.1 .and. 
     1   (ntype_2d.eq.rsc .or. ntype_2d.eq.(ntmax+zcc))) THEN
!         zmncc(:,:,m1+joff) = 0
!         IF (logz) logz = .false.
!         IF (logr .and. lconm1) logz = .true.
         logr = .true.;  logz = .true.
         CALL convert_asym (rmnsc, zmncc)
      END IF


      ALLOCATE (work1(ns*nzeta,12), stat=i)
      IF (i .ne. 0) STOP 'Allocation error in VMEC2000 totzsps'


!     EXTRAPOLATION OF M=0 MODES FOR LAMBDA
!
      IF (jlam(m0) .gt. 1) lmncc(1,:,m0+joff) = lmncc(2,:,m0+joff)

!
!     COMPUTE R, Z, AND LAMBDA IN REAL SPACE
!     NOTE: LU = d(Lam)/du, LV = -d(Lam)/dv
!
         m = m_2d
         mparity = MOD(m,2)
         mj = m+joff
         work1 = 0
         j1 = jmin1(m)
!
!        INVERSE TRANSFORM IN N-ZETA, FOR FIXED M
!
         n = n_2d
            ni = n+ioff
            DO k = 1, nzeta
               l = ns*(k-1)
               j1l = j1+l;  nsl = ns+l
               IF (logr)
     1         work1(j1l:nsl,1) = work1(j1l:nsl,1) 
     2                          + rmnsc(j1:ns,ni,mj)*cosnv(k,n)
               IF (logz)
     1         work1(j1l:nsl,6) = work1(j1l:nsl,6)
     2                          + zmncc(j1:ns,ni,mj)*cosnv(k,n)
               IF (logl)
     1         work1(j1l:nsl,10) = work1(j1l:nsl,10) 
     2                          + lmncc(j1:ns,ni,mj)*cosnv(k,n)

               IF (.not.lthreed) CYCLE
     
               IF (logr) THEN          
               work1(j1l:nsl,2) = work1(j1l:nsl,2) 
     1                          + rmncs(j1:ns,ni,mj)*sinnv(k,n)
               work1(j1l:nsl,3) = work1(j1l:nsl,3) 
     1                          + rmnsc(j1:ns,ni,mj)*sinnvn(k,n)
               work1(j1l:nsl,4) = work1(j1l:nsl,4) 
     1                          + rmncs(j1:ns,ni,mj)*cosnvn(k,n)
               END IF

               IF (logz) THEN
               work1(j1l:nsl,5) = work1(j1l:nsl,5) 
     1                          + zmnss(j1:ns,ni,mj)*sinnv(k,n)
               work1(j1l:nsl,7) = work1(j1l:nsl,7) 
     1                          + zmnss(j1:ns,ni,mj)*cosnvn(k,n)
               work1(j1l:nsl,8) = work1(j1l:nsl,8) 
     1                          + zmncc(j1:ns,ni,mj)*sinnvn(k,n)
               END IF

               IF (logl) THEN
               work1(j1l:nsl,9) = work1(j1l:nsl,9) 
     1                          + lmnss(j1:ns,ni,mj)*sinnv(k,n)
               work1(j1l:nsl,11) = work1(j1l:nsl,11)
     1                          + lmnss(j1:ns,ni,mj)*cosnvn(k,n)
               work1(j1l:nsl,12) = work1(j1l:nsl,12)
     1                          + lmncc(j1:ns,ni,mj)*sinnvn(k,n)
               END IF

            END DO
!
!        INVERSE TRANSFORM IN M-THETA, FOR ALL RADIAL, ZETA VALUES
!
         DO i = 1, ntheta2
            cosmux = xmpq(m,1)*cosmu(i,m)
            sinmux = xmpq(m,1)*sinmu(i,m)
            IF (logr) THEN
            r11(:,i,mparity) = r11(:,i,mparity) + work1(:,1)*
     1            sinmu(i,m)
            ru1(:,i,mparity) = ru1(:,i,mparity) + work1(:,1)*
     1            cosmum(i,m)
            rcn1(:,i,mparity) = rcn1(:,i,mparity) + work1(:,1)*
     1            sinmux
            END IF

            IF (logz) THEN
            z11(:,i,mparity) = z11(:,i,mparity) + work1(:,6)*
     1            cosmu(i,m)
            zu1(:,i,mparity) = zu1(:,i,mparity) + work1(:,6)*
     1            sinmum(i,m)
            zcn1(:,i,mparity) = zcn1(:,i,mparity) + work1(:,6)*
     1            cosmux
            END IF
 
            IF (logl)
     1      lu1(:,i,mparity) = lu1(:,i,mparity) + work1(:,10)*
     2            sinmum(i,m)

            IF (.not.lthreed) CYCLE

            IF (logr) THEN
            r11(:,i,mparity) = r11(:,i,mparity) + work1(:,2)*
     1            cosmu(i,m)
            ru1(:,i,mparity) = ru1(:,i,mparity) + work1(:,2)*
     1            sinmum(i,m)
            rcn1(:,i,mparity) = rcn1(:,i,mparity) + work1(:,2)*
     1            cosmux
            rv1(:,i,mparity) = rv1(:,i,mparity) + work1(:,3)*
     1            sinmu(i,m) + work1(:,4)*cosmu(i,m)
            END IF

            IF (logz) THEN
            z11(:,i,mparity) = z11(:,i,mparity) + work1(:,5)*
     1            sinmu(i,m)
            zu1(:,i,mparity) = zu1(:,i,mparity) + work1(:,5)*
     1            cosmum(i,m)
            zcn1(:,i,mparity) = zcn1(:,i,mparity) + work1(:,5)*
     1            sinmux
            zv1(:,i,mparity) = zv1(:,i,mparity) + work1(:,7)*
     1            sinmu(i,m) + work1(:,8)*cosmu(i,m)
            END IF

            IF (logl) THEN
            lu1(:,i,mparity) = lu1(:,i,mparity) + work1(:,9)*
     1            cosmum(i,m)
            lv1(:,i,mparity) = lv1(:,i,mparity) - (work1(:,11)*
     1            sinmu(i,m) + work1(:,12)*cosmu(i,m))
            END IF

         END DO
      
      DEALLOCATE (work1)

      END SUBROUTINE totzspa_hess


      SUBROUTINE convert_asym(rmnsc, zmncc)
      USE vmec_main, p5_2 => cp5
      USE precon2d, ONLY: ictrl_prec2d
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1), INTENT(inout) :: 
     1                                           rmnsc, zmncc
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), DIMENSION(ns,0:ntor) :: temp
C-----------------------------------------------
!
!     CONVERT FROM INTERNAL REPRESENTATION TO RMNSC, ZMNCC FOURIER FORM
!
      IF (.not.lconm1) RETURN

      temp(:,:) = rmnsc(:,:,1)
      rmnsc(:,:,1) = temp(:,:) + zmncc(:,:,1)
      zmncc(:,:,1) = temp(:,:) - zmncc(:,:,1)

      END SUBROUTINE convert_asym

      END MODULE totzsp_mod
