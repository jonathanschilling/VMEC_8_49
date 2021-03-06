      SUBROUTINE bcovar (lu, lv, lmnsc00)
      USE vmec_main, fpsi => bvco, p5 => cp5
      USE vmec_params, ONLY: ns4, signgs, pdamp, lamscale
      USE realspace
      USE vforces, r12 => armn_o, ru12 => azmn_e, gsqrt => azmn_o,
     1             rs => bzmn_e, zs => brmn_e, zu12 => armn_e,
     2             bsubu_e => clmn_e, bsubv_e => blmn_e, 
     3             bsubu_o => clmn_o, bsubv_o => blmn_o,
     4             bsq => bzmn_o, phipog => brmn_o
      USE vsvd, ONLY: phifac, phifsave, imovephi
      USE xstuff, ONLY: xc
      USE precon2d, ONLY: ictrl_prec2d, lHess_exact,
     1                    ctor_prec2d
      USE fbal
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), DIMENSION(nrzt,0:1), INTENT(inout) :: lu, lv
      REAL(rprec), DIMENSION(ns), INTENT(inout) :: lmnsc00
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!     GENERALLY, IF TEMPORAL CONVERGENCE IS POOR, TRY TO INCREASE PDAMP (< 1)
!     (STORED IN VMEC_PARAMS)
      REAL(rprec), PARAMETER :: c1p5 = (one + p5)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: l, js, ndim
      REAL(rprec) :: r2
      REAL(rprec) :: arnorm, aznorm, dnorm1, volume, tcon_mul
      REAL(rprec), POINTER, DIMENSION(:) :: luu, luv, lvv, tau
      REAL(rprec) :: curpol_temp
      REAL(rprec), DIMENSION(:), POINTER :: bsupu, bsubuh, 
     1                                      bsupv, bsubvh, r12sq
      LOGICAL :: lctor
!-----------------------------------------------
      ndim = 1+nrzt

!
!     POINTER ALIAS ASSIGNMENTS
!   
      tau => extra1(:,1);  luu => extra2(:,1);  
      luv => extra3(:,1);  lvv => extra4(:,1)

      bsupu => luu;  bsubuh => bsubu_o
      bsupv => luv;  bsubvh => bsubv_o
      r12sq => bsq


!
!     FOR OPTIMIZATION ON CRAY, MUST USE COMPILER DIRECTIVES TO
!     GET VECTORIZATION OF LOOPS INVOLVING (MORE THAN ONE) POINTER!
!
      guu(ndim) = 0;  guv = 0;  gvv = 0  

!
!     COMPUTE METRIC ELEMENTS GIJ ON HALF MESH
!     FIRST, GIJ = EVEN PART (ON FULL MESH), LIJ = ODD PART (ON FULL MESH)
!     THEN, GIJ(HALF) = < GIJ(even)> + SHALF < GIJ(odd) >
!

      r12sq(1:nrzt) = sqrts(1:nrzt)*sqrts(1:nrzt)
      guu(1:nrzt)   = ru(1:nrzt,0)*ru(1:nrzt,0) 
     1              + zu(1:nrzt,0)*zu(1:nrzt,0) + r12sq(1:nrzt)*
     2              ( ru(1:nrzt,1)*ru(1:nrzt,1)
     3            +   zu(1:nrzt,1)*zu(1:nrzt,1))

      luu(1:nrzt)   = (ru(1:nrzt,0)*ru(1:nrzt,1) 
     1              +  zu(1:nrzt,0)*zu(1:nrzt,1))*2
      phipog(1:nrzt)= 2*r1(1:nrzt,0)*r1(1:nrzt,1) 

      IF (lthreed) THEN
         guv(1:nrzt)   = ru(1:nrzt,0)*rv(1:nrzt,0)
     1                 + zu(1:nrzt,0)*zv(1:nrzt,0) + r12sq(1:nrzt)*
     2                 ( ru(1:nrzt,1)*rv(1:nrzt,1) 
     3                 + zu(1:nrzt,1)*zv(1:nrzt,1) )
         luv(1:nrzt)   = ru(1:nrzt,0)*rv(1:nrzt,1) 
     1                 + ru(1:nrzt,1)*rv(1:nrzt,0)
     2                 + zu(1:nrzt,0)*zv(1:nrzt,1) 
     3                 + zu(1:nrzt,1)*zv(1:nrzt,0)
         gvv(1:nrzt)   = rv(1:nrzt,0)*rv(1:nrzt,0)
     1                 + zv(1:nrzt,0)*zv(1:nrzt,0) + r12sq(1:nrzt)*
     2                 ( rv(1:nrzt,1)*rv(1:nrzt,1)
     3                 + zv(1:nrzt,1)*zv(1:nrzt,1) )
         lvv(1:nrzt)   =(rv(1:nrzt,0)*rv(1:nrzt,1) 
     1                 + zv(1:nrzt,0)*zv(1:nrzt,1))*2
      END IF

      r12sq(1:nrzt) = r1(1:nrzt,0)*r1(1:nrzt,0) + r12sq(1:nrzt)*
     1                r1(1:nrzt,1)*r1(1:nrzt,1)                       

!DIR$ IVDEP
      DO l = nrzt, 2, -1
         guu(l) = p5*(guu(l) + guu(l-1) + shalf(l)*(luu(l) + luu(l-1)))
         r12sq(l) = p5*(r12sq(l) + r12sq(l-1) + shalf(l)*             !Comment: r12sq = r12**2
     1                (phipog(l) + phipog(l-1)))                      
      END DO

      IF (lthreed) THEN
!DIR$ IVDEP
         DO l = nrzt, 2, -1
            guv(l) = p5*(guv(l) + guv(l-1) +
     1         shalf(l)*(luv(l) + luv(l-1)))
            gvv(l) = p5*(gvv(l) + gvv(l-1) +
     1         shalf(l)*(lvv(l) + lvv(l-1)))
         END DO
      END IF

      tau(1:nrzt) = gsqrt(1:nrzt)
      gsqrt(1:nrzt) = r12(1:nrzt)*tau(1:nrzt)      
      gsqrt(1:nrzt:ns) = gsqrt(2:nrzt:ns)

      gvv(2:nrzt) = gvv(2:nrzt) + r12sq(2:nrzt)

      phipog = 0
      WHERE (gsqrt(2:ndim) .ne. zero) 
     1       phipog(2:ndim) = one/gsqrt(2:ndim)
      phipog(1:ndim:ns) = 0

      vp(1) = 0;  vp(ns+1) = 0
      DO js = 2, ns
         vp(js) = signgs*SUM(gsqrt(js:nrzt:ns)*wint(js:nrzt:ns))
      END DO
      IF (iter2 .eq. 1) voli = twopi*twopi*hs*SUM(vp(2:ns))

!
!     COMPUTE CONTRA-VARIANT COMPONENTS OF B (Bsupu,v) ON RADIAL HALF-MESH
!     TO ACCOMODATE LRFP=T CASES, THE OVERALL PHIP FACTOR (PRIOR TO v8.46)
!     HAS BEEN REMOVED FROM PHIPOG, SO NOW PHIPOG == 1/GSQRT!
!
!     NOTE: LU = LAMU == d(LAM)/du, LV = -LAMV == -d(LAM)/dv COMING INTO THIS ROUTINE
!     WILL ADD CHIP IN CALL TO ADD_FLUXES. THE NET BSUPU, BSUPV ARE (PHIPOG=1/GSQRT AS NOTED ABOVE):
!
!          BSUPU = PHIPOG*(chip + LAMV),
!          BSUPV = PHIPOG*(phip + LAMU)
!
      lu = lu*lamscale
      lv = lv*lamscale

      DO js=1,ns
         lu(js:nrzt:ns,0)=lu(js:nrzt:ns,0)+phipf(js)
      END DO

      bsupu(2:nrzt) = p5*phipog(2:nrzt)*(lv(2:nrzt,0) + lv(1:nrzt-1,0) 
     1              + shalf(2:nrzt)*(lv(2:nrzt,1) + lv(1:nrzt-1,1)))
      bsupv(2:nrzt) = p5*phipog(2:nrzt)*(lu(2:nrzt,0) + lu(1:nrzt-1,0) 
     1              + shalf(2:nrzt)*(lu(2:nrzt,1) + lu(1:nrzt-1,1)))
      bsupu(1) =0;  bsupu(ndim) = 0
      bsupv(1) =0;  bsupv(ndim) = 0

!
!     UPDATE IOTA EITHER OF TWO WAYS:
!     1)  FOR ictrl_prec2d = 0, SOLVE THE LINEAR ALGEBRAIC EQUATION <Bsubu> = icurv 
!         FOR iotas  
!     2)  FOR ictrl_prec2d > 0, EVOLVE IOTAS IN TIME, USING Force-iota  = <Bsubu> - icurv. 
!
!     NEED TO DO IT WAY (#2) TO EASILY COMPUTE THE HESSIAN ELEMENTS DUE TO LAMBDA-VARIATIONS. 
!     IOTAS IS "STORED" AT LOCATION LAMBDA-SC(0,0) IN XC-ARRAY [USE THIS COMPONENT SO IT 
!     WILL WORK EVEN FOR 2D PLASMA], ALTHOUGH ITS VARIATION IS LIKE THAT OF LV-CS(0,0), 
!     WITH N -> 1 IN THE HESSIAN CALCULATION ROUTINES (Compute_Hessian_Flam_lam, etc.)
!

      IF (ncurr .eq. 1) THEN
        IF (ictrl_prec2d .eq. 2) THEN
            lmnsc00(2:ns) = iotas(2:ns)            !Initialize
         ELSE IF (ictrl_prec2d .ne. 0) THEN
            iotas(2:ns) = lmnsc00(2:ns)            !Evolution
         END IF
      END IF

!     COMPUTE (IF NEEDED) AND ADD CHIP TO BSUPU
      CALL add_fluxes(phipog, bsupu, bsupv, ictrl_prec2d.eq.0)

!
!     COMPUTE LAMBDA FORCE KERNELS (COVARIANT B COMPONENT bsubu,v) ON RADIAL HALF-MESH
!
      bsubuh(1:nrzt)=guu(1:nrzt)*bsupu(1:nrzt)+guv(1:nrzt)*bsupv(1:nrzt)
      bsubvh(1:nrzt)=guv(1:nrzt)*bsupu(1:nrzt)+gvv(1:nrzt)*bsupv(1:nrzt)
      bsubuh(ndim) = 0
      bsubvh(ndim) = 0

!
!     COMPUTE MAGNETIC AND KINETIC PRESSURE ON RADIAL HALF-MESH
!
      bsq(:nrzt) = p5*(bsupu(:nrzt)*bsubuh(:nrzt) 
     1           +     bsupv(:nrzt)*bsubvh(:nrzt))

      wb = hs*ABS(SUM(wint(:nrzt)*gsqrt(:nrzt)*bsq(:nrzt)))

!SPH: MAKE CALL HERE, USE lv FOR SCRATCH "gp" ARRAY
!      CALL an_pressure(lv,pperp_ns)
      pres(2:ns) = mass(2:ns)/vp(2:ns)**gamma
      wp = hs*SUM(vp(2:ns)*pres(2:ns))

!     COMPUTE AND ADD KINETIC PRESSURE TO MAGNETIC PRESSURE
!      bsq(2:nrzt)    = bsq(2:nrzt) + pperp(2:nrzt)
      DO js=2,ns
         bsq(js:nrzt:ns) = bsq(js:nrzt:ns) + pres(js)
      END DO

!WAC-SPH: MODIFY EFFECTIVE CURRENT K = curl(sigma_an*B)
!      phipog(1:nrzt) = phipog(1:nrzt)*sigma_an(1:nrzt)
!      bsubuh(1:nrzt) = bsubuh(1:nrzt)*sigma_an(1:nrzt)
!      bsubvh(1:nrzt) = bsubvh(1:nrzt)*sigma_an(1:nrzt)


!SPH122407-MOVED HERE: COMPUTE LAMBDA FULL MESH FORCES
      lvv = phipog(:ndim)*gvv
      bsubv_e(1:nrzt) = p5*(lvv(1:nrzt)+lvv(2:ndim))*lu(1:nrzt,0)

      lvv = lvv*shalf
      bsubu_e(:nrzt) = guv(:nrzt)*bsupu(:nrzt)   !*sigma_an(:nrzt)
      bsubu_e(ndim) = 0
      bsubv_e(1:nrzt) = bsubv_e(1:nrzt) 
     1            + p5*((lvv(1:nrzt) + lvv(2:ndim))*lu(1:nrzt,1)
     2            +      bsubu_e(1:nrzt) + bsubu_e(2:ndim))

!
!     COMPUTE AVERAGE FORCE BALANCE AND TOROIDAL/POLOIDAL CURRENTS
!
!WAC: UPDATE buco, bvco AFTER pressure called
      CALL calc_fbal(bsubuh, bsubvh)
    
      rbtor0= c1p5*fpsi(2)  - p5*fpsi(3)
      rbtor = c1p5*fpsi(ns) - p5*fpsi(ns-1)
!
!     (SPH:08/19/04)
!     MUST AVOID BREAKING TRI-DIAGONAL RADIAL COUPLING AT EDGE WHEN USING PRECONDITIONER
!     CTOR IS PASSED TO VACUUM TO COMPUTE EDGE BSQVAC, SO IT CAN ONLY DEPEND ON NS, NS-1
!     THUS, CTOR ~ buco(ns) WORKS, WITH REMAINDER A FIXED CONSTANT.
!
!     ALSO, IF USING FAST SWEEP IN COMPUTE_BLOCKS, MUST MAKE CTOR CONSTANT
!     TO AVOID BREAKING SYMMETRY OF A+(ns-1) AND B-(ns) HESSIAN ELEMENTS
!
!     TO GET CORRECT HESSIAN, USE THE CTOR=ctor_prec2d +... ASSIGNMENT
!     FOR ictrl_prec2d.ne.0 (replace ictrl_prec2d.gt.1 with ictrl_prec2d.ne.0 in IF test below)
!
!

!     NEXT COMPUTE COVARIANT BSUBV COMPONENT ~ lvv ON FULL RADIAL MESH BY AVERAGING HALF-MESH METRICS
!     NOTE: EDGE VALUES AT JS=NS DOWN BY 1/2
!     THIS IS NEEDED FOR NUMERICAL STABILITY

      IF (lHess_exact) THEN
         lctor = lfreeb .and. ictrl_prec2d.ne.0      !Yields correct hessian near edge
      ELSE
         lctor = lfreeb .and. ictrl_prec2d.gt.1      !Yields better accuracy in solution
      END IF
      IF (lctor) THEN       
         IF (ictrl_prec2d .eq. 2) 
     1       ctor_prec2d = signgs*twopi*p5*(buco(ns) - buco(ns1))
         ctor = ctor_prec2d + signgs*twopi*buco(ns)
      ELSE
         ctor = signgs*twopi*(c1p5*buco(ns) - p5*buco(ns1))
      END IF

!
!     AVERAGE LAMBDA FORCES ONTO FULL RADIAL MESH
!     USE BLENDING FOR bsubv_e FOR NUMERICAL STABILITY NEAR AXIS
!
      DO l=1,ns
         lvv(l:nrzt:ns) = bdamp(l)
      END DO
      bsubu_e(1:nrzt) = p5*(bsubuh(1:nrzt) + bsubuh(2:ndim))
      bsubv_e(1:nrzt) = bsubv_e(1:nrzt)*lvv(1:nrzt) + p5*(1-lvv(1:nrzt))
     1                *(bsubvh(1:nrzt) + bsubvh(2:ndim))

!
!     COMPUTE R,Z AND LAMBDA PRE-CONDITIONING MATRIX
!     ELEMENTS AND FORCE NORMS: NOTE THAT lu=>czmn, lv=>crmn externally
!     SO THIS STORES bsupv in czmn_e, bsupu in crmn_e
!
      IF (iequi.eq.1) THEN
         lu(:nrzt,0) = bsupv(:nrzt)
         lv(:nrzt,0) = bsupu(:nrzt)
      END IF
 
!
!     COMPUTE PRECONDITIONING (1D) AND SCALING PARAMETERS
!     NO NEED TO RECOMPUTE WHEN 2D-PRECONDITIONER ON
!
      IF ((MOD(iter2-iter1,ns4).eq.0 .and. iequi.eq.0) 
     1        .and. ictrl_prec2d.eq.0) THEN
         phifsave = phifac
         phipog(:nrzt) = phipog(:nrzt)*wint(:nrzt)
         CALL lamcal(phipog, guu, guv, gvv)
         CALL precondn(bsupv,bsq,gsqrt,r12,zs,zu12,zu,zu(1,1),
     1                 z1(1,1),arm,ard,brm,brd,crd,rzu_fac,cos01)
         CALL precondn(bsupv,bsq,gsqrt,r12,rs,ru12,ru,ru(1,1),
     1                 r1(1,1),azm,azd,bzm,bzd,crd,rru_fac,sin01)

         rzu_fac(2:ns-1) = sqrts(2:ns-1)*rzu_fac(2:ns-1)
         rru_fac(2:ns-1) = sqrts(2:ns-1)*rru_fac(2:ns-1)
         frcc_fac(2:ns-1) = one/rzu_fac(2:ns-1);  rzu_fac = rzu_fac/2
         fzsc_fac(2:ns-1) =-one/rru_fac(2:ns-1);  rru_fac = rru_fac/2

         guu(:nrzt) = guu(:nrzt)*r12(:nrzt)**2
         volume = hs*SUM(vp(2:ns))
         r2 = MAX(wb,wp)/volume
         fnorm = one/(SUM(guu(1:nrzt)*wint(1:nrzt))*(r2*r2))
         fnorm1 = one/SUM(xc(1+ns:2*irzloff)**2)                   !Norm for preconditioned R,Z forces
         dnorm1 = one/(nzeta*(ntheta2-1))
         fnormL = one/(dnorm1*SUM(bsubuh**2 + bsubvh**2))
!         r3 = one/(2*r0scale)**2
!         fnorm2 = one/MAX(SUM(xc(2*irzloff+1:3*irzloff)**2),r3/4) !Norm for preconditioned Lambda force

!
!        COMPUTE CONSTRAINT FORCE SCALING FACTOR (TCON)
!        OVERRIDE USER INPUT VALUE HERE
!
         r2 = ns
         tcon0 = MIN(ABS(tcon0), one)                              !!ignore large tcon0 from old-style files
         tcon_mul = tcon0*(1 + r2*(one/60 + r2/(200*120)))

         tcon_mul = tcon_mul/((4*r0scale**2)**2)                   !!Scaling of ard, azd (2*r0scale**2); 
                                                                   !!Scaling of cos**2 in alias (4*r0scale**2)

         DO js = 2, ns-1
           arnorm = SUM(wint(js:nrzt:ns)*ru0(js:nrzt:ns)**2)
           aznorm = SUM(wint(js:nrzt:ns)*zu0(js:nrzt:ns)**2)
           IF (arnorm.eq.zero .or. aznorm.eq.zero)
     1        STOP 'arnorm or aznorm=0 in bcovar'

           tcon(js) = MIN(ABS(ard(js,1)/arnorm),
     1                    ABS(azd(js,1)/aznorm)) * tcon_mul*(32*hs)**2
         END DO
         tcon(ns) = p5*tcon(ns-1)
      ENDIF

!
!     COMPUTE COVARIANT BSUBU,V (EVEN, ODD) ON HALF RADIAL MESH
!     FOR FORCE BALANCE AND RETURN (IEQUI=1)
!
!
!     COMPUTE COVARIANT BSUBU,V (EVEN, ODD) ON HALF RADIAL MESH
!     FOR FORCE BALANCE AND RETURN (IEQUI=1)
!
      IF (iequi .eq. 1) THEN

!         IF (.FALSE.) THEN
         DO js = ns-1,2,-1
            DO l = js, nrzt, ns
               bsubvh(l) = 2*bsubv_e(l) - bsubvh(l+1)
            END DO
         END DO

!     ADJUST <bsubvh> AFTER MESH-BLENDING
         DO js = 2, ns
            curpol_temp = fpsi(js) 
     1                  - SUM(bsubvh(js:nrzt:ns)*wint(js:nrzt:ns))
           DO l = js, nrzt, ns
               bsubvh(l) = bsubvh(l) + curpol_temp
            END DO
         END DO
!         END IF

         bsubu_e(:nrzt) = bsubuh(:nrzt)
         bsubv_e(:nrzt) = bsubvh(:nrzt)

!         DO js = 2, ns
!            WRITE (100, *) 'JS: ', js, 'BSUBU, BSUBV'
!            WRITE (100, '(1p,6e12.4)') 
!     1          bsubu_e(js:nrzt:ns), bsubvh(js:nrzt:ns)
!         END DO

         bsubu_o(:nrzt) = shalf(:nrzt)*bsubu_e(:nrzt)
         bsubv_o(:nrzt) = shalf(:nrzt)*bsubv_e(:nrzt)

         RETURN

      END IF

!     MINUS SIGN => HESSIAN DIAGONALS ARE POSITIVE
      bsubu_e = -bsubu_e
      bsubv_e = -bsubv_e
      bsubu_o(:nrzt)  = sqrts(:nrzt)*bsubu_e(:nrzt)
      bsubv_o(:nrzt)  = sqrts(:nrzt)*bsubv_e(:nrzt)
!
!     STORE LU * LV COMBINATIONS USED IN FORCES
!
!WAC, SPH122407
      lvv(2:nrzt) = gsqrt(2:nrzt)  !*sigma_an(2:nrzt)
      guu(2:nrzt)  = bsupu(2:nrzt)*bsupu(2:nrzt)*lvv(2:nrzt)
      guv(2:nrzt)  = bsupu(2:nrzt)*bsupv(2:nrzt)*lvv(2:nrzt)
      gvv(2:nrzt)  = bsupv(2:nrzt)*bsupv(2:nrzt)*lvv(2:nrzt)
      lv(2:nrzt,0) = bsq(2:nrzt)*tau(2:nrzt)
      lu(2:nrzt,0) = bsq(2:nrzt)*r12(2:nrzt)

      END SUBROUTINE bcovar

