      SUBROUTINE jxbforce(bsupu, bsupv, bsubu, bsubv, bsubs, bsubsu,
     1   bsubsv, gsqrt, bsq, itheta, izeta, ier_flag)
      USE safe_open_mod
      USE vmec_main
      USE vmec_params, ONLY: mscale, nscale, signgs, mnyq, nnyq, 
     1                       successful_term_flag
      USE realspace
!!#undef NETCDF     !!undefine NETCDF if text output desired
#ifdef NETCDF
      USE ezcdf
#endif
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), DIMENSION(ns,nznt), INTENT(in) ::
     1  bsupu, bsupv, bsq, gsqrt
!ONLY CHANGE COMPONENT 1 (MUST KEEP COMPONENT 0 UNCHANGED FOR USE IN WROUT!)
      REAL(rprec), DIMENSION(ns,nznt,0:1), INTENT(inout), TARGET ::
     1  bsubu, bsubv
      REAL(rprec), DIMENSION(ns,nznt), INTENT(inout), TARGET :: bsubs
      REAL(rprec), DIMENSION(ns,nznt), INTENT(out) :: itheta, izeta
      REAL(rprec), DIMENSION(ns,nznt,0:1) :: bsubsu, bsubsv
      INTEGER, INTENT(in) :: ier_flag
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      LOGICAL, PARAMETER :: lprint = .false.      !!Prints out bsubs spectrum to fort.33
      REAL(rprec), PARAMETER :: two=2, p5=0.5_dp, c1p5=1.5_dp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER lk, lz, lt, k, m, js, j, n, injxbout, mparity, nznt1
      INTEGER :: njxbout = jxbout0, info, mtrig, ntrig
      INTEGER, PARAMETER :: ns_skip = 1, nu_skip = 1, nv_skip = 1
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::
     1    bdotj, bsubuv, bsubvu
      REAL(rprec), DIMENSION(:,:,:,:), ALLOCATABLE :: bsubsmn
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE ::   brhomn,
     1     bsubs3, bsubv3, bsubu3, jxb_gradp, jcrossb, sqrtg3,
     2     bsupv3, bsupu3, jsups3, jsupv3, jsupu3, jdotb_sqrtg
      REAL(rprec), POINTER :: bs1(:), bu1(:,:), bv1(:,:)
      REAL(rprec), DIMENSION(:), ALLOCATABLE     :: jperpu, jperpv, 
     2    sqgb2, sqrtg, jp2, jxb, jxb2, bsupu1, bsupv1, bsubu1, bsubv1,
     3    avforce, aminfor, amaxfor, toroidal_angle, phin, pprim
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE   :: bsubua, bsubva
      REAL(rprec) ::
     1    bsubsmn1, bsubsmn2, bsubvmn1, bsubvmn2, bsubumn1, bsubumn2,
     1    bsubsmn3, bsubsmn4, bsubvmn3, bsubvmn4, bsubumn3, bsubumn4,
     2    dnorm1, tcos1, tcos2, tsini1, tsini2, tcosi1, tcosi2, 
     3    tcosm1, tcosm2, tcosn1, tcosn2, tsinm1, tsinm2, tsin1, tsin2,
     4    tsinn1, tsinn2, pprime, tjnorm, ovp, pnorm, brho00(ns)
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::
     1    bsubu_s, bsubu_a, bsubv_s, bsubv_a
      REAL(rprec), DIMENSION(:), ALLOCATABLE ::
     1    bsubs_s, bsubs_a
      CHARACTER(LEN=100) :: jxbout_file
      CHARACTER(LEN=100) :: legend(13)
      LOGICAL :: lprint_flag
!-----------------------------------------------
#ifdef NETCDF
      CHARACTER(LEN=*), PARAMETER ::
     1  vn_legend = 'legend', vn_radial_surfaces = 'radial_surfaces',
     2  vn_poloidal_grid_points = 'poloidal_grid_points',
     3  vn_toroidal_grid_points = 'toroidal_grid_points',
     4  vn_mpol = 'mpol', vn_ntor = 'ntor', vn_phin = 'phin',
     5  vn_toroidal_angle = 'toroidal_angle',
     6  vn_avforce = 'avforce', vn_jdotb = 'surf_av_jdotb',
     7  vn_sqg_bdotj = 'sqrt(g)*jdotb', vn_sqrtg = 'sqrt(g)',
     8  vn_bdotgradv = 'bdotgradv', vn_amaxfor = 'amaxfor',
     9  vn_aminfor = 'aminfor', vn_pprime = 'pprime',
     A  vn_jsupu = 'jsupu', vn_jsupv = 'jsupv',
     B  vn_jsups = 'jsups', vn_bsupu = 'bsupu',
     C  vn_bsupv = 'bsupv', vn_jcrossb = 'jcrossb',
     D  vn_jxb_gradp = 'jxb_gradp', vn_bsubu = 'bsubu',
     E  vn_bsubv = 'bsubv', vn_bsubs = 'bsubs'
!-----------------------------------------------
#endif
      lprint_flag = (ier_flag.eq.successful_term_flag)
      IF (lprint_flag) THEN
#ifdef NETCDF
      jxbout_file = 'jxbout_'//TRIM(input_extension)//'.nc'

      CALL cdf_open(njxbout,jxbout_file,'w',injxbout)
#else
      jxbout_file = 'jxbout_'//TRIM(input_extension)//'.txt'
      CALL safe_open(njxbout, injxbout, jxbout_file, 'replace',
     1    'formatted')
#endif
      IF (injxbout .ne. 0) THEN
         PRINT *,' Error opening JXBOUT file in jxbforce'
         RETURN
      END IF

!
!     PROGRAM FOR COMPUTING LOCAL JXB = grad-p FORCE BALANCE
!
!     Compute u (=theta), v (=zeta) derivatives of B sub s
!
      legend(1) = " S = normalized toroidal flux (0 - 1)"
      IF (lasym) THEN
         legend(2) = " U = VMEC poloidal angle (0 - 2*pi, FULL period)"
      ELSE
         legend(2) = " U = VMEC poloidal angle (0 - pi, HALF a period)"
      END IF
      legend(3) = " V = VMEC (geometric) toroidal angle (0 - 2*pi)"
      legend(4) = " SQRT(g') = |SQRT(g-VMEC)| / VOL':" //
     1  " Cylindrical-to-s,u,v Jacobian normed to volume derivative"
      legend(5) = " VOL = Int_s'=0,s Int_u Int_v |SQRT(g_VMEC)| :" //
     1  " plasma volume  enclosed by surface s'=s"
      legend(6) = " VOL' = d(VOL)/ds: differential volume element"
      legend(7) = " Es = SQRT(g') [grad(U) X grad(V)] : covariant" //
     1   " radial unit vector (based on volume radial coordinate)"
      legend(8) =
     1  " BSUP{U,V} = B DOT GRAD{U,V}:  contravariant components of B"
      legend(9) = " JSUP{U,V} = SQRT(g') J DOT GRAD{U,V}"
      legend(10)=
     1  " J X B = Es DOT [J X B]: covariant component of J X B force"
      legend(11)= " J * B = J DOT B * SQRT(g')"
      legend(12)= " p' = dp/d(VOL): pressure gradient (based on" //
     1  " volume radial coordinate)"
      legend(13)= " <JSUP{U,V}> = Int_u Int_v [JSUP{U,V}]/dV/ds"

#ifndef NETCDF
      WRITE (njxbout,5) (ns1-1)/ns_skip, ntheta3/nu_skip, nzeta/nv_skip,
     1    mpol, ntor, phiedge
 5    FORMAT(/,' Radial surfaces = ',i3, ' Poloidal grid points = ',i3,
     1         ' Toroidal grid points = ',i3,/,
     2         ' Poloidal modes = ',i3,' Toroidal modes = ', i3,
     3         ' Toroidal Flux  = ',1pe12.3)
      WRITE (njxbout, 6) (legend(j), j=1,13)
 6    FORMAT(/,100('-'),/,' LEGEND:',/,100('-'),/,
     1  2(3(a,/),/),5(a,/),/,2(a,/),100('-'),//)
#endif
      ENDIF

      nznt1 = nzeta*ntheta2
      ALLOCATE (avforce(ns),aminfor(ns),amaxfor(ns))
      ALLOCATE (bdotj(ns,nznt), bsubuv(ns,nznt),
     1          bsubvu(ns,nznt), jperpu(nznt), jperpv(nznt), 
     2          sqgb2(nznt), brhomn(0:mnyq,-nnyq:nnyq,0:1),jp2(nznt),
     3          jxb(nznt), jxb2(nznt), bsupu1(nznt), 
     3          bsubua(nznt1,0:1), bsubva(nznt1,0:1), 
     4          bsupv1(nznt), bsubu1(nznt), bsubv1(nznt),
     5          bsubsmn(ns,0:mnyq,-nnyq:nnyq,0:1),
     6          bsubs_s(nznt), bsubs_a(nznt), sqrtg(nznt),
     7          bsubu_s(nznt1,0:1), bsubu_a(nznt1,0:1),
     8          bsubv_s(nznt1,0:1), bsubv_a(nznt1,0:1), stat=j)
      IF (j .ne. 0) STOP 'Allocation error in jxbforce'

!
!     NOTE: bsubuv, bsubvu are used to compute the radial current (should be zero)
!
      bsubsu = 0; bsubsv = 0; bsubuv = 0; bsubvu = 0; bdotj  = 0
      bsubs(1,:) = 0; bsubsmn = 0

      radial: DO js = 1, ns
!
!     Put bsubs on full mesh
!
         IF (js.gt.1 .and. js.lt.ns) THEN     
            bsubs(js,:) = p5*(bsubs(js,:) + bsubs(js+1,:))
         END IF

         bsubu(js,:,1) = bsubu(js,:,1)/shalf(js)
         bsubv(js,:,1) = bsubv(js,:,1)/shalf(js)
         bsubua = 0;   bsubva = 0

!        _s: symmetric in u,v  _a: antisymmetric in u,v on half (ntheta2) interval
         IF (lasym)  THEN
            bs1=>bsubs(js,:); bu1=>bsubu(js,:,:); bv1=>bsubv(js,:,:)
            CALL fsym_fft (bs1, bu1, bv1, bsubs_s, bsubu_s, bsubv_s, 
     1                     bsubs_a, bsubu_a, bsubv_a)
         ELSE
            bsubs_s(:) = bsubs(js,:)
            bsubu_s = bsubu(js,:,:); bsubv_s = bsubv(js,:,:)
         END IF

!
!        FOURIER LOW-PASS FILTER bsubX

         DO m = 0, mpol1
            mparity = MOD(m, 2)
            DO n = 0, ntor
!
!        FOURIER TRANSFORM
!
               dnorm1 = one/r0scale**2
               IF (m .eq. mnyq) dnorm1 = p5*dnorm1
               IF (n.eq.nnyq .and. n.ne.0) dnorm1 = p5*dnorm1
               bsubsmn1 = 0;  bsubsmn2 = 0
               IF (lasym) THEN
                  bsubsmn3 = 0;  bsubsmn4 = 0
               END IF
               bsubumn1 = 0;  bsubumn2 = 0;  bsubvmn1 = 0;  bsubvmn2 = 0
               IF (lasym) THEN
                  bsubumn3 = 0; bsubumn4 = 0; bsubvmn3 = 0; bsubvmn4 = 0
               END IF

               DO k = 1, nzeta
                  lk = k
                  DO j = 1, ntheta2
                     tsini1 = sinmui(j,m)*cosnv(k,n)*dnorm1
                     tsini2 = cosmui(j,m)*sinnv(k,n)*dnorm1
                     tcosi1 = cosmui(j,m)*cosnv(k,n)*dnorm1
                     tcosi2 = sinmui(j,m)*sinnv(k,n)*dnorm1
                     bsubsmn1 = bsubsmn1 + tsini1*bsubs_s(lk)
                     bsubsmn2 = bsubsmn2 + tsini2*bsubs_s(lk)
                     bsubvmn1 = bsubvmn1 + tcosi1*bsubv_s(lk, mparity)
                     bsubvmn2 = bsubvmn2 + tcosi2*bsubv_s(lk, mparity)
                     bsubumn1 = bsubumn1 + tcosi1*bsubu_s(lk, mparity)
                     bsubumn2 = bsubumn2 + tcosi2*bsubu_s(lk, mparity)

                     IF (lasym) THEN
                     bsubsmn3 = bsubsmn3 + tcosi1*bsubs_a(lk)
                     bsubsmn4 = bsubsmn4 + tcosi2*bsubs_a(lk)
                     bsubvmn3 = bsubvmn3 + tsini1*bsubv_a(lk, mparity)
                     bsubvmn4 = bsubvmn4 + tsini2*bsubv_a(lk, mparity)
                     bsubumn3 = bsubumn3 + tsini1*bsubu_a(lk, mparity)
                     bsubumn4 = bsubumn4 + tsini2*bsubu_a(lk, mparity)
                     END IF

                     lk = lk + nzeta

                  END DO
               END DO

!
!              FOURIER INVERSE TRANSFORM
!              Compute on u-v grid (must add symmetric, antisymmetric parts for lasym=T)
! 
               DO k = 1, nzeta
                  lk = k
                  DO j = 1, ntheta2
                     tcos1 = cosmu(j,m)*cosnv(k,n)
                     tcos2 = sinmu(j,m)*sinnv(k,n)
                     bsubua(lk,0) = bsubua(lk,0) + tcos1*bsubumn1 +
     1                  tcos2*bsubumn2
                     bsubva(lk,0) = bsubva(lk,0) + tcos1*bsubvmn1 +
     1                  tcos2*bsubvmn2

                     tcosm1 = cosmum(j,m)*cosnv(k,n)
                     tcosm2 = sinmum(j,m)*sinnv(k,n)
                     bsubsu(js,lk,0) = bsubsu(js,lk,0) +
     1                  tcosm1*bsubsmn1 + tcosm2*bsubsmn2
                     tcosn1 = sinmu(j,m)*sinnvn(k,n)
                     tcosn2 = cosmu(j,m)*cosnvn(k,n)
                     bsubsv(js,lk,0) = bsubsv(js,lk,0) +
     1                  tcosn1*bsubsmn1 + tcosn2*bsubsmn2
                     bsubvu(js,lk) = bsubvu(js,lk) + 
     1                               sinmum(j,m)*cosnv(k,n)*bsubvmn1 +
     2                               cosmum(j,m)*sinnv(k,n)*bsubvmn2
                     bsubuv(js,lk) = bsubuv(js,lk) + 
     1                               cosmu(j,m)*sinnvn(k,n)*bsubumn1 +
     2                               sinmu(j,m)*cosnvn(k,n)*bsubumn2

                     IF (lasym) THEN
                     tsin1 = sinmu(j,m)*cosnv(k,n)
                     tsin2 = cosmu(j,m)*sinnv(k,n)
                     bsubua(lk,1) = bsubua(lk,1) + tsin1*bsubumn3 +
     1                  tsin2*bsubumn4
                     bsubva(lk,1) = bsubva(lk,1) + tsin1*bsubvmn3 +
     1                  tsin2*bsubvmn4

                     tsinm1 = sinmum(j,m)*cosnv(k,n)
                     tsinm2 = cosmum(j,m)*sinnv(k,n)
                     bsubsu(js,lk,1) = bsubsu(js,lk,1) +
     1                   tsinm1*bsubsmn3 + tsinm2*bsubsmn4
                     tsinn1 = cosmu(j,m)*sinnvn(k,n)
                     tsinn2 = sinmu(j,m)*cosnvn(k,n)
                     bsubsv(js,lk,1) = bsubsv(js,lk,1) +
     1                   tsinn1*bsubsmn3 + tsinn2*bsubsmn4
                     bsubvu(js,lk) = bsubvu(js,lk) + 
     1                               cosmum(j,m)*cosnv(k,n)*bsubvmn3 +
     2                               sinmum(j,m)*sinnv(k,n)*bsubvmn4
                     bsubuv(js,lk) = bsubuv(js,lk) + 
     1                               sinmu(j,m)*sinnvn(k,n)*bsubumn3 +
     2                               cosmu(j,m)*cosnvn(k,n)*bsubumn4
                     END IF

                     lk = lk + nzeta

                  END DO
               END DO

!
!              bsubsmn: coefficients of sin(mu)cos(nv), n>=0, cos(mu)sin(nv), n<0 (type=0)
!                                       cos(mu)cos(nv), n>=0, sin(mu)sin(nv), n<0 (type=1, nonzero only for lasym=T)
!
               IF (.not.lprint) CYCLE          !Don't need these except for comparison

               bsubsmn(js,m,n,0) = bsubsmn1
               IF (n .gt. 0) bsubsmn(js,m,-n,0) = bsubsmn2
             
               IF (.not.lasym) CYCLE

               bsubsmn(js,m,n,1) = bsubsmn3
               IF (n .gt. 0) bsubsmn(js,m,-n,0) = bsubsmn4

            END DO
         END DO

         IF (lasym) THEN
!           EXTEND FILTERED bsubu, bsubv TO NTHETA3 MESH AND STORE IN COMPONENT=1
!           NOTE: INDEX 0 - COS(mu-nv) SYMMETRY; 1 - SIN(mu-nv) SYMMETRY
            CALL fext_fft (bsubu(js,:,1), bsubua(:,0), bsubua(:,1))
            CALL fext_fft (bsubv(js,:,1), bsubva(:,0), bsubva(:,1))
         ELSE
            bsubu(js,:,1) = bsubua(:,0)
            bsubv(js,:,1) = bsubva(:,0)
         END IF

      END DO radial

      DEALLOCATE (bsubua, bsubva)

!     EXTEND bsubsu, bsubsv TO NTHETA3 MESH
      IF (lasym) CALL fsym_invfft (bsubsu, bsubsv)

      DEALLOCATE (bsubs_s, bsubs_a, bsubu_s,
     1            bsubu_a, bsubv_s, bsubv_a, stat=lk)

!
!     Compute end point values for bsubs
!
      bsubs(1,:)  = 2*bsubs(2,:)  - bsubs(3,:)
      bsubs(ns,:) = 2*bsubs(ns,:) - bsubs(ns-1,:)
!
!     Now compute currents on the FULL radial mesh
!     Here:
!
!     Itheta = sqrt(g) * Jsupu
!     Izeta  = sqrt(g) * Jsupv
!     Jsupx  = J dot grad(x)                          x=(u,v)
!     jxb    = (J X B) dot (grad-u X grad-v) sqrt(g)  
!     bdotj  = sqrt(g)*J dot B
!     jperpx = (B X gradp) dot grad(x) / |B|**2       x=(u,v)
!     sqgb2  = sqrt(g)*|B|**2
!     sqrtg  = sqrt(g)
!     pprime = dp/dV
!
!     jp2   == |j-perp|**2 = jperpu**2 * guu + 2*jperpu*jperpv*guv + jperv**2 * gvv
!     This was compared to the alternative expression (agreed very well):
!     |j-perp|**2 = |grad-s|**2 * (dp/ds)**2 / |B|**2
!
!     Note: Multiply currents, pressure by 1/mu0 to get in mks units!
!           TWOPI*TWOPI factor incorporated in vp (thru ovp factor below), so V' = (2pi)**2*vp
!
#ifdef NETCDF
      ALLOCATE(
     1     bsubs3(ns,nzeta,ntheta3), bsubv3(ns,nzeta,ntheta3), 
     2     bsubu3(ns,nzeta,ntheta3), jxb_gradp(ns,nzeta,ntheta3), 
     3     jcrossb(ns,nzeta,ntheta3), bsupv3(ns,nzeta,ntheta3), 
     4     bsupu3(ns,nzeta,ntheta3), jsups3(ns,nzeta,ntheta3), 
     5     jsupv3(ns,nzeta,ntheta3), jsupu3(ns,nzeta,ntheta3),
     6     jdotb_sqrtg(ns,nzeta,ntheta3), sqrtg3(ns,nzeta,ntheta3), 
     7     phin(ns), toroidal_angle(nzeta), stat=j)

      bsubs3=0; bsubv3=0; bsubu3=0; jxb_gradp=0
      jcrossb=0 ; bsupv3=0; bsupu3=0; jsups3=0
      jsupv3=0; jsupu3=0; phin=0; phin(ns)=1
      jdotb_sqrtg=0; sqrtg3=0 
#endif 
      ALLOCATE (pprim(ns),stat=j)
      pprim=0

      avforce=0; aminfor=0; amaxfor=0
      dnorm1 = twopi*twopi

      DO js = 2, ns1
         ovp = two/(vp(js+1) + vp(js))/dnorm1
         tjnorm = ovp*signgs
         sqgb2(:nznt) = (gsqrt(js+1,:nznt)*(bsq(js+1,:nznt)- pres(js+1))
     1                +  gsqrt(js,:nznt)  *(bsq(js,:nznt) - pres(js)))
!        TAKE THIS OUT: MAY BE POORLY CONVERGED AT THIS POINT....
!         IF (ANY(sqgb2(:nznt)*signgs .le. zero)) 
!     1       STOP ' SQGB2 <= 0 in JXBFORCE'
         pprime = ohs*(pres(js+1)-pres(js))/mu0              !dp/ds here

!        NOTE: COMPONENT 1 OF bsubX, X=u,v, IS THE FILTERED COVARIANT COMPONENT 
!              COMPONENT 0 IS UNFILTERED 
         jperpu(:nznt) = p5*(bsubv(js+1,:nznt,1) + bsubv(js,:nznt,1))*
     1                       pprime/sqgb2
         jperpv(:nznt) =-p5*(bsubu(js+1,:nznt,1) + bsubu(js,:nznt,1))*
     1                       pprime/sqgb2
         bsubu1(:nznt) = p5*(bsubu(js+1,:nznt,1) + bsubu(js,:nznt,1))
         bsubv1(:nznt) = p5*(bsubv(js+1,:nznt,1) + bsubv(js,:nznt,1))

         jp2(:nznt)=p5*(jperpu**2*(guu(js+1:nrzt:ns) + guu(js:nrzt:ns))
     1          + 2*jperpu*jperpv*(guv(js+1:nrzt:ns) + guv(js:nrzt:ns))
     2          +       jperpv**2*(gvv(js+1:nrzt:ns) + gvv(js:nrzt:ns)))
         itheta(js,:nznt) =  bsubsv(js,:nznt,0) - ohs*
     1                      (bsubv(js+1,:nznt,1) - bsubv(js,:nznt,1))
         izeta(js,:nznt)  = -bsubsu(js,:nznt,0) + ohs*
     1                      (bsubu(js+1,:nznt,1) - bsubu(js,:nznt,1))
         itheta(js,:nznt) = itheta(js,:nznt)/mu0
         izeta(js,:nznt)  = izeta(js,:nznt)/mu0
         sqrtg(:) = p5*(gsqrt(js,:) + gsqrt(js+1,:))
         bsupu1(:nznt) = p5*(bsupu(js+1,:nznt)*gsqrt(js+1,:)
     1                 +     bsupu(js,:nznt)  *gsqrt(js,:)) / sqrtg(:)
         bsupv1(:nznt) = p5*(bsupv(js+1,:nznt)*gsqrt(js+1,:)
     1                 +     bsupv(js,:nznt)  *gsqrt(js,:)) / sqrtg(:)

         jxb(:nznt) = ovp*(itheta(js,:nznt) * bsupv1(:nznt)
     1              -      izeta (js,:nznt) * bsupu1(:nznt))
         bdotj(js,:nznt) = itheta(js,:nznt) * bsubu1(:nznt) +
     1                     izeta (js,:nznt) * bsubv1(:nznt)
         pprime = ovp*pprime
         pnorm = one/(ABS(pprime) + EPSILON(pprime))
         amaxfor(js) = MAXVAL(jxb(1:nznt)-pprime)*pnorm
         aminfor(js) = MINVAL(jxb(1:nznt)-pprime)*pnorm
         avforce(js) = SUM(wint(2:nrzt:ns)*(jxb(:nznt) - pprime))
         amaxfor(js) = 100*MIN(amaxfor(js),9.999_dp)
         aminfor(js) = 100*MAX(aminfor(js),-9.999_dp)
         pprim(js) = pprime
!        Compute <j dot B>, <B sup v> = signgs*phip
!        jpar2 = <j||**2>, jperp2 = <j-perp**2>,  with <...> = flux surface average

         jdotb(js) = dnorm1*tjnorm*SUM(bdotj(js,:nznt)*wint(2:nrzt:ns))
         bdotb(js) = dnorm1*tjnorm*SUM(sqgb2(:nznt)*wint(2:nrzt:ns))
              
         bdotgradv(js) = p5*dnorm1*tjnorm*(phip(js) + phip(js+1))
         jpar2(js) = dnorm1*tjnorm*
     1            SUM(bdotj(js,:nznt)**2*wint(2:nrzt:ns)/sqgb2(:nznt))
         jperp2(js)= dnorm1*tjnorm*
     1            SUM(jp2(:nznt)*wint(2:nrzt:ns)*sqrtg(:nznt))

         IF (MOD(js,ns_skip) .eq. 0 .and. lprint_flag) THEN
#ifdef NETCDF
            phin(js) = phi(js)/phi(ns)
            DO lz = 1, nzeta
               toroidal_angle(lz)=REAL(360*(lz-1),rprec)/nzeta
               DO lt = 1, ntheta3
                  lk = lz + nzeta*(lt-1)
                  jsupu3 (js,lz,lt) = ovp*itheta(js,lk)
                  jsupv3 (js,lz,lt) = ovp*izeta(js,lk) 
                  jsups3 (js,lz,lt) = ovp*(bsubuv(js,lk)
     1                              -      bsubvu(js,lk))/mu0
                  bsupu3 (js,lz,lt) = bsupu1(lk)
                  bsupv3 (js,lz,lt) = bsupv1(lk)
                  jcrossb (js,lz,lt) = jxb(lk)
                  jxb_gradp (js,lz,lt) = (jxb(lk) - pprime)
                  jdotb_sqrtg (js,lz,lt) = ovp*bdotj(js,lk) 
                  sqrtg3(js,lz,lt) = sqrtg(lk)*ovp
                  bsubu3(js,lz,lt) = bsubu(js,lk,1)
                  bsubv3(js,lz,lt) = bsubv(js,lk,1)
                  bsubs3(js,lz,lt) = bsubs(js,lk)
               END DO
            END DO
#else
            WRITE (njxbout, 200) phi(js), avforce(js), jdotb(js),
     1         bdotgradv(js), pprime, one/ovp, 
     2         (twopi**2)*tjnorm*SUM(itheta(js,:)*wint(js:nrzt:ns)),
     3         (twopi**2)*tjnorm*SUM(izeta (js,:)*wint(js:nrzt:ns)),
     4         amaxfor(js), aminfor(js)
            WRITE (njxbout, 90)
            DO lz = 1, nzeta, nv_skip
               WRITE (njxbout, 100) REAL(360*(lz-1),rprec)/nzeta, lz
               DO lt = 1, ntheta3, nu_skip
                  lk = lz + nzeta*(lt - 1)
                  WRITE (njxbout, 110) lt, tjnorm*itheta(js,lk),
     1              tjnorm*izeta(js,lk), ovp*(bsubuv(js,lk) -
     2              bsubvu(js,lk))/mu0, bsupu1(lk), bsupv1(lk),
     3              sqrtg(lk)*ovp, jxb(lk), jxb(lk) - pprime,
     4              ovp*bdotj(js,lk), bsubu(js,lk,1),
     5              bsubv(js,lk,1), bsubs(js,lk)
               END DO
            END DO
#endif
         ENDIF
      END DO

      izeta(1,:nznt) = two*izeta(2,:nznt) - izeta(3,:nznt)           !!Need in wrout
      izeta(ns,:nznt)= two*izeta(ns-1,:nznt) - izeta(ns-2,:nznt)     !!Need in wrout
      jdotb(1) = two*jdotb(2) - jdotb(3)
      jdotb(ns) = two*jdotb(ns-1) - jdotb(ns-2)
      bdotb(1) = two*bdotb(3) - bdotb(2)
      bdotb(ns) = two*bdotb(ns-1) - bdotb(ns-2)
      bdotgradv(1) = two*bdotgradv(2) - bdotgradv(3)
      bdotgradv(ns) = two*bdotgradv(ns-1) - bdotgradv(ns-2)
      jpar2(1)   = 0; jpar2(ns)  = 0; jperp2(1)  = 0; jperp2(ns) = 0
      pprim(1) = 2*pprim(ns-1) - pprim(ns-2)
      pprim(ns) = 2*pprim(ns-1) - pprim(ns-2)

      IF (lprint_flag) THEN
#ifdef NETCDF
      CALL cdf_define(njxbout,vn_legend,legend)
      CALL cdf_define(njxbout,vn_mpol,mpol)
      CALL cdf_define(njxbout,vn_ntor,ntor)
      CALL cdf_define(njxbout,vn_phin,phin)
      CALL cdf_define(njxbout,vn_radial_surfaces,ns)
      CALL cdf_define(njxbout,vn_poloidal_grid_points,ntheta3)
      CALL cdf_define(njxbout,vn_toroidal_grid_points,nzeta)
      CALL cdf_define(njxbout,vn_avforce,avforce)
      CALL cdf_define(njxbout,vn_jdotb,jdotb)
 
      CALL cdf_define(njxbout,vn_sqg_bdotj,jdotb_sqrtg)
      CALL cdf_define(njxbout,vn_sqrtg,sqrtg3)

      CALL cdf_define(njxbout,vn_bdotgradv,bdotgradv)
      CALL cdf_define(njxbout,vn_pprime,pprim)
      CALL cdf_define(njxbout,vn_aminfor,aminfor)
      CALL cdf_define(njxbout,vn_amaxfor,amaxfor)
      CALL cdf_define(njxbout,vn_jsupu,jsupu3)
      CALL cdf_define(njxbout,vn_jsupv,jsupv3)
      CALL cdf_define(njxbout,vn_jsups,jsups3)
      CALL cdf_define(njxbout,vn_bsupu,bsupu3)
      CALL cdf_define(njxbout,vn_bsupv,bsupv3)
      CALL cdf_define(njxbout,vn_jcrossb,jcrossb)
      CALL cdf_define(njxbout,vn_jxb_gradp,jxb_gradp)
      CALL cdf_define(njxbout,vn_bsubu,bsubu3)
      CALL cdf_define(njxbout,vn_bsubv,bsubv3)
      CALL cdf_define(njxbout,vn_bsubs,bsubs3)

      CALL cdf_write(njxbout,vn_legend,legend)
      CALL cdf_write(njxbout,vn_mpol,mpol)
      CALL cdf_write(njxbout,vn_ntor,ntor)
      CALL cdf_write(njxbout,vn_phin,phin)
      CALL cdf_write(njxbout,vn_radial_surfaces,ns)
      CALL cdf_write(njxbout,vn_poloidal_grid_points,ntheta3)
      CALL cdf_write(njxbout,vn_toroidal_grid_points,nzeta)
      CALL cdf_write(njxbout,vn_avforce,avforce)
      CALL cdf_write(njxbout,vn_jdotb,jdotb)

      CALL cdf_write(njxbout,vn_sqg_bdotj,jdotb_sqrtg)
      CALL cdf_write(njxbout,vn_sqrtg,sqrtg3)

      CALL cdf_write(njxbout,vn_bdotgradv,bdotgradv)
      CALL cdf_write(njxbout,vn_pprime,pprim)
      CALL cdf_write(njxbout,vn_aminfor,aminfor)
      CALL cdf_write(njxbout,vn_amaxfor,amaxfor)
      CALL cdf_write(njxbout,vn_jsupu,jsupu3)
      CALL cdf_write(njxbout,vn_jsupv,jsupv3)
      CALL cdf_write(njxbout,vn_jsups,jsups3)
      CALL cdf_write(njxbout,vn_bsupu,bsupu3)
      CALL cdf_write(njxbout,vn_bsupv,bsupv3)
      CALL cdf_write(njxbout,vn_jcrossb,jcrossb)
      CALL cdf_write(njxbout,vn_jxb_gradp,jxb_gradp)
      CALL cdf_write(njxbout,vn_bsubu,bsubu3)
      CALL cdf_write(njxbout,vn_bsubv,bsubv3)
      CALL cdf_write(njxbout,vn_bsubs,bsubs3)
 
      CALL cdf_close(njxbout)

      DEALLOCATE(
     1     bsubs3, bsubv3, bsubu3, jxb_gradp, jcrossb, bsupv3, 
     2     bsupu3, jsups3, jsupv3, jsupu3, jdotb_sqrtg, phin, 
     3     toroidal_angle, sqrtg3, stat=j)

#else
      CLOSE (njxbout)

   90 FORMAT(/"   LU      JSUPU      JSUPV      JSUPS      BSUPU",
     1   "      BSUPV   SQRT(g')     J X B   J X B - p'     J * B",
     2   "      BSUBU      BSUBV      BSUBS   "/)
  100 FORMAT( " TOROIDAL ANGLE (PER PERIOD) = ", f8.3," DEGREES",
     1        " (PLANE #", i3,")")
  110 FORMAT(i5,1p,12e11.3)
  200 FORMAT(/" TOROIDAL FLUX =  ",1p,e12.3,3x,"<J X B - p'> = ",
     1   e12.3,3x,"<J DOT B> = ",e12.3,3x,
     2   "<B DOT GRAD(V)> = ",e12.3,/,
     2   " dp/d(VOL) [p'] = ",e12.3,3x,'d(VOL)/ds    = ',e12.3,3x,
     2   "<JSUPU>   = ",e12.3,3x,"<JSUPV>         = ",e12.3,/,
     3   " MAXIMUM FORCE DEVIATIONS (RELATIVE TO p'): ",sp,0p,f7.2,"%",
     4     3x,f7.2,"%")
#endif
      END IF
      
      DEALLOCATE (jperpu, jperpv, sqgb2, sqrtg, jp2, brhomn, bsubsmn, 
     1    jxb, jxb2, bsupu1, bsupv1, bsubu1, bsubv1, avforce, aminfor, 
     2    amaxfor, pprim, stat=j)
!
!     COMPUTE MERCIER CRITERION
!
      bdotj = mu0*bdotj
      CALL Mercier(gsqrt,bsq,bdotj,iotas,wint,r1,ru,rv,zu,zv,bsubu,
     1             vp,phips,pres,ns,nznt)

      DEALLOCATE (bdotj, bsubuv, bsubvu, stat=j)

      END SUBROUTINE jxbforce
