      SUBROUTINE add_fluxes(overg, bsupu, bsupv, lcurrent)
      USE vmec_main
      USE realspace, ONLY: wint, guu, guv, chip, phip
      USE precon2d, ONLY: ictrl_prec2d
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), DIMENSION(nrzt), INTENT(in)    :: overg
      REAL(rprec), DIMENSION(nrzt), INTENT(inout) :: bsupu, bsupv
      LOGICAL, INTENT(in) :: lcurrent
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(rprec), PARAMETER :: p5=0.5_dp, c1p5=1.5_dp
      REAL(rprec), PARAMETER :: iotaped = 0.10
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: js, l
      REAL(rprec) :: top, bot
!-----------------------------------------------
!
!     ADD MAGNETIC FLUX (CHIP, PHIP) TERMS TO BSUPU=-OVERG*LAM_V, BSUPV=OVERG*LAM_U
!     COMPUTE FLUX FROM ITOR = <B_u>, ITOR(s) = integrated toroidal current (icurv)
!     IF ncurr == 1
!
      IF (.not.lcurrent .or. ncurr.eq.0) GOTO 100

      DO js = 2, ns
         top = icurv(js)
         bot = 0
         DO l = js, nrzt, ns
            top = top - wint(l)*(guu(l)*bsupu(l) + guv(l)*bsupv(l))
            bot = bot + wint(l)*overg(l)*guu(l)
         END DO
         IF (bot .ne. zero) chips(js) = top/bot
         IF (phips(js) .ne. zero) iotas(js) = chips(js)/phips(js)
      END DO

 100  CONTINUE

!     CHANGE THIS FOR lRFP = T  (solve for phips?)
      IF (.not.lcurrent) chips(2:ns) = iotas(2:ns)*phips(2:ns)

      DO js = 2, ns
         chip(js:nrzt:ns) = chips(js)
      END DO

      chipf(2:ns1) = (chips(2:ns1) + chips(3:ns1+1))/2
      chipf(ns)    = 2*chips(ns)-chips(ns1)


!     Do not compute iota too near origin
      IF (lrfp) THEN
         iotaf(1)  = one/(c1p5/iotas(2) - p5/iotas(3))
         iotaf(ns) = one/(c1p5/iotas(ns) - p5/iotas(ns-1))
         DO js = 2, ns-1
            iotaf(js) = 2.0_dp/(one/iotas(js) + one/iotas(js+1))
         END DO

      ELSE
         iotaf(1)  = c1p5*iotas(2) - p5*iotas(3)               !zero gradient near axis
         iotaf(ns) = c1p5*iotas(ns) - p5*iotas(ns-1)
         DO js = 2, ns-1
            iotaf(js) = p5*(iotas(js) + iotas(js+1))
         END DO
      END IF

      bsupu(:nrzt) = bsupu(:nrzt)+chip(:nrzt)*overg(:nrzt)

!     FOR PRECONDITIONER, GIVE IOTA A PEDESTAL VALUE TO AVOID SINGULAR HESSIAN
!     WITHOUT THIS, THE D3D curhole CASES, WHERE iota CHANGES SIGN, DO NOT CONVERGE WELL IF AT ALL
!      IF (ictrl_prec2d .ge. 2) THEN
!         DO js = 2, ns
!           IF (ABS(iotas(js)) .le. iotaped) THEN
!              bsupu(js:nrzt:ns) = bsupu(js:nrzt:ns) +                   &
!     &         (SIGN(iotaped,chips(js)) - chips(js))*overg(js:nrzt:ns)
!           END IF
!         END DO
!      END IF

      END SUBROUTINE add_fluxes
