      MODULE fbal
      USE stel_kinds, ONLY: dp
      REAL(dp), DIMENSION(:), ALLOCATABLE :: rzu_fac, rru_fac,
     1  frcc_fac, fzsc_fac

      CONTAINS

      SUBROUTINE calc_fbal(bsubu, bsubv)
      USE vmec_main, ONLY: buco, bvco, equif, 
     1                     jcurv, jcuru, chipf, vp, pres, 
     2                     phipf, vpphi, presgrad, ohs
      USE vmec_params, ONLY: signgs
      USE vmec_dim, ONLY: ns, nrzt
      USE vmec_input, ONLY: lrfp
      USE realspace, ONLY: wint, phip
C-----------------------------------------------
      REAL(dp), INTENT(in) :: bsubu(1:nrzt), bsubv(1:nrzt)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER  :: js, ns1
C-----------------------------------------------
      DO js = 2, ns
         buco(js) = SUM(bsubu(js:nrzt:ns)*wint(js:nrzt:ns))
         bvco(js) = SUM(bsubv(js:nrzt:ns)*wint(js:nrzt:ns))
      END DO

!     FROM AMPERE'S LAW, JcurX are angle averages of jac*JsupX, so
!                        JcurX = (dV/ds)/twopi**2 <JsupX> where <...> is flux surface average
      ns1 = ns-1
      DO js = 2, ns1
         jcurv(js) = (signgs*ohs)*(buco(js+1) - buco(js))
         jcuru(js) =-(signgs*ohs)*(bvco(js+1) - bvco(js))
!FOR RFP         vpphi(js) = (vp(js+1)/phip(js+1) + vp(js)/phip(js))/2
         vpphi(js) = (vp(js+1) + vp(js))/2
         presgrad(js) = (pres(js+1) - pres(js))*ohs
         equif(js) = (-phipf(js)*jcuru(js) + chipf(js)*jcurv(js))
     1                /vpphi(js) + presgrad(js)
      END DO

      equif(1) = 0
      equif(ns) = 0

!      IF ((fsqr+fsqz+fsql).gt.1.E-08_dp .and. ictrl_prec2d.eq.0) RETURN
!      IF (lRFP) RETURN
     
!      equif(2:ns1) = (-jcuru(2:ns1) + iotaf(2:ns1)*jcurv(2:ns1))
!     1             /vpphi(2:ns1)  + presgrad(2:ns1)

      END SUBROUTINE calc_fbal

      END MODULE fbal
