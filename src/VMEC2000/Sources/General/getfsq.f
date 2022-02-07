      SUBROUTINE getfsq(gcr, gcz, gnormr, gnormz, gnorm, medge)
      USE vmec_main
      USE vmec_params, ONLY: ntmax
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: medge
      REAL(rprec), INTENT(out) :: gnormr, gnormz
      REAL(rprec), INTENT(in)  :: gnorm
      REAL(rprec), DIMENSION(ns,mnsize*ntmax), INTENT(in) :: gcr, gcz
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: jsmax
C-----------------------------------------------
      jsmax = ns1 + medge
      gnormr = gnorm * SUM(gcr(:jsmax,:)**2)
      gnormz = gnorm * SUM(gcz(:jsmax,:)**2)

      END SUBROUTINE getfsq
