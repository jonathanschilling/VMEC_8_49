      MODULE spline_parm
      USE stel_kinds
      IMPLICIT INTEGER*4 (i-n), REAL(rprec) (a-h, o-z)
      INTEGER, PARAMETER :: nw=257, nh=257, nwnh=nw*nh
      INTEGER, PARAMETER :: nh2=2*nh, nwrk=2*(nw+1)*nh
      INTEGER, PARAMETER :: kubicx = 4, kubicy = 4, 
     1           lubicx = nw - kubicx + 1,
     2           lubicy = nh - kubicy + 1,
     3           kujunk = kubicx*kubicy*lubicx*lubicy
      END MODULE spline_parm
