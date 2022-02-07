      SUBROUTINE seva2d(bkx,lx,bky,ly,cs,nx,ny,xl,yl,fs,ier,icalc)
c------------------------------------------------------------------------------
c--  S.Thompson  92/05/18                                                    --
c--    Bicubic spline routines.                                              --
c--    Put together with routines from E.Solano.                             --
c--  SMWolfe     93/12/17                                                    --
c--    Modifed to avoid over-writing the original array.                     --
c--              94/04/28                                                    --
c--    Updated.                                                              --
c------------------------------------------------------------------------------
c  Inputs:
c
c      cs       - array of spline coefficients of dimension (kubicx,
c                 lubicx,kubicy,lubicy) from sets2d.
c
c      bkx, bky - interval coefficients of length lubicx+1 and lubicy+1 from
c                 sets2d.
c
c      lx, ly   - number of terms in bkx and bky from sets2d.
c
c      xl, yl   - the point at which interpolations are desired.
c
c      nx, ny   - grid dimensions
c
c  Outputs:
c
c      fs       - vector containing results depending on icalc:
c                 icalc              fs
c                   1                f
c                   2                fx
c                   3                fy
c                   4                fxy
c                   5                fxx
c                   6                fyy
c
c      ier      - error parameter.
c
c-------------------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE spline_parm
      IMPLICIT NONE
c
      INTEGER :: ier, lx, ly, nx, ny, icalc
      INTEGER :: lef, ibk, jj, mflag, ndummy
c      real*8 cs(kubicx,lubicx,kubicy,lubicy),xl,yl,fs(6),bkx(1),bky(1)
      REAL(rprec) :: xl,yl,fs(6),bkx(lx),bky(ly)
      REAL(rprec) :: cs(kubicx,nx-kubicx+1,kubicy,ny-kubicy+1)
      REAL(rprec) :: ppvalw
c
c  Local Variable Specifications:
c
      REAL(rprec) :: work0(4),work1(4),work2(4), h
      INTEGER, PARAMETER :: n00 = 0, n11 = 1, n22 = 2
c
c  Evaluate function and its partial derivatives at (XL, YL):
c
c
c  First do all the lookup and interpolation stuff.
c  This is the most time consuming part of the evaluation, so
c  don't do more than needed.
c
      CALL interv(bky,ly,yl,lef,mflag)
      CALL interv(bkx,lx,xl,ibk,ndummy)
      h = xl - bkx(ibk)
      DO 41 jj=1,4
         work0(jj) = ppvalw(cs(1,ibk,jj,lef),h,n00)
         IF (icalc.eq.1) GOTO41
         work1(jj) = ppvalw(cs(1,ibk,jj,lef),h,n11)
         IF (icalc.le.4) GOTO41
         work2(jj) = ppvalw(cs(1,ibk,jj,lef),h,n22)
 41   CONTINUE
      h = yl - bky(lef)
      fs(1) = ppvalw(work0,h,n00)
      IF (icalc.eq.1) RETURN
      fs(2) = ppvalw(work1,h,n00)
      IF (icalc.eq.2) RETURN
      fs(3) = ppvalw(work0,h,n11)
      IF (icalc.eq.3) RETURN
      fs(4) = ppvalw(work1,h,n11)
      IF (icalc.eq.4) RETURN
      fs(5) = ppvalw(work2,h,n00)
      IF (icalc.eq.5) RETURN
      fs(6) = ppvalw(work0,h,n22)
C
      RETURN
      END

      SUBROUTINE sets2d(s,cs,x,nx,bkx,lx,y,ny,bky,ly,wk,ier)
c------------------------------------------------------------------------------
c--  S.Thompson  92/05/18                                                    --
c--    Bicubic spline routines.                                              --
c--    Put together with routines from E.Solano.                             --
c--  SMWolfe     93/12/17                                                    --
c--    Modifed to avoid over-writing the original array.                     --
c--              94/04/28                                                    --
c--    Updated.                                                              --
c------------------------------------------------------------------------------
c  Inputs:
c
c      s     - nx by ny array containing the function values at (x,y).
c              This is a 1-d array, k=k=(i-1)*ny+j.
c              Modified ordering k = (j-1)*nx+i
c
c      x, y  - (x,y) location, arrays of length nx and ny.
c
c  Outputs:
c
c      cs    - array of spline coefficients of dimension (kubicx,
c              lubicx,kubicy,lubicy).
c
c      bkx, bky - interval coefficients of length lubicx+1 and lubicy+1.
c
c      lx, ly -   number of terms in bkx and bky.
c
c      ier   - rror parameter.
c
c  Work arrays:
c
c      wk    - of dimension at least nx by ny.
c------------------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      USE spline_parm
c
c      DIMENSION s(1), x(nx), y(ny), wk(nx,ny),
c     .          xknot(kubicx + nw), yknot(kubicy + nh),
c     .          cs(kubicx, lubicx, kubicy, lubicy),
c     .          bkx(lubicx + 1), bky(lubicy + 1)
      INTEGER :: i, j, k, lx, nx, ny, ly, ier
      REAL(rprec) :: s(nx*ny), x(nx), y(ny), wk(nx,ny),
     .          xknot(kubicx + nx), yknot(kubicy + ny),
     .          cs(kubicx, nx-kubicx+1, kubicy, ny-kubicy+1),
     .          bkx(nx-kubicx + 2), bky(ny-kubicy + 2)
c
c  Set up knots:
c
c     WRITE (6,*) x(38),y(42),s(41*nx+38)*x(38)
      CALL eknot (nx, x, kubicx, xknot)		
      CALL eknot (ny, y, kubicy, yknot)			
c
c  Save the original, use the work array
c
      DO 10 I=1,NX
      DO 10 j=1,ny
c        k=(i-1)*ny+j
         k=(j-1)*nx+i
  10     wk(i,j) = s(k)
c
c  Calculate spline coefficients:
c
      CALL spl2bc (x, y, nx, ny,xknot, yknot, wk)	
c
c  Coefficients stored in bkx, bky, and c:
c
      CALL spl2pp (nx, ny, xknot, yknot, wk, bkx, lx, bky, ly, cs)
c
      RETURN
      END

      SUBROUTINE spl2bc(rgrid,zgrid,mw,mh,rknot,zknot,copy)
calculates the b-spline coeficients
      USE stel_kinds, ONLY: rprec
      IMPLICIT INTEGER*4 (i-n), REAL(rprec) (a-h, o-z)
      PARAMETER (nw=257,nh=257,krord=4,kzord=4)
      REAL(rprec) :: rgrid(mw),zgrid(mh)
c      DIMENSION rknot(nw+krord),zknot(nh+kzord),copy(mw,mh)
      REAL(rprec) :: rknot(mw+krord),zknot(mh+kzord),copy(mw,mh)
c------------------------------------------------------------------
c-- change DIMENSION of work2 and work3 from nw to nh            --
c-- to ensure the cases when nh > nw     ll, 93/04/01            --
c------------------------------------------------------------------
      DIMENSION work1(mw,mh),work2(mh),work3(mh,2*krord-1)
      CALL spli2d(rgrid,copy,rknot,mw,krord,mh,work2,work3,work1,iflag)
      IF (iflag.ne.1) PRINT *,' error in first spli2d, iflag=',iflag
      CALL spli2d(zgrid,work1,zknot,mh,kzord,mw,work2,work3,copy,iflag)
      IF (iflag.ne.1) PRINT *,' error in second spli2d, iflag=',iflag
      RETURN
      END

      SUBROUTINE spl2pp(mw,mh,rknot,zknot,copy,breakr,lr,breakz,lz,coef)
c translates to pp representation
      USE stel_kinds, ONLY: rprec
      IMPLICIT INTEGER*4 (i-n), REAL(rprec) (a-h, o-z)
      PARAMETER (nw=257,nh=257,krord=4,kzord=4)
c      PARAMETER  (lr0=nw-krord+1,lz0=nh-kzord+1)
c      DIMENSION rknot(nw+krord),zknot(nh+kzord)
c      DIMENSION copy(mw,mh),coef(krord,lr0,kzord,lz0)
c      DIMENSION breakr(lr0+1),breakz(lz0+1)
c      DIMENSION work4(krord,nw,nh), work5(nh,krord,lr0)
c     *         ,work6(kzord,kzord,nw,krord)
      DIMENSION rknot(mw+krord),zknot(mh+kzord)
      DIMENSION copy(mw,mh),coef(krord,mw-krord+1,kzord,mh-kzord+1)
      DIMENSION breakr(mw-krord+2),breakz(mh-kzord+2)
      DIMENSION work4(krord,nw,nh), work5(mh,krord,mw-krord+1)
     *         ,work6(kzord,kzord,nw,krord)
      EQUIVALENCE (work4,work6)
      CALL bspp2d(rknot,copy,mw,krord,mh,work4,breakr,work5,lr)
      ndum=lr*krord
      CALL bspp2d(zknot,work5,mh,kzord,ndum    ,work6,breakz,coef,lz)
      RETURN
      END

      SUBROUTINE eknot(n,x,k,xk)
c given the ordered data points x(1)<...<x(n), this SUBROUTINE generates
c a knot sequence with not-a-knot end conditions (like BSNAK from IMSL)
c Some of this is discussed in de Boor(1978), page 211.
      USE stel_kinds, ONLY: rprec
      IMPLICIT INTEGER*4 (i-n), REAL(rprec) (a-h, o-z)
        DIMENSION x(n),xk(n+k)
        INTEGER kh
c
        DO i=1,k
        xk(i)=x(1)
        ii=i+n
        xk(ii)= x(n)+1.e-5
        END DO
        kh=k/2
        k2=kh+kh
        IF (k2.eq.k) THEN
c even k, place knots at data points
        DO i=k+1,n
        xk(i)=x(i-kh)
        END DO
        ELSE
c odd k, place knots in between data points
        DO i=k+1,n
        xk(i)=.5*(x(i-kh)+x(i-1-kh))
        END DO
        END IF
        RETURN
        END

      SUBROUTINE spli2d ( tau, gtau, t, n, k, m, work, q, bcoef, iflag )
calls bsplvb, banfac/slv
c  this is an extended version of  splint , for the use in tensor prod-
c  uct interpolation.
c
c   spli2d  produces the b-spline coeff.s  bcoef(j,.)  of the spline of
c   order  k  with knots  t (i), i=1,..., n + k , which takes on the
c   value  gtau (i,j)  at  tau (i), i=1,..., n ; j=1,..., m .
c
c******  i n p u t  ******
c  tau   array of length  n , containing data point abscissae.
c  a s s u m p t i o n . . .  tau  is strictly increasing
c  gtau(.,j)  corresponding array of length  n , containing data point
c        ordinates, j=1,...,m
c  t     knot sequence, of length  n+k
c  n     number of data points and dimension of spline space  s(k,t)
c  k     order of spline
c  m     number of data sets
c
c******  w o r k   a r e a  ******
c  work  a vector of length  n
c
c******  o u t p u t  ******
c  q     array of order  (n,2*k-1), containing the triangular factoriz-
c        ation of the coefficient matrix of the linear system for the b-
c        coefficients of the spline interpolant.
c           the b-coeffs for the interpolant of an additional data set
c        (tau(i),htau(i)), i=1,...,n  with the same data abscissae can
c        be obtained without going through all the calculations in this
c        routine, simply by loading  htau  into  bcoef  and then execut-
c        ing the    call banslv ( q, n, n, 2*k-1, k, bcoef )
c  bcoef the b-coefficients of the interpolant, of length  n
c  iflag an integer indicating success (= 1)  or failure (= 2)
c        the linear system to be solved is (theoretically) invertible if
c        and only if
c              t(i) .lt. tau(i) .lt. tau(i+k),    all i.
c        violation of this condition is certain to lead to  iflag = 2 .
c
c******  m e t h o d  ******
c     the i-th equation of the linear system  a*bcoef = b  for the b-co-
c  effs of the interpolant enforces interpolation at  tau(i), i=1,...,n.
c  hence,  b(i) = gtau(i), all i, and  a  is a band matrix with  2k-1
c  bands (if it is invertible).
c     the matrix  a  is generated row by row and stored, diagonal by di-
c  agonal, in the  c o l u m n s  of the array  q , with the main diag-
c  onal going into column  k .  see comments in the program below.
c     the banded system is then solved by a call to  banfac (which con-
c  structs the triangular factorization for  a  and stores it again in
c   q ), followed by a call to  banslv (which then obtains the solution
c   bcoef  by substitution).
c     banfac  does no pivoting, since the total positivity of the matrix
c  a  makes this unnecessary.
c
c      integer iflag,k,m,n,i,ilp1mx,j,jj,kpkm1,left,np1
c      real*8 bcoef(m,n),gtau(n,m),q(n,7),t(n+k),tau(n),work(n),taui
      USE stel_kinds, ONLY: rprec
      IMPLICIT INTEGER*4 (i-n), REAL(rprec) (a-h, o-z)
      DIMENSION bcoef(m,n),gtau(n,m),q(n,2*k-1),t(n+k),tau(n),work(n)
c
      nnn=1
      np1 = n + 1
      kpkm1 = 2*k - 1
      left = k
c
c  ***   loop over  i  to construct the  n  interpolation equations
      DO 30 i=1,n
         iindex=i
         taui = tau(iindex)
         ilp1mx = MIN(iindex+k,np1)
c        *** zero out all entries in row  i  of  a (in the 2k-1 bands)
         DO 13 j=1,kpkm1
   13       q(iindex,j) = 0.
c        *** find  left  in the closed interval (i,i+k-1) such that
c                t(left) .le. tau(i) .lt. t(left+1)
c        matrix is singular if this is not possible
         left = MAX(left,i)
         IF (taui .lt. t(left))         GOTO 998
   15       IF (taui .lt. t(left+1))    GOTO 16
            left = left + 1
            IF (left .lt. ilp1mx)       GOTO 15
         left = left - 1
         IF (taui .gt. t(left+1))       GOTO 998
c        *** the i-th equation enforces interpolation at taui, hence
c        a(i,j) = b(j,k,t)(taui), all j. only the  k  entries with  j =
c        left-k+1,...,left actually might be nonzero. these  k  numbers
c        are returned, in  work  (used for temp.storage here), by the
c        following
   16    CALL bsplvb ( t, k, nnn, taui, left, work )
c        we therefore want  work(j) = b(left-k+j)(taui) to go into
c        a(i,left-k+j), i.e., into  q(i,left-i+j), since the i-th row of
c        a  is so stored in the i-th row of  q  that the (i,i)-entry of
c        a  goes into the  k-th  entry of  q.
         jj = left - iindex
         DO 29 j=1,k
            jj = jj+1
            q(iindex,jj) = work(j)
   29    CONTINUE
   30    CONTINUE
c
c     ***obtain factorization of  a  , stored again in  q.
      CALL banfac ( q, n, n, kpkm1, k, iflag )
                                        GOTO (40,999), iflag
c     *** solve  a*bcoef = gtau  by backsubstitution
   40 DO 50 j=1,m
         DO 41 i=1,n
   41       work(i) = gtau(i,j)
         CALL banslv ( q, n, n, kpkm1, k, work )
         DO 50 i=1,n
   50    bcoef(j,i) = work(i)
                                        RETURN
  998 iflag = 2
  999 PRINT 699
  699 FORMAT(41h linear system in  splint  not invertible)
                                        RETURN
      END

      SUBROUTINE bspp2d ( t, bcoef, n, k, m, scrtch, break, coef, l )
calls  bsplvb
c  this is an extended version of  bsplpp  for use with tensor products
c
converts the b-representation  t, bcoef(.,j), n, k  of some spline into
c  its pp-representation  break, coef(j,.,.), l, k ; j=1, ..., m  .
c
c******  i n p u t  ******
c  t     knot sequence, of length  n+k
c  bcoef(.,j) b-spline coefficient sequence, of length  n ;j=1,...,m
c  n     length of  bcoef  and  dimension of spline space  s(k,t)
c  k     order of the spline
c  m     number of data sets
c
c******  w o r k   a r e a  ******
c  scrtch   of size  (k,k,m), needed to contain bcoeffs of a piece of
c        the spline and its  k-1  derivatives   for each of the m sets
c
c******  o u t p u t  ******
c  break breakpoint sequence, of length  l+1, contains (in increasing
c        order) the distinct points in the sequence  t(k), ..., t(n+1)
c  coef(mm,.,.)  array of size (k,n), with  coef(mm,i,j) = (i-1)st der-
c        ivative of  mm-th  spline at break(j) from the right, mm=1,.,m
c  l     number of polynomial pieces which make up the spline in the
c        interval  (t(k), t(n+1))
c
c******  m e t h o d  ******
c     for each breakpoint interval, the  k  relevant b-coeffs of the
c  spline are found and then differenced repeatedly to get the b-coeffs
c  of all the derivatives of the spline on that interval. the spline and
c  its first  k-1  derivatives are then evaluated at the left end
c  point of that interval, using  bsplvb  repeatedly to obtain the val-
c  ues of all b-splines of the appropriate order at that point.
c
      USE stel_kinds, ONLY: rprec
      IMPLICIT INTEGER*4 (i-n), REAL(rprec) (a-h, o-z)
      PARAMETER (kmax=4)
      INTEGER k,l,m,n,   i,j,jp1,kmj,left
      REAL(rprec) :: bcoef(n,m), break(*), coef(m,k,*),scrtch(k,k,m),
     *     t(n+k), biatx(kmax)
      REAL(rprec) :: diff,fkmj,sum
c
      n11=1
      n22=2
      l = 0
      break(1) = t(k)
      DO 50 left=k,n
c        find the next nontrivial knot interval.
         IF (t(left+1) .eq. t(left))    GOTO 50
         l = l + 1
         break(l+1) = t(left+1)
         IF (k .gt. 1)                  GOTO 9
         DO 5 mm=1,m
    5       coef(mm,1,l) = bcoef(left,mm)
                                        GOTO 50
c        store the k b-spline coeff.s relevant to current knot interval
c        in  scrtch(.,1) .
    9    DO 10 i=1,k
            DO 10 mm=1,m
   10          scrtch(i,1,mm) = bcoef(left-k+i,mm)
c        for j=1,...,k-1, compute the k-j b-spline coeff.s relevant to
c        current knot interval for the j-th derivative by differencing
c        those for the (j-1)st derivative, and store in scrtch(.,j+1) .
         DO 20 jp1=2,k
            j = jp1 - 1
            kmj = k - j
            fkmj = REAL(kmj, rprec)
            DO 20 i=1,kmj
               diff = (t(left+i) - t(left+i - kmj))/fkmj
               IF (diff .le. 0.)         GOTO 20
               DO 15 mm=1,m
   15             scrtch(i,jp1,mm) =
     *            (scrtch(i+1,j,mm) - scrtch(i,j,mm))/diff
   20          CONTINUE
c        starting with the one b-spline of order 1 not zero at t(left),
c        find the values at t(left) of the j+1 b-splines of order j+1
c        not identically zero there from those of order j, then combine
c        with the b-spline coeff.s found earlier to compute the (k-j)-
c        th derivative at t(left) of the given spline.
         CALL bsplvb ( t, n11, n11, t(left), left, biatx )
         DO 25 mm=1,m
   25       coef(mm,k,l) = scrtch(1  ,k,mm)
         DO 30 jp1=2,k
            CALL bsplvb ( t, jp1, n22, t(left), left, biatx )
            kmj = k+1 - jp1
            DO 30 mm=1,m
               sum = 0.
               DO 28 i=1,jp1
   28             sum = biatx(i)*scrtch(i,kmj,mm) + sum
   30          coef(mm,kmj,l) = sum
   50    CONTINUE
         RETURN
      END

      SUBROUTINE bsplvb ( t, jhigh, index, x, left, biatx )
calculates the value of all possibly nonzero b-splines at  x  of order
c
c               jout  =  max( jhigh , (j+1)*(index-1) )
c
c  with knot sequence  t .
c
c******  i n p u t  ******
c  t.....knot sequence, of length  left + jout  , assumed to be nonde-
c        creasing.  a s s u m p t i o n . . . .
c                       t(left)  .lt.  t(left + 1)   .
c   d i v i s i o n  b y  z e r o  will result if  t(left) = t(left+1)
c  jhigh,
c  index.....integers which determine the order  jout = max(jhigh,
c        (j+1)*(index-1))  of the b-splines whose values at  x  are to
c        be returned.  index  is used to avoid recalculations when seve-
c        ral columns of the triangular array of b-spline values are nee-
c        ded (e.g., in  bvalue  or in  bsplvd ). precisely,
c                     if  index = 1 ,
c        the calculation starts from scratch and the entire triangular
c        array of b-spline values of orders 1,2,...,jhigh  is generated
c        order by order , i.e., column by column .
c                     if  index = 2 ,
c        only the b-spline values of order  j+1, j+2, ..., jout  are ge-
c        nerated, the assumption being that  biatx , j , deltal , deltar
c        are, on entry, as they were on exit at the previous call.
c           in particular, if  jhigh = 0, then  jout = j+1, i.e., just
c        the next column of b-spline values is generated.
c
c  w a r n i n g . . .  the restriction   jout .le. jmax (= 20)  is im-
c        posed arbitrarily by the dimension statement for  deltal  and
c        deltar  below, but is  n o w h e r e  c h e c k e d  for .
c
c  x.....the point at which the b-splines are to be evaluated.
c  left.....an integer chosen (usually) so that
c                  t(left) .le. x .le. t(left+1)  .
c
c******  o u t p u t  ******
c  biatx.....array of length  jout , with  biatx(i)  containing the val-
c        ue at  x  of the polynomial of order  jout  which agrees with
c        the b-spline  b(left-jout+i,jout,t)  on the interval (t(left),
c        t(left+1)) .
c
c******  m e t h o d  ******
c  the recurrence relation
c
c                       x - t(i)              t(i+j+1) - x
c     b(i,j+1)(x)  =  -----------b(i,j)(x) + ---------------b(i+1,j)(x)
c                     t(i+j)-t(i)            t(i+j+1)-t(i+1)
c
c  is used (repeatedly) to generate the (j+1)-vector  b(left-j,j+1)(x),
c  ...,b(left,j+1)(x)  from the j-vector  b(left-j+1,j)(x),...,
c  b(left,j)(x), storing the new values in  biatx  over the old. the
c  facts that
c            b(i,1) = 1  if  t(i) .le. x .lt. t(i+1)
c  and that
c            b(i,j)(x) = 0  unless  t(i) .le. x .lt. t(i+j)
c  are used. the particular organization of the calculations follows al-
c  gorithm  (8)  in chapter x of the text.
c
      USE stel_kinds, ONLY: rprec
      IMPLICIT INTEGER*4 (i-n), REAL(rprec) (a-h, o-z)
      PARAMETER(jmax = 4)
      INTEGER :: index, jhigh, left, i, j=1, jp1
      REAL(rprec) :: x,saved,term
c      real*8 biatx(jhigh),t(1),x,
      DIMENSION deltal(jmax),deltar(jmax)
      DIMENSION biatx(jhigh), t(left+jhigh)
current fortran standard makes it impossible to specify the length of
c  t  and of  biatx  precisely without the introduction of otherwise
c  superfluous additional arguments.
      SAVE deltal,deltar  ! (valid in fortran 77)
c
                                        GOTO (10,20), index
   10 j = 1
      biatx(1) = 1.
      IF (j .ge. jhigh)                 GOTO 99
c
   20    jp1 = j + 1
         deltar(j) = t(left+j) - x
         deltal(j) = x - t(left+1-j)
         saved = 0.
         DO 26 i=1,j
            term = biatx(i)/(deltar(i) + deltal(jp1-i))
            biatx(i) = saved + deltar(i)*term
   26       saved = deltal(jp1-i)*term
         biatx(jp1) = saved
         j = jp1
         IF (j .lt. jhigh)              GOTO 20
c
   99                                   RETURN
      END

      FUNCTION ppvalw (coef, x, jd )
C-----------------------------------------------------------------------        
C  Modified for optimization by S.J. Thompson, 30-Aug-1993
c  Revised to eliminate call to interv by S.M.Wolfe, 17-Dec-1993
c          and to use ASF's for evaluation
c  This routine performs only the innermost guts of the spline evaluation
c  Assumes k=4 (cubic spline only). No other cases considered.
c does not call  interv
calculates value at  x  of  jd-th derivative of pp fct from pp-repr
c
c******  i n p u t  ****** to PPVALU, on which this is based.
c  break, coef, l, k.....forms the pp-representation of the function  f
c        to be evaluated. specifically, the j-th derivative of  f  is
c        given by
c
c     (d**j)f(x) = coef(j+1,i) + h*(coef(j+2,i) + h*( ... (coef(k-1,i) +
c                             + h*coef(k,i)/(k-j-1))/(k-j-2) ... )/2)/1
c
c        with  h = x - break(i),  and
c
c       i  =  max( 1 , max( j ;  break(j) .le. x , 1 .le. j .le. l ) ).
c
c  x.....the point at which to evaluate.
c        as used here, x is the distance from the break, not the absolute 
c        position. 
c  jd.....integer*4 giving the order of the derivative to be evaluat-
c        ed.  a s s u m e d  to be zero or positive.
c
c******  o u t p u t  ******
c  ppvalw.....the value of the (jd)-th derivative of  f  at  x.
c
c******  m e t h o d  ******
c     the interval index  i , appropriate for  x , is found through a
c  call to  interv . the formula above for the  jd-th derivative
c  of  f  is then evaluated (by nested multipication).
c
C-----------------------------------------------------------------------        
C   Variable declarations.
C-----------------------------------------------------------------------        
      USE stel_kinds, ONLY: rprec
      IMPLICIT INTEGER*4 (i-n), REAL(rprec) (a-h, o-z)
      INTEGER, INTENT(in) :: jd
      REAL(rprec) :: ppvalw,x
      REAL(rprec) :: coef(4)
c----------------------------------------------------------------------
c ASF's may be slightly more efficient than the alternative
c----------------------------------------------------------------------
      d2(xx) = coef(4)*xx + coef(3)
      d1(xx) = (coef(4)*xx/2 + coef(3))*xx + coef(2)
      d0(xx) = ((coef(4)*xx/3 + coef(3))*xx/2 + 
     .           coef(2))*xx + coef(1)
C-----------------------------------------------------------------------        
C   Derivatives of order k or higher are identically zero.
C-----------------------------------------------------------------------        
C   Evaluate jd-th derivative of i-th polynomial piece at x .
C-----------------------------------------------------------------------        
      GOTO(1,2,3) jd+1
      ppvalw = 0.
      PRINT *, 'Error (ppvalw): JD must be 0, 1, or 2.'
      PRINT *, 'Execution terminated.'
      RETURN
 1    ppvalw = d0(x)	! k = 4 , jd = 0
      RETURN
 2    ppvalw = d1(x)	! k = 4 , jd = 1
      RETURN
 3    ppvalw = d2(x)	! k = 4 , jd = 2
      RETURN
      END
C
      SUBROUTINE banslv ( a, nrow, n, ndiag, middle, b )
      USE stel_kinds, ONLY: rprec
      IMPLICIT INTEGER*4 (i-n), REAL(rprec) (a-h, o-z)
      DIMENSION a(nrow,ndiag),b(n)
      IF (n .eq. 1)                     GOTO 21
      ilo = middle - 1
      IF (ilo .lt. 1)                   GOTO 21
      DO 19 i=2,n
         jmax = MIN(i-1,ilo)
         DO 19 j=1,jmax
   19       b(i) = b(i) - b(i-j)*a(i,middle-j)
c
   21 ihi = ndiag-middle
      DO 30 i=n,1,-1
         jmax = MIN(n-i,ihi)
         IF (jmax .lt. 1)               GOTO 30
         DO 25 j=1,jmax
   25       b(i) = b(i) - b(i+j)*a(i,middle+j)
   30    b(i) = b(i)/a(i,middle)
                                        RETURN
      END

      SUBROUTINE banfac ( a, nrow, n, ndiag, middle, iflag )
      USE stel_kinds, ONLY: rprec
      IMPLICIT INTEGER*4 (i-n), REAL(rprec) (a-h, o-z)
      DIMENSION a(nrow,ndiag)
      iflag = 1
      ilo = middle - 1
      IF (ilo)                          999,10,19
   10 DO 11 i=1,n
         IF(a(i,1) .eq. 0.)             GOTO 999
   11    CONTINUE
                                        RETURN
   19 ihi = ndiag - middle
      IF (ihi)                          999,20,29
   20 DO 25 i=1,n
         IF (a(i,middle) .eq. 0.)       GOTO 999
         jmax = MIN(ilo,n-i)
         IF (jmax .lt. 1)               GOTO 25
         DO 23 j=1,jmax
   23       a(i+j,middle-j) = a(i+j,middle-j)/a(i,middle)
   25    CONTINUE
                                        RETURN
   29 DO 50 i=1,n
         diag = a(i,middle)
         IF (diag .eq. 0.)              GOTO 999
         jmax = MIN(ilo,n-i)
         IF(jmax .lt. 1)                GOTO 50
         kmax = MIN(ihi,n-i)
         DO 33 j=1,jmax
            mmj = middle-j
            a(i+j,mmj) = a(i+j,mmj)/diag
            DO 33 k=1,kmax
   33          a(i+j,mmj+k) = a(i+j,mmj+k) - a(i+j,mmj)*a(i,middle+k)
   50    CONTINUE
                                        RETURN
  999 iflag = 2
                                        RETURN
      END

      SUBROUTINE interv ( xt, lxt, x, left, mflag )
computes  left = max( i ; 1 .le. i .le. lxt  .and.  xt(i) .le. x )  .
c
c******  i n p u t  ******
c  xt.....a real*8 sequence, of length  lxt , assumed to be nondecreasing
c  lxt.....number of terms in the sequence  xt .
c  x.....the point whose location with respect to the sequence  xt  is
c        to be determined.
c
c******  o u t p u t  ******
c  left, mflag.....both integers, whose value is
c
c   1     -1      if               x .lt.  xt(1)
c   i      0      if   xt(i)  .le. x .lt. xt(i+1)
c  lxt     1      if  xt(lxt) .le. x
c
c        in particular,  mflag = 0 is the 'usual' case.  mflag .ne. 0
c        indicates that  x  lies outside the halfopen interval
c        xt(1) .le. y .lt. xt(lxt) . the asymmetric treatment of the
c        interval is due to the decision to make all pp functions cont-
c        inuous from the right.
c
c******  m e t h o d  ******
c  the program is designed to be efficient in the common situation that
c  it is called repeatedly, with  x  taken from an increasing or decrea-
c  sing sequence. this will happen, e.g., when a pp function is to be
c  graphed. the first guess for  left  is therefore taken to be the val-
c  ue returned at the previous call and stored in the  l o c a l  varia-
c  ble  ilo . a first check ascertains that  ilo .lt. lxt (this is nec-
c  essary since the present call may have nothing to do with the previ-
c  ous call). then, if  xt(ilo) .le. x .lt. xt(ilo+1), we set  left =
c  ilo  and are done after just three comparisons.
c     otherwise, we repeatedly double the difference  istep = ihi - ilo
c  while also moving  ilo  and  ihi  in the direction of  x , until
c                      xt(ilo) .le. x .lt. xt(ihi) ,
c  after which we use bisection to get, in addition, ilo+1 = ihi .
c  left = ilo  is then returned.
c
      USE stel_kinds, ONLY: rprec
      IMPLICIT INTEGER*4 (i-n), REAL(rprec) (a-h, o-z)
      INTEGER left,lxt,mflag,   ihi,ilo,istep,middle
      REAL(rprec) :: x
      DIMENSION xt(lxt)
      DATA ilo /1/
c     save ilo  (a valid fortran statement in the new 1977 standard)
      ihi = ilo + 1
      IF (ihi .lt. lxt)                 GOTO 20
         IF (x .ge. xt(lxt))            GOTO 110
         IF (lxt .le. 1)                GOTO 90
         ilo = lxt - 1
         ihi = lxt
c
   20 IF (x .ge. xt(ihi))               GOTO 40
      IF (x .ge. xt(ilo))               GOTO 100
c
c              **** now x .lt. xt(ilo) . decrease  ilo  to capture  x .
   30 istep = 1
   31    ihi = ilo
         ilo = ihi - istep
         IF (ilo .le. 1)                GOTO 35
         IF (x .ge. xt(ilo))            GOTO 50
         istep = istep*2
                                        GOTO 31
   35 ilo = 1
      IF (x .lt. xt(1))                 GOTO 90
                                        GOTO 50
c              **** now x .ge. xt(ihi) . increase  ihi  to capture  x .
   40 istep = 1
   41    ilo = ihi
         ihi = ilo + istep
         IF (ihi .ge. lxt)              GOTO 45
         IF (x .lt. xt(ihi))            GOTO 50
         istep = istep*2
                                        GOTO 41
   45 IF (x .ge. xt(lxt))               GOTO 110
      ihi = lxt
c
c           **** now xt(ilo) .le. x .lt. xt(ihi) . narrow the interval.
   50 middle = (ilo + ihi)/2
      IF (middle .eq. ilo)              GOTO 100
c     note. it is assumed that middle = ilo in case ihi = ilo+1 .
      IF (x .lt. xt(middle))            GOTO 53
         ilo = middle
                                        GOTO 50
   53    ihi = middle
                                        GOTO 50
c**** set output and return.
   90 mflag = -1
      left = 1
                                        RETURN
  100 mflag = 0
      left = ilo
                                        RETURN
  110 mflag = 1
      left = lxt
                                        RETURN
      END
c
c   This routine is required IF the CVS revision numbers are to 
c   survive an optimization.
c
c
c   $Date: 2005/08/10 15:39:26 $ $Author: hirshman $
c
      SUBROUTINE spline_rev(i)
      character*10 s 
      IF( i .eq. 0) s = 
     .'@(#)$RCSfile: spline.f,v $ $Revision: 1.18 $\000'
      RETURN
      END

