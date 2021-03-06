      MODULE fdjac_mod
      USE stel_kinds

      INTEGER, PARAMETER :: flag_singletask = -1, flag_cleanup = -100

      INTEGER :: m, n, ncnt, max_processors, num_lm_params
      REAL(rprec), DIMENSION(:), POINTER :: xp, wap
      LOGICAL, DIMENSION(:), ALLOCATABLE :: flip
      REAL(rprec) :: eps

      CONTAINS
 
      SUBROUTINE fdjac2_mp(fcn, m, n, x, fvec, fjac, ldfjac,
     1           iflag, ncnt, epsfcn, fnorm_min, x_min, wa)
      USE stel_kinds
      USE mpi_params
      IMPLICIT NONE
#if defined(MPI_OPT)
      include 'mpif.h'                                       !mpi stuff
#endif
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: m, n, ldfjac, ncnt
      INTEGER :: iflag
      REAL(rprec), INTENT(in) :: epsfcn
      REAL(rprec) :: x(n)
      REAL(rprec), INTENT(in) :: fvec(m)
      REAL(rprec) :: fjac(ldfjac,n)
      REAL(rprec) :: fnorm_min, x_min(n), wa(m)
      EXTERNAL fcn
#if defined(MPI_OPT)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, j, iproc_min, num_split, joff
      INTEGER :: iflag_array(1)
      REAL(rprec) :: eps, epsmch, dpmpar, temp, enorm, cur_norm
      REAL(rprec), PARAMETER :: one = 1, zero = 0
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: buffer, h,
     1     fnorm_array                                       !mpi stuff
      REAL(rprec), ALLOCATABLE :: fjac_overflow(:,:), fnorm_overflow(:)
      INTEGER :: istat, ierr_flag, ikey
      LOGICAL :: lflip
      LOGICAL, ALLOCATABLE     :: flip_overflow(:)
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      EXTERNAL dpmpar, enorm
C-----------------------------------------------
c
c     SUBROUTINE fdjac2
c
c     this SUBROUTINE computes a forward-difference approximation
c     to the m by n jacobian matrix ASSOCIATED with a specified
c     problem of m functions in n variables.
c
c     the SUBROUTINE statement is
c
c       SUBROUTINE fdjac2(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa)
c
c     WHERE
c
c       fcn is the name of the user-supplied SUBROUTINE which
c         calculates the functions. fcn must be declared
c         in an EXTERNAL statement in the user calling
c         program, and should be written as follows.
c
c         SUBROUTINE fcn(m,n,x,fvec,iflag)
c         INTEGER m,n,iflag
c         REAL(rprec) x(n),fvec(m)
c         ----------
c         calculate the functions at x and
c         RETURN this vector in fvec.
c         ----------
c         RETURN
c         END
c
c         the value of iflag should not be changed by fcn unless
c         the user wants to terminate execution of fdjac2.
c         in this CASE set iflag to a negative INTEGER.
c
c       m is a positive INTEGER input variable set to the number
c         of functions.
c
c       n is a positive INTEGER input variable set to the number
c         of variables. n must not exceed m.
c
c       x is an input array of length n.
c
c       fvec is an input array of length m which must contain the
c         functions evaluated at x.
c
c       fjac is an output m by n array which contains the
c         approximation to the jacobian matrix evaluated at x.
c
c       ldfjac is a positive INTEGER input variable not less than m
c         which specifies the leading DIMENSION of the array fjac.
c
c       iflag is an INTEGER variable which can be used to terminate
c         the execution of fdjac2. see description of fcn.
c
c       epsfcn is an input variable used in determining a suitable
c         step length for the forward-difference approximation. this
c         approximation assumes that the relative errors in the
c         functions are of the order of epsfcn. If epsfcn is less
c         than the machine precision, it is assumed that the relative
c         errors in the functions are of the order of the machine
c         precision.
c
c       wa is a work array of length m.
c
c     subprograms called
c
c       user-supplied ...... fcn
c
c       MINpack-supplied ... dpmpar
c
c       fortran-supplied ... abs,max,sqrt
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
      ALLOCATE (fnorm_array(n), h(n), buffer(m), stat=istat)
      IF (istat .ne. 0) STOP 'Allocation error in fdjac2_mp'
c
c     epsmch is the machine precision.
c
      epsmch = dpmpar(1)
c
      eps = SQRT(MAX(epsfcn,epsmch))
      cur_norm = enorm(m,fvec)

      h(1:n) = eps*ABS(x(1:n))
      WHERE (h .eq. zero) h=eps
      WHERE (flip) h = -h
      ierr_flag = 0


!     NOTE: if numprocs > n, num_split = 0 and we go to the overflow 
!           ("left-over") loop for num_split=mod(n,numprocs) == n
      num_split = n/numprocs
      MPI_COMM_WORKERS = MPI_COMM_WORLD
      joff = 1

      DO i = 1, num_split
         j    = myid + joff
         temp = x(j)
         x(j) = temp + h(j)

         iflag = j
         CALL fcn(m, n, x, wa, iflag, ncnt)

         x(j) = temp
         buffer(:m) = (wa(:m) - fvec(:m))/h(j)

!
!        STORE FNORM OF PERTURBED STATE (X + H)
!
         temp = enorm(m,wa)
         IF (temp > cur_norm) lflip = .not. flip(j)

!        WRITE ITERATION, PROCESSOR, CHISQ TO SCREEN
         WRITE (6, '(2x,i6,8x,i3,7x,1es12.4)') ncnt+j,
     1         myid+1, temp**2
         CALL flush(6)
      
         CALL MPI_ALLGATHER(buffer, m, MPI_REAL8, 
     1        fjac(1,joff), m, MPI_REAL8, MPI_COMM_WORLD, ierr_mpi)
         CALL MPI_ALLGATHER(temp, 1, MPI_REAL8, 
     1        fnorm_array(joff), 1, MPI_REAL8, MPI_COMM_WORLD, ierr_mpi)
         CALL MPI_ALLGATHER(lflip, 1, MPI_LOGICAL, 
     1        flip(joff), 1, MPI_LOGICAL, MPI_COMM_WORLD, ierr_mpi)

         joff = joff + numprocs
 
         IF (ierr_mpi .ne. 0) GO TO 100
      END DO

!
!     account for left-over variables OR special case when numprocs>n: 
!     must redefine the MPIC_COMM_WORDERS group to consist of MOD(n,numprocs) only
!
      num_split = MOD(n, numprocs)          !overflow...
      IF (num_split .ne. 0) THEN
      
         ALLOCATE(fjac_overflow(m,numprocs), fnorm_overflow(numprocs),
     1            flip_overflow(numprocs), stat = istat)
         IF (istat .ne. 0) THEN
            PRINT *,' istat = ',istat,' in fdjac_mod, myid: ', myid
            STOP
         END IF
         j    = myid + joff
         IF (j .le. n) THEN
            temp = x(j)
            x(j) = temp + h(j)
            ikey = WORKER_SPLIT_KEY+1
         ELSE
            temp = 0
            ikey = MPI_UNDEFINED
         END IF

!     Set new WORKERS communicator to include first MOD(n,numprocs) only
         CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, ikey, worker_id, 
     1                       MPI_COMM_WORKERS, ierr_mpi)
         IF (ierr_mpi .ne. 0) 
     1       STOP 'MPI_COMM_SPLIT error in fdjac_mp_queue'


         IF (j .le. n) THEN            !(MPI_COMM_WORKERS ONLY ACCESS HERE)
            iflag = j
            CALL fcn(m, n, x, wa, iflag, ncnt)
            x(j) = temp
            buffer(:m) = (wa(:m) - fvec(:m))/h(j)
            temp = enorm(m,wa)
            IF (temp > cur_norm) lflip = .not. flip(j)

            WRITE (6, '(2x,i6,8x,i3,7x,1es12.4)') ncnt+j,
     1             myid+1, temp**2
      
            CALL MPI_COMM_FREE(MPI_COMM_WORKERS, ierr_mpi)
         ELSE                          !(ALL OTHER PROCESSORS GO HERE)
            buffer = 0; lflip = .true.
         END IF

         CALL MPI_ALLGATHER(buffer, m, MPI_REAL8, 
     1        fjac_overflow, m, MPI_REAL8, MPI_COMM_WORLD, ierr_mpi)
         CALL MPI_ALLGATHER(temp, 1, MPI_REAL8, 
     1        fnorm_overflow, 1, MPI_REAL8, MPI_COMM_WORLD, ierr_mpi)
         CALL MPI_ALLGATHER(lflip, 1, MPI_LOGICAL, 
     1        flip_overflow, 1, MPI_LOGICAL, MPI_COMM_WORLD, ierr_mpi)

         fjac(:,joff:n)      = fjac_overflow(:,1:num_split)
         fnorm_array(joff:n) = fnorm_overflow(1:num_split)
         flip(joff:n)        = flip_overflow(1:num_split)

         DEALLOCATE(fjac_overflow, fnorm_overflow, flip_overflow)
 
      END IF
      
!
!     Find processor with minimum fnorm_min value and broadcast wa (=fvec_min), x_min, fnorm_min
!
      fnorm_min   = MINVAL(fnorm_array)
      iflag_array = MINLOC(fnorm_array)
      iproc_min = iflag_array(1)
      IF (iproc_min .le. 0 .or. iproc_min .gt. n) THEN
         PRINT *,'IPROC_MIN=',iproc_min,' out of range in fdjac2_mp'
         STOP
      END IF
      wa(:) = fjac(:, iproc_min)*h(iproc_min) + fvec(:)
      x_min(:) = x(:)
      x_min(iproc_min) = x(iproc_min) + h(iproc_min)

      DEALLOCATE (h, fnorm_array, buffer)
!
!     Do any special cleanup now for iflag = flag_cleanup
!     Barrier appears in fcn call
!
      ikey = flag_cleanup
      MPI_COMM_WORKERS = MPI_COMM_WORLD
      CALL fcn(m, n, x, wa, ikey, ncnt)

      RETURN

 100  CONTINUE
      WRITE (6, *) ' MPI_BCAST error in FDJAC2_MP: IERR=', ierr_mpi,
     1             ' PROCESSOR: ',myid
#endif
      END SUBROUTINE fdjac2_mp


      SUBROUTINE fdjac2_mp_queue(fcn, m, n, x, fvec, fjac, ldfjac,
     1           iflag, ncnt, epsfcn, fnorm_min, x_min, wa)
      USE stel_kinds
      USE mpi_params
      IMPLICIT NONE
#if defined(MPI_OPT)
      include 'mpif.h'                                       !mpi stuff
#endif
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: m, n, ldfjac, ncnt
      INTEGER :: iflag
      REAL(rprec), INTENT(in) :: epsfcn
      REAL(rprec) :: x(n)
      REAL(rprec), INTENT(in) :: fvec(m)
      REAL(rprec) :: fjac(ldfjac,n)
      REAL(rprec) :: fnorm_min, x_min(n), wa(m)
      EXTERNAL fcn
#if defined(MPI_OPT)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, j, iproc_min
      INTEGER :: iflag_array(1)
      REAL(rprec) :: eps, epsmch, dpmpar, temp, enorm, cur_norm
      REAL(rprec), PARAMETER :: one = 1, zero = 0
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: buffer, h,
     1     fnorm_array                                       !mpi stuff
      INTEGER :: status(MPI_STATUS_size)                     !mpi stuff
      INTEGER :: numsent, sender, istat, ikey                !mpi stuff
      INTEGER :: anstype, column, ierr_flag                  !mpi stuff
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      EXTERNAL dpmpar, enorm
C-----------------------------------------------
c
c     SUBROUTINE fdjac2
c
c     this SUBROUTINE computes a forward-difference approximation
c     to the m by n jacobian matrix ASSOCIATED with a specified
c     problem of m functions in n variables.
c
c     the SUBROUTINE statement is
c
c       SUBROUTINE fdjac2(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa)
c
c     WHERE
c
c       fcn is the name of the user-supplied SUBROUTINE which
c         calculates the functions. fcn must be declared
c         in an EXTERNAL statement in the user calling
c         program, and should be written as follows.
c
c         SUBROUTINE fcn(m,n,x,fvec,iflag)
c         INTEGER m,n,iflag
c         REAL(rprec) x(n),fvec(m)
c         ----------
c         calculate the functions at x and
c         RETURN this vector in fvec.
c         ----------
c         RETURN
c         END
c
c         the value of iflag should not be changed by fcn unless
c         the user wants to terminate execution of fdjac2.
c         in this case set iflag to a negative integer.
c
c       m is a positive INTEGER input variable set to the number
c         of functions.
c
c       n is a positive INTEGER input variable set to the number
c         of variables. n must not exceed m.
c
c       x is an input array of length n.
c
c       fvec is an input array of length m which must contain the
c         functions evaluated at x.
c
c       fjac is an output m by n array which contains the
c         approximation to the jacobian matrix evaluated at x.
c
c       ldfjac is a positive INTEGER input variable not less than m
c         which specifies the leading DIMENSION of the array fjac.
c
c       iflag is an INTEGER variable which can be used to terminate
c         the execution of fdjac2. see description of fcn.
c
c       epsfcn is an input variable used in determining a suitable
c         step length for the forward-difference approximation. this
c         approximation assumes that the relative errors in the
c         functions are of the order of epsfcn. If epsfcn is less
c         than the machine precision, it is assumed that the relative
c         errors in the functions are of the order of the machine
c         precision.
c
c       wa is a work array of length m.
c
c     subprograms called
c
c       user-supplied ...... fcn
c
c       MINpack-supplied ... dpmpar
c
c       fortran-supplied ... abs,max,sqrt
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
      ALLOCATE (buffer(m), fnorm_array(n), h(n), stat=istat)
      IF (istat .ne. 0) STOP 'Allocation error in fdjac2_mp'
      ierr_flag = 0
c
c     epsmch is the machine precision.
c
      epsmch = dpmpar(1)
c
      eps = SQRT(MAX(epsfcn,epsmch))

c     Set up WORKER COMMUNICATOR
      IF (myid .eq. master) THEN
         ikey = MPI_UNDEFINED
      ELSE
         ikey = WORKER_SPLIT_KEY+1
      END IF
      CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, ikey, worker_id, 
     1                    MPI_COMM_WORKERS, ierr_mpi)
      IF (ierr_mpi .ne. 0) STOP 'MPI_COMM_SPLIT error in fdjac_mp_queue'

c
c     Begin parallel MPI master/worker version of the function call.  A
c     bank queue or master/worker MPI algorithm is used here to achieve
c     parallel load balancing in the somewhat uneven work load involved
c     in calculating the Jacobian function needed in the Levenberg-Marquardt
c     optimization.  This model is based on the example given in Chapt. 5 of
c     "The LAM Companion to Using MPI" by Zdzislaw Meglicki (online see:
c     http://carpanta.dc.fi.udc.es/docs/mpi/mpi-lam/mpi.html).  These
c     modifications were made by D.A. Spong 8/23/2000.
c
c     ****Master portion of the code****
c


      IF (myid .eq. master) THEN
         numsent = 0    !numsent is a counter used to track how many
                        !jobs have been sent to workers
         cur_norm = enorm(m,fvec)
c
c        calculate the displacements
c
         h(1:n) = eps*ABS(x(1:n))
         WHERE (h .eq. zero) h=eps
         WHERE (flip) h = -h

c     Send forward difference displacements from master to each
c           worker process. Tag these with the column number (j)
c
         DO j = 1,MIN(numprocs-1,n)
            numsent = numsent+1
            temp = x(numsent)
            x(numsent) = temp + h(numsent)
            CALL MPI_SEND(x, n, MPI_REAL8, j,
     1                  numsent, MPI_COMM_WORLD, ierr_mpi)
            IF (ierr_mpi .ne. 0) STOP 'MPI_SEND error(1) in fdjac2'
            x(numsent) = temp
         END DO          !j = 1,MIN(numprocs-1,n)
c
c      Looping through the columns, collect answers from the workers.
c      As answers are received, new uncalculated columns are sent
c      out to these same workers.
c
         DO j = 1,n
            CALL MPI_RECV(wa, m, MPI_REAL8,
     1           MPI_ANY_SOURCE, MPI_ANY_TAG,
     2           MPI_COMM_WORLD, status, ierr_mpi)
            IF (ierr_mpi .ne. 0) STOP 'MPI_RECV error(1) in fdjac2'
            sender     = status(MPI_SOURCE)
            anstype    = status(MPI_TAG)       ! column is tag value
            IF (anstype .gt. n) STOP 'ANSTYPE > N IN FDJAC2'

            fjac(:m,anstype) = (wa(:m) - fvec(:m))/h(anstype)
!
!           STORE FNORM OF PERTURBED STATE (X + H)
!
            temp = enorm(m,wa)
            fnorm_array(anstype) = temp
            IF (temp > cur_norm) flip(anstype) = .not. flip(anstype)

!           WRITE ITERATION, PROCESSOR, CHISQ TO SCREEN
            WRITE (6, '(2x,i6,8x,i3,7x,1es12.4)') ncnt+anstype,
     1             sender, temp**2

c
c           If more columns are left, then send another column to the worker(sender)
c           that just sent in an answer
c
            IF (numsent .lt. n) THEN
               numsent = numsent+1
               temp = x(numsent)
               x(numsent) = temp + h(numsent)

               CALL MPI_SEND(x, n, MPI_REAL8,
     1                       sender, numsent, MPI_COMM_WORLD, ierr_mpi)
               IF (ierr_mpi .ne. 0) STOP 'MPI_SEND error(2) in fdjac2'
               x(numsent) = temp

            ELSE                ! Tell workers that there is no more work to do

               CALL MPI_SEND(MPI_BOTTOM, 0, MPI_REAL8,
     1                       sender, 0, MPI_COMM_WORLD, ierr_mpi)
               IF (ierr_mpi .ne. 0) STOP 'MPI_end error(3) in fdjac2'
            ENDIF      ! IF( myid .eq. master ) THEN
         END DO     ! DO j = 1,n
c
c     ****Worker portion of the code****
c        Skip this when processor id exceeds work to be done
c
      ELSE IF (myid .le. n) THEN        ! IF( myid .ne. master )
c
c        Otherwise accept the next available column, check the tag,
c        and if the tag is non-zero call subroutine fcn.
c        If the tag is zero, there are no more columns
c        and worker skips to the end.
c
 90      CALL MPI_RECV(x, n, MPI_REAL8, master,
     1                 MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr_mpi)
         IF (ierr_mpi .ne. 0) STOP 'MPI_RECV error(2) in fdjac2'

         column = status(MPI_TAG)
         IF (column .ne. 0) THEN
            iflag = column
c           Call the chisq fcn for the portion of displacement vector which
c           was just received. Note that WA stores the local fvec_min array
            CALL fcn(m, n, x, wa, iflag, ncnt)

            IF (iflag.ne.0 .and. ierr_flag.eq.0) ierr_flag = iflag
c
c           Send this function evaluation back to the master process tagged
c           with the column number so the master knows where to put it
c
            CALL MPI_SEND(wa, m, MPI_REAL8, master,
     1                    column, MPI_COMM_WORLD, ierr_mpi)
            IF (ierr_mpi .ne. 0) STOP 'MPI_SEND error(4) in fdjac2'
            GO TO 90    !Return to 90 and check if master process has sent any more jobs
         END IF
      ENDIF       ! IF( myid .ne. master )

!
!     Broadcast the fjac matrix
!
      DO j=1,n
         IF (myid .eq. master) buffer(:m) = fjac(:m,j)
         CALL MPI_BCAST(buffer, m, MPI_REAL8, master,
     1        MPI_COMM_WORLD, ierr_mpi)
         IF (ierr_mpi .ne. 0) GO TO 100
         IF (myid .ne. master) fjac(:m,j) = buffer(:m)
      END DO

!
!     Find processor with minimum fnorm_min value and broadcast wa (=fvec_min), x_min, fnorm_min
!
      IF (myid .eq. master) THEN
         fnorm_min   = MINVAL(fnorm_array)
         iflag_array = MINLOC(fnorm_array)
         iproc_min = iflag_array(1)
         IF (iproc_min .le. 0 .or. iproc_min .gt. n) THEN
            PRINT *,'IPROC_MIN=',iproc_min,' out of range in fdjac2_mp'
            STOP
         END IF
         wa(:) = fjac(:, iproc_min)*h(iproc_min) + fvec(:)
         x_min(:) = x(:)
         x_min(iproc_min) = x(iproc_min) + h(iproc_min)
      END IF

      CALL MPI_BCAST(fnorm_min, 1, MPI_REAL8, master,
     1               MPI_COMM_WORLD, ierr_mpi)
      IF (ierr_mpi .ne. 0) GOTO 100
      CALL MPI_BCAST(wa, m, MPI_REAL8, master, MPI_COMM_WORLD,ierr_mpi)
      IF (ierr_mpi .ne. 0) GOTO 100
      CALL MPI_BCAST(x_min, n, MPI_REAL8, master, MPI_COMM_WORLD, 
     1               ierr_mpi)
      IF (ierr_mpi .ne. 0) GOTO 100
c
c     Original serial version of the Jacobian calculation:
c
c      DO j = 1, n
c         temp = x(j)
c         h = eps*ABS(temp)
c         IF (h .eq. zero) h = eps
c         x(j) = temp + h
c         CALL fcn (m, n, x, wa, iflag)
c         IF (iflag .lt. 0) EXIT
c         x(j) = temp
c         fjac(:m,j) = (wa - fvec)/h
c      END DO

      DEALLOCATE (h, buffer, fnorm_array)
!
!     Do any special cleanup now for iflag = flag_cleanup
!     Barrier appears in fcn call
!
      column = flag_cleanup
      CALL fcn(m, n, x, wa, column, ncnt)

!
!     Reassign initial x value to all processors and perform error handling
!     This is necessary because other processors had x + h in them
!
      CALL MPI_BCAST(x, n, MPI_REAL8, master, MPI_COMM_WORLD, ierr_mpi)
      IF (ierr_mpi .ne. 0) GOTO 100

      IF (ierr_flag .ne. 0) THEN
         iflag = ierr_flag
         CALL MPI_BCAST(iflag, 1, MPI_INTEGER, myid,
     1        MPI_COMM_WORLD, ierr_mpi)
         IF (ierr_mpi .ne. 0) GOTO 100
      END IF

      IF (myid .ne. master)
     1   CALL MPI_COMM_FREE(MPI_COMM_WORKERS, ierr_mpi)      

      RETURN

 100  CONTINUE
      WRITE (6, *) ' MPI_BCAST error in FDJAC2_MP: IERR=', ierr_mpi,
     1             ' PROCESSOR: ',myid
#endif
      END SUBROUTINE fdjac2_mp_queue


      SUBROUTINE fdjac2(fcn, m1, n1, x, fvec, fjac, ldfjac, iflag,
     1    ncnt1, epsfcn, wa, time, fnorm_min, x_min, fvec_min)
      USE stel_kinds
#if !defined(MPI_OPT)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: m1, n1, ldfjac, ncnt1
      INTEGER, TARGET :: iflag
      REAL(rprec) ::  epsfcn, time
      REAL(rprec), TARGET :: x(n), wa(m)
      REAL(rprec), INTENT(in) :: fvec(m)
      REAL(rprec), INTENT(out) :: fjac(ldfjac,n)
      REAL(rprec), INTENT(out) :: fnorm_min, x_min(n), fvec_min(m)
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: zero = 0
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: j, istat, iread, ic1, ic2, irate, count_max
      REAL(rprec) :: eps, epsmch, h, dpmpar, temp, cur_norm
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
      EXTERNAL fcn, dpmpar, multiprocess, fdjac_parallel
      REAL(rprec), EXTERNAL :: enorm
C-----------------------------------------------
c
c     SUBROUTINE fdjac2
c
c     this SUBROUTINE computes a forward-difference approximation
c     to the m by n jacobian matrix ASSOCIATED with a specified
c     problem of m functions in n variables.
c
c     Here
c
c       fcn is the name of the user-supplied SUBROUTINE which
c         calculates the functions. fcn must be declared
c         in an external statement in the user calling
c         program (see LMDIF1 for documentation), and should be written as follows:
c
c         SUBROUTINE fcn(m,n,x,fvec,iflag,ncnt)
c         INTEGER m,n,iflag
c         REAL(rprec) x(n),fvec(m)
c         ----------
c         calculate the functions at x and
c         RETURN this vector in fvec.
c         ----------
c         RETURN
c         END
c
c       fjac is an output m by n array which CONTAINS the
c         approximation to the jacobian matrix evaluated at x.
c
c       ldfjac is a positive INTEGER input variable not less than m
c         which specifies the leading DIMENSION of the array fjac.
c
c       iflag is an INTEGER variable which can be used to terminate
c         the execution of fdjac2. see description of fcn.
c
c       epsfcn is an input variable used in determining a suitable
c         step length for the forward-difference approximation. this
c         approximation assumes that the relative errors in the
c         functions are of the order of epsfcn. IF epsfcn is less
c         than the machine precision, it is assumed that the relative
c         errors in the functions are of the order of the machine
c         precision.
c
c       wa is a work array of length m.
c
c     subprograms called
c
c       user-supplied ...... fcn
c
c       MINpack-supplied ... dpmpar
c
c       fortran-supplied ... ABS,max,sqrt
c
c     argonne national laboratory. MINpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********

c
c     epsmch is the machine precision.
c
      epsmch = dpmpar(1)
c
      eps = SQRT(MAX(epsfcn,epsmch))
!
!     Load fdjac_mod values. pointers will automatically update targets...
!     Prepare for multi-processing...
!
      m = m1
      n = n1
      ncnt = ncnt1
      xp => x
      wap => wa

!     Find min chisq = fnorm**2 state for this jacobian evaluation
!     (Do NOT retain from previous evaluations, or could get into a non-converging loop...)

      fnorm_min = HUGE(fnorm_min)

      CALL SYSTEM_CLOCK(ic1, irate)
      CALL multiprocess(n, max_processors, fdjac_parallel, fcn)
      CALL SYSTEM_CLOCK(ic2, irate, count_max)
      IF (ic2 .lt. ic1) ic2 = ic2 + count_max

c     DO j = 1, n
c         temp = x(j)
c         h = eps*ABS(temp)
c         IF (h .eq. zero) h = eps
c         x(j) = temp + h
c         CALL fcn (m, n, x, wa, iflag, ncnt)
c         IF (iflag .lt. 0) EXIT
c         x(j) = temp
c         fjac(:m,j) = (wa - fvec)/h
c      END DO


      cur_norm = enorm(m,fvec)      ! where are we now?

      DO j = 1, n

        READ (j+1000, iostat=iread) istat, iflag, h, temp
        IF (iread .ne. 0) THEN
           WRITE (6, *) 'Error reading from file fort.', j+1000,
     1     ' in fdjac2: IOSTAT = ', iread
           iflag = -14
        ELSE IF (j .ne. istat) THEN
           WRITE (6, *) 'Wrong value for INDEX j READ in fdjac2'
           iflag = -14
        END IF

        IF (iflag .ne. 0) EXIT

#if defined(CRAY)
        DO k = 1, m
           READ (j+1000) wa(k)
        END DO
#else
        READ (j+1000) wa
#endif
        fjac(:m,j) = (wa - fvec)/h

        IF( temp > cur_norm) flip(j) = .not. flip(j)  ! flip for next time

        IF( temp < fnorm_min) THEN
           fnorm_min = temp
           fvec_min = wa
#if defined(CRAY)
           DO k = 1, n
              READ (j+1000) x_min(k)
           END DO
#else
           READ (j+1000) x_min
#endif
        END IF

        CLOSE (j+1000, status='delete')                        !!Needed to run correctly in multi-tasking...

      END DO

!
!     Do any special cleanup now for iflag = flag_cleanup
!
      iflag = flag_cleanup
      CALL fcn(m, n, x, wa, iflag, ncnt)

      time = time + REAL(ic2 - ic1)/REAL(irate)                !!Time in multi-process CALL
#endif
      END SUBROUTINE fdjac2


      END MODULE fdjac_mod

