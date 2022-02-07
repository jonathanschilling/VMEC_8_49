!     SPH011108: replace all INTEGER(iprec) with INTEGER
!============================================================================
      MODULE ip_beamline
!--------------------------------------------------------------------------
!  Module to define data TYPEs for Faraday Rotation
!-------------------------------------------------------------------------
      USE stel_kinds
      USE stel_constants
      USE v3_utilities       ! assert_eq is here
      IMPLICIT NONE


!==========================================================================
      TYPE ip_beam
!-----------------------------------------------------------------------
! declares variables relating to the interferometer/polarimeter beam 
! path ("q") through the plasma.
! 
!   q_unit(3)   : unit vector in the direction of the path
!   q0vec(3)    : Cartesian vector of path starting point 
!   qfvec(3)    : Cartesian vector of path end point
!   q0          : path distance at path start point (often assigned to be 0.0)
!   qdist       : path length (distance btwn q0 and qf)
!   wavelength  : wavelength of inferometer/polarimeter beam
!  B_ratio(i,j,k) : precomputed B/I (mag field/current) for points along the beam path
!               :  i = index to specify the "step point" along the path
!               :  j = index to specify the desired "coil collection" 
!               :  k = index (1:3) to specify X,Y, or Z component of the B field
! s_name        : character variable for "short name" of the int/polarimeter beam
! l_name        : character variable for "long name" of the int/polarimeter beam
!------------------------------------------------------------------------------------- 
      REAL(rprec), DIMENSION(3)              :: q_unit
      REAL(rprec), DIMENSION(3)              :: q0vec
      REAL(rprec), DIMENSION(3)              :: qfvec
      REAL(rprec)                            :: q0
      REAL(rprec)                            :: qdist
      REAL(rprec)                            :: wavelength
      REAL(rprec), DIMENSION(:,:,:), POINTER :: B_ratio => null()
      CHARACTER (len=30)                     :: s_name                                 
      CHARACTER (len=80)                     :: l_name
      END TYPE ip_beam




!===============================================================================
      TYPE int_pol
!---------------------------------------------------------------------------------
! declares variables relating to the microwave/laser beam of the interferometer
!
! nbeams     : integer variable to store the number of individual microwave/laser beams
! ipbeam     : ip_beam derived TYPE variable to store info for each individual beam
! s_name     : character variable for "short name" of the interferometer/polarimeter array
! l_name     : character variable for "long name" of the interferometer/polarimeter array
!
!--------------------------------------------------------------------------------
      INTEGER                              :: nbeams
      TYPE(ip_beam), DIMENSION(:), POINTER :: ipbeam => null() 
      CHARACTER (len=30)                   :: s_name                                 
      CHARACTER (len=80)                   :: l_name
      END TYPE int_pol


!===============================================================================
      TYPE int_pol_coll
!---------------------------------------------------------------------------------
! declares variables relating to a COLLECTION of interferometer/polarimeter arrays
! n_ip       : number of int/pol units contained in the collection
! ip         : int_pol TYPE variable to store info about each int/pol unit in the collection
! s_name     : character variable for "short name" of the interferometer/polarimeter array
! l_name     : character variable for "long name" of the interferometer/polarimeter array
!
!--------------------------------------------------------------------------------
      INTEGER                              :: n_ip
      TYPE(int_pol), DIMENSION(:), POINTER :: ip => null() 
      CHARACTER (len=30)                   :: s_name                                 
      CHARACTER (len=80)                   :: l_name
      END TYPE int_pol_coll
     



      INTERFACE ASSIGNMENT (=)
         MODULE PROCEDURE ipbeam_to_ipbeam, ipbeam_a_to_ipbeam_a,              &
     &                    int_pol_to_int_pol, int_pol_a_to_int_pol_a,          &
     &                    int_pol_coll_to_int_pol_coll
      END INTERFACE




!----------------------------------------------------------------------------
!  subroutines for use in module beamline
!---------------------------------------------------------------------------
      CONTAINS

!===============================================================================
      SUBROUTINE int_pol_coll_construct( this, n_ip, ip_array,                &
     &                                   s_name, l_name)

      IMPLICIT NONE

!  Declare Input Arguments 
      TYPE (int_pol_coll), INTENT (inout)             :: this
      INTEGER, INTENT(in)                             :: n_ip
      TYPE (int_pol), DIMENSION(n_ip), INTENT(in)     :: ip_array
      CHARACTER (len=*), INTENT(in)                   :: s_name
      CHARACTER (len=*), INTENT(in)                   :: l_name

!  Declare local variables
      INTEGER :: arr_size, ier1
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'int_pol_coll_construct: '

!  Start of executable code

!      WRITE(*,*) ' Executing int_pol_coll_construct'


      CALL int_pol_coll_destroy(this)


!........verify that n_ip matches actual ip_array size......!
      arr_size = SIZE(ip_array)
      if ( arr_size .NE. n_ip ) then
        write(*,*) sub_name, 'ERROR.  n_ip NOT = to array size'
      else
!...........allocate memory for the ip pointer array...............!
        ALLOCATE( this % ip(n_ip) )

        this % n_ip   = n_ip
        this % s_name = TRIM(s_name)
        this % l_name = TRIM(l_name)

!.......we actually need to "manually" construct each ip() array element SEPARATELY
!.......based on the sizes of the individual elements in array ip_array().........!
        this % ip  = ip_array
      end if



      END SUBROUTINE int_pol_coll_construct


!===============================================================================
      SUBROUTINE int_pol_construct(this, nbeams, beam_array, s_name,           &
     &                             l_name)

      IMPLICIT NONE

!  Declare Input Arguments 
      TYPE (int_pol), INTENT (inout)                  :: this
      INTEGER, INTENT(in)                             :: nbeams
      TYPE (ip_beam), DIMENSION(nbeams), INTENT(in)   :: beam_array
      CHARACTER (len=*), INTENT(in)                   :: s_name
      CHARACTER (len=*), INTENT(in)                   :: l_name

!  Declare local variables
      INTEGER :: arr_size, ier1
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'int_pol_construct: '

!  Start of executable code

!      WRITE(*,*) ' Executing int_pol_construct'


      CALL int_pol_destroy(this)


!........verify that nbeams matches actual beam_array size......!
      arr_size = SIZE(beam_array)
      if ( arr_size .NE. nbeams ) then
        write(*,*) sub_name, 'ERROR.  nbeams NOT = to array size'
      else
!...........allocate memory for the ipbeam pointer array...............!
        ALLOCATE( this % ipbeam(nbeams) )

        this % nbeams = nbeams
        this % s_name = TRIM(s_name)
        this % l_name = TRIM(l_name)

!.......we actually need to "manually" construct each ipbeam() array element SEPARATELY
!.......based on the sizes of the individual elements in array beam_array().........!
        this % ipbeam = beam_array
      end if



      END SUBROUTINE int_pol_construct






!===============================================================================
      SUBROUTINE ip_beam_construct(this, q0vec, qfvec, q0,wavelength,          &
     &                             B_ratio, s_name, l_name)

      USE math_utilities
      IMPLICIT NONE

!  Declare Arguments 
      TYPE (ip_beam), INTENT (inout)          :: this


      REAL(rprec), DIMENSION(3), INTENT(in)     :: q0vec
      REAL(rprec), DIMENSION(3), INTENT(in)     :: qfvec
      REAL(rprec), INTENT(in)                   :: q0
      REAL(rprec), INTENT(in)                   :: wavelength
      REAL(rprec), DIMENSION(:,:,:), INTENT(in) :: B_ratio
      CHARACTER (len=*), INTENT(in)             :: s_name
      CHARACTER (len=*), INTENT(in)             :: l_name

!  Declare local variables
      INTEGER, DIMENSION(3) :: dim
      INTEGER :: nd1, nd2, ier1
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'ip_beam_construct: '

      REAL(rprec), DIMENSION(3)       :: q_unit
      REAL(rprec)                     :: qdist

!  Start of executable code

!      WRITE(*,*) ' Executing ip_beam_construct'

!  The eq_param_var type contains pointers to arrays. Therefore,
!  to deallocate space, need to destroy this
      CALL ip_beam_destroy(this)

 
!       q_unit = 0.0
!       qdist = 0.0
      q_unit = ( qfvec - q0vec)/ MAGNITUDE(qfvec - q0vec)
      qdist  = DIST(q0vec, qfvec)

      this % q_unit      = q_unit
      this % q0vec       = q0vec
      this % qfvec       = qfvec
      this % q0          = q0
      this % qdist       = qdist
      this % wavelength  = wavelength
      this % s_name      = TRIM(s_name)
      this % l_name      = TRIM(l_name)

!......allocate an array with the same size & shape of the input array.......!      
      dim(1:3) = SHAPE(B_ratio)
      ALLOCATE( this%B_ratio( dim(1), dim(2), dim(3) ) )

      this % B_ratio = B_ratio

      END SUBROUTINE ip_beam_construct


!==============================================================================
      SUBROUTINE ip_beam_destroy(this)
      IMPLICIT NONE

!  Declare Arguments 
      TYPE (ip_beam), INTENT(inout) :: this

!  Declare local variables
      INTEGER :: ier1
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'ip_beam_destroy: '

!  Start of executable code

!  Get rid of all components
      
      this % q_unit     = zero
      this % q0vec      = zero
      this % qfvec      = zero
      this % q0         = zero
      this % qdist      = zero
      this % wavelength = zero
      this % s_name = ''
      this % l_name = ''

!.........deallocate the B_ratio array (and nullify the pointer itself as well).......!
      IF ( ASSOCIATED( this % B_ratio) ) THEN
         DEALLOCATE(this % B_ratio,STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'dealloc int')
      ENDIF

      END SUBROUTINE ip_beam_destroy


!==================================================================================
      SUBROUTINE ip_beam_destroy_a(this)
!-------------------------------------------------------------------------------
!  FUNCTION: Destroys an array of ip_beam TYPE variables
!-------------------------------------------------------------------------------
!  Declare Arguments 
      TYPE (ip_beam), DIMENSION(:), INTENT(inout) :: this
      INTEGER :: i, nsize
      
      nsize = SIZE(this)
      DO i = 1, nsize
         CALL ip_beam_destroy( this(i) )
      END DO
      
      END SUBROUTINE ip_beam_destroy_a




!==============================================================================
      SUBROUTINE int_pol_destroy(this)
      IMPLICIT NONE

!  Declare Arguments 
      TYPE(int_pol), INTENT(inout) :: this

!  Declare local variables
      INTEGER :: ier1, arr_size
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'int_pol_destroy: '

!  Start of executable code

      this % nbeams = zero
      this % s_name = ''
      this % l_name = ''

!!.........deallocate the ipbeam array (and nullify the pointer itself as well).......!
!      IF ( ASSOCIATED( this % ipbeam) ) THEN
!         DEALLOCATE(this % ipbeam,STAT=ier1)
!         CALL assert_eq(0,ier1,sub_name // 'dealloc int')
!      ENDIF




!......deallocate the ipbeam array element by element to avoid memory leaks....!
      IF (ASSOCIATED(this % ipbeam)) THEN
         arr_size = SIZE(this % ipbeam)
         CALL ip_beam_destroy_a(this % ipbeam(1:arr_size))         
         DEALLOCATE(this % ipbeam, STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'dealloc int')
      END IF


      END SUBROUTINE int_pol_destroy




!==================================================================================
      SUBROUTINE int_pol_destroy_a(this)
!-------------------------------------------------------------------------------
!  FUNCTION: Destroys an array of int_pol TYPE variables
!-------------------------------------------------------------------------------
!  Declare Arguments 
      TYPE (int_pol), DIMENSION(:), INTENT(inout) :: this
      INTEGER :: i, nsize
      
      nsize = SIZE(this)
      DO i = 1, nsize
         CALL int_pol_destroy( this(i) )
      END DO
      
      END SUBROUTINE int_pol_destroy_a



!==============================================================================
      SUBROUTINE int_pol_coll_destroy(this)
      IMPLICIT NONE

!  Declare Arguments 
      TYPE(int_pol_coll), INTENT(inout) :: this

!  Declare local variables
      INTEGER :: ier1, arr_size
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'int_pol_coll_destroy: '

!  Start of executable code

      this % n_ip   = zero
      this % s_name = ''
      this % l_name = ''

!!.........deallocate the ip array (and nullify the pointer itself as well).......!
!      IF ( ASSOCIATED( this % ip ) ) THEN
!         DEALLOCATE(this % ip, STAT=ier1)
!         CALL assert_eq(0,ier1,sub_name // 'dealloc int')
!      ENDIF




!......deallocate the ip array element by element to avoid memory leaks....!
      IF (ASSOCIATED(this % ip)) THEN
         arr_size = SIZE(this % ip)
         CALL int_pol_destroy_a(this % ip(1:arr_size))         
         DEALLOCATE(this % ip, STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'dealloc int')
      END IF


      END SUBROUTINE int_pol_coll_destroy




!=====================================================================================
      SUBROUTINE ipbeam_to_ipbeam(left,right)
!---------------------------------------------------------------------------------
!
! FUNCTION:  overloads the '=' operator to execute  left = right   for variables of TYPE ip_beam
!----------------------------------------------------------------------------------
      IMPLICIT NONE
      
!  Declare Arguments 
      TYPE (ip_beam), INTENT(out) :: left
      TYPE (ip_beam), INTENT(in) :: right
      
!  Declare local variables
      INTEGER :: size1, size2, size3
      LOGICAL, DIMENSION(5) :: lassert
      LOGICAL :: assert
      CHARACTER(len =*), PARAMETER :: subname = 'ipbeam_to_ipbeam: '

!  Start of executable code
!  Check to see if the 'use as allocatable array, pointers are pointing
!  to the same location (i.e tests against  A = A case)
      lassert(1) = ASSOCIATED(left % B_Ratio, right % B_ratio)

      assert = lassert(1)

      IF(assert .eqv. .TRUE.) THEN
!...............write error message and exit subroutine....................!
         write(*,*) subname, 'WARNING: left beam equal to right beam' 
      ELSE

!.....the destroy routines are giving me trouble, so destroy "manually".  JS 9/26/05
!!........destroy the current "left" beam and reassign it values from "right".....!
!        call ip_beam_destroy(left)

!.........copy all Non-pointer variables.................................!
        left % q_unit     = right % q_unit
        left % q0vec      = right % q0vec
        left % qfvec      = right % qfvec
        left % q0         = right % q0
        left % qdist      = right % qdist
        left % wavelength = right % wavelength
        left % s_name     = right % s_name
        left % l_name     = right % l_name


!...........copy Three dimensional pointer array......................................!
        size1 = SIZE( right % B_ratio , 1)
        size2 = SIZE( right % B_ratio , 2)
        size3 = SIZE( right % B_ratio , 3)
        IF(ASSOCIATED(left % B_ratio)) DEALLOCATE(left % B_ratio)
        ALLOCATE(left % B_ratio(size1, size2, size3) )
        left % B_ratio = right % B_ratio

      END IF 
         
      END SUBROUTINE ipbeam_to_ipbeam



!===============================================================================
      SUBROUTINE ipbeam_a_to_ipbeam_a(left,right)
!---------------------------------------------------------------------------------------------
!  FUNCTION:  overloading of the "=" operator to excute  left = right assignment for TYPE ip_beam
!
!  NOTE:  arrays must be ONE dimensional and have equal sizes.  IF left is left(a:b) and right
!         is right(c:d)  then left(a+1) = right(c+1) and so on  (assuming  b-a = d-c )
!
!----------------------------------------------------------------------------------------
      IMPLICIT NONE
!      
!  Declare Arguments 
      TYPE (ip_beam), DIMENSION(:), INTENT (out) :: left
      TYPE (ip_beam), DIMENSION(:), INTENT (in) :: right
      
!  Declare temporary variables
      INTEGER :: i, nleft, nright, low_b, high_b, b_diff
      CHARACTER(len =*), PARAMETER :: subname = 'ipbeam_a_to_ipbeam_a'
         
!  Start of executable code
      nleft = SIZE(left)
      nright = SIZE(right)
      IF(nleft .ne. nright) THEN
         write(*,*) subname, 'WARNING: nleft not equal to nright' 
      END IF

!.......account for possibility that array indices are offset ( e.g.  0:4 vs 1:5)
      low_b  = LBOUND(left, 1)
      high_b  = UBOUND(left, 1)
      b_diff = LBOUND(right, 1) - low_b

      do i= low_b, high_b
        left(i) = right(i + b_diff)
      end do

      END SUBROUTINE ipbeam_a_to_ipbeam_a


!=====================================================================================
      SUBROUTINE int_pol_to_int_pol(left,right)
!---------------------------------------------------------------------------------
!
! FUNCTION:  overloads the '=' operator to execute  left = right   for variables of TYPE int_pol
!----------------------------------------------------------------------------------
      IMPLICIT NONE
      
!  Declare Arguments 
      TYPE (int_pol), INTENT(out) :: left
      TYPE (int_pol), INTENT(in) :: right
      
!  Declare local variables
      INTEGER :: r_size, tempsize
      LOGICAL, DIMENSION(5) :: lassert
      LOGICAL :: assert
      CHARACTER(len =*), PARAMETER :: subname = 'int_pol_to_int_pol: '

!  Start of executable code
!  Check to see if the 'use as allocatable array, pointers are pointing
!  to the same location (i.e tests against  A = A case)

      lassert(1) = ASSOCIATED(left%ipbeam, right%ipbeam)
      assert = lassert(1) 

      IF(assert .eqv. .TRUE.) THEN
!...............write error message and exit subroutine....................!
         write(*,*) subname,'WARNING: left ipbeam same as right ipbeam' 
      ELSE

!.....the destroy routines are giving me trouble, so destroy "manually".  JS 9/26/05
!!........destroy the current "left" beam and reassign it values from "right".....!
!        call int_pol_destroy(left)

!.........copy all Non-pointer variables.................................!
        left % nbeams  = right % nbeams
        left % s_name  = right % s_name
        left % l_name  = right % l_name


!...........now deal with the one dimensional pointer array...............................!

!....deallocate the existing left array, so that left & right sizes will match...........!
        IF(ASSOCIATED(left % ipbeam) ) DEALLOCATE(left % ipbeam)

!...........allocate "left" to the same size as "right"............................!
        r_size = SIZE( right % ipbeam )
        ALLOCATE(left % ipbeam(r_size) )

        left % ipbeam = right % ipbeam
      END IF 
         
      END SUBROUTINE int_pol_to_int_pol


!===============================================================================
      SUBROUTINE int_pol_a_to_int_pol_a(left,right)
!---------------------------------------------------------------------------------------------
!  FUNCTION:  overloading of the "=" operator to excute   left = right  assignment for TYPE ip_beam
!
!  NOTE:  arrays must be ONE dimensional and have equal sizes.  IF left is left(a:b) and right
!         is right(c:d)  then left(a+1) = right(c+1) and so on  (assuming  b-a = d-c )
!
!----------------------------------------------------------------------------------------
      IMPLICIT NONE
!      
!  Declare Arguments 
      TYPE (int_pol), DIMENSION(:), INTENT (out) :: left
      TYPE (int_pol), DIMENSION(:), INTENT (in) :: right
      
!  Declare temporary variables
      INTEGER :: i, nleft, nright, low_b, high_b, b_diff
      CHARACTER(len =*), PARAMETER :: subname = 'int_pol_a_to_int_pol_a'
         
!  Start of executable code
      nleft = SIZE(left)
      nright = SIZE(right)

      IF(nleft .ne. nright) THEN
         write(*,*) subname, 'WARNING: nleft not equal to nright' 

      ELSE
        low_b  = LBOUND(left, 1)
        high_b  = UBOUND(left, 1)
        b_diff = LBOUND(right, 1) - low_b

        do i= low_b, high_b
          left(i) = right(i + b_diff)
        end do
      END IF

      END SUBROUTINE int_pol_a_to_int_pol_a



!=====================================================================================
      SUBROUTINE int_pol_coll_to_int_pol_coll(left,right)
!---------------------------------------------------------------------------------
!
! FUNCTION:  overloads the '=' operator to execute  left = right   for variables of TYPE int_pol
!----------------------------------------------------------------------------------
      IMPLICIT NONE
      
!  Declare Arguments 
      TYPE (int_pol_coll), INTENT(out) :: left
      TYPE (int_pol_coll), INTENT(in) :: right
      
!  Declare local variables
      INTEGER :: r_size, tempsize
      LOGICAL, DIMENSION(5) :: lassert
      LOGICAL :: assert
      CHARACTER(len =*), PARAMETER ::                                          &
     &                   subname = 'int_pol_coll_to_int_pol_coll: '

!  Start of executable code
!  Check to see if the 'use as allocatable array, pointers are pointing
!  to the same location (i.e tests against  A = A case)

      lassert(1) = ASSOCIATED(left%ip, right%ip)
      assert = lassert(1) 

      IF(assert .eqv. .TRUE.) THEN
!...............write error message and exit subroutine....................!
         write(*,*) subname,'WARNING: left ip same as right ip' 
      ELSE

!.....the destroy routines are giving me trouble, so destroy "manually".  JS 9/26/05
!!........destroy the current "left" beam and reassign it values from "right".....!
!        call int_pol_coll_destroy(left)

!.........copy all Non-pointer variables.................................!
        left % n_ip    = right % n_ip
        left % s_name  = right % s_name
        left % l_name  = right % l_name


!...........now deal with the one dimensional pointer array...............................!

!....deallocate the existing left array, so that left & right sizes will match...........!
        IF(ASSOCIATED(left % ip) ) DEALLOCATE(left % ip)

!...........allocate "left" to the same size as "right"............................!
        r_size = SIZE( right % ip )
        ALLOCATE(left % ip(r_size) )

        left % ip = right % ip
      END IF 
         
      END SUBROUTINE int_pol_coll_to_int_pol_coll


      END MODULE ip_beamline


