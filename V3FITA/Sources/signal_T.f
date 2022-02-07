!*******************************************************************************
!  File signal_T.f
!  Contains module signal_T
!  Defines derived-types: signal_desc, signal_data

!*******************************************************************************
!  MODULE signal_T
!    (Signal Type Definition, for the V3FIT code)
! SECTION I.    VARIABLE DECLARATIONS
! SECTION II.   DERIVED-TYPE DECLARATIONS
! SECTION III.  INTERFACE BLOCKS
! SECTION IV.   CONSTRUCTION SUBROUTINES
! SECTION V.    DESTRUCTION SUBROUTINES
! SECTION VI.   ASSIGNMENT SUBROUTINES
! SECTION VII.  OUTPUT SUBROUTINES

! SECTION XII.  AUXILIARY SUBROUTINES
! SECTION XIII. DEBUGGING SUBROUTINES
! SECTION XV.   DUPLICATE CODING FOR TESTING
! SECTION XVI.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
      MODULE signal_T

!*******************************************************************************
! SECTION I. VARIABLE DECLARATIONS
!*******************************************************************************

!-------------------------------------------------------------------------------
!  Type declarations - lengths of reals, integers, and complexes.
!  Frequently used mathematical constants, lots of extra precision.
!-------------------------------------------------------------------------------
      USE stel_kinds, only : rprec, cprec
      USE stel_constants, only : pi, twopi, one, zero
     
!-------------------------------------------------------------------------------
!  Use Statements for other structures, V3 Utilities
!-------------------------------------------------------------------------------

      USE v3_utilities
      USE diagnostic_T
      USE geometric_T
      USE coosig_T

!-------------------------------------------------------------------------------
!  Implicit None comes after USE statements, before other declarations
!-------------------------------------------------------------------------------
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Make type declarations and constants Private, so there are no conflicts.
!-------------------------------------------------------------------------------
      PRIVATE rprec, cprec, pi, twopi, one, zero

!-------------------------------------------------------------------------------
!  Lengths of Character Variables
!-------------------------------------------------------------------------------
      INTEGER, PARAMETER, PRIVATE :: type_len=10      
      INTEGER, PARAMETER, PRIVATE :: sn_len=30      
      INTEGER, PARAMETER, PRIVATE :: ln_len=80      
      INTEGER, PARAMETER, PRIVATE :: units_len=30      

!*******************************************************************************
! SECTION II. DERIVED-TYPE DECLARATIONS
! 1)   Signal Description:
!         signal_desc  
!     Type of signal specified by  % s_type.
!     Allowable values of s_type:
!       diagnostic     -  implemented in module diagnostic_T
!                             Diagnostic Types are: mddc, ipsl, sxrch
!       geometric      - implemented in module geometric_T
!                             Only geometric type at this time is edge_limit
!       coosig [Combination Of Other SIGnals] - implemented in module coosig_T
!
!  A signal_desc is supposed to hold all the information necessary to describe
!  the signal. For an s_type='diagnostic' signal, this is mostly a pointer to a 
!  diagnostic_desc. For an s_type='geometric' signal, this is mostly a pointer 
!  to a geometric_desc. And for an s_type='coosig' signal, this is mostly a
!  pointer to a coosig_desc
!
! 2)   Signal Data:
!         signal_data
!     Type of signal data specified by  % sd_type.
!     Allowable values of sd_type:
!       'observe'        - based on diagnostic data, observation
!       'model'          - computed from the model
!
!  A signal_data is supposed to hold the computed signal. 
!  For an sd_type='observe', the signal_data is computed from diagnostic
!  observations. (And in the case of magnetic diagnostics, the signal computation
!  is a pretty trivial function of the magnetic diagnostic data). 
!  For an sd_type='model', the signal_data is computed from the signal
!   description and the equilibrium state (model).
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Declare type signal_desc
!    Common to all s_types                                                     
!       s_type          character, type of signal
!       s_name          character, short name of signal
!       l_name          character, long name of signal
!       units           character, physical units that the signal is measured in
!       diag            type diagnostic_desc pointer. Description of diagnostic.
!         ^ USE AS POINTER
!           Used for s_type = 'diagnostic'
!       geom            type geometric_desc pointer. Description of geometric.
!         ^ USE AS POINTER
!         Used for s_type = 'geometric'
!       coosig            type coosig_desc pointer. Description of a combination
!                         of isgnals.
!         ^ USE AS POINTER
!         Used for s_type = 'coosig'
!-------------------------------------------------------------------------------
      TYPE signal_desc
         CHARACTER (len=type_len)         :: s_type
         CHARACTER (len=sn_len)           :: s_name                                 
         CHARACTER (len=ln_len)           :: l_name
         CHARACTER (len=units_len)        :: units                                 
         TYPE (diagnostic_desc), POINTER  :: diag => null()
         TYPE (geometric_desc), POINTER   :: geom => null()
         TYPE (coosig_desc), POINTER      :: coosig => null()
       END TYPE signal_desc
!-------------------------------------------------------------------------------
!  Declare type signal_data
!    sd_type          character, type of signal data
!                        Either 'model' or 'observe'
!    desc             pointer to signal_desc structure
!      ^ USE AS POINTER
!    data             real array
!      ^ USE AS ALLOCATABLE ARRAY
!    sigma            real array, should be same length as data
!      ^ USE AS ALLOCATABLE ARRAY
!-------------------------------------------------------------------------------
      TYPE signal_data
         CHARACTER (len=type_len)           :: sd_type
         TYPE (signal_desc), POINTER        :: desc => null()
         REAL(rprec), DIMENSION(:), POINTER :: data => null()
         REAL(rprec), DIMENSION(:), POINTER :: sigma => null()
         REAL(rprec)                        :: weight
      END TYPE signal_data

!*******************************************************************************
! SECTION III. INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Assignment for structures
!-------------------------------------------------------------------------------
      INTERFACE ASSIGNMENT (=)
         MODULE PROCEDURE signal_desc_assign,                                  &
     &                    signal_data_assign, signal_data_assign_a
      END INTERFACE

!-------------------------------------------------------------------------------
!  Generic construct
!-------------------------------------------------------------------------------
      INTERFACE signal_construct
         MODULE PROCEDURE signal_desc_construct_d,                             &
     &                    signal_desc_construct_g,                             &
     &                    signal_desc_construct_c,                             &
     &                    signal_data_construct,                               &
     &                    signal_data_construct_sds
         END INTERFACE

!-------------------------------------------------------------------------------
!  Generic destroy
!-------------------------------------------------------------------------------
      INTERFACE signal_destroy
         MODULE PROCEDURE signal_desc_destroy,                                 &
     &                    signal_data_destroy
      END INTERFACE

!-------------------------------------------------------------------------------
!  Generic write
!-------------------------------------------------------------------------------
      INTERFACE signal_write
         MODULE PROCEDURE signal_desc_write,                                   &
     &                    signal_data_write
      END INTERFACE

!-------------------------------------------------------------------------------
!  Interface block for testing goes here. 
!-------------------------------------------------------------------------------

      CONTAINS
!*******************************************************************************
! SECTION IV. CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Construct a signal_desc
!
!  For s_type = 'diagnostic' 
!-------------------------------------------------------------------------------
      SUBROUTINE signal_desc_construct_d(this,s_type,s_name,l_name,            &
     &   units,diagnostic)

      IMPLICIT NONE

!  Declare Arguments 
      TYPE (signal_desc), INTENT(inout)          :: this
      CHARACTER (len=*), INTENT(in)              :: s_type
      CHARACTER (len=*), INTENT(in)              :: s_name
      CHARACTER (len=*), INTENT(in)              :: l_name
      CHARACTER (len=*), INTENT(in)              :: units
      TYPE (diagnostic_desc), INTENT(in), TARGET :: diagnostic

!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'signal_desc_construct_d: '

!  Start of executable code

!  Assignments for all s_type's
      this % s_name = TRIM(ADJUSTL(s_name))
      this % l_name = TRIM(ADJUSTL(l_name))
      this % units = TRIM(ADJUSTL(units))
      
!  This subroutine assumes s_type is diagnostic
      this % s_type = 'diagnostic'
      this % diag => diagnostic
      
      IF (TRIM(ADJUSTL(s_type)) .ne. 'diagnostic') THEN
         CALL err_warn(sub_name // 'expected s_type "diagnostic": ',           &
     &      char=s_type)
      ENDIF
      
      END SUBROUTINE signal_desc_construct_d
!-------------------------------------------------------------------------------
!  Construct a signal_desc
!
!  For s_type = 'geometric' 
!-------------------------------------------------------------------------------
      SUBROUTINE signal_desc_construct_g(this,s_type,s_name,l_name,            &
     &   units,geometric)

      IMPLICIT NONE

!  Declare Arguments 
      TYPE (signal_desc), INTENT(inout)          :: this
      CHARACTER (len=*), INTENT(in)              :: s_type
      CHARACTER (len=*), INTENT(in)              :: s_name
      CHARACTER (len=*), INTENT(in)              :: l_name
      CHARACTER (len=*), INTENT(in)              :: units
      TYPE (geometric_desc), INTENT(in), TARGET  :: geometric

!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'signal_desc_construct_g: '

!  Start of executable code

!  Assignments for all s_type's
      this % s_name = TRIM(ADJUSTL(s_name))
      this % l_name = TRIM(ADJUSTL(l_name))
      this % units = TRIM(ADJUSTL(units))
      
!  This subroutine assumes s_type is diagnostic
      this % s_type = 'geometric'
      this % geom => geometric
      
      IF (TRIM(ADJUSTL(s_type)) .ne. 'geometric') THEN
         CALL err_warn(sub_name // 'expected s_type "geometric": ',           &
     &      char=s_type)
      ENDIF
      
      END SUBROUTINE signal_desc_construct_g
!-------------------------------------------------------------------------------
!  Construct a signal_desc
!
!  For s_type = 'coosig' 
!-------------------------------------------------------------------------------
      SUBROUTINE signal_desc_construct_c(this,s_type,s_name,l_name,            &
     &   units,coosig)

      IMPLICIT NONE

!  Declare Arguments 
      TYPE (signal_desc), INTENT(inout)          :: this
      CHARACTER (len=*), INTENT(in)              :: s_type
      CHARACTER (len=*), INTENT(in)              :: s_name
      CHARACTER (len=*), INTENT(in)              :: l_name
      CHARACTER (len=*), INTENT(in)              :: units
      TYPE (coosig_desc), INTENT(in), TARGET     :: coosig

!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'signal_desc_construct_c: '

!  Start of executable code

!  Assignments for all s_type's
      this % s_name = TRIM(ADJUSTL(s_name))
      this % l_name = TRIM(ADJUSTL(l_name))
      this % units = TRIM(ADJUSTL(units))
      
!  This subroutine assumes s_type is coosig
      this % s_type = 'coosig'
      this % coosig => coosig
      
      IF (TRIM(ADJUSTL(s_type)) .ne. 'coosig') THEN
         CALL err_warn(sub_name // 'expected s_type "coosig": ',           &
     &      char=s_type)
      ENDIF
      
      END SUBROUTINE signal_desc_construct_c

!-------------------------------------------------------------------------------
!  Construct a signal_data, array of data
!-------------------------------------------------------------------------------
      SUBROUTINE signal_data_construct(this,desc,sd_type,data,sigma,           &
     &   weight)

      IMPLICIT NONE

!  Declare Arguments 
      TYPE (signal_data), INTENT(inout)                    :: this
      TYPE (signal_desc), TARGET, INTENT(in)               :: desc
      CHARACTER(len=*), INTENT(in)                         :: sd_type
      REAL(rprec), DIMENSION(:), INTENT(in)                :: data
      REAL(rprec), DIMENSION(:), INTENT(in)                :: sigma
      REAL(rprec), OPTIONAL, INTENT(in)                    :: weight

!  Declare local variables
      INTEGER              :: ier1, ier2, n_data, n_sigma
      CHARACTER(len=type_len)     :: temp_sd_type
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'signal_data_construct: '

!  Start of executable code

!  Destroy existing arrays
      CALL signal_data_destroy(this)
            
!  Pointer assignment
      this % desc => desc

!  Different coding, depending on sd_type
!   Case for model includes e, E, for compatibility with old code, where
!   'equilibrium' was used instead of 'model'
      temp_sd_type = ADJUSTL(sd_type)
      SELECT CASE (temp_sd_type(1:1))
      CASE ('m','M','e','E')
         this % sd_type = 'model'
      CASE ('o','O')
         this % sd_type = 'observe'

      CASE DEFAULT
         CALL err_fatal(sub_name // 'unrecognized sd_type: ',                  &
     &      char=sd_type)
      END SELECT ! Different coding depending on sd_type

!  Take length of data from data argument
      n_data = SIZE(data)
      n_sigma = SIZE(sigma)

!  Error exit if the size of the data is not the same as the size of the sigma.
      CALL assert_eq(n_data,n_sigma,sub_name // 'n_data & n_sigma')
      
      ALLOCATE(this % data(1:n_data),STAT=ier1)
      ALLOCATE(this % sigma(1:n_data),STAT=ier2)
      CALL assert_eq(0,ier1,ier2,sub_name // 'alloc')

      this % data = data
      this % sigma = sigma
      
      IF(PRESENT(weight)) THEN
         this % weight = weight
      ELSE
         this % weight = one
      ENDIF
      
      END SUBROUTINE signal_data_construct

!-------------------------------------------------------------------------------
!  Construct a signal_data, scalar data, sigma (sds)
!  This also takes care of the case where the data and/or sigma are missing
!  When data and/or sigma are missing, it is probably a model_compute
!  subroutine that will allocate and fill in the data.
!  (Could write a signal_data_change_data subroutine - JDH 2011-07-07)
!-------------------------------------------------------------------------------
      SUBROUTINE signal_data_construct_sds(this,desc,sd_type,data,sigma,       &
     &   weight)

      IMPLICIT NONE

!  Declare Arguments 
      TYPE (signal_data), INTENT(inout)                    :: this
      TYPE (signal_desc), TARGET, INTENT(in)               :: desc
      CHARACTER(len=*), INTENT(in)                         :: sd_type
      REAL(rprec), INTENT(in), OPTIONAL                    :: data
      REAL(rprec), INTENT(in), OPTIONAL                    :: sigma
      REAL(rprec), OPTIONAL, INTENT(in)                    :: weight

!  Declare local variables
      REAL(rprec), DIMENSION(1)      :: array_data
      REAL(rprec), DIMENSION(1)      :: array_sigma

!  Start of executable code
      IF(PRESENT(data)) THEN
         array_data = data
      ELSE
         array_data = zero
      ENDIF
      IF(PRESENT(sigma)) THEN
         array_sigma = sigma
      ELSE
         array_sigma = zero
      ENDIF
      IF(PRESENT(weight)) THEN
         CALL signal_data_construct(this,desc,sd_type,array_data,              &
     &   array_sigma,weight)
      ELSE
         CALL signal_data_construct(this,desc,sd_type,array_data,              &
     &      array_sigma)
      ENDIF

      END SUBROUTINE signal_data_construct_sds
      
!*******************************************************************************
! SECTION V. DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Destroy a signal_desc
!-------------------------------------------------------------------------------
      SUBROUTINE signal_desc_destroy(this)
      IMPLICIT NONE

!  Declare Arguments 
      TYPE (signal_desc), INTENT(inout) :: this

!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'signal_desc_destroy: '

!  Start of executable code

!  Components common to all s_type's
      this % s_name = ' '
      this % l_name = ' '
      this % units = ' '

!  Different coding, depending on s_type
      SELECT CASE (TRIM(ADJUSTL(this % s_type)))
      CASE ('diagnostic')
         this % s_type = ' '
         this % diag => null()
!  ? Is this the behaviour desired, where the pointer to the diagnostic_desc
!  is merely nulled, and the target diagnostic_desc is untouched ?
!  This could lead to a memory leak. 2007-06-12
!  Same consideration for geometric. 2009-01-21
      CASE ('geometric')
         this % s_type = ' '
         this % geom => null()
!  Same consideration for coosig. 2011-07-06
      CASE ('coosig')
         this % s_type = ' '
         this % coosig => null()

      CASE DEFAULT
         CALL err_fatal(sub_name // 'unrecognized s_type: ',                   &
     &      char=this % s_type)
      END SELECT ! Different coding depending on s_type

      END SUBROUTINE signal_desc_destroy

!-------------------------------------------------------------------------------
!  Destroy a signal_data
!-------------------------------------------------------------------------------
      SUBROUTINE signal_data_destroy(this)
      IMPLICIT NONE

!  Declare Arguments 
      TYPE (signal_data), INTENT(inout) :: this

!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'signal_data_destroy: '
      INTEGER :: ier1, ier2

!  Start of executable code
!  Initialize STAT variables
      ier1 = 0
      ier2 = 0

!  Get rid of all components
      this % desc => null()
      this % sd_type = ' '
      IF (ASSOCIATED(this % data)) DEALLOCATE(this % data,STAT=ier1)
      IF (ASSOCIATED(this % sigma)) DEALLOCATE(this % sigma,STAT=ier2)
      CALL assert_eq(0,ier1,ier2,sub_name // 'deallocate error')
      this % weight = zero

      END SUBROUTINE signal_data_destroy

!*******************************************************************************
! SECTION VI. ASSIGNMENT SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Assignment for signal_desc
!-------------------------------------------------------------------------------
      SUBROUTINE signal_desc_assign(left,right)

!  2007-06-12. Would this be the same as generic assignment?

      IMPLICIT NONE
      
!  Declare Arguments 
      TYPE (signal_desc), INTENT (inout) :: left
      TYPE (signal_desc), INTENT (in) :: right
      
!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'signal_desc_assign: '
         
!  Start of executable code

!  Pointer Assignment for diagnostic_desc, geometric_desc, and coosig_desc
      left % diag => right % diag
      left % geom => right % geom
      left % coosig => right % coosig

      left % s_type = right % s_type
      left % s_name = right % s_name
      left % l_name = right % l_name
      left % units = right % units
         
      END SUBROUTINE signal_desc_assign

!-------------------------------------------------------------------------------
!  Assignment for signal_data
!-------------------------------------------------------------------------------
      SUBROUTINE signal_data_assign(left,right)

      IMPLICIT NONE
      
!  Declare Arguments 
      TYPE (signal_data), INTENT (inout) :: left
      TYPE (signal_data), INTENT (in) :: right
      
!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'signal_data_assign: '
      CHARACTER (len=*), PARAMETER :: err_mess1 =                              &
     & 'left-right data or sigma are the same?. FIX IT'
      INTEGER :: ier1, ier2
         
!  Start of executable code

!  Check to see if the 'use as allocatable array' pointers are pointing
!  to the same location
      CALL assert(.not.ASSOCIATED(left % data,right % data),                   &
     &            .not.ASSOCIATED(left % sigma,right % sigma),                 &
     &   sub_name // err_mess1)

!  Destroy left
      CALL signal_data_destroy(left)
      
      left % sd_type = right % sd_type

!  Pointer Assignment for desc
      left % desc => right % desc

!  Allocate space for the 'use as allocatable array' pointer variables
      ALLOCATE(left % data(1:SIZE(right % data)),STAT=ier1)
      ALLOCATE(left % sigma(1:SIZE(right % sigma)),STAT=ier2)
      CALL assert_eq(0,ier1,ier2,sub_name // 'Allocation error')

      left % data = right % data
      left % sigma = right % sigma
      left % weight = right % weight
         
      END SUBROUTINE signal_data_assign

!-------------------------------------------------------------------------------
!  Assignment for array of signal_data
!-------------------------------------------------------------------------------
      SUBROUTINE signal_data_assign_a(left,right)

      IMPLICIT NONE
      
!  Declare Arguments 
      TYPE (signal_data), DIMENSION(:), INTENT (inout) :: left
      TYPE (signal_data), DIMENSION(:), INTENT (in) :: right
      
!  Declare local variables
      INTEGER :: n_left, n_right, i
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'signal_data_assign_a: '
      CHARACTER (len=*), PARAMETER :: err_mess1 =                              &
     & 'left-right array lengths are not the same.'
         
!  Start of executable code

      n_left = SIZE(left)
      n_right = SIZE(right)
      CALL assert_eq(n_left,n_right,sub_name // err_mess1)
      DO i = 1,n_left
         left(i) = right(i)
      END DO
         
      END SUBROUTINE signal_data_assign_a

!*******************************************************************************
! SECTION VII.  OUTPUT SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Write out the contents of a signal_desc
!-------------------------------------------------------------------------------

      SUBROUTINE signal_desc_write(this,identifier,unit,verbose)
      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_desc), INTENT (in) :: this
      CHARACTER (len=*), INTENT(in), OPTIONAL :: identifier
      INTEGER, INTENT(in), OPTIONAL :: unit
      INTEGER, INTENT(in), OPTIONAL :: verbose
!  identifier   character variable, also written out
!  unit         I/O unit number to write to
!  verbose      integer, to specify verbosity level of write

!  Declare local variables and constants
      INTEGER :: iv_default = 1
      INTEGER :: iv
      INTEGER :: iou_default = 6
      INTEGER :: iou
      CHARACTER (len=60) :: id

!  Declare Format array
      CHARACTER(len=*), PARAMETER, DIMENSION(9) :: fmt1 = (/                   &
     & '(" start signal_desc write, called with id = ",a)    ',                &
     & '(" s_type = ",a)                                     ',                &
     & '(" s_name = ",a)                                     ',                &
     & '(" l_name = ",a)                                     ',                &
     & '(" units = ",a)                                      ',                &
     & '(" diagnostic_desc s_name = ",a)                     ',                &
     & '(" geometric_desc name = ",a)                        ',                &
     & '(" coosig_desc name = ",a)                           ',                &
     & '(" end signal_desc write, called with id = ",a)      '                 &
     &  /) 

!  start of executable code
!  Check for arguments present
      IF (PRESENT(identifier)) THEN
         id = identifier
      ELSE
         id = ' '
      END IF

      IF (PRESENT(unit)) THEN
         iou = unit
      ELSE
         iou = iou_default
      END IF

      IF (PRESENT(verbose)) THEN
         iv = verbose
      ELSE
         iv = iv_default
      END IF
      
!  Select Case of Verbosity Level
!  Will need Select Case on s_type
      SELECT CASE(iv)
      CASE( :0)  ! VERY Terse
         WRITE(iou,*) this % s_type
         WRITE(iou,*) this % s_name
         WRITE(iou,*) this % l_name
         WRITE(iou,*) this % units
         SELECT CASE (this % s_type)
         CASE('diagnostic')
            WRITE(iou,*) this % diag % s_name
         CASE('geometric')
            WRITE(iou,*) this % geom % name
         CASE('coosig')
            WRITE(iou,*) this % coosig % name
         CASE DEFAULT
            WRITE(iou,*) 'bad this % s_type = ', this % s_type
         END SELECT
      
      CASE(1:)    ! Default, more verbose
         WRITE(iou,fmt1(1)) id
         WRITE(iou,fmt1(2)) this % s_type
         WRITE(iou,fmt1(3)) this % s_name
         WRITE(iou,fmt1(4)) this % l_name
         WRITE(iou,fmt1(5)) this % units
         SELECT CASE (this % s_type)
         CASE('diagnostic')
            WRITE(iou,fmt1(6)) this % diag % s_name
         CASE('geometric')
            WRITE(iou,fmt1(7)) this % geom % name
         CASE('coosig')
            WRITE(iou,fmt1(8)) this % coosig % name
         CASE DEFAULT
            WRITE(iou,*) 'bad this % s_type = ', this % s_type
         END SELECT
         WRITE(iou,fmt1(9)) id
      
      END SELECT

      END SUBROUTINE signal_desc_write

!-------------------------------------------------------------------------------
!  Write out the contents of a signal_data
!-------------------------------------------------------------------------------

      SUBROUTINE signal_data_write(this,identifier,unit,verbose)
      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_data), INTENT (in) :: this
      CHARACTER (len=*), INTENT(in), OPTIONAL :: identifier
      INTEGER, INTENT(in), OPTIONAL :: unit
      INTEGER, INTENT(in), OPTIONAL :: verbose
!  identifier   character variable, also written out
!  unit         I/O unit number to write to
!  verbose      integer, to specify verbosity level of write

!  Declare local variables and constants
      INTEGER      :: iv_default = 1
      INTEGER      :: iv
      INTEGER      :: iou_default = 6
      INTEGER      :: iou
      CHARACTER (len=60)  :: id
      INTEGER      :: i, n_data

!  Declare Format array
      CHARACTER(len=*), PARAMETER, DIMENSION(5) :: fmt1 = (/                   &
     & '(/" start signal_data write, called with id = ",a)           ',        &
     & '(" signal_desc s_name = ",a,"| signal_data sd_type = ",a)    ',        &
     & '(" number of data points = ",i4, " weight = ",es10.3)        ',        &
     & '(" index   data:    sigma: ",/,(1x,i4,2(3x,es12.5)))         ',        &
     & '(" end signal_data write, called with id = ",a)              '         &
     &  /) 

!  start of executable code
!  Check for arguments present
      IF (PRESENT(identifier)) THEN
         id = identifier
      ELSE
         id = ' '
      END IF

      IF (PRESENT(unit)) THEN
         iou = unit
      ELSE
         iou = iou_default
      END IF

      IF (PRESENT(verbose)) THEN
         iv = verbose
      ELSE
         iv = iv_default
      END IF

!  Select Case of Verbosity Level
      SELECT CASE(iv)
      CASE( :0)  ! VERY Terse
         WRITE(iou,*) this % desc % s_name
         WRITE(iou,*) this % sd_type
         WRITE(iou,*) this % data
         WRITE(iou,*) this % sigma
         WRITE(iou,*) this % weight
      
      CASE(1:)    ! Default, more verbose
         n_data = SIZE(this % data)
         WRITE(iou,fmt1(1)) id
         WRITE(iou,fmt1(2)) this % desc % s_name, this % sd_type 
         WRITE(iou,fmt1(3)) n_data, this % weight
         WRITE(iou,fmt1(4)) (i,this % data(i),this % sigma(i),                 &
     &      i=1,n_data)
         WRITE(iou,fmt1(5)) id
      
      END SELECT

      END SUBROUTINE signal_data_write

!-------------------------------------------------------------------------------
!  Write out the model and observed signal data (arrays)
!-------------------------------------------------------------------------------

      SUBROUTINE signal_data_moa_write(sdm,sdo,identifier,unit,verbose)
      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_data), DIMENSION(:), INTENT (in) :: sdm, sdo
      CHARACTER (len=*), INTENT(in), OPTIONAL :: identifier
      INTEGER, INTENT(in), OPTIONAL :: unit
      INTEGER, INTENT(in), OPTIONAL :: verbose
!  identifier   character variable, also written out
!  unit         I/O unit number to write to
!  verbose      integer, to specify verbosity level of write

!  Declare local variables and constants
      INTEGER      :: iv_default = 1
      INTEGER      :: iv
      INTEGER      :: iou_default = 6
      INTEGER      :: iou
      CHARACTER (len=60)  :: id
      INTEGER      :: i, n_data, nsdm, nsdo, n_data_model, icount_geom
      INTEGER      :: it, nterms
      REAL(rprec)  :: ssigma, evec, sweight
      CHARACTER (len=type_len)     :: dgc_type
      CHARACTER(len=*), PARAMETER  :: sub_name =                               &
     &  'signal_data_moa_write: '

!  Declare Format array ||| 2011-01-21 Changed to format statements.

! ------------------------------------------------------------------------------
!  Start of Executable Code

!  Check for arguments present
      IF (PRESENT(identifier)) THEN
         id = identifier
      ELSE
         id = ' '
      END IF

      IF (PRESENT(unit)) THEN
         iou = unit
      ELSE
         iou = iou_default
      END IF

      IF (PRESENT(verbose)) THEN
         iv = verbose
      ELSE
         iv = iv_default
      END IF

!  Check that lengths of sdm and sdo are the same
      nsdm = SIZE(sdm)
      nsdo = SIZE(sdo)
      CALL assert_eq(nsdm,nsdo,sub_name // 'sdm, sdo unequal lengths')

! For now, ignore verbosity level
! Write header information
      WRITE(iou,1000) id
      WRITE(iou,'(/a)') '[diagnostic mddc - Model(2) is plasma signal]'
      WRITE(iou,2000) nsdm
      WRITE(iou,3000) 

!  Loop over arrays
      icount_geom = 0
      DO i = 1,nsdm
         CALL assert(sdm(i) % sd_type .eq. 'model',                            &
     &      sub_name // 'sd_type bad - model')
         CALL assert(sdo(i) % sd_type .eq. 'observe',                          &
     &      sub_name // 'sd_type bad - observe')
         CALL assert(sdm(i) % desc%s_name .eq. sdo(i) % desc%s_name,           &
     &      sub_name // 's_name disagreement')
         n_data_model = SIZE(sdm(i) % data)
!  ssigma coding copied from subroutine gsq_evaluate_g2
!  Problem with zero ssigma should have been detected earlier.
         ssigma = SQRT(sdm(i) % sigma(1) ** 2 + sdo(i) % sigma(1) ** 2)
         sweight = sdo(i) % weight
!  Weight comes from sdo
         evec = SQRT(sweight) * (sdo(i) % data(1) - sdm(i) % data(1))          &
     &      / ssigma
         IF (sdm(i) % desc % s_type .eq. 'diagnostic') THEN
            dgc_type = sdm(i) % desc % diag % d_type
         ELSEIF (sdm(i) % desc % s_type .eq. 'geometric') THEN
            dgc_type = sdm(i) % desc % geom % g_type
            icount_geom = icount_geom + 1
         ELSEIF (sdm(i) % desc % s_type .eq. 'coosig') THEN
            dgc_type = sdm(i) % desc % coosig % comb_type
         ELSE
            dgc_type = 'xxxx'
         ENDIF
         IF (n_data_model .ge. 2) THEN
            WRITE(iou,4000) i, sdm(i) % desc % s_type, dgc_type,               &
     &         sdm(i) %desc % s_name, sdo(i) % data(1),                        &
     &         sdm(i) % data(1), ssigma, sweight, evec,                        &
     &         sdm(i) % data(2)
         ELSE
            WRITE(iou,4000) i, sdm(i) % desc % s_type, dgc_type,               &
     &         sdm(i) %desc % s_name, sdo(i) % data(1),                        &
     &         sdm(i) % data(1), ssigma, sweight, evec
         ENDIF
      END DO

!  Extra information about Geometric Signals
      IF (icount_geom .gt. 0) THEN
         WRITE(iou,4100)
         DO i = 1,nsdm
            IF (sdm(i) % desc % s_type .eq. 'geometric') THEN
               WRITE(iou,4200) i, sdm(i) % desc % s_name,                      &
     &            sdm(i) % data(2), sdm(i) % data(3), sdm(i) % data(4)
            ELSE
               CYCLE
            ENDIF
         END DO
      ENDIF

!  Extra information about COOSIG Signals
      DO i = 1,nsdm
         IF (sdm(i) % desc % s_type .eq. 'coosig') THEN
            WRITE(iou,4300) i, sdm(i) % desc % coosig % comb_type
            nterms = SIZE(sdm(i) % desc % coosig % index)
            DO it = 1,nterms
               WRITE(iou,4400) it, sdm(i) % desc % coosig % index(it),         &
     &            sdm(i) % desc % coosig % a(it)
            END DO
         ENDIF
      END DO

      WRITE(iou,5000) id

1000  FORMAT(/" *** START signal_data_moa_write, called with id = ",a)
2000  FORMAT(i4," Signals"/)
3000  FORMAT(3x,"#",t15,"type",t30,"s_name",t54,"Observe",t68,"Model",         &
     &   t81,"Sigma",t94,"weight",t106,"evec",t118,"Model(2)")
4000  FORMAT(i4,2x,2(a10,1x),a20,2(2x,es12.5),3(2x,es10.3),2x,es12.5)
4100  FORMAT(/"Geometrics - Position on s=1 surface of Max iso_fun",/,
     &  3x,"#",t15,"s_name",t30,"R",t44,"Phi",t58,"Z")
4200  FORMAT(i4,10x,1(a10,1x),3(2x,es12.5))
4300  FORMAT(/"Signal ",i4," is a Combination Of Other SIGnals "               &
     & "(coosig) of type '",a,"'"/                                             &
     &   "     term #    Other signal #        Coefficient")
4400  FORMAT(t5,i4,t22,i4,t36,es12.5)
5000  FORMAT(/" *** END signal_data_moa_write, called with id = ",a) 

      END SUBROUTINE signal_data_moa_write

!-------------------------------------------------------------------------------
!  Write out the model signal data (array)
!-------------------------------------------------------------------------------

      SUBROUTINE signal_data_ma_write(sdm,identifier,unit,verbose)

!  JDH First version 2011-03-15. Based on signal_data_moa_write
      IMPLICIT NONE

!  Declare Arguments
      TYPE (signal_data), DIMENSION(:), INTENT (in) :: sdm
      CHARACTER (len=*), INTENT(in), OPTIONAL :: identifier
      INTEGER, INTENT(in), OPTIONAL :: unit
      INTEGER, INTENT(in), OPTIONAL :: verbose
!  identifier   character variable, also written out
!  unit         I/O unit number to write to
!  verbose      integer, to specify verbosity level of write

!  Declare local variables and constants
      INTEGER      :: iv_default = 1
      INTEGER      :: iv
      INTEGER      :: iou_default = 6
      INTEGER      :: iou
      CHARACTER (len=60)  :: id
      INTEGER      :: i, n_data, nsdm,n_data_model, icount_geom
      INTEGER      :: it, nterms
      REAL(rprec)  :: ssigma, evec, sweight
      CHARACTER (len=type_len)     :: dgc_type
      CHARACTER(len=*), PARAMETER  :: sub_name =                               &
     &  'signal_data_ma_write: '

! ------------------------------------------------------------------------------
!  Start of Executable Code

      nsdm = SIZE(sdm)

!  Check for arguments present
      IF (PRESENT(identifier)) THEN
         id = identifier
      ELSE
         id = ' '
      END IF

      IF (PRESENT(unit)) THEN
         iou = unit
      ELSE
         iou = iou_default
      END IF

      IF (PRESENT(verbose)) THEN
         iv = verbose
      ELSE
         iv = iv_default
      END IF

! For now, ignore verbosity level
! Write header information
      WRITE(iou,1000) id
      WRITE(iou,'(/a)') '[diagnostic mddc - Model(2) is plasma signal]'
      WRITE(iou,2000) nsdm
      WRITE(iou,3000) 

!  Loop over arrays
      icount_geom = 0
      DO i = 1,nsdm
         CALL assert(sdm(i) % sd_type .eq. 'model',                            &
     &      sub_name // 'sd_type bad - model')
         n_data_model = SIZE(sdm(i) % data)
         IF (sdm(i) % desc % s_type .eq. 'diagnostic') THEN
            dgc_type = sdm(i) % desc % diag % d_type
         ELSEIF (sdm(i) % desc % s_type .eq. 'geometric') THEN
            dgc_type = sdm(i) % desc % geom % g_type
            icount_geom = icount_geom + 1
         ELSEIF (sdm(i) % desc % s_type .eq. 'coosig') THEN
            dgc_type = sdm(i) % desc % coosig % comb_type
         ELSE
            dgc_type = 'xxxx'
         ENDIF
         IF (n_data_model .ge. 2) THEN
            WRITE(iou,4000) i, sdm(i) % desc % s_type, dgc_type,                &
     &         sdm(i) %desc % s_name, sdm(i) % data(1),                        &
     &         sdm(i) % data(2)
         ELSE
            WRITE(iou,4000) i, sdm(i) % desc % s_type, dgc_type,                &
     &         sdm(i) %desc % s_name, sdm(i) % data(1)
         ENDIF
      END DO

!  Extra information about Geometric Signals
      IF (icount_geom .gt. 0) THEN
         WRITE(iou,4100)
         DO i = 1,nsdm
            IF (sdm(i) % desc % s_type .eq. 'geometric') THEN
               WRITE(iou,4200) i, sdm(i) % desc % s_name,                      &
     &            sdm(i) % data(2), sdm(i) % data(3), sdm(i) % data(4)
            ELSE
               CYCLE
            ENDIF
         END DO
      ENDIF

!  Extra information about COOSIG Signals
      DO i = 1,nsdm
         IF (sdm(i) % desc % s_type .eq. 'coosig') THEN
            WRITE(iou,4300) i, sdm(i) % desc % coosig % comb_type
            nterms = SIZE(sdm(i) % desc % coosig % index)
            DO it = 1,nterms
               WRITE(iou,4400) it, sdm(i) % desc % coosig % index(it),         &
     &            sdm(i) % desc % coosig % a(it)
            END DO
         ENDIF
      END DO
 
      WRITE(iou,5000) id

1000  FORMAT(/" *** START signal_data_ma_write, called with id = ",a)
2000  FORMAT(i4," Signals"/)
3000  FORMAT(3x,"#",t15,"type",t30,"s_name",t54,"Model",t68,"Model(2)")
4000  FORMAT(i4,2x,2(a10,1x),a20,2(2x,es12.5))
4100  FORMAT(/'Geometrics - Position on s=1 surface of Max iso_fun',/,
     &  3x,"#",t15,"s_name",t30,"R",t44,"Phi",t58,"Z")
4200  FORMAT(i4,9x,1(a10,1x),3(2x,es12.5))
4300  FORMAT(/"Signal ",i4," is a Combination Of Other SIGnals, type ",        &
     &   a/"     term #    Other signal #      Coefficient")
4400  FORMAT(5x,i4,10x,i4,10x,es12.5)
5000  FORMAT(/" *** END signal_data_ma_write, called with id = ",a) 

      END SUBROUTINE signal_data_ma_write

!*******************************************************************************
! SECTION XVI.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
!  
!  JDH 08-23-04. 
!     Modifying diagnostic_T.f to get signal_T.f

!  JDH 08-24-2004
!     First thorough cleanup of signals
!     Needs precomp coding
!
!  JDH 08-26-2004
!     Added signal_mrf coding
!
!  JDH 09-02-2004
!     Added flux_factor to signal_desc. Added code_name, code_verson, date_run, 
!     field_coils_id, and extcur_mg to signal_mrf.
!
!  JDH 09-10-2004
!     Deleted flux_factor from signal_desc. (Flux Factor in diagnostic_desc is
!     good enough). 
!
!  JDH 12-02-2004
!     Fixed move of extcur_mg in signal_mrf_construct
!
!  JDH 12-11-2004
!     Eliminated 'pointer' attribute of mrf component of signal_desc. Removed
!     desc component of signal_mrf. Added l_mrf_def component to signal_desc.
!     Added subroutine signal_desc_assign.
!
!  JDH 07-04-2005
!     Added subroutine signal_mrf_write. Fixed sub_name in signal_mrf_assign.
!
!  JDH 09-18-2006
!     Modified SUBROUTINE signal_data_write
!
!  JDH 09-20-2006
!     Removed references to sd_type 'equilibrium'. Now use 'model'.
!
!  JDH 2007-06-12
!     Eliminated signal_mrf. Changed only s_type to 'diagnostic'.
!
!  JDH 2007-07-05
!     Added subroutine for assignment of a 1-d array of signal_data
!     Fortran 2003 feature allowing ALLOCATABLE components would be NICE
!
!  JDH 2008-01-19 - Added => null() to pointer declarations
!
!  JDH 2008-01-21
!    SPH Jan 2008 eliminated iprec. Completed elimination.
!    Initialized STAT variables in signal_data_destroy
!
!  JDH 2009-01-21
!    Added s_type = 'geometric'
!
!  JDH 2010-06-10
!    Added subroutine signal_data_moa_write
!
!  JDH 2011-01-18
!    Added weight component to signal_data type. 
!
!  JDH 2011-03-16
!    Added subroutine to write a model data array, signal_data_ma_write
!
!  JDH 2011-07-06
!    Added s_type = 'coosig', Combination Of Other SIGnals
!    Changed name from signal_data_construct_s to signal_data_construct_sds, to
!    emphasize that it is the data and sigma that is scalar. Made arguments
!    data and sigma optional (Easier to use for model data creation)
!
!  JDH 2011-07-11
!    Added printout  of extra information on geometrics and coosig in
!    subroutines signal_data_moa_write and signal_data_ma_write

      END MODULE signal_T
