!*******************************************************************************
! Combination Of Other SIGnals - COOSIG
!  File coosig_T.f
!  Contains module coosig_T
!  Defines derived-types: coosig_desc

!*******************************************************************************
!  MODULE coosig_T
!    (coosig Type Definition, for the V3FIT code)
! SECTION I.    VARIABLE DECLARATIONS
! SECTION II.   DERIVED-TYPE DECLARATIONS
! SECTION III.  INTERFACE BLOCKS
! SECTION IV.   CONSTRUCTION SUBROUTINES
! SECTION V.    DESTRUCTION SUBROUTINES
! SECTION VI.   ASSIGNMENT SUBROUTINES
! SECTION VII.  OUTPUT SUBROUTINES

! SECTION XVI.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
      MODULE coosig_T

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
      INTEGER, PARAMETER, PRIVATE :: name_len=30      
      INTEGER, PARAMETER, PRIVATE :: units_len=30      

!*******************************************************************************
! SECTION II. DERIVED-TYPE DECLARATIONS
!   coosig  [Combination Of SIGnals] Description:
!       coosig_desc  
!     Type of coosig specified by  % comb_type.
!     Allowable values of comb_type:
!       'sum'
!       'max'
!       'min'
!
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Declare type coosig_desc
!   Common to all comb_types
!       comb_type          character, type of combination
!          sum                s = Sum_over_i[a(i) * signal(index(i))]
!          max                s = Max[a(i) * signal(index(i))]
!          min                s = Min[a(i) * signal(index(i))]
!       name               character, name of coosig
!       units              character, physical units that the data is measured in
!       sigma_default      real, default value of the uncertainty in the data
!       n_sig              number of signals to combine
!       index              array of indices - to specify signals to combine
!       a                  array of real coefficients, for combination
!-------------------------------------------------------------------------------
      TYPE coosig_desc
         CHARACTER (len=type_len)           :: comb_type
         CHARACTER (len=name_len)           :: name
         CHARACTER (len=units_len)          :: units
         REAL(rprec)                        :: sigma_default
         INTEGER, DIMENSION(:), POINTER     :: index => null()
         REAL(rprec), DIMENSION(:), POINTER :: a => null()
      END TYPE coosig_desc

!*******************************************************************************
! SECTION III. INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Assignment for structures
!-------------------------------------------------------------------------------
      INTERFACE ASSIGNMENT (=)
         MODULE PROCEDURE coosig_desc_assign, coosig_desc_assign_a
      END INTERFACE

!-------------------------------------------------------------------------------
!  Generic construct
!-------------------------------------------------------------------------------
      INTERFACE coosig_construct
         MODULE PROCEDURE coosig_desc_construct
      END INTERFACE

!-------------------------------------------------------------------------------
!  Generic destroy
!-------------------------------------------------------------------------------
      INTERFACE coosig_destroy
         MODULE PROCEDURE coosig_desc_destroy
      END INTERFACE

!-------------------------------------------------------------------------------
!  Generic write
!-------------------------------------------------------------------------------
      INTERFACE coosig_write
         MODULE PROCEDURE coosig_desc_write
      END INTERFACE

!-------------------------------------------------------------------------------
!  Interface block for testing goes here. 
!-------------------------------------------------------------------------------

      CONTAINS
!*******************************************************************************
! SECTION IV. CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Construct a coosig_desc
!-------------------------------------------------------------------------------
      SUBROUTINE coosig_desc_construct(this,comb_type,name,                    &
     &   units,sigma_default,index,a)

      IMPLICIT NONE

!  Declare Arguments 
      TYPE (coosig_desc), INTENT(inout)          :: this
      CHARACTER (len=*), INTENT(in)              :: comb_type
      CHARACTER (len=*), INTENT(in)              :: name
      CHARACTER (len=*), INTENT(in)              :: units
      REAL(rprec), INTENT(in)                    :: sigma_default
      INTEGER, DIMENSION(:), INTENT(in)          :: index
      REAL(rprec), DIMENSION(:), INTENT(in)      :: a

!  Declare local variables
      INTEGER :: n_index, n_a, ier1, ier2
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'coosig_desc_construct: '

!  Start of executable code

!  Scalar assignments
      this % name = TRIM(ADJUSTL(name))
      this % units = TRIM(ADJUSTL(units))
      this % sigma_default = sigma_default

!  Different coding, depending on comb_type
      SELECT CASE (TRIM(ADJUSTL(comb_type)))
      CASE ('sum')
         this % comb_type = 'sum'
      CASE ('max')
         this % comb_type = 'max'
      CASE ('min')
         this % comb_type = 'min'

      CASE DEFAULT
         CALL err_fatal(sub_name // 'unrecognized comb_type: ',                &
     &      char=comb_type)
      END SELECT ! Different coding depending on comb_type
      
!  Take length of arrays from arguments
      n_index = SIZE(index)
      n_a = SIZE(a)

!  Error exit if the size of the data is not the same as the size of the sigma.
      CALL assert_eq(n_index,n_a,sub_name // 'n_index & n_a')
      
      ALLOCATE(this % index(1:n_index),STAT=ier1)
      ALLOCATE(this % a(1:n_a),STAT=ier2)
      CALL assert_eq(0,ier1,ier2,sub_name // 'alloc')

      this % index = index
      this % a = a

      END SUBROUTINE coosig_desc_construct

!*******************************************************************************
! SECTION V. DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Destroy a coosig_desc
!-------------------------------------------------------------------------------
      SUBROUTINE coosig_desc_destroy(this)
      IMPLICIT NONE

!  Declare Arguments 
      TYPE (coosig_desc), INTENT(inout) :: this

!  Declare local variables
      INTEGER :: ier1, ier2
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'coosig_desc_destroy: '

!  Start of executable code

!  Get rid of all components
      this % comb_type = ' '
      this % name = ' '
      this % units = ' '
      this % sigma_default = zero
      IF (ASSOCIATED(this % index)) DEALLOCATE(this % index,STAT=ier1)
      IF (ASSOCIATED(this % a)) DEALLOCATE(this % a,STAT=ier2)
      CALL assert_eq(0,ier1,ier2,sub_name // 'deallocate error')

      END SUBROUTINE coosig_desc_destroy

!*******************************************************************************
! SECTION VI. ASSIGNMENT SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Assignment for coosig_desc
!-------------------------------------------------------------------------------
      SUBROUTINE coosig_desc_assign(left,right)

!  Can't get by with intrinsic assignment, because of the 'use as allocatable
!  array" pointers index and a.

      IMPLICIT NONE
      
!  Declare Arguments 
      TYPE (coosig_desc), INTENT (inout) :: left
      TYPE (coosig_desc), INTENT (in) :: right
      
!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'coosig_desc_assign: '
      CHARACTER (len=*), PARAMETER :: err_mess1 =                              &
     & 'left-right index or a are the same?. FIX IT'
      INTEGER :: ier1, ier2
         
!  Start of executable code

!  Check to see if the 'use as allocatable array' pointers are pointing
!  to the same location
      CALL assert(.not.ASSOCIATED(left % index,right % index),                 &
     &            .not.ASSOCIATED(left % a,right % a),                         &
     &   sub_name // err_mess1)

!  Destroy the left
      CALL coosig_desc_destroy(left)
      
!  Allocate space for the 'use as allocatable array' pointer variables
      ALLOCATE(left % index(1:SIZE(right % index)),STAT=ier1)
      ALLOCATE(left % a(1:SIZE(right % a)),STAT=ier2)
      CALL assert_eq(0,ier1,ier2,sub_name // 'Allocation error')

      left % comb_type = right % comb_type
      left % name = right % name
      left % units = right % units
      left % sigma_default = right % sigma_default
      left % index = right % index
      left % a = right % a
         
      END SUBROUTINE coosig_desc_assign
!-------------------------------------------------------------------------------
!  Assignment for array of coosig_desc
!-------------------------------------------------------------------------------
      SUBROUTINE coosig_desc_assign_a(left,right)

      IMPLICIT NONE
      
!  Declare Arguments 
      TYPE (coosig_desc), DIMENSION(:), INTENT (inout) :: left
      TYPE (coosig_desc), DIMENSION(:), INTENT (in) :: right
      
!  Declare local variables
      INTEGER :: n_left, n_right, i
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'coosig_desc_assign_a: '
      CHARACTER (len=*), PARAMETER :: err_mess1 =                              &
     & 'left-right array lengths are not the same.'
         
!  Start of executable code

      n_left = SIZE(left)
      n_right = SIZE(right)
      CALL assert_eq(n_left,n_right,sub_name // err_mess1)
      DO i = 1,n_left
         left(i) = right(i)
      END DO
         
      END SUBROUTINE coosig_desc_assign_a
          
!*******************************************************************************
! SECTION VII.  OUTPUT SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Write out the contents of a coosig_desc
!-------------------------------------------------------------------------------

      SUBROUTINE coosig_desc_write(this,identifier,unit,verbose)
      IMPLICIT NONE

!  Declare Arguments
      TYPE (coosig_desc), INTENT (in) :: this
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
      INTEGER :: n_signals
      CHARACTER (len=60) :: id
      INTEGER      :: i

!  Declare Format array
      CHARACTER(len=*), PARAMETER, DIMENSION(7) :: fmt1 = (/                   &
     & '(" start coosig_desc write, called with id = ",a)   ',                 &
     & '(" comb_type = ",a)                                 ',                 &
     & '(" name = ",a)                                      ',                 &
     & '(" units = ",a)                                     ',                 &
     & '(" sigma_default = ",es10.3)                        ',                 &
     & '(" i   index     a       ",/,(1x,2(i4,3x),es12.5))  ',                 &
     & '(" end coosig_desc write, called with id = ",a)     '                  &
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
         WRITE(iou,*) this % comb_type
         WRITE(iou,*) this % name
         WRITE(iou,*) this % units
         WRITE(iou,*) this % sigma_default
      
      CASE(1:)    ! Default, more verbose
         n_signals = SIZE(this % index)
         WRITE(iou,fmt1(1)) id
         WRITE(iou,fmt1(2)) this % comb_type
         WRITE(iou,fmt1(3)) this % name
         WRITE(iou,fmt1(4)) this % units
         WRITE(iou,fmt1(5)) this % sigma_default
         WRITE(iou,fmt1(6)) (i,this % index(i),this % a(i),                    &
     &      i=1,n_signals)
         WRITE(iou,fmt1(7)) id
      
      END SELECT

      END SUBROUTINE coosig_desc_write

!*******************************************************************************
! SECTION XVI.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
!
!  JDH 2011-07-06.  First version of coosig_T. Copied and edited from
!      geometric_T, used some examples from signal_T
!
!  JDH 2011-09-07.  
!      Modified coosig_desc_write - eliminate n_data in favor of n_signals
!
!  ---------- Below - comments from diagnostic_T  ------------------------------

      END MODULE coosig_T
