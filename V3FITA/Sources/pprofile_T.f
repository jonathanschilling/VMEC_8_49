!*******************************************************************************
!  File pprofile_T.f
!  Contains module pprofile_T
!  Defines derived-types: pprofile

!  Module is part of the V3FIT stellarator equilibrium reconstruction code.


!*******************************************************************************
!  MODULE pprofile_T
! SECTION I.      VARIABLE DECLARATIONS
! SECTION II.     DERIVED TYPE DECLARATIONS
! SECTION III.    INTERFACE BLOCKS
! SECTION IV.     CONSTRUCTION SUBROUTINES
! SECTION V.      DESTRUCTION SUBROUTINES
! SECTION VI.     ASSIGNMENT SUBROUTINES
! SECTION VII.    EVALUATION SUBROUTINES
! SECTION VIII.   PUT SUBROUTINES
! SECTION IX.     GET SUBROUTINES
! SECTION X.      WRITE SUBROUTINES
! SECTION XI.     COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
      MODULE pprofile_T

!*******************************************************************************
! SECTION I. VARIABLE DECLARATIONS
!*******************************************************************************

!-------------------------------------------------------------------------------
!  Type declarations - lengths of reals, integers, and complexes.
!-------------------------------------------------------------------------------
      USE stel_kinds

!-------------------------------------------------------------------------------
!  Frequently used mathematical constants, lots of extra precision.
!-------------------------------------------------------------------------------
      USE stel_constants

!-------------------------------------------------------------------------------
!  Use Statements for  V3 Utilities
!-------------------------------------------------------------------------------
      USE v3_utilities
      
      IMPLICIT NONE
!-------------------------------------------------------------------------------
!  Lengths of Character Variables
!-------------------------------------------------------------------------------
      INTEGER, PARAMETER, PRIVATE :: p_type_len=20      
!-------------------------------------------------------------------------------
!  Array Bounds for the components
!-------------------------------------------------------------------------------
      INTEGER, PARAMETER, PRIVATE :: ilb_b = 0     
      INTEGER, PARAMETER, PRIVATE :: iub_b = 20    
      INTEGER, PARAMETER, PRIVATE :: iub_asf = 101    
!-------------------------------------------------------------------------------
!  Other Dependencies
!    1) pprofile_construct uses to_lower, from LIBSTELL
!-------------------------------------------------------------------------------

!*******************************************************************************
! SECTION II. DERIVED TYPE (STRUCTURE) DECLARATIONS
!
!   Derived Type for Parameterized Profiles
!      Variables of type pprofile will be used for quantities such as the
!      electron density, ion density, electron temperature
!      They are assumed to be functions of the radial coordinate, which varies 
!      from 0. to 1. The derived type will carry the parameters and the type
!      of parameterization. A module function will evaluate the profile as
!      a function of s.
!
!      The code is partly modelled after the code in VMEC2000 file 
!      profile_functions.f, which implements the current, pressure, and iota
!      profiles.
!
!      One difference from the VMEC2000 profiles is that a pprofile will need
!      to return a value for s arguments larger than 1, indicating that the
!      point is outside the plasma.
!
!      I call the parameters b (not ac, am, ai, as in VMEC2000) to remind the
!      user that these parameterizations are DIFFERENT from the VMEC2000
!      parameterizations.
!
!*******************************************************************************
!
!  
!-------------------------------------------------------------------------------
!  Declare type pprofile
!
!  p_type      Type of parameterization, character variable
!  b           Array of parameters, dimensioned (0:20)
!  as          Array of spline s values (independent variables)
!  af          Array of spline function values (dependent variables)
!-------------------------------------------------------------------------------
      TYPE pprofile
!  Variables for a parameterized profile

         CHARACTER (len = p_type_len)                :: p_type
         REAL(rprec), DIMENSION(ilb_b:iub_b)         :: b              ! (0:20)
         REAL(rprec), DIMENSION(iub_asf)             :: as, af
      END TYPE pprofile

!*******************************************************************************
! SECTION III. INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Generic get 
!-------------------------------------------------------------------------------
      INTERFACE pprofile_get
         MODULE PROCEDURE pprofile_get_ch, pprofile_get_real
      END INTERFACE

      CONTAINS
!*******************************************************************************
! SECTION IV. CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Construct an pprofile
!
!-------------------------------------------------------------------------------
      SUBROUTINE pprofile_construct(this, p_type_arg, b_arg,                   &
     &   as_arg, af_arg)

      IMPLICIT NONE

!  Declare Arguments 
      TYPE (pprofile), INTENT (inout)                 :: this
      CHARACTER (len=*), INTENT(in)                   :: p_type_arg
      REAL(rprec), DIMENSION(:), INTENT(in)           :: b_arg
      REAL(rprec), DIMENSION(:), OPTIONAL, INTENT(in) :: as_arg, af_arg

!  Declare local variables
      INTEGER :: n_b_copy, ilb_b_arg, n_as_copy, n_af_copy
      CHARACTER (len = p_type_len)                :: p_type_lc
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'pprofile_construct: '

!  Start of executable code

!  Convert the type to lower case, and truncate
      p_type_lc = p_type_arg
      CALL tolower(p_type_lc)
      SELECT CASE(TRIM(p_type_lc))
      
         CASE ('two_power','power_series','none','cubic_spline',               &
     &      'akima_spline')
            this % p_type = p_type_lc
      
         CASE DEFAULT
            this % p_type = 'none'
            WRITE(*,*) sub_name, 'Unrecognized p_type:', p_type_lc
            WRITE(*,*) ' *** CHECK YOUR INPUT ***'            
         
      END SELECT

! Store the coefficients, being careful with the indices.
      n_b_copy = MIN(SIZE(b_arg),iub_b - ilb_b + 1)
      ilb_b_arg = LBOUND(b_arg,1)
      this % b = zero
      this % b(ilb_b:ilb_b + n_b_copy) =                                       &
     &    b_arg(ilb_b_arg: ilb_b_arg + n_b_copy)
      
      this % as = zero
      this % af = zero
      IF (PRESENT(as_arg)) THEN
         n_as_copy = MIN(SIZE(as_arg),iub_asf)
         this % as(1:n_as_copy) = as_arg(1:n_as_copy)         
      ENDIF
      IF (PRESENT(af_arg)) THEN
         n_af_copy = MIN(SIZE(af_arg),iub_asf)
         this % af(1:n_af_copy) = af_arg(1:n_af_copy)         
      ENDIF
      
      END SUBROUTINE pprofile_construct

!*******************************************************************************
! SECTION V. DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Destroy an p_profile
!-------------------------------------------------------------------------------
      SUBROUTINE pprofile_destroy(this)
      IMPLICIT NONE

!  Declare Arguments 
      TYPE (pprofile), INTENT(inout) :: this

!  Start of executable code
!  Get rid of all components

      this % p_type   = ''
      this % b        = zero
      this % as       = zero
      this % af       = zero

      END SUBROUTINE pprofile_destroy

!*******************************************************************************
! SECTION VI. ASSIGNMENT SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Assignment for pprofile - Default OK as no Pointers
!-------------------------------------------------------------------------------

!*******************************************************************************
! SECTION VII. EVALUATION FUNCTION
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Evaluate a pprofile at an argument value
!-------------------------------------------------------------------------------

      FUNCTION pprofile_evaluate(this,s_arg)

      IMPLICIT NONE

!  Declare Arguments 
      TYPE (pprofile), INTENT (in)               :: this
      REAL(rprec), INTENT(in)                    :: s_arg
      REAL(rprec)                                :: pprofile_evaluate

!  Declare local variables
      REAL(rprec) :: s_use, s_01
      LOGICAL     :: l_01
      INTEGER     :: i, i_max, iflag
      REAL(rprec) :: pprofile_value
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'pprofile_evaluate: '

!  Start of executable code

!  Define argument and logical variable for s between zero and one
      s_01 = MAX(zero,MIN(one,s_arg))
      l_01 = (s_arg .ge. zero) .and. (s_arg .le. one)

      SELECT CASE(TRIM(this % p_type))
      
         CASE ('none')
            pprofile_value = -1.E49_rprec
         
         CASE ('two_power')
            pprofile_value = this % b(0)
            IF (l_01) THEN
               pprofile_value = pprofile_value + this % b(1) *                 &
     &            ((1 - s_arg ** this % b(2)) ** this % b(3))
            ENDIF
      
         CASE ('power_series')
         !  use s_01, so that value is constant outside
            pprofile_value = zero
            DO i = iub_b, ilb_b, -1      ! all of b array, backwards
               pprofile_value = (pprofile_value + this % b(i)) * s_01
            END DO

         CASE ('cubic_spline')
!  Cubic splines come from LIBSTELL/miscel/cubic_spline.f
!  As of 2011-10-21, this subroutine is used in VMEC profile_functions.f
            i_max = MINLOC(this % as(2:),dim=1)
            IF (i_max .lt. 4) THEN
               WRITE(*,*) 'pprofile:cubic spline: too few as values'
               WRITE(*,*) 'i_max, as = ', i_max
               WRITE(*,*) this % as
               STOP
            ENDIF
            s_use = MIN(this % as(i_max - 1), MAX(s_arg,this % as(1)))
            CALL spline_cubic(s_use,pprofile_value,this % as,this % af,        &
     &         i_max,iflag)
            IF (iflag < 0) THEN
               IF(iflag == -1) THEN
                  STOP 'pprofile: outside value from spline_cubic'
               ELSEIF(iflag == -2) THEN
                  STOP 'pprofile:  decreasing s values in spline_cubic'
               ELSE
                  STOP 'pprofile: unknown error from spline_cubic'
               ENDIF
            ENDIF

         CASE ('akima_spline')
!  Akima splines come from LIBSTELL/miscel/spline_akima.f
!  As of 2011-10-21, this subroutine is used in VMEC profile_functions.f
            i_max = MINLOC(this % as(2:),dim=1)
            IF (i_max .lt. 4) THEN
               WRITE(*,*) 'pprofile:akima spline: too few as values'
               WRITE(*,*) 'i_max, as = ', i_max
               WRITE(*,*) this % as
               STOP
            ENDIF
            s_use = MIN(this % as(i_max - 1), MAX(s_arg,this % as(1)))
            CALL spline_akima(s_use,pprofile_value,this % as,this % af,        &
     &         i_max,iflag)
            IF (iflag < 0) THEN
               STOP 'pprofile: bad value from spline_akima requested'
            ENDIF

         CASE DEFAULT
            WRITE(*,*) sub_name, 'Unrecognized p_type:', this % p_type
            WRITE(*,*) ' SHOULD HAVE BEEN CORRECTLY SET IN ',                  &
     &         'pprofile_construct'
         
      END SELECT
      
      pprofile_evaluate = pprofile_value
      
      END FUNCTION pprofile_evaluate

!*******************************************************************************
! SECTION VIII. PUT SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Put a value into a pprofile
!-------------------------------------------------------------------------------
      SUBROUTINE pprofile_put(this,id,index,value)
      IMPLICIT NONE

!  Declare Arguments 
      TYPE (pprofile), INTENT(inout) :: this
      CHARACTER(len=*), INTENT(in) :: id
      INTEGER, INTENT(in) :: index
      REAL(rprec), INTENT(in) :: value

!  Start of executable code

      SELECT CASE(TRIM(id))
      
      CASE('b')
         IF((index .ge. ilb_b) .and. (index .le. iub_b)) THEN
            this % b(index) = value
         ENDIF
      
      CASE('as')
         IF((index .ge. 1) .and. (index .le. iub_asf)) THEN
            this % as(index) = value
         ENDIF

      CASE('af')
         IF((index .ge. 1) .and. (index .le. iub_asf)) THEN
            this % af(index) = value
         ENDIF

      CASE DEFAULT
         WRITE(*,*) ' BAD ID ARGUMENT, pprofile_put. id=', id
      
      END SELECT   

      END SUBROUTINE pprofile_put
!*******************************************************************************
! SECTION IX. GET SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Get a character value from a pprofile
!    CALL pprofile_get(pp_ne,'p_type',p_type_return)
!-------------------------------------------------------------------------------
      SUBROUTINE pprofile_get_ch(this,id,ch_value)
      IMPLICIT NONE

!  Declare Arguments 
      TYPE (pprofile), INTENT(in) :: this
      CHARACTER(len=*), INTENT(in) :: id
      CHARACTER(len=*), INTENT(inout) :: ch_value

!  Start of executable code

!  The only character information carried is the p_type. Return it

      ch_value = this % p_type

      END SUBROUTINE pprofile_get_ch
!-------------------------------------------------------------------------------
!  Get a real value from a pprofile internal array
!    CALL pprofile_get(pp_ne,'b',i,b_i)
!-------------------------------------------------------------------------------
      SUBROUTINE pprofile_get_real(this,id,index,value)
      IMPLICIT NONE

!  Declare Arguments 
      TYPE (pprofile), INTENT(in) :: this
      CHARACTER(len=*), INTENT(in) :: id
      INTEGER, INTENT(in) :: index
      REAL(rprec), INTENT(out) :: value

!  Start of executable code

      value = -2.E49_rprec     ! indicator that indices out of range
      SELECT CASE(TRIM(id))
      
      CASE('b')
         IF((index .ge. ilb_b) .and. (index .le. iub_b)) THEN
            value = this % b(index)
         ENDIF
      
      CASE('as')
         IF((index .ge. 1) .and. (index .le. iub_asf)) THEN
            value = this % as(index)
         ENDIF

      CASE('af')
         IF((index .ge. 1) .and. (index .le. iub_asf)) THEN
            value = this % af(index)
         ENDIF

      CASE DEFAULT
         WRITE(*,*) ' BAD ID ARGUMENT, pprofile_get. id=', id
      
      END SELECT   

      END SUBROUTINE pprofile_get_real

*******************************************************************************
! SECTION X. WRITE SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Write pprofile
!    CALL pprofile_write(pp_ne,'pp_ne',iou)
!-------------------------------------------------------------------------------
      SUBROUTINE pprofile_write(this,identifier,unit)
      IMPLICIT NONE

!  Declare Arguments 
      TYPE (pprofile), INTENT(in) :: this
      CHARACTER (len=*), INTENT(in), OPTIONAL :: identifier
      INTEGER, INTENT(in), OPTIONAL :: unit
!  identifier   character variable, also written out
!  unit         I/O unit number to write to

!  Declare local variables and constants
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'pprofile_write: '
      INTEGER :: iou_default = 6
      INTEGER :: iou
      INTEGER :: i, i_max
      CHARACTER (len=60) :: id

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

!  Actual write. Could do a select-case on p_type.
      WRITE(iou,1100) id
      WRITE(iou,1200) this % p_type
      
      SELECT CASE(TRIM(this % p_type))
      
         CASE ('none')
         
         CASE ('two_power')
            WRITE(iou,1210)
            WRITE(iou,1300) this % b(0:3)
      
         CASE ('power_series')
            WRITE(iou,1220)
            WRITE(iou,1300) this % b(ilb_b:iub_b)
            
         CASE ('cubic_spline','akima_spline')
            i_max = MINLOC(this % as(2:),dim=1) ! last value used for splines
            WRITE(iou,1230)
            WRITE(iou,1231) (i, this % as(i), this % af(i),i=1,i_max)

         CASE DEFAULT
            WRITE(*,*) sub_name, 'Unrecognized p_type:', this % p_type
            WRITE(*,*) ' SHOULD HAVE BEEN CORRECTLY SET IN ',                  &
     &         'pprofile_construct'
         
      END SELECT

!      WRITE(iou,1400) id

1100  FORMAT(/' Parameterized Profile Write: id = ',a)
1200  FORMAT('   pp_type = ',a)
1210  FORMAT(' b_0 + Th(s)Th(1-s)(b_1 (1 - s ** b_2) ** b_3).',                &
     &   '   b(0:3) = ')
1220  FORMAT(' Th(s)Th(1-s)[Sum_0_n b_i s** i].   b(0:n) = ')
1230  FORMAT(' i       as(i)           af(i)')
1231  FORMAT(1x,i3,2x,es15.8,2x,es15.8)
1300  FORMAT(4(2x,es15.8))
1400  FORMAT(" END Parameterized Profile Write: id = ",a)

      END SUBROUTINE pprofile_write

!*******************************************************************************
! SECTION XI.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
!
!  To Do (2011-10-21)
! 1) add other p_type
!
!  2011-10-21 JDH
!    Add coding for cubic aplines and Akima splines
!
!  2011-10-19 JDH
!    Some code changes
!
!  2011-07-26 JDH
!    First version, based on module density_T
!    Wait on use until have refactored model to have appropriate methods.
           
      END MODULE pprofile_T
