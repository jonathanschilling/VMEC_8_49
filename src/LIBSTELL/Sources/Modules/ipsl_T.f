!*******************************************************************************
!  File ipsl_T.f
!  Contains module ipsl_T
!  Defines derived-types: ipsl_desc
!  A type of Diagnostic - Interferometry Polarimetry Straight Line
!   (A microwave/laser diagnostic that computes line-integrated B field & density)
!
!*******************************************************************************
!  MODULE ipsl_T
!    (IPSL Type Definition, for the V3FIT code)
! SECTION I.    VARIABLE DECLARATIONS
! SECTION II.   DERIVED-TYPE DECLARATIONS
! SECTION III.  INTERFACE BLOCKS
! SECTION IV.   CONSTRUCTION SUBROUTINES
! SECTION V.    DESTRUCTION SUBROUTINES
! SECTION VI.   ASSIGNMENT SUBROUTINES
! SECTION VII.  OUTPUT SUBROUTINES

! SECTION XVI.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
      MODULE ipsl_T

!*******************************************************************************
! SECTION I. VARIABLE DECLARATIONS
!*******************************************************************************

!-------------------------------------------------------------------------------
!  Type declarations - lengths of reals, integers, and complexes.
!  Frequently used mathematical constants, lots of extra precision.
!-------------------------------------------------------------------------------
      USE stel_kinds, only : rprec, iprec, cprec
      USE stel_constants, only : pi, twopi, one, zero
     
!-------------------------------------------------------------------------------
!  Use Statements for other structures, V3 Utilities
!-------------------------------------------------------------------------------
      USE ip_beamline
      USE v3_utilities

!-------------------------------------------------------------------------------
!  Implicit None comes after USE statements, before other declarations
!-------------------------------------------------------------------------------
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Make type declarations and constants Private, so there are no conflicts.
!-------------------------------------------------------------------------------
      PRIVATE rprec, iprec, cprec, pi, twopi, one, zero

!-------------------------------------------------------------------------------
!  Lengths of Character Variables
!-------------------------------------------------------------------------------
      INTEGER(iprec), PARAMETER, PRIVATE :: type_len=10      
      INTEGER(iprec), PARAMETER, PRIVATE :: sn_len=30      
      INTEGER(iprec), PARAMETER, PRIVATE :: ln_len=80      
      INTEGER(iprec), PARAMETER, PRIVATE :: units_len=30      

!*******************************************************************************
! SECTION II. DERIVED-TYPE DECLARATIONS
! 1)  IPSL Description:
!       ipsl_desc  
!     Type of diagnostic specified by  % d_type = 'ipsl'.
!
!
!*******************************************************************************

!-------------------------------------------------------------------------------
!  Declare type ipsl_desc
!       s_name          character, short name of individual beam diagnostic
!       l_name          character, long name of individual beam diagnostic
!       ip_sname        character, short name of source device producing beam
!       ip_lname        character, long name of source device producing beam
!       units           character, physical units that the data is measured in
!       sigma_default   real, default value of the uncertainty in the data
!       ipsl_type       character, keyword from the diagnostic_dot file
!       l_ipbeam_def    logical, definition status of the ipbeam component     
!       ipbeam          type ip_beam, description of int/pol beam
!!-------------------------------------------------------------------------------
      TYPE ipsl_desc
         CHARACTER (len=sn_len)         :: s_name                                 
         CHARACTER (len=ln_len)         :: l_name
         CHARACTER (len=sn_len)         :: ip_sname                                 
         CHARACTER (len=ln_len)         :: ip_lname
         CHARACTER (len=units_len)      :: units                                 
         CHARACTER (len=30)             :: ipsl_type                                 
         LOGICAL                        :: l_ipbeam_def
         REAL(rprec)                    :: sigma_default
         TYPE (ip_beam)                 :: ipbeam 
      END TYPE ipsl_desc

!*******************************************************************************
! SECTION III. INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Assignment for structures
!-------------------------------------------------------------------------------
      INTERFACE ASSIGNMENT (=)
         MODULE PROCEDURE ipsl_desc_assign  
      END INTERFACE

!-------------------------------------------------------------------------------
!  Generic construct
!-------------------------------------------------------------------------------
      INTERFACE ipsl_construct
         MODULE PROCEDURE ipsl_desc_construct
      END INTERFACE

!-------------------------------------------------------------------------------
!  Generic destroy
!-------------------------------------------------------------------------------
      INTERFACE ipsl_destroy
         MODULE PROCEDURE ipsl_desc_destroy
      END INTERFACE

!-------------------------------------------------------------------------------
!  Generic write
!-------------------------------------------------------------------------------
      INTERFACE ipsl_write
         MODULE PROCEDURE ipsl_desc_write
      END INTERFACE


      CONTAINS
!*******************************************************************************
! SECTION IV. CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Construct a ipsl_desc
!
!  For d_type = 'ipsl' (inteferometer/polarimeter straight line)
!-------------------------------------------------------------------------------
      SUBROUTINE ipsl_desc_construct(this,s_name,l_name, ip_sname,             &
     &   ip_lname, units, sigma_default,ipsl_type,ipbeam)

!  NB. 
!  The argument is assigned to the 'this' component. Do NOT call this
!  subroutine with this % ipsl as the ipsl argument.

      IMPLICIT NONE

!  Declare Arguments 
      TYPE (ipsl_desc), INTENT(inout)            :: this
      CHARACTER (len=*), INTENT(in)              :: s_name
      CHARACTER (len=*), INTENT(in)              :: l_name
      CHARACTER (len=*), INTENT(in)              :: ip_sname
      CHARACTER (len=*), INTENT(in)              :: ip_lname
      CHARACTER (len=*), INTENT(in)              :: units
      CHARACTER (len=*), INTENT(in)              :: ipsl_type
      REAL(rprec), INTENT(in)                    :: sigma_default
!      TYPE (ip_beam), INTENT(in), TARGET         :: ipbeam ! Why a TARGET ? 2007-06-11
      TYPE (ip_beam), INTENT(in)                 :: ipbeam 

!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'ipsl_desc_construct: '

!  Start of executable code

!  Destroy the ipbeam component
      CALL ip_beam_destroy(this % ipbeam)

!  Scalar assignments
      this % s_name = TRIM(ADJUSTL(s_name))
      this % l_name = TRIM(ADJUSTL(l_name))
      this % ip_sname = TRIM(ADJUSTL(ip_sname))
      this % ip_lname = TRIM(ADJUSTL(ip_lname))
      this % units = TRIM(ADJUSTL(units))
      this % ipsl_type =  TRIM(ADJUSTL(ipsl_type))
      this % l_ipbeam_def = .TRUE.


!  Derived Type Assignments
      this % ipbeam = ipbeam
      
      END SUBROUTINE ipsl_desc_construct


!*******************************************************************************
! SECTION V. DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Destroy a ipsl_desc
!-------------------------------------------------------------------------------
      SUBROUTINE ipsl_desc_destroy(this)
      IMPLICIT NONE

!  Declare Arguments 
      TYPE (ipsl_desc), INTENT(inout) :: this

!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'ipsl_desc_destroy: '

!  Start of executable code

!  Destroy scalar components
      this % s_name = ' '
      this % l_name = ' '
      this % ip_sname = ' '
      this % ip_lname = ' '
      this % units = ' '
      this % ipsl_type = ' '
      this % sigma_default = zero
      this % l_ipbeam_def = .FALSE.

!  Destroy Derived Type
      CALL ip_beam_destroy(this % ipbeam)

      END SUBROUTINE ipsl_desc_destroy


!*******************************************************************************
! SECTION VI. ASSIGNMENT SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Assignment for ipsl_desc
!-------------------------------------------------------------------------------
      SUBROUTINE ipsl_desc_assign(left,right)


      IMPLICIT NONE
      
!  Declare Arguments 
      TYPE (ipsl_desc), INTENT (inout) :: left
      TYPE (ipsl_desc), INTENT (in) :: right
      
!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'ipsl_desc_assign: '
         
!  Start of executable code
      left % s_name = right % s_name
      left % l_name = right % l_name
      left % ip_sname = right % ip_sname
      left % ip_lname = right % ip_lname
      left % units = right % units
      left % ipsl_type = right % ipsl_type
      left % l_ipbeam_def = right % l_ipbeam_def
      left % sigma_default = right % sigma_default
      left % ipbeam = right % ipbeam
         
      END SUBROUTINE ipsl_desc_assign
!-------------------------------------------------------------------------------

          
!*******************************************************************************
! SECTION VII.  OUTPUT SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Write out the contents  of a ipsl_desc
!-------------------------------------------------------------------------------

      SUBROUTINE ipsl_desc_write(this,identifier,unit,verbose)
      IMPLICIT NONE

!  Declare Arguments
      TYPE (ipsl_desc), INTENT (in) :: this
      CHARACTER (len=*), INTENT(in), OPTIONAL :: identifier
      INTEGER(iprec), INTENT(in), OPTIONAL :: unit
      INTEGER(iprec), INTENT(in), OPTIONAL :: verbose
!  identifier   character variable, also written out
!  unit         I/O unit number to write to
!  verbose      integer, to specify verbosity level of write

!  Declare local variables and constants
      INTEGER(iprec) :: iv_default = 1
      INTEGER(iprec) :: iv
      INTEGER(iprec) :: iou_default = 6
      INTEGER(iprec) :: iou
      CHARACTER (len=60) :: id

!  Declare Format array
      CHARACTER(len=*), PARAMETER, DIMENSION(11) :: fmt1 = (/                  &
     & '(" start ipsl_desc write, called with id = ",a)      ',                &
     & '(" s_name = ",a)                                     ',                &
     & '(" l_name = ",a)                                     ',                &
     & '(" ip_sname = ",a)                                   ',                &
     & '(" ip_lname = ",a)                                   ',                &
     & '(" units = ",a)                                      ',                &
     & '(" l_ipbeam_def = ",L1)                              ',                &
     & '(" ipsl_type = ",a)                                  ',                &
     & '(" ip_beam s_name = ",a)                             ',                &
     & '(" sigma_default = ",es12.5)                         ',                &
     & '(" end ipsl_desc write, called with id = ",a)        '                 &
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
         WRITE(iou,*) this % s_name
         WRITE(iou,*) this % l_name
         WRITE(iou,*) this % ip_sname
         WRITE(iou,*) this % ip_lname
         WRITE(iou,*) this % units
         WRITE(iou,*) this % l_ipbeam_def
         WRITE(iou,*) this % ipsl_type
         WRITE(iou,*) this % ipbeam % s_name
         WRITE(iou,*) this % sigma_default
      
      CASE(1:)    ! Default, more verbose
         WRITE(iou,fmt1(1)) id
         WRITE(iou,fmt1(2)) this % s_name
         WRITE(iou,fmt1(3)) this % l_name
         WRITE(iou,fmt1(4)) this % ip_sname
         WRITE(iou,fmt1(5)) this % ip_lname
         WRITE(iou,fmt1(6)) this % units
         WRITE(iou,fmt1(7)) this % l_ipbeam_def
         WRITE(iou,fmt1(8)) this % ipsl_type
         WRITE(iou,fmt1(9)) this % ipbeam % s_name
         WRITE(iou,fmt1(10)) this % sigma_default
         WRITE(iou,fmt1(11)) id
      
      END SELECT

      END SUBROUTINE ipsl_desc_write



!*******************************************************************************
! SECTION XVI.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
!
!  JMS 2007-06-14. First version of ipsl_T. Copied and edited from mddc_T
!
!  


       
      END MODULE ipsl_T
