!*******************************************************************************
!  File ipsl_cdf.f
!  Contains the module ipsl_cdf
!    Module for defining variables and writing netCDF files, and reading
!    netCDF files with the derived types ipsl_desc and ipsl_data
!    (from the ipsl_T module).
!    
!    Information about the  EZcdf module is at:
!       http://w3.pppl.gov/NTCC/EZcdf/
!
!-------------------------------------------------------------------------------
!   DEPENDENCIES
!-------------------------------------------------------------------------------
!
!    This module uses the following modules:
!       stel_kinds
!       stel_constants
!       ipsl_T
!       bsc
!       bsc_cdf
!       ezcdf
!       v3_utilities
!
!-------------------------------------------------------------------------------
!   CHANGE HISTORY
!-------------------------------------------------------------------------------
!
!  See Section X at the end of the module.
!
!-------------------------------------------------------------------------------
!   USAGE
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!   COMMENTS
!
! 1) All cdf_define calls must be completed before the first cdf_write call.
!-------------------------------------------------------------------------------
!
!*******************************************************************************

!*******************************************************************************
!  MODULE ipsl_cdf
!    
! SECTION I. VARIABLE DECLARATIONS
! SECTION II. INTERFACE BLOCKS
! SECTION III. DEFINE SUBROUTINES
! SECTION IV. WRITE SUBROUTINES
! SECTION V. READ SUBROUTINES
! SECTION VI. AUXILIARY FUNCTIONS AND SUBROUTINES

! SECTION X.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************

      MODULE ipsl_cdf

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
!  Modules to USE
!-------------------------------------------------------------------------------
      USE ipsl_T
      USE ip_beamline
      USE intpol_cdf
      USE ezcdf
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
!  Variable Names for netCDF. Make them Private.
!-------------------------------------------------------------------------------

      CHARACTER (LEN=*), PRIVATE, PARAMETER ::                                 &
     &  vn_s_name = 'ipsl_desc_s_name',                                        &            
     &  vn_l_name = 'ipsl_desc_l_name',                                        &            
     &  vn_ip_sname = 'ipsl_desc_ip_sname',                                        &            
     &  vn_ip_lname = 'ipsl_desc_ip_lname',                                        &            
     &  vn_units = 'ipsl_desc_units',                                          &            
     &  vn_sigma_default = 'ipsl_desc_sigma_default',                          &            
     &  vn_l_ipbeam_def = 'ipsl_desc_l_ipbeam_def',                            &            
     &  vn_ipsl_type = 'ipsl_desc_ipsl_type'

      CHARACTER (LEN=64), PRIVATE ::                                           &
     &  vn_s_name_use,                                                         &
     &  vn_l_name_use,                                                         &
     &  vn_ip_sname_use,                                                         &
     &  vn_ip_lname_use,                                                         &
     &  vn_units_use,                                                          &         
     &  vn_sigma_default_use,                                                  &            
     &  vn_l_ipbeam_def_use,                                                   &            
     &  vn_ipsl_type_use



!-------------------------------------------------------------------------------
!  Lengths of Character Variables
!-------------------------------------------------------------------------------
      INTEGER(iprec), PARAMETER, PRIVATE :: type_len=10      
      INTEGER(iprec), PARAMETER, PRIVATE :: sn_len=30      
      INTEGER(iprec), PARAMETER, PRIVATE :: ln_len=80      
      INTEGER(iprec), PARAMETER, PRIVATE :: units_len=30      

!-------------------------------------------------------------------------------
!  Generic Define
!-------------------------------------------------------------------------------
      INTERFACE ipsl_cdf_define
         MODULE PROCEDURE ipsl_cdf_define_desc
      END INTERFACE

!-------------------------------------------------------------------------------
!  Generic Write
!-------------------------------------------------------------------------------
      INTERFACE ipsl_cdf_write
         MODULE PROCEDURE ipsl_cdf_write_desc
       END INTERFACE

!-------------------------------------------------------------------------------
!  Generic Read
!-------------------------------------------------------------------------------
      INTERFACE ipsl_cdf_read
         MODULE PROCEDURE ipsl_cdf_read_desc
       END INTERFACE
!-------------------------------------------------------------------------------

 

      CONTAINS
!*******************************************************************************
! SECTION III. DEFINE SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------
      SUBROUTINE ipsl_cdf_define_desc(this,iou,prefix)
!  Subroutine to do the appropriate netCDF definition calls for a ipsl_desc
!  
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      TYPE (ipsl_desc), INTENT (in)            :: this
      INTEGER, INTENT(in)                      :: iou
      CHARACTER (len=*), INTENT(in), OPTIONAL  :: prefix

!  this        ipsl_desc - this is what gets defined and written.
!  iou         i/o unit number of the netCDF file. 
!  prefix      character -  prefixed to variable names, so that netCDF
!                  doesn't have problems with repeated, identical names
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'ipsl_cdf_define_desc: '
      CHARACTER(len=32) :: prefix_use

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

! Define the prefix to actually use
      IF (PRESENT(prefix)) THEN
         prefix_use = TRIM(ADJUSTL(prefix))
      ELSE
         prefix_use = ' '
      ENDIF

! Define all vn_--_use variable names
      CALL ipsl_cdf_defvn_desc(prefix_use)
         
! Scalar Components
      CALL cdf_define(iou, TRIM(vn_s_name_use), this % s_name)
      CALL cdf_define(iou, TRIM(vn_l_name_use), this % l_name)
      CALL cdf_define(iou, TRIM(vn_ip_sname_use), this % ip_sname)
      CALL cdf_define(iou, TRIM(vn_ip_lname_use), this % ip_lname)
      CALL cdf_define(iou, TRIM(vn_units_use), this % units)
      CALL cdf_define(iou, TRIM(vn_sigma_default_use),                         &
     &   this % sigma_default)
      CALL cdf_define(iou, TRIM(vn_l_ipbeam_def_use),                          &
     &      this % l_ipbeam_def)
      CALL cdf_define(iou, TRIM(vn_ipsl_type_use),                             &
     &      this % ipsl_type)

!  define ip_beam
      CALL intpol_cdf_define(this % ipbeam,iou,prefix)

     
      RETURN
      
      END SUBROUTINE ipsl_cdf_define_desc



!*******************************************************************************
! SECTION III. WRITE SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------
      SUBROUTINE ipsl_cdf_write_desc(this,iou,prefix)
!  Subroutine to do the appropriate netCDF write calls for a ipsl_desc
!  
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      TYPE (ipsl_desc), INTENT (in)        :: this
      INTEGER, INTENT (in)                       :: iou
      CHARACTER (len=*), INTENT (in), OPTIONAL   :: prefix

!  this        ipsl_desc - this is what gets written.
!  iou         i/o unit number of the netCDF file. 
!  prefix      character -  prefixed to variable names, so that netCDF
!                  doesn't have problems with repeated, identical names
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'ipsl_cdf_write_desc: '
      CHARACTER(len=32) :: prefix_use

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
! Define the prefix to actually use
      IF (PRESENT(prefix)) THEN
         prefix_use = TRIM(ADJUSTL(prefix))
      ELSE
         prefix_use = ' '
      ENDIF

! Define all vn_--_use variable names
      CALL ipsl_cdf_defvn_desc(prefix_use)

! Scalar Components
      CALL cdf_write(iou, TRIM(vn_s_name_use), this % s_name)
      CALL cdf_write(iou, TRIM(vn_l_name_use), this % l_name)
      CALL cdf_write(iou, TRIM(vn_ip_sname_use), this % ip_sname)
      CALL cdf_write(iou, TRIM(vn_ip_lname_use), this % ip_lname)
      CALL cdf_write(iou, TRIM(vn_units_use), this % units)
      CALL cdf_write(iou, TRIM(vn_sigma_default_use),                         &
     &   this % sigma_default)
      CALL cdf_write(iou, TRIM(vn_l_ipbeam_def_use),                          &
     &      this % l_ipbeam_def)
      CALL cdf_write(iou, TRIM(vn_ipsl_type_use),                             &
     &      this % ipsl_type)

!  ip_beam part
      CALL intpol_cdf_write(this % ipbeam,iou,prefix)
 
      
      RETURN
      
      END SUBROUTINE ipsl_cdf_write_desc


!*******************************************************************************
! SECTION IV. READ SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------
      SUBROUTINE ipsl_cdf_read_desc(this,iou,prefix)
!  Subroutine to do the appropriate netCDF read calls for a ipsl_desc
!  
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      TYPE (ipsl_desc), INTENT (inout)        :: this
      INTEGER, INTENT (in)                       :: iou
      CHARACTER (len=*), INTENT (in), OPTIONAL   :: prefix

!  this        ipsl_desc - this is what gets defined and written.
!  iou         i/o unit number of the netCDF file. 
!  prefix      character -  prefixed to variable names, so that netCDF
!                  doesn't have problems with repeated, identical names
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'ipsl_cdf_read_desc: '
      CHARACTER(len=32) :: prefix_use
      INTEGER(iprec), DIMENSION(3) :: dimlens
      
      CHARACTER (len=sn_len)    :: s_name                                 
      CHARACTER (len=ln_len)    :: l_name
      CHARACTER (len=sn_len)    :: ip_sname                                 
      CHARACTER (len=ln_len)    :: ip_lname
      CHARACTER (len=units_len) :: units                                 
      CHARACTER (len=30)        :: ipsl_type                                 
      LOGICAL                   :: l_ipbeam_def
      REAL(rprec)               :: sigma_default
      TYPE (ip_beam)            :: ipbeam

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

! Define the prefix to actually use
      IF (PRESENT(prefix)) THEN
         prefix_use = TRIM(ADJUSTL(prefix))
      ELSE
         prefix_use = ' '
      ENDIF

! Define all vn_--_use variable names
      CALL ipsl_cdf_defvn_desc(prefix_use)
         
! Read Scalar pieces
! Note: Read in to variables local to this subroutine.
      CALL cdf_read(iou, TRIM(vn_s_name_use), s_name)
      CALL cdf_read(iou, TRIM(vn_l_name_use), l_name)
      CALL cdf_read(iou, TRIM(vn_ip_sname_use), ip_sname)
      CALL cdf_read(iou, TRIM(vn_ip_lname_use), ip_lname)
      CALL cdf_read(iou, TRIM(vn_units_use), units)
      CALL cdf_read(iou, TRIM(vn_sigma_default_use), sigma_default)
      CALL cdf_read(iou, TRIM(vn_l_ipbeam_def_use),l_ipbeam_def)
      CALL cdf_read(iou, TRIM(vn_ipsl_type_use),ipsl_type)

!  Read Derived Type
      CALL intpol_cdf_read(ipbeam,iou,prefix)

      

! Create the ipsl_desc, this
      CALL ipsl_desc_construct(this,s_name,l_name,ip_sname,ip_lname,           &
     &   units, sigma_default,ipsl_type,ipbeam)

!  Destroy the local int/pol beamline ipbeam  to avoid memory leakage
      CALL ip_beam_destroy(ipbeam)
      
      RETURN
      
      END SUBROUTINE ipsl_cdf_read_desc
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------


!*******************************************************************************
! SECTION IV. AUXILIARY FUNCTIONS AND SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------
      SUBROUTINE ipsl_cdf_defvn_desc(prefix_use)
!  Subroutine to do define the character variable names for a ipsl_desc,
!  using the prefix. All the vn_ variables are module variables, and so do not
!   need to be declared here

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER (len=*), INTENT (in)   :: prefix_use

!  prefix_use      character -  prefixed to variable names, so that netCDF
!                  doesn't have problems with repeated, identical names

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'ipsl_cdf_defvn_desc: '

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

! Define all variable names
      vn_s_name_use = ipsl_cdf_mknam(prefix_use,vn_s_name)
      vn_l_name_use = ipsl_cdf_mknam(prefix_use,vn_l_name)
      vn_ip_sname_use = ipsl_cdf_mknam(prefix_use,vn_ip_sname)
      vn_ip_lname_use = ipsl_cdf_mknam(prefix_use,vn_ip_lname)
      vn_units_use = ipsl_cdf_mknam(prefix_use,vn_units)
      vn_sigma_default_use = ipsl_cdf_mknam(prefix_use,                        &
     &   vn_sigma_default)
      vn_l_ipbeam_def_use = ipsl_cdf_mknam(prefix_use,                         &
     &   vn_l_ipbeam_def)
      vn_ipsl_type_use = ipsl_cdf_mknam(prefix_use,                            &
     &   vn_ipsl_type)
      
      RETURN
      
      END SUBROUTINE ipsl_cdf_defvn_desc

!========================================================
      FUNCTION ipsl_cdf_mknam(c1,c2)
! A simple function to help in the generation of names

!-----------------------------------------------
!   F u n c t i o n   N a m e
!-----------------------------------------------
      CHARACTER(LEN=64) ipsl_cdf_mknam

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER(LEN=*), INTENT (in) :: c1,c2

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
      IF (LEN_TRIM(c1) .eq. 0) THEN
         ipsl_cdf_mknam = TRIM(c2)
      ELSE
         ipsl_cdf_mknam = ADJUSTL(TRIM(c1) // '_' // TRIM(c2))
      ENDIF
      
      RETURN
       
      END FUNCTION ipsl_cdf_mknam

!-----------------------------------------------
!-----------------------------------------------
!*******************************************************************************
! SECTION X.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
!
! JMS 2007-06-14. First version of ipsl_cdf.
!    Based on diagnostic_cdf
!


      END MODULE ipsl_cdf
