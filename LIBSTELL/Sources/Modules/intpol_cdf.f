      MODULE intpol_cdf
!--------------------------------------------------------------------------
! Module to Generically implement netCDF I/O for TYPE "int_pol"
!
!----------------------------------------------------------------------------

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

!      USE bsc
      USE ip_beamline  ! module containing int_pol TYPE definition
      USE ezcdf
      USE v3_utilities

!-------------------------------------------------------------------------------
!  Implicit None comes after USE statements, before other declarations
!-------------------------------------------------------------------------------
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Lengths of Character Variables
!-------------------------------------------------------------------------------
      INTEGER(iprec), PARAMETER, PRIVATE :: type_len=10      
      INTEGER(iprec), PARAMETER, PRIVATE :: sn_len=30      
      INTEGER(iprec), PARAMETER, PRIVATE :: ln_len=80      
      INTEGER(iprec), PARAMETER, PRIVATE :: units_len=30      
!SPH010408
      INTEGER, PRIVATE, PARAMETER        :: izero=0

!-------------------------------------------------------------------------------
!  Variable Names for netCDF. Make them Private.
!-------------------------------------------------------------------------------

!.......variable names for TYPE int_pol_coll..................................!
      CHARACTER (LEN=*), PRIVATE, PARAMETER ::                                 &
     &  vn_ipc_n_ip   = 'ip_coll_n_ip',                                       &
     &  vn_ipc_ip     = 'ip_coll_ip',                                         &
     &  vn_ipc_s_name = 'ip_coll_s_name',                                      &
     &  vn_ipc_l_name = 'ip_coll_l_name'

     
      CHARACTER (LEN=64), PRIVATE ::                                            &
     &  vn_ipc_n_ip_use,                                                        &
     &  vn_ipc_s_name_use,                                                      &
     &  vn_ipc_l_name_use                                                         


!.......variable names for TYPE int_pol.......................................!
      CHARACTER (LEN=*), PRIVATE, PARAMETER ::                                 &
     &  vn_nbeams = 'int_pol_nbeams',                                          &
     &  vn_ipbeam = 'int_pol_ipbeam',                                          &
     &  vn_s_name = 'int_pol_s_name',                                          &
     &  vn_l_name = 'int_pol_l_name'

     
      CHARACTER (LEN=64), PRIVATE ::                                           &
     &  vn_nbeams_use,                                                         &
     &  vn_s_name_use,                                                         &
     &  vn_l_name_use                                                         


!.......variable names for TYPE ip_beam.......................................!
!....note that q_unit & qdist are fns of q0vec & qfvec & hence don't need CDF I/O.  JS 9/20/05..!
      CHARACTER (LEN=*), PRIVATE, PARAMETER ::                                 &
!     &  vn_q_unit    = 'ip_beam_q_unit',                                      &
     &  vn_q0vec      = 'ip_beam_q0vec',                                       &
     &  vn_qfvec      = 'ip_beam_qfvec',                                       &
     &  vn_q0         = 'ip_beam_q0',                                          &
!     &  vn_qdist     = 'ip_beam_qdist',                                       &
     &  vn_wavelength = 'ip_beam_wavelength',                                  &
     &  vn_B_ratio    = 'ip_beam_B_ratio',                                     &
     &  vn_ipb_s_name = 'ip_beam_s_name',                                      &
     &  vn_ipb_l_name = 'ip_beam_l_name'

     
      CHARACTER (LEN=64), PRIVATE ::                                           &
!     &  vn_q_unit_use,                                                         &
     &  vn_q0vec_use,                                                          &
     &  vn_qfvec_use,                                                          &
     &  vn_q0_use,                                                             &         
!     &  vn_qdist_use,                                                          &
     &  vn_wavelength_use,                                                     &
     &  vn_B_ratio_use,                                                        &
     &  vn_ipb_s_name_use,                                                     &
     &  vn_ipb_l_name_use                                                         






!*******************************************************************************
! SECTION II. INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Generic Define
!-------------------------------------------------------------------------------
      INTERFACE intpol_cdf_define
         MODULE PROCEDURE intpol_cdf_define_int_pol,                               &
     &                    intpol_cdf_define_int_pol_coll,                          &
     &                    intpol_cdf_define_ip_beam
      END INTERFACE

!-------------------------------------------------------------------------------
!  Generic Write
!-------------------------------------------------------------------------------
      INTERFACE intpol_cdf_write
        MODULE PROCEDURE intpol_cdf_write_int_pol,                                &
     &                    intpol_cdf_write_int_pol_coll,                          &
     &                    intpol_cdf_write_ip_beam
      END INTERFACE

!-------------------------------------------------------------------------------
!  Generic Read
!-------------------------------------------------------------------------------
      INTERFACE intpol_cdf_read
        MODULE PROCEDURE intpol_cdf_read_int_pol,                                 &
     &                    intpol_cdf_read_int_pol_coll,                           &
     &                    intpol_cdf_read_ip_beam
      END INTERFACE


!-------------------------------------------------------------------------------

      CONTAINS
!*******************************************************************************
! SECTION III. DEFINE SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
      SUBROUTINE intpol_cdf_define_int_pol_coll(this,iou,prefix)
!  Subroutine to do the appropriate netCDF definition calls for a int_pol_coll

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      TYPE (int_pol_coll), INTENT (in)         :: this
      INTEGER, INTENT(in)                      :: iou
      CHARACTER (len=*), INTENT(in), OPTIONAL  :: prefix

!  this        int_pol_coll - this is what gets defined and written.
!  iou         i/o unit number of the netCDF file. 
!  prefix      character -  prefixed to variable names, so that netCDF
!                  doesn't have problems with repeated, identical names
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER                     :: i, n_ip
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'intpol_cdf_define_int_pol_coll: '
      CHARACTER(len=12)  :: mychar
      CHARACTER(len=32)  :: prefix_use, ip_prefix

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

!      write(*,*) ' now entering ', sub_name
! Define the prefix to actually use
      IF (PRESENT(prefix)) THEN
         prefix_use = TRIM(ADJUSTL(prefix))
      ELSE
         prefix_use = ' '
      ENDIF

! Define all vn_--_use variable names
      CALL intpol_cdf_defvn_int_pol_coll(prefix_use)

! Define Components

      CALL cdf_define(iou, TRIM(vn_ipc_n_ip_use),  this % n_ip)
      CALL cdf_define(iou, TRIM(vn_ipc_s_name_use), this % s_name)
      CALL cdf_define(iou, TRIM(vn_ipc_l_name_use), this % l_name)

!.......define each ip array element using the define_int_pol routine.......!
      n_ip = this % n_ip
!      write(*,*) sub_name, 'n_ip = ', n_ip
      
!......***temporarily forced ONE iteration for debugging.  JS 10/7/05...!
      do i = 1,n_ip
!      do i = 1,1

!.......pass the array index info using the optional "prefix" character variable...!
!.....(this ensures that the "vn" names for each array element are unique).........!
        write(mychar, '(i8)' ) i
        if ( PRESENT(prefix) ) then
          ip_prefix = TRIM(prefix_use)//'_ip'//TRIM(ADJUSTL(mychar))
        else
          ip_prefix = 'ip'//TRIM(ADJUSTL(mychar))
        end if

!        write(*,*) sub_name, 'prefix_use = ', prefix_use
!        write(*,*) sub_name, 'ip_prefix = ', ip_prefix
        CALL intpol_cdf_define_int_pol(this%ip(i), iou, ip_prefix)
      end do
      
!      write(*,*) ' now exiting ', sub_name

      RETURN
      
      END SUBROUTINE intpol_cdf_define_int_pol_coll

!========================================================================
      SUBROUTINE intpol_cdf_defvn_int_pol_coll(prefix_use)
!  Subroutine to do define the character variable names for a int_pol,
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
     &  'intpol_cdf_defvn_int_pol_coll: '

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
!      write(*,*) ' now entering ', sub_name

! Define all variable names
      vn_ipc_n_ip_use = intpol_cdf_mknam(prefix_use,vn_ipc_n_ip)
      vn_ipc_s_name_use = intpol_cdf_mknam(prefix_use,vn_ipc_s_name)
      vn_ipc_l_name_use = intpol_cdf_mknam(prefix_use,vn_ipc_l_name)

!      write(*,*) ' now exiting ', sub_name

      RETURN
      
      END SUBROUTINE intpol_cdf_defvn_int_pol_coll




!=============================================================
      SUBROUTINE intpol_cdf_read_int_pol_coll(this,iou,prefix)
!---------------------------------------------------------------------
!  Subroutine to do the appropriate netCDF read calls for int_pol
!
!
!---------------------------------------------------------------------
!   D u m m y   A r g u m e n t s
!---------------------------------------------------------------------
      TYPE (int_pol_coll), INTENT (inout)        :: this
      INTEGER, INTENT (in)                       :: iou
      CHARACTER (len=*), INTENT (in), OPTIONAL   :: prefix

!  this        int_pol - this is what gets defined and written.
!  iou         i/o unit number of the netCDF file. 
!  prefix      character -  prefixed to variable names, so that netCDF
!                  doesn't have problems with repeated, identical names
!-------------------------------------------------------------------
!   L o c a l   V a r i a b l e s
!-------------------------------------------------------------------
      CHARACTER(len=32)            :: prefix_use
      INTEGER(iprec), DIMENSION(3) :: DIM1, DIM2
      
      INTEGER                                    :: n_ip
!SPH010408      INTEGER(iprec)                             :: n_ip
      TYPE(int_pol), DIMENSION(:), ALLOCATABLE   :: ip
      CHARACTER (len=30)                         :: s_name                                 
      CHARACTER (len=80)                         :: l_name


! Declare local variables
      INTEGER(iprec)           :: i, arr_size
      INTEGER                  :: ierr1, ierr2, ierr3
!SPH010408      INTEGER(iprec)           :: ierr1, ierr2, ierr3
      CHARACTER(len=12)        :: mychar
      CHARACTER(len=32)        :: ip_prefix
      CHARACTER(len=*), PARAMETER  :: sub_name =                               &
     &  'intpol_cdf_read_int_pol_coll: '

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
!      write(*,*) 'now entering ', sub_name

! Define the prefix to actually use
      IF (PRESENT(prefix)) THEN
         prefix_use = TRIM(ADJUSTL(prefix))
      ELSE
         prefix_use = ' '
      ENDIF

!........Define all vn_--_use variable names..........................!
      CALL intpol_cdf_defvn_int_pol_coll(prefix_use)
         
! Read Components
! Note: Read in to variables local to this subroutine.
! Arrays require inquiry regarding size, and allocation before actual reading.
      CALL cdf_read(iou, TRIM(vn_ipc_n_ip_use), n_ip)
      CALL cdf_read(iou, TRIM(vn_ipc_s_name_use), s_name)
      CALL cdf_read(iou, TRIM(vn_ipc_l_name_use), l_name)

!      write(*,*) sub_name, 'READ IN n_ip = ', n_ip 
!      write(*,*) sub_name, 'READ IN s_name = ', s_name
!      write(*,*) sub_name, 'READ IN l_name = ', l_name

!......now read in the ip derived TYPE array, one element at a time.......!

!......use the read-in "n_ip" value to correctly size the allocated array....!
      if ( n_ip > 0 ) then
        arr_size = n_ip
      else
        arr_size = 0
        write(*,*) sub_name//'Error in reading ip array size'
        write(*,*) sub_name//'read in n_ip = ', n_ip
      end if

      ALLOCATE( ip(arr_size),STAT=ierr1 )
      CALL assert_eq(izero,ierr1,sub_name // 'ip array')

!.......read in each ip element, using the array index to produce a unique "vn name"...!
      do i = 1, arr_size

        write(mychar, '(i8)' ) i
        if ( PRESENT(prefix) ) then
          ip_prefix = TRIM(prefix_use)//'_ip'//TRIM(ADJUSTL(mychar))
        else
          ip_prefix = 'ip'//TRIM(ADJUSTL(mychar))
        end if

        CALL intpol_cdf_read_int_pol(ip(i),iou, ip_prefix)
      end do


!........construct the int_pol derived TYPE using the now-filled local variables...!
      CALL int_pol_coll_construct(this, n_ip, ip, s_name, l_name)



!.......deallocate local allocatable array..........................................!
      DEALLOCATE(ip,STAT=ierr1)

!      DIM1(1:3) = SHAPE(this%ip(1)%ipbeam(1)%B_ratio)
!      write(*,*) sub_name, 'DIM1 = ', DIM1
!      write(*,*) sub_name, 'B_ratio = ', this%ip(1)%ipbeam(1)%B_ratio
!      write(*,*) 'Now exiting ', sub_name
      RETURN
      
      END SUBROUTINE intpol_cdf_read_int_pol_coll


!==================================================================================
      SUBROUTINE intpol_cdf_write_int_pol_coll(this,iou,prefix)
!  Subroutine to do the appropriate netCDF write calls for a int_pol_coll
!
!------------------------------------------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------------------------------------
      TYPE (int_pol_coll), INTENT (in)         :: this
      INTEGER, INTENT(in)                      :: iou
      CHARACTER (len=*), INTENT(in), OPTIONAL  :: prefix

!  this        int_pol_coll - this is what gets defined and written.
!  iou         i/o unit number of the netCDF file. 
!  prefix      character -  prefixed to variable names, so that netCDF
!                  doesn't have problems with repeated, identical names
!-----------------------------------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------------------------------
      INTEGER(iprec)           :: i, arr_size, ip_size
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'intpol_cdf_write_int_pol_coll: '
      CHARACTER(len=12)        :: mychar
      CHARACTER(len=32)        :: ip_prefix
      CHARACTER(len=32)        :: prefix_use

!----------------------------------------------------------------------
!  Start of Executable Code
!---------------------------------------------------------------------
!      write(*,*) 'Now entering ', sub_name

! Define the prefix to actually use
      IF (PRESENT(prefix)) THEN
         prefix_use = TRIM(ADJUSTL(prefix))
      ELSE
         prefix_use = ' '
      ENDIF

! Define all vn_--_use variable names
      CALL intpol_cdf_defvn_int_pol_coll(prefix_use)
         
! Write Components
      CALL cdf_write(iou, TRIM(vn_ipc_n_ip_use), this % n_ip)
      CALL cdf_write(iou, TRIM(vn_ipc_s_name_use), this % s_name)
      CALL cdf_write(iou, TRIM(vn_ipc_l_name_use), this % l_name)

      ip_size = SIZE( this %ip)
!      write(*,*) sub_name, 'ip_size = ', ip_size
!      write(*,*) sub_name, 'n_ip = ', this % n_ip

!........use variable "n_ip" to determine the size of the ip() array...........!
      if ( ip_size > 0 .AND. ip_size .EQ. this%n_ip) then
        arr_size = ip_size
      else
        arr_size = 0
        write(*,*) sub_name//'Error in ip array size'
      end if

!      write(*,*) sub_name, 'arr_size = ', arr_size

!.......read in each ip element, using the array index to produce a unique "vn name"...!
!......***temporarily forced ONE iteration for debugging.  JS 10/7/05...!
      do i = 1, arr_size
!      do i = 1, 1

        write(mychar, '(i8)' ) i
        if ( PRESENT(prefix) ) then
          ip_prefix = TRIM(prefix_use)//'_ip'//TRIM(ADJUSTL(mychar))
        else
          ip_prefix = 'ip'//TRIM(ADJUSTL(mychar))
        end if

        CALL intpol_cdf_write_int_pol(this%ip(i), iou, ip_prefix)
      end do

      
      RETURN
      
      END SUBROUTINE intpol_cdf_write_int_pol_coll






!-------------------------------------------------------------------------------
      SUBROUTINE intpol_cdf_define_int_pol(this,iou,prefix)
!  Subroutine to do the appropriate netCDF definition calls for a int_pol

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      TYPE (int_pol), INTENT (in)              :: this
      INTEGER, INTENT(in)                      :: iou
      CHARACTER (len=*), INTENT(in), OPTIONAL  :: prefix

!  this        int_pol - this is what gets defined and written.
!  iou         i/o unit number of the netCDF file. 
!  prefix      character -  prefixed to variable names, so that netCDF
!                  doesn't have problems with repeated, identical names
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER                     :: i, nbeams
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'intpol_cdf_define_int_pol: '
      CHARACTER(len=12)  :: mychar
      CHARACTER(len=32)  :: prefix_use, ipb_prefix

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

!      write(*,*) ' now entering ', sub_name

! Define the prefix to actually use
      IF (PRESENT(prefix)) THEN
         prefix_use = TRIM(ADJUSTL(prefix))
      ELSE
         prefix_use = ' '
      ENDIF

! Define all vn_--_use variable names
      CALL intpol_cdf_defvn_int_pol(prefix_use)

! Define Components

      CALL cdf_define(iou, TRIM(vn_nbeams_use),  this % nbeams)
      CALL cdf_define(iou, TRIM(vn_s_name_use), this % s_name)
      CALL cdf_define(iou, TRIM(vn_l_name_use), this % l_name)

!.......define each int_pol array element using the define_int_pol routine.......!
      nbeams = this % nbeams
!      write(*,*) sub_name, 'nbeams = ', nbeams
      do i = 1,nbeams

!.......pass the array index info using the optional "prefix" character variable...!
!.....(this ensures that the "vn" names for each array element are unique).........!
        write(mychar, '(i8)' ) i
        if ( PRESENT(prefix) ) then
          ipb_prefix = TRIM(prefix_use)//'_ipb'//TRIM(ADJUSTL(mychar))
        else
          ipb_prefix = 'ipb'//TRIM(ADJUSTL(mychar))
        end if

!        write(*,*) sub_name, 'prefix_use = ', prefix_use
!        write(*,*) sub_name, 'ipb_prefix = ', ipb_prefix
        CALL intpol_cdf_define_ip_beam(this%ipbeam(i), iou, ipb_prefix)
      end do
      
!      write(*,*) ' now exiting ', sub_name

      RETURN
      
      END SUBROUTINE intpol_cdf_define_int_pol



!========================================================================
      SUBROUTINE intpol_cdf_defvn_int_pol(prefix_use)
!  Subroutine to do define the character variable names for a int_pol,
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
     &  'intpol_cdf_defvn_int_pol: '

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
!      write(*,*) ' now entering ', sub_name

! Define all variable names
      vn_nbeams_use = intpol_cdf_mknam(prefix_use,vn_nbeams)
      vn_s_name_use = intpol_cdf_mknam(prefix_use,vn_s_name)
      vn_l_name_use = intpol_cdf_mknam(prefix_use,vn_l_name)

!      write(*,*) ' now exiting ', sub_name

      RETURN
      
      END SUBROUTINE intpol_cdf_defvn_int_pol


!=============================================================
      SUBROUTINE intpol_cdf_read_int_pol(this,iou,prefix)
!---------------------------------------------------------------------
!  Subroutine to do the appropriate netCDF read calls for int_pol
!
!
!---------------------------------------------------------------------
!   D u m m y   A r g u m e n t s
!---------------------------------------------------------------------
      TYPE (int_pol), INTENT (inout)             :: this
      INTEGER, INTENT (in)                       :: iou
      CHARACTER (len=*), INTENT (in), OPTIONAL   :: prefix

!  this        int_pol - this is what gets defined and written.
!  iou         i/o unit number of the netCDF file. 
!  prefix      character -  prefixed to variable names, so that netCDF
!                  doesn't have problems with repeated, identical names
!-------------------------------------------------------------------
!   L o c a l   V a r i a b l e s
!-------------------------------------------------------------------
      CHARACTER(len=32)            :: prefix_use
      INTEGER(iprec), DIMENSION(3) :: DIM1, DIM2
      
      INTEGER                                    :: nbeams      
!SPH010408      INTEGER(iprec)                             :: nbeams
      TYPE(ip_beam), DIMENSION(:), ALLOCATABLE   :: ipbeam
      CHARACTER (len=30)                         :: s_name                                 
      CHARACTER (len=80)                         :: l_name


! Declare local variables
      INTEGER(iprec)           :: i, arr_size
      INTEGER                  :: ierr1, ierr2, ierr3
!SPH010408      INTEGER(iprec)           :: ierr1, ierr2, ierr3
      CHARACTER(len=12)        :: mychar
      CHARACTER(len=32)        :: ipb_prefix
      CHARACTER(len=*), PARAMETER  :: sub_name =                               &
     &  'intpol_cdf_read_int_pol: '

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
      write(*,*) ' now entering ', sub_name

! Define the prefix to actually use
      IF (PRESENT(prefix)) THEN
         prefix_use = TRIM(ADJUSTL(prefix))
      ELSE
         prefix_use = ' '
      ENDIF

!........Define all vn_--_use variable names..........................!
      CALL intpol_cdf_defvn_int_pol(prefix_use)
         
! Read Components
! Note: Read in to variables local to this subroutine.
! Arrays require inquiry regarding size, and allocation before actual reading.
      CALL cdf_read(iou, TRIM(vn_nbeams_use), nbeams)
      CALL cdf_read(iou, TRIM(vn_s_name_use), s_name)
      CALL cdf_read(iou, TRIM(vn_l_name_use), l_name)

!      write(*,*) sub_name, 'READ IN nbeams = ', nbeams 
!      write(*,*) sub_name, 'READ IN s_name = ', s_name
!      write(*,*) sub_name, 'READ IN l_name = ', l_name

!......now read in the ipbeam derived TYPE array, one element at a time.......!

!......use the read-in "nbeams" value to correctly size the allocated array....!
      if ( nbeams > 0 ) then
        arr_size = nbeams
      else
        arr_size = 0
        write(*,*) sub_name//'Error in reading ipbeam array size'
        write(*,*) sub_name//'read in nbeams = ', nbeams
      end if

      ALLOCATE( ipbeam(arr_size),STAT=ierr1 )
      CALL assert_eq(izero,ierr1,sub_name // 'ipbeam array')

!.......read in each ipbeam element, using the array index to produce a unique "vn name"...!
      do i = 1, arr_size

        write(mychar, '(i8)' ) i
        if ( PRESENT(prefix) ) then
          ipb_prefix = TRIM(prefix_use)//'_ipb'//TRIM(ADJUSTL(mychar))
        else
          ipb_prefix = 'ipb'//TRIM(ADJUSTL(mychar))
        end if

        CALL intpol_cdf_read_ip_beam(ipbeam(i),iou, ipb_prefix)
      end do

!      DIM1(1:3) = SHAPE(ipbeam(1)%B_ratio)
!      write(*,*) sub_name, 'DIM1 = ', DIM1
!      write(*,*) sub_name, 'nbeams = ', nbeams

!........construct the int_pol derived TYPE using the now-filled local variables...!
      CALL int_pol_construct(this, nbeams, ipbeam, s_name, l_name)

!.......deallocate local allocatable array..........................................!
      DEALLOCATE(ipbeam,STAT=ierr1)

      RETURN
      
      END SUBROUTINE intpol_cdf_read_int_pol




!==================================================================================
      SUBROUTINE intpol_cdf_write_int_pol(this,iou,prefix)
!  Subroutine to do the appropriate netCDF write calls for a int_pol
!
!------------------------------------------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------------------------------------
      TYPE (int_pol), INTENT (in)              :: this
      INTEGER, INTENT(in)                      :: iou
      CHARACTER (len=*), INTENT(in), OPTIONAL  :: prefix

!  this        int_pol - this is what gets defined and written.
!  iou         i/o unit number of the netCDF file. 
!  prefix      character -  prefixed to variable names, so that netCDF
!                  doesn't have problems with repeated, identical names
!-----------------------------------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------------------------------
      INTEGER(iprec)           :: i, arr_size, ipb_size
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'intpol_cdf_write_int_pol: '
      CHARACTER(len=12)        :: mychar
      CHARACTER(len=32)        :: ipb_prefix
      CHARACTER(len=32)        :: prefix_use

!----------------------------------------------------------------------
!  Start of Executable Code
!---------------------------------------------------------------------
!      write(*,*) ' now entering ', sub_name

! Define the prefix to actually use
      IF (PRESENT(prefix)) THEN
         prefix_use = TRIM(ADJUSTL(prefix))
      ELSE
         prefix_use = ' '
      ENDIF

! Define all vn_--_use variable names
      CALL intpol_cdf_defvn_int_pol(prefix_use)
         
! Write Components
      CALL cdf_write(iou, TRIM(vn_nbeams_use), this % nbeams)
      CALL cdf_write(iou, TRIM(vn_s_name_use), this % s_name)
      CALL cdf_write(iou, TRIM(vn_l_name_use), this % l_name)

      ipb_size = SIZE(this % ipbeam)

!........use variable "nbeams" to determine the size of the ipbeam() array...........!
      if ( ipb_size > 0 .AND. ipb_size .EQ. this%nbeams) then
        arr_size = ipb_size
      else
        arr_size = 0
        write(*,*) sub_name//'Error in ipbeam array size'
      end if


!.......read in each ipbeam element, using the array index to produce a unique "vn name"...!
      do i = 1, arr_size

        write(mychar, '(i8)' ) i
        if ( PRESENT(prefix) ) then
          ipb_prefix = TRIM(prefix_use)//'_ipb'//TRIM(ADJUSTL(mychar))
        else
          ipb_prefix = 'ipb'//TRIM(ADJUSTL(mychar))
        end if

        CALL intpol_cdf_write_ip_beam(this%ipbeam(i), iou, ipb_prefix)
      end do

      
      RETURN
      
      END SUBROUTINE intpol_cdf_write_int_pol


!-------------------------------------------------------------------------------
      SUBROUTINE intpol_cdf_define_ip_beam(this,iou,prefix)
!  Subroutine to do the appropriate netCDF definition calls for a ip_beam

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      TYPE (ip_beam), INTENT (in)              :: this
      INTEGER, INTENT(in)                      :: iou
      CHARACTER (len=*), INTENT(in), OPTIONAL  :: prefix

!  this        ip_beam - this is what gets defined and written.
!  iou         i/o unit number of the netCDF file. 
!  prefix      character -  prefixed to variable names, so that netCDF
!                  doesn't have problems with repeated, identical names
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'intpol_cdf_define_ip_beam: '
      CHARACTER(len=32) :: prefix_use

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

!      write(*,*) 'now entering ', sub_name

! Define the prefix to actually use
      IF (PRESENT(prefix)) THEN
         prefix_use = TRIM(ADJUSTL(prefix))
      ELSE
         prefix_use = ' '
      ENDIF

! Define all vn_--_use variable names
      CALL intpol_cdf_defvn_ip_beam(prefix_use)
         
! Define Components

      CALL cdf_define(iou, TRIM(vn_q0vec_use),  this % q0vec)
      CALL cdf_define(iou, TRIM(vn_qfvec_use), this % qfvec)
      CALL cdf_define(iou, TRIM(vn_q0_use), this % q0)


      CALL cdf_define(iou, TRIM(vn_wavelength_use), this % wavelength)
      CALL cdf_define(iou, TRIM(vn_B_ratio_use), this % B_ratio)
      CALL cdf_define(iou, TRIM(vn_ipb_s_name_use), this % s_name)
      CALL cdf_define(iou, TRIM(vn_ipb_l_name_use), this % l_name)

!.......qdist & q_unit are functions of q0vec & qfvec & hence don't get CDF I/O...!
!      CALL cdf_define(iou, TRIM(vn_q_unit_use), this % q_unit)
!      CALL cdf_define(iou, TRIM(vn_qdist_use), this % qdist)

!      write(*,*) 'now exiting ', sub_name
      
      RETURN
      
      END SUBROUTINE intpol_cdf_define_ip_beam




!========================================================================
      SUBROUTINE intpol_cdf_defvn_ip_beam(prefix_use)
!  Subroutine to do define the character variable names for a ip_beam,
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
     &  'intpol_cdf_defvn_ip_beam: '

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
!      write(*,*) 'now entering ', sub_name

! Define all variable names


      vn_q0vec_use = intpol_cdf_mknam(prefix_use,vn_q0vec)
      vn_qfvec_use = intpol_cdf_mknam(prefix_use,vn_qfvec)
      vn_q0_use = intpol_cdf_mknam(prefix_use,vn_q0)

      vn_wavelength_use = intpol_cdf_mknam(prefix_use,vn_wavelength)
      vn_B_ratio_use = intpol_cdf_mknam(prefix_use,vn_B_ratio)
      vn_ipb_s_name_use = intpol_cdf_mknam(prefix_use,vn_ipb_s_name)
      vn_ipb_l_name_use = intpol_cdf_mknam(prefix_use,vn_ipb_l_name)

!.......qdist & q_unit are functions of q0vec & qfvec & hence don't get CDF I/O...!
!      vn_qdist_use = intpol_cdf_mknam(prefix_use,vn_qdist)
!      vn_q_unit_use = intpol_cdf_mknam(prefix_use,vn_q_unit)
      
      RETURN
      
      END SUBROUTINE intpol_cdf_defvn_ip_beam

!====================================================================
      SUBROUTINE intpol_cdf_read_ip_beam(this,iou,prefix)
!---------------------------------------------------------------------
!  Subroutine to do the appropriate netCDF read calls for ip_beam
!
!---------------------------------------------------------------------
!   D u m m y   A r g u m e n t s
!---------------------------------------------------------------------
      TYPE (ip_beam), INTENT (inout)             :: this
      INTEGER, INTENT (in)                       :: iou
      CHARACTER (len=*), INTENT (in), OPTIONAL   :: prefix

!  this        ip_beam - this is what gets defined and written.
!  iou         i/o unit number of the netCDF file. 
!  prefix      character -  prefixed to variable names, so that netCDF
!                  doesn't have problems with repeated, identical names
!-------------------------------------------------------------------
!   L o c a l   V a r i a b l e s
!-------------------------------------------------------------------
      CHARACTER(len=32)           :: prefix_use
      INTEGER, DIMENSION(3)       :: dimlens, DIM1
!SPH010408      INTEGER(iprec), DIMENSION(3) :: dimlens, DIM1
      

      REAL(rprec), DIMENSION(3)   :: q_unit
      REAL(rprec), DIMENSION(3)   :: q0vec
      REAL(rprec), DIMENSION(3)   :: qfvec
      REAL(rprec)                 :: q0
      REAL(rprec)                 :: qdist
      REAL(rprec)                 :: wavelength

      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE   :: B_ratio

      CHARACTER (len=30)          :: s_name                                 
      CHARACTER (len=80)          :: l_name


! Declare local variables
      INTEGER(iprec)           :: q_size, bsize1, bsize2, bsize3    
      INTEGER                  :: ierr1, ierr2, ierr3
!SPH010408      INTEGER(iprec)           :: ierr1, ierr2, ierr3

      CHARACTER(len=*), PARAMETER  :: sub_name =                               &
     &  'intpol_cdf_read_ip_beam: '

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

!      write(*,*) 'now entering ', sub_name


! Define the prefix to actually use
      IF (PRESENT(prefix)) THEN
         prefix_use = TRIM(ADJUSTL(prefix))
      ELSE
         prefix_use = ' '
      ENDIF

! Define all vn_--_use variable names
      CALL intpol_cdf_defvn_ip_beam(prefix_use)
         
! Read Components
! Note: Read in to variables local to this subroutine.
! Arrays require inquiry regarding size, and allocation before actual reading.




!.......qdist not read in, since it is CALCULATED in ip_beam_construct!.....!
!      CALL cdf_inquire(iou, TRIM(vn_q_unit_use),dimlens)
!      q_size = dimlens(1)
!      CALL assert_eq(0,dimlens(2),dimlens(3),                                  &
!     &   sub_name // 'Unexpected q_unit dimensions')
!      CALL assert_eq(0, (q_size - 3),                                          &
!     &   sub_name // 'Unexpected q_unit length')
!      CALL cdf_read(iou, TRIM(vn_q_unit_use), q_unit)


      CALL cdf_inquire(iou, TRIM(vn_q0vec_use),dimlens)
      q_size = dimlens(1)
      CALL assert_eq(izero,dimlens(2),dimlens(3),                              &
     &   sub_name // 'Unexpected q0vec dimensions')
      CALL assert_eq(izero, (q_size - 3),                                      &
     &   sub_name // 'Unexpected q0vec length')
      CALL cdf_read(iou, TRIM(vn_q0vec_use), q0vec)


      CALL cdf_inquire(iou, TRIM(vn_qfvec_use),dimlens)
      q_size = dimlens(1)
      CALL assert_eq(izero,dimlens(2),dimlens(3),                               &
     &   sub_name // 'Unexpected qfvec dimensions')
      CALL assert_eq(izero, (q_size - 3),                                       &
     &   sub_name // 'Unexpected qfvec length')
      CALL cdf_read(iou, TRIM(vn_qfvec_use), qfvec)


      CALL cdf_read(iou, TRIM(vn_q0_use), q0)

!.......qdist not read in, since it is CALCULATED in ip_beam_construct!
!      CALL cdf_read(iou, TRIM(vn_qdist_use),qdist)

      CALL cdf_read(iou, TRIM(vn_wavelength_use),wavelength)


      CALL cdf_inquire(iou, TRIM(vn_B_ratio_use),dimlens)
      bsize1 = dimlens(1)
      bsize2 = dimlens(2)                 
      bsize3 = dimlens(3)                 
      ALLOCATE( B_ratio(bsize1, bsize2, bsize3),STAT=ierr1 )
      CALL assert_eq(izero,ierr1,sub_name // 'B_ratio')
      CALL cdf_read(iou, TRIM(vn_B_ratio_use), B_ratio)



      CALL cdf_read(iou, TRIM(vn_ipb_s_name_use), s_name)
      CALL cdf_read(iou, TRIM(vn_ipb_l_name_use), l_name)




! Create the ip_beam,"this", using the now-filled local variables
      CALL ip_beam_construct(this, q0vec, qfvec, q0, wavelength,                &
     &                             B_ratio, s_name, l_name)
     
      DIM1(1:3) = SHAPE(this%B_ratio)

!.......deallocate local allocatable array..........................................!
      DEALLOCATE(B_ratio,STAT=ierr1)

      RETURN
      
      END SUBROUTINE intpol_cdf_read_ip_beam




!==================================================================================
      SUBROUTINE intpol_cdf_write_ip_beam(this,iou,prefix)
!  Subroutine to do the appropriate netCDF write calls for a ip_beam
!
!------------------------------------------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------------------------------------
      TYPE (ip_beam), INTENT (in)           :: this
      INTEGER, INTENT(in)                      :: iou
      CHARACTER (len=*), INTENT(in), OPTIONAL  :: prefix

!  this        ip_beam - this is what gets defined and written.
!  iou         i/o unit number of the netCDF file. 
!  prefix      character -  prefixed to variable names, so that netCDF
!                  doesn't have problems with repeated, identical names
!-----------------------------------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------------------------------
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'intpol_cdf_write_ip_beam: '
      CHARACTER(len=32) :: prefix_use

!----------------------------------------------------------------------
!  Start of Executable Code
!---------------------------------------------------------------------
!      write(*,*) 'now entering ', sub_name


! Define the prefix to actually use
      IF (PRESENT(prefix)) THEN
         prefix_use = TRIM(ADJUSTL(prefix))
      ELSE
         prefix_use = ' '
      ENDIF

! Define all vn_--_use variable names
      CALL intpol_cdf_defvn_ip_beam(prefix_use)
         
! Write Components

      CALL cdf_write(iou, TRIM(vn_q0vec_use), this % q0vec)
      CALL cdf_write(iou, TRIM(vn_qfvec_use), this % qfvec)
      CALL cdf_write(iou, TRIM(vn_q0_use), this % q0)
      CALL cdf_write(iou, TRIM(vn_wavelength_use), this % wavelength)
      CALL cdf_write(iou, TRIM(vn_B_ratio_use), this % B_ratio)
      CALL cdf_write(iou, TRIM(vn_ipb_s_name_use), this % s_name)
      CALL cdf_write(iou, TRIM(vn_ipb_l_name_use), this % l_name)

!....note that q_unit & qdist are fns of q0vec & qfvec & hence don't need CDF I/O.  JS 9/20/05..!
!      CALL cdf_write(iou, TRIM(vn_q_unit_use), this % q_unit)
!      CALL cdf_write(iou, TRIM(vn_qdist_use), this % qdist)

      
      RETURN
      
      END SUBROUTINE intpol_cdf_write_ip_beam



!==========================================================================
      FUNCTION intpol_cdf_mknam(c1,c2)
!--------------------------------------------------------------
! A simple function to help in the generation of names.  Specifically: it concatenates
!    input strings c1 and c2 with a '_' between them
!-------------------------------------------------------------

!-----------------------------------------------
!   F u n c t i o n   N a m e
!-----------------------------------------------
      CHARACTER(LEN=64) intpol_cdf_mknam

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER(LEN=*), INTENT (in) :: c1,c2

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
      IF (LEN_TRIM(c1) .eq. 0) THEN
         intpol_cdf_mknam = TRIM(c2)
      ELSE
         intpol_cdf_mknam = ADJUSTL(TRIM(c1) // '_' // TRIM(c2))
      ENDIF
      
      RETURN
       
      END FUNCTION intpol_cdf_mknam


      END module intpol_cdf

