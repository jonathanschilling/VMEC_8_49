!*******************************************************************************
!  File diagnostic_dot.f
!  Contains module diagnostic_dot
!    Module for opening and reading a 'diagnostic.' file, and then placing the coil
!    information into a bsc_coilcoll.
!  
!-------------------------------------------------------------------------------
!   DEPENDENCIES
!-------------------------------------------------------------------------------
!
!    This module uses the following modules:
!       stel_kinds
!       stel_constants
!       bsc
!       safe_open_mod
!       mddc_T
!
!-------------------------------------------------------------------------------
!   CHANGE HISTORY
!-------------------------------------------------------------------------------
!
!    JDH. 11.15.2002
!    JDH 11.21.2002 - cleaned up determination of keyword. Added 'end'
!    JDH 12.04.2002 
!          Now using modules st_type and safe_open_mod
!    JDH 12.05.2002 
!          Replaced st_type with stel_kinds and stel_constants
!     JDH 01.20.2003
!          Changed raux from 1 / (2 Pi) to 1, for flux_loop_circular
!     JDH 03.28.2003
!          Added rogowski coils. Changed 'end' to 'end_of_file'
!     JDH 09-10-2004
!          Added coding for writing of diagnostic_desc, in addition to bsc_coils
!     JDH 12-13-2004
!          Eliminated pointer assignments of coils to diagnostic_desc's. No
!          longer needed because of changes in diagnostic_T.
!  JDH 2007-06-11
!     Replaced diagnostic_T stuff with mddc_T
!  JDH 2008-05-16
!     Eliminated (iprec)
!
!  JDH 2012-01-12
!     No longer use raux. Added b_rogowski, i_rogowski, f_rogowski
!     so that signal has units of B (Tesla), I (Ampere), or flux (Weber)
!
!-------------------------------------------------------------------------------
!   USAGE
!-------------------------------------------------------------------------------
!     Usage: Here is a sample of code that shows how to use this module
!
!      USE bsc_T
!      USE diagnostic_dot
!      TYPE(bsc_coilcoll) :: diagnostic_coils
!      TYPE(mddc_desc), DIMENSION(1000) :: mddc_descriptions_a
! ...
! !  Read in diagnostic coils 
!      CALL diagnostic_dot_read(file_name_diagnostic_dot,diagnostic_coils[,
!       &    mddc_descriptions_a])
!
!-------------------------------------------------------------------------------
!   FILE STRUCTURE - diagnostic. file
!-------------------------------------------------------------------------------
!
!  diagnostic. File Structure:
! line 1: Character - identifier
!
!      (0 or more lines, no keyword present. All ignored by parser.)
!    keyword
!      Character - coil identifier
!      lines containing data - number depends on keyword
!
!      (0 or more lines, no keyword present. All ignored by parser.)
!    keyword
!      Character - coil identifier
!      lines containing data - number depends on keyword
!
!    ...
!
!    end
!      (0 or more lines, ignored by parser.)!    
!    
!-------------------------------------------------------------------------------
!   COMMENTS
!-------------------------------------------------------------------------------
!
! 1)  Defined keywords are:
!       dd_keyword(1) = 'flux_loop'
!       dd_keyword(2) = 'flux_loop_circular'
!       dd_keyword(3) = 'magnetic_probe' 
!       dd_keyword(4) = 'magnetic_probe_tokamak' 
!       dd_keyword(5) = 'rogowski' 
!       dd_keyword(6) = 'b_rogowski' 
!       dd_keyword(7) = 'i_rogowski' 
!       dd_keyword(8) = 'f_rogowski' 
!       dd_keyword(9) = 'end_of_file' 
!
!*******************************************************************************

!*******************************************************************************
!  MODULE diagnostic_dot
!    
! SECTION I. VARIABLE DECLARATIONS
! SECTION II. INTERFACE BLOCKS
! SECTION III. PUBLICLY ACCESSIBLE SUBROUTINES
! SECTION IV. PARSING SUBROUTINES
! SECTION V. AUXILIARY SUBROUTINES
!*******************************************************************************
      MODULE diagnostic_dot

!*******************************************************************************
! SECTION I. VARIABLE DECLARATIONS
!*******************************************************************************

      USE stel_kinds
      USE stel_constants
      USE bsc_T
      USE mddc_T
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Variables local to the module. Make them private.
!-------------------------------------------------------------------------------
      INTEGER ::  dd_iou
      INTEGER, PARAMETER :: n_xnode=10000
      INTEGER :: i_line
      CHARACTER(len=200) :: line
      LOGICAL :: dd_l_done
      TYPE (bsc_coil), SAVE,TARGET :: coil_temp
      TYPE (mddc_desc), SAVE :: mddc_desc_temp
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: xnode

      PRIVATE dd_iou, n_xnode, i_line 
      PRIVATE line, dd_l_done, coil_temp, mddc_desc_temp, xnode

!  coil_temp       A temporary bsc_coil. Lahey Fujitsu compiler tells me it MUST have
!                     the SAVE attribute, since 'component initialization was specified
!                     in a module'. JDH 11.18.02
!  mddc_desc_temp  A temporary mddc_desc. Constructed in the individual parsing
!                    subroutines, and assigned (if necessary) in the SELECT CASE 
!                    construct in diagnostic_dot_read  
!  dd_iou          i/o unit number for diagnostic_dot file
!  dd_l_done       logical - true if done reading file
!  i_line          integer - line counter for file
!  line            character - the last line read from the file
!  n_xnode         integer - size for allocation of xnode array.
!  xnode           real array, (3,n_xnode) - positions of coil nodes
          
!*******************************************************************************
! SECTION II. INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!  
!-------------------------------------------------------------------------------

      CONTAINS
!*******************************************************************************
! SECTION III. PUBLICLY ACCESSIBLE SUBROUTINES
!*******************************************************************************

      SUBROUTINE diagnostic_dot_read (diagnostic_file,mddc_desc_a,             &
     &   n_diagn_c)
      
!  Subroutine to parse the large structure of the diagnostic_file
!  Individual coils are parsed by the diagnostic_dot_parse_xxxx subroutines

      USE safe_open_mod
      IMPLICIT NONE
      
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER*(*), INTENT(in)            :: diagnostic_file
!  2012-01-12 JDH Eliminate d_coils - information is in mddc_desc_a
!      TYPE (bsc_coilcoll), INTENT(inout)   :: d_coils
      TYPE (mddc_desc), DIMENSION(:), INTENT(inout) ::                         &
     &   mddc_desc_a
      INTEGER, INTENT(inout) :: n_diagn_c

!  diagnostic_file     Name of the diagnostic. file.
!  d_coils             A bsc_coilcoll - for storing the diagnostic coil information
!  mddc_desc_a         An array of diagnostic_desc, optional.
!  n_diagn_c           Number of diagnostic coils
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: dd_iou_0=32
      INTEGER, PARAMETER :: len_key = 30
      INTEGER, PARAMETER :: n_dd_keyword = 9

!  dd_iou_0        INTEGER - starting I/O unit number, for safe_open
!  len_key         integer - length of the dd_keyword character variables
!  n_dd_keyword    integer - number of dd_keyword s.

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(len=len_key), DIMENSION(1:n_dd_keyword) :: dd_keyword
      INTEGER, DIMENSION(1:n_dd_keyword) :: n_char_key, i_find_key

      INTEGER               :: istat
      INTEGER, DIMENSION(1) :: loc_of_max
      INTEGER        :: i, icount_coils
      CHARACTER (len=30)    :: s_name                                 
      CHARACTER (len=80)    :: l_name
      CHARACTER (len=len_key) :: this_keyword
!  Make sure that len for this_keyword is the same as len for dd_keyword.

!  dd_keyword    character array - keywords for various diagnostic coils types
!  n_char_key    integer array - lengths of the dd_keyword strings
!  i_find_key    integer array - lengths of dd_keyword strings contained in line
!  istat         integer - status variable
!  loc_of_max    Integer array of length 1, to store result from maxloc call.

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

!  Initialize the keywords and keyword lengths
      dd_keyword(1) = 'flux_loop'
      dd_keyword(2) = 'flux_loop_circular'
      dd_keyword(3) = 'magnetic_probe' 
      dd_keyword(4) = 'magnetic_probe_tokamak' 
      dd_keyword(5) = 'rogowski'
      dd_keyword(6) = 'b_rogowski'
      dd_keyword(7) = 'i_rogowski'
      dd_keyword(8) = 'f_rogowski'
      dd_keyword(9) = 'end_of_file'
      n_char_key(1:n_dd_keyword) = len_trim(dd_keyword(1:n_dd_keyword))

!  Allocate space for the xnode array
      IF (ALLOCATED(xnode)) THEN
         WRITE(*,*) 'CAUTION: xnode allocated on'
         WRITE(*,*) '   entry to diagnostic_dot_read'
         DEALLOCATE(xnode)
      END IF
      ALLOCATE(xnode(3,1:n_xnode))
!
!  Open up the 'diagnostic.' file
      dd_iou = dd_iou_0
      CALL safe_open(dd_iou, istat,                                            &
     &     TRIM(diagnostic_file), 'old', 'formatted')
      IF (istat .ne. 0) STOP 'Error opening diagnostic. file'

!  Read the first line
      READ (dd_iou, '(a)', iostat=istat) line
      IF (istat .ne. 0)                                                        &   
     &   STOP 'Error in First line of diagnostic. file '

!  Initilize the bsc_coilcoll d_coil
!  2012-01-12 JDH ELiminate d_coils. Information in mddc_desc_a
!      s_name = TRIM(line)
!      l_name = ''
!      CALL bsc_destroy(d_coils)
!      IF (PRESENT(mddc_desc_a)) THEN
!         CALL bsc_construct(d_coils,s_name,l_name,                             &
!     &      ncoil_init = SIZE(mddc_desc_a))
!      ELSE
!         CALL bsc_construct(d_coils,s_name,l_name)
!      END IF
      icount_coils = 0

!  Read the second line, and get ready to enter the parsing loop
      READ (dd_iou, '(a)', iostat=istat) line
      IF (istat .ne. 0)                                                        &   
     &   STOP 'Error in Second line of diagnostic. file '
      i_line = 2
      dd_l_done = .false.

!  Infinite Loop
!  Character variable line should be defined on entry
      DO 
!  Check for keyword.
!  Trying to identify the LONGEST keyword contained in line
         i_find_key(1:n_dd_keyword) = 0
         DO i = 1,n_dd_keyword
            istat = INDEX(line,dd_keyword(i))
            IF (istat .ne. 0) THEN
               i_find_key(i) = n_char_key(i)
            END IF
         END DO
         IF (maxval(i_find_key(1:n_dd_keyword)) .gt. 0) THEN
            loc_of_max = maxloc(i_find_key(1:n_dd_keyword))
            this_keyword =                                                     &
     &         dd_keyword(loc_of_max(1))
         ELSE
            this_keyword = ''
         END IF

!  Branch on the keyword
         SELECT CASE (this_keyword)
         
            CASE DEFAULT ! Did  not find a keyword. Read another line.
               READ (dd_iou, '(a)', iostat=istat) line
               i_line = i_line + 1
               IF (istat .ne. 0) THEN
                  dd_l_done = .true.
               END IF
            
            CASE ('flux_loop')
               CALL diagnostic_dot_parse_fl
               icount_coils = icount_coils + 1
!               CALL bsc_append(d_coils,coil_temp)
               mddc_desc_a(icount_coils) = mddc_desc_temp
               
            CASE ('flux_loop_circular')
               CALL diagnostic_dot_parse_flc
               icount_coils = icount_coils + 1
!               CALL bsc_append(d_coils,coil_temp)
               mddc_desc_a(icount_coils) = mddc_desc_temp
                           
            CASE ('magnetic_probe')
               CALL diagnostic_dot_parse_mp
               icount_coils = icount_coils + 1
!               CALL bsc_append(d_coils,coil_temp)
               mddc_desc_a(icount_coils) = mddc_desc_temp
            
            CASE ('magnetic_probe_tokamak')
               CALL diagnostic_dot_parse_mptok
               icount_coils = icount_coils + 1
!               CALL bsc_append(d_coils,coil_temp)
               mddc_desc_a(icount_coils) = mddc_desc_temp

            CASE ('rogowski','b_rogowski','i_rogowski','f_rogowski')
               !  Pass this_keyword, so that can get right units
               CALL diagnostic_dot_parse_rogow(this_keyword)
               icount_coils = icount_coils + 1
!               CALL bsc_append(d_coils,coil_temp)
               mddc_desc_a(icount_coils) = mddc_desc_temp
               
            CASE ('end_of_file')
               dd_l_done = .true.

         END SELECT
         
!  Check for too many additions to the mddc_desc_a array
         IF  (icount_coils .gt. SIZE(mddc_desc_a)) THEN 
            STOP 'diagnostic_dot: Too many coils for mddc_desc_a '
         ENDIF

         
         IF (dd_l_done) EXIT ! the infinite loop
         
      END DO

!  Clean up
      CALL bsc_destroy(coil_temp)
      CALL mddc_destroy(mddc_desc_temp)
      DEALLOCATE(xnode)

!  Close the 'diagnostic.' file
      CLOSE(dd_iou)

!  Coil Count
      n_diagn_c = icount_coils
      
      RETURN
      END SUBROUTINE diagnostic_dot_read

!*******************************************************************************
! SECTION IV. PARSING SUBROUTINES
!*******************************************************************************

!-------------------------------------------------------------------------------
!     Expected behaviour for the parsing routines:
!  On entry:
!     - character variable 'line' is defined
!     - logical variable dd_l_done is .false.
!     - integer variable i_line is set
!  During Execution:
!     - integer variable i_line is kept current
!     - reads lines from file, looks for appropriate line structure
!     - presence of inappropriate line indicates end of coil
!     - constructs the bsc_coil coil_temp
!        Use the component % raux to carry response function factors
!  On exit:
!     - 'line' contains the first inappropriate line
!     - End of File or other problem: dd_l_done set to .true.
!     
!-------------------------------------------------------------------------------

      SUBROUTINE diagnostic_dot_parse_fl
!  Subroutine to parse a flux loop
!      dd_keyword(1) = 'flux_loop'

      IMPLICIT NONE

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER            :: istat
      INTEGER     :: inode, nnode
      REAL(rprec)        :: avelsq, closelsq
      CHARACTER (len=30) :: s_name                                 
      CHARACTER (len=80) :: l_name
      
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

!  Parse beyond the keyword in the first line

!  Read the second line, identification
      READ (dd_iou, '(a)', iostat=istat) line
      i_line = i_line + 1
      IF (istat .ne. 0) THEN 
         WRITE(*,*) 'Error reading ID for flux_loop, line ', i_line
         STOP
      END IF

      s_name = TRIM(line)
      l_name = ''
      inode = 1

!  Infinite loop to read node positions
      DO

!  read a new line
         READ (dd_iou, '(a)', iostat=istat) line
         i_line = i_line + 1
         IF (istat .ne. 0) THEN
            dd_l_done = .true.
            line = ''
            EXIT ! Infinite loop
         END IF

!  Reread the line, looking for three reals. Nonzero istat indicates problem
!  Assume that problem indicates end of coil.
         READ(line, *, iostat=istat) xnode(1:3,inode)
         IF (istat .ne. 0) EXIT ! infinite loop

!  Successful read of  more xnode values. Increment inode
         inode = inode + 1

         IF (inode .gt. n_xnode) THEN
            WRITE(*,*) ' inode bigger than n_xnode.', inode, n_xnode
            STOP ' Make n_xnode bigger'
         END IF
      
      END DO

!  Prepare for creation of bsc_coil
      nnode = inode - 1

!  Check to see if coil is closed.
      IF (nnode .gt. 3) THEN
!  First, compute average length ^2
         avelsq = SUM((xnode(1,2:nnode) - xnode(1,1:nnode-1) ** 2)) +          &
     &            SUM((xnode(2,2:nnode) - xnode(2,1:nnode-1) ** 2)) +          &
     &            SUM((xnode(3,2:nnode) - xnode(3,1:nnode-1) ** 2)) 

!  Compute the length ^2 of the last-first segment
         avelsq = avelsq / nnode
         closelsq = SUM((xnode(1:3,1) - xnode(1:3,nnode)) ** 2)

         IF (closelsq /  avelsq .lt. 1.e-12) THEN
            nnode = nnode - 1
         END IF
      END IF

!  Construct the coil. Use a current of one, since it is a diagnostic coil
!  Use raux (response function factor) of 1.
!  2012-01-12 JDH - no longer use raux, use flux_factor
      CALL bsc_construct(coil_temp,'fil_loop',s_name, l_name, one, 
     &   xnode(1:3,1:nnode))

!  Construct the mddc_desc. Note - optional mrf argment is missing
!  flux_loop - units are Tesla meter^2, flux_factor is one
      CALL mddc_construct(mddc_desc_temp,s_name,l_name,'T m2',                 &
     &   one,'flux_loop',coil_temp,flux_factor = one)

      RETURN
      END SUBROUTINE diagnostic_dot_parse_fl

!-----------------------------------------------
!-----------------------------------------------

      SUBROUTINE diagnostic_dot_parse_flc
!  Subroutine to parse a circular flux loop
!      dd_keyword(2) = 'flux_loop_circular'

      IMPLICIT NONE

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER                   :: istat
      REAL(rprec)               :: this_rcirc
      REAL(rprec), DIMENSION(3) :: this_xcent, this_enhat
      CHARACTER (len=30)        :: s_name                                 
      CHARACTER (len=80)        :: l_name
      
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

!  Parse beyond the keyword in the first line

!  Read the second line, identification
      READ (dd_iou, '(a)', iostat=istat) line
      i_line = i_line + 1
      IF (istat .ne. 0) THEN  
         WRITE(*,*) 'Error reading ID, line ', i_line
         STOP 'flux_loop_circular'
      END IF

      s_name = TRIM(line)
      l_name = ''

!  Read third line, radius of the circular flux loop
      READ (dd_iou, '(a)', iostat=istat) line
      i_line = i_line + 1
      IF (istat .ne. 0) THEN
         WRITE(*,*) ' Expected radius in line ', i_line
         STOP
      END IF
         
!  Reread the line, looking for one real. 
      READ(line, *, iostat=istat) this_rcirc
      IF (istat .ne. 0) THEN
         WRITE(*,*) ' Expected radius(2) in line ', i_line
         STOP
      END IF

!  Read fourth line, center of the circular flux loop
      READ (dd_iou, '(a)', iostat=istat) line
      i_line = i_line + 1
      IF (istat .ne. 0) THEN
         WRITE(*,*) ' Expected xcent in line ', i_line
         STOP
      END IF
         
!  Reread the line, looking for three reals.
      READ(line, *, iostat=istat) this_xcent(1:3)
      IF (istat .ne. 0) THEN
         WRITE(*,*) ' Expected xcent(2) in line ', i_line
         STOP
      END IF

!  Read fifth line, vector normal to circular flux loop
      READ (dd_iou, '(a)', iostat=istat) line
      i_line = i_line + 1
      IF (istat .ne. 0) THEN
         WRITE(*,*)  ' Expected enhat in line ', i_line
         STOP
      END IF
         
!  Reread the line, looking for three reals. 
      READ(line, *, iostat=istat) this_enhat(1:3)
      IF (istat .ne. 0) THEN
         WRITE(*,*) ' Expected enhat(2) in line ', i_line
         STOP
      END IF

!  Read sixth line, next keyword
      READ (dd_iou, '(a)', iostat=istat) line
      i_line = i_line + 1
      IF (istat .ne. 0) THEN
         dd_l_done = .true.
      END IF

!  Construct the coil. Use a current of one, since it is a diagnostic coil
!  Use raux (response function factor) of 1. / (2 Pi). (D3D convention)
!  JDh 1.20.03. Change raux from 1/(2 Pi) to 1. Change in Convention.
!  To get D3D tst case to work, will need to change the signals.
!  2012-01-12 JDH - no longer use raux, use flux_factor
      CALL bsc_construct(coil_temp,'fil_circ',s_name, l_name, one,             & 
     &   rcirc = this_rcirc, xcent = this_xcent(1:3),                          &
     &   enhat = this_enhat(1:3))

!  Construct the mddc_desc. Note - optional mrf argment is missing
!  flux_loop_circular - units are Tesla meter^2, flux_factor is one
      CALL mddc_construct(mddc_desc_temp,s_name,l_name,'T m2',                 &
     &   one,'flux_loop_circular',coil_temp,flux_factor = one)

      RETURN
      END SUBROUTINE diagnostic_dot_parse_flc
      
!-----------------------------------------------
!-----------------------------------------------

      SUBROUTINE diagnostic_dot_parse_mp
!  Subroutine to parse a magnetic probe
!      dd_keyword(3) = 'magnetic_probe' 

      IMPLICIT NONE
      
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER                   :: istat
      REAL(rprec)               :: this_rcirc, flux_factor
      REAL(rprec), DIMENSION(3) :: this_xcent, this_enhat
      CHARACTER (len=30)        :: s_name                                 
      CHARACTER (len=80)        :: l_name
      
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

!  Parse beyond the keyword in the first line

!  Read the second line, identification
      READ (dd_iou, '(a)', iostat=istat) line
      i_line = i_line + 1
      IF (istat .ne. 0) THEN  
         WRITE(*,*) 'Error reading ID, line ', i_line
         STOP ' magnetic_probe'
      END IF

      s_name = TRIM(line)
      l_name = ''

!  Read third line, radius of the circular flux loop
      READ (dd_iou, '(a)', iostat=istat) line
      i_line = i_line + 1
      IF (istat .ne. 0) THEN
         WRITE(*,*) ' Expected radius in line ', i_line
         STOP
      END IF
         
!  Reread the line, looking for one real. 
      READ(line, *, iostat=istat) this_rcirc
      IF (istat .ne. 0) THEN
         WRITE(*,*) ' Expected radius(2) in line ', i_line
         STOP
      END IF

!  Read fourth line, center of the circular flux loop
      READ (dd_iou, '(a)', iostat=istat) line
      i_line = i_line + 1
      IF (istat .ne. 0) THEN
         WRITE(*,*) ' Expected xcent in line ', i_line
         STOP
      END IF
         
!  Reread the line, looking for three reals.
      READ(line, *, iostat=istat) this_xcent(1:3)
      IF (istat .ne. 0) THEN
         WRITE(*,*) ' Expected xcent(2) in line ', i_line
         STOP
      END IF

!  Read fifth line, vector normal to circular flux loop
      READ (dd_iou, '(a)', iostat=istat) line
      i_line = i_line + 1
      IF (istat .ne. 0) THEN
         WRITE(*,*) ' Expected enhat in line ', i_line
         STOP
      END IF
         
!  Reread the line, looking for three reals. 
      READ(line, *, iostat=istat) this_enhat(1:3)
      IF (istat .ne. 0) THEN
         WRITE(*,*) ' Expected enhat(2) in line ', i_line
         STOP
      END IF

!  Read sixth line, next keyword
      READ (dd_iou, '(a)', iostat=istat) line
      i_line = i_line + 1
      IF (istat .ne. 0) THEN
         dd_l_done = .true.
      END IF

!  Construct the coil. Use a current of one, since it is a diagnostic coil
!  Use raux (response function factor) of 1. / ( Pi R^2).
!  2012-01-12 JDH - no longer use raux, use flux_factor
      flux_factor = one / (pi * this_rcirc * this_rcirc) 
      CALL bsc_construct(coil_temp,'fil_circ',s_name, l_name, one,             & 
     &   rcirc = this_rcirc, xcent = this_xcent(1:3),                          &
     &   enhat = this_enhat(1:3) )

!  Construct the mddc_desc. Note - optional mrf argment is missing
!  magnetic_probe - units are Tesla, flux_factor is one over area
      CALL mddc_construct(mddc_desc_temp,s_name,l_name,'T',                    &
     &   one,'magnetic_probe',coil_temp,flux_factor = flux_factor)

      RETURN
      END SUBROUTINE diagnostic_dot_parse_mp

!-----------------------------------------------
!-----------------------------------------------

      SUBROUTINE diagnostic_dot_parse_mptok
!  Subroutine to parse a magnetic probe - tokamak parameterization
!      dd_keyword(4) = 'magnetic_probe_tokamak' 

!  Using the parameterization very similar to the one that Lang Lao
!  used in the first version of v3rfun

!  Variables in data line:
!  this_rcirc, capR, phi, z, amp
!  phi and amp are angles, in degrees.

      IMPLICIT NONE
      
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER                   :: istat
      REAL(rprec)               :: this_rcirc, flux_factor
      REAL(rprec), DIMENSION(3) :: this_xcent, this_enhat
      REAL(rprec)               :: capR, phi, z, amp 
      CHARACTER (len=30)        :: s_name                                 
      CHARACTER (len=80)        :: l_name
      
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

!  Parse beyond the keyword in the first line

!  Read the second line, identification
      READ (dd_iou, '(a)', iostat=istat) line
      i_line = i_line + 1
      IF (istat .ne. 0) THEN  
         WRITE(*,*) 'Error reading ID, line ', i_line
         STOP ' magnetic_probe'
      END IF

      s_name = TRIM(line)
      l_name = ''

!  Read third line, Reals.
      READ (dd_iou, '(a)', iostat=istat) line
      i_line = i_line + 1
      IF (istat .ne. 0) THEN
         WRITE(*,*) ' Expected MP Tokamak specification ', i_line
         STOP
      END IF
         
!  Reread the line, looking for 5 reals. 
      READ(line, *, iostat=istat) this_rcirc, capR, phi, z, amp
      IF (istat .ne. 0) THEN
         WRITE(*,*) ' Expected MP Tokamak in line ', i_line
         STOP
      END IF

!  Read sixth line, next keyword
      READ (dd_iou, '(a)', iostat=istat) line
      i_line = i_line + 1
      IF (istat .ne. 0) THEN
         dd_l_done = .true.
      END IF

!  Translate specification to bsc_coil stuff
      this_xcent(1) = capR * cos(phi * degree)
      this_xcent(2) = capR * sin(phi * degree)
      this_xcent(3) = z
      this_enhat(1) = cos(amp * degree) * cos(phi * degree)
      this_enhat(2) = cos(amp * degree) * sin(phi * degree)
      this_enhat(3) = sin(amp * degree)

!  Construct the coil. Use a current of one, since it is a diagnostic coil
!  Use raux (response function factor) of 1. / ( Pi R^2).
!  2012-01-12 JDH - no longer use raux, use flux_factor
      flux_factor = one / (pi * this_rcirc * this_rcirc) 
      CALL bsc_construct(coil_temp,'fil_circ',s_name, l_name, one,             & 
     &   rcirc = this_rcirc, xcent = this_xcent(1:3),                          &
     &   enhat = this_enhat(1:3) )

!  Construct the mddc_desc. Note - optional mrf argment is missing
!  magnetic_probe_tokamak - units are Tesla, flux_factor is one over area
      CALL mddc_construct(mddc_desc_temp,s_name,l_name,'T',                    &
     &   one,'magnetic_probe_tokamak',coil_temp,                               &
     &   flux_factor = flux_factor)

      RETURN
      END SUBROUTINE diagnostic_dot_parse_mptok

!-----------------------------------------------
!-----------------------------------------------

      SUBROUTINE diagnostic_dot_parse_rogow(this_keyword)

!  Subroutine to parse a rogowski coil
!      dd_keyword(1) = 'rogowski', 'b_rogowski', 'i_rogowski', 'f_rogowski',
      IMPLICIT NONE

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER(len=*), INTENT(in) :: this_keyword
!
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER            :: istat
      INTEGER     :: inode, nnode
      REAL(rprec)        :: avelsq, closelsq
      CHARACTER (len=30) :: s_name                                 
      CHARACTER (len=80) :: l_name
      REAL(rprec)        :: anturns, xsarea
      REAL(rprec)        :: flux_factor
      CHARACTER(len=10)  :: units
      
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

!  Parse beyond the keyword in the first line

!  Read the second line, identification
      READ (dd_iou, '(a)', iostat=istat) line
      i_line = i_line + 1
      IF (istat .ne. 0) THEN 
         WRITE(*,*) 'Error reading ID for rogowski, line ', i_line
         STOP
      END IF

      s_name = TRIM(line)
      l_name = ''

! Read Third Line, number of turns and cross-sectional area
      READ (dd_iou, '(a)', iostat=istat) line
      i_line = i_line + 1
      IF (istat .ne. 0) THEN
         WRITE(*,*) ' Expected anturns and xsarea in line ', i_line
         STOP
      END IF
         
!  Reread the line, looking for two reals. 
      READ(line, *, iostat=istat) anturns, xsarea
      IF (istat .ne. 0) THEN
         WRITE(*,*) ' Expected anturns, xsarea in line ', i_line
         STOP
      END IF

!  Check for faulty values of anturns and xsarea
      IF (anturns * xsarea .eq. zero) THEN
         WRITE(*,*) ' anturns * xsarea cant be zero', anturns, xsarea
         STOP 
      END IF

!  Infinite loop to read node positions
      inode = 1
      DO

!  read a new line
         READ (dd_iou, '(a)', iostat=istat) line
         i_line = i_line + 1
         IF (istat .ne. 0) THEN
            dd_l_done = .true.
            line = ''
            EXIT ! Infinite loop
         END IF

!  Reread the line, looking for three reals. Nonzero istat indicates problem
!  Assume that problem indicates end of coil.
         READ(line, *, iostat=istat) xnode(1:3,inode)
         IF (istat .ne. 0) EXIT ! infinite loop

!  Successful read of  more xnode values. Increment inode
         inode = inode + 1

         IF (inode .gt. n_xnode) THEN
            WRITE(*,*) ' inode bigger than n_xnode.', inode, n_xnode
            STOP ' Make n_xnode bigger'
         END IF
      
      END DO

!  Prepare for creation of bsc_coil
      nnode = inode - 1

!  Problem if only one node read in
      IF (nnode .le. 1) THEN
         WRITE(*,*) ' Rogowski: nnode .le. 1', nnode
         STOP ' In diagnostic_dot_parse_rogow'
      END IF
         
!  Construct the coil. Use a current of one, since it is a diagnostic coil
!  2012-01-12 - No longer use raux, taken care of with flux_factor

      CALL bsc_construct(coil_temp,'fil_rogo',s_name, l_name, one,             &
     &   xnode(1:3,1:nnode),anturns = anturns,xsarea = xsarea)
     

      SELECT CASE(this_keyword(1:2))

         CASE('b_','ro')  ! Signal is Integral(B dl) / Integral(dl)
            flux_factor = one / (anturns * xsarea)
            units = 'T'
         
         CASE('i_')   !  Signal is Integral(B dl) / mu0
            flux_factor = one/(coil_temp % ave_n_area * 2.d-7 * twopi)
            units = 'A'
         
         CASE('f_')   !  Signal is N A Integral(B dl) / Integral(dl)
                      !  N - total number of turns (anturns)
                      !  A - cross sectional area (xsarea)
            flux_factor = one
            units = 'T m2'

         CASE DEFAULT
            WRITE(*,*) 'CASE DEFAULT ERROR: diagnostic_dot_parse_rogow'
         
         END SELECT
            
!  Construct the mddc_desc. Note - optional mrf argment is missing
      CALL mddc_construct(mddc_desc_temp,s_name,l_name,TRIM(units),            &
     &   one,TRIM(this_keyword),coil_temp,flux_factor = flux_factor)

      RETURN
      END SUBROUTINE diagnostic_dot_parse_rogow

!*******************************************************************************
! SECTION V. AUXILIARY SUBROUTINES
!*******************************************************************************

      END MODULE diagnostic_dot
