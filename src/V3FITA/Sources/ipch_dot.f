!*******************************************************************************
!  File ipch_dot.f
!  Contains module ipch_dot
!    Module for opening and reading a 'ipch.' file, and then placing the
!     interferometry-polarimetry chord information into an array of ipch_desc
!!  
!-------------------------------------------------------------------------------
!   DEPENDENCIES
!-------------------------------------------------------------------------------
!
!    This module uses the following modules:
!       stel_kinds
!       stel_constants
!       ipch_T
!       safe_open_mod
!
!-------------------------------------------------------------------------------
!   CHANGE HISTORY
!-------------------------------------------------------------------------------
!
!  Based on module sxrch_dot
!  First version JDH 2012-03-15
!
!-------------------------------------------------------------------------------
!   USAGE
!-------------------------------------------------------------------------------
!     Usage: Here is a sample of code that shows how to use this module
!
!      USE ipch_T
!      USE ipch_dot
!      TYPE(ipch_desc), DIMENSION(1000) :: ipch_descriptions_a
! ...
! !  Read in interferometry-polarimetry chords
!      CALL ipch_dot_read(file_name_ipch_dot,ipch_descriptions_a)
!
!-------------------------------------------------------------------------------
!   FILE STRUCTURE - ipch. file
!-------------------------------------------------------------------------------
!
!  ipch. File Structure:
! line 1: Character - identifier
!
!      (0 or more lines, no keyword present. All ignored by parser.)
!    keyword
!      lines containing data - number depends on keyword
!
!      (0 or more lines, no keyword present. All ignored by parser.)
!    keyword
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
!       ipchd_keyword(1) = 'ip_chord_XYZ'
!       ipchd_keyword(2) = 'ip_chord_RPhiDegZ'
!       ipchd_keyword(3) = 'end_of_file' 
!
!*******************************************************************************

!*******************************************************************************
!  MODULE ipch_dot
!    
! SECTION I. VARIABLE DECLARATIONS
! SECTION II. INTERFACE BLOCKS
! SECTION III. PUBLICLY ACCESSIBLE SUBROUTINES
! SECTION IV. PARSING SUBROUTINES
! SECTION V. AUXILIARY SUBROUTINES
!*******************************************************************************
      MODULE ipch_dot

!*******************************************************************************
! SECTION I. VARIABLE DECLARATIONS
!*******************************************************************************

      USE stel_kinds
      USE stel_constants
      USE ipch_T
      IMPLICIT NONE

!-------------------------------------------------------------------------------
!  Variables local to the module. Make them private.
!-------------------------------------------------------------------------------
      INTEGER ::  ipchd_iou
      INTEGER :: i_line
      CHARACTER(len=200) :: line
      LOGICAL :: ipchd_l_done
      TYPE (ipch_desc), SAVE :: ipch_desc_temp

      PRIVATE ipchd_iou,i_line 
      PRIVATE line, ipchd_l_done, ipch_desc_temp

!  ipch_desc_temp    A temporary ipch_desc. SAVE attribute - see note in 
!                         diagostic_dot about Lahey Fujitsu compiler
!  ipchd_iou          i/o unit number for ipch_dot file
!  ipchd_l_done       logical - true if done reading file
!  i_line             integer - line counter for file
!  line               character - the last line read from the file
          
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

      SUBROUTINE ipch_dot_read (ipch_file,ipch_desc_arr,n_chords)
      
!  Subroutine to parse the large structure of the ipch_dot_file
!  Individual ip chords are parsed by the ipch_dot_parse_xxxx subroutines

      USE safe_open_mod
      IMPLICIT NONE
      
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER*(*), INTENT(in)            :: ipch_file
      TYPE (ipch_desc), DIMENSION(:), INTENT(inout) :: ipch_desc_arr
      INTEGER, INTENT(inout)  :: n_chords


!  ipch_file            Name of the ipch. file.
!  ipch_desc_arr        Array of ipch_desc. Assumed to have been allocated
!                        before call to ipch_dot_read
!  n_chords              Number of chords read from the file

!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: ipchd_iou_0=32
      INTEGER, PARAMETER :: len_key = 30
      INTEGER, PARAMETER :: n_ipchd_keyword = 3

!  ipchd_iou_0            INTEGER - starting I/O unit number, for safe_open
!  len_key                 integer - length of the ipchd_keyword character variables
!  n_ipchd_keyword        integer - number of ipchd_keyword s.

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(len=len_key), DIMENSION(1:n_ipchd_keyword) ::                 &
     &   ipchd_keyword
      INTEGER, DIMENSION(1:n_ipchd_keyword) :: n_char_key, i_find_key

      INTEGER               :: istat
      INTEGER, DIMENSION(1) :: loc_of_max
      INTEGER               :: i, icount_chords
      INTEGER               :: n_ipch_desc_arr
      CHARACTER (len=30)    :: s_name                                 
      CHARACTER (len=80)    :: l_name
      CHARACTER (len=len_key) :: this_keyword
!  Make sure that len for this_keyword is the same as len for ipchd_keyword.

!  ipchd_keyword    character array - keywords for various ipch input types
!  n_char_key        integer array - lengths of the ipchd_keyword strings
!  i_find_key        integer array - lengths of ipchd_keyword strings contained in line
!  istat             integer - status variable
!  loc_of_max        Integer array of length 1, to store result from maxloc call.

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

!  Initialize the keywords and keyword lengths
      ipchd_keyword(1) = 'ip_chord_XYZ'
      ipchd_keyword(2) = 'ip_chord_RPhiDegZ'
      ipchd_keyword(3) = 'end_of_file'
      n_char_key(1:n_ipchd_keyword) =                                         &
     &   len_trim(ipchd_keyword(1:n_ipchd_keyword))

!  Find the length of the ipch_desc_arr array
      n_ipch_desc_arr = SIZE(ipch_desc_arr)

!  Open up the 'ipch.' file
      ipchd_iou = ipchd_iou_0
      CALL safe_open(ipchd_iou, istat,                                        &
     &     TRIM(ipch_file), 'old', 'formatted')
      IF (istat .ne. 0) STOP 'Error opening ipch. file'

!  Read the first line
      READ (ipchd_iou, '(a)', iostat=istat) line
      IF (istat .ne. 0)                                                        &   
     &   STOP 'Error in First line of ipch. file '


!  Read the second line, and get ready to enter the parsing loop
      READ (ipchd_iou, '(a)', iostat=istat) line
      IF (istat .ne. 0)                                                        &   
     &   STOP 'Error in Second line of ipch. file '
      i_line = 2
      ipchd_l_done = .false.

      icount_chords = 0

!  Infinite Loop
!  Character variable line should be defined on entry
      DO 
!  Check for keyword.
!  Trying to identify the LONGEST keyword contained in line
         i_find_key(1:n_ipchd_keyword) = 0
         DO i = 1,n_ipchd_keyword
            istat = INDEX(line,ipchd_keyword(i))
            IF (istat .ne. 0) THEN
               i_find_key(i) = n_char_key(i)
            END IF
         END DO
         IF (maxval(i_find_key(1:n_ipchd_keyword)) .gt. 0) THEN
            loc_of_max = maxloc(i_find_key(1:n_ipchd_keyword))
            this_keyword =                                                     &
     &         ipchd_keyword(loc_of_max(1))
         ELSE
            this_keyword = ''
         END IF

!  Branch on the keyword
         SELECT CASE (this_keyword)
         
            CASE DEFAULT ! Did  not find a keyword. Read another line.
               READ (ipchd_iou, '(a)', iostat=istat) line
               i_line = i_line + 1
               IF (istat .ne. 0) THEN
                  ipchd_l_done = .true.
               END IF
            
            CASE ('ip_chord_XYZ')
               CALL ipch_dot_parse_chord('XYZ')
               icount_chords = icount_chords + 1
               IF (icount_chords .gt. n_ipch_desc_arr) THEN
                  WRITE(*,*) 'TOO MANY ipCH. Increase dimensions'
               ELSE
                  ipch_desc_arr(icount_chords) = ipch_desc_temp
               ENDIF
               
            CASE ('ip_chord_RPhiDegZ')
               CALL ipch_dot_parse_chord('RPHiDegZ')
               icount_chords = icount_chords + 1
               IF (icount_chords .gt. n_ipch_desc_arr) THEN
                  WRITE(*,*) 'TOO MANY ipCH. Increase dimensions'
               ELSE
                  ipch_desc_arr(icount_chords) = ipch_desc_temp
               ENDIF
                           
            CASE ('end_of_file')
               ipchd_l_done = .true.

         END SELECT
         
!  Check for too many additions to the ipch_desc_arr array
         IF  (icount_chords .gt. SIZE(ipch_desc_arr)) THEN
            STOP 'ipch_dot: Too many coils for ipch_desc_arr '
         ENDIF
        
         IF (ipchd_l_done) EXIT ! the infinite loop
         
      END DO

!  Clean up
      CALL ipch_desc_destroy(ipch_desc_temp)

!  Close the 'ipch.' file
      CLOSE(ipchd_iou)
      
      n_chords = icount_chords
      
      RETURN
      END SUBROUTINE ipch_dot_read

!*******************************************************************************
! SECTION IV. PARSING SUBROUTINES
!*******************************************************************************

!-------------------------------------------------------------------------------
!     Expected behaviour for the parsing routines:
!  On entry:
!     - character variable 'line' is defined
!     - logical variable ipchd_l_done is .false.
!     - integer variable i_line is set
!  During Execution:
!     - integer variable i_line is kept current
!     - reads lines from file, looks for appropriate line structure
!     - presence of inappropriate line indicates end of chord
!     - constructs ipch_desc_temp
!  On exit:
!     - 'line' contains the first inappropriate line
!     - End of File or other problem: ipchd_l_done set to .true.
!     
!-------------------------------------------------------------------------------

      SUBROUTINE ipch_dot_parse_chord(coordinate_type)
!  Subroutine to parse an ipch chord (redundant chord there)
!      ipchd_keyword(1) = 'ip_chord_XYZ' or 'ip_chord_RPhiDegZ'

      IMPLICIT NONE

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER (len=*) :: coordinate_type
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER                     :: istat
      INTEGER                     :: inode, nnode
      REAL(rprec)                 :: phi_rad
      REAL(rprec), DIMENSION(3)   :: x_start, x_end, xcart_i, xcart_f
      CHARACTER (len=30)          :: chord_name                               
      
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

!  Parse beyond the keyword in the first line

!  Read the second line, chord identification
!-----------------------------------------------
      READ (ipchd_iou, '(a)', iostat=istat) line
      i_line = i_line + 1
      IF (istat .ne. 0) THEN 
         WRITE(*,*) 'Error reading ID for ip chord, line ', i_line
         STOP
      END IF

      chord_name = TRIM(line)

!  Read the third line, start position
!-----------------------------------------------
      READ (ipchd_iou, '(a)', iostat=istat) line
      i_line = i_line + 1
      IF (istat .ne. 0) THEN 
         WRITE(*,*) 'Error ip start position line', i_line
         STOP
      END IF

!  Reread the line, looking for three reals
         READ(line, *, iostat=istat) x_start(1:3)
      IF (istat .ne. 0) THEN 
         WRITE(*,*) 'Error ip start position', i_line
         STOP
      END IF

!  Read the fourth line, end position
!-----------------------------------------------
      READ (ipchd_iou, '(a)', iostat=istat) line
      i_line = i_line + 1
      IF (istat .ne. 0) THEN 
         WRITE(*,*) 'Error ip end position line', i_line
         STOP
      END IF

!  Reread the line, looking for three reals
         READ(line, *, iostat=istat) x_end(1:3)
      IF (istat .ne. 0) THEN 
         WRITE(*,*) 'Error ip end position', i_line
         STOP
      END IF

!  Done reading information for this ipch. Generate the ipch
!-----------------------------------------------

      IF (TRIM(coordinate_type) .eq. 'RPHiDegZ') THEN
         phi_rad = x_start(2) * pi / 180.
         xcart_i(1) = x_start(1) * COS(phi_rad)
         xcart_i(2) = x_start(1) * SIN(phi_rad)
         xcart_i(3) = x_start(3)
         phi_rad = x_end(2) * pi / 180.
         xcart_f(1) = x_end(1) * COS(phi_rad)
         xcart_f(2) = x_end(1) * SIN(phi_rad)
         xcart_f(3) = x_end(3)
      ELSE
         xcart_i = x_start
         xcart_f = x_end
      ENDIF
      
      CALL ipch_desc_construct(ipch_desc_temp,chord_name,                    &
     &    xcart_i,xcart_f)

!  Read the next line, in preparation for returning to main parsing loop
!-----------------------------------------------
      READ (ipchd_iou, '(a)', iostat=istat) line
      i_line = i_line + 1
      IF (istat .ne. 0) THEN 
         ipchd_l_done = .true.
      END IF

      RETURN
      END SUBROUTINE ipch_dot_parse_chord


!*******************************************************************************
! SECTION V. AUXILIARY SUBROUTINES
!*******************************************************************************

      END MODULE ipch_dot
