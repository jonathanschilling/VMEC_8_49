!*******************************************************************************
!  File eq_interface.f
!  Contains module eq_interface

!  Module is part of the V3FIT stellarator equilibrium reconstruction code.
!  It deals with the interface with the equilibrium code.
!  Right now (6/2004), the VMEC code is the only equilibrium code.

!  NB. 2011-07-28
!    This file now contains some functions OUTSIDE the eq_interface module

!*******************************************************************************
!  MODULE eq_mod
!    (Equilibrium)
! SECTION I.      VARIABLE DECLARATIONS
! SECTION II.     INTERFACE BLOCKS
! SECTION III.    EQUILIBRIUM INITIALIZATION SUBROUTINES
! SECTION IV.     EQUILIBRIUM STEP ROUTINES
! SECTION V.      EQUILIBRIUM CHANGE SUBROUTINES
! SECTION VI.     EQUILIBRIUM GET INFORMATION ROUTINES
! SECTION VII.    EQUILIBRIUM AUXILIARY CALCULATION ROUTINES
! SECTION VIII.   EQUILIBRIUM WRITE ROUTINES
! SECTION IX.     FUNCTIONS & SUBROUTINES OUTSIDE THE MODULE
! SECTION X.      COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
      MODULE eq_interface

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
!  Equilibrium Derived Types
!-------------------------------------------------------------------------------
      USE eq_T     
!-------------------------------------------------------------------------------
!  V3 Utilities
!-------------------------------------------------------------------------------
      USE v3_utilities
      
!-------------------------------------------------------------------------------
!  V3 Global
!-------------------------------------------------------------------------------
      USE v3f_global, ONLY: l_zero_xcdot
!  JDH 2008-08-04
!     The USE v3f_global is a KLUDGE. This is an UGLY way to get control 
!     information to the eq_interface. Use for testing, and try to GET RID OF IT

      IMPLICIT NONE
      
!-------------------------------------------------------------------------------
!  Variables for SAVING and RESTORING the VMEC 2D preconditioning variable values
!-------------------------------------------------------------------------------
      CHARACTER(len=10)  :: precon_type_save   ! module vmec_input
      REAL(rprec) :: prec2d_threshold_save     ! module vmec_input
      INTEGER :: ictrl_prec2d_save             ! module precon2d
      LOGICAL :: l_comp_prec2D_save            ! module precon2d
      LOGICAL :: lqmr_save                     ! module gmres_mod (added 2011-10-17)
      
!*******************************************************************************
! SECTION II.  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!SPH  Interface to external runvmec routine
!-------------------------------------------------------------------------------
      INTERFACE
         SUBROUTINE runvmec(ictrl_array, input_file, lscreen,                  &
     &                      reset_file_name)
         INTEGER, INTENT(inout), TARGET :: ictrl_array(5)
         LOGICAL, INTENT(in)            :: lscreen
         CHARACTER(LEN=*), INTENT(in)   :: input_file
         CHARACTER(LEN=*), OPTIONAL     :: reset_file_name
         END SUBROUTINE runvmec

      END INTERFACE
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
      CONTAINS
          
!*******************************************************************************
! SECTION III.  EQUILIBRIUM INITIALIZATION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Initialize the equilibrium solver - from a file
!
!    Reads information in from the file 
!    Initializes the equilibrium solver
!    Constructs the equilibrium structures
!
!    For now, specific to the VMEC code.
!-------------------------------------------------------------------------------
      SUBROUTINE eq_init_file(filename,state)

!-------------------------------------------------------------------------------
!  USE statements
!-------------------------------------------------------------------------------

! VMEC input variables
      USE vmec_input, ONLY: lasym, lfreeb, lforbal, lrfp,                      &
     &   mgrid_file, nfp,                                                      &
     &   ncurr, nvacskip, mpol, ntor, ntheta, nzeta, ns_array,                 &
     &   pcurr_type, piota_type, pmass_type,                                   &
     &   niter, nstep, delt, ftol_array, gamma, tcon0,                         &
     &   ac, ac_aux_s, ac_aux_f, ai, ai_aux_s, ai_aux_f,                       &
     &   am, am_aux_s, am_aux_f, bloat, rbc, zbs, zbc, rbs,                    &     
     &   extcur, curtor, phiedge, pres_scale, l_v3fit

! VMEC dimensions
!  10-21-04 Check with Steve Hirshman, are there any gotchas with this
      USE vmec_dim, ONLY: mns

! VMEC flags for communicating with runvmec
      USE vmec_params, ONLY: restart_flag, readin_flag, timestep_flag,         &
     &   more_iter_flag, norm_term_flag, version => version_

! VMEC fsq - for monitoring convergence
      USE vmec_main, ONLY: fsqr, fsqz, fsql

! VMEC xstuff - state arrays
      USE xstuff, ONLY: xc_vmec_1 => xc, xcdot_vmec_1 => xcdot

! VMEC nextcur - number of external currents
      USE mgrid_mod, ONLY: nextcur

!  VMEC_HISTORY
      USE vmec_history

!-------------------------------------------------------------------------------
!  Declarations
!-------------------------------------------------------------------------------
!  Declare Arguments 
      CHARACTER(len=*), INTENT(inout) :: filename
      TYPE (eq_state), INTENT (inout) :: state

!  Declare local variables
      TYPE (eq_param_fix) :: fixp
      TYPE (eq_param_var) :: varp
      CHARACTER(len=*), PARAMETER :: code = "VMEC2000"
      CHARACTER(len=30) :: s_id = " "
      CHARACTER(len=80) :: l_id = " "
      INTEGER :: ns_index
      REAL(rprec) :: fsq_max

      INTEGER :: ictrl_array(5), ier_flag
      INTEGER :: rbc_bounds(4) ! order lower_1, upper_1, lower_2, upper_2
      LOGICAL, PARAMETER :: lscreen     = .false.
      INTEGER :: index_dat, index_end
      CHARACTER(len=100) :: input_extension
      CHARACTER(len=100) ::wout_filename
      
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'eq_init_file: '

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------

!  Initialize VMEC: read in the file filename to load module vmec_input variables
!  AND initialize xc_vmec state array. Flags = reset(1) + initialize(2) + 
!  timestep(4) = 7.

      ictrl_array = 0
      ictrl_array(1) = restart_flag + readin_flag + timestep_flag
      ictrl_array(3) = 1           ! Run one step to initialize xc and xcdot

!  JDH 11-30-04. For now, limit VMEC to use only the first value in the 
!  ns_array. (i.e., no multigrid, just a single grid: ns_index = 1)
      ns_index = 1
      ictrl_array(4) = ns_index        

!  Set the l_v3fit logical
      l_v3fit = .true.

!  Set vmec_history printing off, so that VMEC doesn't control
!  Set vmec history integers - indicates where called from
      CALL vmec_history_print_flag_off
      CALL eq_history_set(i1 = 0,i2 = 0) 

      CALL runvmec(ictrl_array, filename, lscreen)
      ier_flag = ictrl_array(2)
         
!  Check for abnormal terminations from runvmec
      IF ((ier_flag .ne. norm_term_flag) .and.                                 &
     &    (ier_flag .ne. more_iter_flag)) THEN
             CALL err_fatal(sub_name // 'Bad ier_flag', int = ier_flag)
      ENDIF

!      CALL assert_eq(0,ier_flag,sub_name // 'Reading vmec input data')
!      IF (ier_flag .ne. 0) STOP 'Error reading in vmec input data!'

      fsq_max= MAX(fsqr,fsqz,fsql)
      WRITE(*,*) 'In eq_init_file, fsq_max = ', fsq_max

!  Construct the structures
      CALL eq_param_fix_construct(fixp,                                        &
     &   lasym, lfreeb, lforbal, lrfp, mgrid_file, nfp, ncurr,                 &
     &   nvacskip, mpol, ntor, ntheta, nzeta, mns, ns_array, nextcur,          &
     &   pcurr_type, piota_type, pmass_type,                                   &
     &   niter, nstep, delt, ftol_array, gamma, tcon0)
     
!  2011-01-03 Need to communicate bounds for rbc, etc  correctly
!  Assume rbc, zbc, rbs, zbs are all dimensioned the same.
!  order assumed is lower_1, upper_1, lower_2, upper_2
      rbc_bounds(1) = LBOUND(rbc,1)
      rbc_bounds(2) = UBOUND(rbc,1)
      rbc_bounds(3) = LBOUND(rbc,2)
      rbc_bounds(4) = UBOUND(rbc,2)
      CALL eq_param_var_construct(varp,ns_index, rbc_bounds,                   &
     &   ac, ac_aux_s, ac_aux_f, ai, ai_aux_s, ai_aux_f,                       &
     &   am, am_aux_s, am_aux_f, bloat, rbc, zbs, zbc, rbs,                    &     
     &   extcur, curtor, phiedge, pres_scale)
     
!  JDH 2011-0228. Change below, so that VMEC filenames have input extension
!  Coding for input_extension based on code in runvmec
!      CALL eq_state_construct(state,code,version,s_id,l_id,                    &
!     &   fixp,varp,xc_vmec_1,xcdot_vmec_1,' ',fsq_max)
      index_dat = INDEX(filename,'input.')
      index_end = LEN_TRIM(filename)
      IF (index_dat .gt. 0) THEN
         input_extension  = filename(index_dat+6:index_end)
      ELSE
         input_extension = filename(1:index_end)
      END IF
      wout_filename = 'wout_' // TRIM(input_extension) // '.nc'
      
      CALL eq_state_construct(state,code,version,s_id,l_id,                    &
     &   fixp,varp,xc_vmec_1,xcdot_vmec_1,wout_filename,fsq_max)

      RETURN
      END SUBROUTINE eq_init_file

!-------------------------------------------------------------------------------
!  Initialize the equilibrium solver - from equilibrium structures
!
!    Reads information in from the structures
!    Initializes the equilibrium solver
!
!    For now, specific to the VMEC code.
!-------------------------------------------------------------------------------
      SUBROUTINE eq_init_structure(state)

! Note to SH. I think you'll want some USE statements here, to get access to the
! VMEC variables

!  Declare Arguments 
      TYPE (eq_state), INTENT (inout) :: state

!  Declare local variables

!  Start of executable code

!  Initialize VMEC, using the information in fixp, varp, and state

      RETURN
      END SUBROUTINE eq_init_structure

!*******************************************************************************
! SECTION IV.  EQUILIBRIUM STEP SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Step the equilibrium solver
!    Takes as input the number of iterations numiter and the eq_state state
!    The state that is iterated is the  * internal * state of the equilibrium solver
!    The eq_state argument state contains the state on output.
!    The equilibrium solver MUST have been initialized already.
!
!    For now, specific to the VMEC code.
!-------------------------------------------------------------------------------

      SUBROUTINE eq_step(numiter,state,lconverged,iter_es,                     &
     &   ier_flag_vmec,vmec_flag_name,delta_iterc)

!  JDH 2009-09-10. Changed err_fatal to err_warn for various VMEC flags
!  JDH 2009-07-16. Revised and cleaned up coding

!-------------------------------------------------------------------------------
!  USE statements
!-------------------------------------------------------------------------------

! VMEC xstuff - state arrays
      USE xstuff, ONLY: xc_vmec_2 => xc, xcdot_vmec_2 => xcdot

! VMEC flags for communicating with runvmec
      USE vmec_params

! VMEC fsq - for monitoring convergence
!      iterc - vmec cumulative internal iteration counter
!      SPH013108: add iter1 here
      USE vmec_main, ONLY: fsqr, fsqz, fsql, iter1, iter2, iterc

!  VMEC history - iteration counter
!      USE vmec_history, ONLY: iter_ha
!     ^  JDH 2010-08-03. Use iter2 from vmec_main for iter_es argument

!  IOU for runlog file
      USE v3f_global, ONLY: iou_runlog

!  VMEC inout state. Error recovery is performed directly with out altering the 
!  orginal state.
      USE vmec_input, ONLY: delt

!-------------------------------------------------------------------------------
!  Declarations
!-------------------------------------------------------------------------------
!  numiter         number of iterations for vmec to perform
!  state           Type eq_state, the equilibrium state to iterate
!  lconverged      Logical, whether or not equilibrium has converged
!  iter_es         integer, equilibrium solver internal iteration counter.
!                  (07-24-06 - use iter_ha from vmec_history)
!                  (2010-08-03 - use iter2 from vmec_main)
!                  (2012-06-20 - use iterc from vmec_main)
!  ier_flag_vmec   integer ier_flag returned from runvmec
!  vmec_flag_name  name of the ier_flag from VMEC
!  delta_iterc     Change in VMEC iterc counter, due to runvmec call
!  recovIndex      Integer counter, to keep track of number of recovery attempts
!  delt_local_store  Value of delt on entry. Store, for restoration after
!                    error recovery
!  
!  Declare Arguments 
      INTEGER, INTENT(inout) :: numiter
      TYPE (eq_state), INTENT (inout) :: state
      LOGICAL, INTENT(out), OPTIONAL :: lconverged
      INTEGER, INTENT(out), OPTIONAL :: iter_es
      INTEGER, INTENT(out), OPTIONAL :: ier_flag_vmec
      CHARACTER	(len=*), INTENT(inout), OPTIONAL :: vmec_flag_name
      INTEGER, INTENT(out), OPTIONAL :: delta_iterc     

!  Declare local variables
      INTEGER :: ictrl_array(5), ier_flag, ier1, iterc_before
      LOGICAL, PARAMETER :: lscreen     = .false.
      LOGICAL :: lconverged_local
      CHARACTER	(len=80) :: vmec_flag_name_local
      INTEGER :: index_end
      CHARACTER	(len=100) :: input_filename
      INTEGER :: recovIndex = 0
      REAL(rprec) :: delt_local_store
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'eq_step: '

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!  Check consistency of state variables and VMI variables
!-------------------------------------------------------------------------------
      CALL assert(ASSOCIATED(state % xc),ASSOCIATED(state % xcdot),            &
     &   sub_name // 'state variables not associated')
      CALL assert_eq(SIZE(state % xc),SIZE(xc_vmec_2),                         &
     &   sub_name // 'xc and state % xc not the same size')
      CALL assert_eq(SIZE(state % xcdot),SIZE(xcdot_vmec_2),                   &
     &   sub_name // 'xcdot and state % xcdot not the same size')
      CALL assert(state % fixp % delt .eq. delt,                               &
     &   sub_name // 'state % fixp % delt .ne. vmec_input:delt')

!-------------------------------------------------------------------------------
!  Error recovery - save initial values in _local_store variables
!-------------------------------------------------------------------------------

      delt_local_store = delt

!-------------------------------------------------------------------------------
!  Run vmec until successful or an unrecoverable error has occured.
!-------------------------------------------------------------------------------

!  Internal VMEC state is unchanged before call to runvmec.
      DO

!-------------------------------------------------------------------------------
!  Setup call to runvmec
!-------------------------------------------------------------------------------
!  Use ns_index from the state
         ictrl_array = 0
         ictrl_array(1) = timestep_flag + output_flag + reset_jacdt_flag
         ictrl_array(3) = numiter
         ictrl_array(4) = state % varp % ns_index
         iterc_before = iterc
         iter1        = iter2
      
!  Construct input_filename from state % wout_filename
!  Assumes wout_filename is  = 'wout_' // TRIM(input_extension) // '.nc'

         index_end = LEN_TRIM(state % wout_filename)
         input_filename =                                                      &
     &      'input.'//state % wout_filename(6:index_end - 3)
      
         IF (l_zero_xcdot) CALL eq_change_vmi_zero_xcdot
!  l_zero_xcdot gets here via a USE v3f_global
!  For testing 2008-08-04
         CALL runvmec(ictrl_array, input_filename, lscreen)
         CALL eq_aux_12_undefine(state)

!-------------------------------------------------------------------------------
!  Check error flags for recovery
!-------------------------------------------------------------------------------
!  Set error flag. (Flags declared in vmec_params.f)
!
!  Check error flags for error recovery modes.
!  If a recoverable error has occured, rerun vmec with altered eq_state. Error
!  recovery is attempted a maximum of 2 times.
!  If a non-recoverable error has occured, do not rerun vmec and report error.
!  If vmec terminated sucessfully continue.
!  
!  When altering a fixed parameter, change the value directly in the 
!  VMI (VMEC Internal State) and in the eq_state, to keep them consistent.
!  (Value should have been stored in a _local_store variable, for later
!  restoration.)

         ier_flag = ictrl_array(2)
         SELECT CASE (ier_flag)
         CASE (norm_term_flag)
            vmec_flag_name_local = 'norm_term_flag'       ! = 0
            EXIT
         CASE (bad_jacobian_flag)
            vmec_flag_name_local = 'bad_jacobian_flag'    ! = 1
            EXIT
         CASE (more_iter_flag)
            vmec_flag_name_local = 'more_iter_flag'       ! = 2
            EXIT
         CASE (jac75_flag) ! Jacobian failed decrease DELT
            IF (recovIndex .le. 2) THEN
               WRITE(*,*) sub_name, ' delt change: old, new:',                 &
     &            delt, delt / 2.
               WRITE(iou_runlog,*) sub_name, ' delt change: old, new:',        &
     &            delt, delt / 2.
               delt = delt / 2.0                          ! Note - no EXIT
               state % fixp % delt = delt                 ! Keep state consistent
                                                          ! with VMI
            ELSE
               vmec_flag_name_local = 'jac75_flag'        ! = 4
               EXIT
            ENDIF
         CASE (input_error_flag)
            vmec_flag_name_local = 'input_error_flag'     ! = 5
            EXIT
         CASE (phiedge_error_flag)
            vmec_flag_name_local = 'phiedge_error_flag'   ! = 7
            EXIT
         CASE (ns_error_flag)
            vmec_flag_name_local = 'ns_error_flag'        ! = 8
            EXIT
         CASE (misc_error_flag)
            vmec_flag_name_local = 'misc_error_flag'      ! = 9
            EXIT
         CASE (successful_term_flag)
            vmec_flag_name_local = 'successful_term_flag' ! = 11
            EXIT
         CASE DEFAULT
            vmec_flag_name_local = 'default'
            EXIT
         END SELECT

!-------------------------------------------------------------------------------
!  Reset vmec state with altered recovery variable.
!-------------------------------------------------------------------------------
!  'state' remains unaltered until VMEC variables are copied back into 'state' 
!  after errors are reported.
         WRITE (*,*) 'Attempting error recovery'
         WRITE (iou_runlog,*) 'Attempting error recovery'
!  Reset vmec xc variables.
         CALL eq_change_vmi_cp_xc(state)
         recovIndex = recovIndex + 1
      END DO

!-------------------------------------------------------------------------------
!  Report errors
!-------------------------------------------------------------------------------

      IF (PRESENT(vmec_flag_name)) vmec_flag_name =                            &
     &      vmec_flag_name_local

      IF (PRESENT(ier_flag_vmec)) ier_flag_vmec = ier_flag
         
!  Check for abnormal terminations from runvmec. Only warn
      IF ((ier_flag .ne. norm_term_flag) .and.                                 &
     &       (ier_flag .ne. successful_term_flag) .and.                        &
     &       (ier_flag .ne. more_iter_flag)) THEN
          CALL err_warn(sub_name // 'Bad ier_flag ' //                         &
     &            vmec_flag_name_local, int = ier_flag)
      ENDIF
      
!-------------------------------------------------------------------------------
!  Copy the VMEC variables into the state variables
!-------------------------------------------------------------------------------
!  First, need to check to see if sizes have changed
!  (They would have changed if ns_index was increased)
      IF (SIZE(state % xc) .ne. SIZE(xc_vmec_2)) THEN
         DEALLOCATE(state % xc,STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'dealloc xc')
         ALLOCATE(state % xc(1:SIZE(xc_vmec_2)),STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'alloc xc')
      ENDIF
      IF (SIZE(state % xcdot) .ne. SIZE(xcdot_vmec_2)) THEN
         DEALLOCATE(state % xcdot,STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'dealloc xcdot')
         ALLOCATE(state % xcdot(1:SIZE(xcdot_vmec_2)),STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'alloc xcdot')
      ENDIF

      state % xc = xc_vmec_2
      state % xcdot = xcdot_vmec_2
      state % fsq_max = MAX(fsqr,fsqz,fsql)
!  JDH 2011-0228. Comment out line below. Correct wout_filename should
!    have been put into the eq_state in eq_init_file.
!  Need to get correct wout filename.
!      state % wout_filename = 'wout_.nc'

!-------------------------------------------------------------------------------
!  Error recovery - restore initial values from _local_store variables
!-------------------------------------------------------------------------------
      delt = delt_local_store
      state % fixp % delt = delt_local_store

!-------------------------------------------------------------------------------
!  Logic to see if have converged
!-------------------------------------------------------------------------------
      lconverged_local = .false.
      IF (state%fsq_max .le.                                                   &
     &   state%fixp%ftol_array(state%varp%ns_index)) THEN
         lconverged_local = .true.
      ENDIF
      IF (PRESENT(lconverged)) lconverged = lconverged_local
      
!-------------------------------------------------------------------------------
!  Diagnostic Printout
!-------------------------------------------------------------------------------
      WRITE(*,1000) iterc, iterc - iterc_before,                               &
     &   vmec_flag_name_local, ier_flag, lconverged_local
      WRITE(iou_runlog,1000) iterc, iterc - iterc_before,                      &
     &   vmec_flag_name_local, ier_flag, lconverged_local
1000  FORMAT('   eq_step:', i7,1x,i7,1x,a20,1x,i3,1x,l1)

!  Print out more information when have not converged
      IF (.not. lconverged_local) THEN
         WRITE(*,*) ' lconverged_local = ', lconverged_local,                  &
     &      ' in eq_step'
         WRITE(*,*) ' state%fsq_max = ',state%fsq_max
         WRITE(*,*) ' state%fixp%ftol_array(state%varp%ns_index) = ',          &
     &      state%fixp%ftol_array(state%varp%ns_index)
         WRITE(*,*) ' vmec_flag_name_local = ', vmec_flag_name_local
         WRITE(*,*) ' ier_flag = ', ier_flag
         WRITE(iou_runlog,*) ' lconverged_local = ', lconverged_local,         &
     &      ' in eq_step'
         WRITE(iou_runlog,*) ' state%fsq_max = ',state%fsq_max
         WRITE(iou_runlog,*)                                                   &
     &      ' state%fixp%ftol_array(state%varp%ns_index) = ',                  &
     &      state%fixp%ftol_array(state%varp%ns_index)
         WRITE(iou_runlog,*) ' vmec_flag_name_local = ',                       &
     &      vmec_flag_name_local
         WRITE(iou_runlog,*) ' ier_flag = ', ier_flag
      ENDIF

!-------------------------------------------------------------------------------
!  Clean up loose ends
!-------------------------------------------------------------------------------
         
      IF (PRESENT(delta_iterc)) delta_iterc = iterc - iterc_before

!  Logicals for auxiliarys defined in eq_state set to false
      state % l_def_aux1 = .false.
      state % l_def_aux2 = .false.

!  Cumulative Internal iteration counter
!      iter_es = iter_ha JDH 2010-08-03
      iter_es = iterc
      
      RETURN
      
      END SUBROUTINE eq_step
!*******************************************************************************
! SECTION V.  EQUILIBRIUM CHANGE SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Change the internal state of the equilibrium solver
!
!    For now, specific to the VMEC code.
!-------------------------------------------------------------------------------

      SUBROUTINE eq_change_vmi_varp(state)

!-------------------------------------------------------------------------------
!  USE statements
!-------------------------------------------------------------------------------

! VMEC xstuff - state arrays
      USE xstuff, ONLY: xc_vmec_3 => xc, xcdot_vmec_3 => xcdot,                &
     &   xstore_vmec_3 => xstore

! VMEC input variables. mgrid_file, nfp to reread for extcur change
      USE vmec_input, ONLY:                                                    &
     &   ac, ac_aux_s, ac_aux_f, ai, ai_aux_s, ai_aux_f,                       &
     &   am, am_aux_s, am_aux_f, bloat, rbc, zbs, zbc, rbs,                    &     
     &   extcur, curtor, phiedge, pres_scale,                                  &
     &   mgrid_file, nfp, lfreeb

! VMEC  - icurv array, for checking. currv - for curtor
!   irzloff - for call to profil3d
! VMEC ivac, fsqr, fsqz for ensuring that extcur changes are absorbed
      USE vmec_main, ONLY: icurv, currv, irzloff, ivac, fsqr, fsqz

!  mu0 travels through vparams, and vmec_main to get to readin
      USE stel_constants, ONLY: mu0

!  nv is in vacmod0, used as argument to read_mgrid
      USE vacmod0, ONLY: nv

!  subroutine read_mgrid is in the mgrid_mod module
!  Need mgrid_path_old to get around a test in read_mgrid
      USE mgrid_mod, ONLY: read_mgrid, mgrid_path_old

!-------------------------------------------------------------------------------
!  Declarations
!-------------------------------------------------------------------------------
!  state       Type eq_state. Change the internal state of vmec to match the varp
!  
!  Declare Arguments 
      TYPE (eq_state), INTENT (inout) :: state

!  Declare local variables
!  lextcur_diff      logical to determine if extcur has changed.
!  lrbc_diff         logical to determine if rbc has changed.
!  lzbs_diff         logical to determine if zbs has changed.
!  lzbc_diff         logical to determine if zbc has changed.
!  lrbs_diff         logical to determine if rbs has changed.
!  l_bc_mask         logical array, same bounds as rbc, etc.
!  lscreen_arg       argument for read_mgrid
!  ier_flag_arg      argument for read_mgrid
!  l_reset_xc        argument for profil1d and profil3d
!  ilb1              lower bound of rbc array, first dimension
!  iub1              upper bound of rbc array, first dimension
!  ilb2              lower bound of rbc array, second dimension
!  iub2              upper bound of rbc array, second dimension

      LOGICAL            :: l_reset_xc_3d     = .false.
      LOGICAL, PARAMETER :: l_reset_1d     = .false.
      LOGICAL            :: lextcur_diff, lscreen_arg
      LOGICAL, DIMENSION(:,:), ALLOCATABLE :: l_bc_mask
      LOGICAL, DIMENSION(:), ALLOCATABLE :: l_extcur_mask
      LOGICAL            :: lrbc_diff, lzbs_diff, lzbc_diff, lrbs_diff
      LOGICAL            :: lphiedge_diff
      INTEGER            :: ier_flag_arg, ier1
      INTEGER            :: ilb1, iub1, ilb2, iub2
      REAL(rprec)        :: phiedge_ratio
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'eq_change_vmi_varp: '

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------

!  The auxilliary structures are no longer consistent
      CALL eq_aux_12_undefine(state)

!-------------------------------------------------------------------------------
!  Copy the variable parameters into the VMEC variables
!  Sizes were checked on construction, so don't bother here
!-------------------------------------------------------------------------------
 
      ac = state % varp % ac
      ac_aux_s = state % varp % ac_aux_s
      ac_aux_f = state % varp % ac_aux_f
      ai = state % varp % ai
      ai_aux_s = state % varp % ai_aux_s
      ai_aux_f = state % varp % ai_aux_f
      am = state % varp % am
      am_aux_s = state % varp % am_aux_s
      am_aux_f = state % varp % am_aux_f
      bloat = state % varp % bloat
      curtor = state % varp % curtor
      pres_scale = state % varp % pres_scale

!  See if phiedge has changed, set lphiedge_diff, phiedge_ratio
      lphiedge_diff = .false.
      IF (phiedge .ne. state % varp % phiedge) THEN
         lphiedge_diff = .true.
         phiedge_ratio = one
         IF (phiedge .ne. 0._rprec) THEN
            phiedge_ratio = state % varp % phiedge / phiedge
         ENDIF
         phiedge = state % varp % phiedge
      ENDIF
            
!  See if extcur has changed, set lextcur_diff
      ilb1 = LBOUND(extcur,1)
      iub1 = UBOUND(extcur,1)
      ALLOCATE(l_extcur_mask(ilb1:iub1),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc l_bextcur_mask')

      l_extcur_mask = extcur .ne. state % varp % extcur
      lextcur_diff = any(l_extcur_mask)
      IF (lextcur_diff) extcur = state % varp % extcur

!  See if boundary coefficients have changed, set lrbc_diff, etc.
      ilb1 = LBOUND(rbc,1)
      iub1 = UBOUND(rbc,1)
      ilb2 = LBOUND(rbc,2)
      iub2 = UBOUND(rbc,2)
      ALLOCATE(l_bc_mask(ilb1:iub1,ilb2:iub2),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc l_bc_mask')

      l_bc_mask = rbc .ne. state % varp % rbc
      lrbc_diff = any(l_bc_mask)
      IF (lrbc_diff) rbc = state % varp % rbc

      l_bc_mask = zbs .ne. state % varp % zbs
      lzbs_diff = any(l_bc_mask)
      IF (lzbs_diff) zbs = state % varp % zbs

      l_bc_mask = zbc .ne. state % varp % zbc
      lzbc_diff = any(l_bc_mask)
      IF (lzbc_diff) zbc = state % varp % zbc

      l_bc_mask = rbs .ne. state % varp % rbs
      lrbs_diff = any(l_bc_mask)
      IF (lrbs_diff) rbs = state % varp % rbs
      
!-------------------------------------------------------------------------------
!  Make sure vmec "absorbs" these changes.
!-------------------------------------------------------------------------------
!  Execute relevant statements from readin
!  (For curtor, gets put into currv)
      currv = mu0*curtor              !Convert to Internal units

!  Extra logic and coding for change of phiedge
!  JDH 2010-11-24. Added when VMEC lrfp coding added
!  SPH email --- The reason for the change is that in the new code, instead
!  of phip*(1+lambda) we have phip + lambda, and lambda needs to be scaled
!  by this factor to keep the initial lambda force small.---
!  JDH 2010-11-29. Added coding to als0 rescale xstore
      IF (lphiedge_diff) THEN
         xc_vmec_3(1+2*irzloff:3*irzloff) = phiedge_ratio *                    &
     &      xc_vmec_3(1+2*irzloff:3*irzloff)
         xstore_vmec_3(1+2*irzloff:3*irzloff) = phiedge_ratio *                &
     &      xstore_vmec_3(1+2*irzloff:3*irzloff)
      ENDIF

!  More 'absorption' in profil1d
      CALL profil1d (xc_vmec_3, xcdot_vmec_3, l_reset_1d)

!  2011-01-04 Added boundary coefficients to recon_params
!  Coding to "absorb" boundary coefficient changes into VMEC
!    eq_change_vmi_init_bnd_coeff - coding from readin
!    fsqr, fsqz = one - prevents immediate return
!    profil3d calls - actual changes to xc, xstore
!  ELSE clause insures profil3d gets called
      IF ((.not. lfreeb) .and. (lrbc_diff .or. lzbs_diff .or.                  &
     &   lzbc_diff .or. lrbs_diff)) THEN
         CALL eq_change_vmi_init_bnd_coeff
         fsqr = one
         fsqz = one
         l_reset_xc_3d = .true.
         CALL profil3d(xc_vmec_3(1), xc_vmec_3(1+irzloff),                     &
     &      l_reset_xc_3d)
         CALL profil3d(xstore_vmec_3(1), xstore_vmec_3(1+irzloff),             &
     &      l_reset_xc_3d)
      ELSE
         l_reset_xc_3d = .false.
         CALL profil3d(xc_vmec_3(1), xc_vmec_3(1+irzloff),                     &
     &      l_reset_xc_3d)
      ENDIF

!  Have to reread the mgrid file if the external currents have been changed
!  Have to change mgrid_path_old (module variable in mgrid_mod) to avoid a trap
!  in read_mgrid, that detects if the mgrid file has been read already.
      IF (lextcur_diff) THEN
         mgrid_path_old = ' '
         lscreen_arg = .FALSE.
         ier_flag_arg = 0
         CALL read_mgrid(mgrid_file,extcur,nv, nfp,lscreen_arg,                &
     &      ier_flag_arg)
         CALL assert_eq(0,ier_flag_arg,sub_name // 'bad read_mgrid')

!  Coding added 2007-12-21 JDH
!  Also have to ensure that changed current is used in the first VMEC iteration
!  See VMEC coding in subroutine funct3d.f
!  Want to ensure that coding inside IF statements at lines
!  154 - IF (lfreeb .and. iter2.gt.1 .and. iequi.eq.0) THEN
!    and
!  157  - IF (ivac .ge. 0) THEN
!     and
!  166  - IF (ivac .le. 2) ivacskip = 0
!  gets executed.
!  JDH 2008-01-10
!  Try to prevent execution of 
!  155  - IF ((fsqr + fsqz) .le. 1.e-3_dp) ivac = ivac + 1
!  JDH 2008-01-14 - Did not work without fsqr=fsqz=1.

         ivac = 1
         fsqr = one
         fsqz = one
      ENDIF

      RETURN
      
      END SUBROUTINE eq_change_vmi_varp
      
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

      SUBROUTINE eq_change_vmi_init_bnd_coeff
!
!  Subroutine to (re)initialize various VMEC quantities, when the boundary
!  coefficients (rbc, zbs, zbc, rbs) have been changed.
!  Coding mostly duplicated from VMEC's readin.f
!  Called from eq_change_vmi_varp
!  
!  First version, JDH 2011-01-04
!
      USE stel_kinds
      USE stel_constants
!  NB: vmec_dim and vmec_input are USEd in vmec_main
      USE vmec_input, ONLY: lasym, lfreeb, rbs, zbc, zbs, rbc, mpol,           &
     &   ntor, mfilter_fbdy, nfilter_fbdy
      USE vmec_dim, ONLY: mpol1, ntor1
      USE vmec_main, ONLY: rmn_bdy, zmn_bdy, lthreed, lconm1
      USE vmec_params, ONLY: signgs, rcc, rss, rsc, rcs, zsc, zcs,             &
     &   zcc, zss
!  vmec_main USEd vacmod which USEd vparams
      USE vparams, ONLY: p5 => cp5
      IMPLICIT NONE

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec) :: delta, trc, tzc, rtest, ztest
      REAL(rprec), DIMENSION(:,:), POINTER ::                                  &
     &  rbcc, rbss, rbcs, rbsc, zbcs, zbsc, zbcc, zbss
      REAL(rprec), ALLOCATABLE :: temp(:)
      INTEGER :: m, n, ioff, joff, mj, ni, isgn

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------

      IF (lasym) THEN
!
!       CONVERT TO REPRESENTATION WITH RBS(m=1) = ZBC(m=1)
!
      delta = ATAN( (rbs(0,1) - zbc(0,1))/
     1           (ABS(rbc(0,1)) + ABS(zbs(0,1))) )
      IF (delta .ne. zero) THEN
        DO m = 0,mpol1
          DO n = -ntor,ntor
            trc = rbc(n,m)*COS(m*delta) + rbs(n,m)*SIN(m*delta)
            rbs(n,m) = rbs(n,m)*COS(m*delta) - rbc(n,m)*SIN(m*delta)
            rbc(n,m) = trc
            tzc = zbc(n,m)*COS(m*delta) + zbs(n,m)*SIN(m*delta)
            zbs(n,m) = zbs(n,m)*COS(m*delta) - zbc(n,m)*SIN(m*delta)
            zbc(n,m) = tzc
          ENDDO
        ENDDO
      ENDIF

      ENDIF

!
!     CONVERT TO INTERNAL REPRESENTATION OF MODES
!
!     R = RBCC*COS(M*U)*COS(N*V) + RBSS*SIN(M*U)*SIN(N*V)
!         + RBCS*COS(M*U)*SIN(N*V) + RBSC*SIN(M*U)*COS(N*V)
!     Z = ZBCS*COS(M*U)*SIN(N*V) + ZBSC*SIN(M*U)*COS(N*V)
!         + ZBCC*COS(M*U)*COS(N*V) + ZBSS*SIN(M*U)*SIN(N*V)
!
!
!     POINTER ASSIGNMENTS (NOTE: INDICES START AT 1, NOT 0, FOR POINTERS, EVEN THOUGH
!                          THEY START AT ZERO FOR RMN_BDY)
!     ARRAY STACKING ORDER DETERMINED HERE
!
      rbcc => rmn_bdy(:,:,rcc)
      zbsc => zmn_bdy(:,:,zsc)
      IF (lthreed) THEN
         rbss => rmn_bdy(:,:,rss)
         zbcs => zmn_bdy(:,:,zcs)
      END IF

      IF (lasym) THEN
         rbsc => rmn_bdy(:,:,rsc)
         zbcc => zmn_bdy(:,:,zcc)
         IF (lthreed) THEN
            rbcs => rmn_bdy(:,:,rcs)
            zbss => zmn_bdy(:,:,zss)
         END IF
      ENDIF

      rmn_bdy = 0;  zmn_bdy = 0

      ioff = LBOUND(rbcc,1)
      joff = LBOUND(rbcc,2)

      DO m=0,mpol1
         mj = m+joff
         IF (lfreeb .and. 
     1       (mfilter_fbdy.gt.1 .and. m.gt.mfilter_fbdy)) CYCLE
         DO n=-ntor,ntor
            IF (lfreeb .and. 
     1         (nfilter_fbdy.gt.0 .and. ABS(n).gt.nfilter_fbdy)) CYCLE
            ni = ABS(n) + ioff
            IF (n .eq. 0) THEN
               isgn = 0
            ELSE IF (n .gt. 0) THEN
               isgn = 1
            ELSE
               isgn = -1
            END IF
            rbcc(ni,mj) = rbcc(ni,mj) + rbc(n,m)
            zbsc(ni,mj) = zbsc(ni,mj) + zbs(n,m)
            IF (m .eq. 0) zbsc(ni,mj) = 0

            IF (lthreed) THEN
               rbss(ni,mj) = rbss(ni,mj) + isgn*rbc(n,m)
               zbcs(ni,mj) = zbcs(ni,mj) - isgn*zbs(n,m)
               IF (m .eq. 0) rbss(ni,mj) = 0
            END IF

            IF (lasym) THEN
               rbsc(ni,mj) = rbsc(ni,mj) + rbs(n,m)
               zbcc(ni,mj) = zbcc(ni,mj) + zbc(n,m)
               IF (m .eq. 0) rbsc(ni,mj) = 0
               IF (lthreed) THEN
               rbcs(ni,mj) = rbcs(ni,mj) - isgn*rbs(n,m)
               zbss(ni,mj) = zbss(ni,mj) + isgn*zbc(n,m)
               IF (m .eq. 0) zbss(ni,mj) = 0
               END IF
            END IF

!  JDH commented out IFs below, only result is to WRITE(threed
!            IF (ier_flag_init .ne. norm_term_flag) CYCLE
!            trc = ABS(rbc(n,m)) + ABS(rbs(n,m))
!     1          + ABS(zbc(n,m)) + ABS(zbs(n,m))
!            IF (m .eq. 0) THEN
!               IF (n .lt. 0) CYCLE
!               IF (trc.eq.zero .and. ABS(raxis_cc(n)).eq.zero .and.
!     1             ABS(zaxis_cs(n)).eq.zero) CYCLE
!               WRITE (nthreed,195) n, m, rbc(n,m), rbs(n,m),
!     1                   zbc(n,m), zbs(n,m), raxis_cc(n), raxis_cs(n),
!     2                   zaxis_cc(n), zaxis_cs(n)
!            ELSE
!               IF (trc .eq. zero) CYCLE
!               WRITE (nthreed,195) n, m, rbc(n,m), rbs(n,m),
!     1                   zbc(n,m), zbs(n,m)
!            END IF
         END DO
      END DO
 195  FORMAT(i5,i4,1p,8e12.4)

!
!     CHECK SIGN OF JACOBIAN (SHOULD BE SAME AS SIGNGS)
!
      m = 1
      mj = m+joff
      rtest = SUM(rbcc(1:ntor1,mj))
      ztest = SUM(zbsc(1:ntor1,mj))
      signgs = one
      IF (rtest*ztest .gt. zero) signgs = -one

!
!     CONVERT TO INTERNAL FORM FOR (CONSTRAINED) m=1 MODES
!

      IF (lconm1 .and. (lthreed .or. lasym)) THEN
         ALLOCATE (temp(SIZE(rbcc,1)))
         IF (lthreed) THEN
            mj = 1+joff
            temp = rbss(:,mj)
            rbss(:,mj) = p5*(temp(:) + zbcs(:,mj))
            zbcs(:,mj) = p5*(temp(:) - zbcs(:,mj))
         END IF
         IF (lasym) THEN
            mj = 1+joff
            temp = rbsc(:,mj)
            rbsc(:,mj) = p5*(temp(:) + zbcc(:,mj))
            zbcc(:,mj) = p5*(temp(:) - zbcc(:,mj))
         END IF
         IF (ALLOCATED(temp)) DEALLOCATE (temp)
      END IF
 
      END SUBROUTINE eq_change_vmi_init_bnd_coeff

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

      SUBROUTINE eq_change_vmi_zero_xcdot
!  Subroutine to set the VMEC xcdot array to zero

!-------------------------------------------------------------------------------
!  USE statements
!-------------------------------------------------------------------------------

! VMEC xstuff - state arrays
      USE xstuff, ONLY: xcdot_vmec_4 => xcdot

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------

      xcdot_vmec_4 = zero
      
      END SUBROUTINE eq_change_vmi_zero_xcdot

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

      SUBROUTINE eq_change_vmi_cp_xc(state)
!  Subroutine to copy the xc array from an eq_state to the VMEC internal state

! JDH 2012-04-24
!  Add phiedge change, so that phiedge is consistent with the xc array
!  (otherwise, may mess stuff up with coding in eq_change_vmi_varp)
!-------------------------------------------------------------------------------
!  USE statements
!-------------------------------------------------------------------------------

! VMEC xstuff - state arrays
      USE xstuff, ONLY: xc_vmec_5 => xc
      USE vmec_input, ONLY: phiedge

!-------------------------------------------------------------------------------
!  Declarations
!-------------------------------------------------------------------------------
!  state       Type eq_state.
!  
!  Declare Arguments 
      TYPE (eq_state), INTENT (in) :: state

!  Declare local variables
!  n_xc_vmec      size of the VMEC internal xc array
!  n_xc_state     size of the state xc array
      INTEGER            :: n_xc_vmec, n_xc_state
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'eq_change_cp_xc: '

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------

! Get the array sizes
      n_xc_vmec = SIZE(xc_vmec_5)
      n_xc_state = SIZE(state % xc)
      
      IF (n_xc_vmec .eq. n_xc_state) THEN
         xc_vmec_5 = state % xc
         phiedge = state % varp % phiedge
      ELSE
         WRITE(*,*) 'WARNING from ', sub_name
         WRITE(*,*) ' xc arrays of different sizes'
         WRITE(*,*) ' n_xc_vmec, n_xc_state = ', n_xc_vmec, n_xc_state
         WRITE(*,*) 'END WARNING'
      ENDIF
            
      END SUBROUTINE eq_change_vmi_cp_xc

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

      SUBROUTINE eq_change_vmi_reset_x(state)
!  Subroutine to reset x* variables from an eq_state to the VMEC internal state
!  appropriate for restarting VEC itertions
!     xc_VMEC = xc_STATE
!     xstore_VMEC = xc_STATE
!     xcdot_VMEC = zero

!-------------------------------------------------------------------------------
!  USE statements
!-------------------------------------------------------------------------------

! VMEC xstuff - state arrays
      USE xstuff, ONLY: xc_vmec_6 => xc, xstore_vmec_6 => xstore,              &
     &   xcdot_vmec_6 => xcdot

!-------------------------------------------------------------------------------
!  Declarations
!-------------------------------------------------------------------------------
!  state       Type eq_state.
!  
!  Declare Arguments 
      TYPE (eq_state), INTENT (in) :: state

!  Declare local variables
!  n_xc_vmec      size of the VMEC internal xc array
!  n_xc_state     size of the state xc array
      INTEGER            :: n_xc_vmec, n_xc_state
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'eq_change_reset_x: '

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------

! Get the array sizes
      n_xc_vmec = SIZE(xc_vmec_6)
      n_xc_state = SIZE(state % xc)
      
      IF (n_xc_vmec .eq. n_xc_state) THEN
         xc_vmec_6 = state % xc
         xstore_vmec_6 = state % xc
         xcdot_vmec_6 = zero
      ELSE
         WRITE(*,*) 'WARNING from ', sub_name
         WRITE(*,*) ' xc arrays of different sizes'
         WRITE(*,*) ' n_xc_vmec, n_xc_state = ', n_xc_vmec, n_xc_state
         WRITE(*,*) 'END WARNING'
      ENDIF
            
      END SUBROUTINE eq_change_vmi_reset_x
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

      SUBROUTINE eq_change_vmi_precon2d_save
!  Subroutine to save 2d preconditioning variables from the VMEC internal state

!-------------------------------------------------------------------------------
!  USE statements
!-------------------------------------------------------------------------------

! VMEC xstuff 
       USE precon2d, ONLY: ictrl_prec2d, l_comp_prec2d
       USE vmec_input, ONLY: precon_type, prec2d_threshold
       USE gmres_mod, ONLY: lqmr

!-------------------------------------------------------------------------------
!  Declarations
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------

! save the variables. The _save variables are module variables in eq_interface.
      ictrl_prec2d_save = ictrl_prec2d
      precon_type_save = precon_type
      prec2d_threshold_save = prec2d_threshold
      l_comp_prec2d_save = l_comp_prec2d
      lqmr_save = lqmr
            
      END SUBROUTINE eq_change_vmi_precon2d_save

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

      SUBROUTINE eq_change_vmi_precon2d_restore
!  Subroutine to restore 2d preconditioning variables to the VMEC internal state
!  Should be called AFTER a call to eq_change_vmi_precon2d_save

!-------------------------------------------------------------------------------
!  USE statements
!-------------------------------------------------------------------------------

! VMEC xstuff 
       USE precon2d, ONLY: ictrl_prec2d, l_comp_prec2d
       USE vmec_input, ONLY: precon_type, prec2d_threshold
       USE gmres_mod, ONLY: lqmr

!-------------------------------------------------------------------------------
!  Declarations
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------

! Restore the variables. The _save variables are module variables in eq_interface.
      ictrl_prec2d = ictrl_prec2d_save
      precon_type = precon_type_save
      prec2d_threshold= prec2d_threshold_save
      l_comp_prec2d = l_comp_prec2d_save
      lqmr = lqmr_save
            
      END SUBROUTINE eq_change_vmi_precon2d_restore

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

      SUBROUTINE eq_change_vmi_precon2d_off_same_blocks
!  Subroutine to reset 2d preconditioning variables in the VMEC internal state
!   (1)  Turns off the preconditionning, by setting ictrl_prec2d = 0
!   (2)  Sets l_comp_prec2d = .false, so that the preconditioning blocks
!        will NOT be recomputed the next time VMEC2000 goes through the 
!        prec2d_threshold

!-------------------------------------------------------------------------------
!  USE statements
!-------------------------------------------------------------------------------

! VMEC xstuff 
       USE precon2d, ONLY: ictrl_prec2d, l_comp_prec2d
       USE vmec_input, ONLY: precon_type, prec2d_threshold
       USE gmres_mod, ONLY: lqmr

!-------------------------------------------------------------------------------
!  Declarations
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------

! Restore the variables. The _save variables are module variables in eq_interface.
      ictrl_prec2d = 0
      lqmr = .false.
      l_comp_prec2d = .false.
            
      END SUBROUTINE eq_change_vmi_precon2d_off_same_blocks
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

      SUBROUTINE eq_change_vmi_precon2d_off_new_blocks
!  Subroutine to reset 2d preconditioning variables in the VMEC internal state
!   (1)  Turns off the preconditionning, by setting ictrl_prec2d = 0
!   (2)  Sets l_comp_prec2d = .true., so that the preconditioning blocks
!        WILL be recomputed the next time VMEC2000 goes through the 
!        prec2d_threshold

!-------------------------------------------------------------------------------
!  USE statements
!-------------------------------------------------------------------------------

! VMEC xstuff 
       USE precon2d, ONLY: ictrl_prec2d, l_comp_prec2d
       USE vmec_input, ONLY: precon_type, prec2d_threshold
       USE gmres_mod, ONLY: lqmr

!-------------------------------------------------------------------------------
!  Declarations
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------

! Restore the variables. The _save variables are module variables in eq_interface.
      ictrl_prec2d = 0
      lqmr = .false.
      l_comp_prec2d = .true.
            
      END SUBROUTINE eq_change_vmi_precon2d_off_new_blocks

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

      SUBROUTINE eq_change_vmi_precon2d_stop
!  Subroutine to STOP 2d preconditioning in the VMEC internal state

!-------------------------------------------------------------------------------
!  USE statements
!-------------------------------------------------------------------------------

! VMEC xstuff 
       USE precon2d, ONLY: ictrl_prec2d
       USE vmec_input, ONLY: precon_type, prec2d_threshold
       USE gmres_mod, ONLY: lqmr

!-------------------------------------------------------------------------------
!  Declarations
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------

! Turn off the current preconditioning
      ictrl_prec2d = 0
      lqmr = .false.

! Make sure that preconditioning does NOT get turned back on
      precon_type = 'NONE'; prec2d_threshold = 1.E-30_dp
                  
      END SUBROUTINE eq_change_vmi_precon2d_stop

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

!*******************************************************************************
! SECTION VI.  EQUILIBRIUM GET INFORMATION SUBROUTINES
!*******************************************************************************
!  See FUNCTION eq_interface_get_pressure_s(s_arg), outside the module

!*******************************************************************************
! SECTION VII.  EQUILIBRIUM AUXILIARY CALCULATION SUBROUTINES
!*******************************************************************************



!*******************************************************************************
! SECTION VIII.  EQUILIBRIUM WRITE SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Print out history information
!    Module vmec_history, in VMEC2000/Sources/TimeStep store some FTOL 
!    information after each iteration.
!    This subroutine calls the vmec_hisotry subroutine that prints the
!    information out.
!
!    Specific to the VMEC code.
!-------------------------------------------------------------------------------
      SUBROUTINE eq_history_print

      USE vmec_history

!  Start of executable code

!  write out the vmec_history
      CALL vmec_history_print_flag_on
      CALL vmec_history_print

      RETURN
      END SUBROUTINE eq_history_print

!-------------------------------------------------------------------------------
!  Set integers for vmec_history
!    Use i1 for reconstruction iteration number
!    Use i2 to specify which reconstruction parameter, when calculating jacobian
!
!    Specific to the VMEC code.
!-------------------------------------------------------------------------------
      SUBROUTINE eq_history_set(i1,i2)

      USE vmec_history

!  Declare Arguments 
      INTEGER, OPTIONAL :: i1, i2

!  Start of executable code

      IF (PRESENT(i1)) THEN
         IF (PRESENT(i2)) THEN
            CALL vmec_history_set(i1,i2)  !  both i1 and i2
         ELSE
            CALL vmec_history_set(i1)     !  i1 only
         ENDIF
      ELSE
         IF (PRESENT(i2)) THEN
            CALL vmec_history_set(i2=i2)  !  i2 only
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE eq_history_set

!-------------------------------------------------------------------------------
!  Get integers for vmec_history
!
!    Specific to the VMEC code.
!-------------------------------------------------------------------------------
      SUBROUTINE eq_history_get(i1,i2)

      USE vmec_history

!  Declare Arguments 
      INTEGER :: i1, i2

!  Start of executable code
      CALL vmec_history_get(i1,i2)

      RETURN
      END SUBROUTINE eq_history_get

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

      SUBROUTINE eq_interface_nli_write(input_file_end,istat)
!  Write out a VMEC namelist-input file
!  Modified JDH 2012-06-01

      USE vmec_input
      USE mgrid_mod, ONLY:  curlabel, nextcur
      USE safe_open_mod
      USE write_array_generic
      USE vacfield_mod, ONLY: write_invac
      USE date_and_computer
      USE vmec_params, ONLY: version_
      IMPLICIT NONE

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER(LEN=*) :: input_file_end
      INTEGER, INTENT(inout) :: istat
!  input_file_end    full path to vmec namelist input file, to be written
!  istat             output error indicator

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, DIMENSION(1) :: ins
      INTEGER :: iftol, m, n, i, k, ns0, ns1, iunit,                           &
     &    nextcur_max, iout=23
      CHARACTER(LEN=10) :: date0, time0, zone0
      CHARACTER(LEN=50) :: dateloc, Version
      INTEGER :: imon
      CHARACTER(LEN=80), PARAMETER ::                                          &
     &   banner =                                                              &
     &'! THIS IS WRITTEN BY V3FIT USING VMEC2000, VERSION '

!-----------------------------------------------
!  Start of Executable Code

!  Define dateloc variable, for write at end of file
!   GetComputerInfo is in module date_and_computer
!   DATE_AND_TIME is Fotran intrinsic function
      CALL GetComputerInfo
      CALL DATE_AND_TIME(date0,time0,zone0)
      READ(date0(5:6),'(i2)') imon
      WRITE(dateloc,100) months(imon), date0(7:8), date0(1:4),                 &
     &  time0(1:2),time0(3:4),time0(5:6)
 100  FORMAT('! DATE = ',a3,' ',a2,',',a4,' ',' TIME = ',2(a2,':'),a2)

!  Open the file to write to
      iunit = iout+2
      CALL safe_open(iunit, istat, input_file_end,                             & 
     &	'replace', 'formatted')
      IF (istat .ne. 0) THEN
         istat = -3
         RETURN
      END IF

!  find print out size for NS_ARRAY
      ins = MAXLOC(ns_array)
      ns0 = ns_array(ins(1))

!  find print out size for FTOL_ARRAY
      iftol =0 
      DO WHILE(ftol_array(iftol+1).ne.zero .and. iftol.lt.100)
         iftol = iftol + 1
      END DO

!  find print out size for EXTCUR, onnly when nextcur is zero to start.
!  (NB - redefines nextcur, module variable in mgrid_mod)
      IF (nextcur .eq. 0) THEN
         DO nextcur = SIZE(extcur), 1, -1
            IF (extcur(nextcur) .ne. zero) EXIT
         END DO
         nextcur = MAX(0, nextcur)
      END IF

!
!     Start writing out formatted VMEC data with &INDATA heading
!     so that it can be read in under the indata NAMELIST
!
      WRITE (iunit, '(a7)') '&INDATA'

      WRITE(iunit,*) 
      WRITE(iunit,'("!---- Free / Fixed Boundary ---------------")')
      WRITE(iunit, 1007) 'LFREEB = ', lfreeb
      IF (lfreeb) THEN
         WRITE(iunit,'(2x,3a)') "MGRID_FILE = '", TRIM(mgrid_file), "'"
         IF (nextcur .gt. 0) then
            DO i=1,nextcur 
               WRITE(iunit,'(a,i2.2,a,1p,e22.14,2a)')                          &
     &		"  EXTCUR(",i,") = ", extcur(i),"	  !",TRIM(curlabel(i))
        ENDDO
      ENDIF
      ENDIF

      WRITE(iunit,*) 
      WRITE(iunit,'("!---- Runtime Parameters ---------------")')
      WRITE(iunit,'(2x,a,1p,e10.2)') 'DELT = ', delt
      WRITE(iunit,'(2x,a,i6)') 'NITER = ', niter
      WRITE(iunit,'(2x,a,i4)') 'NSTEP = ', nstep
      WRITE(iunit,'(2x,a,1p,e10.2)') 'TCON0 = ', tcon0
      WRITE(iunit,1007) 'LFORBAL = ', lforbal
      WRITE(iunit, '(2x,a,i4)') 'NVACSKIP = ', nvacskip
      WRITE(iunit,'(2x,a)') 'NS_ARRAY = '
      WRITE(iunit,980) (ns_array(i),i=1,ins(1))
      WRITE(iunit,'(2x,a)') 'FTOL_ARRAY ='
      WRITE(iunit,995) (ftol_array(i),i=1,iftol )

      WRITE(iunit,*) 
      WRITE(iunit,'("!---- Preconditioning Parameters -------------")')
      WRITE(iunit,'(2x,3a)') "PRECON_TYPE = '", TRIM(precon_type),"'" 
      WRITE(iunit,'(2x,a,1p,e14.6)') "PREC2D_THRESHOLD = ",                    &
     &                                prec2d_threshold

      WRITE(iunit,*) 
      WRITE(iunit,'("!---- Grid Parameters -------------")')
      WRITE(iunit,1007) 'LASYM = ', lasym
      WRITE(iunit,'(2x,a,i4)') 'NFP = ', nfp
      WRITE(iunit,'(2(2x,a,i4))') 'MPOL = ', mpol, 'NTOR = ', ntor
      IF (ntheta .gt. 0) WRITE(iunit,'(2x,a,i4)') 'NTHETA = ',ntheta
      IF (nzeta .gt. 0)  WRITE(iunit,'(2x,a,i4)') 'NZETA  = ',nzeta
      WRITE (iunit, 996) phiedge

      WRITE(iunit,*) 
      WRITE(iunit,'("!---- Pressure Parameters -------------")')
      WRITE(iunit, 994) gamma
      WRITE(iunit, '(2x,a,1p,e14.6)') 'BLOAT = ', bloat
      WRITE(iunit, 997) pres_scale
      WRITE(iunit, '(2x,3a)') "PMASS_TYPE = '", pmass_type, "'"      
      WRITE(iunit, '(2x,a,/,2x,3(1pe22.14))')'AM = ', am
      SELECT CASE(TRIM(pmass_type))
         CASE ('Akima_spline','cubic_spline')
            CALL write_array(iunit,'AM_AUX_S',am_aux_s,SIZE(am_aux_s))
            CALL write_array(iunit,'AM_AUX_F',am_aux_f,SIZE(am_aux_s))
      END SELECT

      WRITE(iunit,*) 
      WRITE(iunit,'("!---- Current / Iota Parameters -------------")')
      WRITE(iunit, '(2x,a,i4)') 'NCURR = ', ncurr
      WRITE(iunit, 998) curtor
      WRITE(iunit, '(2x,3a)') "PCURR_TYPE = '", pcurr_type, "'"      
      WRITE(iunit, '(2x,a,/,2x,3(1pe22.14))')'AC = ', ac
      SELECT CASE(TRIM(pcurr_type))
         CASE ('Akima_spline_Ip','Akima_spline_I',                             &
     &           'cubic_spline_Ip','cubic_spline_I')
            CALL write_array(iunit,'AC_AUX_S',ac_aux_s,SIZE(am_aux_s))
            CALL write_array(iunit,'AC_AUX_F',ac_aux_f,SIZE(am_aux_s))
      END SELECT
      WRITE(iunit, '(2x,3a)') "PIOTA_TYPE = '", piota_type, "'"      
      WRITE(iunit, '(2x,a,/,2x,3(1pe22.14))')'AI = ', ai
      SELECT CASE(TRIM(piota_type))
         CASE ('Akima_spline','cubic_spline')
            CALL write_array(iunit,'AI_AUX_S',ai_aux_s,SIZE(am_aux_s))
            CALL write_array(iunit,'AI_AUX_F',ai_aux_f,SIZE(am_aux_s))
      END SELECT

      WRITE(iunit,*) 
      WRITE(iunit,'("!---- Axis, Boundary -------------")')
      IF (mfilter_fbdy.gt.0)                                                   &
     &    WRITE (iunit, '(2x,a,i4)') 'mfilter_fbdy  = ',mfilter_fbdy
      IF (nfilter_fbdy.gt.0)                                                   & 
     &    WRITE (iunit, '(2x,a,i4)') 'nfilter_fbdy  = ',nfilter_fbdy

      CALL write_rbzb (iunit, istat)
      IF (istat .ne. 0) THEN
         istat = -5
         RETURN
      ENDIF
!  End of Namelist
       WRITE (iunit, 1009) 

! Information after end of namelist       
      CALL GetComputerInfo
      Version = TRIM(ADJUSTL(version_))
      WRITE(iunit,*) 
      WRITE(iunit,'(a,x,a,/,a,x,a,/,a,/,a,x,a)')                               &
     &  TRIM(banner),TRIM(Version),                                            &
     &     '! COMPUTER: ', TRIM(computer),                                     &
     &  TRIM(dateloc)
      CLOSE (iunit)
      RETURN

  980 FORMAT(2x,16i5)
  981 FORMAT(2x,a,10(i5))
  990 FORMAT(2x,a,6(1p,e22.14))
  991 FORMAT(1x,a,5(1p,e122.14))
  993 FORMAT(2(2x,a,1p,e22.14))
  994 FORMAT(2x,'GAMMA = ',1p,e14.6)
  995 FORMAT(2x,1p,4e14.6)
  996 FORMAT(2x,'PHIEDGE = ',1p,e22.14)
  997 FORMAT(2x,'PRES_SCALE = ',1p,e22.14)
  998 FORMAT(2x,'CURTOR = ',1p,e22.14)
 1000 FORMAT(5(1p,e22.14))
 1007 FORMAT(2x,a,L1)
 1008 FORMAT(2x,a,L1,2x,a,i2)
 1009 FORMAT("/",a)     

      END  SUBROUTINE eq_interface_nli_write

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

      SUBROUTINE write_rbzb(iunit, istat)
      USE vmec_input
      USE write_array_generic
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: iunit, istat, m, n
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      CHARACTER(LEN=*), DIMENSION(7), PARAMETER :: outfmt =                    &
     &   (/ "(2(a6,i1,a1,i1,a4,1p,e22.14,3x))",                                &
     &      "(2(a6,i1,a1,i2,a4,1p,e21.14,3x))",                                &
     &      "(2(a6,i2,a1,i1,a4,1p,e21.14,3x))",                                &
     &      "(2(a6,i2,a1,i2,a4,1p,e21.14,3x))",                                &
     &      "(2(a6,i3,a1,i1,a4,1p,e21.14,3x))",                                &
     &      "(2(a6,i3,a1,i2,a4,1p,e21.14,3x))",                                &
     &      "(2(a6,i ,a1,i ,a4,1p,e21.14,3x))" /)
      CHARACTER(len=LEN(outfmt(1))) :: outcfmt
      CALL write_array(iunit, 'RAXIS_CC', raxis_cc, ntor+1, low_index=0)
      CALL write_array(iunit, 'ZAXIS_CS', zaxis_cs, ntor+1, low_index=0)
      CALL write_array(iunit, 'RAXIS_CS', raxis_cs, ntor+1, low_index=0)
      CALL write_array(iunit, 'ZAXIS_CC', zaxis_cc, ntor+1, low_index=0)
      DO m = 0, mpol-1
         DO n = -ntor, ntor
            IF ((rbc(n,m).ne.zero) .or. (zbs(n,m).ne.zero)) THEN
!
!     handle formatting for up to 2 digit n by 3 digit m.
!     while this is probably overkill, we at least have to handle up
!     thru 2 digit n by 2 digit m to handle known cases.
!
               outcfmt = outfmt(1)
               IF( m > 9) outcfmt = outfmt(2)
               IF(( n>-10 .and. n<0 ) .or. (n > 9 .and. n < 100)) THEN
                  outcfmt = outfmt(3)
                  IF( m > 9) outcfmt = outfmt(4)
               ELSE IF( n>-100 .and. n< -9 ) THEN
                  outcfmt = outfmt(5)
                  IF( m > 9) outcfmt = outfmt(6)
               ENDIF

               IF( n>= 100 .or. n <= -100 .or. m >=100) THEN
                  outcfmt = outfmt(7)
               ENDIF
             IF(lasym) then
               WRITE (iunit, outcfmt, err=2000)                                &
     &            '  RBC(', n, ',', m, ') = ', rbc(n,m),                       &
     &            '  ZBS(', n, ',', m, ') = ', zbs(n,m),                       &
     &            '  RBS(', n, ',', m, ') = ', rbs(n,m),                       &
     &            '  ZBC(', n, ',', m, ') = ', zbc(n,m)
             ELSE
               WRITE (iunit, outcfmt, err=2000)                                &
     &            '  RBC(', n, ',', m, ') = ', rbc(n,m),                       &
     &            '  ZBS(', n, ',', m, ') = ', zbs(n,m)
             ENDIF
            ENDIF
         END DO
      END DO

 100  FORMAT(a,(1p,4e22.14))

      istat = 0
      RETURN

 2000 istat = -5

      END SUBROUTINE write_rbzb
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

      END MODULE eq_interface

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

!*******************************************************************************
! SECTION IX.     FUNCTIONS & SUBROUTINES OUTSIDE THE MODULE
!*******************************************************************************
!-------------------------------------------------------------------------------
! The function eq_interface_get_pressure_s needs to be OUTSIDE the
! eq_interface module because:
! 1) it is called by eq_get_pressure_s from module eq_T
! 2) Module eq_T is USEd by eq_interface
! 3) This was the first way I could get it to compile and link
!-------------------------------------------------------------------------------
!  Get pressure of the VMEC internal state
!-------------------------------------------------------------------------------
      FUNCTION eq_interface_get_pressure_s(s_arg)
!  FUNCTION to return the pressure of the VMEC internal state, given a value of
!  the radial coordinate s.
      USE stel_kinds
      USE stel_constants, ONLY: mu0
      IMPLICIT NONE
      
!  Declare Arguments
      REAL(rprec)                     :: eq_interface_get_pressure_s
      REAL(rprec), INTENT (in)        :: s_arg

!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'eq_interface_get_pressure_s: '

! Interface to external function pmass
      INTERFACE
         FUNCTION pmass (xx)
         USE stel_kinds
         REAL(rprec) :: xx, pmass
         END FUNCTION pmass
      END INTERFACE
               
!  Start of executable code.

!  One problems here (JDH 2011-07-27)
!  1) VMEC does not really know the pressure, it know the mass.
!        Resolve by ASSUMING the mass is the pressure.
!        To Do - add coding to check gamma, emit warning if mass != pressure

!  Divide by factor of mu0, to convert back to pascals
!  [See multiplication by mu0 in pmass, and comment there:
!  "NOTE: On entry, am is in pascals. pmass internal units are
!   mu0*pascals (B**2 units)"  Quote copied 2011-07-28]

      eq_interface_get_pressure_s = pmass(s_arg) / mu0
      
      RETURN
         
      END FUNCTION eq_interface_get_pressure_s

!*******************************************************************************
! SECTION X.      COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
!  
!  JDH 09-21-2004
!     First version of module
!
!  JDH 11-02-2004
!     First version of eq_step
!
!  JDH 11-26-2004
!     Modified for revised eq_T module
!
!  JDH 11-28-2004
!     Added nextcur access (in module mgrid_mod). Set state logical variables 
!     in eq_step.
!
!  JDH 12-06-2004
!    ns_index issues, added specification of ictrl_array(4) to eq_step, Change
!    size of xc and xcdot if needed.
!
!  JDH 06-28-2006
!    Added iter_es argument to eq_step, equilibrium solver iteration counter
!
!  JDH 07-18-06
!    Added subroutine eq_history_print
!
!  JDH 07-20-06
!    Modified eq_step, so that don't copy varp items into vmec at start.
!    Add subroutine eq_change_vmi_varp
!
!  JDH 07-24-06
!    Modified eq_step a bit more: only step, don't copy from state to vmec
!    internal variables. Changes to eq_change_vmi_varp: change ac also.
!
!  JDH 11-28-2006
!    Changed eq_change_vmi_varp. Changes pres_scale and am (I hope)
!
!  JDH 12-27-2006
!    Changed reset_jdt_flag to reset_jacdt_flag, consistency with latest runvmec.f
!
!  JDH 12-29-2006
!    Added coding to eq_change_vmi_varp for phiedge and extcur
!
!  JDH 2007-12-21
!    Changed test value for extcur change detection (lextcur_diff)
!    from 1.e-9 to 1.e-20, in subroutine eq_change_vmi_varp
!
!     SPH011408: Replace all INTEGER(iprec) with INTEGER
!
!  JDH 2008-01-14 (started 2007-12-21)
!    Additional changes so that extcur change "takes" in VMEC
!
!  JDH 2008-08-04
!    Added subroutines eq_change_vmi_zero_xcdot and eq_change_vmi_cp_xc(state)
!
!  JDH 2009-03-24
!    Changes from Joan Jimenez via Steve Hirshman. Peculiarity of compiler on
!    origin2000. Complaints about repeated name aliases in USE statements in
!    different module subroutines.
!
!  JDH 2009-05-15
!    Set the l_v3fit logical in eq_init_file. Add pcurr_type and pmass_type
!    variables to eq_interface_nli_write.
!
!  JDH 2009-08-27
!    Significant changes to eq_step. Actual coding changes written 2009-07-16
!    Revised and cleaned up code.
!
!  JDH 2009-11-16
!    Slight changes to coding
!
!  JDH 2010-11-24. 
!    Change to eq_change_vmi_varp, extra logic with VMEC lrfp coding. when
!    phiedge changes.
!
!  JDH 2010-12-06. 
!    Changed argument sequence for calls to eq_param_fix_construct and 
!    eq_param_var_construct, in eq_init_file. Corresponding changes in
!    eq_T.
!    Added more _var variables to be copied in eq_change_vmi_varp
!
!  JDH 2011-01-04
!    Changed argument sequence for eq_param_var_construct, in eq_init_file.
!    Corresponding changes in eq_T. Working on adding boundary coefficients
!    (rbc, zbs, zbc, rbs) as reconstruction parameters.
!    Right now, no coding in eq_change_vmi_varp to "absorb" changes in
!    boundary coefficients.
!
!  JDH 2011-01-11
!    Significant rewrite of eq_change_vmi_varp
!
!  JDH 2011-02-07 (Actually added 2011-06-22)
!    Added eq_change_vmi_precon2d_save, _restore, and _reset. Anticipating
!    experimentation with control of 2d preconditioning.
!
!  JDH 2011-02-21 (Actually added 2011-06-22)
!    Modify eq_change_vmi_precon2d*, to add l_comp_prec2d logical
!
!  JDH 2011-06-22
!    Earlier, had added change so that VMEC file names have correct extensions
!
!  JDH 2011-07-28
!    Added eq_interface_get_pressure_s, outside the module
!
!  JDH 2011-09-07
!    Undefine the aux structures as appropriate
!
!  JDH 2011-10-17
!    Add lqmr to variables controlled for eq_change_vmi_precon2d_ subroutines
!
!  JDH 2012-04-24
!    In eq_change_vmi_cp_xc, added change of phiedge, to keep phiedge consistent
!    with the xc array.
!
!  MC and JDH, 2012-07-09. Modified eq_step to add error recovery possibilities.
         
