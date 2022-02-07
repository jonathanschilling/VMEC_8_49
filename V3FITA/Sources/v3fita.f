!******************************************************************************
!  File v3fit.f
!*******************************************************************************
!**********************************************************************
!**                                                                  **
!**   MAIN PROGRAM:  Equilibrium Reconstruction                      **
!**                                                                  **
!**                                                                  **
!**     SUBPROGRAM DESCRIPTION:                                      **
!**          v3fit uses diagnostic data to reconstruct a VMEC        **
!**	     equilibrium that matches the data                           **
!**                                                                  **
!**     CALLING ARGUMENTS: This is a stand alone code                **
!**          xv3fit arg1                                             **
!**          arg1: optional input namelist file name                 **
!**                                                                  **
!**     Author: James D. Hanson, with help from:                     **
!**             Steve Hirshman                                       **
!**             Ed Lazarus                                           **
!**             Lang Lao                                             **
!**             Steve Knowlton                                       **
!**                                                                  **
!**     REFERENCES:                                                  **
!**          (1)                                                     **
!**                                                                  **
!**********************************************************************
!-------------------------------------------------------------------------------
!   DEPENDENCIES
!-------------------------------------------------------------------------------
!
!  This file uses the following modules:
!    stel_kinds
!       located in LIBSTELL/Sources/Modules
!
!    stel_constants
!       located in LIBSTELL/Sources/Modules
!
!    ezcdf
!       located in LIBSTELL/Sources/Ezcdf
!
!    read_wout_mod
!       located in LIBSTELL/Sources/Modules
!
!    safe_open_mod 
!       located in LIBSTELL/Sources/Modules
!
!    v3f_global
!       located in V3FITA/Sources
!       storage for global variables
!
!    recon_param_T
!       located in V3FITA/Sources
!       Type declarations for reconstruction parameters
!
!-------------------------------------------------------------------------------
!   CHANGE HISTORY
!-------------------------------------------------------------------------------
!
!  See Section V, at end of file
!
!-------------------------------------------------------------------------------
!   USAGE
!-------------------------------------------------------------------------------
!
!  Executing 'xv3fit -h' will printout a help message
!
!    INPUT FILES
!
!  There are many input files read in by xv3fit:
! 1) A namelist input file
!      The first argument on the execute line is the name of this file
!      If there is no argument on the execute line, the default name of the
!         namelist input file is v3fit.in
!      For example, 'xv3fita my_nl_file' reads namelist input from a file named
!         my_nl_file
!
! 2) A LIST file
!
! 3) Many mdsig.nc files o one for each diagnostic
!
! 4) A VMEC namelist input file
!
!   OUTPUT FILES
!      
!-------------------------------------------------------------------------------
!   COMMENTS
!-------------------------------------------------------------------------------
!
!  The module v3_global contains the declarations for the globally accessible
!   variables.
!
!*******************************************************************************
!  CODE V3FITA
!    
! SECTION I.	Main Program
! SECTION II.	Initialization Subroutines
!    initialize_read   - called from Main
!       init_command_line       - parse the command line
!       init_main_nli           - read in the main namelist input file
!       init_signal             - read in the signal information
!                                 (9/04 - reads in mdsig files)
!       init_equilibrium        - read in the initial mhd equilibrium state
!    initialize_compute - called from Main
! SECTION III.	TASK SUBROUTINES
!
!
! SECTION IV. Cleanup Subroutines 
!    cleanup           - called from Main (NYI 07-01-2005)
!    cleanup_vmec      - deallocate vmec arrays and write out files
!        (NYI 07-01-2005)
!    help               - (NYI 07-01-2005)
! SECTION V.	SUBROUTINES FOR TESTNG
! SECTION VI.	COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************

!*******************************************************************************
! SECTION I. MAIN PROGRAM
!*******************************************************************************

      PROGRAM v3fit
!**     MAIN PROGRAM:  Equilibrium reconstruction                    **
!**                                                                  **
      USE stel_kinds
      USE stel_constants
!  USE modules that contain variables that should be persistent
      USE v3f_global
      USE gsq_mod
      IMPLICIT NONE

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: i1=0, iero

      REAL(rprec)        :: time_on, time_off
      CHARACTER(len=100) :: v3f_nli_filename_out
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
      CALL second0(time_on)

      CALL init_command_line

      CALL init_main_nli
      
!  vmec_nli_filename from v3f_global
      CALL init_v3f_nli_out(vmec_nli_filename,v3f_nli_filename_out)
      
      WRITE(*,*) ' *** TASK IS ', my_task
      
      SELECT CASE(TRIM(my_task))
      
      CASE('reconstruct_a1')
         CALL init_signal
         CALL init_equilibrium
         CALL task_reconstruct_a1
         CALL eq_interface_nli_write(v3f_nli_filename_out,iero)

      CASE('gsq_jac_test')
         CALL init_signal
         CALL init_equilibrium
         CALL task_gsq_jac_test

      CASE('vmec_v3post') !  Run VMEC, compute signals, mimicking v3post

         CALL init_signal
         CALL init_equilibrium
         CALL task_vmec_v3post

      CASE('vmec_monitor')
         CALL init_signal
         CALL init_equilibrium
         CALL task_vmec_monitor

      CASE('vmec_change_rp')
         CALL init_signal
         CALL init_equilibrium
         CALL task_vmec_change_rp

      CASE('vmec_changep')
         CALL init_signal
         CALL init_equilibrium
         CALL task_vmec_changep
         CALL eq_interface_nli_write(v3f_nli_filename_out,iero)

      CASE('vmec_change_one_p')
         CALL init_signal
         CALL init_equilibrium
         CALL task_vmec_change_one_p
         CALL eq_interface_nli_write(v3f_nli_filename_out,iero)

      CASE('gsq_on_grid')
         CALL init_signal
         CALL init_equilibrium
         CALL task_gsq_on_grid

      CASE('test') !  Generic test - free to use
         CALL task_test
      
      CASE('sxr') ! SXR_chord testing  GJH 2010-01-05
         CALL init_signal
         CALL init_equilibrium
         CALL task_vmec_v3post
                 
      CASE('add_noise_a1') ! Adds noise to observed signals, and reconstructs
         CALL init_signal
         CALL init_equilibrium
         CALL task_add_noise_a1

      CASE('faraday_rotation') !  Run VMEC, compute signals, mimicking v3post, then run faraday_rotation
         CALL init_signal
         CALL init_equilibrium
         CALL task_vmec_v3post
         CALL AJAX_FARADAY_INIT(1,0)

      CASE('f_rotation_and_cm') !  Run VMEC, etc , then run faraday_rotation AND Cotton-Mouton
         CALL init_signal
         CALL init_equilibrium
         CALL task_vmec_v3post
         CALL AJAX_FARADAY_INIT(4,0)

      CASE('faraday_only_no_vmec') !  run faraday rotation stuff ONLY (ie no VMEC, V3POST, etc)
         CALL init_signal_ipsl
         CALL AJAX_FARADAY_INIT(1,1)

      CASE('fr_cm_only_no_vmec') !  run Faraday & Cotton-Mouton stuff ONLY (ie no VMEC, V3POST, etc)
         CALL init_signal_ipsl
         CALL AJAX_FARADAY_INIT(4,1)

      CASE('bfield_grid')
!         CALL init_signal
!         CALL init_equilibrium
         CALL task_bfield_grid

!  2009-03-05 JDH - Deleted tasks 'ray_tracing' and 'rtrace_circle', along with
!     corresponding subroutines task_ray_tracing and task_rtrace_circle

      CASE DEFAULT
         WRITE(*,*) 'V3FITA MAIN, CASE SELECT - DEFAULT'
         WRITE(*,*) 'my_task = ', my_task

      END SELECT

!-----------------------------------------------
!  End of Executable Code
!  Need to write this (10/18/04: SPH)
!-----------------------------------------------
!      CALL cleanup

      CALL second0(time_off)

      WRITE (*,'(/,a,f10.2)') 'V3FIT RUN TIME (s): ',                          &
     &       time_off-time_on
      WRITE (iou_recout,'(/,a,f10.2)') 'V3FIT RUN  TIME (s): ',                &
     &       time_off-time_on
      WRITE (iou_runlog,'(/,a,f10.2)') 'V3FIT RUN  TIME (s): ',                &
     &       time_off-time_on

      END PROGRAM v3fit

!*******************************************************************************
! SECTION II.	Initialization Subroutines
!*******************************************************************************
!
!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------
!
      SUBROUTINE init_command_line

      USE stel_kinds
      USE stel_constants
      USE v3f_global, only: num_cl_args, command_line_args,                    &
     &   main_nli_filename, main_nli_filename_default, l_debug
      USE v3_utilities                                     ! subroutine err_warn
      
      IMPLICIT NONE

!  Declare local variables
      INTEGER           :: i
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'init_command_line: '

!  Declare Format array for help message'
      CHARACTER(len=*), PARAMETER, DIMENSION(7) :: fmt1 = (/                   &
     & '                                                             ',        &
     & 'Help (-h) message for program V3FIT.                         ',        &
     & '  Source of message is subroutine init_command_line          ',        &
     & '  The single command line argument is the file name of the   ',        &
     & '  main namelist input file. If there is no argument, the code',        &
     & '  uses v3fit.in as the main namelist input file name.        ',        &
     & 'End of Help (-h) message for program V3FIT.                  '         &
     &  /) 

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
!----------------------------------------------------------------------
!  Parse the execute line.
!    Fortran 95 does NOT have intrinsics for getting command line arguments
!    The subroutine getcarg (LIBSTELL/Miscel/getcarg.f) has versions of this
!    functionality for WIN32, LINUX, VMS, CRAY machines.
!
!    The equivalent Fortran 2003 intrinsics will be 
!    command_argument_count and get_command_argument
!    See Section 18.0.2 of Metcalf, Reid, & Cohen, Fortran 95/2003 Explained
!----------------------------------------------------------------------
      CALL getcarg(1, command_line_args(1), num_cl_args)
      DO i = 2, MIN(num_cl_args,SIZE(command_line_args))
         CALL getcarg(i, command_line_args(i), num_cl_args)
      END DO
      
!  Choose action based on the number of arguments
!  (Will need a more sophisticated parsing algorithm if the 
!  number of command line argments gets much larger)
      SELECT CASE(num_cl_args)
      CASE(0)
         main_nli_filename = main_nli_filename_default
      
      CASE(1)
         main_nli_filename = TRIM(ADJUSTL(command_line_args(1)))

      CASE DEFAULT
         main_nli_filename = TRIM(ADJUSTL(command_line_args(1)))
         CALL err_warn(sub_name // ' Expected 0 or 1 cl arguments',            &
     &  int=num_cl_args)

      END SELECT
      
!  Look for help flag, print out help message, and quit
      IF (TRIM(ADJUSTL(command_line_args(1))) .eq. '-h') THEN
         DO i = 1,SIZE(fmt1)
            WRITE(*,*) fmt1(i)
         END DO
         STOP '-h flag causes termination of program'
      END IF
      
!  Look for debug flag. Print out indicator, and set global logical
!    (debug flag uses default main NLI file)
      l_debug = .FALSE.
      IF (TRIM(ADJUSTL(command_line_args(1))) .eq. '-d') THEN
         main_nli_filename = main_nli_filename_default
         l_debug = .TRUE.
         CALL err_warn(sub_name // ' Saw -d argument. l_debug set')
      END IF

      RETURN
      END SUBROUTINE init_command_line
!
!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------
!
      SUBROUTINE init_main_nli
!---------------------------------------------------------------------
! FUNCTION: reads in NML variables from V3FIT namelist file
!
! **Updates: added intpol_filename to namelist.  JS 6/27/07
!--------------------------------------------------------------------      
      USE stel_kinds
      USE stel_constants
      USE v3f_global
!SPH010408
      USE safe_open_mod
      USE v3_utilities                                  ! subroutine assert_eq
      USE recon_param_T
      USE pprofile_T
      USE model_T
      
      IMPLICIT NONE

!----------------------------------------------------------------------
! L O C A L Variable Declarations
!----------------------------------------------------------------------
      INTEGER   :: i, istat, ier1, ier2, ier3, ios1, l_rl_fn, l_so_fn,         &
     &   j, imin, imax
      CHARACTER(LEN=300)                 :: mnli_filename = ''

      TYPE(pprofile) :: pp_ne_temp, pp_sxrem_temp, pp_te_temp
      TYPE(eq_state) :: eqstate_temp

      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'init_main_nli: '

! NLI variables
!----------------------------------------------------------------------
!  Definition of rparam array (reconstruction parameters)
!  NB, only dimensioned 10. Fix when necessary.
!     JDH 2007-12-20. Necessary, increased to 100.
!  Don't bother with inputing values. Will get from the state.
      INTEGER                            :: n_rc = 0
      CHARACTER (len=20), DIMENSION(100) :: rc_type = ''
      INTEGER, DIMENSION(100)            :: rc_index = 0
      REAL(rprec), DIMENSION(100)        :: rc_value = zero
      INTEGER                            :: n_rp = 0
      CHARACTER (len=20), DIMENSION(100) :: rp_type = ''
      INTEGER, DIMENSION(100)            :: rp_index = 0
      INTEGER, DIMENSION(100)            :: rp_index2 = 0
      REAL(rprec), DIMENSION(100)        :: rp_vrnc = zero

!  Signal Data Observed (sdo) sigma and weight specification
!  Declare length from v3f_global, parameter na_sdata_o_default
      INTEGER, DIMENSION(na_sdata_o_default) :: sdo_s_spec_imin = 0
      INTEGER, DIMENSION(na_sdata_o_default) :: sdo_s_spec_imax = 0
      INTEGER, DIMENSION(na_sdata_o_default) :: sdo_w_spec_imin = 0
      INTEGER, DIMENSION(na_sdata_o_default) :: sdo_w_spec_imax = 0
      REAL(rprec), DIMENSION(na_sdata_o_default) ::                            & 
     &   sdo_s_spec_floor = zero
      REAL(rprec), DIMENSION(na_sdata_o_default) ::                            &
     &   sdo_s_spec_fraction = zero
      REAL(rprec), DIMENSION(na_sdata_o_default) ::                            & 
     &   sdo_w_spec_weight = zero

!  Model Specification Variables, to be put in to model_a (v3f_global)
      CHARACTER(LEN=20)             :: model_sxrem_type='none'
      CHARACTER(LEN=20)             :: model_te_type='none'
      REAL(rprec)                   :: ne_pp_unit=1.e18
      REAL(rprec)                   :: ne_min=1.e15
      REAL(rprec)                   :: e_pressure_fraction=0.5

! Parameterized Profile Input variables 
!-------------------------------------------------------------------------------
!  electron density ne
      CHARACTER(LEN=20)             :: pp_ne_ptype='none'
      REAL(rprec), DIMENSION(0:20)  :: pp_ne_b
      REAL(rprec), DIMENSION(101)   :: pp_ne_as, pp_ne_af
!  Soft x-ray emissivity
      CHARACTER(LEN=20)             :: pp_sxrem_ptype='none'
      REAL(rprec), DIMENSION(0:20)  :: pp_sxrem_b
      REAL(rprec), DIMENSION(101)   :: pp_sxrem_as, pp_sxrem_af
!  electron temperature
      CHARACTER(LEN=20)             :: pp_te_ptype='none'
      REAL(rprec), DIMENSION(0:20)  :: pp_te_b
      REAL(rprec), DIMENSION(101)   :: pp_te_as, pp_te_af

!----------------------------------------------------------------------
!  Namelist v3fit_main_nli definition

!          File name variables
!  main_nli_filename     File name for main namelist input (v3f_global)
!  mdsig_list_filename   File name for list of MDSIG files (v3f_global)
!  intpol_filename       File name for interferometry/polarimetry (v3f_global)
!  vmec_nli_filename     File name for VMEC namelist input (v3f_global)
!  sxr_nli_filename		 file holding sxr camera information   GJH 2009-08-26 (v3f_global)
!   ^ Eliminate sxr_nli_filename JDH 2011-10-18
!  sxrch_dot_filename	 file holding sxr chord information   JDH 2011-10-18 (v3f_global)
!  thscte_dot_filename	 file holding Thomson Scattering information   JDH 2011-10-23 (v3f_global)
!  ipch_dot_filenam      file holding interferometry-polarimetry information   JDH 2012-03-15 (v3f_global)
!
!       Array allocation sizes, numbers to allocate
!  na_d_desc             d_desc array, diagnostic_desc (v3f_global) 
!  na_g_desc             g_desc array, geometric_desc (v3f_global)
!  na_c_desc             c_desc array, coosig_desc (v3f_global)
!  na_s_desc             s_desc array, signal_desc (v3f_global)
!  na_rparam             rparam array, recon_param (v3f_global)
!
!       Task specification, work variables
!  my_task               Character: to specify task in MAIN select case (v3f_global)
!  i_work                Array of integers, work space (v3f_global)
!  r_work                Array of reals, work space (v3f_global)
!  c_work                Array of character variables, work space (v3f_global)
!
!      Model profile specification
!  pp_ne_ptype           model electron density profile, parameterized profile type
!  pp_ne_b               array of b_coefficients for the electron density profile
!  pp_ne_as              array of as_coefficients electron density splines
!  pp_ne_af              array of af_coefficients electron density splines
!  pp_sxrem_ptype        model soft x-ray emissivity profile, parameterized profile type
!  pp_sxrem_b            array of b_coefficients for the soft x-ray emissivity profile
!  pp_sxrem_as           array of as_coefficients soft x-ray emissivity splines
!  pp_sxrem_af           array of af_coefficients soft x-ray emissivity splines
!  pp_te_ptype           model electron temperature profile, parameterized profile type
!  pp_te_b               array of b_coefficients for the e-temperature profile
!  pp_te_as              array of as_coefficients e-temperature splines
!  pp_te_af              array of af_coefficients e-temperature splines
!
!      Model Specification Variables
!  model_sxrem_type     Character variable to specify how the Soft x-ray 
!                       emissivity will be calculated by the model.
!      'none'             Not part of the model.
!      'pp_sxrem'         Calculated from the parameterized profile pp_sxrem
!      'pp_ne_vmec_p'     Calculated from the pp_ne and the VMEC pressure
!  model_te_type        Character variable to specify how the e-temperature 
!                       will be calculated by the model.
!      'none'             Not part of the model.
!      'pp_te'            Calculated from the parameterized profile pp_te
!      'pp_ne_vmec_p'     Calculated from the pp_ne and the VMEC pressure
!  ne_pp_unit           Scaling/conversion factor for the density. 
!                       ne * ne_pp_unit is assumed to have units of m ** -3.
!  ne_min               Minimum electron density, m ** -3.
!  e_pressure_fraction  Electron pressure fraction of the total pressure
!                       Used when model_te_type = 'pp_ne_vmec_p'
!
!        Reconstruction constraints
!  n_rc                  Number of reconstruction constraints
!  rc_type               Array of reconstruction constraint types
!  rc_index              Array of reconstruction constraint indices
!  rc_value              Array of reconstruction constraint values
!
!       Reconstruction parameter specification
!  n_rp                  Number of reconstruction parameters
!  rp_type               Array of reconstruction parameter types
!  rp_index              Array of reconstruction parameter indices
!  rp_index2             Array of reconstruction parameter 2nd indices
!  rp_vrnc               Array of reconstruction parameter variances
!
!        Signal Data specification
!  n_sdata_o             Number of signal_data observations (v3f_global)
!  iw_sdo_verbose        Integer to control write out of signal_data
!                           observations. (-1, no write, otherwise, verbosity)
!  sdo_data_a            Array of data for observed signals (v3f_global)
!  sdo_sigma_a           Array of sigmas for observed signals (v3f_global)
!  sdo_weight_a          Array of weight for observed signals (v3f_global)
!
!        Signal Data sigma and weight specification
!  sdo_s_spec_imin       Array of minimum index values, for specifying sigmas
!  sdo_s_spec_imax       Array of maximum index values, for specifying sigmas
!  sdo_s_spec_floor      Array of floor values, for specifying sigmas
!  sdo_s_spec_fraction   Array of fractional values, for specifying sigmas
!      floor > 0, sigma = max(floor,fraction * sdo)
!      floor < 0, sigma = sqrt(floor ** 2 + (fraction * sdo) ** 2 )
!  sdo_w_spec_imin       Array of minimum index values, for specifying weights
!  sdo_w_spec_imax       Array of maximum index values, for specifying weights
!  sdo_w_spec_weight     Array of weight values, for specifying weights
!
!      Geometric information
!  r_major_radius        For circular torus geometric (v3g_global)
!  a_minor_radius        For circular torus geometric (v3g_global)
!
!     lif - Limiter_Iso Function - all declared in v3f_global
!       first dimension of arrays is dimensioned na_lif
!
!  n_lif                 number of limiter_iso functions
!  n_phis_lif            array of number of phi values (:)
!  lif_arz               array of r-z polynomial coefficients, (:,0:4,0:4)
!  lif_rc                array of r offset values (:)
!  lif_zc                array of z offset values (:)
!  lif_sigma             array of sigma values (:)
!  lif_phi_degree        array of phi values (:,na_phis_lif)
!  lif_on_edge           array of logical values (:)
!                           true  - edge of plasma is on limiter
!                           false - edge of plasma is within limiter
!
!     COOSIG - Combination Of Other SIGnals NLI Variables
!  n_coosig         Number of new signals of signal_type coosig
!  coosig_indices   2d Integer Array of indices (to list of signals)
!                   for combinations. First index is for which new signal,
!                   second index is for terms in the combination
!  coosig_coeff     2d real array of coefficients for the combination
!                   (a's in the type)
!  coosig_type      Character array of combination types
!                     'sum' - linear combination
!                     'max' - maximum of wieghted signals
!                     'min' - minimum of weighted signals
!  coosig_name      Character array of coosig names
!  coosig_units     Character array of coosig units
!
!    Reconstruction step control NLI variables  
!  nrstep         max number of reconstruction steps to perform
!  dg2_stop       Stopping criterion on change in g^2
!  cut_svd        cutoff value for relative Singular Values
!  cut_eff        cutoff value for expected step efficiency
!  cut_marg_eff   cutoff value for expected marginal step efficiency
!  cut_delta_a    cutoff value for expected step size
!  cut_dg2        cutoff value for expected change in g^2
!  astep_max      maximum allowable normalized step size
!  step_type      character specification of reconstruction step type
!                   'sl' - straight line
!                   'seg' - segmented
!                   'lm' - Levenberg-Marquardt
!----------------------------------------------------------------------
      NAMELIST/v3fit_main_nli/  main_nli_filename, mdsig_list_filename,        &
     &   intpol_filename, vmec_nli_filename, geometric_filename,               &
     &   sxrch_dot_filename, thscte_dot_filename, ipch_dot_filename,           &
     &   na_d_desc, na_g_desc, na_c_desc, na_s_desc, na_rparam,                &
     &   my_task, i_work, c_work, r_work,                                      &

     &   pp_ne_ptype, pp_ne_b, pp_ne_as, pp_ne_af,                             &
     &   pp_sxrem_ptype, pp_sxrem_b, pp_sxrem_as, pp_sxrem_af,                 &
     &   model_sxrem_type,                                                     &
     &   pp_te_ptype, pp_te_b, pp_te_as, pp_te_af,                             &
     &   model_te_type,                                                        &
     &   ne_pp_unit, ne_min, e_pressure_fraction,                              &
     &   n_rc, rc_type, rc_index, rc_value,                                    &
     &   n_rp, rp_type, rp_index, rp_index2, rp_vrnc,                          &
     &   n_sdata_o, iw_sdo_verbose, sdo_data_a, sdo_sigma_a,                   &
     &   r_major_radius, a_minor_radius,                                       &
     &   n_lif, n_phi_lif, lif_arz, lif_rc, lif_zc, lif_sigma,                 & 
     &   lif_phi_degree, lif_on_edge,                                          &
     &   n_coosig, coosig_indices, coosig_coeff, coosig_type,                  &
     &   coosig_name, coosig_units,                                            &
     &   l_zero_xcdot,                                                         & 
     &   sdo_s_spec_imin, sdo_s_spec_imax, sdo_s_spec_floor,                   & 
     &   sdo_s_spec_fraction,                                                  & 
     &   sdo_w_spec_imin, sdo_w_spec_imax, sdo_w_spec_weight,                  &
     &   nrstep, dg2_stop, cut_svd, cut_eff, cut_marg_eff, cut_delta_a,        &
     &   cut_dg2, astep_max, step_type

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

!.........initialize name of (optionally input) ipsl file from V3RFUN..........!
      intpol_filename = ''

!----------------------------------------------------------------------
!-- open input/output units                                          --
!----------------------------------------------------------------------
      CALL safe_open(iou_mnli, istat, TRIM(main_nli_filename),                 &
     &               'old', 'formatted')
      CALL assert_eq(0,istat,sub_name //                                       &
     &   ' Safe_open of main_nli_filename failed')

!----------------------------------------------------------------------
!-- read inputs from namelist                                        --
!----------------------------------------------------------------------
!      READ (iou_mnli, nml=v3fit_main_nli, iostat = ios1)
!      CALL assert_eq(0,ios1,sub_name // ' error reading nli from ' //          &
!     &   main_nli_filename)
! Comment out above, and read in below, so that get better error messages from
! run-time system
      READ (iou_mnli, nml=v3fit_main_nli)
      CLOSE(iou_mnli, iostat=ios1)
      CALL assert_eq(0,ios1,sub_name // ' error closing ' //                   &
     &   main_nli_filename)
      WRITE(*,*) ' *** V3FIT NLI READ FROM ',                                  &
     &    TRIM(ADJUSTL(main_nli_filename))

!----------------------------------------------------------------------
!-- Open the runlog file, recout file (REConstruction OUTput)        --
!----------------------------------------------------------------------
      mnli_filename = TRIM(ADJUSTL(main_nli_filename))
      IF ((mnli_filename(1:4) .eq. 'v3fi') .and.
     &    (mnli_filename(6:6) .eq. '.')          ) THEN
         runlog_filename = 'runlog.' // TRIM(mnli_filename(7:))
         recout_filename = 'recout.' // TRIM(mnli_filename(7:))
      ELSE
         runlog_filename = 'runlog.' // TRIM(mnli_filename) 
         recout_filename = 'recout.' // TRIM(mnli_filename) 
      ENDIF
      CALL safe_open(iou_runlog, istat, TRIM(runlog_filename),                 &
     &               'replace', 'formatted', delim_in='none')
      CALL assert_eq(0,istat,sub_name //                                       &
     &   ' Safe_open of runlog_filename failed')
           CALL safe_open(iou_recout, istat, TRIM(recout_filename),            &
     &               'replace', 'formatted', delim_in='none')
      CALL assert_eq(0,istat,sub_name //                                       &
     &   ' Safe_open of recout_filename failed')

!----------------------------------------------------------------------
!-- Initial Writes to the runlog and recout files                    --
!----------------------------------------------------------------------

      WRITE(iou_runlog,*) 'V3FITA RUN'
      WRITE(iou_runlog,*) '  Namelist Input from file ',                       &
     &   TRIM(main_nli_filename)
      WRITE(iou_runlog,*) '  my_task is ', my_task
      
      INQUIRE(iou_runlog, RECL=i)
      WRITE(iou_runlog,*) " This file's record length is ", i

      WRITE(iou_recout,*) 'V3FITA RUN'
      WRITE(iou_recout,*) '  Namelist Input from file ',                       &
     &  TRIM(main_nli_filename)
      WRITE(iou_recout,*) '  my_task is ', my_task
      
!----------------------------------------------------------------------
!-- Allocate Space for Arrays in v3f_global                          --
!----------------------------------------------------------------------
      na_d_desc = MAX(na_d_desc,na_d_desc_default)
      ALLOCATE(d_desc(na_d_desc),STAT=ier1)
      na_g_desc = MAX(na_g_desc,na_g_desc_default)
      ALLOCATE(g_desc(na_g_desc),STAT=ier2)
      CALL assert_eq(0,ier1,ier2,sub_name // 'alloc 1')
      na_c_desc = MAX(na_c_desc,na_c_desc_default)
      ALLOCATE(c_desc(na_c_desc),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc 1.5')
      na_s_desc = MAX(na_s_desc,na_s_desc_default)
      ALLOCATE(s_desc(na_s_desc),STAT=ier1)
      na_rparam = MAX(na_rparam,na_rparam_default)
      ALLOCATE(rparam(na_rparam),STAT=ier2)
      CALL assert_eq(0,ier1,ier2,sub_name // 'alloc 2')

      na_sdata_o = MAX(na_sdata_o,na_sdata_o_default)
      ALLOCATE(sd_observe(na_sdata_o),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc 3')

!----------------------------------------------------------------------
!-- Initialize recon_param derived types.                            --
!----------------------------------------------------------------------
!  This might better be done in a separate subroutine, or through a separate
!  file, but for now I don't want to proliferate too many input files. JDH 08-21-06
!  Values initialized to zero - will get values later from the state.
!    2007-12-20 JDH. Bad coding practice below. Fudge - change from 10 to 100
      n_rparam = Max(Min(n_rp,100),0)
      DO i = 1,n_rparam
         CALL recon_param_construct(rparam(i),rp_type(i),zero,                 &
     &      rp_index(i),rp_index2(i),rp_vrnc(i))
      END DO

!----------------------------------------------------------------------
!-- Initialize recon_cnstrnts derived types.                         --
!----------------------------------------------------------------------
      CALL recon_cnstrnts_construct(rcnstrnts)
!    2007-12-20 JDH. Bad coding practice below. Fudge - change from 10 to 100
      n_rc = Max(Min(n_rc,100),0)
      DO i = 1,n_rc
         CALL recon_cnstrnts_append(rcnstrnts,rc_type(i),rc_value(i),          &
     &      rc_index(i))
      END DO

!----------------------------------------------------------------------
!-- Easy specification of observed signal data sigmas                --
!----------------------------------------------------------------------
      DO j = 1,SIZE(sdo_s_spec_imin)
         IF (sdo_s_spec_imin(j) .le. 0) EXIT
         imin = MAX(MIN(sdo_s_spec_imin(j),n_sdata_o),1)
         imax = MIN(MAX(sdo_s_spec_imax(j),imin),n_sdata_o)
         DO i = imin,imax
            IF (sdo_sigma_a(i) .eq. 0) THEN
               IF (sdo_s_spec_floor(j) .ge. 0) THEN
                  sdo_sigma_a(i) = MAX(sdo_s_spec_floor(j),                    & 
     &               ABS(sdo_s_spec_fraction(j) * sdo_data_a(i)))
               ELSE
                  sdo_sigma_a(i) = SQRT(sdo_s_spec_floor(j) ** 2 +             & 
     &               (sdo_s_spec_fraction(j) * sdo_data_a(i)) ** 2)
               ENDIF
            ENDIF
         END DO
      END DO
!----------------------------------------------------------------------
!-- Easy specification of observed signal data weights               --
!----------------------------------------------------------------------
      sdo_weight_a = one
      DO j = 1,SIZE(sdo_w_spec_imin)
         IF (sdo_w_spec_imin(j) .le. 0) EXIT
         imin = MAX(MIN(sdo_w_spec_imin(j),n_sdata_o),1)
         imax = MIN(MAX(sdo_w_spec_imax(j),imin),n_sdata_o)
         DO i = imin,imax
            sdo_weight_a(i) = sdo_w_spec_weight(j)
         END DO
      END DO

!----------------------------------------------------------------------
!-- Initialize observed signal data                                  --
!----------------------------------------------------------------------
!  Matching with signal descriptions, and creation of signal_data types
!  is done in init_signal

!----------------------------------------------------------------------
!-- Initialize parameterized profiles and specification variables in the model
!----------------------------------------------------------------------
!  model_a is from v3f_global. pp_ne_temp is local
!  JDH 2011-10-25 - actually construct the model_a

      CALL pprofile_construct(pp_ne_temp,pp_ne_ptype,pp_ne_b,pp_ne_as,         &
     &   pp_ne_af)
!      model_a % pp_ne = pp_ne_temp
      CALL pprofile_write(pp_ne_temp,'init_nli pp_ne',iou_recout)

      CALL pprofile_construct(pp_sxrem_temp,pp_sxrem_ptype,pp_sxrem_b,         &
     &   pp_sxrem_as,pp_sxrem_af)
!      model_a % pp_sxrem = pp_sxrem_temp
      CALL pprofile_write(pp_sxrem_temp,'init_nli pp_sxrem',iou_recout)

      CALL pprofile_construct(pp_te_temp,pp_te_ptype,pp_te_b,                  &
     &   pp_te_as,pp_te_af)
!      model_a % pp_te = pp_te_temp
      CALL pprofile_write(pp_te_temp,'init_nli pp_te',iou_recout)
     
!      model_a % model_sxrem_type = model_sxrem_type
!      model_a % model_te_type = model_te_type

      WRITE(iou_recout,*) 
      WRITE(iou_recout,*) '  model_sxrem_type is ',                            &
     &   TRIM(model_sxrem_type)
      WRITE(iou_recout,*) '  model_te_type is ', TRIM(model_te_type)
      WRITE(iou_recout,1000) '  ne_pp_unit is ', ne_pp_unit
      WRITE(iou_recout,1000) '  ne_min is ', ne_min
      WRITE(iou_recout,1000) '  e_pressure_fraction is ',                         &
     &   e_pressure_fraction
1000  FORMAT(a,2x,es12.5)

!  Create the model. The eq_state is a temporary variable. The initialization 
!  with an actual state is done with the call
!      CALL eq_init_file(vmec_nli_filename,model_a % eqstate)
!  in subroutine init_equilibrium

      CALL model_construct(model_a,eqstate = eqstate_temp,                     &
     &   pp_ne = pp_ne_temp,                                                   &
     &   pp_sxrem = pp_sxrem_temp,                                             &
     &   pp_te = pp_te_temp,                                                   &
     &   model_sxrem_type = model_sxrem_type,                                  &
     &   model_te_type = model_te_type,                                        &
     &   ne_pp_unit = ne_pp_unit,                                              &
     &   ne_min = ne_min,                                                      &
     &   e_pressure_fraction = e_pressure_fraction)

      CALL pprofile_destroy(pp_ne_temp)
      CALL pprofile_destroy(pp_sxrem_temp)
      CALL pprofile_destroy(pp_te_temp)

      RETURN
      END SUBROUTINE init_main_nli

!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------
!
      SUBROUTINE init_v3f_nli_out(name_in,name_out)
!--------------------------------------------------------------------------
!  Subroutine to initialize generate the name for the output NLI file
!  JDH 2011-09-02 - First Version
!
!  The input V3FIT NLI file name may be a complete path
!  We want the output file name to NOT overwrite the input file name
!  So, we get rid of anything that is directory specification
!  and we append an identifier

      USE stel_kinds
      
      IMPLICIT NONE

!----------------------------------------------------------------------
! A R G U M E N T  Variable Declarations
!----------------------------------------------------------------------
!  name_in       Input file name - may have complete path
!  name_out      output file name - to use
      CHARACTER(len=*), INTENT(in)      :: name_in
      CHARACTER(len=*), INTENT(inout)   :: name_out
      
!----------------------------------------------------------------------
! L O C A L Variable Declarations
!----------------------------------------------------------------------
!  suffix        Suffix to append to generate the name_out
!  len_suffix    Length of the suffix string
!  len_in        Length of name_in
!  len_out       Length of name_out

      CHARACTER(len=7)             :: suffix = '_v3fout'
      CHARACTER(len=LEN(name_in))  :: name_in_local
      CHARACTER(len=LEN(name_out)) :: name_out_temp
      INTEGER                      :: len_suffix, len_in, len_out
      INTEGER                      :: last_slash, len_trim_out_temp
      INTEGER                      :: len_save
      
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

      len_suffix = LEN(suffix)
      len_in = LEN(name_in)
      len_out = LEN(name_out)
      
      name_in_local = name_in
      last_slash = INDEX(name_in_local,'/',.TRUE.)
      
!  Eliminate any directories in the path.
      name_out_temp = name_in_local(last_slash + 1:)
      
      len_trim_out_temp = LEN_TRIM(name_out_temp)
      len_save = MIN(len_trim_out_temp,len_out - len_suffix)
      
      name_out = name_out_temp(1:len_save) // suffix

      RETURN
      END SUBROUTINE init_v3f_nli_out
!
!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------
!
      SUBROUTINE init_signal
!--------------------------------------------------------------------------
!  Subroutine to initialize the signals
!  Reads in mdsig files, written by V3RFUN - done in init_signal_mddc
!
! *modified to initialize ipsl signals as well.  JMS 6/27/07
!  JDH 2008-01-20 - Clean up, eliminate some extraneous writes.
!  GJH 2010-01-15 - Added SXR initialization - done in init_signal_sxrch
!---------------------------------------------------------------------------
!      USE ezcdf
!      USE bsc_cdf
      USE diagnostic_T
      USE diagnostic_cdf
      USE signal_T
      USE signal_cdf
!SPH010408
      USE safe_open_mod
      USE v3_utilities
      USE v3f_global
      USE gsq_mod
      
      IMPLICIT NONE

!----------------------------------------------------------------------
! L O C A L Variable Declarations
!----------------------------------------------------------------------
!  mdsig_filename_pp     _pp means plus path. Path from mdsig_list_path
!  mdsig_list_path       path to mdsig files. Extracted from mdsig_list_filename
!  isig                  Counter for signals
!  idiag                 Counter for diagnostics
!  imddc                 Counter for mddc's
!  idc                   Counter read from mdsig_list files
      
! GJH 2010-01-21 INTEGER imddc declaration ?? JDH 2010-06-11 IS THIS USED?

      INTEGER              :: isig, idiag, imddc
      INTEGER              :: i, istat, ic, mem, na_min, idc
      INTEGER              :: ier1, ilastslash
      TYPE (bsc_coil)             :: temp_coil
      CHARACTER(len=300)          :: mdsig_filename, mdsig_filename_pp,        &
     &   mdsig_list_path
      CHARACTER(len=80)           :: c_ident, s_ident
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'init_signal: '

!----------------------------------------------------------------------
!      Variables from v3f_global, used in this subroutine:
!  mdsig_list_filename
!  iou_mdsig_list
!  sxrch_dot_filename
!  ipch_dot_filename
!----------------------------------------------------------------------

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
      WRITE(*,*) ' *** IN SUBROUTINE ', sub_name
      isig = 0

      mdsig_list_filename = ADJUSTL(mdsig_list_filename)
      intpol_filename     = ADJUSTL(intpol_filename)
      sxrch_dot_filename  = ADJUSTL(sxrch_dot_filename) 
      ipch_dot_filename  = ADJUSTL(ipch_dot_filename) 
!----------------------------------------------------------------------
!--  Diagnostics Section                                             --
!----------------------------------------------------------------------
      idiag = 0
!----------------------------------------------------------------------
!--  mddc (Magnetic Diagnostics) SubSection                          --
!----------------------------------------------------------------------

      CALL init_signal_mddc(isig,idiag)
!----------------------------------------------------------------------
!--  END mddc (Magnetic Diagnostics) SubSection                      --
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!--  sxrch (Soft X-Ray Diagnostics) SubSection  ! GJH 2010-01-15     --
!----------------------------------------------------------------------

      IF (TRIM(sxrch_dot_filename) .NE. '') THEN
         CALL init_signal_sxrch(isig,idiag)    
      ENDIF

!----------------------------------------------------------------------
!--  END sxrch (Soft X-Ray Diagnostics) SubSection                   --
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!--  ipch (Interferometry-Polarimetry CHord) SubSection              --
!----------------------------------------------------------------------

      IF (TRIM(ipch_dot_filename) .NE. '') THEN
         CALL init_signal_ipch(isig,idiag)    
      ENDIF
!----------------------------------------------------------------------
!--  END ipch (Interferometry-Polarimetry CHord) SubSection          --
!----------------------------------------------------------------------

!----------------------------------------------------------------------
!--  thscte (Thomson Scattering Te) SubSection                       --
!----------------------------------------------------------------------

      IF (TRIM(thscte_dot_filename) .NE. '') THEN
         CALL init_signal_thscte(isig,idiag)    
      ENDIF
!----------------------------------------------------------------------
!--  END thscte (Thomson Scattering Te) SubSection                   --
!----------------------------------------------------------------------

!      IPSL Deprecated 2012-03-15
!----------------------------------------------------------------------
!--  ipsl (Interferometry/Polarimetry Straight Line) SubSection      --
!----------------------------------------------------------------------
      IF ( TRIM(intpol_filename) .NE. '' ) THEN
        CALL init_signal_ipsl
      END IF

!----------------------------------------------------------------------
!--  END Diagnostics Section  –  START Geometric Section             --
!----------------------------------------------------------------------

      CALL init_signal_geometric
 
!----------------------------------------------------------------------
!--  END Geometric Section                                           --
!-–  START Coosig (Combination Of Other Signals) Section             --
!----------------------------------------------------------------------

      CALL init_signal_coosig
  
!----------------------------------------------------------------------
!-- Initialize observed signal data                                  --
!----------------------------------------------------------------------
!  data and sigma values for observed signals were NLI in main NLI
!  Communicated to here in arrays in v3f_global
!  Here we match with signal descriptions, and create signal_data

!  Assumption when there are geometrics
!    1) NLI sdo have been input appropriately, enough for diagnostics and geometrics
!    2) sdo_sigmas are important
!    3) sdo_data should be zero
!  Test and enforce some of these assumptions

      IF (n_sdata_o .GT. 0) THEN
      
         IF (n_sdata_o .ne. n_s_desc) THEN   ! Problem with #s of signals, data
            WRITE(*,*) 'In ', sub_name, 
     &       ' Problem - n_sdata_o .ne. n_s_desc' ,n_sdata_o, n_s_desc
            DO ic = 1, n_s_desc
               WRITE(s_ident,'("ic=",i3)') ic
               CALL signal_write(s_desc(ic),s_ident)
            END DO
         CALL assert_eq(n_sdata_o,n_s_desc,sub_name //                         &
     &       ' problem with n_sdata_o, n_s_desc ') 
         ENDIF
         
         DO ic = 1,n_s_desc
            IF (s_desc(ic) % s_type .eq. 'geometric') THEN
               sdo_data_a(ic) = zero
            ENDIF
            CALL signal_construct(sd_observe(ic),s_desc(ic),                   &
     &         'observation',sdo_data_a(ic),sdo_sigma_a(ic),                   &
     &          sdo_weight_a(ic))
            IF (iw_sdo_verbose .ge. 0) THEN
               WRITE(c_ident,*) ' ic = ',ic
               CALL signal_write(sd_observe(ic),c_ident,                       &
     &            verbose = iw_sdo_verbose)
            END IF
         END DO
      ELSE
         WRITE(*,*) 'No observed signal data initialized'
      ENDIF

!----------------------------------------------------------------------
!-- Initialize gsq_mod weights                                --
!----------------------------------------------------------------------
!  10-03-06. Allocation and assignment now done in gsq_initialize

      RETURN
      END SUBROUTINE init_signal
!
!----------------------------------------------------------------------

!======================================================================
      SUBROUTINE init_signal_sxrch(isig,idiag)
!----------------------------------------------------------------------
! Initializes the sxr_chord diagnostic signals
! stores information in the global d_desc and s_desc arrays
!----------------------------------------------------------------------
      
      USE sxrch_T
      USE v3f_global
!      USE sxr_routines
      USE sxrch_dot
      USE signal_T
      
      IMPLICIT NONE
!----------------------------------------------------------------------
! A R G U M E N T  Variable Declarations
!----------------------------------------------------------------------
!  isig                  Counter for signals
!  idiag                 Counter for diagnostics
!----------------------------------------------------------------------
        INTEGER, INTENT(inout)              :: isig, idiag
!----------------------------------------------------------------------
! other definitions
!  sxr_chords            holds the sxr chord info until transferred to
!                        signal_desc
!  n_chords              number of sxr_chords
!  ic                    generic counter
!  mem                   memory
!  isigin, idiagin       holders for input isig and idiag
!----------------------------------------------------------------------
      
        TYPE (sxrch_desc),DIMENSION(:),ALLOCATABLE :: sxr_chords
        INTEGER :: na_sxr_chords, n_chords
        INTEGER :: ier1
        INTEGER :: ic
        INTEGER :: mem
        INTEGER :: isigin, idiagin
        
        CHARACTER(LEN=*),PARAMETER :: d_type='sxrch'
        CHARACTER(LEN=30)           :: s_name
        CHARACTER(LEN=100)          :: l_name
        CHARACTER(LEN=100)           :: units
        CHARACTER(LEN=*),PARAMETER :: sub_name='init_signal_sxrch'
        
        REAL (rprec), PARAMETER           :: sigma_default=0.1
        
!----------------------------------------------------------------------
!                       Start of Executable Code
!----------------------------------------------------------------------

!----------------------------------------------------------------------
! this is where-am-I stuff and can be eliminated later
!----------------------------------------------------------------------
      WRITE(*,*) ' *** IN SUBROUTINE ', sub_name
      isigin=isig
      idiagin=idiag
!      WRITE(*,*) '----------------------------------'
!      WRITE(*,*)'Camera file - ',TRIM(sxr_nli_filename)
 
!----------------------------------------------------------------------
! Allocate the sxr_chords array. Read in the sxrch_dot file
!----------------------------------------------------------------------

      na_sxr_chords = SIZE(d_desc) - idiag
      ALLOCATE(sxr_chords(na_sxr_chords),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc sxr_chords')
      
      CALL sxrch_dot_read(sxrch_dot_filename,sxr_chords,n_chords)
!      CALL init_sxr_chords(sxr_nli_filename,sxr_chords)   ! from sxr_routines
!      n_chords=SIZE(sxr_chords)
        
!----------------------------------------------------------------------------        
! construct the d_desc arrays to hold the diagnostic information
! and store the sxrch info in them
!----------------------------------------------------------------------------                
      DO ic=1,n_chords ! for all sxr chords
         idiag=idiag+1
         s_name=sxr_chords(ic) % chord_name
         l_name= ''
         units='arb'
         CALL diagnostic_desc_construct_sxrch(d_desc(idiag),d_type,            &
     &      s_name,l_name, units,sigma_default,sxr_chords(ic))

!         WRITE(96,*) ic, s_name, sxr_chords(ic) % Ri,                          &
!     &      sxr_chords(ic) % Zi,                                               &
!     &      sxr_chords(ic) % Rf, sxr_chords(ic) % Zf            !  DEBUG  
!
      END DO

!----------------------------------------------------------------------------
! set the global variables for:
! the number of sxr diagnostic chords
! the number of sxr signals
! and the total number of diagnostics
!----------------------------------------------------------------------------        
      n_d_desc_sxrch = n_chords
      n_s_desc_sxrch = n_chords
      n_d_desc      = n_d_desc + n_d_desc_sxrch
        
      CALL assert_eq(idiag,n_d_desc,                                           &
     &   sub_name // 'idiag ne n_d_desc') 
!----------------------------------------------------------------------------
! Construct the signal corresponding to sxrch diagnostic
!----------------------------------------------------------------------------
!  for all of the chords loop over
       
      DO ic= n_s_desc+1, n_s_desc + n_chords
         isig = ic
         idiag = ic
         CALL signal_construct(s_desc(isig),'diagnostic',                      &
     &      d_desc(idiag) % s_name,d_desc(idiag) % l_name,                     &
     &      d_desc(idiag) % units,d_desc(idiag))
      END DO
      n_s_desc=n_s_desc+n_s_desc_sxrch
      
!  Clean Up
      DEALLOCATE(sxr_chords)
      
      CALL getmem(mem)
      WRITE(*,*) sub_name,' Current Memory Usage (kB) is ',mem
     
      END SUBROUTINE init_signal_sxrch
      
!======================================================================
      SUBROUTINE init_signal_ipch(isig,idiag)
!----------------------------------------------------------------------
! Initializes the ip_chord diagnostic signals
! stores information in the global d_desc and s_desc arrays
!----------------------------------------------------------------------
      
      USE ipch_T
      USE v3f_global
!      USE ip_routines
      USE ipch_dot
      USE signal_T
      
      IMPLICIT NONE
!----------------------------------------------------------------------
! A R G U M E N T  Variable Declarations
!----------------------------------------------------------------------
!  isig                  Counter for signals
!  idiag                 Counter for diagnostics
!----------------------------------------------------------------------
        INTEGER, INTENT(inout)              :: isig, idiag
!----------------------------------------------------------------------
! other definitions
!  ip_chords            holds the ip chord info until transferred to
!                        signal_desc
!  n_chords              number of ip_chords
!  ic                    generic counter
!  mem                   memory
!  isigin, idiagin       holders for input isig and idiag
!----------------------------------------------------------------------
      
        TYPE (ipch_desc),DIMENSION(:),ALLOCATABLE :: ip_chords
        INTEGER :: na_ip_chords, n_chords
        INTEGER :: ier1
        INTEGER :: ic
        INTEGER :: mem
        INTEGER :: isigin, idiagin
        
        CHARACTER(LEN=*),PARAMETER :: d_type='ipch'
        CHARACTER(LEN=30)           :: s_name
        CHARACTER(LEN=100)          :: l_name
        CHARACTER(LEN=100)           :: units
        CHARACTER(LEN=*),PARAMETER :: sub_name='init_signal_ipch'
        
        REAL (rprec), PARAMETER           :: sigma_default=0.1
        
!----------------------------------------------------------------------
!                       Start of Executable Code
!----------------------------------------------------------------------

!----------------------------------------------------------------------
! this is where-am-I stuff and can be eliminated later
!----------------------------------------------------------------------
      WRITE(*,*) ' *** IN SUBROUTINE ', sub_name
      isigin=isig
      idiagin=idiag
 
!----------------------------------------------------------------------
! Allocate the ip_chords array. Read in the ipch_dot file
!----------------------------------------------------------------------

      na_ip_chords = SIZE(d_desc) - idiag
      ALLOCATE(ip_chords(na_ip_chords),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc ip_chords')
      
      CALL ipch_dot_read(ipch_dot_filename,ip_chords,n_chords)
        
!----------------------------------------------------------------------------        
! construct the d_desc arrays to hold the diagnostic information
! and store the ipch info in them
!----------------------------------------------------------------------------                
      DO ic=1,n_chords ! for all ip chords
         idiag=idiag+1
         s_name=ip_chords(ic) % chord_name
         l_name= ''
         units='arb'
         CALL diagnostic_desc_construct_ipch(d_desc(idiag),d_type,             &
     &      s_name,l_name, units,sigma_default,ip_chords(ic))

      END DO

!----------------------------------------------------------------------------
! set the global variables for:
! the number of ip diagnostic chords
! the number of ip signals
! and the total number of diagnostics
!----------------------------------------------------------------------------        
      n_d_desc_ipch = n_chords
      n_s_desc_ipch = n_chords
      n_d_desc      = n_d_desc + n_d_desc_ipch
        
      CALL assert_eq(idiag,n_d_desc,                                           &
     &   sub_name // 'idiag ne n_d_desc') 
!----------------------------------------------------------------------------
! Construct the signal corresponding to ipch diagnostic
!----------------------------------------------------------------------------
!  for all of the chords loop over
       
      DO ic= n_s_desc+1, n_s_desc + n_chords
         isig = ic
         idiag = ic
         CALL signal_construct(s_desc(isig),'diagnostic',                      &
     &      d_desc(idiag) % s_name,d_desc(idiag) % l_name,                     &
     &      d_desc(idiag) % units,d_desc(idiag))
      END DO
      n_s_desc=n_s_desc+n_s_desc_ipch
      
!  Clean Up
      DEALLOCATE(ip_chords)
      
      CALL getmem(mem)
      WRITE(*,*) sub_name,' Current Memory Usage (kB) is ',mem
     
      END SUBROUTINE init_signal_ipch

!======================================================================
      SUBROUTINE init_signal_thscte(isig,idiag)
!----------------------------------------------------------------------
! Initializes the Thomson Scattering diagnostic signals
! stores information in the global d_desc and s_desc arrays
!----------------------------------------------------------------------
      
      USE thscte_T
      USE v3f_global
!      USE sxr_routines
      USE thscte_dot
      USE signal_T
      
      IMPLICIT NONE
!----------------------------------------------------------------------
! A R G U M E N T  Variable Declarations
!----------------------------------------------------------------------
!  isig                  Counter for signals
!  idiag                 Counter for diagnostics
!----------------------------------------------------------------------
        INTEGER, INTENT(inout)              :: isig, idiag
!----------------------------------------------------------------------
! other definitions
!  thscte_pts            holds the sxr chord info until transferred to
!                        signal_desc
!  n_chords              number of thscte_pts
!  ic                    generic counter
!  mem                   memory
!  isigin, idiagin       holders for input isig and idiag
!----------------------------------------------------------------------
      
        TYPE (thscte_desc),DIMENSION(:),ALLOCATABLE :: thscte_pts
        INTEGER :: na_thscte_pts, n_chords
        INTEGER :: ier1
        INTEGER :: ic
        INTEGER :: mem
        INTEGER :: isigin, idiagin
        
        CHARACTER(LEN=*),PARAMETER :: d_type='thscte'
        CHARACTER(LEN=30)           :: s_name
        CHARACTER(LEN=100)          :: l_name
        CHARACTER(LEN=100)           :: units
        CHARACTER(LEN=*),PARAMETER :: sub_name='init_signal_thscte'
        
        REAL (rprec), PARAMETER           :: sigma_default=0.1
        
!----------------------------------------------------------------------
!                       Start of Executable Code
!----------------------------------------------------------------------

!----------------------------------------------------------------------
! this is where-am-I stuff and can be eliminated later
!----------------------------------------------------------------------
      WRITE(*,*) ' *** IN SUBROUTINE ', sub_name
      isigin = isig
      idiagin = idiag
 
!----------------------------------------------------------------------
! Allocate the thscte_pts array. Read in the thscte_dot file
!----------------------------------------------------------------------

      na_thscte_pts = SIZE(d_desc) - idiag
      ALLOCATE(thscte_pts(na_thscte_pts),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc thscte_pts')
      
      CALL thscte_dot_read(thscte_dot_filename,thscte_pts,n_chords)
        
!----------------------------------------------------------------------------        
! construct the d_desc arrays to hold the diagnostic information
! and store the thscte info in them
!----------------------------------------------------------------------------                
      DO ic=1,n_chords ! for all sxr chords
         idiag = idiag+1
         s_name = thscte_pts(ic) % chord_name
         l_name= ''
         units='arb'
         CALL diagnostic_desc_cnstrct_thscte(d_desc(idiag),d_type,             &
     &      s_name,l_name, units,sigma_default,thscte_pts(ic))

      END DO

!----------------------------------------------------------------------------
! set the global variables for:
! the number of sxr diagnostic chords
! the number of sxr signals
! and the total number of diagnostics
!----------------------------------------------------------------------------        
      n_d_desc_thscte = n_chords
      n_s_desc_thscte = n_chords
      n_d_desc      = n_d_desc + n_d_desc_thscte
        
      CALL assert_eq(idiag,n_d_desc,                                           &
     &   sub_name // 'idiag ne n_d_desc') 
!----------------------------------------------------------------------------
! Construct the signal corresponding to thscte diagnostic
!----------------------------------------------------------------------------
!  for all of the chords loop over
       
      DO ic= n_s_desc+1, n_s_desc + n_chords
         isig = ic
         idiag = ic
         CALL signal_construct(s_desc(isig),'diagnostic',                      &
     &      d_desc(idiag) % s_name,d_desc(idiag) % l_name,                     &
     &      d_desc(idiag) % units,d_desc(idiag))
      END DO
      n_s_desc=n_s_desc+n_s_desc_thscte
      
!  Clean Up
      DEALLOCATE(thscte_pts)
      
      CALL getmem(mem)
      WRITE(*,*) sub_name,' Current Memory Usage (kB) is ',mem
     
      END SUBROUTINE init_signal_thscte

!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------
      SUBROUTINE init_signal_mddc(isig,idiag)
!--------------------------------------------------------------------------
!  Subroutine to initialize the mddc diagnostic signals
!  Reads in mdsig files, written by V3POST 
!
!!  JDH 2009-03-13 - Separated out of init_signal
!---------------------------------------------------------------------------
      USE diagnostic_T
      USE diagnostic_cdf
      USE signal_T
      USE signal_cdf
!SPH010408
      USE safe_open_mod
      USE v3_utilities
      USE v3f_global
      USE gsq_mod
      
      IMPLICIT NONE

!----------------------------------------------------------------------
! A R G U M E N T  Variable Declarations
!----------------------------------------------------------------------
!  isig                  Counter for signals
!  idiag                 Counter for diagnostics

      INTEGER, INTENT(inout)              :: isig, idiag

!----------------------------------------------------------------------
! L O C A L Variable Declarations
!----------------------------------------------------------------------
!  mdsig_filename_pp     _pp means plus path. Path from mdsig_list_path
!  mdsig_list_path       path to mdsig files. Extracted from mdsig_list_filename
!  imddc                 Counter for mddc's
!  idc                   Counter read from mdsig_list files

      INTEGER              :: imddc
      INTEGER              :: i, istat, ic, mem, na_min, idc
      INTEGER              :: ier1, ilastslash
      TYPE (bsc_coil)             :: temp_coil
      CHARACTER(len=300)          :: mdsig_filename, mdsig_filename_pp,        &
     &   mdsig_list_path
      CHARACTER(len=80)           :: c_ident
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'init_signal_mddc: '

!----------------------------------------------------------------------
!      Variables from v3f_global, used in this subroutine:
!  mdsig_list_filename
!  iou_mdsig_list
!----------------------------------------------------------------------

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
      WRITE(*,*) ' *** IN SUBROUTINE ', sub_name
!----------------------------------------------------------------------
!--  mddc (Magnetic Diagnostics)                                     --
!----------------------------------------------------------------------

!  Check file name, to see if need to be here.
      IF ( TRIM(mdsig_list_filename) .EQ. '' ) RETURN

!-- open input/output units                                          --
      CALL safe_open(iou_mdsig_list, istat, TRIM(mdsig_list_filename),         &
     &               'old', 'formatted')
      CALL assert_eq(0,istat,sub_name //                                       &
     &   ' Safe_open of mdsig_list_filename failed')

!-- Extract the mdsig_list_path                                      --
      ilastslash = INDEX(TRIM(mdsig_list_filename),'/',back=.true.)
      mdsig_list_path = mdsig_list_filename(1:ilastslash)

!-- Loop over files in the list file                                 --

      imddc = 0
      na_min = MIN(na_d_desc,na_s_desc)
      DO  ! loop over magnetic diagnostic signals
         imddc = imddc + 1
         idiag = idiag + 1
         READ(iou_mdsig_list,*,iostat=istat) idc, mdsig_filename
         IF (istat .ne. 0) EXIT ! Loop
         
!         WRITE(*,*) 'imddc, idiag, mdsig_filename', 
!     &         imddc, idiag, mdsig_filename
         
         IF (idiag .gt. na_min) THEN
            CALL err_fatal(sub_name // ' Too many mdsig files',                &
     &         int=imddc)
         END IF
         IF (imddc .ne. idc) CALL err_warn(sub_name //                         &
     &      'imddc .ne. idc',int=idc,char=mdsig_filename)

         mdsig_filename_pp = TRIM(mdsig_list_path) //                          &
     &      ADJUSTL(mdsig_filename)
         CALL cdf_open(iou_mdsig, TRIM(mdsig_filename_pp), 'r',                &
     &       istat)
         CALL assert_eq(0,istat,sub_name // 'cdf open error')

!  Netcdf read of diagnostic_desc
!  Arrays of diagnostic_desc (d_desc) are in module v3f_global
         CALL diagnostic_cdf_read(d_desc(idiag),iou_mdsig)
         n_d_desc_mddc = imddc         ! these two should be taken out of
         n_d_desc = idiag              ! loop GJH 2010-01-21

         CALL cdf_close(iou_mdsig)

      END DO  !  loop over magnetic diagnostic signals

      WRITE(*,'(a,i4,a)') 'Read ',n_d_desc_mddc,                               &
     &   ' magnetic diagnostic signal (mdsig) files'

!...........Construct the signal corresponding to mddc diagnostic............................!
!.......NOTE: modified to constuct s_desc mddc signals *AFTER* the d_desc array has been.....!
!.......completely constructed (to make adding new signals more transparent).  JS 6/28/07....!

      n_s_desc = n_d_desc
      DO isig= 1, n_s_desc
        idiag = isig
        CALL signal_construct(s_desc(isig),'diagnostic',                       &
     &      d_desc(idiag) % s_name,d_desc(idiag) % l_name,                     &
     &      d_desc(idiag) % units,d_desc(idiag))
      END DO

      CALL getmem(mem)
      WRITE(*,*) 'Current Memory Usage (kB) is ',mem
      
!  Check that the r-z-phi grids are all the same. (Assert_eq doesn't check reals,
!  so don't bother with rmin, rmax, etc.)
!  Note assumption that signals start at 1.
      IF (n_s_desc .ge. 2) THEN
         
         CALL assert_eq(d_desc(1:n_s_desc) % mddc % mrf % ir,                  &
     &      sub_name // 'bad ir')
         CALL assert_eq(d_desc(1:n_s_desc) % mddc % mrf % jz,                  &
     &      sub_name // 'bad jz')
         CALL assert_eq(d_desc(1:n_s_desc) % mddc % mrf % kp,                  &
     &      sub_name // 'bad kp')
         CALL assert_eq(d_desc(1:n_s_desc) % mddc % mrf % kp_store,            &
     &      sub_name // 'bad kp_store')
         CALL assert_eq(d_desc(1:n_s_desc) % mddc % mrf %                      &
     &       n_field_periods,sub_name // 'bad n_field_periods')
     
         DO i = 2,n_s_desc
            CALL assert((d_desc(1) % mddc % mrf % lstell_sym .EQV.             &
     &         d_desc(i) % mddc % mrf % lstell_sym),                           &
     &         sub_name // 'bad lstell_sym')
         END DO
      END IF
!----------------------------------------------------------------------
!--  END mddc (Magnetic Diagnostics)                                 --
!----------------------------------------------------------------------

      RETURN
      END SUBROUTINE init_signal_mddc
!
!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------
      SUBROUTINE init_signal_ipsl
!--------------------------------------------------------------------------------
!
! Function:  Initialization routine for interferometery/polarimetry.  
!               (At the moment, all it does is fill the global ipcoll_g TYPE variable
!                using data from the input ip_*.nc file that was output by V3RFUN and then
!                converts in into an array of TYPE diagnostic_desc)
!
! INPUTS:   (none)
!
! OUTPUTS:  (none)
!
! global inputs: n_d_desc, n_s_desc
!
! global outputs: ipcoll_g, ipsl_cdf_infile, n_d_desc_g, n_s_desc_g, d_desc_g, s_desc_g
!
!  created June, 2007 by J. Shields
!  **updated June 28,2007 to fill existing global variables s_desc & d_desc
! 
! CAVEAT: Current version of routine assumes only two types of diagnostics (mddc & ipsl)
!         so it may need to be modified when more types are added.
! 
!----------------------------------------------------------------------------------
      USE stel_kinds
      USE stel_constants
      USE eq_T         ! module containing VMEC output info
      USE ip_beamline  ! module containing int_pol TYPE definition
      USE ip_global    ! module containing interfermeter/polarimeter global variables
      USE v3f_global, ONLY: intpol_filename, n_d_desc, n_s_desc,               &
     &                    n_d_desc_mddc, n_s_desc_mddc,                        &
     &                    n_d_desc_ipsl , n_s_desc_ipsl, d_desc, s_desc,       &
     &                    n_sdata_o, my_task
      USE faraday_mod, ONLY: read_intpol_data_from_cdf
      IMPLICIT NONE


!........local variables........................................................!

      INTEGER     :: counter = 0
      INTEGER     :: ilastslash
      INTEGER     :: i, j, k, idiag, isig
      INTEGER     :: idiag_ip_min            ! index of 1st ipsl element in d_desc
      INTEGER     :: isig_ip_min             ! index of 1st ipsl element in s_desc
      INTEGER     :: iflag, istat
      INTEGER     :: n_ip_units
      INTEGER     :: sum_nbeams

      CHARACTER(len=80) :: mychar1, filestr
      CHARACTER(len=80) :: mychar2, message
      CHARACTER(len=80) :: frcm_outfile
      CHARACTER(len=300) :: intpol_path
      CHARACTER(len=*), PARAMETER :: subname = 'init_signal_ipsl: '

      write(*,*) '===================================================='    
      write(*,*) ' '
      write(*,*) 'Hi.  Now in ', subname

!      write(*,*) subname, 'intpol file = ', TRIM(intpol_filename)

      write(*,*) subname, 'input n_d_desc = ', n_d_desc
      write(*,*) subname, 'input n_s_desc = ', n_s_desc

!..........**DEBUGGING: hard-wire global variable n_sdata_o = 0 to allow
!...........simple testing of task_ipsl_signals_test.  (obviously this hard-wire will
!...........need to be DELETED at some point...). JS 6/29/07
!      n_sdata_o = 0

      intpol_filename = ADJUSTL(intpol_filename)

!........read in int/polarimeter path info from input .nc file...............!
      call read_intpol_data_from_cdf(intpol_filename, ipcoll_g)


!........verify that the "intpol collection" variable ipcoll_g was read in correctly...!
      do i = 1,1

        write(*,*) '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
        write(*,*) subname, '**** READ-IN intpol element = ', i
        write(*,*) 'READ IN nbeams = ', ipcoll_g%ip(i)%nbeams   
        write(*,*) 'READ IN q0vec 1 = ', ipcoll_g%ip(i)%ipbeam(1)%q0vec   
        write(*,*) 'READ IN qfvec 1 = ', ipcoll_g%ip(i)%ipbeam(1)% qfvec   
        write(*,*) 'READ IN q0vec 2 = ', ipcoll_g%ip(i)%ipbeam(2)%q0vec   
        write(*,*) 'READ IN qfvec 2 = ', ipcoll_g%ip(i)%ipbeam(2)% qfvec   

        write(*,*) 'READ IN beam1 q0 = ', ipcoll_g%ip(i)%ipbeam(1)%q0   
        write(*,*) 'READ IN wavelength = ',                                    &
     &             ipcoll_g%ip(i)%ipbeam(1)% wavelength
!        write(*,*) 'READ IN B_ratio = ', ipcoll_g%ip(i)%ipbeam(1)% B_ratio
        write(*,*) 'READ IN s_name 1= ',                                       &
     &               ipcoll_g%ip(i)%ipbeam(1)% s_name
        write(*,*) 'READ IN l_name 1= ',                                       &
     &              ipcoll_g%ip(i)%ipbeam(1)% l_name
        write(*,*) 'READ IN s_name 2= ',                                       &
     &              ipcoll_g%ip(i)%ipbeam(2)% s_name
        write(*,*) 'READ IN l_name 2= ',                                       &
     &              ipcoll_g%ip(i)%ipbeam(2)% l_name
      end do

!..............allocate the global d_desc_g & s_desc_g arrays......................!
      n_ip_units = ipcoll_g%n_ip

      sum_nbeams = 0
      do i = 1, n_ip_units
        sum_nbeams = sum_nbeams + ipcoll_g % ip(i) % nbeams
      end do


      SELECT CASE(TRIM(my_task))
      
      CASE('reconstruct_a1')

!..........compute the total # of ipsl diagnostics..............!
        n_d_desc_ipsl = 2 * sum_nbeams
        n_s_desc_ipsl = n_d_desc_ipsl


!        n_d_desc_g = n_ipsl_diags
!        n_s_desc_g = n_ipsl_diags
!        ALLOCATE( d_desc_g(n_d_desc_g) )
!        ALLOCATE( s_desc_g(n_s_desc_g) )

        write(*,*) subname, 'n_ip_units = ', n_ip_units
        write(*,*) subname, 'n_d_desc_g = ', n_d_desc_g


!..........convert info into a (global) "diagnostic_desc" TYPE variable. (this entails treating each
!..........microwave beam as a completely separate diagnostic)...............!
        idiag_ip_min = n_d_desc + 1
        isig_ip_min  = n_s_desc + 1

!.........initialize diagnostic index to last filled mddc index...................!
!        idiag = 0 
        idiag = n_d_desc 


        do i = 1, n_ip_units
          do j = 1, ipcoll_g % ip(i) % nbeams
            do k = 1,2
              idiag = idiag + 1
              d_desc(idiag)%d_type = 'ipsl'
              d_desc(idiag)%ipsl%ip_sname = TRIM(ipcoll_g%ip(i)%s_name)
              d_desc(idiag)%ipsl%ip_lname = TRIM(ipcoll_g%ip(i)%l_name)
              d_desc(idiag)%ipsl%ipbeam   = ipcoll_g%ip(i)%ipbeam(j)

              SELECT CASE (k)
              CASE(1)
                d_desc(idiag)%ipsl%ipsl_type = 'interf'              
              CASE(2)
                d_desc(idiag)%ipsl%ipsl_type = 'polarim'              
              CASE DEFAULT
                write(*,*)  subname, 'error in CASE selection'
              END SELECT

!..............define the "long name" of the diagnostic to be a concatentation of the ip_type, .....!
!..............and the short names of the "ip unit" (ie the source) and the individual beam...........!
              d_desc(idiag)%l_name =                                           &
     &                          TRIM(d_desc(idiag)%ipsl%ipsl_type) //          &
     &                   '_' // TRIM(ipcoll_g%ip(i)%s_name) //                 &
     &                   '_' // TRIM(ipcoll_g%ip(i)%ipbeam(j)%s_name)
            end do
          end do
        end do


        write(*,*) subname, 'final idiag after loop = ', idiag

!...........modify n_d_desc to include ipsl diagnostics..................!
        n_d_desc = n_d_desc + n_d_desc_ipsl
        n_s_desc = n_s_desc + n_s_desc_ipsl



!..............now construct the signal description array.......................!
        do i= isig_ip_min, n_s_desc
           CALL signal_construct( s_desc(i),'diagnostic',                      &
     &        d_desc(i) % s_name,d_desc(i) % l_name,                           &
     &        d_desc(i) % units,d_desc(i) )
        end do


!        n_d_desc = 2
!        n_s_desc = 2

        write(*,*) subname, 'output n_d_desc = ', n_d_desc
        write(*,*) subname, 'output n_s_desc = ', n_s_desc

      CASE DEFAULT
         WRITE(*,*) subname, 'CASE SELECT - DEFAULT'
         WRITE(*,*) 'my_task = ', my_task
      END SELECT

      write(*,*) 'Now exiting ', subname


      RETURN
      END SUBROUTINE init_signal_ipsl
!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------
!
      SUBROUTINE init_signal_geometric
!--------------------------------------------------------------------------
!  Subroutine to initialize the geometric signals
!  First Version - 2009-02-01
!  2009-04-13 - Change to new parameterization of limiter_iso
!  SUBJECT TO CHANGE
!---------------------------------------------------------------------------
!      USE ezcdf
!      USE bsc_cdf
      USE signal_T
      USE geometric_T
      USE edge_limit_T
      USE limiter_iso_T
      USE safe_open_mod
      USE v3_utilities
      USE v3f_global
      
      IMPLICIT NONE

!----------------------------------------------------------------------
! L O C A L Variable Declarations
!----------------------------------------------------------------------
      TYPE (geometric_desc)     :: geom_local
      TYPE (edge_limit_desc)    :: el_desc_local
      TYPE (limiter_iso)        :: lim_iso_local

      INTEGER                            :: numin_loc
      REAL(rprec), DIMENSION(na_phi_lif) :: vgrid_loc     

      CHARACTER(len=20) :: el_name, el_units, el_type, g_name, g_units
      REAL(rprec)       :: el_sigma_default, g_sigma_default
      LOGICAL           :: el_l_on_edge

      INTEGER :: isig, igeom, i, j, n_phi_loc
      
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'init_signal_geometric: '

!----------------------------------------------------------------------
!      Variables from v3f_global, used in this subroutine:
!  geometric_filename
!  n_s_desc
!  r_major_radius
!  a_minor_radius
!  n_g_desc
!  g_desc
!  lif variables
!----------------------------------------------------------------------

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
      WRITE(*,*) ' *** IN SUBROUTINE ', sub_name

!  Assumption - n_s_desc (from v3f_global) is the number of signal
!  descriptions so far. 
!  In the subroutine we increment this counter as we add geometric signals

!  First, check to see if wish to use NLI variables r_major_radius and
!  a_minor_radius to define a circular torus.
!  NB. Coded to OVERWRITE the first lif variables.
!  NB. 2009-04-13 - May decide to eliminate this option. It is hardwired
!        for CTH.

      IF ((r_major_radius .gt. zero) .AND.                                     &
     &    (a_minor_radius .gt. zero)) THEN
         n_lif = MAX(1,n_lif)
         lif_arz(1,:,:) = zero
         lif_arz(1,0,0) = - a_minor_radius ** 2
         lif_arz(1,2,0) = one
         lif_arz(1,0,2) = one
         lif_rc(1) = r_major_radius
         lif_zc(1) = zero
         lif_sigma(1) = 0.001_rprec
         n_phi_lif(1) = 4
         lif_phi_degree(1,1) = zero 
         lif_phi_degree(1,2) = 18. 
         lif_phi_degree(1,3) = 36.  
         lif_phi_degree(1,4) = 54.
         lif_on_edge(1) = .true.
      ENDIF

      numin_loc = 20     ! Plucked from thin air - minimum number of points on s=1
!  Loop over the lif_ variables
      DO i = 1,n_lif

! vgrid_loc - phi values (radians) at which to evaluate limiter_iso function
         CALL assert((n_phi_lif(i) .LE. na_phi_lif),sub_name //                &
     &     'n_phi_lif too big','WARNING')
         CALL assert((n_phi_lif(i) .GE. 1),sub_name //                         &
     &     'n_phi_lif too small','WARNING')
         n_phi_loc = MAX(MIN(na_phi_lif,n_phi_lif(i)),1)
         DO j = 1,n_phi_loc
            vgrid_loc(j) = lif_phi_degree(i,j) * twopi / 360.
         END DO
      
! Construct a local limiter_iso
!      SUBROUTINE limiter_iso_construct(this,arz,rc,zc,numin,vgrid)
         CALL limiter_iso_construct(lim_iso_local,lif_arz(i,:,:),              &
     &      lif_rc(i),lif_zc(i),numin_loc,                                     &
     &      vgrid_loc(1:n_phi_loc))
         
!  construct the edge_limit
!      SUBROUTINE edge_limit_desc_construct(this,name,edge_limit_type,          &
!     &   units,sigma_default,l_on_edge,lim_iso)
         WRITE(el_name,'("edge_lim_",i1)') i
!         el_name = 'in isg'
         el_units = 'meter'
         el_sigma_default = lif_sigma(i)
         el_type = 'iso_fun'
         el_l_on_edge = lif_on_edge(i)
         CALL edge_limit_construct(el_desc_local,el_name,el_type,              &
     &      el_units,el_sigma_default,el_l_on_edge,lim_iso_local)

!  construct the geometric
!      SUBROUTINE geometric_desc_construct(this,g_type,name,                    &
!     &   units,sigma_default,el_desc)
         g_name = el_name
         g_units = el_units
         g_sigma_default = el_sigma_default
         CALL geometric_construct(g_desc(i),'edge_limit',g_name,               &
     &      g_units,g_sigma_default,el_desc_local)
      END DO
      n_g_desc = n_lif

!  Now all the geometrics have been constructed. Construct the corresponding
!  signal descriptions
      DO isig = n_s_desc + 1,n_s_desc + n_g_desc
         igeom = isig - n_s_desc
        CALL signal_construct(s_desc(isig),'geometric',                        &
     &      g_desc(igeom) % name,g_desc(igeom) % name,                         &
     &      g_desc(igeom) % units,g_desc(igeom))
      END DO

!  increment the number of signals
      n_s_desc = n_s_desc + n_g_desc

      RETURN
      END SUBROUTINE init_signal_geometric

!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------
!
      SUBROUTINE init_signal_coosig
!--------------------------------------------------------------------------
!  Subroutine to initialize the coosig signals
!  First Version - 2011-07-08
!  SUBJECT TO CHANGE
!---------------------------------------------------------------------------
      USE signal_T
      USE coosig_T
      USE safe_open_mod
      USE v3_utilities
      USE v3f_global
      
      IMPLICIT NONE

!----------------------------------------------------------------------
! L O C A L Variable Declarations
!----------------------------------------------------------------------
      TYPE (coosig_desc)        :: cd_local

      INTEGER                            :: numin_loc
      REAL(rprec), DIMENSION(na_phi_lif) :: vgrid_loc     

!      CHARACTER(len=20) :: el_name, el_units, el_type, g_name, g_units
!      REAL(rprec)       :: el_sigma_default, g_sigma_default
!      LOGICAL           :: el_l_on_edge

!      INTEGER :: isig, igeom, i, j, n_phi_loc

      INTEGER :: i, j, ic, n_terms, isig
      
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'init_signal_coosig: '


!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
      WRITE(*,*) ' *** IN SUBROUTINE ', sub_name

!  Assumption - n_s_desc (from v3f_global) is the number of signal
!  descriptions so far. 
!  In the subroutine we increment this counter as we add coosig signals

!  Loop over the coosig_ variables
      ic = 0  ! counter for coosigs
      DO i = 1,n_coosig

!  Find out how many terms there are in this coosig
         n_terms = 0
         DO j = 1,na_coef
            IF (coosig_indices(i,j) .gt. 0) THEN
               n_terms = n_terms + 1
            ELSE
               EXIT ! j-index loop
            ENDIF
         END DO
         
         IF (n_terms .le. 0) THEN
            CALL err_warn(sub_name // 'no terms in coosig ',int=i)
            CYCLE ! i-index loop
         ENDIF
         
!  Define the coosig_desc. Note the different indices ic and i (NLI)
!  Note that the c_desc array, declared in v3f_global, is allocated
!  in init_main_nli.
      
         ic = ic + 1
         CALL coosig_construct(c_desc(ic),coosig_type(i),                      &
     &      coosig_name(i),coosig_units(i),zero,
     &      coosig_indices(i,1:n_terms),coosig_coeff(i,1:n_terms))
      END DO
      n_s_desc_coosig = ic
               
!  Now all the coosigs have been constructed. Construct the corresponding
!  signal descriptions
      DO isig = n_s_desc + 1,n_s_desc + n_s_desc_coosig
         ic = isig - n_s_desc
        CALL signal_construct(s_desc(isig),'coosig',                           &
     &      c_desc(ic) % name,c_desc(ic) % name,                               &
     &      c_desc(ic) % units,c_desc(ic))
      END DO

!  increment the number of signals
      n_s_desc = n_s_desc + n_s_desc_coosig

      RETURN
      END SUBROUTINE init_signal_coosig

!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------
!
      SUBROUTINE init_equilibrium
      USE v3f_global
      USE eq_interface

      IMPLICIT NONE

!----------------------------------------------------------------------
! L O C A L Variable Declarations
!----------------------------------------------------------------------
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'init_equilibrium: '
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
      CALL eq_init_file(vmec_nli_filename,model_a % eqstate)
      
      WRITE(*,*) ' *** EXITING ', sub_name

      RETURN
      END SUBROUTINE init_equilibrium
!! 
!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------

!*******************************************************************************
! SECTION III.   TASK SUBROUTINES
!*******************************************************************************
!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------
!
      SUBROUTINE task_reconstruct_a1
!  Task subroutine, Reconstruction, Algorithm 1
! 
!  2011-03-02 JDH
!    Added choice of step type - straight-line or segmented
!    Added coding to retry a step if it converged, but had higher g^2.
!
!  2011-02-28 JDH
!    Move reconstruction step control variables to v3f_global
!
! 2009-07-10 JDH
!    Add loop over step trys - to make code more robust. Also moved
!    recon_param print to subroutine in recon_param_T
!
! 2008-06-04  JDH 
!    Added more sophisticated method of choosing number of SV
!    to use. Necessitated revising r_work assumptions
!  10-07-06. First version of task subroutine
!  12-28-2006. Significant revisions to print out.
!  2007-06-22. Rearrange variables a bit, revise printout at end,
!    Add simple convergence criterion

      USE v3f_global, ONLY: nrstep, dg2_stop, cut_svd, cut_eff,                &
     &   cut_marg_eff, cut_delta_a, cut_dg2, astep_max,                        &
     &   step_type, model_a, i_work, r_work, s_desc
      USE signal_mc
      USE recon_param_T
      USE recon_param_model
      USE gsq_mod
      USE mmaw
      IMPLICIT NONE
!----------------------------------------------------------------------
! A R G U M E N T Declarations
!----------------------------------------------------------------------

!----------------------------------------------------------------------
! L O C A L Variable Declarations
!----------------------------------------------------------------------

!      Derived Types
!  rp_old       type Recon_Param, vector, old values
!  rp_new       type Recon_Param, vector, new values
!  state_vec    type eq_state, pointer to component of model_a
!  sd_model     type signal_data, array, model values

!      Save variables, for use in printout at end of subroutine
!  save_g2      vector to save gsq values
!  save_rpvra   two-d array of reals, save reconstruction parameter value
!               step history (Recon_Param % Value Real Array)

!      Other local variables
!  rpvra_da     vector of reals, delta of normalized parameters
!  rprva_dr     vector of reals, delta of unnormalized parameters
!  temp_npvec   vector of length gsq_np, for parameter error printout
!  nrs          integer, number of actual reconstruction steps taken
!  irstep       integer, index for reconstruction step loop

      REAL(rprec), DIMENSION(5) :: step_cntrl_array

      TYPE(recon_param), DIMENSION(:), ALLOCATABLE :: rp_old, rp_new
      TYPE(eq_state), POINTER :: state_vec => null()
      TYPE(eq_state) :: eq_state_0
      TYPE(signal_data), DIMENSION(:), ALLOCATABLE :: sd_model

      REAL(rprec), DIMENSION(:), ALLOCATABLE    :: save_g2, save_g2_exp
      REAL(rprec), DIMENSION(:), ALLOCATABLE    :: save_da
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE  :: save_rpvra
      INTEGER, DIMENSION(:), ALLOCATABLE :: save_nsv

      INTEGER       :: irstep, irp, i, kuse, nrs
      INTEGER       :: numiter, isig, ier1, ier2, ier3
      INTEGER       :: iter_v

      INTEGER       :: i_try, n_try = 5
      REAL(rprec)  :: try_factor
      REAL(rprec), DIMENSION(5) :: try_factor_a =                              &
     &    (/ 1., .5, .25, .0625, .015625 /)

      REAL(rprec), DIMENSION(:), ALLOCATABLE :: rpvra_da, rpvra_dr,            &
     &   temp_npvec
      REAL(rprec) :: sum_da, sum_dr, this_da, this_dr, net_dr, net_da,         &
     &   g2_diff, fracstep_itry, astep_try, astep_i_try_1, astep_use,          &
     &   g2_exp

      LOGICAL              :: lconv, l_step_succeed
      CHARACTER(len=25)     :: a_string
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'task_reconstruct_a1: '
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
      WRITE(*,*) ' *** IN SUBROUTINE ', sub_name
      WRITE(iou_runlog,*) ' *** IN SUBROUTINE ', sub_name

      IF (cut_svd .eq. zero) THEN
!  Old - deprecated code. Still using i_work and r_work
!  Eliminated even older code, only svd_cut and dg2_step (r_work(1:2)) 2011-0228
         nrstep = i_work(1)
         dg2_stop = r_work(1)
         cut_svd = r_work(4)
         cut_eff = r_work(5)
         cut_marg_eff = r_work(6)
         cut_delta_a = r_work(7)
         cut_dg2 = r_work(8)
         astep_max = r_work(9)
         step_cntrl_array = r_work(4:8)
      ELSE
!  New coding. Use reconstruction step control variables - NLI
         step_cntrl_array(1) = cut_svd
         step_cntrl_array(2) = cut_eff
         step_cntrl_array(3) = cut_marg_eff
         step_cntrl_array(4) = cut_delta_a
         step_cntrl_array(5) = cut_dg2
      ENDIF

!  point to equilibrium state
      state_vec => model_a % eqstate

!  Allocate space for signal data stuff.
!  Allocate sd_model to be number of observed signals
!  Note that sd_observe is allocated larger.

      CALL assert_eq(n_sdata_o,n_s_desc,sub_name //                            &
     &       ' problem with n_sdata_o, n_s_desc ') 
      ALLOCATE(sd_model(n_s_desc),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc sd_model')

!  Create the sd_model - makes sure that the signal_desc is properly pointed to
      DO isig = 1,n_s_desc
         CALL signal_construct(sd_model(isig),s_desc(isig),'model')
      END DO

!  Iterate VMEC to equilibrium      
      numiter = state_vec % fixp % niter
      CALL eq_history_set(i1 = 0)   !  Reconstruction step zero
      CALL eq_step(numiter,state_vec,lconverged = lconv,                       &
     &   iter_es = iter_v)
      WRITE(*,*) ' Initial VMEC equilibration'
      WRITE(iou_runlog,*) ' Initial VMEC equilibration'
      WRITE(*,*)  '   lconverged = ', lconv, ' iter_ha = ', iter_v
      WRITE(iou_runlog,*)  '   lconverged = ', lconv,                          &
     &   ' iter_ha = ', iter_v
      CALL assert(lconv,sub_name // ' no initial convergence')
      eq_state_0 = state_vec

! Compute model signals
      CALL signal_mc_model_compute(sd_model,model_a)
 
!  Allocate space for reconstruction parameters and arrays
      ALLOCATE(rp_old(n_rparam),rp_new(n_rparam),STAT=ier1)
      ALLOCATE(save_rpvra(0:nrstep,n_rparam),STAT=ier2)
      ALLOCATE(rpvra_da(n_rparam),rpvra_dr(n_rparam),STAT=ier3)
      CALL assert_eq(0,ier1,ier2,ier3,sub_name // 'alloc rp_')
      ALLOCATE(save_g2(0:nrstep),STAT=ier1)
      ALLOCATE(temp_npvec(gsq_nrp),STAT=ier2)
      ALLOCATE(save_g2_exp(1:nrstep),STAT=ier3)
      CALL assert_eq(0,ier1,ier2,ier3,sub_name // 'alloc save_g2')
      ALLOCATE(save_nsv(1:nrstep),STAT=ier1)
      ALLOCATE(save_da(1:nrstep),STAT=ier2)
      CALL assert_eq(0,ier1,ier2,sub_name // 'alloc save_nsv')

!  Warnings for Fixed-Boundary Reconstructions
!   (Don't warn for edge_limit signal - may want just for information)
!   (VMEC doesn't use extcur for fixed boundary)
      IF ( .NOT. state_vec % fixp % lfreeb) THEN
         DO i = 1,n_rparam
            IF (rparam(i) % p_type .EQ. 'extcur')  THEN
               CALL err_warn(sub_name // 'Fixed-Boundary ' //                  &
     &           ' recon_param Problem',rparam(i) % p_type,int=i,              &
     &            log=state_vec % fixp % lfreeb)
            ENDIF
         END DO
      ENDIF
      
!  Get the initial reconstruction parameters, g2
      rp_old = rparam(1:n_rparam)
      CALL recon_param_model_get(rp_old,model_a)
      save_rpvra(0,:) = rp_old(:) % value
      CALL gsq_evaluate_g2(sd_model,sd_observe(1:n_s_desc))
      save_g2(0) = gsq_g2

!  Initial Print Out
      WRITE(*,*)
      WRITE(*,*) ' ** RECONSTRUCTION STEP 0'
      WRITE(*,*) '  g2 = ', save_g2(0)
      WRITE(iou_runlog,*)
      WRITE(iou_runlog,*) ' ** RECONSTRUCTION STEP 0'
      WRITE(iou_runlog,*) '  g2 = ', save_g2(0)
      CALL recon_param_write(rp_old,                                           &
     &   ' Initial Reconstruction Parameters',iou_runlog)
      
!  Loop over the reconstruction steps
      sum_dr = 0.
      sum_da = 0.
      DO irstep = 1,nrstep
         CALL eq_history_set(i1 = irstep)   !  Reconstruction step #
         WRITE(*,*)
         WRITE(*,*) ' RECONSTRUCTION STEP ', irstep 
         WRITE(iou_runlog,*)
         WRITE(iou_runlog,*) ' RECONSTRUCTION STEP ', irstep 

!  Evaluate and write out the Jacobian
         CALL gsq_evaluate_jac(model_a,sd_model,sd_observe(1:n_s_desc),        &
     &      rp_old,rcnstrnts)
         WRITE(a_string,'(a,i2)') ' In Loop, irstep = ',irstep
         CALL gsq_write_jac(a_string,iou_runlog,verbose = 1)

!  Save the gsq_evector values, so that can compute expected changes
!  to gsq (in gsq_evaluate_delta_a_x)
         CALL gsq_save_evector
     
!  Compute # SV to use for step. Note that gsq_evaluate_jac has called gsq_evaluate_g2.
         CALL gsq_evaluate_kuse(step_cntrl_array,kuse)

! Loop over attempts to take a step
         astep_i_try_1 = astep_max
         DO i_try = 1,n_try
            l_step_succeed = .false.
            try_factor = try_factor_a(i_try)
            astep_try = astep_i_try_1 * try_factor
            
            SELECT CASE (TRIM(step_type))
               CASE ('seg')
                  CALL gsq_evaluate_delta_a_seg(astep_try,kuse,                &
     &               rpvra_da,rpvra_dr,astep_use,g2_exp)
               CASE ('lm')  ! Levenberg-Marquardt - not yet implemented
                  CALL gsq_evaluate_delta_a_lm(astep_try,kuse,                 &
     &               rpvra_da,rpvra_dr,astep_use,g2_exp)
               CASE DEFAULT  !  straight line
                  CALL gsq_evaluate_delta_a_sl(astep_try,kuse,                 &
     &               rpvra_da,rpvra_dr,astep_use,g2_exp)
            END SELECT 

            IF (i_try .eq. 1) astep_i_try_1 = astep_use
            rp_new = rp_old
            CALL recon_param_change_value(rp_new,delta_value=rpvra_dr)
            save_rpvra(irstep,:) = rp_new(:) % value
            save_nsv(irstep) = kuse
            save_g2_exp(irstep) = g2_exp
            WRITE(*,*) '    i_try ', i_try, ' Expected g2 = ',                 &
     &         save_g2_exp(irstep)
            WRITE(iou_runlog,*) '    i_try ', i_try,                           &
     &         ' Expected g2 = ', save_g2_exp(irstep)
         
!  Print out information about parameters
            CALL recon_param_write_ond(rp_old,rp_new,unit=iou_runlog)

!  Restore eq_state_0 xc vector to the VMEC Internal State
!  JDH 2011-06-21 Change eq_change_vmi_cp_xc to eq_change_vmi_reset_xc
!  JDH 2011-06-24 - Back to eq_change_vmi_cp_xc
            CALL eq_change_vmi_cp_xc(eq_state_0)
!            CALL eq_change_vmi_reset_x(eq_state_0)
            
!  Change the parameters
            CALL recon_param_model_put(rp_new,model_a,rcnstrnts,'VMI')

!  Turn off VMEC 2d preconditioning
            CALL eq_change_vmi_precon2d_off_new_blocks

!  Iterate VMEC to equilibrium at new parameter values     
            CALL eq_history_set(i2 = - i_try)   !  Try #
            numiter = state_vec % fixp % niter
            CALL eq_step(numiter,state_vec,lconverged = lconv,                 &
     &         iter_es = iter_v)
     
            IF (.not. lconv) THEN
                  WRITE(*,*) 'VMEC did not converge.'
                  WRITE(iou_runlog,*) 'VMEC did not converge.'
               IF (i_try .LT. n_try) THEN
                  WRITE(*,*) 'Changing step size factor to ',                  &
     &               try_factor_a(i_try + 1)
                  WRITE(iou_runlog,*)                                          &
     &               'Changing step size factor to ',                          &
     &                try_factor_a(i_try + 1)
               ELSE
                  WRITE(*,*) 'This was last try '
                  WRITE(iou_runlog,*) 'This was last try '                  
               ENDIF
               WRITE(*,*)
               WRITE(iou_runlog,*)
               CYCLE ! the loop over step size trys
            ENDIF

!  Recompute model signals and gsq
            CALL signal_mc_model_compute(sd_model,model_a)
            CALL gsq_evaluate_g2(sd_model,sd_observe(1:n_s_desc))

!  The step converged. But was it a useful step? Check it.
!    EXIT loop over step size trys if a good step, otherwise decrease
!    step size to try
            IF (gsq_g2 .le. save_g2(irstep - 1)) THEN
               l_step_succeed = .true.
               EXIT ! loop over step size trys
            ELSE
               WRITE(*,*) 'VMEC Converged, but g2 increased: ',                &
     &            gsq_g2
               WRITE(iou_runlog,*) 'VMEC Converged, but g2',                   &
     &            'increased: ',  gsq_g2
               IF (i_try .LT. n_try) THEN
                  WRITE(*,*) 'Changing step size factor to ',                  &
     &               try_factor_a(i_try + 1)
                  WRITE(iou_runlog,*)                                          &
     &               'Changing step size factor to ',                          &
     &                try_factor_a(i_try + 1)
               ELSE
                  WRITE(*,*) 'This was last try '
                  WRITE(iou_runlog,*) 'This was last try '                                    
               ENDIF
               WRITE(*,*)
               WRITE(iou_runlog,*)
            ENDIF           
               
        END DO ! loop over step size trys

!  Unsuccessful step
      IF (.not. l_step_succeed) THEN
         WRITE(*,*)
         WRITE(*,*) ' Result of Step - DID NOT SUCCEED '
         WRITE(*,*) '  g2 old      '
         WRITE(*,'(3(4x,es11.4))') save_g2(irstep-1)
         WRITE(iou_runlog,*)
         WRITE(iou_runlog,*) ' Result of Step - DID NOT SUCCEED '
         WRITE(iou_runlog,*) '  g2 old      '
         WRITE(iou_runlog,'(3(4x,es11.4))') save_g2(irstep-1)
         EXIT  ! Loop over reconstruction steps
      ENDIF

!  Save information about step
         eq_state_0 = state_vec

!  JDH 2011-10-20. Fix this_da and this_dr computation. No need for try_factor.
!         this_da = try_factor *                                                &
!     &      SQRT(DOT_PRODUCT(rpvra_da(:),rpvra_da(:)))
         this_da = SQRT(DOT_PRODUCT(rpvra_da(:),rpvra_da(:)))
         save_da(irstep) = this_da
!         this_dr = try_factor *                                                & 
!     &      SQRT(DOT_PRODUCT(rpvra_dr(:),rpvra_dr(:)))
         this_dr = SQRT(DOT_PRODUCT(rpvra_dr(:),rpvra_dr(:)))

         sum_da = sum_da + this_da
         sum_dr = sum_dr + this_dr
!         WRITE(*,'(2(a,es11.4))') '  length da = ', this_da,                  &
!     &      ' length dr = ', this_dr
         WRITE(iou_runlog,*) '  length da = ', this_da,                        &
     &      ' length dr = ', this_dr

         save_g2(irstep) = gsq_g2

!  Print Out Step Result
         WRITE(*,*)
         WRITE(*,*) ' Result of Step '
         WRITE(*,*) '  g2 old      g2 new     g2 diff'
         g2_diff = save_g2(irstep) - save_g2(irstep-1)
         WRITE(*,'(3(4x,es11.4))') save_g2(irstep-1), save_g2(irstep),         &
     &      g2_diff
         WRITE(iou_runlog,*)
         WRITE(iou_runlog,*) ' Result of Step '
         WRITE(iou_runlog,*) '  g2 old      g2 new     g2 diff'
         WRITE(iou_runlog,'(3(4x,es11.4))') save_g2(irstep-1),                 &
     &      save_g2(irstep), g2_diff

        rp_old = rp_new

!  This is where an exit test should go
!  Convergence Tests
         IF ((g2_diff .lt. zero) .and. (g2_diff .gt. - dg2_stop)) THEN
            WRITE(*,*) 
            WRITE(*,*) 'dg2_stop = ', dg2_stop
            WRITE(*,*) 'g2 decreasing slowly - reconstructed'
            WRITE(iou_runlog,*) 
            WRITE(iou_runlog,*) 'dg2_stop = ', dg2_stop
            WRITE(iou_runlog,*)                                                &
     &            'g2 decreasing slowly - reconstructed'
            EXIT
         ENDIF
         
         IF (save_g2(irstep) .lt. dg2_stop) THEN
            WRITE(*,*) 
            WRITE(*,*) 'dg2_stop = ', dg2_stop
            WRITE(*,*) 'g2 small - reconstructed'
            WRITE(iou_runlog,*) 
            WRITE(iou_runlog,*) 'dg2_stop = ', dg2_stop
            WRITE(iou_runlog,*) 'g2 small - reconstructed'
            EXIT
         ENDIF

         IF (irstep .ge. nrstep) THEN
            WRITE(*,*) 
            WRITE(*,*) 'dg2_stop = ', dg2_stop
            WRITE(*,*) 'nrstep = ', nrstep
            WRITE(*,*) irstep,' steps completed - may be reconstructed'
            WRITE(iou_runlog,*) 
            WRITE(iou_runlog,*) 'dg2_stop = ', dg2_stop
            WRITE(iou_runlog,*) 'nrstep = ', nrstep
            WRITE(iou_runlog,*) irstep,                                        &
     &          ' steps completed - may be reconstructed'
         ENDIF

      END DO

!  Print out final information
      nrs = MIN(irstep,nrstep)
      IF (.not. l_step_succeed) nrs = nrs - 1
      WRITE(*,*)
      WRITE(*,*) ' RECONSTRUCTION FINISHED '
      net_dr = SQRT(DOT_PRODUCT(save_rpvra(nrs,:) - save_rpvra(0,:),           &
     &   save_rpvra(nrs,:) - save_rpvra(0,:)))
      net_da = SQRT(DOT_PRODUCT(                                               &
     &   (save_rpvra(nrs,:) - save_rpvra(0,:)) / gsq_ppi(:),                   &
     &   (save_rpvra(nrs,:) - save_rpvra(0,:)) / gsq_ppi(:)))
      WRITE(*,*) '  net delta parameters (unnormalized) = ', net_dr
      WRITE(*,*) '  sum delta parameters (unnormalized) = ', sum_dr
      WRITE(*,*) '  net delta parameters (normalized) = ', net_da
      WRITE(*,*) '  sum delta parameters (normalized) = ', sum_da
      WRITE(*,*) '  g2: start, end, diff  = ', save_g2(0),                     &
     &   save_g2(nrs), save_g2(nrs) - save_g2(0)

      WRITE(iou_runlog,*)
      WRITE(iou_runlog,*) ' RECONSTRUCTION FINISHED '
      WRITE(iou_runlog,*) '  net delta parameters (unnorm) = ', net_dr
      WRITE(iou_runlog,*) '  sum delta parameters (unnorm) = ', sum_dr
      WRITE(iou_runlog,*) '  net delta parameters (normal) = ', net_da
      WRITE(iou_runlog,*) '  sum delta parameters (normal) = ', sum_da
      WRITE(iou_runlog,*) '  g2: start, end, diff  = ', save_g2(0),            &
     &   save_g2(nrs), save_g2(nrs) - save_g2(0)

      WRITE(iou_recout,*)
      WRITE(iou_recout,*) ' RECONSTRUCTION FINISHED '
      WRITE(iou_recout,*) '  net delta parameters (unnorm) = ', net_dr
      WRITE(iou_recout,*) '  sum delta parameters (unnorm) = ', sum_dr
      WRITE(iou_recout,*) '  net delta parameters (normal) = ', net_da
      WRITE(iou_recout,*) '  sum delta parameters (normal) = ', sum_da
      WRITE(iou_recout,*) '  g2: start, end, diff  = ', save_g2(0),            &
     &   save_g2(nrs), save_g2(nrs) - save_g2(0)

      WRITE(*,2000)
      WRITE(*,2100) 0, save_g2(0)

      WRITE(iou_runlog,2000)
      WRITE(iou_runlog,2100) 0, save_g2(0)

      WRITE(iou_recout,2000)
      WRITE(iou_recout,2100) 0, save_g2(0)
      DO i = 1,nrs
         WRITE(*,2100) i, save_g2(i), save_g2_exp(i), save_nsv(i),             &
     &      save_da(i)
         WRITE(iou_runlog,2100) i, save_g2(i), save_g2_exp(i),                 &
     &      save_nsv(i), save_da(i)
         WRITE(iou_recout,2100) i, save_g2(i), save_g2_exp(i),                 &
     &      save_nsv(i), save_da(i)
      END DO      
2000  FORMAT(//' Reconstruction Steps'/                                        &
     &   'irstep',t10,'g^2',t24,'g^2',t35,'# SV',t43,                          &
     &   'step size'/                                                          &
     &   t24,'expected',t35,'used',t43,'(normalized)')
2100  FORMAT(2x,i3,2x,es11.4,2x,es11.4,3x,i3,3x,es11.4)

!  Write out the confidence limit information
!      CALL gsq_write_conf(rp_new,'end reconstruct_a1',6)
      CALL gsq_write_conf(rp_new,'end reconstruct_a1',iou_recout)

!  Write out the most redundant parameters
!  Commented out JDH 2010-06-10. Not very useful.
!      CALL gsq_write_mr(rp_new,'end reconstruct_a1',6)
!      CALL gsq_write_mr(rp_new,'end reconstruct_a1',iou_recout)

!  Write out the model and observed signals
      CALL signal_data_moa_write(sd_model,sd_observe(1:n_s_desc),
     &   'End of reconstruct_a1',iou_recout)

!  Write out the model profiles
      CALL model_write_profiles(model_a,' Finished',iou_recout)

!  write out the vmec_history
      CALL eq_history_print

!  If the argument is present - return the irstep variable
!      IF (PRESENT(irstep_arg)) THEN
!         irstep_arg = irstep
!      ENDIF

      RETURN
      END SUBROUTINE task_reconstruct_a1
!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------
!
      SUBROUTINE task_gsq_jac_test
!  Task subroutine, to test the gsq jacobian calculation

!  10-03-06. First version of task subroutine

      USE v3f_global
      USE signal_mc
      USE recon_param_T
      USE recon_param_model
      USE gsq_mod
      USE mmaw
      
      IMPLICIT NONE

!----------------------------------------------------------------------
! L O C A L Variable Declarations
!----------------------------------------------------------------------
      TYPE(eq_state), POINTER :: state_vec => null()
      TYPE(signal_data), DIMENSION(:), ALLOCATABLE :: sd_model
      INTEGER       :: numiter, isig, ier1
      INTEGER       :: iter_v
      LOGICAL              :: lconv, lsuccess
      CHARACTER(len=*), PARAMETER :: fmt1 =                                    &
     &  '(i4,1x,i3,3(es20.13,2x))'
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'task_gsq_jac_test: '
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
!  point to equilibrium state
      state_vec => model_a % eqstate

!  Allocate space for signal data stuff.
!  Allocate sd_model to be number of observed signals
!  Note that sd_observe is allocated larger.

      ALLOCATE(sd_model(n_s_desc),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc')

!  Iterate VMEC to equilibrium      
!      numiter = state_vec % fixp % niter
!      CALL eq_step(numiter,state_vec,lconverged = lconv,                       &
!     &   iter_es = iter_v)
!      WRITE(*,*)  ' gsq_jac_test:   lconverged = ', lconv,                     &
!     &   ' iter_ha = ', iter_v
!      DO isig = 1,n_s_desc
!         CALL signal_mc_model_compute(sd_model(isig),s_desc(isig),             &
!     &      model_a)
!      END DO
            
!  Evaluate Jacobian
      CALL gsq_evaluate_jac(model_a,sd_model,sd_observe(1:n_s_desc),           &
     &   rparam(1:n_rparam),rcnstrnts)

!  Write out jacobian and SVD information. Arrays are from gsq_mod
      mmaw_iou = 50
      CALL mmaw_wreal('gsqjacobianobserve', gsq_jacobian_observe)
      CALL mmaw_wreal('gsqjacobianmodel', gsq_jacobian_model)
      CALL mmaw_wreal('gsqjacobian', gsq_jacobian)
      CALL mmaw_wreal('gsqjacobiannorm', gsq_jacobian_norm)
      CALL mmaw_wreal('gsqjsvdu', gsq_jsvd_u)
      CALL mmaw_wreal('gsqjsvdvt', gsq_jsvd_vt)
      CALL mmaw_wreal('gsqjsvdw', gsq_jsvd_w)

!      WRITE(*,*)
!      WRITE(*,*) 'gsq_jacobian_observe', gsq_jacobian_observe
!      WRITE(*,*)
!      WRITE(*,*) 'gsq_jacobian_model', gsq_jacobian_model
!      WRITE(*,*)
!      WRITE(*,*) 'gsq_jacobian', gsq_jacobian
!      WRITE(*,*)
!      WRITE(*,*) 'gsq_jacobian_norm', gsq_jacobian_norm
!      WRITE(*,*)
!      WRITE(*,*) 'gsq_jsvd_u', gsq_jsvd_u
!      WRITE(*,*)
!      WRITE(*,*) 'gsq_jsvd_vt', gsq_jsvd_vt
!      WRITE(*,*)
!      WRITE(*,*) 'gsq_jsvd_w', gsq_jsvd_w
!      WRITE(*,*)

!  write out the vmec_history
      CALL eq_history_print

      RETURN
      END SUBROUTINE task_gsq_jac_test
!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------
!
      SUBROUTINE task_gsq_on_grid
!  Task subroutine, to evaluate the gsq function on a (2d) grid of parameter
!   values 

!  09-20-06. First version of task subroutine

      USE v3f_global
      USE signal_mc
      USE recon_param_T
      USE recon_param_model
      USE gsq_mod
      USE mmaw
      
      IMPLICIT NONE

!----------------------------------------------------------------------
! L O C A L Variable Declarations
!----------------------------------------------------------------------
      TYPE(eq_state), POINTER :: state_vec => null()
      TYPE(signal_data), DIMENSION(:), ALLOCATABLE :: sd_model
      INTEGER       :: n1, n2, i1, i2, numiter, isig, ier1
      REAL(rprec)          :: p1min, p1max, p2min, p2max, p1, p2
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE  ::  gsq
      INTEGER       :: iter_v
      LOGICAL              :: lconv, lsuccess
      CHARACTER(len=30)    :: p1_str, p2_str, gsq_str
      CHARACTER(len=*), PARAMETER :: fmt1 =                                    &
     &  '(i4,1x,i3,3(es20.13,2x))'
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'task_gsq_on_grid: '
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
!  point to equilibrium state
      state_vec => model_a % eqstate

!  Allocate space for signal data stuff.
!  Allocate sd_model to be number of observed signals
!  Note that sd_observe is allocated larger.

      ALLOCATE(sd_model(n_s_desc),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc')

!  Create the sd_model - makes sure that the signal_desc is properly pointed to
      DO isig = 1,n_s_desc
         CALL signal_construct(sd_model(isig),s_desc(isig),'model')
      END DO
      
      n1 = max(i_work(1),1)
      n2 = max(i_work(2),1)
      p1min = r_work(1)
      p1max = r_work(2)
      p2min = r_work(3)
      p2max = r_work(4)
      
      WRITE(49,*) 'gsq on grid'
      WRITE(49,*) 'n1, n2 = ', n1, n2
      WRITE(49,*) 'p1min, max = ', p1min, p1max
      WRITE(49,*) 'p2min, max = ', p2min, p2max
      WRITE(49,*) 
      WRITE(49,*) 'i1, i2, p1, p2, gsq'

!  Allocate space for gsq
      ALLOCATE(gsq(n1,n2),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc gsq')
      
      DO i1 = 1,n1
         p1 = p1min + (i1 - 1) * (p1max - p1min) / max(n1-1,1)
         rparam(1) % value = p1
         DO i2 = 1,n2
            p2 = p2min + (i2 - 1) * (p2max - p2min) / max(n2-1,1)
            rparam(2) % value = p2
            CALL recon_param_model_put(rparam(1:2),model_a,rcnstrnts,          &
     &         'VMI')
            
            numiter = state_vec % fixp % niter
            CALL eq_step(numiter,state_vec,lconverged = lconv,                 &
     &         iter_es = iter_v)
            WRITE(*,*) i1, i2, '   lconverged = ', lconv,                      &
     &         ' iter_ha = ', iter_v
            CALL signal_mc_model_compute(sd_model,model_a)
            
            CALL gsq_evaluate_g2(sd_model,sd_observe(1:n_s_desc))
            gsq(i1,i2) = gsq_g2
            
            WRITE(49,fmt1) i1, i2, p1, p2, gsq(i1,i2)
         END DO
      END DO

      WRITE(49,*)
      WRITE(49,*) ' rewrite gsq for Mma'
      DO i1 = 1,n1
         p1 = p1min + (i1 - 1) * (p1max - p1min) / max(n1-1,1)
         CALL mmaw_r2s(p1,12,p1_str)
         DO i2 = 1,n2
            p2 = p2min + (i2 - 1) * (p2max - p2min) / max(n2-1,1)
            CALL mmaw_r2s(p2,12,p2_str)
            CALL mmaw_r2s(gsq(i1,i2),12,gsq_str)
            WRITE(49,'("{",a,",",a,",",a,"}")') TRIM(p1_str),                  &
     &         TRIM(p2_str), TRIM(gsq_str)
         END DO
      END DO
      WRITE(49,*) 'end gsq on grid'

!  write out the vmec_history
      CALL eq_history_print

      RETURN
      END SUBROUTINE task_gsq_on_grid

!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------
!
      SUBROUTINE task_vmec_changep
!  Task subroutine, to start to run VMEC from a Namelist input file, and then
!  change a parameter and continue running

!  08-21-06. Added coding so that also recompute signals just one iteration
!   after change of parameter.

!  09-11-06. Changed coding to use recon_param
!  11-29-2006. Commented out the one-step stuff
!  2008-05-21 JDH
!    Pick a direction in parameter space, and go as far as possible in that direction.

!  2011-01-05 Rewritten, simplified
      USE v3f_global
      USE signal_mc
      USE recon_param_T
      USE recon_param_model
      USE gsq_mod
      
      IMPLICIT NONE

!----------------------------------------------------------------------
! L O C A L Variable Declarations
!----------------------------------------------------------------------
      TYPE(eq_state), POINTER :: state_vec => null()
      TYPE(signal_data), DIMENSION(:), ALLOCATABLE :: sd_model
      TYPE(recon_param), DIMENSION(:), ALLOCATABLE :: rp_old, rp_new
      INTEGER       :: numiter, isig, ier1, ier2, i1, irp
      INTEGER       :: iter_v
      LOGICAL              :: lconv, lsuccess
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: rpvra_dr

      CHARACTER(len=*), PARAMETER :: fmt1 =                                    &
     &  '(i5,1x,a30,1x,2(es20.13,2x))'
      CHARACTER(len=*), PARAMETER :: fmt2 =                                    &
     &  '(i5,1x,a30,1x,es11.4)'  
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'task_vmec_changep_2: '
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
!  point to equilibrium state
      state_vec => model_a % eqstate

!  Allocate space for signal data stuff.
!  Allocate sd_model to be number of observed signals
!  Note that sd_observe is allocated larger.

      CALL assert_eq(n_sdata_o,n_s_desc,sub_name //                            &
     &       ' problem with n_sdata_o, n_s_desc ') 
      ALLOCATE(sd_model(n_s_desc),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc sd_model')

!  Create the sd_model - makes sure that the signal_desc is properly pointed to
      DO isig = 1,n_s_desc
         CALL signal_construct(sd_model(isig),s_desc(isig),'model')
      END DO

!  Iterate VMEC to equilibrium      
      numiter = state_vec % fixp % niter
      CALL eq_history_set(i1 = 0)   !  Reconstruction step zero
      CALL eq_step(numiter,state_vec,lconverged = lconv,                       &
     &   iter_es = iter_v)
      WRITE(*,*) ' Initial VMEC equilibration'
      WRITE(iou_runlog,*) ' Initial VMEC equilibration'
      WRITE(*,*)  '   lconverged = ', lconv, ' iter_ha = ', iter_v
      WRITE(iou_runlog,*)  '   lconverged = ', lconv,                          &
     &   ' iter_ha = ', iter_v
      CALL assert(lconv,sub_name // ' no initial convergence')

!-----------------------------------------------
!  Compute the signals. Write result to iou_recout
      CALL signal_mc_model_compute(sd_model(isig),model_a)
      CALL signal_data_ma_write(sd_model(1:n_s_desc),unit=iou_recout)
      
!-----------------------------------------------
 
!  Allocate space for reconstruction parameters and arrays
      ALLOCATE(rp_old(n_rparam),rp_new(n_rparam),STAT=ier1)
      ALLOCATE(rpvra_dr(n_rparam),STAT=ier2)
      CALL assert_eq(0,ier1,ier2,sub_name // 'alloc rp_')

!  Get the initial reconstruction parameters, g2
      rp_old = rparam(1:n_rparam)
      CALL recon_param_model_get(rp_old,model_a)

!  Change the reconstruction parameters in rp_new
      rp_new = rp_old
      DO irp = 1, n_rparam
         rpvra_dr(irp) = r_work(irp)
      END DO
      CALL recon_param_change_value(rp_new,delta_value = rpvra_dr)

!  Print out information about parameters
      CALL recon_param_write_ond(rp_old,rp_new,unit=iou_recout)

!  Change the parameters in VMEC
      CALL recon_param_model_put(rp_new,model_a,rcnstrnts,'VMI')

!  Iterate VMEC to equilibrium at new parameter values     
      CALL eq_history_set(i1 = 1)   
      numiter = state_vec % fixp % niter
      CALL eq_step(numiter,state_vec,lconverged = lconv,                       &
     &         iter_es = iter_v)
     
!-----------------------------------------------
!  Compute the signals. Write result to iou_recout
      CALL signal_mc_model_compute(sd_model(isig),model_a)
      CALL signal_data_ma_write(sd_model(1:n_s_desc),unit=iou_recout)

!-----------------------------------------------
!  write out the vmec_history
      CALL eq_history_print

      RETURN
      END SUBROUTINE task_vmec_changep

!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------
!
      SUBROUTINE task_vmec_change_one_p
!  Task subroutine, to start to run VMEC from a Namelist input file, and then
!  change JUST ONE parameter and continue running, multiple times.

!  2011-01-12. First version, based on task_vmec_changep

      USE v3f_global
      USE signal_mc
      USE recon_param_T
      USE recon_param_model
      USE gsq_mod
      
      IMPLICIT NONE

!----------------------------------------------------------------------
! L O C A L Variable Declarations
!----------------------------------------------------------------------
      TYPE(eq_state), POINTER :: state_vec => null()
      TYPE(signal_data), DIMENSION(:), ALLOCATABLE :: sd_model
      TYPE(recon_param), DIMENSION(:), ALLOCATABLE :: rp_old, rp_new
      INTEGER       :: numiter, isig, ier1, ier2, i1, j, jmax, ii, jj
      INTEGER       :: iub1
      INTEGER       :: iter_v
      INTEGER, PARAMETER       :: ioul=37
      LOGICAL       :: lconv, lsuccess
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: rpvra_dr
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: sd1_loop_array,              &
     &  sd2_loop_array

      CHARACTER(len=*), PARAMETER :: fmt1 =                                    &
     &  '(i5,1x,a30,1x,2(es20.13,2x))'
      CHARACTER(len=*), PARAMETER :: fmt2 =                                    &
     &  '(i5,1x,a30,1x,es11.4)'  
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'task_vmec_change_one_p: '
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
!  point to equilibrium state
      state_vec => model_a % eqstate

!  Allocate space for signal data stuff.
!  Allocate sd_model to be number of observed signals
!  Note that sd_observe is allocated larger.

      CALL assert_eq(n_sdata_o,n_s_desc,sub_name //                            &
     &       ' problem with n_sdata_o, n_s_desc ') 
      ALLOCATE(sd_model(n_s_desc),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc sd_model')

!  Create the sd_model - makes sure that the signal_desc is properly pointed to
      DO isig = 1,n_s_desc
         CALL signal_construct(sd_model(isig),s_desc(isig),'model')
      END DO

!  Should only be ONE reconstruction parameter (Warning only)
      CALL assert_eq(1,n_rparam,sub_name // 'n_rparam should be 1','W') 

!  Allocate space for reconstruction parameters and arrays
      ALLOCATE(rp_old(n_rparam),rp_new(n_rparam),STAT=ier1)
      ALLOCATE(rpvra_dr(n_rparam),STAT=ier2)
      CALL assert_eq(0,ier1,ier2,sub_name // 'alloc rp_')

!  Get the initial reconstruction parameters, g2
      rp_old = rparam(1:n_rparam)
      CALL recon_param_model_get(rp_old,model_a)

!  Allocate space for signal - loop arrays
      iub1 = UBOUND(r_work,1)
      ALLOCATE(sd1_loop_array(n_s_desc,0:iub1),STAT=ier1)
      ALLOCATE(sd2_loop_array(n_s_desc,0:iub1),STAT=ier2)
      CALL assert_eq(0,ier1,ier2,sub_name // 'alloc sd_loop')

!  Loop over successive changes to first reconstruction parameter
      j = 0
      DO

!  Iterate VMEC to equilibrium      
         numiter = state_vec % fixp % niter
         CALL eq_history_set(i1 = j)   !  Identify loop iteration
         CALL eq_step(numiter,state_vec,lconverged = lconv,                    &
     &      iter_es = iter_v)
         WRITE(*,*) ' VMEC equilibration #, converged, iter_ha ', j,           &
     &      lconv, iter_v
         CALL assert(lconv,sub_name // ' no initial convergence')

! Compute model signals
      CALL signal_mc_model_compute(sd_model,model_a)
 
      DO isig = 1,n_s_desc
         sd1_loop_array(isig,j) = sd_model(isig) % data(1)
         sd2_loop_array(isig,j) = sd_model(isig) % data(2)
      END DO
      CALL signal_data_ma_write(sd_model(1:n_s_desc),unit=iou_recout)
!-----------------------------------------------
 
!  Change the reconstruction parameters in rp_new
         j = j + 1
         IF (r_work(j) .eq. zero) THEN
            jmax = j - 1
            EXIT                                  ! the main DO loop
         ENDIF
         rp_new = rp_old
         rpvra_dr = zero
         rpvra_dr(1) = r_work(j)      ! only change first recon_param
         CALL recon_param_change_value(rp_new,delta_value = rpvra_dr)

!  Print out information about parameters
         CALL recon_param_write_ond(rp_old,rp_new,unit=iou_recout)

!  Change the parameters in VMEC
         CALL recon_param_model_put(rp_new,model_a,rcnstrnts,'VMI')
         WRITE(*,*) 'Recon Param type, index, index2, value ',                 &
     &      rp_new(1) % p_type, rp_new(1) % index, rp_new(1) % index2,         &
     &      rp_new(1) % value
         WRITE(iou_recout,*) 'Recon Param type, index, index2, value ',        &
     &      rp_new(1) % p_type, rp_new(1) % index, rp_new(1) % index2,         &
     &      rp_new(1) % value
         rp_old = rp_new

      END DO  ! Loop over successive changes to first reconstruction parameter

!-----------------------------------------------
!  write out the sdx_loop arrays

      WRITE(ioul,*) 'OUTPUT from task_vmec_change_one_p'
      WRITE(ioul,*) 'First, signal data(1), second - signal data(2)'
      WRITE(ioul,2000) r_work(1:jmax)
      DO ii = 1, n_s_desc      
         WRITE(ioul,2100) ii, (sd1_loop_array(ii,jj), jj = 0, jmax)
      END DO
      WRITE(ioul,'(//)')
      WRITE(ioul,2000) r_work(1:jmax)
      DO ii = 1, n_s_desc  
         WRITE(ioul,2100) ii, (sd2_loop_array(ii,jj), jj = 0, jmax)
      END DO
      
2000  FORMAT(/,'r_work(1:jmax)',t21,10(1x,es16.9))
2100  FORMAT(i4,10(1x,es16.9))

!-----------------------------------------------
!  write out the vmec_history
      CALL eq_history_print

      RETURN
      END SUBROUTINE task_vmec_change_one_p
!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------
!
      SUBROUTINE task_vmec_change_rp
!  Task subroutine, to start to run VMEC from a Namelist input file, and then
!  change reconstruction parameters and continue running, multiple times.

!  2011-06-07. First version, based on task_vmec_change_one_p

!  All the reconstruction parameters will be changed. Note that the value of a
!  reconstruction parameter specified in the VMEC NLI file is NOT USED.

!  irp - index to the reconstruction parameters
!  nrp - number of reconstruction parameters

!  Significance of the iwork and rwork arrays:
!    i_work
!  1       nstep             number of steps (Different values of recon parameters)
!  2       rp_var_type(1)    integer to specify type of rp variation
!  ...     ...                 0 - linear
!  ...     ...                 1 - cosine
!  ...     ...                 2 - sine
!  nrp+1   rp_var_type(nrp)  

!    r_work
!  1       p_a(1)       Parameter specification a for first recon param
!  2       p_b(1)       Parameter specification b for first recon param
!  ...                     Define dx = (istep - 1) / (nstep - 1): 0 <= dx <= 1
!  ...                     linear, rp_var_type(irp) = 0
!  ...                        p = p_a + dx * (p_b - p_a) 
!  ...                     cosine, rp_var_type(irp) = 1
!  ...                        p = p_a + p_b * cos(dx * twopi)
!  ...                     sine, rp_var_type(irp) = 2
!  ...                        p = p_a + p_b * sin(dx * twopi)
!  2 * nrp  p_b(nrp)

      USE v3f_global
      USE signal_mc
      USE recon_param_T
      USE recon_param_model
      USE gsq_mod
      
      IMPLICIT NONE

!----------------------------------------------------------------------
! L O C A L Variable Declarations
!----------------------------------------------------------------------
      TYPE(eq_state), POINTER :: state_vec => null()
      TYPE(signal_data), DIMENSION(:), ALLOCATABLE :: sd_model
      TYPE(recon_param), DIMENSION(:), ALLOCATABLE :: rpT
      
      INTEGER :: ier1, ier2, iter_v, numiter, isig
      INTEGER :: istep, nstep, irp, nrp
      INTEGER, DIMENSION(:), ALLOCATABLE :: rp_var_type
      
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: p_a, p_b, rp_vals
      
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: sd1_save_array,              &
     &  sd2_save_array, rp_vals_save_array

      REAL(rprec)   :: dx
      LOGICAL       :: lconv

      CHARACTER(len=*), PARAMETER :: fmt1 =                                    &
     &  '(i5,3x,i2,3x,2(es20.13,2x))'
       CHARACTER(len=*), PARAMETER :: sub_name =                               &
     &  'task_vmec_change_rp: '
!----------------------------------------------------------------------
! L O C A L Variables
!----------------------------------------------------------------------

!      Derived Types
!  rpT          type Recon_Param, vector, old values
!  state_vec    type eq_state, pointer to component of model_a
!  sd_model     type signal_data, array, model values

!      Local Scalars
!  ier1         integer, error flag
!  istep        index for step loop
!  nstep        number of steps (different set of recon param values)
!  irp          index for loops over reconstruction parameters
!  nrp          number of reconstruction parameters

!      Local 1-d arrays
!  rp_var_type  integer, array of type of recon param variation
!  p_a          real array, recon parameter specification a
!  p_b          real array, recon parameter specification b
!  rp_vals      real array, values of reconstruction parameters

!      Local 2-d arrays
!  sd1_save_array        signal data (1) : n_s_desc X nstep
!  sd1_save_array        signal data (2) : n_s_desc X nstep
!  rp_vals_save_array    recon param values : nrp X nstep

!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
!-----------------------------------------------
!  Initialization and Allocation

      WRITE(iou_runlog,*) ' task_vmec_change_rp:  nrp = ', n_rparam
      WRITE(*,*) ' task_vmec_change_rp:  nrp = ', n_rparam

!  Convert from work arrays to local variables
      nrp = n_rparam                              ! n_rparam from v3f_global
      nstep = MAX(1,i_work(1))                    ! i_work from v3f_global
      ALLOCATE(rp_var_type(n_rparam),STAT=ier1) 
      CALL assert_eq(0,ier1,sub_name // 'alloc rp_var_type')
      ALLOCATE(p_a(n_rparam),p_b(n_rparam),STAT=ier1)   
      CALL assert_eq(0,ier1,sub_name // 'alloc p_a, p_b')
      DO irp = 1,nrp
         rp_var_type(irp) = i_work(irp + 1)
         IF ((rp_var_type(irp) .lt. 0) .or. (rp_var_type(irp) .gt. 2))         &
     &      CALL err_fatal(sub_name // 'rp_var_type out of range ' //          &
     &      ' index  is',int=irp)
         p_a(irp) = r_work(2 * irp - 1)
         p_b(irp) = r_work(2 * irp)
      END DO
      
!  point to equilibrium state. model_a from v3f_global
      state_vec => model_a % eqstate

!  Allocate space for model signal data stuff. n_s_desc from v3f_global.
      ALLOCATE(sd_model(n_s_desc),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc sd_model')

!  Create the sd_model - makes sure that the signal_desc is properly pointed to
      DO isig = 1,n_s_desc
         CALL signal_construct(sd_model(isig),s_desc(isig),'model')
      END DO

!  Should at least one reconstruction parameter
      CALL assert(n_rparam .ge. 1,sub_name // 'n_rparam .lt. 1') 

!  Allocate space for reconstruction parameters and arrays
      ALLOCATE(rpT(n_rparam),STAT=ier1)
      ALLOCATE(rp_vals(n_rparam),STAT=ier2)
      CALL assert_eq(0,ier1,ier2,sub_name // 'alloc rp_')

!  Allocate space for local 2-d arrays
      ALLOCATE(sd1_save_array(n_s_desc,nstep),                                 &
     &   sd2_save_array(n_s_desc,nstep) ,STAT=ier1)
      ALLOCATE(rp_vals_save_array(nrp,nstep),STAT=ier2)
      CALL assert_eq(0,ier1,ier2,sub_name // 'alloc 2d')

!  Get the reconstruction parameters from the state
      rpT = rparam(1:n_rparam)
      CALL recon_param_model_get(rpT,model_a)

!-----------------------------------------------
!  Printout before the step loop
      WRITE(iou_runlog,*) 
      WRITE(*,*) 
      WRITE(iou_runlog,*) '   irp      var_type   p_a   p_b'
      WRITE(*,*) '   irp      var_type   p_a   p_b'
      DO irp = 1,nrp
         WRITE(iou_runlog,fmt1) irp, rp_var_type(irp), p_a(irp),               &
     &      p_b(irp)
         WRITE(*,fmt1) irp, rp_var_type(irp), p_a(irp),                        &
     &      p_b(irp)
      END DO
      WRITE(iou_recout,*) 
      WRITE(*,*) 

!-----------------------------------------------
!  Step Loop
!-----------------------------------------------
      DO istep = 1,nstep
         WRITE(iou_runlog,'(//)')
         WRITE(iou_runlog,*) '*************** istep = ', istep
         WRITE(*,*) 'istep = ', istep
!  Define the reconstruction parameter values
         DO irp = 1,nrp
            IF (nstep .le. 1) THEN
               dx = zero
            ELSE
               dx = (istep - 1.) / (nstep - 1.)
            ENDIF
            
            SELECT CASE(rp_var_type(irp))
            
            CASE(0) !  linear
               rp_vals(irp) = p_a(irp) + dx * (p_b(irp) - p_a(irp))
            
            CASE(1) ! Cosine
               rp_vals(irp) = p_a(irp) + p_b(irp) * COS(dx * twopi)

            CASE(2) ! Sine
               rp_vals(irp) = p_a(irp) + p_b(irp) * SIN(dx * twopi)
            
            CASE DEFAULT
               CALL err_fatal(sub_name // 'rp_var_type out of range '//        &
     &         ' index  is',int=irp)

            END SELECT
         END DO
         rp_vals_save_array(:,istep) = rp_vals(:)
        
!  Set the new reconstruction parameter values. Change the parameters in VMEC.
         CALL recon_param_change_value(rpT,new_value = rp_vals)
         CALL recon_param_model_put(rpT,model_a,rcnstrnts,'VMI')

!  Print out the reconstruction parameters
         CALL recon_param_write(rpT,unit=iou_runlog)
           
!  Iterate VMEC to equilibrium      
         numiter = state_vec % fixp % niter
         CALL eq_history_set(i1 = istep)   !  Identify step
         CALL eq_step(numiter,state_vec,lconverged = lconv,                    &
     &      iter_es = iter_v)
         WRITE(*,*) ' VMEC equilibration #, converged, iter_ha ',istep,        &
     &      lconv, iter_v
         CALL assert(lconv,sub_name // ' no initial convergence')

!  Compute the signals. 
         CALL signal_mc_model_compute(sd_model,model_a)
         DO isig = 1,n_s_desc
            sd1_save_array(isig,istep) = sd_model(isig) % data(1)
            sd2_save_array(isig,istep) = sd_model(isig) % data(2)
         END DO
         
         CALL signal_data_ma_write(sd_model(1:n_s_desc),                       &
     &      unit=iou_runlog)
     
!-----------------------------------------------
!  End of the Step Loop
      END DO
     
!-----------------------------------------------
!  Write to the recout file
      WRITE(iou_recout,'(//)')
      WRITE(iou_recout,*) 'task_vmec_change_rp.  nstep = ', nstep
      WRITE(iou_recout,'(/)')
      WRITE(iou_recout,*) '  *******  Reconstruction Parameters'
      DO irp = 1,nrp
         WRITE(iou_recout,'()')
         WRITE(iou_recout,*) 'irp, rp_var_type:', irp, rp_var_type(irp)
         WRITE(iou_recout,*) 'p_a, p_b:', p_a(irp), p_b(irp)
         WRITE(iou_recout,*) 'Type, indices: ', rpT(irp) % p_type,             & 
     &      rpT(irp) % index, rpT(irp) % index2
         WRITE(iou_recout,*) 'Values:'
         WRITE(iou_recout,*) (rp_vals_save_array(irp,istep),                   & 
     &      istep = 1,nstep)
      END DO
      
      WRITE(iou_recout,'(//)')
      WRITE(iou_recout,*) '  *******  Signal Data (1)'
      WRITE(iou_recout,*)                                                      &
     &   '     For magnetic diagnostics, (1) is the total signal'
      DO isig = 1,n_s_desc
         WRITE(iou_recout,'()')
         WRITE(iou_recout,*) 'isig = ', isig,'  (1) Values:'
         WRITE(iou_recout,*) (sd1_save_array(isig,istep),                      & 
     &      istep = 1,nstep)
      END DO
      
      WRITE(iou_recout,'(//)')
      WRITE(iou_recout,*) '  *******  Signal Data (2)'
      WRITE(iou_recout,*)                                                      &
     &   '     For magnetic diagnostics, (2) is the plasma signal'
      DO isig = 1,n_s_desc
         WRITE(iou_recout,'()')
         WRITE(iou_recout,*) 'isig = ', isig,'  (2) Values:'
         WRITE(iou_recout,*) (sd2_save_array(isig,istep),                      & 
     &      istep = 1,nstep)
      END DO

      WRITE(iou_recout,*) 'Done with task_vmec_change_rp. '

!-----------------------------------------------
!  Write out the vmec_history
      CALL eq_history_print

      RETURN
      END SUBROUTINE task_vmec_change_rp

!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------
!
      SUBROUTINE task_vmec_monitor
!  Task to run vmec, stopping to compute signals along the way.

      USE v3f_global
      USE signal_mc
      USE eq_interface
      
       IMPLICIT NONE

!  i_work(1)     From v3f_global. Used to specify number of iterations between
!                signals computations.

!----------------------------------------------------------------------
! L O C A L Variable Declarations
!----------------------------------------------------------------------
!  iter_reset   Integer, to save number f iterations when VMEC resets its
!               internal counter (which it does when incrementing ns_index)
!  iter_v       Integer, vmec internal iteration counter, comes from eq_step

      TYPE(eq_state), POINTER :: state_vec => null()
      TYPE(signal_data), DIMENSION(:), ALLOCATABLE :: sd_model
      INTEGER  :: i, numiter, isig, j, iter_reset, iter_v, imax, ier1
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: psignals
      LOGICAL :: lconv, lsuccess
      CHARACTER(len=*), PARAMETER :: fmt1 =                                    &
     &  '(i5,1x,i5,1x,es11.4,2x,20(es20.13,2x))'
      CHARACTER(len=400) :: signal_names = ' '
      CHARACTER(len=20) :: char_num
      CHARACTER(len=80) :: sys_command, new_filename
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'task_vmec_monitor: '
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
      ALLOCATE(psignals(n_s_desc))
      state_vec => model_a % eqstate

!  Allocate space for signal data stuff.
      ALLOCATE(sd_model(n_s_desc),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc sd_model')

!  Create the sd_model - makes sure that the signal_desc is properly pointed to
      DO isig = 1,n_s_desc
         CALL signal_construct(sd_model(isig),s_desc(isig),'model')
      END DO
      
      DO i = 1,n_s_desc
         signal_names = TRIM(signal_names) // ' ' //                           &
     &                  TRIM(s_desc(i) % s_name)
      END DO
      WRITE(47,*) 'ns_index, iter_cumul, fsq, ',TRIM(signal_names)
     
      numiter = i_work(1) 
      imax = state_vec % fixp % niter / numiter + 1
      iter_reset = 0
      DO i = 1,imax
         CALL eq_step(numiter,state_vec,lconverged = lconv,                    &
     &      iter_es = iter_v)
!         WRITE(*,*) 'iter_cumul, fsq_max, lconv:  ', iter_reset + iter_v,      &
!     &      state_vec % fsq_max, ' ', lconv
     
!  Compute the signals, temporarily store in psignals
         CALL signal_mc_model_compute(sd_model,model_a)
         DO isig = 1,n_s_desc
            psignals(isig) = sd_model(isig) % data(2)
         END DO

!  Write signals to text file
         WRITE(47,fmt1) state_vec % varp % ns_index,                           &
     &      iter_reset + iter_v, state_vec % fsq_max,                          &
     &      psignals(1:n_s_desc)

!  Copy the wout file to a unique name, if i_work(2) is 1. (Kludge, 06-30-06)
!  wout_.nc is the name when vmec_input:input_extension is blank
!  (Coding in VMEC subroutine wrout)
      IF (i_work(2) .eq. 1) THEN
         WRITE(char_num,'(i8)') iter_reset + iter_v
         new_filename = 'wout_'// TRIM(ADJUSTL(char_num))                      &
     &      // '.nc'
!  Construct a system call that will copy the wout file.
!  Note that the subroutine system is NOT standard Fortran.
         WRITE(sys_command,"('cp wout_.nc ',a)") TRIM(new_filename)
         CALL system(sys_command)
      ENDIF
         
!  If have converged at this value of ns_array, increment ns_index
!  If can't increment ns_index succesfully, then should be done.
         IF (lconv) THEN
            WRITE(*,*) 'Converged at this ns value'
            CALL eq_ns_index_inc(state_vec,lsuccess)
            IF (lsuccess) THEN
               WRITE(*,*) 'ns_index incremented'
               iter_reset = iter_reset + iter_v
            ELSE
               WRITE(*,*) 'ns_index not successfully incremented'
               WRITE(*,*) 'Assume done with largest ns'
               EXIT
            ENDIF
         END IF
      END DO
      
      WRITE(*,*) 'Completed calls to eq_step from task_vmec_monitor'

      RETURN
      END SUBROUTINE task_vmec_monitor

!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------
!
      SUBROUTINE task_vmec_v3post
!  Task subroutine, to first run VMEC from a Namlist input file, and then
!  compute the signals, mimicking V3POST

      USE v3f_global
      USE signal_mc
      USE eq_interface
      
      IMPLICIT NONE

!----------------------------------------------------------------------
! L O C A L Variable Declarations
!----------------------------------------------------------------------
      TYPE(eq_state), POINTER :: state_vec => null()
      TYPE(signal_data), DIMENSION(:), ALLOCATABLE :: sd_model
      INTEGER  :: numiter, isig, iter_v, ier1
      LOGICAL :: lconv, lsuccess
      CHARACTER(len=*), PARAMETER :: fmt1 =                                    &
     &  '(i5,1x,a30,1x,2(es20.13,2x))'
      CHARACTER(len=*), PARAMETER :: fmt2 =                                    &
     &  '(i5,1x,a30,1x,es11.4)'  !GJH added 2010-01-26
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'task_vmec_v3post: '
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
      state_vec => model_a % eqstate

!  Allocate sd_model 
      ALLOCATE(sd_model(n_s_desc),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc sd_model')

!  Create the sd_model - makes sure that the signal_desc is properly pointed to
      DO isig = 1,n_s_desc
         CALL signal_construct(sd_model(isig),s_desc(isig),'model')
      END DO

      numiter = state_vec % fixp % niter
!  Loop over ns_index values
      DO
         WRITE(*,*) ' ns_index = ',state_vec % varp % ns_index
         CALL eq_step(numiter,state_vec,lconverged = lconv,                    &
     &      iter_es = iter_v)
         IF (lconv) THEN
            !  VMEC converged at this ns_index value
            WRITE(*,*) '   converged, numiter, fsq_max, iter_v',               &
     &           numiter, state_vec % fsq_max, iter_v
            CALL eq_ns_index_inc(state_vec,lsuccess)
            IF (lsuccess) THEN
               ! ns_index successfully incremented, run vmec some more
               CYCLE
            ELSE
               ! ns_index not incremented, assume that done
               EXIT
            ENDIF
         ELSE
            WRITE(*,*) 'VMEC DID NOT CONVERGE'
         END IF
      END DO

!  Now, compute the signals. Write results
      CALL signal_mc_model_compute(sd_model,model_a)
      CALL signal_data_ma_write(sd_model(1:n_s_desc),unit=iou_recout)
!-----------------------------------------------

!  write out the vmec_history
      CALL eq_history_print

      RETURN
      END SUBROUTINE task_vmec_v3post

!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------
!
!
      SUBROUTINE task_add_noise_a1
!  Task subroutine, add noise to Reconstruction, Algorithm 1

!  2010-03-02  GJH - fsigma was input but not used. 
!                    Previously unit_noise*signal_sigma was added to data
!                    Now unit_noise*signal_sigma*fsigma is added to data
!  2007-07-03. First version - based on task_reconstruct_a1
!  2008-06-12. Changed _work indices, so don't collide with reconstruct_a1 values

      USE v3f_global
!  ONLY model_a, i_work, r_work + others
      USE signal_mc
      USE recon_param_T
      USE recon_param_model
      USE gsq_mod
      USE mmaw
      
      IMPLICIT NONE

!----------------------------------------------------------------------
! L O C A L Variable Declarations
!----------------------------------------------------------------------

!      
!      Namelist Input - data from v3f_global

!   ** Used by add_noise_a1 **
!  nrecon       i_work(15): number of reconstructions to do
!  fsigma       r_work(15): fraction sigma - scaling factor for 
!                  determining noise level to add
!  l_rseed      i_work(16): 0 (default) - do not use a random seed
!                          non-zero - use a random number seed

!      Derived Types
!  rp_original    Type recon_param, vector, original values - to save
!  rp_done        Type recon_param, vector, values at end of reconstruction
!  sdo_original   Type signal_data, vector, original values - to save

!      Save variables, for use in printout at end of subroutine
!  save_g2      vector to save gsq values (index - reconstruction number)
!  save_rp      two-d array of reals, save reconstruction parameter values
!               (1st index - reconstruction number, 2nd index - parameter)
!  save_g2_atexact  vector to save gsq values, evaluated at the exact parameters

!      Other local variables
!  sig_noise    vector of reals, signal noise levels
!  irstep_used  integer, argument to task_reconstruct_a1. Number of 
!                  reconstruction steps used
!               GJH - 2010-03-02 note that irstep_used is not being 
!                     used in this routine

      INTEGER  ::  nrecon, l_rseed
      REAL(rprec)     :: fsigma

      TYPE(recon_param), DIMENSION(:), ALLOCATABLE :: rp_original
      TYPE(recon_param), DIMENSION(:), ALLOCATABLE :: rp_done
      TYPE(signal_data), DIMENSION(:), ALLOCATABLE :: sdo_original

      REAL(rprec), DIMENSION(:), ALLOCATABLE    :: save_g2,                    &
     &   save_g2_atexact
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE  :: save_rp

      REAL(rprec), DIMENSION(:), ALLOCATABLE :: sig_noise
      CHARACTER(len=30)    :: rp1_str, rp2_str, gsq_str
      INTEGER :: nsdo, nrp, irecon, isdo, ier1, ier2, ier3,                    &
     &  irstep_used=0

      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'task_add_noise_a1: '
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
      WRITE(*,*) ' *** IN SUBROUTINE ', sub_name

      nrecon = i_work(15)
      fsigma = r_work(15)
      l_rseed = i_work(16)

!  Randomize seed for random number generation (or not)
      IF (l_rseed .ne. 0) THEN
         WRITE(*,*) sub_name,'Randomizing seed'
         CALL random_seed
      ENDIF
      
! Type signal_data initialization. Allocate, assign sufficient
!  sdo_original is local
!  sd_observe is from v3f_global
      nsdo = n_sdata_o
      ALLOCATE(sdo_original(nsdo),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc sdo_original')
!      DO isdo = 1,nsdo
!         CALL signal_construct(sdo_original(isdo),                             &
!     &      sd_observe(isdo) % desc,sd_observe(isdo) % sd_type,                &
!     &      sd_observe(isdo) % data,sd_observe(isdo) % sigma)
!      END DO
      sdo_original = sd_observe(1:nsdo)

!  Type recon_param initialization. Allocate and assign.
!  No construction necessary, as has no pointer components.
!  rp_original and rp_done are local
!  rparam is from v3f_global
      nrp = n_rparam
      ALLOCATE(rp_original(nrp),STAT=ier1)
      ALLOCATE(rp_done(nrp),STAT=ier2)
      CALL assert_eq(0,ier1,ier2,sub_name // 'alloc rp')
      rp_original = rparam(1:nrp)
      rp_done = rparam(1:nrp)
      CALL recon_param_model_get(rp_original,model_a)

!  Allocate other variables
      ALLOCATE(sig_noise(nsdo),STAT=ier1)
      ALLOCATE(save_g2(nrecon),STAT=ier2)
      ALLOCATE(save_g2_atexact(nrecon),STAT=ier3)
      CALL assert_eq(0,ier1,ier2,ier3,sub_name // 'alloc various 1')
      ALLOCATE(save_rp(nrecon,nrp),STAT=ier1)
      CALL assert_eq(0,ier1,sub_name // 'alloc various 2')

!  Loop over number of reconstructions
      DO irecon = 1,nrecon

!  Generate and apply noise to observations
         CALL gauss_rand(nsdo,sig_noise(1:nsdo))
         sig_noise(1:nsdo)=sig_noise(1:nsdo)*fsigma  ! GJH 2010-03-02 
                                                     ! added fsigma product
         save_g2_atexact(irecon) = DOT_PRODUCT(sig_noise(1:nsdo),              &
     &      sig_noise(1:nsdo))
     
         DO isdo = 1,nsdo
            sig_noise(isdo) = sig_noise(isdo) *                                &
     &         sdo_original(isdo) % sigma(1)  ! GJH 2010-03-02 added fsigma 
            sd_observe(isdo) % data(1) = sdo_original(isdo) % data(1)          &
     &         + sig_noise(isdo)
         END DO

!  Start at original reconstruction parameter values. model_a is in v3f_global
         CALL recon_param_model_put(rp_original,model_a,rcnstrnts,             &
     &      'VMI')

!  Do the reconstruction
         CALL task_reconstruct_a1

!  Save g^2, reconstructed parameter values
         save_g2(irecon) = gsq_g2
         CALL recon_param_model_get(rp_done,model_a)
         save_rp(irecon,1:nrp) = rp_done(1:nrp) % value

!  Write to standard output, and a separate file
         WRITE(*,*) 
         WRITE(*,*) '@@@@@@@@@@@@@@@'
         WRITE(*,*) 'irecon, g^2, g^2 at exact,/ rp'
         WRITE(*,*) irecon,  save_g2(irecon),
     &      save_g2_atexact(irecon)
         WRITE(*,*) save_rp(irecon,1:nrp)
         WRITE(*,*) '@@@@@@@@@@@@@@@ '
         WRITE(*,*) 
         WRITE(63,*)  irecon, irstep_used, save_g2(irecon),                    &
     &       save_g2_atexact(irecon), save_rp(irecon,1:nrp)
      END DO

!  Rewrite so that easy for Mma input - only for 2 parameters
      IF(nrp .eq. 2) THEN
         DO irecon = 1,nrecon
            CALL mmaw_r2s(save_g2(irecon),12,gsq_str)
            CALL mmaw_r2s(save_rp(irecon,1),12,rp1_str)
            CALL mmaw_r2s(save_rp(irecon,2),12,rp2_str)
            WRITE(63,'("{",a,",",a,",",a,"},")') TRIM(rp1_str),                &
     &         TRIM(rp2_str), TRIM(gsq_str)
         END DO
      ENDIF
          
      RETURN
      END SUBROUTINE task_add_noise_a1
               
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
      SUBROUTINE task_test
!  Task subroutine, generic test
!  Developers should feel free to use this for any testing that they wish.
!  If this task will be used by others, change the task name to something 
!  other than test

!  2008-02-22 JDH Test most_redundant subroutine

      USE stel_kinds
      USE v3_utilities

      IMPLICIT NONE

      REAL(rprec), DIMENSION(6,5) :: atest
      REAL(rprec), DIMENSION(5) :: svprod_array
      INTEGER, DIMENSION(5) :: ncol_array, j_col_elim
      INTEGER :: i,j
      
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'task_test: '
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------
      WRITE(*,*) 'most_redundant test, random 6 by 5'

      CALL random_number(atest)
      atest = 1. + atest
      DO i = 1,6
          WRITE(*,1100) i,(atest(i,j),j=1,5)
      END DO
1100  FORMAT(2x,i2,5(3x,es12.3))
      WRITE(*,*) 
      
      CALL most_redundant(atest,ncol_array,svprod_array,j_col_elim)
      WRITE(*,*) 'i   col-remain         svprod             j_col_elim'
1000  FORMAT(i2,4x,i2,6x,es12.3,4x,i2)
      WRITE(*,1000) (i, ncol_array(i), svprod_array(i), j_col_elim(i),         &
     &   i=1,5)
      WRITE(*,*) 
      WRITE(*,*) 
      
      WRITE(*,*) 'most_redundant test, Column 2 down by factor to 10.'
      DO i = 1,6
         atest(i,2) = 1.e-1 * atest(i,2)
      END DO
      DO i = 1,6
          WRITE(*,1100) i,(atest(i,j),j=1,5)
      END DO

      CALL most_redundant(atest,ncol_array,svprod_array,j_col_elim)
      WRITE(*,*) 'i   col-remain         svprod             j_col_elim'
      WRITE(*,1000) (i, ncol_array(i), svprod_array(i), j_col_elim(i),         &
     &   i=1,5)
      WRITE(*,*) 
      WRITE(*,*) 
      
      WRITE(*,*) 'most_redundant test, Column 2 down by factor 100'
      DO i = 1,6
         atest(i,2) = 1.e-1 * atest(i,2)
      END DO
      DO i = 1,6
          WRITE(*,1100) i,(atest(i,j),j=1,5)
      END DO

      CALL most_redundant(atest,ncol_array,svprod_array,j_col_elim)
      WRITE(*,*) 'i   col-remain         svprod             j_col_elim'
      WRITE(*,1000) (i, ncol_array(i), svprod_array(i), j_col_elim(i),         &
     &   i=1,5)
      WRITE(*,*) 
      WRITE(*,*) 

      WRITE(*,*) 'most_redundant test, Column 2 down by factor 1000'
      DO i = 1,6
         atest(i,2) = 1.e-1 * atest(i,2)
      END DO
      DO i = 1,6
          WRITE(*,1100) i,(atest(i,j),j=1,5)
      END DO

      CALL most_redundant(atest,ncol_array,svprod_array,j_col_elim)
      WRITE(*,*) 'i   col-remain         svprod             j_col_elim'
      WRITE(*,1000) (i, ncol_array(i), svprod_array(i), j_col_elim(i),         &
     &   i=1,5)
      WRITE(*,*) 
      WRITE(*,*) 

      WRITE(*,*) 'most_redundant test, 1,2,3,4,5'
      CALL random_number(atest)
      atest = 1. + atest
      DO i = 1,6
         atest(i,2) = 1.e-1 * atest(i,2)
         atest(i,3) = 1.e-3 * atest(i,3)
         atest(i,4) = 1.e-6 * atest(i,4)
         atest(i,5) = 1.e-10 * atest(i,5)
      END DO
      DO i = 1,6
          WRITE(*,1100) i,(atest(i,j),j=1,5)
      END DO

      CALL most_redundant(atest,ncol_array,svprod_array,j_col_elim)
      WRITE(*,*) 'i   col-remain         svprod             j_col_elim'
      WRITE(*,1000) (i, ncol_array(i), svprod_array(i), j_col_elim(i),         &
     &   i=1,5)
      WRITE(*,*) 
      WRITE(*,*) 
      
      WRITE(*,*) 'most_redundant test, 5,4,3,2,1'
      CALL random_number(atest)
      atest = 1. + atest
      DO i = 1,6
         atest(i,4) = 1.e-1 * atest(i,4)
         atest(i,3) = 1.e-3 * atest(i,3)
         atest(i,2) = 1.e-6 * atest(i,2)
         atest(i,1) = 1.e-10 * atest(i,1)
      END DO
      DO i = 1,6
          WRITE(*,1100) i,(atest(i,j),j=1,5)
      END DO

      CALL most_redundant(atest,ncol_array,svprod_array,j_col_elim)
      WRITE(*,*) 'i   col-remain         svprod             j_col_elim'
      WRITE(*,1000) (i, ncol_array(i), svprod_array(i), j_col_elim(i),         &
     &   i=1,5)
      WRITE(*,*) 
      WRITE(*,*) 
      
      WRITE(*,*) 'most_redundant test, 3,5 zero'
      CALL random_number(atest)
      atest = 1. + atest
      DO i = 1,6
         atest(i,3) = 0.
         atest(i,5) = 0.
      END DO
      DO i = 1,6
          WRITE(*,1100) i,(atest(i,j),j=1,5)
      END DO

      CALL most_redundant(atest,ncol_array,svprod_array,j_col_elim)
      WRITE(*,*) 'i   col-remain         svprod             j_col_elim'
      WRITE(*,1000) (i, ncol_array(i), svprod_array(i), j_col_elim(i),         &
     &   i=1,5)
      WRITE(*,*) 
      WRITE(*,*) 
      
      RETURN
      END SUBROUTINE task_test
 
!===============================================================================
      SUBROUTINE task_bfield_grid
!-------------------------------------------------------------------------------
!  Task subroutine to invoke the bfield_grid_mod module, which computes a 2D grid
!      of B field values.
!
!  Created 10/23/06 by J. Shields
!-------------------------------------------------------------------------------
!      USE v3f_global
!      USE signal_mc
!      USE eq_interface
      USE bfield_grid_mod
      
      IMPLICIT NONE

!----------------------------------------------------------------------
! L O C A L Variable Declarations
!----------------------------------------------------------------------
      INTEGER  :: i,j,k,m
      REAL(rprec)     :: p_area_min, p_area_max

      call B_GRID_DRIVER

      RETURN
      END SUBROUTINE task_bfield_grid

!  2009-03-05 JDH - Deleted tasks 'ray_tracing' and 'rtrace_circle', along with
!     corresponding subroutines task_ray_tracing and task_rtrace_circle

!*******************************************************************************
! SECTION V.	SUBROUTINES FOR TESTNG
!*******************************************************************************

      SUBROUTINE getmem(mem)
!  Subroutine to get the memory used by the executing code
!  4 August 2004   James D. Hanson
!  24 September 2004 JDH - Change algorithm, so works on logjam, also.
!  2009-08-20 JDH. Changed ps commmand, so that works on Solaris machines also
!  2011-08-26 JDH. Fixed problem when executable was a full-path-name.
!  Last Revised: 2011-08-26
!-------------------------------------------------------------------------------
! The subroutine uses non-standard intrinsic functions to determine the
! memory size of the executing code.
! The current version works for:
! 1)   the IBM XLF compiler on a Macintosh running OS X
! 2)   Lahey-Fujitsu 6.2 on a Linux machine (logjam)
!-------------------------------------------------------------------------------

      IMPLICIT NONE
      
!  Declare arguments
      INTEGER :: mem

!    Declare local variables
!  iou             I/O unit number, for temporary file.
!  odd_file_name   file name of temporary file
!  sys_command     character variable to hold system call
!  executable      character, to hold the name of this executable code
!  line            character variable that holds a line of the result of the ps command
!  istat           integer status variable
!  line_number     integer variable to hold the line number where the executable occurs
!  i               do loop counter
!  ilast_slash     index to the location of the last '/' in the executable name

      INTEGER :: iou = 99
      CHARACTER(len=80) :: odd_file_name = 'w3rt7fkv'
      CHARACTER(len=80) :: sys_command
      CHARACTER(len=80) :: executable
      CHARACTER(len=80) :: line
      INTEGER :: istat, line_number, i, ilast_slash

#if defined(WIN32)
      mem = 0
      RETURN
#endif

!  Start of executable code

!  Find the name of this executable
      CALL getarg(0,executable)
      executable = TRIM(ADJUSTL(executable))
!      WRITE(*,*) ' Executable is ', executable

!  The executable name may be a full path. Toss out all but the last
!  component.
      ilast_slash = INDEX(executable,'/',.TRUE.)
      executable = executable(ilast_slash + 1:)
!      WRITE(*,*) ' Last component of executable is ', executable

!  Construct a system call that will generate output that contains
!  the memory. 
      WRITE(sys_command,"('ps -o vsz -o comm > ',a)")                          &
     &   TRIM(odd_file_name)

!      WRITE(*,*) sys_command
      
!  Execute the command
      CALL system(sys_command)
      
!      WRITE(*,*) 'Executed system call'
      
!  Open the file that the result of the ps command got written to
      OPEN(iou,file=odd_file_name,iostat=istat)
      
!      WRITE(*,*) 'Opened file ', TRIM(odd_file_name)
!      WRITE(*,*) 'iostat is ', istat

      IF (istat .ne. 0) THEN
         WRITE(*,*) 'Error opening file in getmem. iostat = ',istat
         mem = 0
         RETURN
      END IF
      
! Read through the file, looking for the line that contains the executable
      line_number = 0
      DO i = 1,1000
         READ (iou, '(a)', iostat=istat) line
         IF (istat .ne. 0) THEN
            mem = 0
            WRITE(*,*) 'getmem: istat .ne. 0, i = ',i
            WRITE(*,*) '  Could not find executable: ',                        & 
     &         TRIM(ADJUSTL(executable))
            RETURN
         END IF
         istat = INDEX(line,TRIM(ADJUSTL(executable)))
         IF (istat .ne. 0) THEN
            line_number = i
            EXIT
         END IF
      END DO
      
      IF (line_number .eq. 0) THEN
         WRITE(*,*) 'getmem: Could not find executable ',i
         mem = 0
         RETURN
      END IF
      
! Reread the line, and detect the memory size
 !     READ(line,*)  i, mem     
      READ(line,*)   mem     
      
!  Delete the odd_file_name file
      CLOSE(iou,iostat = istat,status='delete')

      IF (istat .ne. 0) THEN
         WRITE(*,*) 'Error closing file in getmem. iostat = ',istat
      END IF
      
      RETURN
      END SUBROUTINE getmem
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
      SUBROUTINE gauss_rand(n,x)
!  Subroutine to fill an array with gaussian random variables of
!  unit variance. Pass size of array, so that don't have to 
!  fuss with interfaces

      USE stel_kinds

      INTEGER, INTENT(in) :: n
      REAL(rprec), DIMENSION(n), INTENT(inout) :: x

      INTEGER :: np1o2, i, j
      REAL(rprec), DIMENSION(2) :: v
      REAL(rprec) :: rsq, fac

!  start of executable code
      np1o2 = (n + 1) / 2
      j = 0
      DO i = 1,np1o2
         j = j + 1
100      CALL RANDOM_NUMBER(v)
         v = 2. * v - 1.
         rsq = v(1) * v(1) + v(2) * v(2)
         IF ((rsq .gt. 1.) .or. (rsq .eq. 0.)) GO TO 100
         fac = SQRT(-2. * LOG(rsq) / rsq)
         x(j) = v(1) * fac
         j = j + 1
         IF (j .le. n) THEN
            x(j) = v(2) * fac
         ENDIF
      END DO
      RETURN
      END SUBROUTINE gauss_rand
!----------------------------------------------------------------------
!*******************************************************************************
!----------------------------------------------------------------------
!
!*******************************************************************************
!  Start of test subroutine 1
!*******************************************************************************

      SUBROUTINE test_1
!  JDH 07-05-2005. This version is used for some memory testing.
      
      USE stel_kinds

      IMPLICIT NONE

!----------------------------------------------------------------------
! L O C A L Variable Declarations
!----------------------------------------------------------------------
      REAL(rprec), DIMENSION(:), ALLOCATABLE, TARGET :: a, c 
      REAL(rprec), DIMENSION(:,:,:), POINTER :: b1 => null()
      REAL(rprec), DIMENSION(:,:,:), POINTER :: b2 => null()
      REAL(rprec), DIMENSION(:,:,:), POINTER :: b3 => null()
      REAL(rprec), DIMENSION(:,:,:), POINTER :: d1 => null()
      REAL(rprec), DIMENSION(:,:,:), POINTER :: d2 => null()
      REAL(rprec), DIMENSION(:,:,:), POINTER :: d3 => null()
      
      
      REAL(rprec), DIMENSION(:), POINTER :: p1 => null()
      REAL(rprec), DIMENSION(:), POINTER :: p2 => null()
      INTEGER :: i, istat, mb
      INTEGER :: mem, mem_old, expkb
      INTEGER :: i1, i2, i3, is
      
!-----------------------------------------------
!  Start of Executable Code
!-----------------------------------------------

      WRITE(*,*) 'Starting test subroutine test_1'
      
!      CALL system("ps -O vsz")
      
      CALL getmem(mem)
      WRITE(*,*) 'In test_1, result from getmem is', mem
      mem_old = mem
      
      WRITE(*,*) '#  i1  i2  i3  size  exp(kB)  mem-mem_old   mem'
      
!  One, single dimension, 101 size
      i1 = 3 * 32 * 101 * 101
      ALLOCATE(a(i1)) 
      expkb = i1 * 8 / 1024
      CALL getmem(mem)
      WRITE(*,*) 1, i1, 0, 0, i1, expkb, mem - mem_old, mem
      mem_old = mem

!  three, three dimension, 101 size
      i1 = 101
      i2 = 101
      i3 = 32
      ALLOCATE(b1(i1,i2,i3),b2(i1,i2,i3),b3(i1,i2,i3))
      is = 3 * i1 * i2 * i3
      expkb = is * 8 / 1024
      CALL getmem(mem)
      WRITE(*,*) 3, i1, i2, i3, is, expkb, mem - mem_old, mem
      mem_old = mem
      
!  One, single dimension, 128 size
      i1 = 3 * 32 * 128 * 128
      ALLOCATE(c(i1)) 
      expkb = i1 * 8 / 1024
      CALL getmem(mem)
      WRITE(*,*) 1, i1, 0, 0, i1, expkb, mem - mem_old, mem
      mem_old = mem

!  three, three dimension, 128 size
      i1 = 128
      i2 = 128
      i3 = 32
      ALLOCATE(d1(i1,i2,i3),d2(i1,i2,i3),d3(i1,i2,i3))
      is = 3 * i1 * i2 * i3
      expkb = is * 8 / 1024
      CALL getmem(mem)
      WRITE(*,*) 3, i1, i2, i3, is, expkb, mem - mem_old, mem
      mem_old = mem
      
      WRITE(*,*) ' now do deallocations'
      DEALLOCATE(a)
      CALL getmem(mem)
      WRITE(*,*) 'a', mem - mem_old, mem
      mem_old = mem
      
      DEALLOCATE(b1, b2, b3)
      CALL getmem(mem)
      WRITE(*,*) 'b', mem - mem_old, mem
      mem_old = mem

      DEALLOCATE(c)
      CALL getmem(mem)
      WRITE(*,*) 'c', mem - mem_old, mem
      mem_old = mem

      DEALLOCATE(d1, d2, d3)
      CALL getmem(mem)
      WRITE(*,*) 'd', mem - mem_old, mem
      mem_old = mem
 
 
      CALL getmem(mem)
      
      WRITE(*,*) 'In test_1, result from getmem is', mem         
 
      STOP
      
      END

!===================================================================================
      SUBROUTINE AJAX_FARADAY_INIT(integral_mode, in_dataset)
!--------------------------------------------------------------------------------
!
! Function:  Initializes the Faraday Rotation (aka interferometry/polarimetry) 
!             module by loading the current equilibrium into AJAX and calling the
!             top-level Faraday rotation subroutine.
!
! INPUTS:   integral_mode       ! integer toggle for integral of interest (FR, CM, interferometry, etc)
!                                       0 = INITIALIZATION_ONLY     ! no integration performed 
!                                       1 = F_ROTATION_ONLY         ! polarimetry integration only
!                                       2 = PHASE_SHIFT_ONLY        ! interferometry integration only
!                                       3 =  F_ROTATION_AND_PSHIFT
!                                       4 = FR_AND_COTTON_MOUTON 
!
!           in_dataset        !  specifies what type of data should be used to compute the
!                                B field info, etc
!                                       0 = V3FIT_IN_MEMORY     ! "model a" version of VMEC info
!                                       1 = GENERIC_WOUT_FILE   ! std VMEC output file
!                                       2 = HOULBERG_DATAFILE   ! test file that came w. TRACK package
!
!
! created Aug, 2005 by J. Shields
!
!----------------------------------------------------------------------------------
      USE stel_kinds
      USE stel_constants
      USE v3read_wout  ! loads an equilibrium from a wout file into AJAX
      USE faraday_mod  ! module integrating over Faraday rotation angle
      IMPLICIT NONE

!.........dummy variables...................................................................!
      INTEGER, INTENT(IN) :: integral_mode   ! specifies which integrals code should compute
      INTEGER, INTENT(IN) :: in_dataset   ! specifies where VMEC info is taken from


!.........local variables.......................................................................!
      INTEGER :: mysize = 99
      INTEGER :: counter = 0
      INTEGER :: myint, istat =1, i
      INTEGER :: input_dataset               ! specifies where VMEC info is taken from
      INTEGER :: input_dataset_default       ! default location to take VMEC info from
      INTEGER :: selected_integrals          ! local toggle to specify integrals of interest
      INTEGER :: selected_integrals_default
      INTEGER, PARAMETER :: V3FIT_IN_MEMORY = 0
      INTEGER, PARAMETER :: GENERIC_WOUT_FILE = 1
      INTEGER, PARAMETER :: HOULBERG_DATAFILE = 2
      INTEGER, PARAMETER :: INITIALIZATION_ONLY   = 0
      INTEGER, PARAMETER :: F_ROTATION_ONLY       = 1
      INTEGER, PARAMETER :: PHASE_SHIFT_ONLY      = 2
      INTEGER, PARAMETER :: F_ROTATION_AND_PSHIFT = 3
      INTEGER, PARAMETER :: FR_AND_COTTON_MOUTON  = 4

      CHARACTER(len=*), PARAMETER :: in_woutfile = 'wout_.nc'     
      CHARACTER(len=*), PARAMETER :: subname = 'AJAX_FARADAY_INIT: '     

      write(*,*) subname, 'passed dataset toggle = ', in_dataset


!      input_dataset_default = HOULBERG_DATAFILE
!      input_dataset_default = GENERIC_WOUT_FILE
      input_dataset_default = V3FIT_IN_MEMORY

      selected_integrals_default = F_ROTATION_ONLY

!.............check input integral_mode for bad input............!
      SELECT CASE(integral_mode)
      CASE(0)
         selected_integrals = selected_integrals_default
      CASE(1:4) 
        selected_integrals = integral_mode
      CASE DEFAULT 
        write(*,*) subname, "Unknown integral mode.  Using default FR"
        selected_integrals = selected_integrals_default
      END SELECT

!...........optional stuff apparently ONLY works in "explicit interfaces"!
!      if( PRESENT(in_dataset) ) then
!        input_dataset = in_dataset
!      else
!        input_dataset = input_dataset_default
!      end if

      SELECT CASE(in_dataset)
      CASE(0:2) 
        input_dataset = in_dataset
      CASE DEFAULT 
        write(*,*) subname, "Unknown dataset mode.  Using default."
        input_dataset = input_dataset_default
      END SELECT



!.......link "model_a" values with the AJAX RZLAM & MAGFLUX subroutines.  JS 3/7/05....!
      if ( input_dataset == HOULBERG_DATAFILE) then
         call JS_TRACK_AJAX_DR
      else if ( input_dataset == GENERIC_WOUT_FILE) then        
        call READEQ_FROM_WOUT(in_woutfile)
      else if ( input_dataset == V3FIT_IN_MEMORY) then        
        call JS_READEQ
      else
         write(*,*) subname, 'input_dataset LOGIC ERROR'
      end if

      if ( integral_mode .NE. INITIALIZATION_ONLY) then
        call  FARADAY_ROTATION(selected_integrals)
      end if

      RETURN
      END SUBROUTINE AJAX_FARADAY_INIT

!===================================================================================


!*******************************************************************************
! SECTION VI.	COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************

!  JDH 11-9-04 Added coding in init_signal, to check that all the magnetic 
!     diagnostic mrf's have the same r-z-phi grid. (If not, rerun V3RFUN)
!
!  JDH 07-05-2005. Added getmem2 subroutine. Totally frustrated with trying to
!     understand why memory behaves the way it does.
!
!  JMS 12-14-05.  Added supplementary modules for Interferometry/Polarimetry, as well
!     as a "driver" routine AJAX_FARADAY_INIT.  A logical toggle WANT_F_ROTATION has
!     been added to MAIN() to optionally enable/disable the Faraday Rotation routines.
!     In addition, I have added a temporary "kludge" routine initial_compute_JS_debug
!     that (optionally) replaces the default initial_compute in order to speed debugging.
!
!  JDH 06-23-2006. Added my_task character variable, SELECT CASE in MAIN.
!  
!  JDH 06-28-2006. Added i_work, r_work, and c_work. Added task vmec_monitor, similar
!     to previous initialize_compute 
!
!  JDH 09-11-06. Using recon_param s in task vmec_changep.
!
!  JDH 09-22-06. Added task gsq_on_grid, compute the g^2 function on a 
!     two dimensional grid of parameter values.
!
!  JDH 10-03-06. Add task gsq_jac_test.
!
!  JDH 10-06-06. Eliminate mmaw_r2s, added module mmaw.
!
!  JDH 2007-06-24
!    NB - Also look for change comments in individual subroutines
!    Changes to task_reconstruct_a1, printout modified, added simple
!    convergence test
!
!  JDH 2007-06-29
!    Cleaned up some continuations. Call to gsq_write_conf in reconstruct_a1.
!
!  2007-07-03 JDH
!    Added task_test - a generic test task, to be used by anybody
!    Added task_add_noise_a1 - adds noise to observed signal data, and reconstructs
!    Added subroutine gauss_rand - gaussian distribution, for task_add_noise_a1
!
! JMS 2007-08-07 moved VMEC_B_INIT call (needed for ipsl stuff) from reconstruct_a1 to
!                 module signal_model
!
!  JDH 2007-10-06
!    Modifications for reconstruction constraints.
!
!  JDH 2007-12-20
!    Increase dimensions of local arrays in init_main_nli
!
!  JDH 2008-01-19 - Added => null() to pointer declarations
!
!  JDH 2008-01-21
!    Eliminated task_magd_vmec_v3post. Too many compiler errors in mag_study_mod.f
!    Eliminated optional argument for task_reconstruct_a1 - not used
!    Added iou_rcl printouts to task_reconstruct_a1
!
!  JDH 2009-01-16
!    Cleaned up SPH's fix of command line parsing in subroutine init_commmand_line
!
!  JDH 2009-03-03
!    Change from signal_model to signal_mc
!
!  JDH 2009-03-05 - Deleted tasks 'ray_tracing' and 'rtrace_circle', along with
!     corresponding subroutines task_ray_tracing and task_rtrace_circle
!
!  JDH 2009-03-14 - Fix up some logic when only signal is geometric
!  JDH 2009-04-07 - geometric_filename now initialized in v3f_global
!
!  JDH 2009-04-13 - Changed parameterization of limiter_iso function (lif_).
!      NLI now can access full range of parameters.
!
!  JDH 2011-01-18 - Added easy specification of sdo sigmas and weights in
!     init_main_nli
!
!  JDH 2011-03-16
!    Improved writing of signals.
!
!  JDH 2011-06-07
!   Added task vmec_change_rp
!
!  JDH 2011-06-21
!    Slight code change to deal with preconditioning better, in reconstruct_a1
!
!  JDH 2011-08-01
!    Refactor sxrc -> sxrch
!
!  JDH 2011-08-26
!    Eliminated getmem2 - no longer used.
!    slight modification to getmem - prints out executable name when
!    returns from error.
!    Also fixed problem with full-path-names for executables.
!
!  JDH 2011-10-23
!    Missing lots of changes in the past 2 months. Added parameterized profiles, sxrch_dot
!
!  JDH 2011-10-23
!    Added pp_te and some thscte_dot file stuff
!
