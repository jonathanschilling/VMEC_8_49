!*******************************************************************************
!  File model_T.f
!  Contains module model_T
!  Defines derived-types: model

!*******************************************************************************
!  MODULE model_T
!    (Model Type Definition, for the V3FIT code)
! SECTION I.    Variable Declarations
! SECTION II.   Derived Type Declarations
! SECTION III.  Interface Blocks
! SECTION IV.   Construction Subroutine
! SECTION V.    Destruction Subroutines
! SECTION VI.   Assignment Subroutines
! SECTION VII.  Get Subroutines
! SECTION VIII. WRITE Subroutines

! SECTION XVI.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
      MODULE model_T

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
      USE bsc_T
      USE v3_utilities
      USE eq_T
!      USE density_T
      USE pprofile_T

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
      INTEGER(iprec), PARAMETER, PRIVATE :: type_len=20

!*******************************************************************************
! SECTION II. DERIVED-TYPE DECLARATIONS
! 1)   Model
!         model  
!     A model will contain an equilibrium state, and any other information 
!     that is needed to compute model signals. For example, if the experiment
!     has interferometry and polarimetry diagnsotics, then there will need
!     to be a model density profile.
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Declare type model
!         ---- Equilibrium State ----
!    eqstate         Type eq_state
!
!         ---- Parameterized Profiles (type pprofile)  -----
!    pp_ne          electron density
!    pp_sxrem       Soft x-ray emissivity
!    pp_te          Electron temperature
!
!         ----  Other Information ----
!  model_sxrem_type    Character variable to specify how the Soft x-ray 
!                      emissivity will be calculated by the model.
!      'none'             Not part of the model.
!      'pp_sxrem'         Calculated from the parameterized profile pp_sxrem
!      'pp_ne_vmec_p'     Calculated from the pp_ne and the VMEC pressure
!
!  model_te_type    Character variable to specify how the electron temperature 
!                   will be calculated by the model.
!      'none'             Not part of the model.
!      'pp_te'            Calculated from the parameterized profile pp_te
!      'pp_ne_vmec_p'     Calculated from the pp_ne and the VMEC pressure
!
!  ne_pp_unit        Scaling/conversion factor for the density. ne * ne_pp_unit
!                    is assumed to have units of m ** -3.
!
!  ne_min            Minimum electron density, m ** -3.
!
!  e_pressure_fraction    Electron pressure fraction of the total pressure
!                         Used when model_te_type = 'pp_ne_vmec_p'
!-------------------------------------------------------------------------------
      TYPE model
         TYPE (eq_state)          :: eqstate
         TYPE (pprofile)          :: pp_ne
         TYPE (pprofile)          :: pp_sxrem
         TYPE (pprofile)          :: pp_te
         CHARACTER(LEN=type_len)  :: model_sxrem_type
         CHARACTER(LEN=type_len)  :: model_te_type
         REAL(rprec)              :: ne_pp_unit
         REAL(rprec)              :: ne_min
         REAL(rprec)              :: e_pressure_fraction
      END TYPE model

!*******************************************************************************
! SECTION III. INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Assignment for structures
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!  Generic construct
!-------------------------------------------------------------------------------
 !-------------------------------------------------------------------------------
!  Generic destroy
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!  Generic write
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!  Interface block for testing goes here. 
!-------------------------------------------------------------------------------

      CONTAINS
!*******************************************************************************
! SECTION IV. CONSTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Construct a model
!-------------------------------------------------------------------------------
      SUBROUTINE model_construct(this,eqstate, pp_ne, pp_sxrem,                &
     &   pp_te, model_sxrem_type, model_te_type,ne_pp_unit,ne_min,             &
     &   e_pressure_fraction)

      IMPLICIT NONE

!  Declare Arguments 
      TYPE (model), INTENT(inout)              :: this
      TYPE (eq_state), INTENT(in), OPTIONAL    :: eqstate
      TYPE (pprofile), INTENT(in), OPTIONAL    :: pp_ne
      TYPE (pprofile), INTENT(in), OPTIONAL    :: pp_sxrem
      TYPE (pprofile), INTENT(in), OPTIONAL    :: pp_te
      CHARACTER(LEN=20), INTENT(in), OPTIONAL  :: model_sxrem_type
      CHARACTER(LEN=20), INTENT(in), OPTIONAL  :: model_te_type
      REAL(rprec), INTENT(in), OPTIONAL        :: ne_pp_unit
      REAL(rprec), INTENT(in), OPTIONAL        :: ne_min
      REAL(rprec), INTENT(in), OPTIONAL        :: e_pressure_fraction
      
!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'model_construct: '

!  Start of executable code

!  Copy the OPTIONAL arguments
      IF (PRESENT(eqstate)) THEN
         this % eqstate = eqstate
      ENDIF

      IF (PRESENT(pp_ne)) THEN
         this % pp_ne = pp_ne
      ENDIF

      IF (PRESENT(pp_sxrem)) THEN
         this % pp_sxrem = pp_sxrem
      ENDIF

      IF (PRESENT(pp_te)) THEN
         this % pp_te = pp_te
      ENDIF

      IF (PRESENT(model_sxrem_type)) THEN
         this % model_sxrem_type = model_sxrem_type
      ENDIF

      IF (PRESENT(model_te_type)) THEN
         this % model_te_type = model_te_type
      ENDIF
 
      IF (PRESENT(ne_pp_unit)) THEN
         this % ne_pp_unit = ne_pp_unit
      ENDIF
     
      IF (PRESENT(ne_min)) THEN
         this % ne_min = ne_min
      ENDIF
     
      IF (PRESENT(e_pressure_fraction)) THEN
         this % e_pressure_fraction = e_pressure_fraction
      ENDIF
     
      END SUBROUTINE model_construct

      
!*******************************************************************************
! SECTION V. DESTRUCTION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Destroy a model
!-------------------------------------------------------------------------------
      SUBROUTINE model_destroy(this)
      IMPLICIT NONE

!  Declare Arguments 
      TYPE (model), INTENT(inout) :: this

!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'model_destroy: '

!  Start of executable code

!  As of 12-14-2004, there isn't an eq_destroy for a state.
!      CALL eq_destroy(this % eqstate)

!.........destroy the parameterized profiles
      CALL pprofile_destroy( this % pp_ne)
      CALL pprofile_destroy( this % pp_sxrem)
      CALL pprofile_destroy( this % pp_te)
      
      this % model_sxrem_type = 'none'
      this % model_te_type = 'none'

      this % ne_pp_unit = zero
      this % ne_min = zero
      this % e_pressure_fraction = zero

      END SUBROUTINE model_destroy

!*******************************************************************************
! SECTION VI. ASSIGNMENT SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Assignment for model
!-------------------------------------------------------------------------------
      SUBROUTINE model_assign(left,right)

      IMPLICIT NONE
      
!  Declare Arguments 
      TYPE (model), INTENT (inout) :: left
      TYPE (model), INTENT (in) :: right
      
!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'model_assign: '
         
!  Start of executable code

      left % eqstate = right % eqstate
      left % pp_ne = right % pp_ne
      left % pp_sxrem = right % pp_sxrem
      left % pp_te = right % pp_te
      left % model_sxrem_type = right % model_sxrem_type
      left % model_te_type = right % model_te_type
      left % ne_pp_unit = right % ne_pp_unit
      left % ne_min = right % ne_min
      left % e_pressure_fraction = right % e_pressure_fraction
         
      END SUBROUTINE model_assign

!*******************************************************************************
! SECTION VII. GET SUBROUTINES
!*******************************************************************************
! Information needed about anything in the model should be accessed through 
! these subroutines. (That's the goal. JDH 2011-07-14)
!
!  Need to get the magnetic diagnostics working through here.
! 
!  For profiles, the _xcyl functions call the _s functions.
!-------------------------------------------------------------------------------
!  Plasma surface information (Pass through to eq)
!-------------------------------------------------------------------------------
      SUBROUTINE model_get_seq1_rz(this,phi,r_at_seq1,z_at_seq1)
!  Subroutine to return the R and Z positions of points on the
!  plasma surface (s = 1 surface), on a constant_phi plane.
      IMPLICIT NONE
      
!  Declare Arguments 
      TYPE (model), INTENT (inout)             :: this
      REAL(rprec), INTENT (in)                 :: phi
      REAL(rprec), DIMENSION(:), INTENT(inout) :: r_at_seq1, z_at_seq1
      
!  Start of executable code. Pass through to the eq.
      CALL eq_get_seq1_rz(this % eqstate,phi,r_at_seq1,z_at_seq1)
         
      END SUBROUTINE model_get_seq1_rz
!-------------------------------------------------------------------------------
!  Electron Density ne
!
!     Given Cylindrical Coordinates XCYL
!-------------------------------------------------------------------------------
      FUNCTION model_get_ne_xcyl(this,xcyl)
!  FUNCTION to return the electron density, given a model and 
!  Cylindrical coordinates xcyl(3) [ = R, Phi, Z]
      IMPLICIT NONE
      
!  Declare Arguments 
      REAL(rprec) :: model_get_ne_xcyl
      TYPE (model), INTENT (inout)             :: this
      REAL (rprec), INTENT (in)                :: xcyl(3)

!  Declare local variables
      REAL(rprec) :: c_flux(3)               ! VMEC coordinates (s,u,v)
      REAL(rprec) :: s                       ! VMEC radial coordinate
      REAL(rprec) :: ne                      ! electron density
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'model_get_ne_s: '
               
!  Start of executable code.

!  Get the radial coordinate s
      c_flux = (/zero, zero, zero /)
      CALL eq_get_flux(this % eqstate,xcyl,c_flux)
      s = c_flux(1)

      ne = model_get_ne_s(this,s)       
      model_get_ne_xcyl = ne
      RETURN
         
      END FUNCTION model_get_ne_xcyl
!-------------------------------------------------------------------------------
!  Electron Density ne
!
!     Given VMEC radial coordinate S
!-------------------------------------------------------------------------------
      FUNCTION model_get_ne_s(this,s_arg)
!  FUNCTION to return the electron density in m ** -3, given a model 
!  and a value of the radial coordinate s.
!  Note the conversion factor ne_pp_unit.
      IMPLICIT NONE
      
!  Declare Arguments 
      REAL(rprec) :: ne                      ! electron density
      TYPE (model), INTENT (inout)             :: this
      REAL (rprec), INTENT (in)                :: s_arg

!  Declare local variables
      REAL(rprec) :: model_get_ne_s
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'model_get_ne_s: '
               
!  Start of executable code.

      model_get_ne_s = MAX(this % ne_min,this % ne_pp_unit *                   &
     &   pprofile_evaluate(this % pp_ne,s_arg))
      RETURN
         
      END FUNCTION model_get_ne_s

!-------------------------------------------------------------------------------
!  Soft x-ray emissivity
!
!     Given Cylindrical Coordinates XCYL
!-------------------------------------------------------------------------------
      FUNCTION model_get_sxrem_xcyl(this,xcyl)
!  FUNCTION to return the soft x-ray emissivity, given a model and 
!  Cylindrical coordinates xcyl(3) [ = R, Phi, Z]
      IMPLICIT NONE
      
!  Declare Arguments 
      REAL(rprec) :: model_get_sxrem_xcyl
      TYPE (model), INTENT (inout)             :: this
      REAL (rprec), INTENT (in)                :: xcyl(3)

!  Declare local variables
      REAL(rprec) :: c_flux(3)               ! VMEC coordinates (s,u,v)
      REAL(rprec) :: s                       ! VMEC radial coordinate
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'model_get_sxrem_xcyl: '
               
!  Start of executable code.

!  Get the radial coordinate s
      c_flux = (/zero, zero, zero /)
      CALL eq_get_flux(this % eqstate,xcyl,c_flux)
      s = c_flux(1)

      model_get_sxrem_xcyl = model_get_sxrem_s(this,s)

      RETURN         
      END FUNCTION model_get_sxrem_xcyl

!-------------------------------------------------------------------------------
!  Soft x-ray emissivity
!
!     Given VMEC radial coordinate S
!-------------------------------------------------------------------------------
      FUNCTION model_get_sxrem_s(this,s_arg)
!  FUNCTION to return the soft x-ray emissivity, given a model and 
!  and a value of the radial coordinate s.
      IMPLICIT NONE
      
!  Declare Arguments 
      REAL(rprec) :: model_get_sxrem_s
      TYPE (model), INTENT (inout)             :: this
      REAL (rprec), INTENT (in)                :: s_arg

!  Declare local variables
      REAL(rprec) :: ne 
      REAL(rprec) :: pressure
      REAL(rprec) :: emissivity
      REAL(rprec) :: ne_norm = 1.E18         ! Normalizing constant, per Hartwell routine
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'model_get_ne_s: '
               
!  Start of executable code.
      SELECT CASE(TRIM(this % model_sxrem_type))
      
      CASE('pp_sxrem')
         emissivity = pprofile_evaluate(this % pp_sxrem,s_arg)     

      CASE('pp_ne_vmec_p')
         IF (s_arg .le. one) THEN
            ne = model_get_ne_s(this,s_arg)
            pressure = MAX(0.,eq_get_pressure_s(this % eqstate,s_arg))
            emissivity = (ne / ne_norm) ** (3/2) * SQRT(pressure)   ! per comments         
         ELSE
            emissivity = zero
         ENDIF
         
      CASE ('none')
         emissivity = -2.E49_rprec
         
      CASE DEFAULT
         WRITE(*,*) ' IN model_get_sxrem_s. model_sxrem_type is ',             &
     &      TRIM(this % model_sxrem_type)
         WRITE(*,*) '  Needs to be pp_sxrem or pp_ne_vmec_p'
         STOP

      END SELECT
       
      model_get_sxrem_s = emissivity

      RETURN
      END FUNCTION model_get_sxrem_s

!-------------------------------------------------------------------------------
!  Line integral for Interferometry - ne
!
!     Given Cylindrical Coordinates XCYL
!-------------------------------------------------------------------------------
      FUNCTION model_get_ip_i_xcyl(this,xcyl)
!  FUNCTION to return the interferometry integrand, given a model and 
!  Cylindrical coordinates xcyl(3) [ = R, Phi, Z]
      IMPLICIT NONE
      
!  Declare Arguments 
      REAL(rprec) :: model_get_ip_i_xcyl
      TYPE (model), INTENT (inout)             :: this
      REAL (rprec), INTENT (in)                :: xcyl(3)

!  Declare local variables
      REAL(rprec) :: c_flux(3)               ! VMEC coordinates (s,u,v)
      REAL(rprec) :: s                       ! VMEC radial coordinate
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'model_get_ip_i_xcyl: '
               
!  Start of executable code.

!  Get the radial coordinate s
      c_flux = (/zero, zero, zero /)
      CALL eq_get_flux(this % eqstate,xcyl,c_flux)
      s = c_flux(1)

      model_get_ip_i_xcyl = model_get_ip_i_s(this,s)

      RETURN         
      END FUNCTION model_get_ip_i_xcyl

!-------------------------------------------------------------------------------
!  Line integral for Interferometry - ne
!
!     Given VMEC radial coordinate S
!-------------------------------------------------------------------------------
      FUNCTION model_get_ip_i_s(this,s_arg)
!  FUNCTION to return the soft x-ray emissivity, given a model and 
!  and a value of the radial coordinate s.
      IMPLICIT NONE
      
!  Declare Arguments 
      REAL(rprec) :: model_get_ip_i_s
      TYPE (model), INTENT (inout)             :: this
      REAL (rprec), INTENT (in)                :: s_arg

!  Declare local variables
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'model_get_ip_i_s: '
               
!  Start of executable code.
       
      model_get_ip_i_s = model_get_ne_s(this,s_arg)

      RETURN
      END FUNCTION model_get_ip_i_s
!-------------------------------------------------------------------------------
!  Electron Temperature 
!
!     Given Cylindrical Coordinates XCYL
!-------------------------------------------------------------------------------
      FUNCTION model_get_te_xcyl(this,xcyl)
!  FUNCTION to return the electron temperature in eV, given a model and 
!  Cylindrical coordinates xcyl(3) [ = R, Phi, Z]
      IMPLICIT NONE
      
!  Declare Arguments 
      REAL(rprec) :: model_get_te_xcyl
      TYPE (model), INTENT (inout)             :: this
      REAL (rprec), INTENT (in)                :: xcyl(3)

!  Declare local variables
      REAL(rprec) :: c_flux(3)               ! VMEC coordinates (s,u,v)
      REAL(rprec) :: s                       ! VMEC radial coordinate
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'model_get_te_xcyl: '
               
!  Start of executable code.

!  Get the radial coordinate s
      c_flux = (/zero, zero, zero /)
      CALL eq_get_flux(this % eqstate,xcyl,c_flux)
      s = c_flux(1)
       
      model_get_te_xcyl = model_get_te_s(this,s)
      RETURN
         
      END FUNCTION model_get_te_xcyl
!-------------------------------------------------------------------------------
!  Electron Temperature 
!
!     Given VMEC radial coordinate S
!-------------------------------------------------------------------------------
      FUNCTION model_get_te_s(this,s_arg)
!  FUNCTION to return the electron temperature in eV, given a model and 
!  and a value of the radial coordinate s.
      IMPLICIT NONE
      
!  Declare Arguments 
      REAL(rprec) :: model_get_te_s
      TYPE (model), INTENT (inout)             :: this
      REAL (rprec), INTENT (in)                :: s_arg

!  Declare Parameters
      REAL(rprec), PARAMETER :: eV_per_Joule = one / 1.602e-19

!  Declare local variables
      REAL(rprec) :: ne                      ! electron density (m ** -3)
      REAL(rprec) :: te                      ! electron temperature (eV)
      REAL(rprec) :: pressure                ! pressure
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'model_get_te_s: '
               
!  Start of executable code.

!  Get the radial coordinate s
      SELECT CASE(TRIM(this % model_te_type))
      
      CASE('pp_te')
         te = pprofile_evaluate(this % pp_te,s_arg)     

      CASE('pp_ne_vmec_p')
         ne = model_get_ne_s(this,s_arg)
         pressure = MAX(0.,eq_get_pressure_s(this % eqstate,s_arg))
         te = eV_per_Joule * this % e_pressure_fraction * pressure / ne         
         
      CASE ('none')
         te = -2.E49_rprec
         
      CASE DEFAULT
         WRITE(*,*) ' IN model_get_te_s. model_te_type is ',                   &
     &      TRIM(this % model_te_type)
         WRITE(*,*) '  Needs to be pp_te or pp_ne_vmec_p'
         STOP

      END SELECT
       
      model_get_te_s = te
      RETURN
         
      END FUNCTION model_get_te_s
!*******************************************************************************
! SECTION VIII.  WRITE SUBROUTINES
!*******************************************************************************
      SUBROUTINE model_write_profiles(a_model,identifier,unit)
      IMPLICIT NONE

!  Declare Arguments
      TYPE (model), INTENT (inout) :: a_model
      CHARACTER (len=*), INTENT(in), OPTIONAL :: identifier
      INTEGER, INTENT(in), OPTIONAL :: unit
!  identifier   character variable, also written out
!  unit         I/O unit number to write to

!  Declare local variables and constants
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'model_write_profiles: '
      INTEGER :: iv_default = 1
      INTEGER :: iv
      INTEGER :: iou_default = 6
      INTEGER :: iou
      INTEGER :: i, nrow 
      CHARACTER (len=60) :: id
      REAL(rprec) :: s, p, ne, te, sxrem

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

!  Write Header
      WRITE(iou,900) id

      CALL pprofile_write(a_model % pp_ne,sub_name // 'ne',iou)
      CALL pprofile_write(a_model % pp_sxrem,sub_name // 'sxrem',iou)
      CALL pprofile_write(a_model % pp_te,sub_name // 'te',iou)

      WRITE(iou,1000) 
      WRITE(iou,1010) 
      WRITE(iou,1020) TRIM(a_model % model_te_type),                           &
     &   TRIM(a_model % model_sxrem_type)

      nrow = a_model % eqstate % fixp % ns_array(                              &
     &   a_model % eqstate % varp % ns_index) + 2
      DO i = 1,nrow + 2
         s = (i - one) / nrow
         p = eq_get_pressure_s(a_model % eqstate,s)
         ne = model_get_ne_s(a_model,s)
         te = model_get_te_s(a_model,s)
         sxrem = model_get_sxrem_s(a_model,s)
         WRITE(iou,1100) s, p, ne, te, sxrem
      END DO

      WRITE(iou,1200) id
      
      RETURN

900   FORMAT(/' *** Model Write Profiles. id = ',a)
1000  FORMAT(/t5,'S',t16,'PRESSURE',t30,'e DENSITY',t46,                       &
     &   'e TEMPERATURE',t63,'SXR EMISSIVITY')
1010  FORMAT(t18,'(Pascal)',t32,'(m-3)',t48,'(eV)',t65,'( - )')
1020  FORMAT(t18,'VMEC',t32,'pp_ne',t48,a12,t65,a12/80('-'))
1100  FORMAT(2x,es9.2,t16,es10.3,t30,es10.3,t46,es10.3,t63,es10.3)
1200  FORMAT(/' *** END Model Write Profiles. id = ',a)

      END SUBROUTINE model_write_profiles
!*******************************************************************************
! SECTION XVI.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
!  
!  JDH 2011-07-14
!    Added model_get_seq1_rz
!
!  JDH 12-14-04. 
!     First version, modified signal_T
!  JMS 7-24-07.
!     added density to the model
!
!  JDH 2011-07-14
!    Add model_get_seq1_rz
!
!  JDH 2011-07-28
!    Add model_get_sxrem_s, model_get_ne_s, model_get_pressure_s
!
!  JDH 2011-08-07
!    Correct typo in emissivity, * SQRT(pressure) instead of + SQRT(pressure)
!
!  JDH 2011-09-06
!    model_get_szr_em_s: add MAX(0.,-) for density and pressure
!
!  2011-09-08 JDH
!    Replaced model_get_sxrem_s with _xcyl. Eliminated
!    model_get_pressure_s.
!
!  2011-10-19 JDH
!    Convert from e_density to pprofile pp_ne.
!
!  2011-10-23 JDH
!    Added model_sxrem_type and pp_sxrem
!
!  2011-10-23 JDH
!    Added model_te_type and pp_te, FUNCTION model_get_te_xcyl
!
!  2011-10-25 JDH
!    Added ne_pp_unit, so that ne returned by model_get is in m ** -3.
!    Added e_pressure_fraction. Rearranged the get functions.
!    Added model_write_profiles.

      END MODULE model_T
      