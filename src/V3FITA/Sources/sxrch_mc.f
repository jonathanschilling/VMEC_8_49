!*******************************************************************************
!  File sxrch_mc.f
!  Contains module sxrch_mc
!           SUBROUTINE sxrch_mc_model_compute
!           FUNCTION emissivity_model
!
!  Module is part of the V3FIT stellarator equilibrium reconstruction code.
!  The postfix _mc in the name indicates Model Compute
!  The module contains the subroutine sxrch_mc_model_compute
!  Both an sxrch_desc and a model are needed for this computation.

!*******************************************************************************
!  MODULE sxrch_mc
!    (sxrch - Soft X-Ray Model Compute Calculations)
! SECTION I.     VARIABLE DECLARATIONS
! SECTION II.    INTERFACE BLOCKS
! SECTION III.   MAIN COMPUTATIONAL SUBROUTINE
! SECTION IV.    _OLD (SUPERCEDED) CODE
! SECTION X.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
      MODULE sxrch_mc

!*******************************************************************************
! SECTION I. VARIABLE DECLARATIONS
!*******************************************************************************

!-------------------------------------------------------------------------------
!  Type declarations, constants, utilities
!-------------------------------------------------------------------------------
      USE stel_kinds
      USE v3f_global
      USE vmec_utils
      USE read_wout_mod

!-------------------------------------------------------------------------------
!  sxrch Derived Type
!-------------------------------------------------------------------------------
      USE sxrch_T
      
!-------------------------------------------------------------------------------
!  Model Derived Types
!-------------------------------------------------------------------------------
      USE model_T

!------------------------------------------------------------------------------
      IMPLICIT NONE
      
!*******************************************************************************
! SECTION II.  INTERFACE BLOCKS
!*******************************************************************************

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
      CONTAINS
          
!*******************************************************************************
! SECTION III.  MAIN COMPUTATIONAL SUBROUTINE
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Compute an sxrch signal
!
!    Information comes from the sxrch_desc and the model
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!  Actual computation of the model signal is in this subroutine
!    s_type = diagnostic
!      d_type = sxrch
!        Soft X-Ray Chord Diagnostic
!        signal_model_compute_sxrch 
!
!-------------------------------------------------------------------------------

      SUBROUTINE sxrch_mc_model_compute(a_sxrch,a_model,mod_signal,            &
     &   mod_sigma)

!-------------------------------------------------------------------------------
! ARGUMENTS
! a_sxrch        - type sxrch_desc - holds sxr chord info
! a_model       - type model - holds eq_state and density model information
! mod_signal    - output of the generated signal
! mod_sigma     - output sigma      
!------------------------------------------------------------------------------- 
      TYPE (sxrch_desc), INTENT (inout)     :: a_sxrch
      TYPE (model), INTENT (inout), TARGET :: a_model
      REAL(rprec), POINTER, DIMENSION(:) :: mod_signal, mod_sigma

!------------------------------------------------------------------------------- 
! local variables
!------------------------------------------------------------------------------- 
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'sxrch_mc_model_compute'
     
!-------------------------------------------------------------------------------      
! cyl2flx declarations  - used for the subroutine cyl2flx call
!------------------------------------------------------------------------------- 
      REAL(rprec), DIMENSION(3) :: r_cyl(1:3)=0.0
      REAL(rprec), DIMENSION(3) :: xcart, dxcart
      INTEGER :: points
       
!------------------------------------------------------------------------------- 
! other declarations
!
!-------------------------------------------------------------------------------        
      INTEGER       :: istat1
      REAL(rprec)   :: dx               ! integration step size
      INTEGER       :: ipoints          ! number of integration points
      INTEGER       :: istep
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: sxrem_array

!-------------------------------------------------------------------------------
!  START OF EXECUTABLE CODE
!-------------------------------------------------------------------------------         
!  Allocate data and sigma
         IF (ASSOCIATED(mod_signal)) THEN
           DEALLOCATE(mod_signal,mod_sigma, stat=istat1)
           CALL assert_eq(0,istat1,sub_name //                                 &
     &                     'mod_signal, sigma dealloc')
         ENDIF
         ALLOCATE(mod_signal(1),mod_sigma(1),stat=istat1)
         CALL assert_eq(0,istat1,sub_name // 'mod_signal, sigma alloc')

!  Allocate space for sxrem_array
         ipoints=100
         ALLOCATE(sxrem_array(ipoints),stat=istat1)
         CALL assert_eq(0,istat1,sub_name // 'alloc sxrem_array')

         dxcart = (a_sxrch % xcart_f - a_sxrch % xcart_i)/ipoints
         dx = SQRT(dxcart(1) ** 2 + dxcart(2) ** 2 + dxcart(3) ** 2)
         DO istep=1,ipoints
            xcart = a_sxrch % xcart_i + istep * dxcart
            r_cyl(1) = SQRT(xcart(1) ** 2 + xcart(2) ** 2)
            r_cyl(2) = ATAN2(xcart(2),xcart(1))
            r_cyl(3) = xcart(3)
            sxrem_array(istep) =                                               &
     &         model_get_sxrem_xcyl(a_model,r_cyl) * dx
         ENDDO
         mod_signal(1)=SUM(sxrem_array)
         DEALLOCATE(sxrem_array,stat=istat1)
         CALL assert_eq(0,istat1,sub_name // 'dealloc sxrem_array')

! make up a sigma         
!         mod_sigma(1)=0.5   Changed to zero JDH 2011-10-17
         mod_sigma(1)=0.0
         
      END SUBROUTINE sxrch_mc_model_compute

!*******************************************************************************
! SECTION X.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
!
!  JDH 2009-03-03. First version of module. Based on coding from previous
!    module signal_model
!
!  GJH 2010-01-20. SXR version of coding based on module mddc_mc
!
!  JDH 2011-07-28
!    Refactoring - move emissivity function to the model. Test.
!    Moved to _old - sxrch_mc_model_compute_old and emissivity_model_old
!
!  JDH 2011-08-01
!    Refactor sxrc -> sxrch
!
!  JDH 2011-09-08
!    Remove old code. Use eq_get_flux instead of cyl2flx.
!
!  JDH 2011-10-17
!     Change to cartesian coordinates in sxrch
!
!  JDH 2011-10-25
!     Some refactoring, changed s_array to sxrem_array.

      END MODULE sxrch_mc
