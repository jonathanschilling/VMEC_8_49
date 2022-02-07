!*******************************************************************************
!  File thscte_mc.f
!  Contains module thscte_mc
!           SUBROUTINE thscte_mc_model_compute
!
!  Module is part of the V3FIT stellarator equilibrium reconstruction code.
!  The postfix _mc in the name indicates Model Compute
!  The module contains the subroutine thscte_mc_model_compute
!  Both an thscte_desc and a model are needed for this computation.

!*******************************************************************************
!  MODULE thscte_mc
!    (thscte - Thomson Scattering Te Model Compute Calculations)
! SECTION I.     VARIABLE DECLARATIONS
! SECTION II.    INTERFACE BLOCKS
! SECTION III.   MAIN COMPUTATIONAL SUBROUTINE
! SECTION IV.    _OLD (SUPERCEDED) CODE
! SECTION X.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
      MODULE thscte_mc

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
!  thscte Derived Type
!-------------------------------------------------------------------------------
      USE thscte_T
      
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
!  Compute an thscte signal
!
!    Information comes from the thscte_desc and the model
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!  Actual computation of the model signal is in this subroutine
!    s_type = diagnostic
!      d_type = thscte
!        Thomson scattering Te 
!        signal_model_compute_thscte 
!
!-------------------------------------------------------------------------------

      SUBROUTINE thscte_mc_model_compute(a_thscte,a_model,mod_signal,          &
     &   mod_sigma)

!-------------------------------------------------------------------------------
! ARGUMENTS
! a_thscte       - type thscte_desc - holds Thomson Scattering info
! a_model       - type model - holds eq_state and te profile information
! mod_signal    - output of the generated signal
! mod_sigma     - output sigma      
!------------------------------------------------------------------------------- 
      TYPE (thscte_desc), INTENT (inout)     :: a_thscte
      TYPE (model), INTENT (inout), TARGET   :: a_model
      REAL(rprec), POINTER, DIMENSION(:)     :: mod_signal, mod_sigma

!------------------------------------------------------------------------------- 
! Local Variables
!------------------------------------------------------------------------------- 
      INTEGER       :: istat1
      REAL(rprec), DIMENSION(3) :: r_cyl(1:3)=0.0
      REAL(rprec), DIMENSION(3) :: xcart
      REAL(rprec)               :: te

      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'thscte_mc_model_compute'

!-------------------------------------------------------------------------------
!  START OF EXECUTABLE CODE
!-------------------------------------------------------------------------------
         
!  Allocate data and sigma
      IF (ASSOCIATED(mod_signal)) THEN
         DEALLOCATE(mod_signal,mod_sigma, stat=istat1)
         CALL assert_eq(0,istat1,sub_name //                                   &
     &                     'mod_signal, sigma dealloc')
      ENDIF
      ALLOCATE(mod_signal(1),mod_sigma(1),stat=istat1)
      CALL assert_eq(0,istat1,sub_name // 'mod_signal, sigma alloc')

      xcart = a_thscte % xcart
      r_cyl(1) = SQRT(xcart(1) ** 2 + xcart(2) ** 2)
      r_cyl(2) = ATAN2(xcart(2),xcart(1))
      r_cyl(3) = xcart(3)
      te = model_get_te_xcyl(a_model,r_cyl)
      mod_signal(1) = te

! make up a sigma         
         mod_sigma(1)=0.
         
      END SUBROUTINE thscte_mc_model_compute

!*******************************************************************************
! SECTION X.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
!
!  JDH 2011-10-23
!    First version for thscte_mc. Based on sxrch_mc

      END MODULE thscte_mc
