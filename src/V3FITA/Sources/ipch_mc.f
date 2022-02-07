!*******************************************************************************
!  File ipch_mc.f
!  Contains module ipch_mc
!           SUBROUTINE ipch_mc_model_compute
!           FUNCTION emissivity_model
!
!  Module is part of the V3FIT stellarator equilibrium reconstruction code.
!  The postfix _mc in the name indicates Model Compute
!  The module contains the subroutine ipch_mc_model_compute
!  Both an ipch_desc and a model are needed for this computation.

!*******************************************************************************
!  MODULE ipch_mc
!    (ipch - Interferometry-Polarimetry Model Compute Calculations)
! SECTION I.     VARIABLE DECLARATIONS
! SECTION II.    INTERFACE BLOCKS
! SECTION III.   MAIN COMPUTATIONAL SUBROUTINE
! SECTION IV.    _OLD (SUPERCEDED) CODE
! SECTION X.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
      MODULE ipch_mc

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
!  ipch Derived Type
!-------------------------------------------------------------------------------
      USE ipch_T
      
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
!  Compute an ipch signal
!
!    Information comes from the ipch_desc and the model
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!  Actual computation of the model signal is in this subroutine
!    s_type = diagnostic
!      d_type = ipch
!        Interferometry-Polarimetry Diagnostic
!        signal_model_compute_ipch 
!
!-------------------------------------------------------------------------------

      SUBROUTINE ipch_mc_model_compute(a_ipch,a_model,mod_signal,              &
     &   mod_sigma)

!-------------------------------------------------------------------------------
! ARGUMENTS
! a_ipch        - type ipch_desc - holds ip chord info
! a_model       - type model - holds eq_state and density model information
! mod_signal    - output of the generated signal
! mod_sigma     - output sigma      
!------------------------------------------------------------------------------- 
      TYPE (ipch_desc), INTENT (inout)     :: a_ipch
      TYPE (model), INTENT (inout), TARGET :: a_model
      REAL(rprec), POINTER, DIMENSION(:) :: mod_signal, mod_sigma

!------------------------------------------------------------------------------- 
! local variables
!------------------------------------------------------------------------------- 
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'ipch_mc_model_compute'
     
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
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: ipem_array

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

!  Allocate space for ipem_array
         ipoints=200
         ALLOCATE(ipem_array(ipoints),stat=istat1)
         CALL assert_eq(0,istat1,sub_name // 'alloc ipem_array')

         dxcart = (a_ipch % xcart_f - a_ipch % xcart_i)/ipoints
         dx = SQRT(dxcart(1) ** 2 + dxcart(2) ** 2 + dxcart(3) ** 2)
         DO istep=1,ipoints
            xcart = a_ipch % xcart_i + istep * dxcart
            r_cyl(1) = SQRT(xcart(1) ** 2 + xcart(2) ** 2)
            r_cyl(2) = ATAN2(xcart(2),xcart(1))
            r_cyl(3) = xcart(3)
            SELECT CASE(a_ipch % ip_type)
            
            CASE DEFAULT
               ipem_array(istep) =                                             &
     &            model_get_ip_i_xcyl(a_model,r_cyl) * dx
            END SELECT
         ENDDO
         mod_signal(1)=SUM(ipem_array)
         DEALLOCATE(ipem_array,stat=istat1)
         CALL assert_eq(0,istat1,sub_name // 'dealloc ipem_array')

! make up a sigma         
!         mod_sigma(1)=0.5   Changed to zero JDH 2011-10-17
         mod_sigma(1)=0.0
         
      END SUBROUTINE ipch_mc_model_compute

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
!
!  JDH 2012-03-15
!    First version for ip - based on sxr

      END MODULE ipch_mc
