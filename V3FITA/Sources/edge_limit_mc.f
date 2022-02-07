!*******************************************************************************
!  File edge_limit_mc.f
!  Contains module edge_limit_mc

!  Module is part of the V3FIT stellarator equilibrium reconstruction code.
!  The postfix _mc in the name indicates Model Compute
!  The module contains the subroutine edge_limit_mc_model_compute
!  Both an edge_limit_desc and a model are needed for this computation.

!*******************************************************************************
!  MODULE edge_limit_mc
!    (edge_limit - Model Compute Calculations)
! SECTION I.     VARIABLE DECLARATIONS
! SECTION II.    INTERFACE BLOCKS
! SECTION III.   MAIN COMPUTATIONAL SUBROUTINE
!*******************************************************************************
      MODULE edge_limit_mc

!*******************************************************************************
! SECTION I. VARIABLE DECLARATIONS
!*******************************************************************************

!-------------------------------------------------------------------------------
!  Type declarations, constants, utilities
!-------------------------------------------------------------------------------
      USE stel_kinds
      USE stel_constants
      USE v3_utilities

!-------------------------------------------------------------------------------
!  MDDC Derived Types
!-------------------------------------------------------------------------------
      USE edge_limit_T
      
!-------------------------------------------------------------------------------
!  Model Derived Types
!-------------------------------------------------------------------------------
      USE model_T

!-------------------------------------------------------------------------------
!  Module for limiter_iso function
!-------------------------------------------------------------------------------
      USE limiter_iso_T

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
!  Compute an edge_limit signal
!
!    Information comes from the edge_limit_desc and the model
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!    s_type = geometric
!      g_type = edge_limit
!        edge limit
!        signal_model_compute_edge_limit
!           % data(1) - limiter function - zero on limiter, negative inside plasma
!           % data(2) - R position of maximum limiter function
!           % data(3) - phi position of maximum limiter function
!           % data(4) - Z position of maximum limiter function
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

      SUBROUTINE edge_limit_mc_model_compute(a_el_desc,a_model,            &
     &   arr_data,arr_sigma)
!-------------------------------------------------------------------------------
!  Declarations
!-------------------------------------------------------------------------------
!  Declare Arguments 
      TYPE (edge_limit_desc), INTENT (inout), TARGET :: a_el_desc
      TYPE (model), INTENT (inout), TARGET :: a_model
      REAL(rprec), POINTER, DIMENSION(:) :: arr_data, arr_sigma
!  Note - convention is that arr_data and arr_sigma
!    1) are deallocated in this subroutine, and then
!    2) are allocated in this subroutine

!  Declare local variables
      TYPE (eq_state), POINTER             :: state => null()
      TYPE (edge_limit_desc), POINTER      :: my_el_desc => null()
      TYPE (limiter_iso), POINTER          :: my_iso => null()
      
      REAL(rprec) :: fval_max, fval
      REAL(rprec), DIMENSION(3) :: rpz_at_max, rpz_edge
      INTEGER :: n_data
      INTEGER :: ist, jumax, kvmax, j, k, ier1, istat1

      REAL(rprec), ALLOCATABLE, DIMENSION(:)  :: r_at_seq1, z_at_seq1
      REAL(rprec) :: phi
      INTEGER :: nu

      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'edge_limit_mc_model_compute: '

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------

!  Point state to the model equilibrium state
      state => a_model % eqstate
! 
!  Point to edge_limit
      my_el_desc => a_el_desc

      SELECT CASE (TRIM(ADJUSTL(my_el_desc % edge_limit_type)))
      CASE ('iso_fun')
         my_iso => my_el_desc % lim_iso
      
!  Allocate data and sigma
         n_data = 4
         IF (ASSOCIATED(arr_data)) THEN
            DEALLOCATE(arr_data,arr_sigma, stat=istat1)
            CALL assert_eq(0,istat1,sub_name // 'arr_data, s dealloc')
         ENDIF
         ALLOCATE(arr_data(n_data),arr_sigma(n_data),stat=istat1)
         CALL assert_eq(0,istat1,sub_name // 'arr_data, sigma alloc')

!  Determine number of poloidal points and allocate R and Z arrays
!  NOTE - access to state variables here.
         nu = MAX(my_iso % numin,state % fixp % ns_array(                      &
     &      state % varp % ns_index),2 * state % fixp % mpol + 6)
         ALLOCATE(r_at_seq1(nu),z_at_seq1(nu),stat=istat1)
         CALL assert_eq(0,istat1,sub_name // 'arr_data, r_at_seq1')

!----------------------------------------------------------------------
!-- Loop over toroidal planes                                        --
!----------------------------------------------------------------------
!  loop around the s = 1 surface
!  Evaluate limiter function for each point
!  take maximum, store R, Phi, and Z at location of maximum

         fval_max = -1.e10
         rpz_at_max = (/ zero, zero, zero /)
      
         jumax = nu
         kvmax = SIZE(my_iso % vgrid,1)

         TOR_PLANE: DO k = 1, kvmax
!  Get R and Z values at s=1, toroidal angle phi
            phi = my_iso % vgrid(k)
            CALL model_get_seq1_rz(a_model,phi,r_at_seq1,z_at_seq1)

            DO j = 1, jumax   ! loop over poloidal positions
               rpz_edge(1) = r_at_seq1(j)
               rpz_edge(2) = phi
               rpz_edge(3) = z_at_seq1(j)
               CALL limiter_iso_f_cyl(my_iso,rpz_edge,fval)
               IF (fval .gt. fval_max) THEN
                  fval_max = fval
                  rpz_at_max = rpz_edge
               ENDIF
            END DO

         END DO TOR_PLANE
         
         IF (my_el_desc % l_on_edge) THEN
            arr_data(1) = fval_max
         ELSE
            arr_data(1) = MAX(fval_max,zero)
         ENDIF
         arr_data(2:4) = rpz_at_max
         arr_sigma(1:4) = 1.E-10

      CASE ('polygon')
         CALL err_fatal(sub_name // 'This type NYI: ',                         &
     &      char=my_el_desc % edge_limit_type)

      CASE DEFAULT
         CALL err_fatal(sub_name // 'unrecognized edge_limit_type: ',          &
     &      char=my_el_desc % edge_limit_type)

      END SELECT ! Different coding depending on edge_limit_type

!  Deallocate the _at_seq1 arrays
      IF (ALLOCATED(r_at_seq1)) THEN
         DEALLOCATE(r_at_seq1,z_at_seq1,STAT=ier1)
         CALL assert_eq(0,ier1,sub_name // 'deallocate error')
      END IF
       
      RETURN
      END SUBROUTINE edge_limit_mc_model_compute

!*******************************************************************************
! SECTION IV.  SUPERCEDED CODE
!*******************************************************************************
!-------------------------------------------------------------------------------
!  JDH 2011-07-21
!    New version of edge_limit_mc_model_compute - more logical access to 
!    model data. Coding below is unused, and should be deleted soon.
!       SUBROUTINE edge_limit_mc_model_compute_old
!  JDH 2011-07-27
!    Eliminated edge_limit_mc_model_compute_old and edge_limit_mc_aux_rz_one_old
!-------------------------------------------------------------------------------

!*******************************************************************************
! SECTION X.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
!
!  JDH 2009-03-03. First version of module. Based on coding from previous
!    module signal_model
!
!  2009-03-14 JDH
!     Fixed so that aux1 is DEFINED, if need be.
!
!  2011-07-21 JDH
!     Added testing for model_get_seq1_rz
!
!  2011-07-22 JDH
!     Added edge_limit_mc_model_compute_new - test compiling, slow transition
!        11:00 AM
!     Moved edge_limit_mc_model_compute -> edge_limit_mc_model_compute_old
!           edge_limit_mc_aux_rz_one -> edge_limit_mc_aux_rz_one_old
!           edge_limit_mc_model_compute_new -> edge_limit_mc_model_compute
!
!  2011-07-27 JDH
!    Eliminate:
!      edge_limit_mc_model_compute_old
!      edge_limit_mc_aux_rz_one_old
!
      END MODULE edge_limit_mc
