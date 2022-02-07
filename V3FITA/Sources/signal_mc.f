!*******************************************************************************
!  File signal_mc.f
!  Contains module signal_mc

!  Module is part of the V3FIT stellarator equilibrium reconstruction code.
!  The postfix _mc in the name indicates Model Compute
!  The module contains the subroutine signal_mc_model_compute
!  Both an signal_desc and a model are needed for this computation.

!*******************************************************************************
!  MODULE signal_mc
!    (Signal - Model Calculations)
! SECTION I.     VARIABLE DECLARATIONS
! SECTION II.    INTERFACE BLOCKS
! SECTION III.   MAIN COMPUTATIONAL SUBROUTINE
!*******************************************************************************
      MODULE signal_mc

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
!  Signal Derived Types
!-------------------------------------------------------------------------------
      USE signal_T
      
!-------------------------------------------------------------------------------
!  Model Derived Types
!-------------------------------------------------------------------------------
      USE model_T

!-------------------------------------------------------------------------------
!  Extra Modules, for IPSL stuff
!-------------------------------------------------------------------------------
      USE read_wout_mod, ONLY: read_wout_file      
      USE b_transform_vmec, ONLY: VMEC_B_INIT

!-------------------------------------------------------------------------------
!  Other _mc modules
!-------------------------------------------------------------------------------
      USE mddc_mc
      USE ipsl_mc
      USE sxrch_mc           !GJH 2010-01-20
      USE ipch_mc            !JDH 2012-03-15
      USE thscte_mc
      USE edge_limit_mc
!  Note that there is no coosig_mc module. Coosig's are computed here.
      
      IMPLICIT NONE
      
!*******************************************************************************
! SECTION II.  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!  signal_mc_model_compute
!-------------------------------------------------------------------------------
      INTERFACE signal_mc_model_compute
         MODULE PROCEDURE signal_mc_model_compute,                             &
     &      signal_mc_model_compute_a
      END INTERFACE

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
      CONTAINS
          
!*******************************************************************************
! SECTION III.  MAIN COMPUTATIONAL SUBROUTINE
!*******************************************************************************
!-------------------------------------------------------------------------------
!  Compute a model signal
!
!    Information comes from the signal and the equilibrium
!-------------------------------------------------------------------------------
!      SUBROUTINE signal_mc_model_compute(sdata,sdesc,a_model)
! JDH 2011-07-07. Comment out above - get sdesc from sdata pointer
      SUBROUTINE signal_mc_model_compute(sdata,a_model)

!-------------------------------------------------------------------------------
!  Declarations
!-------------------------------------------------------------------------------
!  Declare Arguments 
      TYPE (signal_data), INTENT (inout) :: sdata
      TYPE (model), INTENT (inout) :: a_model

!  Declare local variables
      TYPE (signal_desc), POINTER :: sdesc
      INTEGER :: ierr
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'signal_mc_model_compute: '

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!  Things common to all s_type
!-------------------------------------------------------------------------------

!  Deallocate data and sigma
      IF (ASSOCIATED(sdata % data)) THEN
         DEALLOCATE(sdata % data,sdata % sigma)
      ENDIF

!  Set signal_desc pointer
      sdesc => sdata % desc

!  Set sd_type
      sdata % sd_type = 'model'

!-------------------------------------------------------------------------------
!  Different coding, depending on s_type
!-------------------------------------------------------------------------------
      SELECT CASE (TRIM(ADJUSTL(sdesc % s_type)))
      CASE ('diagnostic')
         SELECT CASE (TRIM(ADJUSTL(sdesc % diag % d_type)))
         CASE ('mddc')
             CALL mddc_mc_model_compute(sdesc % diag % mddc,a_model,           &
     &               sdata % data,sdata % sigma)
         CASE ('ipsl')
             CALL read_wout_file('wout_.nc', ierr)
             if (ierr .ne. 0 ) write(*,*) sub_name // 'read_wout error'
             call VMEC_B_INIT(1)
             CALL ipsl_mc_model_compute(sdesc % diag % ipsl,a_model,           &
     &               sdata % data,sdata % sigma)
         CASE ('sxrch')
             CALL sxrch_mc_model_compute(sdesc % diag % sxrch,a_model,         &
     &               sdata % data, sdata % sigma)
         CASE ('ipch')
             CALL ipch_mc_model_compute(sdesc % diag % ipch,a_model,           &
     &               sdata % data, sdata % sigma)
         CASE ('thscte')
             CALL thscte_mc_model_compute(sdesc % diag % thscte,               &
     &               a_model, sdata % data, sdata % sigma)
         CASE DEFAULT
            CALL err_fatal(sub_name // 'unrecognized d_type: ',                &
     &         sdesc % diag % d_type)
         END SELECT ! Select Case on d_type
      
      CASE ('geometric')
         SELECT CASE (TRIM(ADJUSTL(sdesc % geom % g_type)))
         CASE ('edge_limit')
             CALL edge_limit_mc_model_compute(sdesc % geom % el_desc,          &
     &               a_model,sdata % data,sdata % sigma)
         CASE DEFAULT
            CALL err_fatal(sub_name // 'unrecognized g_type: ',                &
     &         sdesc % geom % g_type)
         END SELECT ! Select Case on g_type
      
      CASE ('coosig')
!  Combination Of Other SIGnals
!    A coosig needs to know the other signals, which are specified by indices
!    in an array of signal_data. This subroutine only knows about a single 
!    signal, so the coosig can't be computed. Do nothing.
         
      CASE DEFAULT
         CALL err_fatal(sub_name // 'unrecognized s_type: ',                   &
     &      char=sdesc % s_type)
      END SELECT ! Different coding depending on d_type

      RETURN
      END SUBROUTINE signal_mc_model_compute

!-------------------------------------------------------------------------------
!  Compute an array of model signals
!
!    Information comes from the signals and the equilibrium
!-------------------------------------------------------------------------------
      SUBROUTINE signal_mc_model_compute_a(sdata,a_model)

!-------------------------------------------------------------------------------
!  Declarations
!-------------------------------------------------------------------------------
!  Declare Arguments 
      TYPE (signal_data), DIMENSION(:), INTENT (inout) :: sdata
      TYPE (model), INTENT (inout) :: a_model

!  Declare local variables
      TYPE (signal_desc), POINTER :: sdesc
      TYPE (coosig_desc), POINTER :: this_coosig
      INTEGER :: ierr, n_sdata, n_data, isig, ic, ind, nc
      REAL(rprec) :: tdata, tsigma
      CHARACTER(len=*), PARAMETER :: sub_name =                                &
     &  'signal_mc_model_compute_a: '

!-------------------------------------------------------------------------------
!  Start of executable code
!-------------------------------------------------------------------------------

!  Get size of sdata array
      n_sdata = SIZE(sdata)

!  First loop over arrays
      DO isig = 1,n_sdata
         CALL signal_mc_model_compute(sdata(isig),a_model)
      END DO

!  Second loop over arrays, to get the coosig (Combinations Of Other SIGnals)
!  signals.
      DO isig = 1,n_sdata
         sdesc => sdata(isig) % desc
         IF (TRIM(ADJUSTL(sdesc % s_type)) .ne. 'coosig') CYCLE

!  Here we have a coosig. Take care of the combination here, as this subroutine
!  knows the other signals

!  Take care of allocation of data and sigma
         IF (ASSOCIATED(sdata(isig) % data)) THEN
            n_data = SIZE(sdata(isig) % data)
            IF (n_data .ne. 1) THEN
               DEALLOCATE(sdata(isig) % data,sdata(isig) % sigma)
               ALLOCATE(sdata(isig) % data(1),sdata(isig) % sigma(1))
            ENDIF
         ELSE
            ALLOCATE(sdata(isig) % data(1),sdata(isig) % sigma(1))
         ENDIF

!  Set coosig pointer
         this_coosig => sdesc % coosig

!  Set sd_type
         sdata(isig) % sd_type = 'model'

!  Form the combination. Put the select case on the combination type within the
!  loop over the signal to combine (index ic)
         tdata = zero
         IF (TRIM(ADJUSTL(this_coosig % comb_type)) .eq. 'max')                &
     &      tdata = -1.e99_rprec
         IF (TRIM(ADJUSTL(this_coosig % comb_type)) .eq. 'min')                &
     &      tdata = 1.e99_rprec
         tsigma = zero
         nc = SIZE(this_coosig % index)
         DO ic = 1, nc
            ind = this_coosig % index(ic)
            CALL assert(ind .ge. 1,ind .le. n_sdata, ind .ne. isig,                 &
     &         sub_name // 'combining signals, bad index')
            SELECT CASE(TRIM(ADJUSTL(this_coosig % comb_type)))
               CASE ('sum')
                  tdata = tdata + this_coosig % a(ic) *                        &
     &               sdata(ind) % data(1)
                  tsigma = tsigma + this_coosig % a(ic) *                      &
     &               sdata(ind) % sigma(1)
               CASE ('max')
                  tdata = MAX(tdata,this_coosig % a(ic) *                      &
     &               sdata(ind) % data(1))
                  tsigma = MAX(tsigma,this_coosig % a(ic) *                    &
     &               sdata(ind) % sigma(1))
               CASE ('min')
                  tdata = MIN(tdata,this_coosig % a(ic) *                      &
     &               sdata(ind) % data(1))
                  tsigma = MIN(tsigma,this_coosig % a(ic) *                    &
     &               sdata(ind) % sigma(1))

               CASE DEFAULT
                  CALL err_fatal(sub_name //'unrecognized comb_type: ',        &
     &               char=this_coosig % comb_type)
            END SELECT ! Different coding depending on comb_type
         END DO ! of the loop over signal to combine
         
         sdata(isig) % data(1) = tdata
         sdata(isig) % sigma(1) = tsigma
      END DO ! of the second loop over the arrays

      RETURN
      END SUBROUTINE signal_mc_model_compute_a

!*******************************************************************************
! SECTION X.  COMMENTS FOR DIFFERENT REVISIONS
!*******************************************************************************
!  
!  JDH 11-28-2004
!     First version of module
!
!  JDH 12-14-2004
!     Moved from signal_eq to signal_model
!
!  JDH 2007-06-12
!     Changed from signal_mrf to mddc_mrf located in diagnostic_desc
!
!  JDH 2007-10-05
!     Slight change to error trapping in compute_mddc
!
!  JDH 2008-01-19 - Added => null() to pointer declarations
!
!  JDH 2009-01-30. Added subroutines signal_model_compute_edge_limit and
!     signal_model_aux_rz_one.
!
!  JDH 2009-03-03. Revised structure, removed stuff, changed name
!    from signal_model to signal_mc
!
!  JDH 2010-06-11. Added GJH coding for sxrch.
!
!  JDH 2011-07-07. Removed sdesc argument from signal_mc_model_compute. Now the
!    signal description comes from the pointer in the signal data.
!    Added coding for coosig (Combination Of Other SIGnals)
!    Added signal_mc_model_compute_a - coosig's computed in here.
!
!  JDH 2011-08-01
!    Refactor sxrc -> sxrch
!
!  JDH 2012-03-15
!    Added coding for ipch
!

      END MODULE signal_mc
