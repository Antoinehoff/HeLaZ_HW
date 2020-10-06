SUBROUTINE stepon
  !   Advance one time step, (num_step=4 for Runge Kutta 4 scheme)

  USE basic
  USE time_integration
  USE fields, ONLY: moments_e, moments_i, phi
  USE array , ONLY: moments_rhs_e, moments_rhs_i, Sepj, Sipj
  USE grid
  USE advance_field_routine, ONLY: advance_time_level, advance_field
  USE model
  USE utility, ONLY: checkfield
  use prec_const
  !USE convolution, ONLY: filter_Sapj, Orszag_filter
  IMPLICIT NONE

  INTEGER :: num_step, ip,ij


   DO num_step=1,ntimelevel ! eg RK4 compute successively k1, k2, k3, k4

      ! Compute right hand side of moments hierarchy equation
      CALL moments_eq_rhs
      ! Advance from updatetlevel to updatetlevel+1 (according to num. scheme)
      CALL advance_time_level
      ! Update the moments with the hierarchy RHS (step by step)
      CALL advance_field(moments_e(1,1,:,:,:),moments_rhs_e(1,1,:,:,:))
      CALL advance_field(moments_i(1,1,:,:,:),moments_rhs_i(1,1,:,:,:))
      ! Update electrostatic potential
      CALL poisson
      ! Update nonlinear term
      IF (NON_LIN) THEN
        CALL compute_Sapj(Sepj,Sipj)
      ENDIF

      CALL checkfield_all()
   END DO

   CONTAINS

      SUBROUTINE checkfield_all ! Check all the fields for inf or nan
        IF(.NOT.nlend) THEN
           nlend=nlend .or. checkfield(phi,' phi')
           DO ip=ips_e,ipe_e
             DO ij=ijs_e,ije_e
              nlend=nlend .or. checkfield(moments_e(ip,ij,:,:,updatetlevel),' moments_e')
             ENDDO
           ENDDO

           DO ip=ips_i,ipe_i
             DO ij=ijs_i,ije_i
              nlend=nlend .or. checkfield(moments_i(ip,ij,:,:,updatetlevel),' moments_i')
             ENDDO
           ENDDO

           DO ip=ips_i,ipe_i
            DO ij=ijs_i,ije_i
              nlend=nlend .or. checkfield(Sepj(ip,ij,:,:),' Sepj')
              nlend=nlend .or. checkfield(Sipj(ip,ij,:,:),' Sipj')
            ENDDO
          ENDDO
        ENDIF
      END SUBROUTINE checkfield_all

END SUBROUTINE stepon
