SUBROUTINE moments_eq_rhs

  USE basic
  USE time_integration
  USE array
  USE fields
  USE grid
  USE model
  USE prec_const
  IMPLICIT NONE

  REAL(dp)    :: kr, kz, kperp2
!write(*,*) '----------------------------------------'

      ! Loop on kspace
  DO ikr = ikrs,ikre
    DO ikz = ikzs,ikze
      kr     = krarray(ikr)   ! Poloidal wavevector
      kz     = kzarray(ikz)   ! Toroidal wavevector
      kperp2 = kr**2 + kz**2  ! perpendicular wavevector

      !!!!!!!!! Zeta RHS !!!!!!!!!
      moments_rhs_e(1,1,ikr,ikz,updatetlevel) = &
      + lambdaD * (phi(ikr,ikz)-moments_i(1,1,ikr,ikz,updatetlevel))&
      - nu * kperp2**2 * moments_e(1,1,ikr,ikz,updatetlevel)

      !!!!!!!!! n RHS !!!!!!!!!
      moments_rhs_i(1,1,ikr,ikz,updatetlevel) = &
      + lambdaD * (phi(ikr,ikz)-moments_i(1,1,ikr,ikz,updatetlevel))&
      - nu * kperp2**2 * moments_i(1,1,ikr,ikz,updatetlevel)&
      - eta_n * imagu * kz * phi(ikr,ikz)

      !!!! Non linear term
      IF (NON_LIN) THEN
        ! Zeta
        moments_rhs_e(1,1,ikr,ikz,updatetlevel) = &
            moments_rhs_e(1,1,ikr,ikz,updatetlevel) - Sepj(1,1,ikr,ikz)
        ! n
        moments_rhs_i(1,1,ikr,ikz,updatetlevel) = &
            moments_rhs_i(1,1,ikr,ikz,updatetlevel) - Sipj(1,1,ikr,ikz)
      ENDIF

    END DO
  END DO

  ! Modified Hasegawa Wakatni (f(z=0)=0 to substract z-average of the fields )
  IF (MHW) THEN
    ! Zeta
    moments_rhs_e(1,1,:,ikz_0,updatetlevel) = &
        moments_rhs_e(1,1,:,ikz_0,updatetlevel) - lambdaD * &
        (phi(:,ikz_0)-moments_i(1,1,:,ikz_0,updatetlevel))
    ! n
    moments_rhs_i(1,1,:,ikz_0,updatetlevel) = &
        moments_rhs_i(1,1,:,ikz_0,updatetlevel) - lambdaD * &
        (phi(:,ikz_0)-moments_i(1,1,:,ikz_0,updatetlevel))
  ENDIF

END SUBROUTINE moments_eq_rhs
