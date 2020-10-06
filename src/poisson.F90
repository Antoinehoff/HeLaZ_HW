SUBROUTINE poisson
  ! Solve poisson equation to get phi

  USE time_integration, ONLY: updatetlevel
  USE array
  USE fields
  USE grid

  USE prec_const
  IMPLICIT NONE

  REAL(dp)    :: kr, kz, kperp2

  DO ikr=ikrs,ikre
    DO ikz=ikzs,ikze
      kr     = krarray(ikr)
      kz     = kzarray(ikz)
      kperp2 = kr**2 + kz**2

      IF ( kperp2 .NE. 0) THEN
        phi(ikr,ikz) = -moments_e(1,1,ikr,ikz,updatetlevel)/kperp2
      ELSE
        phi(ikr,ikz) = 0._dp
      ENDIF

    END DO
  END DO

END SUBROUTINE poisson
