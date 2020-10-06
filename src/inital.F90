!******************************************************************************!
!!!!!! initialize the moments and load/build coeff tables
!******************************************************************************!
SUBROUTINE inital

  USE basic
  USE prec_const

  implicit none
  !!!!!! Set the moments arrays Nepj, Nipj !!!!!!
  CALL init_moments

  IF (RESTART) THEN
    CALL load_cp
  ENDIF

END SUBROUTINE inital
!******************************************************************************!

!******************************************************************************!
!!!!!!! Initialize the moments randomly
!******************************************************************************!
SUBROUTINE init_moments
  USE basic
  USE grid
  USE fields
  USE initial_par
  USE time_integration
  USE prec_const
  USE model, ONLY: NON_LIN
  IMPLICIT NONE

  REAL(dp) ::kr, kz, sigma, gain
  INTEGER, DIMENSION(12) :: iseedarr
  REAL(dp) :: noise

  ! Seed random number generator
  iseedarr(:)=iseed
  CALL RANDOM_SEED(PUT=iseedarr)

  CALL set_updatetlevel(1)

  sigma = 2._dp
  gain  = 0.5_dp
  !**** Gaussian initialization (Hakim 2017)
    DO ikr=1,nkr
      kr = krarray(ikr)
      DO ikz=1,nkz
        kz = kzarray(ikz)
        moments_i( 1,1, ikr,ikz, :) = gain*sigma/SQRT2 * exp(-(kr**2+kz**2)*sigma**2/4._dp)
        phi(ikr,ikz)                = moments_i( 1,1, ikr,ikz, 1)
        moments_e( 1,1, ikr,ikz, :) = -(kr**2+kz**2) * phi(ikr,ikz)
      END DO
    END DO

    CALL poisson

    IF (NON_LIN) THEN
      CALL compute_Sapj
    ENDIF
END SUBROUTINE init_moments
!______________________________________________________________________________!

!______________________________________________________________________________!
!!!!!!!! Load stored state to continue a previous simulation
!______________________________________________________________________________!
SUBROUTINE load_cp
  USE basic
  USE futils,          ONLY: openf, closef, getarr, getatt
  USE grid
  USE fields
  USE diagnostics_par
  USE time_integration
  IMPLICIT NONE

  WRITE(rstfile,'(a,a3)') TRIM(rstfile0),'.h5'

  WRITE(*,'(3x,a)') "Resume from previous run"

  CALL openf(rstfile, fidrst)
  CALL getatt(fidrst, '/Basic', 'cstep', cstep)
  CALL getatt(fidrst, '/Basic', 'time', time)
  CALL getatt(fidrst, '/Basic', 'jobnum', jobnum)
  jobnum = jobnum+1
  CALL getatt(fidrst, '/Basic', 'iframe2d',iframe2d)
  CALL getatt(fidrst, '/Basic', 'iframe5d',iframe5d)
  iframe2d = iframe2d-1; iframe5d = iframe5d-1

  ! Read state of system from restart file
  CALL getarr(fidrst, '/Basic/moments_e', moments_e(ips_e:ipe_e,ijs_e:ije_e,ikrs:ikre,ikzs:ikze,1),ionode=0)
  CALL getarr(fidrst, '/Basic/moments_i', moments_i(ips_i:ipe_i,ijs_i:ije_i,ikrs:ikre,ikzs:ikze,1),ionode=0)
  CALL closef(fidrst)

  WRITE(*,'(3x,a)') "Reading from restart file "//TRIM(rstfile)//" completed!"

  ! Update auxiliary arrays based on moments initialization
  CALL poisson
  CALL compute_Sapj

END SUBROUTINE load_cp
