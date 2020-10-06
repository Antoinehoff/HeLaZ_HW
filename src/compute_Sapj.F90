SUBROUTINE compute_Sapj

  USE array, ONLY : dnjs, Sepj, Sipj
  USE basic
  USE fourier, ONLY : convolve_2D_F2F
  USE fields!, ONLY : phi, moments_e, moments_i
  USE grid
  USE model
  USE prec_const
  USE time_integration!, ONLY : updatetlevel
  IMPLICIT NONE

  COMPLEX(dp), DIMENSION(Nkr,Nkz) :: F_, G_, CONV
  REAL(dp):: kr, kz
  
  !!!!!!!!!!!!!!!!!!!! Zeta non linear term computation (Sepj)!!!!!!!!!!
  !! First convolution (term dr phi dz g)
  DO ikr = 1,Nkr ! Loop over kr
    DO ikz = 1,Nkz ! Loop over kz
      kr     = krarray(ikr)
      kz     = kzarray(ikz)
      ! First convolution term
      F_(ikr,ikz) = imagu*kr * phi(ikr,ikz)
      ! Second convolution term
      G_(ikr,ikz) = imagu*kz * moments_e(1,1,ikr,ikz,updatetlevel)
    ENDDO
  ENDDO

  CALL convolve_2D_F2F( F_, G_, CONV ) ! Convolve and go back to Fourier space
  Sepj(1,1,:,:) = CONV ! Add it to Sepj

  !! Second convolution (term -dz phi dr g)
  DO ikr = 1,Nkr ! Loop over kr
    DO ikz = 1,Nkz ! Loop over kz
      kr     = krarray(ikr)
      kz     = kzarray(ikz)
      ! First convolution term
      F_(ikr,ikz) = imagu*kz * phi(ikr,ikz)
      ! Second convolution term
      G_(ikr,ikz) = imagu*kr * moments_e(1,1,ikr,ikz,updatetlevel)
    ENDDO
  ENDDO

  CALL convolve_2D_F2F( F_, G_, CONV ) ! Convolve and go back to Fourier space
  Sepj(1,1,:,:) = (Sepj(1,1,:,:) - CONV)*Lr*Lz ! Add it to Sepj
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!! Density non linear term computation (Sipj)!!!!!!!!!!
  !! First convolution (term dr phi dz g)
  DO ikr = 1,Nkr ! Loop over kr
    DO ikz = 1,Nkz ! Loop over kz
      kr     = krarray(ikr)
      kz     = kzarray(ikz)
      ! First convolution term
      F_(ikr,ikz) = imagu*kr  * phi(ikr,ikz)
      ! Second convolution term
      G_(ikr,ikz) = imagu*kz  * moments_i(1,1,ikr,ikz,updatetlevel)
    ENDDO
  ENDDO

  CALL convolve_2D_F2F( F_, G_, CONV ) ! Convolve and go back to Fourier space
  Sipj(1,1,:,:) = CONV ! Add it to Sepj

  !! second convolution (term -dz phi dr g)
  DO ikr = 1,Nkr ! Loop over kr
    DO ikz = 1,Nkz ! Loop over kz
      kr     = krarray(ikr)
      kz     = kzarray(ikz)
      ! First convolution term
      F_(ikr,ikz) = imagu*kz * phi(ikr,ikz)
      ! Second convolution term
      G_(ikr,ikz) = imagu*kr * moments_i(1,1,ikr,ikz,updatetlevel)
    ENDDO
  ENDDO

  CALL convolve_2D_F2F( F_, G_, CONV ) ! Convolve and go back to Fourier space
  Sipj(1,1,:,:) = (Sipj(1,1,:,:) - CONV)*Lr*Lz ! Add it to Sepj
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE compute_Sapj
