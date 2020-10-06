MODULE fourier
  USE prec_const
  USE MKL_DFTI
  USE grid
  implicit none

  PRIVATE
    ! MKL DFTI descripto for FFTs
  type(DFTI_DESCRIPTOR), POINTER :: Desc_fft_r
  type(DFTI_DESCRIPTOR), POINTER :: Desc_ifft_r
  type(DFTI_DESCRIPTOR), POINTER :: Desc_fft2_r,  Desc_fft2_z
  type(DFTI_DESCRIPTOR), POINTER :: Desc_ifft2_r, Desc_ifft2_z

  COMPLEX(dp), DIMENSION(:,:), ALLOCATABLE :: tmp_2D    ! Auxiliary fields
  COMPLEX(dp), DIMENSION(:)  , ALLOCATABLE :: tmp_kr_cpx
  COMPLEX(dp), DIMENSION(:)  , ALLOCATABLE :: tmp_kz_cpx
  COMPLEX(dp), DIMENSION(:)  , ALLOCATABLE :: tmp_r_cpx
  COMPLEX(dp), DIMENSION(:)  , ALLOCATABLE :: tmp_z_cpx
  REAL(dp),    DIMENSION(:)  , ALLOCATABLE :: tmp_r_real
  REAL(dp),    DIMENSION(:)  , ALLOCATABLE :: tmp_z_real

  PUBLIC :: fft_r2cc
  PUBLIC :: ifft_cc2r
  PUBLIC :: fft2_r2cc
  PUBLIC :: ifft2_cc2r
  PUBLIC :: convolve_2D_F2F
  PUBLIC :: set_descriptors, free_descriptors

  CONTAINS

  SUBROUTINE fft_r2cc(fx_in, Fk_out)
    IMPLICIT NONE

    REAL(dp),    DIMENSION(Nr),  INTENT(IN) :: fx_in
    COMPLEX(dp), DIMENSION(Nkr), INTENT(OUT):: Fk_out

    INTEGER :: Status ! Mkl status variable

    Status = DftiComputeForward  ( Desc_fft_r, fx_in, Fk_out)

  END SUBROUTINE fft_r2cc

  SUBROUTINE ifft_cc2r(Fk_in, fx_out)
    IMPLICIT NONE

    COMPLEX(dp), DIMENSION(Nkr),  INTENT(IN) :: Fk_in
    REAL(dp),    DIMENSION(Nr),   INTENT(OUT):: fx_out

    INTEGER :: Status ! Mkl status variable

    Status = DftiComputeBackward ( Desc_ifft_r, Fk_in, fx_out)

  END SUBROUTINE ifft_cc2r

  SUBROUTINE fft2_r2cc( ffx_, FFk_ )
      IMPLICIT NONE
      REAL(dp),    DIMENSION(Nr,Nz), INTENT(IN)   :: ffx_
      COMPLEX(dp), DIMENSION(Nkr,Nkz), INTENT(OUT):: FFk_
      INTEGER :: Status ! Mkl status variable

      !!! 2D Forward FFT ________________________!
      DO ir = 1,Nr
        tmp_z_real = ffx_(ir,:)
        Status = DftiComputeForward( Desc_fft2_z, tmp_z_real, tmp_kz_cpx)
        tmp_2D(ir,:) = tmp_kz_cpx
      ENDDO
      DO iz = 1,Nkz
        tmp_kr_cpx = tmp_2D(:,iz)
        Status = DftiComputeForward( Desc_fft2_r, tmp_kr_cpx)
        FFk_(:,iz) = tmp_r_cpx
      ENDDO

  END SUBROUTINE fft2_r2cc

  SUBROUTINE ifft2_cc2r( FFk_, ffx_ )
      IMPLICIT NONE
      COMPLEX(dp), DIMENSION(Nkr,Nkz), INTENT(IN)  :: FFk_
      REAL(dp),    DIMENSION(Nr,Nz), INTENT(OUT) :: ffx_
      INTEGER :: Status ! Mkl status variable

      !!! 2D Backward FFT ________________________!
      DO ikz = 1,Nkz
        tmp_r_cpx = FFk_(:,ikz)
        Status = DftiComputeBackward( Desc_ifft2_r, tmp_r_cpx)
        tmp_2D(:,ikz) = tmp_r_cpx
      ENDDO
      DO ir = 1,Nr
        tmp_kz_cpx = tmp_2D(ir,:)
        Status = DftiComputeBackward( Desc_ifft2_z, tmp_kz_cpx, tmp_z_real)
        ffx_(ir,:) = tmp_z_real
      ENDDO
      IF (Status .NE. 0) THEN
        WRITE(*,*) 'Warning : an error occured in iFFT2 computations...'
      ENDIF

  END SUBROUTINE ifft2_cc2r

  !!! Convolution 2D Fourier to Fourier
  !   - Compute the convolution using the convolution theorem and MKL
  SUBROUTINE convolve_2D_F2F( F_2D, G_2D, C_2D )
    IMPLICIT NONE
    COMPLEX(dp), DIMENSION(Nkr,Nkz), INTENT(IN)  :: F_2D, G_2D  ! input fields
    COMPLEX(dp), DIMENSION(Nkr,Nkz), INTENT(OUT) :: C_2D  ! output convolutioned field
    REAL(dp), DIMENSION(Nr,Nz) :: ff, gg ! iFFT of input fields
    REAL(dp), DIMENSION(Nr,Nz) :: ffgg  ! will receive the product of f*g in Real

    !!! Backward FFT of F -> f _______________________________________________________!
    CALL ifft2_cc2r(F_2D,ff)
    !!! Backward FFT of G -> g _______________________________________________________!
    CALL ifft2_cc2r(G_2D,gg)
    !!! Product c = f*g ______________________________________________________________!
    ffgg = ff*gg
    !!! Forward FFT of fg -> Conv.  Real -> Complex conjugate ________________________!
    CALL fft2_r2cc(ffgg,C_2D)

    !!! Anti-Aliasing 2/3 rule____________________!
    DO ikz = 1,Nkz
      DO ikr = 1,Nkr
        C_2D(ikr,ikz) = C_2D(ikr,ikz) * AA_r(ikr)*AA_z(ikz);
      ENDDO
    ENDDO

END SUBROUTINE convolve_2D_F2F

SUBROUTINE check_status(Status)
  IMPLICIT NONE
  INTEGER :: Status ! Mkl status variable
  IF (Status .NE. 0) THEN
    WRITE(*,*) 'Warning : an error occured in FFT2 computations...'
  ENDIF
END SUBROUTINE check_status

SUBROUTINE set_descriptors
  IMPLICIT NONE
  INTEGER :: Status
  WRITE(*,*) 'Set descriptors'

  Status = DftiCreateDescriptor( Desc_fft_r, DFTI_DOUBLE, DFTI_REAL, 1, Nkr)
  Status = DftiSetValue        ( Desc_fft_r, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
  Status = DftiCommitDescriptor( Desc_fft_r)

  Status = DftiCreateDescriptor( Desc_ifft_r, DFTI_DOUBLE, DFTI_REAL, 1, Nr)
  Status = DftiSetValue        ( Desc_ifft_r, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
  Status = DftiSetValue        ( Desc_ifft_r, DFTI_BACKWARD_SCALE, 1._dp/REAL(Nr,dp))
  Status = DftiCommitDescriptor( Desc_ifft_r)

  Status = DftiCreateDescriptor( Desc_fft2_z, DFTI_DOUBLE, DFTI_REAL, 1, Nz)
  Status = DftiSetValue        ( Desc_fft2_z, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
  Status = DftiSetValue        ( Desc_fft2_z, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
  Status = DftiCommitDescriptor( Desc_fft2_z)

  Status = DftiCreateDescriptor( Desc_fft2_r, DFTI_DOUBLE, DFTI_COMPLEX, 1, Nr)
  Status = DftiSetValue        ( Desc_fft2_r, DFTI_PLACEMENT, DFTI_INPLACE)
  Status = DftiSetValue        ( Desc_fft2_r, DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
  Status = DftiCommitDescriptor( Desc_fft2_r)

  Status = DftiCreateDescriptor( Desc_ifft2_r, DFTI_DOUBLE, DFTI_COMPLEX, 1, Nr)
  Status = DftiSetValue        ( Desc_ifft2_r, DFTI_PLACEMENT, DFTI_INPLACE)
  Status = DftiSetValue        ( Desc_ifft2_r, DFTI_BACKWARD_SCALE, 1._dp/REAL(Nr,dp))
  Status = DftiCommitDescriptor( Desc_ifft2_r)

  Status = DftiCreateDescriptor( Desc_ifft2_z, DFTI_DOUBLE, DFTI_REAL, 1, Nz)
  Status = DftiSetValue        ( Desc_ifft2_z, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
  Status = DftiSetValue        ( Desc_ifft2_z, DFTI_BACKWARD_SCALE, 1._dp/REAL(Nz,dp))
  Status = DftiCommitDescriptor( Desc_ifft2_z)
  CALL check_status(Status)

  ALLOCATE(tmp_2D(1:Nr,1:Nz));
  ALLOCATE(tmp_kr_cpx(1:Nkr));  ALLOCATE(tmp_kz_cpx(1:Nkz))
  ALLOCATE(tmp_r_cpx(1:Nr));    ALLOCATE(tmp_z_cpx(1:Nz))
  ALLOCATE(tmp_r_real(1:Nr));   ALLOCATE(tmp_z_real(1:Nkz))
END SUBROUTINE set_descriptors

SUBROUTINE free_descriptors
  IMPLICIT NONE
  INTEGER :: Status
  WRITE(*,*) 'Free descriptors'
  Status = DftiFreeDescriptor  ( Desc_fft_r)
  Status = DftiFreeDescriptor  ( Desc_ifft_r)
  Status = DftiFreeDescriptor  ( Desc_fft2_r)
  Status = DftiFreeDescriptor  ( Desc_fft2_z)
  Status = DftiFreeDescriptor  ( Desc_ifft2_r)
  Status = DftiFreeDescriptor  ( Desc_ifft2_z)
  CALL check_status(Status)

END SUBROUTINE free_descriptors

END MODULE fourier
