SUBROUTINE memory
  ! Allocate arrays (done dynamically otherwise size is unknown)

  USE array
  USE basic
  USE fields
  USE grid
  USE time_integration
  USE model, ONLY: CO, NON_LIN

  USE prec_const
  IMPLICIT NONE

  ! Moments and moments rhs
  CALL allocate_array(     moments_e, ips_e,ipe_e, ijs_e,ije_e, ikrs,ikre, ikzs,ikze, 1,ntimelevel )
  CALL allocate_array(     moments_i, ips_i,ipe_i, ijs_i,ije_i, ikrs,ikre, ikzs,ikze, 1,ntimelevel )
  CALL allocate_array( moments_rhs_e, ips_e,ipe_e, ijs_e,ije_e, ikrs,ikre, ikzs,ikze, 1,ntimelevel )
  CALL allocate_array( moments_rhs_i, ips_i,ipe_i, ijs_i,ije_i, ikrs,ikre, ikzs,ikze, 1,ntimelevel )

  ! Electrostatic potential
  CALL allocate_array(phi, ikrs,ikre, ikzs,ikze)

  ! Non linear terms and dnjs table
  CALL allocate_array( Sepj, ips_e,ipe_e, ijs_e,ije_e, ikrs,ikre, ikzs,ikze )
  CALL allocate_array( Sipj, ips_i,ipe_i, ijs_i,ije_i, ikrs,ikre, ikzs,ikze )
  CALL allocate_array( dnjs, 1,maxj+1, 1,maxj+1, 1,maxj+1)

END SUBROUTINE memory
