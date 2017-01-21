module multipole_parameters
  
  ! multipole moments
  
  ! unpolarized multipole moments, previously read from file multipoles
  ! public: d0, q0, o0, h0
  
  use data_types
  
  implicit none
  
  private
  
  ! au to Debye constant
  
  ! Ang to au constant
  real(dp), parameter :: ang_to_au = 1.88972666351031921149_dp
  real(dp), parameter :: a0_len = 1.0_dp/ang_to_au    !bohr to other length, now (angstrom, A)
  real(dp), parameter :: au_to_debye = 2.5417709_dp
  real(dp), parameter :: ea0_dip = a0_len!au_to_debye
  
  
  ! components of the unpolarized dipole, quadrupole, ...
!  real(dp), parameter :: d0_1 = -0.72981_dp
!  real(dp), parameter :: q0_1 = 1.95532_dp, q0_2 = -1.85867_dp, q0_3 = -0.09665_dp
  real(dp), parameter :: d0_1 = -0.76219201_dp 
  real(dp), parameter :: q0_1 = 1.950317739166_dp, q0_2 = -1.840618450315_dp, q0_3 = -0.109699288851_dp 
  real(dp), parameter :: o0_1 = -3.27190_dp, o0_2 = 1.36606_dp, o0_3 = 1.90585_dp
  real(dp), parameter :: h0_1 = -0.94903_dp, h0_2 = -3.38490_dp, h0_3 = 4.33393_dp, &
                         h0_4 = 4.09835_dp,  h0_5 = -0.71345_dp, h0_6 = -3.62048_dp

  
  real(dp), public, parameter, dimension(3) :: d0 = ea0_dip * [0.0_dp, 0.0_dp, d0_1]
  
  real(dp), public, parameter, dimension(3,3) :: q0 = ea0_dip*a0_len * &
       reshape([q0_1, 0.0_dp, 0.0_dp, &
       0.0_dp, q0_2, 0.0_dp, &
       0.0_dp, 0.0_dp, q0_3], shape(q0), order = [2,1])
  
  real(dp), public, parameter, dimension(3,3,3) :: o0 = ea0_dip*a0_len**2 * &
       reshape([0.0_dp, 0.0_dp, o0_1,  &       ! (1,1,1), (1,1,2), (1,1,3)
       0.0_dp, 0.0_dp, 0.0_dp, &       ! (1,2,1), (1,2,2), (1,2,3)
       o0_1, 0.0_dp, 0.0_dp,  &       ! (1,3,1), (1,3,2), (1,3,3)
       0.0_dp, 0.0_dp, 0.0_dp, &       ! (2,1,1), (2,1,2), (2,1,3)
       0.0_dp, 0.0_dp, o0_2,  &       ! (2,2,1), (2,2,2), (2,2,3)
       0.0_dp, o0_2, 0.0_dp,  &       ! (2,3,1), (2,3,2), (2,3,3)
       o0_1, 0.0_dp, 0.0_dp,  &       ! (3,1,1), (3,1,2), (3,1,3)
       0.0_dp, o0_2, 0.0_dp,  &       ! (3,2,1), (3,2,2), (3,2,3)
       0.0_dp, 0.0_dp, o0_3], &       ! (3,3,1), (3,3,2), (3,3,3)
       shape(o0), order = [3,2,1])
  
  real(dp), public, parameter, dimension(3,3,3,3) :: h0 = ea0_dip*a0_len**3 * &
       reshape([h0_1, 0.0_dp, 0.0_dp, 0.0_dp, h0_2, 0.0_dp, 0.0_dp, 0.0_dp, h0_3,  &  ! (1,1,1,1), ... (1,1,3,3)
       0.0_dp, h0_2, 0.0_dp, h0_2, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &  ! (1,2,1,1), ... (1,2,3,3)
       0.0_dp, 0.0_dp, h0_3, 0.0_dp, 0.0_dp, 0.0_dp, h0_3, 0.0_dp, 0.0_dp, &  ! (1,3,1,1), ... (1,3,3,3)
       0.0_dp, h0_2, 0.0_dp, h0_2, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &  ! (2,1,1,1), ... (2,1,3,3)
       h0_2, 0.0_dp, 0.0_dp, 0.0_dp, h0_4, 0.0_dp, 0.0_dp, 0.0_dp, h0_5,  &  ! (2,2,1,1), ... (2,2,3,3)
       0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, h0_5, 0.0_dp, h0_5, 0.0_dp, &  ! (2,3,1,1), ... (2,3,3,3)
       0.0_dp, 0.0_dp, h0_3, 0.0_dp, 0.0_dp, 0.0_dp, h0_3, 0.0_dp, 0.0_dp, &  ! (3,1,1,1), ... (3,1,3,3)
       0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, h0_5, 0.0_dp, h0_5, 0.0_dp, &  ! (3,2,1,1), ... (3,2,3,3)
       h0_3, 0.0_dp, 0.0_dp, 0.0_dp, h0_5, 0.0_dp, 0.0_dp, 0.0_dp, h0_6],  &  ! (3,3,1,1), ... (3,3,3,3)
       shape(h0), order = [4,3,2,1])
  
  
end module multipole_parameters
