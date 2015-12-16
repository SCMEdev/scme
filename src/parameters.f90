module parameters

  ! parameters in param.pot, multipoles and polar
  
  use data_types
  
  implicit none
  
  ! from param.pot
  
  real(dp), parameter :: tt_b = 4.40_dp              ! b parameter for Tang_Toennies
  real(dp), parameter :: coreInt_c1 = -1.44350_dp, coreInt_c2 = 20169.17400_dp, coreInt_c3 = -3.46581_dp, &
       coreInt_c4 = 1.0_dp, coreInt_c5_r = -1.0_dp
  real(dp), parameter :: rho_a1 = 1.0250853530d-01, rho_a2 = -1.7246118649d-04, rho_a3 = 1.0219555575d-07,&
       rho_a4 = -2.6087710711d-11, rho_a5 = 3.0605430553d-15, rho_a6 = -1.3290133918d-19
  integer, parameter :: rho_n = 5
  
  real(dp), dimension(3), parameter :: rcut = (/ 3.0, 3.0, 3.0/), rskin = (/ 3.0_dp, 3.0_dp, 3.0_dp/), skndpt = rskin-rcut
  real(dp), parameter :: Ro_0 = 0.9572_dp, theta_0 = 104.52_dp
  
  integer, parameter :: num_cells = 0 ! NC
  
end module parameters
