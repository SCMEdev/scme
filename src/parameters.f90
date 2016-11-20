module parameters

  ! parameters in param.pot, multipoles and polar
  
  use data_types
  
  implicit none
  
  ! from param.pot
  
!  real(dp), parameter :: tt_b = 4.40_dp              ! b parameter for Tang_Toennies

!JÖ changed the form to more readable:
  real(dp), parameter :: tt_b         =  2.8811729_dp
  real(dp), parameter :: coreInt_c1   = -1.44350_dp
  real(dp), parameter :: coreInt_c2   = 20169.17400_dp
  real(dp), parameter :: coreInt_c3   = -3.46581_dp
  real(dp), parameter :: coreInt_c4   =  1.0_dp
  real(dp), parameter :: coreInt_c5_r = -1.0_dp
  real(dp), parameter :: rho_a1 =  1.0250853530e-01_dp !JÖ changed exponential from d-02 to e-02_dp 
  real(dp), parameter :: rho_a2 = -1.7246118649e-04_dp !JÖ
  real(dp), parameter :: rho_a3 =  1.0219555575e-07_dp !JÖ
  real(dp), parameter :: rho_a4 = -2.6087710711e-11_dp !JÖ
  real(dp), parameter :: rho_a5 =  3.0605430553e-15_dp !JÖ
  real(dp), parameter :: rho_a6 = -1.3290133918e-19_dp !JÖ
  integer, parameter :: rho_n = 5
  
  real(dp), dimension(3), parameter :: rcut = (/ 3.0_dp, 3.0_dp, 3.0_dp/) !JÖ added _dp to 3.0
  real(dp), dimension(3), parameter :: rskin = (/ 3.0_dp, 3.0_dp, 3.0_dp/)
  real(dp), dimension(3), parameter :: skndpt = rskin-rcut
  real(dp), parameter :: Ro_0 = 0.9572_dp, theta_0 = 104.52_dp
  
  integer, parameter :: num_cells = 0 ! NC
  
end module parameters
