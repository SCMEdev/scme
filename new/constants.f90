module constants
use data_types
implicit none
  real(dp), parameter :: Eh_J        = 4.35974434e-18_dp        !! CODATA 2010
  real(dp), parameter :: Na          = 6.02214129e+23_dp        !! CODATA 2010
  real(dp), parameter :: kcal_J      = 4184.0_dp                
  real(dp), parameter :: Eh_kcalmol  = Eh_J*Na/kcal_J;        
  real(dp), parameter :: Bohr_A      = 0.52917721092_dp         !! CODATA 2010
  real(dp), parameter :: A_Bohr      = 1.0_dp
  real(dp), parameter :: c0          = 299792458.0_dp           !! m/s CODATA 2010
  real(dp), parameter :: ea0         = 8.47835326e-30_dp        !! C*m CODATA 2010
  real(dp), parameter :: D_au        = (1.0_dp/c0)*1.0e-21_dp/ea0  !! e * Bohr
  real(dp), parameter :: D           = D_au*Bohr_A;           !! e * A
  real(dp), parameter :: h_Js        = 6.62606957e-34_dp        !! J*s CODATA 2010
  real(dp), parameter :: hbar_Js     = 1.054571726e-34_dp       !! J*s CODATA 2010
  real(dp), parameter :: Eh_cm1      = 1.0e-2_dp*Eh_J/(c0*h_Js) !! cm-1
  real(dp), parameter :: cm1_kcalmol = Eh_kcalmol/Eh_cm1;
  real(dp), parameter :: cm1_eV      = 1.239842573148e-4_dp
  real(dp), parameter :: kB          = 1.3806488e-23_dp         !! JK-1 CODATA 2010
  real(dp), parameter :: e           =  1.602176565e-19_dp      !! C CODATA 2010
  real(dp), parameter :: E_cc        = 1.0e-7_dp*(c0*e*c0*e)/1.0e-10_dp !! in J !!interaction energy of 2 unit charges 1A apart
  real(dp), parameter :: CHARGECON   = dsqrt(E_cc*Na/kcal_J)  !! from NIST web site
  real(dp), parameter :: H_mass      = 1.00782503207_dp
  real(dp), parameter :: O_mass      = 15.99491461956_dp
end module

!! The _dp (or d0) is important
