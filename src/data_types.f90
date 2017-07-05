module data_types
  
  implicit none
  
  private
  public dp, tt_b, num_cells,A_a0, a0_A, cm1_eV, coulomb_k !, ea0_Deb, eA_Deb, xyz, hho, kk1,kk2,
  
  integer, parameter :: dp = selected_real_kind(15, 307) ! double precision ~ 8 byte float

!scme.f90 >>
    real(dp), parameter :: pi = 3.14159265358979324_dp
    real(dp), parameter :: kk1       =            1.0_dp! *2.5417709_dp
    real(dp), parameter :: kk2       =            1.0_dp! *1.88972666351031921149_dp
    real(dp), parameter :: coulomb_k = 14.39975841_dp   ! /4.803206799_dp**2 ! Eh_eV*a0_A/ea0_deb**2 ish
!scme.f90 <<

!constants.f90  (or constants.h) >>
  real(dp), parameter :: cm1_eV      = 1.2398419739e-4_dp !Codata 2014    1.239842573148e-4_dp 
  real(dp), parameter :: a0_A        = 0.52917721067_dp   !Codata 2014      0.52917721092_dp   !! CODATA 2010 0.529 177 210 67(12) later
  real(dp), parameter :: A_a0        = 1.0_dp/a0_A
!constants.f90  (or constants.h) <<

!parameters.f90 >>
  real(dp), parameter :: tt_b         =  2.8811729_dp
  integer, parameter :: num_cells = 0 ! NC
!parameters.f90 <<
  

!position type
!  type h2o
!   real(dp), dimension(3) :: h1, h2, o
!  end type
!
!type xyz 
!  real(dp) :: x(3)
!end type
!
!type hho
!  type(xyz) :: a(3)
!end type


!force type  
!  type h2o_forces
!    real(dp), dimension(3) :: h1, h2, o
!  endtype
  
end module data_types

!  integer, parameter :: sp = selected_real_kind(6,37) ! single precision, ~ 4 byte float
!  integer, parameter :: medint = selected_int_kind(9) ! double precision ~ 4 byte float
!  integer, parameter :: longint = selected_int_kind(15) ! 15 digits, ~ 8 byte integer



!!constants.f90  (or constants.h) >>
!  real(dp), parameter :: cm1_eV      = 1.239842573148e-4_dp
!  real(dp), parameter :: a0_A        = 0.52917721092_dp         !! CODATA 2010
!  real(dp), parameter :: A_a0        = 1.0_dp/a0_A
!
!
!  real(dp), parameter :: Eh_J        = 4.35974434e-18_dp        !! CODATA 2010
!  real(dp), parameter :: Na          = 6.02214129e+23_dp        !! CODATA 2010
!  real(dp), parameter :: kcal_J      = 4184.0_dp                
!  real(dp), parameter :: Eh_kcalmol  = Eh_J*Na/kcal_J;        
!  real(dp), parameter :: c0          = 299792458.0_dp           !! m/s CODATA 2010
!  real(dp), parameter :: ea0         = 8.47835326e-30_dp        !! C*m CODATA 2010
!  real(dp), parameter :: Deb_ea0     = (1.0_dp/c0)*1.0e-21_dp/ea0  !! Debye to e * Bohr
!  real(dp), parameter :: ea0_Deb     = 1.0_dp/Deb_ea0               !kk1 
!  real(dp), parameter :: Deb_eA      = Deb_ea0*a0_A;           !! Debye to e * A
!  real(dp), parameter :: eA_Deb      = ea0_Deb*A_a0           !! Debye to e * A
!  real(dp), parameter :: h_Js        = 6.62606957e-34_dp        !! J*s CODATA 2010
!  real(dp), parameter :: hbar_Js     = 1.054571726e-34_dp       !! J*s CODATA 2010
!  real(dp), parameter :: Eh_cm1      = 1.0e-2_dp*Eh_J/(c0*h_Js) !! cm-1
!  real(dp), parameter :: cm1_kcalmol = Eh_kcalmol/Eh_cm1;
!  real(dp), parameter :: kB          = 1.3806488e-23_dp         !! JK-1 CODATA 2010
!  real(dp), parameter :: e           =  1.602176565e-19_dp      !! C CODATA 2010
!  real(dp), parameter :: E_cc        = 1.0e3_dp*(c0*e*c0*e)!1.0e-7_dp*(c0*e*c0*e)/1.0e-10_dp !! in J !!interaction energy of 2 unit charges 1A apart
!  real(dp), parameter :: CHARGECON   = dsqrt(E_cc*Na/kcal_J)  !! from NIST web site
!  real(dp), parameter :: H_mass      = 1.00782503207_dp
!  real(dp), parameter :: O_mass      = 15.99491461956_dp
!!constants.f90  (or constants.h) <<



!!parameters.f90 >>
!  real(dp), parameter :: coreInt_c1   = -1.44350_dp
!  real(dp), parameter :: coreInt_c2   = 20169.17400_dp
!  real(dp), parameter :: coreInt_c3   = -3.46581_dp
!  real(dp), parameter :: coreInt_c4   =  1.0_dp
!  real(dp), parameter :: coreInt_c5_r = -1.0_dp
!  real(dp), parameter :: rho_a1 =  1.0250853530e-01_dp !JÖ changed exponential from d-02 to e-02_dp 
!  real(dp), parameter :: rho_a2 = -1.7246118649e-04_dp !JÖ
!  real(dp), parameter :: rho_a3 =  1.0219555575e-07_dp !JÖ
!  real(dp), parameter :: rho_a4 = -2.6087710711e-11_dp !JÖ
!  real(dp), parameter :: rho_a5 =  3.0605430553e-15_dp !JÖ
!  real(dp), parameter :: rho_a6 = -1.3290133918e-19_dp !JÖ
!  integer, parameter :: rho_n = 5
!  
!  real(dp), dimension(3), parameter :: rcut = (/ 3.0_dp, 3.0_dp, 3.0_dp/) !JÖ added _dp to 3.0
!  real(dp), dimension(3), parameter :: rskin = (/ 3.0_dp, 3.0_dp, 3.0_dp/)
!  real(dp), dimension(3), parameter :: skndpt = rskin-rcut
!  real(dp), parameter :: Ro_0 = 0.9572_dp, theta_0 = 104.52_dp
!!parameters.f90 <<
  
