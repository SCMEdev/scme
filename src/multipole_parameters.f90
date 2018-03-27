module multipole_parameters
  
  ! multipole moments
  
  ! unpolarized multipole moments, previously read from file multipoles
  ! public: d0, q0, o0, h0
  
  use data_types, only:dp, a0_A
  
  implicit none
  
  private
  
  ! au to Debye constant
  
  ! Ang to au constant
  !real(dp), parameter :: ang_to_au = 1.88972666351031921149_dp
  !real(dp), parameter :: au_to_debye = 2.5417709_dp
  !real(dp), parameter :: a0_len  = a0_A!1.0_dp/ang_to_au    !bohr to other length, now (angstrom, A)
  !real(dp), parameter :: ea0_dip = a0_A!au_to_debye
  
  
  ! components of the unpolarized dipole, quadrupole, ...
!  real(dp), parameter :: d0_1 = -0.72981_dp
!  real(dp), parameter :: q0_1 = 1.95532_dp, q0_2 = -1.85867_dp, q0_3 = -0.09665_dp
  real(dp), parameter :: d0_1 = -0.76219201_dp 
   real(dp), parameter :: q0_1 = 1.950317739166_dp, q0_2 = -1.840618450315_dp, q0_3 = -0.109699288851_dp 
   
  real(dp), parameter :: o0_1 = -3.27190_dp, o0_2 = 1.36606_dp, o0_3 = 1.90585_dp
  
  real(dp), parameter :: h0_1 = -0.94903_dp, h0_2 = -3.38490_dp, h0_3 = 4.33393_dp, &
                         h0_4 = 4.09835_dp,  h0_5 = -0.71345_dp, h0_6 = -3.62048_dp
  ! jo from gaussian:
  !jö traceless:
  !real(dp), parameter :: q0_1 = 1.941232948358537_dp ,&
  !                       q0_2 = -1.840976227504731_dp,& 
  !                       q0_3 = -0.100256704235994_dp 
  !tracefull
  !real(dp), parameter :: q0_1 = -4.964547074614011_dp/1.5_dp ,&
  !                       q0_2 = -8.746756250477278_dp/1.5_dp,& 
  !                       q0_3 = -7.006148149636545_dp/1.5_dp 
                         
  
  ! Compressed rigid tensors
  real(dp), public, parameter :: q1_00(3)  =  a0_A**1 *[0.0_dp, 0.0_dp, d0_1]
  real(dp), public, parameter :: q2_00(6)  =  a0_A**2 *[q0_1,0.0_dp, 0.0_dp, q0_2,0.0_dp,q0_3]
  real(dp), public, parameter :: q3_00(10) =  a0_A**3 *[0.0_dp,0.0_dp,o0_1,0.0_dp,0.0_dp,0.0_dp,0.0_dp,o0_2,0.0_dp,o0_3]
  real(dp), public, parameter :: q4_00(15) =  a0_A**4 *[h0_1,0.0_dp,0.0_dp,h0_2,0.0_dp,h0_3,0.0_dp,0.0_dp,0.0_dp,0.0_dp,h0_4,0.0_dp,h0_5,0.0_dp,h0_6]
                                                        !11   12      13      22   23     33
                                                        !111   112    113   122    123    133    222    223  233   333
                                                        !1111  1112   1113  1122  1123  1133  1222   1223  1233   1333   2222  2223  2233  2333  3333
  !real(dp), parameter :: qn_00(35) = [0.0_dp,q1_00,q2_00,q3_00,q4_00]
  
  
  !Full rigid tensors
  real(dp), public, parameter, dimension(3) :: d0 = a0_A**1 * [0.0_dp, 0.0_dp, d0_1]
  
  real(dp), public, parameter, dimension(3,3) :: q0 = a0_A**2 * &
       reshape([q0_1, 0.0_dp, 0.0_dp, &
       0.0_dp, q0_2, 0.0_dp, &
       0.0_dp, 0.0_dp, q0_3], shape(q0), order = [2,1])
  
  real(dp), public, parameter, dimension(3,3,3) :: o0 = a0_A**3 * &
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
  
  real(dp), public, parameter, dimension(3,3,3,3) :: h0 = a0_A**4 * &
       reshape([h0_1, 0.0_dp, 0.0_dp, 0.0_dp, h0_2, 0.0_dp, 0.0_dp, 0.0_dp, h0_3,  &  ! (1,1,1,1), ... (1,1,3,3)
       !        1111  1112    1113    1121    1122  1123    1131    1132    1133
       0.0_dp, h0_2, 0.0_dp, h0_2, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &  ! (1,2,1,1), ... (1,2,3,3)
       0.0_dp, 0.0_dp, h0_3, 0.0_dp, 0.0_dp, 0.0_dp, h0_3, 0.0_dp, 0.0_dp, &  ! (1,3,1,1), ... (1,3,3,3)
       0.0_dp, h0_2, 0.0_dp, h0_2, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &  ! (2,1,1,1), ... (2,1,3,3)
       h0_2, 0.0_dp, 0.0_dp, 0.0_dp, h0_4, 0.0_dp, 0.0_dp, 0.0_dp, h0_5,  &  ! (2,2,1,1), ... (2,2,3,3)
       !2211 2212    2213    2221    2222  2223    2231    2232    2233
       0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, h0_5, 0.0_dp, h0_5, 0.0_dp, &  ! (2,3,1,1), ... (2,3,3,3)
       0.0_dp, 0.0_dp, h0_3, 0.0_dp, 0.0_dp, 0.0_dp, h0_3, 0.0_dp, 0.0_dp, &  ! (3,1,1,1), ... (3,1,3,3)
       0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, h0_5, 0.0_dp, h0_5, 0.0_dp, &  ! (3,2,1,1), ... (3,2,3,3)
       h0_3, 0.0_dp, 0.0_dp, 0.0_dp, h0_5, 0.0_dp, 0.0_dp, 0.0_dp, h0_6],  &  ! (3,3,1,1), ... (3,3,3,3)
       shape(h0), order = [4,3,2,1])
  
  
end module multipole_parameters
