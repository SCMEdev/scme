module polariz_parameters
  
  ! multipole moment polarizabilities, previously read from file multipoles
  ! public: dd0, dq0, hp0,qq0
  
  use data_types, only:dp
  
  implicit none
  
  private
  
  ! au to Debye constant
  !real(dp), parameter :: au_to_debye = 2.5417709_dp
  ! Ang to au constant
  real(dp), parameter :: ang_to_au = 1.88972666351031921149_dp
  
  ! components of the polarizabilities
!  real(dp), parameter :: x_x = 10.31146_dp
!  real(dp), parameter :: y_y = 9.54890_dp
!  real(dp), parameter :: z_z = 9.90656_dp
  real(dp), parameter :: x_x = 9.932447_dp
  real(dp), parameter :: y_y = 9.438219_dp
  real(dp), parameter :: z_z = 9.638610_dp
  
!  real(dp), parameter :: x_xz = -8.42037_dp
!  real(dp), parameter :: y_yz = -1.33400_dp
!  real(dp), parameter :: z_xx = -2.91254_dp
!  real(dp), parameter :: z_yy =  4.72407_dp
!  real(dp), parameter :: z_zz = -1.81153_dp
  real(dp), parameter :: x_xz = -6.602274_dp
  real(dp), parameter :: y_yz = -2.613041_dp
  real(dp), parameter :: z_xx = -1.016789_dp
  real(dp), parameter :: z_yy = 3.819314_dp
  real(dp), parameter :: z_zz = -2.802525_dp 
 
!  real(dp), parameter :: xx_xx = 12.11907_dp
!  real(dp), parameter :: xx_yy = -6.95326_dp
!  real(dp), parameter :: xx_zz = -5.16582_dp
!  real(dp), parameter :: xy_xy = 7.86225_dp
!  real(dp), parameter :: xz_xz = 11.98862_dp
!  real(dp), parameter :: yy_yy = 11.24741_dp
!  real(dp), parameter :: yy_zz = -4.29415_dp
!  real(dp), parameter :: yz_yz = 6.77226_dp
!  real(dp), parameter :: zz_zz = 9.45997_dp
  real(dp), parameter :: xx_xx = 42.633711_dp / 3.0_dp
  real(dp), parameter :: xx_yy = -23.219259_dp / 3.0_dp
  real(dp), parameter :: xx_zz = -19.414453_dp / 3.0_dp
  real(dp), parameter :: xy_xy = 33.508394_dp / 3.0_dp
  real(dp), parameter :: xz_xz = 38.044926_dp / 3.0_dp
  real(dp), parameter :: yy_yy = 45.676750_dp / 3.0_dp
  real(dp), parameter :: yy_zz = -22.457490_dp / 3.0_dp
  real(dp), parameter :: yz_yz = 33.104536_dp / 3.0_dp
  real(dp), parameter :: zz_zz = 41.871945_dp / 3.0_dp 
 
  
    !COMPRESSED matrices
    real(dp), public, parameter :: p11_data(3,3) = ang_to_au**(-3) * &
    reshape([&
    x_x,  0.0_dp, 0.0_dp, &
    0.0_dp, y_y,  0.0_dp, &
    0.0_dp, 0.0_dp, z_z &
    ],shape(p11_data), order = [2,1])
    
    real(dp), public, parameter :: p12_data(3,6) = ang_to_au**(-4) * &
    reshape([ &
    0.0_dp,0.0_dp,x_xz,0.0_dp,0.0_dp,0.0_dp, &
    ! 111   112    113   122    123    133
    0.0_dp,0.0_dp,0.0_dp,0.0_dp,y_yz,0.0_dp, &
    ! 211  212     213    222    223    233
    z_xx,0.0_dp,0.0_dp,z_yy,0.0_dp,z_zz &
    ! 311  312    313    322    323   333
    ],shape(p12_data), order = [2,1])
  

    ! 111 112 113 122 123 133
    ! 211 212 213 222 223 233
    ! 311 312 313 322 323 333

    !1: 11 12 13 22 23 33
    !2: 11 12 13 22 23 33
    !3: 11 12 13 22 23 33
    
    real(dp), public, parameter :: p21_data(6,3) = transpose(p12_data)
    
    real(dp), public, parameter :: p22_data(6,6) = ang_to_au**(-5) * &
    reshape([ &
    xx_xx, 0.0_dp, 0.0_dp, xx_yy, 0.0_dp, xx_zz, &
    ! 1111   1112   1113     1122   1123   1133 
    0.0_dp, xy_xy, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp,  &
    ! 1112   1212   1213     1222   1223   1233 
    0.0_dp, 0.0_dp, xz_xz, 0.0_dp, 0.0_dp, 0.0_dp,  &
    ! 1113   1213   1313     1322   1323   1333 
    xx_yy, 0.0_dp, 0.0_dp, yy_yy, 0.0_dp, yy_zz,  &
    ! 1122   1222   1322     2222   2223   2233 
    0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, yz_yz, 0.0_dp,  &
    ! 1123   1223   1323     2223   2323   2333 
    xx_zz, 0.0_dp, 0.0_dp, yy_zz, 0.0_dp, zz_zz   &
    ! 1133   1233   1333     2233   2333   3333 
    ],shape(p22_data), order = [2,1])
    

    


  ! quadrupole-quadrupole polarizability
  real(dp), public, parameter, dimension(3,3,3,3) :: qq0 = (ang_to_au)**(-5) * &
       reshape([&
       xx_xx, 0.0_dp, 0.0_dp, 0.0_dp, xx_yy, 0.0_dp, 0.0_dp, 0.0_dp, xx_zz,  &  ! (1,1,1,1), ... (1,1,3,3)
       0.0_dp, xy_xy, 0.0_dp, xy_xy, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &  ! (1,2,1,1), ... (1,2,3,3)
       0.0_dp, 0.0_dp, xz_xz, 0.0_dp, 0.0_dp, 0.0_dp, xz_xz, 0.0_dp, 0.0_dp, &  ! (1,3,1,1), ... (1,3,3,3)
       0.0_dp, xy_xy, 0.0_dp, xy_xy, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &  ! (2,1,1,1), ... (2,1,3,3)
       xx_yy, 0.0_dp, 0.0_dp, 0.0_dp, yy_yy, 0.0_dp, 0.0_dp, 0.0_dp, yy_zz,  &  ! (2,2,1,1), ... (2,2,3,3)
       0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, yz_yz, 0.0_dp, yz_yz, 0.0_dp, &  ! (2,3,1,1), ... (2,3,3,3)
       0.0_dp, 0.0_dp, xz_xz, 0.0_dp, 0.0_dp, 0.0_dp, xz_xz, 0.0_dp, 0.0_dp, &  ! (3,1,1,1), ... (3,1,3,3)
       0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, yz_yz, 0.0_dp, yz_yz, 0.0_dp, &  ! (3,2,1,1), ... (3,2,3,3)
       xx_zz, 0.0_dp, 0.0_dp, 0.0_dp, yy_zz, 0.0_dp, 0.0_dp, 0.0_dp, zz_zz], &  ! (3,3,1,1), ... (3,3,3,3)
       shape(qq0), order = [4,3,2,1])
  
    ! 1111 1112 1113 1122 1123 1133 
    ! 1112 1212 1213 1222 1223 1233 
    ! 1113 1213 1313 1322 1323 1333 
    ! 1122 1222 1322 2222 2223 2233 
    ! 1123 1223 1323 2223 2323 2333 
    ! 1133 1233 1333 2233 2333 3333 

    ! 11: 11 12 13 22 23 33 
    ! 12: 11 12 13 22 23 33 
    ! 13: 11 12 13 22 23 33 
    ! 22: 11 12 13 22 23 33 
    ! 23: 11 12 13 22 23 33 
    ! 33: 11 12 13 22 23 33 
  
    
  ! dipole-dipole polarizability
  real(dp), public, parameter, dimension(3,3) :: dd0 = (ang_to_au)**(-3) * &
       reshape([x_x, 0.0_dp, 0.0_dp, &
       0.0_dp, y_y, 0.0_dp, &
       0.0_dp, 0.0_dp, z_z ], shape(dd0), order = [2,1])
  
  ! dipole-quadrupole polarizability
  real(dp), public, parameter, dimension(3,3,3) :: dq0 = (ang_to_au)**(-4) * &
       reshape([0.0_dp, 0.0_dp, x_xz,  & ! (1,1,1), (1,1,2), (1,1,3)
       0.0_dp, 0.0_dp, 0.0_dp, &          ! (1,2,1), (1,2,2), (1,2,3)
       x_xz, 0.0_dp, 0.0_dp,  &          ! (1,3,1), (1,3,2), (1,3,3)
       0.0_dp, 0.0_dp, 0.0_dp, &          ! (2,1,1), (2,1,2), (2,1,3)
       0.0_dp, 0.0_dp, y_yz,  &          ! (2,2,1), (2,2,2), (2,2,3)
       0.0_dp, y_yz, 0.0_dp,  &          ! (2,3,1), (2,3,2), (2,3,3)
       z_xx, 0.0_dp, 0.0_dp,  &          ! (3,1,1), (3,1,2), (3,1,3)
       0.0_dp, z_yy, 0.0_dp,  &          ! (3,2,1), (3,2,2), (3,2,3)
       0.0_dp, 0.0_dp, z_zz], &          ! (3,3,1), (3,3,2), (3,3,3)
       shape(dq0), order = [3,2,1])
  
  
  
  
  
  
  ! hp0? set to zero in readPolariz
  real(dp), public, parameter, dimension(3,3,3) :: hp0 =  &
       reshape([0.0_dp, 0.0_dp, 0.0_dp, &
       0.0_dp, 0.0_dp, 0.0_dp, &
       0.0_dp, 0.0_dp, 0.0_dp, &
       0.0_dp, 0.0_dp, 0.0_dp, &
       0.0_dp, 0.0_dp, 0.0_dp, &
       0.0_dp, 0.0_dp, 0.0_dp, &
       0.0_dp, 0.0_dp, 0.0_dp, &
       0.0_dp, 0.0_dp, 0.0_dp, &
       0.0_dp, 0.0_dp, 0.0_dp],&
       shape(hp0), order = [3,2,1])
  
end module polariz_parameters
