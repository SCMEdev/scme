! Copyright (C)  2015-2016  SCMEdev
! Licenced under LGPLv3. See LICENCE for details.


#include "mifu.h"

module test_scme
  use mifu_asserts
  use printer_mod, only: printer
  use data_types, only: au_to_debye, coulomb_k, a0_A, A_a0

  ! Include the module we test.
  use scme

  contains

  subroutine test_scme_monomer1()
    implicit none
    ! ----------------------------------------
    integer, parameter :: n_atoms = 3
    integer :: i
    real*8  :: u_tot = 0.0d0
    real*8  :: u_tot_ref = 0.0d0
    real*8,dimension(45000)     :: coordinates
    real*8,dimension(3)         :: trans
    real*8,dimension(3)         :: cell
    real*8,dimension(45000)     :: forces
    real*8,dimension(n_atoms*3) :: ref_forces
    ! ----------------------------------------

    ! ---------------------------------------------------------
    ! Test 1) Test the forces and energy against hardcoded
    !         reference values.
    ! ---------------------------------------------------------
    coordinates(1) =  0.018       ! x H1 on W1
    coordinates(2) = -0.739       ! y H1 on W1
    coordinates(3) =  0.521       ! z H1 on W1
    coordinates(4) = -0.815       ! x H2 on W1
    coordinates(5) = -0.673
    coordinates(6) = -0.592
    coordinates(7) =  0.000       ! x O on W1
    coordinates(8) =  0.000
    coordinates(9) =  0.000

    cell(1) = 35.000
    cell(2) = 35.000
    cell(3) = 35.000

    forces(:) = 0.0

    ! Call the scme function.
    call scme_calculate(n_atoms, coordinates, cell, forces, u_tot)

    ! Check the total energy against reference value.
    u_tot_ref = 1.3254072629995397d0
    MIFU_ASSERT_REAL_EQUAL(u_tot, u_tot_ref, 1.0e-10)

    ! Check the forces.
    ref_forces(1) =  1.2532152102008445d0
    ref_forces(2) = -2.6108158941929331d0
    ref_forces(3) =  3.3666330806012423d0
    ref_forces(4) =  3.2864124801108985d0
    ref_forces(5) =  3.7924992888626039d0
    ref_forces(6) =  1.6604054612584063d0
    ref_forces(7) = -4.5396276903117423d0
    ref_forces(8) = -1.1816833946696708d0
    ref_forces(9) = -5.0270385418596488d0

    ! Check.
    do i=1,9
       MIFU_ASSERT_REAL_EQUAL(ref_forces(i), forces(i), 1.0e-10)
    end do

    ! ---------------------------------------------------------
    ! Test 2) Test that at translation of the coordinates does
    !         not alter the forces and energies.
    ! ---------------------------------------------------------

    ! Translate the molecule and veify that the energy and forces do not change.
    trans(1) =  2.d0
    trans(2) = -1.d0
    trans(3) =  3.d0

    do i=0,2
       coordinates(i*3+1:i*3+3) = coordinates(i*3+1:i*3+3) + trans(1:3)
    end do

    ! Call the scme function to get the new forces and energy.
    call scme_calculate(n_atoms, coordinates, cell, forces, u_tot)

    ! Check that the energy and forces has not changed.
    MIFU_ASSERT_REAL_EQUAL(u_tot, u_tot_ref, 1.0e-10)
    do i=1,9
       MIFU_ASSERT_REAL_EQUAL(ref_forces(i), forces(i), 1.0e-10)
    end do

    ! ---------------------------------------------------------
    ! Test 3) Test that a long translation (way out of the box)
    !         does not alter the forces and energies.
    ! ---------------------------------------------------------

    ! Translate the molecule and veify that the energy and forces do not change.
    trans(1) =  1202.d0
    trans(2) =  1201.d0
    trans(3) = -2203.d0

    do i=0,2
       coordinates(i*3+1:i*3+3) = coordinates(i*3+1:i*3+3) + trans(1:3)
    end do

    ! Call the scme function to get the new forces and energy.
    call scme_calculate(n_atoms, coordinates, cell, forces, u_tot)

    ! Check that the energy and forces has not changed.
    MIFU_ASSERT_REAL_EQUAL(u_tot, u_tot_ref, 1.0e-10)
    do i=1,9
       MIFU_ASSERT_REAL_EQUAL(ref_forces(i), forces(i), 1.0e-10)
    end do

  end subroutine test_scme_monomer1


  subroutine test_scme_dimer1()
    implicit none
    ! ----------------------------------------
    integer, parameter :: n_atoms = 6
    integer :: i
    real*8  :: u_tot = 0.0d0
    real*8  :: u_tot_ref = 0.0d0
    real*8,dimension(45000)     :: coordinates
    real*8,dimension(3)         :: cell
    real*8,dimension(45000)     :: forces
    real*8,dimension(n_atoms*3) :: ref_forces
    ! ----------------------------------------

    ! ---------------------------------------------------------
    ! Test 1) Test the forces and energy against hardcoded
    !         reference values.
    ! ---------------------------------------------------------

    coordinates(1)  =  0.018       ! x H1 on W1
    coordinates(2)  = -0.739       ! y H1 on W1
    coordinates(3)  =  0.521       ! z H1 on W1
    coordinates(4)  = -0.815       ! x H2 on W1
    coordinates(5)  = -0.673
    coordinates(6)  = -0.592
    coordinates(7)  =  1.942       ! x H1 on W2
    coordinates(8)  = -1.964       ! y H1 on W2
    coordinates(9)  =  1.223       ! z H1 on W2
    coordinates(10) =  2.319       ! x H2 on W2
    coordinates(11) = -3.450
    coordinates(12) =  1.136
    coordinates(13) =  0.000       ! x O on W1
    coordinates(14) =  0.000
    coordinates(15) =  0.000
    coordinates(16) =  1.675       ! x O on W2
    coordinates(17) = -2.803
    coordinates(18) =  0.911

    cell(1) = 35.000
    cell(2) = 35.000
    cell(3) = 35.000

    forces(:) = 0.0

    ! Call the scme function.
    call scme_calculate(n_atoms, coordinates, cell, forces, u_tot)

    ! Check the total energy against reference value.
    u_tot_ref = 1.2904844927113648d0

    MIFU_ASSERT_REAL_EQUAL(u_tot_ref, u_tot, 1.0e-10)

    ! Check the forces against reference values.
    ref_forces(1)  =  1.3580482852700133d0
    ref_forces(2)  = -2.6108158941929331d0
    ref_forces(3)  =  3.5061422780398157d0
    ref_forces(4)  =  3.2864124801108985d0
    ref_forces(5)  =  3.9382894622077109d0
    ref_forces(6)  =  1.5157008465161304d0
    ref_forces(7)  =  0.80474379761950909d0
    ref_forces(8)  =  0.92790818581970669d0
    ref_forces(9)  =  0.52332582709118647d0
    ref_forces(10) =  0.76242032020453132d0
    ref_forces(11) = -0.30407784368126489d0
    ref_forces(12) =  0.37221067319580348d0
    ref_forces(13) = -4.6406547014988444d0
    ref_forces(14) = -1.4408190402173098d0
    ref_forces(15) = -5.0071366501451626d0
    ref_forces(16) = -1.5747776828267221d0
    ref_forces(17) = -0.51132962120957237d0
    ref_forces(18) = -0.91313228694841819d0

    ! Check.
    do i=1,18
       MIFU_ASSERT_REAL_EQUAL(ref_forces(i), forces(i), 1.0e-10)
    end do

  end subroutine test_scme_dimer1


  subroutine test_scme_perf()
    implicit none
    ! ----------------------------------------
    integer, parameter :: n_atoms = 6
    integer :: i
    real*8  :: u_tot = 0.0d0
    real*8  :: u_tot_ref = 0.0d0
    real*8,dimension(18)     :: coordinates
    real*8,dimension(3)         :: cell
    real*8,dimension(54)     :: forces
    real*8,dimension(n_atoms*3) :: ref_forces
    ! ----------------------------------------

    ! ---------------------------------------------------------
    ! Test 1) Looping calculation of forces and energy to
    !         run timings.
    ! ---------------------------------------------------------

    coordinates(1)  =  0.018       ! x H1 on W1
    coordinates(2)  = -0.739       ! y H1 on W1
    coordinates(3)  =  0.521       ! z H1 on W1
    coordinates(4)  = -0.815       ! x H2 on W1
    coordinates(5)  = -0.673
    coordinates(6)  = -0.592
    coordinates(7)  =  1.942       ! x H1 on W2
    coordinates(8)  = -1.964       ! y H1 on W2
    coordinates(9)  =  1.223       ! z H1 on W2
    coordinates(10) =  2.319       ! x H2 on W2
    coordinates(11) = -3.450
    coordinates(12) =  1.136
    coordinates(13) =  0.000       ! x O on W1
    coordinates(14) =  0.000
    coordinates(15) =  0.000
    coordinates(16) =  1.675       ! x O on W2
    coordinates(17) = -2.803
    coordinates(18) =  0.911

    cell(1) = 35.000
    cell(2) = 35.000
    cell(3) = 35.000

    forces(:) = 0.0

    ! Call the scme function.
    do i=1,1000
       call scme_calculate(n_atoms, coordinates, cell, forces, u_tot)
    end do

  end subroutine test_scme_perf

  subroutine test_scme_cluster_2_perf()
    implicit none
    ! ----------------------------------------
    integer, parameter :: n_atoms = 6
    integer :: i
    real*8  :: u_tot = 0.0d0
    real*8  :: u_tot_ref = 0.0d0
    real*8,dimension(18)     :: coordinates
    real*8,dimension(3)         :: cell
    real*8,dimension(54)     :: forces
    real*8,dimension(n_atoms*3) :: ref_forces
    ! ----------------------------------------

    ! ---------------------------------------------------------
    ! Test 1) Looping calculation of forces and energy to
    !         run timings.
    ! ---------------------------------------------------------

    coordinates(1)  =  0.018       ! x H1 on W1
    coordinates(2)  = -0.739       ! y H1 on W1
    coordinates(3)  =  0.521       ! z H1 on W1
    coordinates(4)  = -0.815       ! x H2 on W1
    coordinates(5)  = -0.673
    coordinates(6)  = -0.592
    coordinates(7)  =  1.942       ! x H1 on W2
    coordinates(8)  = -1.964       ! y H1 on W2
    coordinates(9)  =  1.223       ! z H1 on W2
    coordinates(10) =  2.319       ! x H2 on W2
    coordinates(11) = -3.450
    coordinates(12) =  1.136
    coordinates(13) =  0.000       ! x O on W1
    coordinates(14) =  0.000
    coordinates(15) =  0.000
    coordinates(16) =  1.675       ! x O on W2
    coordinates(17) = -2.803
    coordinates(18) =  0.911

    cell(1) = 35.000
    cell(2) = 35.000
    cell(3) = 35.000

    forces(:) = 0.0

    ! Call the scme function.
    do i=1,1
       call scme_calculate(n_atoms, coordinates, cell, forces, u_tot)
    end do

  end subroutine 

  subroutine test_scme_cluster_6_perf()
    implicit none
    ! ----------------------------------------
    integer, parameter :: n_atoms = 18 !from 6 waters
    integer :: i
    real*8  :: u_tot = 0.0d0
    real*8  :: u_tot_ref = 0.0d0
    real*8,dimension(3,n_atoms) :: coordinates
    real*8,dimension(3)         :: cell
    real*8,dimension(3,n_atoms) :: forces
    real*8,dimension(3,n_atoms) :: ref_forces
    real*8,dimension(3,n_atoms/3) :: dip_perm,dip_ind
    integer p
    ! ----------------------------------------
!this was for the old version of the ra,fa arrays
!    p=1
!           coordinates(p:p+2) = [-2.006698000d+00 , -4.223270000d-01,  2.219847000d+00]!H1w1
!    p=p+3; coordinates(p:p+2) = [-6.010540000d-01 , -5.969720000d-01,  1.553718000d+00]!H2w1
!    p=p+3; coordinates(p:p+2) = [-2.516835000d+00 , -7.667650000d-01, -1.733766000d+00]!H1w2
!    p=p+3; coordinates(p:p+2) = [-1.888941000d+00 , -4.796530000d-01, -3.476240000d-01]!H2w2
!    p=p+3; coordinates(p:p+2) = [-9.898310000d-01 ,  1.592736000d+00, -8.774190000d-01]!H1w3
!    p=p+3; coordinates(p:p+2) = [-9.477200000d-01 ,  1.533567000d+00,  6.252280000d-01]!H2w3
!    p=p+3; coordinates(p:p+2) = [ 1.542224000d+00 , -3.936920000d-01,  1.344373000d+00]!H1w4
!    p=p+3; coordinates(p:p+2) = [ 9.795570000d-01 , -1.522041000d+00,  5.278330000d-01]!H2w4
!    p=p+3; coordinates(p:p+2) = [ 1.470709000d+00 , -5.709330000d-01, -1.277710000d+00]!H1w5
!    p=p+3; coordinates(p:p+2) = [ 6.516100000d-02 , -1.118951000d+00, -1.522886000d+00]!H2w5
!    p=p+3; coordinates(p:p+2) = [ 2.674716000d+00 ,  1.735342000d+00, -2.379950000d-01]!H1w6
!    p=p+3; coordinates(p:p+2) = [ 1.141637000d+00 ,  1.532266000d+00, -1.401210000d-01]!H2w6
!    p=p+3; coordinates(p:p+2) = [-1.502169000d+00 , -1.913590000d-01,  1.434927000d+00]!O w1
!    p=p+3; coordinates(p:p+2) = [-1.744575000d+00 , -3.823480000d-01, -1.309144000d+00]!O w2
!    p=p+3; coordinates(p:p+2) = [-5.604090000d-01 ,  2.017830000d+00, -1.219840000d-01]!O w3
!    p=p+3; coordinates(p:p+2) = [ 9.648030000d-01 , -1.165765000d+00,  1.439987000d+00]!O w4
!    p=p+3; coordinates(p:p+2) = [ 9.747050000d-01 , -1.401503000d+00, -1.335970000d+00]!O w5
!    p=p+3; coordinates(p:p+2) = [ 2.002280000d+00 ,  1.057824000d+00, -1.245020000d-01]!O w6

    p=1
           coordinates(:,p) = [-2.006698000d+00 , -4.223270000d-01,  2.219847000d+00]!H1w1
    p=p+1; coordinates(:,p) = [-6.010540000d-01 , -5.969720000d-01,  1.553718000d+00]!H2w1
    p=p+1; coordinates(:,p) = [-1.502169000d+00 , -1.913590000d-01,  1.434927000d+00]!O w1
    p=p+1; coordinates(:,p) = [-2.516835000d+00 , -7.667650000d-01, -1.733766000d+00]!H1w2
    p=p+1; coordinates(:,p) = [-1.888941000d+00 , -4.796530000d-01, -3.476240000d-01]!H2w2
    p=p+1; coordinates(:,p) = [-1.744575000d+00 , -3.823480000d-01, -1.309144000d+00]!O w2
    p=p+1; coordinates(:,p) = [-9.898310000d-01 ,  1.592736000d+00, -8.774190000d-01]!H1w3
    p=p+1; coordinates(:,p) = [-9.477200000d-01 ,  1.533567000d+00,  6.252280000d-01]!H2w3
    p=p+1; coordinates(:,p) = [-5.604090000d-01 ,  2.017830000d+00, -1.219840000d-01]!O w3
    p=p+1; coordinates(:,p) = [ 1.542224000d+00 , -3.936920000d-01,  1.344373000d+00]!H1w4
    p=p+1; coordinates(:,p) = [ 9.795570000d-01 , -1.522041000d+00,  5.278330000d-01]!H2w4
    p=p+1; coordinates(:,p) = [ 9.648030000d-01 , -1.165765000d+00,  1.439987000d+00]!O w4
    p=p+1; coordinates(:,p) = [ 1.470709000d+00 , -5.709330000d-01, -1.277710000d+00]!H1w5
    p=p+1; coordinates(:,p) = [ 6.516100000d-02 , -1.118951000d+00, -1.522886000d+00]!H2w5
    p=p+1; coordinates(:,p) = [ 9.747050000d-01 , -1.401503000d+00, -1.335970000d+00]!O w5
    p=p+1; coordinates(:,p) = [ 2.674716000d+00 ,  1.735342000d+00, -2.379950000d-01]!H1w6
    p=p+1; coordinates(:,p) = [ 1.141637000d+00 ,  1.532266000d+00, -1.401210000d-01]!H2w6
    p=p+1; coordinates(:,p) = [ 2.002280000d+00 ,  1.057824000d+00, -1.245020000d-01]!O w6

    print*, "# coordinates should be ",6*3*3,". It is: ", p,'atoms*3=',p*3, 'coordinates'
    cell(1) = 35.000
    cell(2) = 35.000
    cell(3) = 35.000

    forces = 0

    ! Call the scme function.
    do i=1,1
       call scme_calculate(n_atoms, coordinates, cell, forces, u_tot)!,dip_perm,dip_ind)
    end do
    
!    call printer(coordinates,'coordinates',2,.false.)
!    call printer(forces,'forces',2,.false.)
!    call printer(u_tot,'u_tot',2,.false.)
!    
!    call printer(dip_perm*A_a0*au_to_debye,'dip_perm',2,.false.)
!    call printer(dip_ind *A_a0*au_to_debye,'dip_ind',2,.false.)
    
  end subroutine 



  subroutine test_scme_cluster_9_perf()
    implicit none
    ! ----------------------------------------
    integer, parameter :: n_atoms = 27 !from 9 waters
    integer :: i
    real*8  :: u_tot = 0.0d0
    real*8  :: u_tot_ref = 0.0d0
    real*8,dimension(n_atoms*3) :: coordinates
    real*8,dimension(3)         :: cell
    real*8,dimension(n_atoms*3) :: forces
    real*8,dimension(n_atoms*3) :: ref_forces
    integer p
    ! ----------------------------------------

    ! ---------------------------------------------------------
    ! Test 1) Looping calculation of forces and energy to
    !         run timings.
    ! ---------------------------------------------------------
    p=1
           coordinates(p:p+2) = [ -1.4736, -1.3180, -1.6898 ]!H
    p=p+3; coordinates(p:p+2) = [ -0.8686,  0.0696, -1.2669 ]!H
    p=p+3; coordinates(p:p+2) = [ -1.0907, -0.8252, -0.9284 ]!O
    p=p+3; coordinates(p:p+2) = [ -0.7162,  3.3599,  0.8471 ]!H
    p=p+3; coordinates(p:p+2) = [ -0.8391,  1.8572,  1.0804 ]!H
    p=p+3; coordinates(p:p+2) = [ -0.1530,  2.5613,  0.8758 ]!O
    p=p+3; coordinates(p:p+2) = [ -0.3571,  2.1947, -1.5612 ]!H
    p=p+3; coordinates(p:p+2) = [ -1.4623,  2.0190, -2.5554 ]!H
    p=p+3; coordinates(p:p+2) = [ -0.6756,  1.5259, -2.2016 ]!O
    p=p+3; coordinates(p:p+2) = [ -1.8516, -0.1730,  0.5943 ]!H
    p=p+3; coordinates(p:p+2) = [ -2.6503,  0.1754,  1.8833 ]!H
    p=p+3; coordinates(p:p+2) = [ -2.1394,  0.5531,  1.1850 ]!O
    p=p+3; coordinates(p:p+2) = [ -1.6640,  4.9992, -0.1046 ]!H
    p=p+3; coordinates(p:p+2) = [ -2.1483,  5.5160,  1.3337 ]!H
    p=p+3; coordinates(p:p+2) = [ -1.9288,  4.7626,  0.8354 ]!O
    p=p+3; coordinates(p:p+2) = [ -3.1531,  1.9088,  0.2640 ]!H
    p=p+3; coordinates(p:p+2) = [ -3.2814,  3.4999,  0.2831 ]!H
    p=p+3; coordinates(p:p+2) = [ -3.6595,  2.6702, -0.1104 ]!O
    p=p+3; coordinates(p:p+2) = [  1.5809,  4.9571, -1.0815 ]!H
    p=p+3; coordinates(p:p+2) = [  0.9649,  3.6761, -0.4809 ]!H
    p=p+3; coordinates(p:p+2) = [  1.1160,  4.1296, -1.3211 ]!O
    p=p+3; coordinates(p:p+2) = [ -3.3448,  2.9950, -1.7281 ]!H
    p=p+3; coordinates(p:p+2) = [ -3.7567,  3.1501, -3.2142 ]!H
    p=p+3; coordinates(p:p+2) = [ -2.9721,  3.2005, -2.6291 ]!O
    p=p+3; coordinates(p:p+2) = [ -0.4045,  4.8032, -1.9164 ]!H
    p=p+3; coordinates(p:p+2) = [ -1.8429,  4.6655, -2.4245 ]!H
    p=p+3; coordinates(p:p+2) = [ -1.2805,  5.2116, -1.8113 ]!O
    
    print*, "things in coordinates should be 81 and it is: ", p+2
    cell(1) = 35.000
    cell(2) = 35.000
    cell(3) = 35.000

    forces(:) = 0.0

    ! Call the scme function.
    do i=1,1
       call scme_calculate(n_atoms, coordinates, cell, forces, u_tot)
    end do
    
    print*, "test_scme_cluster_9_perf() energy is", u_tot
    
    


  end subroutine 

  subroutine test_scme_cluster_15_perf()
    implicit none
    ! ----------------------------------------
    integer, parameter :: n_atoms = 45
    integer :: i
    real*8  :: u_tot = 0.0d0
    real*8  :: u_tot_ref = 0.0d0
    real*8,dimension(n_atoms*3) :: coordinates
    real*8,dimension(3)         :: cell
    real*8,dimension(n_atoms*3) :: forces
    real*8,dimension(n_atoms*3) :: ref_forces
    integer p
    ! ----------------------------------------

    ! ---------------------------------------------------------
    ! Test 1) Looping calculation of forces and energy to
    !         run timings.
    ! ---------------------------------------------------------
    p=1 !the data:
           coordinates(p:p+2) = [ 1.80482, -1.73505, -0.85955]!H
    p=p+3; coordinates(p:p+2) = [ 2.42260, -0.35826, -0.70808]!H
    p=p+3; coordinates(p:p+2) = [-0.04077, -0.15230, -0.04336]!H
    p=p+3; coordinates(p:p+2) = [-0.01324,  0.97891,  0.98093]!H
    p=p+3; coordinates(p:p+2) = [ 4.24637,  1.27773, -0.86698]!H
    p=p+3; coordinates(p:p+2) = [ 3.83542,  0.83590, -2.20574]!H
    p=p+3; coordinates(p:p+2) = [-1.33710,  4.11792, -0.55783]!H
    p=p+3; coordinates(p:p+2) = [-0.15858,  4.05926, -1.54358]!H
    p=p+3; coordinates(p:p+2) = [ 3.70433,  3.38958,  2.75945]!H
    p=p+3; coordinates(p:p+2) = [ 2.87922,  3.94545,  1.55536]!H
    p=p+3; coordinates(p:p+2) = [ 0.48484,  1.78382, -0.82579]!H
    p=p+3; coordinates(p:p+2) = [ 1.83919,  1.86647, -1.67688]!H
    p=p+3; coordinates(p:p+2) = [ 1.36128,  7.31484, -0.25489]!H
    p=p+3; coordinates(p:p+2) = [ 0.50437,  5.95037, -0.40769]!H
    p=p+3; coordinates(p:p+2) = [ 2.15403,  0.94247,  2.28060]!H
    p=p+3; coordinates(p:p+2) = [ 2.89100, -0.27132,  1.62148]!H
    p=p+3; coordinates(p:p+2) = [ 1.30408,  2.61876,  2.42595]!H
    p=p+3; coordinates(p:p+2) = [-0.16326,  2.47088,  2.68855]!H
    p=p+3; coordinates(p:p+2) = [ 2.63809,  5.08587, -0.37526]!H
    p=p+3; coordinates(p:p+2) = [ 2.24243,  3.64913, -0.71816]!H
    p=p+3; coordinates(p:p+2) = [-2.24224,  3.03902,  1.19737]!H
    p=p+3; coordinates(p:p+2) = [-1.71070,  2.17173,  0.06858]!H
    p=p+3; coordinates(p:p+2) = [ 0.72546,  5.83903,  1.80321]!H
    p=p+3; coordinates(p:p+2) = [ 1.40607,  5.02224,  3.03626]!H
    p=p+3; coordinates(p:p+2) = [ 4.51732,  1.44236,  2.67606]!H
    p=p+3; coordinates(p:p+2) = [ 5.39808,  2.08574,  3.73668]!H
    p=p+3; coordinates(p:p+2) = [ 5.55434,  2.62557,  0.92644]!H
    p=p+3; coordinates(p:p+2) = [ 4.66419,  3.33488, -0.23372]!H
    p=p+3; coordinates(p:p+2) = [-1.45378,  3.50547,  3.98412]!H
    p=p+3; coordinates(p:p+2) = [-0.77806,  4.38400,  2.89675]!H
    p=p+3; coordinates(p:p+2) = [ 1.66322, -0.92324, -0.35186]!O
    p=p+3; coordinates(p:p+2) = [-0.45776,  0.78681,  0.11944]!O
    p=p+3; coordinates(p:p+2) = [ 3.53951,  0.77639, -1.32046]!O
    p=p+3; coordinates(p:p+2) = [-0.57027,  4.66525, -0.91234]!O
    p=p+3; coordinates(p:p+2) = [ 2.82467,  3.80518,  2.50783]!O
    p=p+3; coordinates(p:p+2) = [ 1.05912,  2.41662, -1.33542]!O
    p=p+3; coordinates(p:p+2) = [ 1.24000,  6.41845,  0.03330]!O
    p=p+3; coordinates(p:p+2) = [ 2.98116,  0.48514,  2.18258]!O
    p=p+3; coordinates(p:p+2) = [ 0.60944,  1.93689,  2.41997]!O
    p=p+3; coordinates(p:p+2) = [ 2.96096,  4.17215, -0.24367]!O
    p=p+3; coordinates(p:p+2) = [-2.15967,  3.04991,  0.22393]!O
    p=p+3; coordinates(p:p+2) = [ 0.56317,  5.47605,  2.72575]!O
    p=p+3; coordinates(p:p+2) = [ 5.08862,  2.23134,  2.85097]!O
    p=p+3; coordinates(p:p+2) = [ 5.20638,  2.53755, -0.01077]!O
    p=p+3; coordinates(p:p+2) = [-1.42229,  3.63362,  3.02443]!O
    print*, "things in coordinates should be",15*3*3, "and it is: ", p+2
    
    cell(1) = 35.000
    cell(2) = 35.000
    cell(3) = 35.000

    forces(:) = 0.0

    ! Call the scme function.
    call scme_calculate(n_atoms, coordinates, cell, forces, u_tot)
    
  end subroutine 
  
  subroutine test_scme_cluster_75_perf()
    implicit none
    ! ----------------------------------------
    integer, parameter :: n_atoms = 15*5*3 !15cluster 4*repeted 3*atomsPerWater
    integer :: i
    real*8  :: u_tot = 0.0d0
    real*8  :: u_tot_ref = 0.0d0
    real*8,dimension(n_atoms*3) :: coordinates
    real*8,dimension(3)         :: cell
    real*8,dimension(n_atoms*3) :: forces
    real*8,dimension(n_atoms*3) :: ref_forces
    integer p
    ! ----------------------------------------

    ! ---------------------------------------------------------
    ! In this subroutine I repeated the 15 cluster in the 
    ! (0),x,y,z direction. Can still expand in xy,xz,yz,xyz
    ! to make it larger
    ! ---------------------------------------------------------
    p=1 !the data:
           coordinates(p:p+2) = [ 1.80482, -1.73505, -0.85955]!H
    p=p+3; coordinates(p:p+2) = [ 2.42260, -0.35826, -0.70808]!H
    p=p+3; coordinates(p:p+2) = [-0.04077, -0.15230, -0.04336]!H
    p=p+3; coordinates(p:p+2) = [-0.01324,  0.97891,  0.98093]!H
    p=p+3; coordinates(p:p+2) = [ 4.24637,  1.27773, -0.86698]!H
    p=p+3; coordinates(p:p+2) = [ 3.83542,  0.83590, -2.20574]!H
    p=p+3; coordinates(p:p+2) = [-1.33710,  4.11792, -0.55783]!H
    p=p+3; coordinates(p:p+2) = [-0.15858,  4.05926, -1.54358]!H
    p=p+3; coordinates(p:p+2) = [ 3.70433,  3.38958,  2.75945]!H
    p=p+3; coordinates(p:p+2) = [ 2.87922,  3.94545,  1.55536]!H
    p=p+3; coordinates(p:p+2) = [ 0.48484,  1.78382, -0.82579]!H
    p=p+3; coordinates(p:p+2) = [ 1.83919,  1.86647, -1.67688]!H
    p=p+3; coordinates(p:p+2) = [ 1.36128,  7.31484, -0.25489]!H
    p=p+3; coordinates(p:p+2) = [ 0.50437,  5.95037, -0.40769]!H
    p=p+3; coordinates(p:p+2) = [ 2.15403,  0.94247,  2.28060]!H
    p=p+3; coordinates(p:p+2) = [ 2.89100, -0.27132,  1.62148]!H
    p=p+3; coordinates(p:p+2) = [ 1.30408,  2.61876,  2.42595]!H
    p=p+3; coordinates(p:p+2) = [-0.16326,  2.47088,  2.68855]!H
    p=p+3; coordinates(p:p+2) = [ 2.63809,  5.08587, -0.37526]!H
    p=p+3; coordinates(p:p+2) = [ 2.24243,  3.64913, -0.71816]!H
    p=p+3; coordinates(p:p+2) = [-2.24224,  3.03902,  1.19737]!H
    p=p+3; coordinates(p:p+2) = [-1.71070,  2.17173,  0.06858]!H
    p=p+3; coordinates(p:p+2) = [ 0.72546,  5.83903,  1.80321]!H
    p=p+3; coordinates(p:p+2) = [ 1.40607,  5.02224,  3.03626]!H
    p=p+3; coordinates(p:p+2) = [ 4.51732,  1.44236,  2.67606]!H
    p=p+3; coordinates(p:p+2) = [ 5.39808,  2.08574,  3.73668]!H
    p=p+3; coordinates(p:p+2) = [ 5.55434,  2.62557,  0.92644]!H
    p=p+3; coordinates(p:p+2) = [ 4.66419,  3.33488, -0.23372]!H
    p=p+3; coordinates(p:p+2) = [-1.45378,  3.50547,  3.98412]!H
    p=p+3; coordinates(p:p+2) = [-0.77806,  4.38400,  2.89675]!H
    
    p=p+3; coordinates(p:p+2) = [( 1.80482+7.0), -1.73505, -0.85955]!H
    p=p+3; coordinates(p:p+2) = [( 2.42260+7.0), -0.35826, -0.70808]!H
    p=p+3; coordinates(p:p+2) = [(-0.04077+7.0), -0.15230, -0.04336]!H
    p=p+3; coordinates(p:p+2) = [(-0.01324+7.0),  0.97891,  0.98093]!H
    p=p+3; coordinates(p:p+2) = [( 4.24637+7.0),  1.27773, -0.86698]!H
    p=p+3; coordinates(p:p+2) = [( 3.83542+7.0),  0.83590, -2.20574]!H
    p=p+3; coordinates(p:p+2) = [(-1.33710+7.0),  4.11792, -0.55783]!H
    p=p+3; coordinates(p:p+2) = [(-0.15858+7.0),  4.05926, -1.54358]!H
    p=p+3; coordinates(p:p+2) = [( 3.70433+7.0),  3.38958,  2.75945]!H
    p=p+3; coordinates(p:p+2) = [( 2.87922+7.0),  3.94545,  1.55536]!H
    p=p+3; coordinates(p:p+2) = [( 0.48484+7.0),  1.78382, -0.82579]!H
    p=p+3; coordinates(p:p+2) = [( 1.83919+7.0),  1.86647, -1.67688]!H
    p=p+3; coordinates(p:p+2) = [( 1.36128+7.0),  7.31484, -0.25489]!H
    p=p+3; coordinates(p:p+2) = [( 0.50437+7.0),  5.95037, -0.40769]!H
    p=p+3; coordinates(p:p+2) = [( 2.15403+7.0),  0.94247,  2.28060]!H
    p=p+3; coordinates(p:p+2) = [( 2.89100+7.0), -0.27132,  1.62148]!H
    p=p+3; coordinates(p:p+2) = [( 1.30408+7.0),  2.61876,  2.42595]!H
    p=p+3; coordinates(p:p+2) = [(-0.16326+7.0),  2.47088,  2.68855]!H
    p=p+3; coordinates(p:p+2) = [( 2.63809+7.0),  5.08587, -0.37526]!H
    p=p+3; coordinates(p:p+2) = [( 2.24243+7.0),  3.64913, -0.71816]!H
    p=p+3; coordinates(p:p+2) = [(-2.24224+7.0),  3.03902,  1.19737]!H
    p=p+3; coordinates(p:p+2) = [(-1.71070+7.0),  2.17173,  0.06858]!H
    p=p+3; coordinates(p:p+2) = [( 0.72546+7.0),  5.83903,  1.80321]!H
    p=p+3; coordinates(p:p+2) = [( 1.40607+7.0),  5.02224,  3.03626]!H
    p=p+3; coordinates(p:p+2) = [( 4.51732+7.0),  1.44236,  2.67606]!H
    p=p+3; coordinates(p:p+2) = [( 5.39808+7.0),  2.08574,  3.73668]!H
    p=p+3; coordinates(p:p+2) = [( 5.55434+7.0),  2.62557,  0.92644]!H
    p=p+3; coordinates(p:p+2) = [( 4.66419+7.0),  3.33488, -0.23372]!H
    p=p+3; coordinates(p:p+2) = [(-1.45378+7.0),  3.50547,  3.98412]!H
    p=p+3; coordinates(p:p+2) = [(-0.77806+7.0),  4.38400,  2.89675]!H

    p=p+3; coordinates(p:p+2) = [ 1.80482, (-1.73505+7.0), -0.85955]!H
    p=p+3; coordinates(p:p+2) = [ 2.42260, (-0.35826+7.0), -0.70808]!H
    p=p+3; coordinates(p:p+2) = [-0.04077, (-0.15230+7.0), -0.04336]!H
    p=p+3; coordinates(p:p+2) = [-0.01324, ( 0.97891+7.0),  0.98093]!H
    p=p+3; coordinates(p:p+2) = [ 4.24637, ( 1.27773+7.0), -0.86698]!H
    p=p+3; coordinates(p:p+2) = [ 3.83542, ( 0.83590+7.0), -2.20574]!H
    p=p+3; coordinates(p:p+2) = [-1.33710, ( 4.11792+7.0), -0.55783]!H
    p=p+3; coordinates(p:p+2) = [-0.15858, ( 4.05926+7.0), -1.54358]!H
    p=p+3; coordinates(p:p+2) = [ 3.70433, ( 3.38958+7.0),  2.75945]!H
    p=p+3; coordinates(p:p+2) = [ 2.87922, ( 3.94545+7.0),  1.55536]!H
    p=p+3; coordinates(p:p+2) = [ 0.48484, ( 1.78382+7.0), -0.82579]!H
    p=p+3; coordinates(p:p+2) = [ 1.83919, ( 1.86647+7.0), -1.67688]!H
    p=p+3; coordinates(p:p+2) = [ 1.36128, ( 7.31484+7.0), -0.25489]!H
    p=p+3; coordinates(p:p+2) = [ 0.50437, ( 5.95037+7.0), -0.40769]!H
    p=p+3; coordinates(p:p+2) = [ 2.15403, ( 0.94247+7.0),  2.28060]!H
    p=p+3; coordinates(p:p+2) = [ 2.89100, (-0.27132+7.0),  1.62148]!H
    p=p+3; coordinates(p:p+2) = [ 1.30408, ( 2.61876+7.0),  2.42595]!H
    p=p+3; coordinates(p:p+2) = [-0.16326, ( 2.47088+7.0),  2.68855]!H
    p=p+3; coordinates(p:p+2) = [ 2.63809, ( 5.08587+7.0), -0.37526]!H
    p=p+3; coordinates(p:p+2) = [ 2.24243, ( 3.64913+7.0), -0.71816]!H
    p=p+3; coordinates(p:p+2) = [-2.24224, ( 3.03902+7.0),  1.19737]!H
    p=p+3; coordinates(p:p+2) = [-1.71070, ( 2.17173+7.0),  0.06858]!H
    p=p+3; coordinates(p:p+2) = [ 0.72546, ( 5.83903+7.0),  1.80321]!H
    p=p+3; coordinates(p:p+2) = [ 1.40607, ( 5.02224+7.0),  3.03626]!H
    p=p+3; coordinates(p:p+2) = [ 4.51732, ( 1.44236+7.0),  2.67606]!H
    p=p+3; coordinates(p:p+2) = [ 5.39808, ( 2.08574+7.0),  3.73668]!H
    p=p+3; coordinates(p:p+2) = [ 5.55434, ( 2.62557+7.0),  0.92644]!H
    p=p+3; coordinates(p:p+2) = [ 4.66419, ( 3.33488+7.0), -0.23372]!H
    p=p+3; coordinates(p:p+2) = [-1.45378, ( 3.50547+7.0),  3.98412]!H
    p=p+3; coordinates(p:p+2) = [-0.77806, ( 4.38400+7.0),  2.89675]!H

    p=p+3; coordinates(p:p+2) = [ 1.80482, -1.73505, (-0.85955+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 2.42260, -0.35826, (-0.70808+7.0)]!H
    p=p+3; coordinates(p:p+2) = [-0.04077, -0.15230, (-0.04336+7.0)]!H
    p=p+3; coordinates(p:p+2) = [-0.01324,  0.97891, ( 0.98093+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 4.24637,  1.27773, (-0.86698+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 3.83542,  0.83590, (-2.20574+7.0)]!H
    p=p+3; coordinates(p:p+2) = [-1.33710,  4.11792, (-0.55783+7.0)]!H
    p=p+3; coordinates(p:p+2) = [-0.15858,  4.05926, (-1.54358+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 3.70433,  3.38958, ( 2.75945+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 2.87922,  3.94545, ( 1.55536+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 0.48484,  1.78382, (-0.82579+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 1.83919,  1.86647, (-1.67688+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 1.36128,  7.31484, (-0.25489+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 0.50437,  5.95037, (-0.40769+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 2.15403,  0.94247, ( 2.28060+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 2.89100, -0.27132, ( 1.62148+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 1.30408,  2.61876, ( 2.42595+7.0)]!H
    p=p+3; coordinates(p:p+2) = [-0.16326,  2.47088, ( 2.68855+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 2.63809,  5.08587, (-0.37526+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 2.24243,  3.64913, (-0.71816+7.0)]!H
    p=p+3; coordinates(p:p+2) = [-2.24224,  3.03902, ( 1.19737+7.0)]!H
    p=p+3; coordinates(p:p+2) = [-1.71070,  2.17173, ( 0.06858+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 0.72546,  5.83903, ( 1.80321+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 1.40607,  5.02224, ( 3.03626+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 4.51732,  1.44236, ( 2.67606+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 5.39808,  2.08574, ( 3.73668+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 5.55434,  2.62557, ( 0.92644+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 4.66419,  3.33488, (-0.23372+7.0)]!H
    p=p+3; coordinates(p:p+2) = [-1.45378,  3.50547, ( 3.98412+7.0)]!H
    p=p+3; coordinates(p:p+2) = [-0.77806,  4.38400, ( 2.89675+7.0)]!H

    p=p+3; coordinates(p:p+2) = [ 1.80482, (-1.73505+7.0), (-0.85955+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 2.42260, (-0.35826+7.0), (-0.70808+7.0)]!H
    p=p+3; coordinates(p:p+2) = [-0.04077, (-0.15230+7.0), (-0.04336+7.0)]!H
    p=p+3; coordinates(p:p+2) = [-0.01324, ( 0.97891+7.0), ( 0.98093+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 4.24637, ( 1.27773+7.0), (-0.86698+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 3.83542, ( 0.83590+7.0), (-2.20574+7.0)]!H
    p=p+3; coordinates(p:p+2) = [-1.33710, ( 4.11792+7.0), (-0.55783+7.0)]!H
    p=p+3; coordinates(p:p+2) = [-0.15858, ( 4.05926+7.0), (-1.54358+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 3.70433, ( 3.38958+7.0), ( 2.75945+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 2.87922, ( 3.94545+7.0), ( 1.55536+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 0.48484, ( 1.78382+7.0), (-0.82579+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 1.83919, ( 1.86647+7.0), (-1.67688+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 1.36128, ( 7.31484+7.0), (-0.25489+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 0.50437, ( 5.95037+7.0), (-0.40769+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 2.15403, ( 0.94247+7.0), ( 2.28060+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 2.89100, (-0.27132+7.0), ( 1.62148+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 1.30408, ( 2.61876+7.0), ( 2.42595+7.0)]!H
    p=p+3; coordinates(p:p+2) = [-0.16326, ( 2.47088+7.0), ( 2.68855+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 2.63809, ( 5.08587+7.0), (-0.37526+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 2.24243, ( 3.64913+7.0), (-0.71816+7.0)]!H
    p=p+3; coordinates(p:p+2) = [-2.24224, ( 3.03902+7.0), ( 1.19737+7.0)]!H
    p=p+3; coordinates(p:p+2) = [-1.71070, ( 2.17173+7.0), ( 0.06858+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 0.72546, ( 5.83903+7.0), ( 1.80321+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 1.40607, ( 5.02224+7.0), ( 3.03626+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 4.51732, ( 1.44236+7.0), ( 2.67606+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 5.39808, ( 2.08574+7.0), ( 3.73668+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 5.55434, ( 2.62557+7.0), ( 0.92644+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 4.66419, ( 3.33488+7.0), (-0.23372+7.0)]!H
    p=p+3; coordinates(p:p+2) = [-1.45378, ( 3.50547+7.0), ( 3.98412+7.0)]!H
    p=p+3; coordinates(p:p+2) = [-0.77806, ( 4.38400+7.0), ( 2.89675+7.0)]!H    
    
    p=p+3; coordinates(p:p+2) = [ 1.66322, -0.92324, -0.35186]!O 0
    p=p+3; coordinates(p:p+2) = [-0.45776,  0.78681,  0.11944]!O
    p=p+3; coordinates(p:p+2) = [ 3.53951,  0.77639, -1.32046]!O
    p=p+3; coordinates(p:p+2) = [-0.57027,  4.66525, -0.91234]!O
    p=p+3; coordinates(p:p+2) = [ 2.82467,  3.80518,  2.50783]!O
    p=p+3; coordinates(p:p+2) = [ 1.05912,  2.41662, -1.33542]!O
    p=p+3; coordinates(p:p+2) = [ 1.24000,  6.41845,  0.03330]!O
    p=p+3; coordinates(p:p+2) = [ 2.98116,  0.48514,  2.18258]!O
    p=p+3; coordinates(p:p+2) = [ 0.60944,  1.93689,  2.41997]!O
    p=p+3; coordinates(p:p+2) = [ 2.96096,  4.17215, -0.24367]!O
    p=p+3; coordinates(p:p+2) = [-2.15967,  3.04991,  0.22393]!O
    p=p+3; coordinates(p:p+2) = [ 0.56317,  5.47605,  2.72575]!O
    p=p+3; coordinates(p:p+2) = [ 5.08862,  2.23134,  2.85097]!O
    p=p+3; coordinates(p:p+2) = [ 5.20638,  2.53755, -0.01077]!O
    p=p+3; coordinates(p:p+2) = [-1.42229,  3.63362,  3.02443]!O

    p=p+3; coordinates(p:p+2) = [( 1.66322+7.0), -0.92324, -0.35186]!O x
    p=p+3; coordinates(p:p+2) = [(-0.45776+7.0),  0.78681,  0.11944]!O
    p=p+3; coordinates(p:p+2) = [( 3.53951+7.0),  0.77639, -1.32046]!O
    p=p+3; coordinates(p:p+2) = [(-0.57027+7.0),  4.66525, -0.91234]!O
    p=p+3; coordinates(p:p+2) = [( 2.82467+7.0),  3.80518,  2.50783]!O
    p=p+3; coordinates(p:p+2) = [( 1.05912+7.0),  2.41662, -1.33542]!O
    p=p+3; coordinates(p:p+2) = [( 1.24000+7.0),  6.41845,  0.03330]!O
    p=p+3; coordinates(p:p+2) = [( 2.98116+7.0),  0.48514,  2.18258]!O
    p=p+3; coordinates(p:p+2) = [( 0.60944+7.0),  1.93689,  2.41997]!O
    p=p+3; coordinates(p:p+2) = [( 2.96096+7.0),  4.17215, -0.24367]!O
    p=p+3; coordinates(p:p+2) = [(-2.15967+7.0),  3.04991,  0.22393]!O
    p=p+3; coordinates(p:p+2) = [( 0.56317+7.0),  5.47605,  2.72575]!O
    p=p+3; coordinates(p:p+2) = [( 5.08862+7.0),  2.23134,  2.85097]!O
    p=p+3; coordinates(p:p+2) = [( 5.20638+7.0),  2.53755, -0.01077]!O
    p=p+3; coordinates(p:p+2) = [(-1.42229+7.0),  3.63362,  3.02443]!O


    p=p+3; coordinates(p:p+2) = [ 1.66322, (-0.92324+7.0), -0.35186]!O y
    p=p+3; coordinates(p:p+2) = [-0.45776, ( 0.78681+7.0),  0.11944]!O
    p=p+3; coordinates(p:p+2) = [ 3.53951, ( 0.77639+7.0), -1.32046]!O
    p=p+3; coordinates(p:p+2) = [-0.57027, ( 4.66525+7.0), -0.91234]!O
    p=p+3; coordinates(p:p+2) = [ 2.82467, ( 3.80518+7.0),  2.50783]!O
    p=p+3; coordinates(p:p+2) = [ 1.05912, ( 2.41662+7.0), -1.33542]!O
    p=p+3; coordinates(p:p+2) = [ 1.24000, ( 6.41845+7.0),  0.03330]!O
    p=p+3; coordinates(p:p+2) = [ 2.98116, ( 0.48514+7.0),  2.18258]!O
    p=p+3; coordinates(p:p+2) = [ 0.60944, ( 1.93689+7.0),  2.41997]!O
    p=p+3; coordinates(p:p+2) = [ 2.96096, ( 4.17215+7.0), -0.24367]!O
    p=p+3; coordinates(p:p+2) = [-2.15967, ( 3.04991+7.0),  0.22393]!O
    p=p+3; coordinates(p:p+2) = [ 0.56317, ( 5.47605+7.0),  2.72575]!O
    p=p+3; coordinates(p:p+2) = [ 5.08862, ( 2.23134+7.0),  2.85097]!O
    p=p+3; coordinates(p:p+2) = [ 5.20638, ( 2.53755+7.0), -0.01077]!O
    p=p+3; coordinates(p:p+2) = [-1.42229, ( 3.63362+7.0),  3.02443]!O

    p=p+3; coordinates(p:p+2) = [ 1.66322, -0.92324, (-0.35186+7.0)]!O z
    p=p+3; coordinates(p:p+2) = [-0.45776,  0.78681, ( 0.11944+7.0)]!O
    p=p+3; coordinates(p:p+2) = [ 3.53951,  0.77639, (-1.32046+7.0)]!O
    p=p+3; coordinates(p:p+2) = [-0.57027,  4.66525, (-0.91234+7.0)]!O
    p=p+3; coordinates(p:p+2) = [ 2.82467,  3.80518, ( 2.50783+7.0)]!O
    p=p+3; coordinates(p:p+2) = [ 1.05912,  2.41662, (-1.33542+7.0)]!O
    p=p+3; coordinates(p:p+2) = [ 1.24000,  6.41845, ( 0.03330+7.0)]!O
    p=p+3; coordinates(p:p+2) = [ 2.98116,  0.48514, ( 2.18258+7.0)]!O
    p=p+3; coordinates(p:p+2) = [ 0.60944,  1.93689, ( 2.41997+7.0)]!O
    p=p+3; coordinates(p:p+2) = [ 2.96096,  4.17215, (-0.24367+7.0)]!O
    p=p+3; coordinates(p:p+2) = [-2.15967,  3.04991, ( 0.22393+7.0)]!O
    p=p+3; coordinates(p:p+2) = [ 0.56317,  5.47605, ( 2.72575+7.0)]!O
    p=p+3; coordinates(p:p+2) = [ 5.08862,  2.23134, ( 2.85097+7.0)]!O
    p=p+3; coordinates(p:p+2) = [ 5.20638,  2.53755, (-0.01077+7.0)]!O
    p=p+3; coordinates(p:p+2) = [-1.42229,  3.63362, ( 3.02443+7.0)]!O

    p=p+3; coordinates(p:p+2) = [ 1.66322, (-0.92324+7.0), (-0.35186+7.0)]!O yz
    p=p+3; coordinates(p:p+2) = [-0.45776, ( 0.78681+7.0), ( 0.11944+7.0)]!O
    p=p+3; coordinates(p:p+2) = [ 3.53951, ( 0.77639+7.0), (-1.32046+7.0)]!O
    p=p+3; coordinates(p:p+2) = [-0.57027, ( 4.66525+7.0), (-0.91234+7.0)]!O
    p=p+3; coordinates(p:p+2) = [ 2.82467, ( 3.80518+7.0), ( 2.50783+7.0)]!O
    p=p+3; coordinates(p:p+2) = [ 1.05912, ( 2.41662+7.0), (-1.33542+7.0)]!O
    p=p+3; coordinates(p:p+2) = [ 1.24000, ( 6.41845+7.0), ( 0.03330+7.0)]!O
    p=p+3; coordinates(p:p+2) = [ 2.98116, ( 0.48514+7.0), ( 2.18258+7.0)]!O
    p=p+3; coordinates(p:p+2) = [ 0.60944, ( 1.93689+7.0), ( 2.41997+7.0)]!O
    p=p+3; coordinates(p:p+2) = [ 2.96096, ( 4.17215+7.0), (-0.24367+7.0)]!O
    p=p+3; coordinates(p:p+2) = [-2.15967, ( 3.04991+7.0), ( 0.22393+7.0)]!O
    p=p+3; coordinates(p:p+2) = [ 0.56317, ( 5.47605+7.0), ( 2.72575+7.0)]!O
    p=p+3; coordinates(p:p+2) = [ 5.08862, ( 2.23134+7.0), ( 2.85097+7.0)]!O
    p=p+3; coordinates(p:p+2) = [ 5.20638, ( 2.53755+7.0), (-0.01077+7.0)]!O
    p=p+3; coordinates(p:p+2) = [-1.42229, ( 3.63362+7.0), ( 3.02443+7.0)]!O


    print*, "things in coordinates should be",n_atoms*3, "and it is: ", p+2
    
    cell(1) = 35.000
    cell(2) = 35.000
    cell(3) = 35.000

    forces(:) = 0.0

    ! Call the scme function.
    do i=1,1
       call scme_calculate(n_atoms, coordinates, cell, forces, u_tot)
    end do

  end subroutine 
  

  
!    p=p+3; coordinates(p:p+2) = [ 1.80482, -1.73505, -0.85955]!H
!    p=p+3; coordinates(p:p+2) = [ 2.42260, -0.35826, -0.70808]!H
!    p=p+3; coordinates(p:p+2) = [-0.04077, -0.15230, -0.04336]!H
!    p=p+3; coordinates(p:p+2) = [-0.01324,  0.97891,  0.98093]!H
!    p=p+3; coordinates(p:p+2) = [ 4.24637,  1.27773, -0.86698]!H
!    p=p+3; coordinates(p:p+2) = [ 3.83542,  0.83590, -2.20574]!H
!    p=p+3; coordinates(p:p+2) = [-1.33710,  4.11792, -0.55783]!H
!    p=p+3; coordinates(p:p+2) = [-0.15858,  4.05926, -1.54358]!H
!    p=p+3; coordinates(p:p+2) = [ 3.70433,  3.38958,  2.75945]!H
!    p=p+3; coordinates(p:p+2) = [ 2.87922,  3.94545,  1.55536]!H
!    p=p+3; coordinates(p:p+2) = [ 0.48484,  1.78382, -0.82579]!H
!    p=p+3; coordinates(p:p+2) = [ 1.83919,  1.86647, -1.67688]!H
!    p=p+3; coordinates(p:p+2) = [ 1.36128,  7.31484, -0.25489]!H
!    p=p+3; coordinates(p:p+2) = [ 0.50437,  5.95037, -0.40769]!H
!    p=p+3; coordinates(p:p+2) = [ 2.15403,  0.94247,  2.28060]!H
!    p=p+3; coordinates(p:p+2) = [ 2.89100, -0.27132,  1.62148]!H
!    p=p+3; coordinates(p:p+2) = [ 1.30408,  2.61876,  2.42595]!H
!    p=p+3; coordinates(p:p+2) = [-0.16326,  2.47088,  2.68855]!H
!    p=p+3; coordinates(p:p+2) = [ 2.63809,  5.08587, -0.37526]!H
!    p=p+3; coordinates(p:p+2) = [ 2.24243,  3.64913, -0.71816]!H
!    p=p+3; coordinates(p:p+2) = [-2.24224,  3.03902,  1.19737]!H
!    p=p+3; coordinates(p:p+2) = [-1.71070,  2.17173,  0.06858]!H
!    p=p+3; coordinates(p:p+2) = [ 0.72546,  5.83903,  1.80321]!H
!    p=p+3; coordinates(p:p+2) = [ 1.40607,  5.02224,  3.03626]!H
!    p=p+3; coordinates(p:p+2) = [ 4.51732,  1.44236,  2.67606]!H
!    p=p+3; coordinates(p:p+2) = [ 5.39808,  2.08574,  3.73668]!H
!    p=p+3; coordinates(p:p+2) = [ 5.55434,  2.62557,  0.92644]!H
!    p=p+3; coordinates(p:p+2) = [ 4.66419,  3.33488, -0.23372]!H
!    p=p+3; coordinates(p:p+2) = [-1.45378,  3.50547,  3.98412]!H
!    p=p+3; coordinates(p:p+2) = [-0.77806,  4.38400,  2.89675]!H
!    p=p+3; coordinates(p:p+2) = [ 1.66322, -0.92324, -0.35186]!O
!    p=p+3; coordinates(p:p+2) = [-0.45776,  0.78681,  0.11944]!O
!    p=p+3; coordinates(p:p+2) = [ 3.53951,  0.77639, -1.32046]!O
!    p=p+3; coordinates(p:p+2) = [-0.57027,  4.66525, -0.91234]!O
!    p=p+3; coordinates(p:p+2) = [ 2.82467,  3.80518,  2.50783]!O
!    p=p+3; coordinates(p:p+2) = [ 1.05912,  2.41662, -1.33542]!O
!    p=p+3; coordinates(p:p+2) = [ 1.24000,  6.41845,  0.03330]!O
!    p=p+3; coordinates(p:p+2) = [ 2.98116,  0.48514,  2.18258]!O
!    p=p+3; coordinates(p:p+2) = [ 0.60944,  1.93689,  2.41997]!O
!    p=p+3; coordinates(p:p+2) = [ 2.96096,  4.17215, -0.24367]!O
!    p=p+3; coordinates(p:p+2) = [-2.15967,  3.04991,  0.22393]!O
!    p=p+3; coordinates(p:p+2) = [ 0.56317,  5.47605,  2.72575]!O
!    p=p+3; coordinates(p:p+2) = [ 5.08862,  2.23134,  2.85097]!O
!    p=p+3; coordinates(p:p+2) = [ 5.20638,  2.53755, -0.01077]!O
!    p=p+3; coordinates(p:p+2) = [-1.42229,  3.63362,  3.02443]!O
  
  subroutine test_scme_cluster_120_perf()
    implicit none
    ! ----------------------------------------
    integer, parameter :: n_atoms = 15*8*3 !15cluster 4*repeted 3*atomsPerWater
    integer :: i
    real*8  :: u_tot = 0.0d0
    real*8  :: u_tot_ref = 0.0d0
    real*8,dimension(n_atoms*3) :: coordinates
    real*8,dimension(3)         :: cell
    real*8,dimension(n_atoms*3) :: forces
    real*8,dimension(n_atoms*3) :: ref_forces
    integer p
    ! ----------------------------------------

    ! ---------------------------------------------------------
    ! In this subroutine I repeated the 15 cluster in the 
    ! (0),x,y,z direction. Can still expand in xy,xz,yz,xyz
    ! to make it larger
    ! ---------------------------------------------------------
    p=1 !the data:
           coordinates(p:p+2) = [ 1.80482, -1.73505, -0.85955]!H 0
    p=p+3; coordinates(p:p+2) = [ 2.42260, -0.35826, -0.70808]!H
    p=p+3; coordinates(p:p+2) = [-0.04077, -0.15230, -0.04336]!H
    p=p+3; coordinates(p:p+2) = [-0.01324,  0.97891,  0.98093]!H
    p=p+3; coordinates(p:p+2) = [ 4.24637,  1.27773, -0.86698]!H
    p=p+3; coordinates(p:p+2) = [ 3.83542,  0.83590, -2.20574]!H
    p=p+3; coordinates(p:p+2) = [-1.33710,  4.11792, -0.55783]!H
    p=p+3; coordinates(p:p+2) = [-0.15858,  4.05926, -1.54358]!H
    p=p+3; coordinates(p:p+2) = [ 3.70433,  3.38958,  2.75945]!H
    p=p+3; coordinates(p:p+2) = [ 2.87922,  3.94545,  1.55536]!H
    p=p+3; coordinates(p:p+2) = [ 0.48484,  1.78382, -0.82579]!H
    p=p+3; coordinates(p:p+2) = [ 1.83919,  1.86647, -1.67688]!H
    p=p+3; coordinates(p:p+2) = [ 1.36128,  7.31484, -0.25489]!H
    p=p+3; coordinates(p:p+2) = [ 0.50437,  5.95037, -0.40769]!H
    p=p+3; coordinates(p:p+2) = [ 2.15403,  0.94247,  2.28060]!H
    p=p+3; coordinates(p:p+2) = [ 2.89100, -0.27132,  1.62148]!H
    p=p+3; coordinates(p:p+2) = [ 1.30408,  2.61876,  2.42595]!H
    p=p+3; coordinates(p:p+2) = [-0.16326,  2.47088,  2.68855]!H
    p=p+3; coordinates(p:p+2) = [ 2.63809,  5.08587, -0.37526]!H
    p=p+3; coordinates(p:p+2) = [ 2.24243,  3.64913, -0.71816]!H
    p=p+3; coordinates(p:p+2) = [-2.24224,  3.03902,  1.19737]!H
    p=p+3; coordinates(p:p+2) = [-1.71070,  2.17173,  0.06858]!H
    p=p+3; coordinates(p:p+2) = [ 0.72546,  5.83903,  1.80321]!H
    p=p+3; coordinates(p:p+2) = [ 1.40607,  5.02224,  3.03626]!H
    p=p+3; coordinates(p:p+2) = [ 4.51732,  1.44236,  2.67606]!H
    p=p+3; coordinates(p:p+2) = [ 5.39808,  2.08574,  3.73668]!H
    p=p+3; coordinates(p:p+2) = [ 5.55434,  2.62557,  0.92644]!H
    p=p+3; coordinates(p:p+2) = [ 4.66419,  3.33488, -0.23372]!H
    p=p+3; coordinates(p:p+2) = [-1.45378,  3.50547,  3.98412]!H
    p=p+3; coordinates(p:p+2) = [-0.77806,  4.38400,  2.89675]!H
    
    p=p+3; coordinates(p:p+2) = [( 1.80482+7.0), -1.73505, -0.85955]!H x
    p=p+3; coordinates(p:p+2) = [( 2.42260+7.0), -0.35826, -0.70808]!H
    p=p+3; coordinates(p:p+2) = [(-0.04077+7.0), -0.15230, -0.04336]!H
    p=p+3; coordinates(p:p+2) = [(-0.01324+7.0),  0.97891,  0.98093]!H
    p=p+3; coordinates(p:p+2) = [( 4.24637+7.0),  1.27773, -0.86698]!H
    p=p+3; coordinates(p:p+2) = [( 3.83542+7.0),  0.83590, -2.20574]!H
    p=p+3; coordinates(p:p+2) = [(-1.33710+7.0),  4.11792, -0.55783]!H
    p=p+3; coordinates(p:p+2) = [(-0.15858+7.0),  4.05926, -1.54358]!H
    p=p+3; coordinates(p:p+2) = [( 3.70433+7.0),  3.38958,  2.75945]!H
    p=p+3; coordinates(p:p+2) = [( 2.87922+7.0),  3.94545,  1.55536]!H
    p=p+3; coordinates(p:p+2) = [( 0.48484+7.0),  1.78382, -0.82579]!H
    p=p+3; coordinates(p:p+2) = [( 1.83919+7.0),  1.86647, -1.67688]!H
    p=p+3; coordinates(p:p+2) = [( 1.36128+7.0),  7.31484, -0.25489]!H
    p=p+3; coordinates(p:p+2) = [( 0.50437+7.0),  5.95037, -0.40769]!H
    p=p+3; coordinates(p:p+2) = [( 2.15403+7.0),  0.94247,  2.28060]!H
    p=p+3; coordinates(p:p+2) = [( 2.89100+7.0), -0.27132,  1.62148]!H
    p=p+3; coordinates(p:p+2) = [( 1.30408+7.0),  2.61876,  2.42595]!H
    p=p+3; coordinates(p:p+2) = [(-0.16326+7.0),  2.47088,  2.68855]!H
    p=p+3; coordinates(p:p+2) = [( 2.63809+7.0),  5.08587, -0.37526]!H
    p=p+3; coordinates(p:p+2) = [( 2.24243+7.0),  3.64913, -0.71816]!H
    p=p+3; coordinates(p:p+2) = [(-2.24224+7.0),  3.03902,  1.19737]!H
    p=p+3; coordinates(p:p+2) = [(-1.71070+7.0),  2.17173,  0.06858]!H
    p=p+3; coordinates(p:p+2) = [( 0.72546+7.0),  5.83903,  1.80321]!H
    p=p+3; coordinates(p:p+2) = [( 1.40607+7.0),  5.02224,  3.03626]!H
    p=p+3; coordinates(p:p+2) = [( 4.51732+7.0),  1.44236,  2.67606]!H
    p=p+3; coordinates(p:p+2) = [( 5.39808+7.0),  2.08574,  3.73668]!H
    p=p+3; coordinates(p:p+2) = [( 5.55434+7.0),  2.62557,  0.92644]!H
    p=p+3; coordinates(p:p+2) = [( 4.66419+7.0),  3.33488, -0.23372]!H
    p=p+3; coordinates(p:p+2) = [(-1.45378+7.0),  3.50547,  3.98412]!H
    p=p+3; coordinates(p:p+2) = [(-0.77806+7.0),  4.38400,  2.89675]!H

    p=p+3; coordinates(p:p+2) = [ 1.80482, (-1.73505+7.0), -0.85955]!H y
    p=p+3; coordinates(p:p+2) = [ 2.42260, (-0.35826+7.0), -0.70808]!H
    p=p+3; coordinates(p:p+2) = [-0.04077, (-0.15230+7.0), -0.04336]!H
    p=p+3; coordinates(p:p+2) = [-0.01324, ( 0.97891+7.0),  0.98093]!H
    p=p+3; coordinates(p:p+2) = [ 4.24637, ( 1.27773+7.0), -0.86698]!H
    p=p+3; coordinates(p:p+2) = [ 3.83542, ( 0.83590+7.0), -2.20574]!H
    p=p+3; coordinates(p:p+2) = [-1.33710, ( 4.11792+7.0), -0.55783]!H
    p=p+3; coordinates(p:p+2) = [-0.15858, ( 4.05926+7.0), -1.54358]!H
    p=p+3; coordinates(p:p+2) = [ 3.70433, ( 3.38958+7.0),  2.75945]!H
    p=p+3; coordinates(p:p+2) = [ 2.87922, ( 3.94545+7.0),  1.55536]!H
    p=p+3; coordinates(p:p+2) = [ 0.48484, ( 1.78382+7.0), -0.82579]!H
    p=p+3; coordinates(p:p+2) = [ 1.83919, ( 1.86647+7.0), -1.67688]!H
    p=p+3; coordinates(p:p+2) = [ 1.36128, ( 7.31484+7.0), -0.25489]!H
    p=p+3; coordinates(p:p+2) = [ 0.50437, ( 5.95037+7.0), -0.40769]!H
    p=p+3; coordinates(p:p+2) = [ 2.15403, ( 0.94247+7.0),  2.28060]!H
    p=p+3; coordinates(p:p+2) = [ 2.89100, (-0.27132+7.0),  1.62148]!H
    p=p+3; coordinates(p:p+2) = [ 1.30408, ( 2.61876+7.0),  2.42595]!H
    p=p+3; coordinates(p:p+2) = [-0.16326, ( 2.47088+7.0),  2.68855]!H
    p=p+3; coordinates(p:p+2) = [ 2.63809, ( 5.08587+7.0), -0.37526]!H
    p=p+3; coordinates(p:p+2) = [ 2.24243, ( 3.64913+7.0), -0.71816]!H
    p=p+3; coordinates(p:p+2) = [-2.24224, ( 3.03902+7.0),  1.19737]!H
    p=p+3; coordinates(p:p+2) = [-1.71070, ( 2.17173+7.0),  0.06858]!H
    p=p+3; coordinates(p:p+2) = [ 0.72546, ( 5.83903+7.0),  1.80321]!H
    p=p+3; coordinates(p:p+2) = [ 1.40607, ( 5.02224+7.0),  3.03626]!H
    p=p+3; coordinates(p:p+2) = [ 4.51732, ( 1.44236+7.0),  2.67606]!H
    p=p+3; coordinates(p:p+2) = [ 5.39808, ( 2.08574+7.0),  3.73668]!H
    p=p+3; coordinates(p:p+2) = [ 5.55434, ( 2.62557+7.0),  0.92644]!H
    p=p+3; coordinates(p:p+2) = [ 4.66419, ( 3.33488+7.0), -0.23372]!H
    p=p+3; coordinates(p:p+2) = [-1.45378, ( 3.50547+7.0),  3.98412]!H
    p=p+3; coordinates(p:p+2) = [-0.77806, ( 4.38400+7.0),  2.89675]!H

    p=p+3; coordinates(p:p+2) = [ 1.80482, -1.73505, (-0.85955+7.0)]!H z
    p=p+3; coordinates(p:p+2) = [ 2.42260, -0.35826, (-0.70808+7.0)]!H
    p=p+3; coordinates(p:p+2) = [-0.04077, -0.15230, (-0.04336+7.0)]!H
    p=p+3; coordinates(p:p+2) = [-0.01324,  0.97891, ( 0.98093+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 4.24637,  1.27773, (-0.86698+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 3.83542,  0.83590, (-2.20574+7.0)]!H
    p=p+3; coordinates(p:p+2) = [-1.33710,  4.11792, (-0.55783+7.0)]!H
    p=p+3; coordinates(p:p+2) = [-0.15858,  4.05926, (-1.54358+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 3.70433,  3.38958, ( 2.75945+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 2.87922,  3.94545, ( 1.55536+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 0.48484,  1.78382, (-0.82579+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 1.83919,  1.86647, (-1.67688+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 1.36128,  7.31484, (-0.25489+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 0.50437,  5.95037, (-0.40769+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 2.15403,  0.94247, ( 2.28060+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 2.89100, -0.27132, ( 1.62148+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 1.30408,  2.61876, ( 2.42595+7.0)]!H
    p=p+3; coordinates(p:p+2) = [-0.16326,  2.47088, ( 2.68855+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 2.63809,  5.08587, (-0.37526+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 2.24243,  3.64913, (-0.71816+7.0)]!H
    p=p+3; coordinates(p:p+2) = [-2.24224,  3.03902, ( 1.19737+7.0)]!H
    p=p+3; coordinates(p:p+2) = [-1.71070,  2.17173, ( 0.06858+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 0.72546,  5.83903, ( 1.80321+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 1.40607,  5.02224, ( 3.03626+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 4.51732,  1.44236, ( 2.67606+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 5.39808,  2.08574, ( 3.73668+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 5.55434,  2.62557, ( 0.92644+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 4.66419,  3.33488, (-0.23372+7.0)]!H
    p=p+3; coordinates(p:p+2) = [-1.45378,  3.50547, ( 3.98412+7.0)]!H
    p=p+3; coordinates(p:p+2) = [-0.77806,  4.38400, ( 2.89675+7.0)]!H

    p=p+3; coordinates(p:p+2) = [ 1.80482, (-1.73505+7.0), (-0.85955+7.0)]!H yz
    p=p+3; coordinates(p:p+2) = [ 2.42260, (-0.35826+7.0), (-0.70808+7.0)]!H
    p=p+3; coordinates(p:p+2) = [-0.04077, (-0.15230+7.0), (-0.04336+7.0)]!H
    p=p+3; coordinates(p:p+2) = [-0.01324, ( 0.97891+7.0), ( 0.98093+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 4.24637, ( 1.27773+7.0), (-0.86698+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 3.83542, ( 0.83590+7.0), (-2.20574+7.0)]!H
    p=p+3; coordinates(p:p+2) = [-1.33710, ( 4.11792+7.0), (-0.55783+7.0)]!H
    p=p+3; coordinates(p:p+2) = [-0.15858, ( 4.05926+7.0), (-1.54358+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 3.70433, ( 3.38958+7.0), ( 2.75945+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 2.87922, ( 3.94545+7.0), ( 1.55536+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 0.48484, ( 1.78382+7.0), (-0.82579+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 1.83919, ( 1.86647+7.0), (-1.67688+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 1.36128, ( 7.31484+7.0), (-0.25489+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 0.50437, ( 5.95037+7.0), (-0.40769+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 2.15403, ( 0.94247+7.0), ( 2.28060+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 2.89100, (-0.27132+7.0), ( 1.62148+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 1.30408, ( 2.61876+7.0), ( 2.42595+7.0)]!H
    p=p+3; coordinates(p:p+2) = [-0.16326, ( 2.47088+7.0), ( 2.68855+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 2.63809, ( 5.08587+7.0), (-0.37526+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 2.24243, ( 3.64913+7.0), (-0.71816+7.0)]!H
    p=p+3; coordinates(p:p+2) = [-2.24224, ( 3.03902+7.0), ( 1.19737+7.0)]!H
    p=p+3; coordinates(p:p+2) = [-1.71070, ( 2.17173+7.0), ( 0.06858+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 0.72546, ( 5.83903+7.0), ( 1.80321+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 1.40607, ( 5.02224+7.0), ( 3.03626+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 4.51732, ( 1.44236+7.0), ( 2.67606+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 5.39808, ( 2.08574+7.0), ( 3.73668+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 5.55434, ( 2.62557+7.0), ( 0.92644+7.0)]!H
    p=p+3; coordinates(p:p+2) = [ 4.66419, ( 3.33488+7.0), (-0.23372+7.0)]!H
    p=p+3; coordinates(p:p+2) = [-1.45378, ( 3.50547+7.0), ( 3.98412+7.0)]!H
    p=p+3; coordinates(p:p+2) = [-0.77806, ( 4.38400+7.0), ( 2.89675+7.0)]!H  
    
    p=p+3; coordinates(p:p+2) = [( 1.80482+7.0), -1.73505, (-0.85955+7.0)]!H xz
    p=p+3; coordinates(p:p+2) = [( 2.42260+7.0), -0.35826, (-0.70808+7.0)]!H
    p=p+3; coordinates(p:p+2) = [(-0.04077+7.0), -0.15230, (-0.04336+7.0)]!H
    p=p+3; coordinates(p:p+2) = [(-0.01324+7.0),  0.97891, ( 0.98093+7.0)]!H
    p=p+3; coordinates(p:p+2) = [( 4.24637+7.0),  1.27773, (-0.86698+7.0)]!H
    p=p+3; coordinates(p:p+2) = [( 3.83542+7.0),  0.83590, (-2.20574+7.0)]!H
    p=p+3; coordinates(p:p+2) = [(-1.33710+7.0),  4.11792, (-0.55783+7.0)]!H
    p=p+3; coordinates(p:p+2) = [(-0.15858+7.0),  4.05926, (-1.54358+7.0)]!H
    p=p+3; coordinates(p:p+2) = [( 3.70433+7.0),  3.38958, ( 2.75945+7.0)]!H
    p=p+3; coordinates(p:p+2) = [( 2.87922+7.0),  3.94545, ( 1.55536+7.0)]!H
    p=p+3; coordinates(p:p+2) = [( 0.48484+7.0),  1.78382, (-0.82579+7.0)]!H
    p=p+3; coordinates(p:p+2) = [( 1.83919+7.0),  1.86647, (-1.67688+7.0)]!H
    p=p+3; coordinates(p:p+2) = [( 1.36128+7.0),  7.31484, (-0.25489+7.0)]!H
    p=p+3; coordinates(p:p+2) = [( 0.50437+7.0),  5.95037, (-0.40769+7.0)]!H
    p=p+3; coordinates(p:p+2) = [( 2.15403+7.0),  0.94247, ( 2.28060+7.0)]!H
    p=p+3; coordinates(p:p+2) = [( 2.89100+7.0), -0.27132, ( 1.62148+7.0)]!H
    p=p+3; coordinates(p:p+2) = [( 1.30408+7.0),  2.61876, ( 2.42595+7.0)]!H
    p=p+3; coordinates(p:p+2) = [(-0.16326+7.0),  2.47088, ( 2.68855+7.0)]!H
    p=p+3; coordinates(p:p+2) = [( 2.63809+7.0),  5.08587, (-0.37526+7.0)]!H
    p=p+3; coordinates(p:p+2) = [( 2.24243+7.0),  3.64913, (-0.71816+7.0)]!H
    p=p+3; coordinates(p:p+2) = [(-2.24224+7.0),  3.03902, ( 1.19737+7.0)]!H
    p=p+3; coordinates(p:p+2) = [(-1.71070+7.0),  2.17173, ( 0.06858+7.0)]!H
    p=p+3; coordinates(p:p+2) = [( 0.72546+7.0),  5.83903, ( 1.80321+7.0)]!H
    p=p+3; coordinates(p:p+2) = [( 1.40607+7.0),  5.02224, ( 3.03626+7.0)]!H
    p=p+3; coordinates(p:p+2) = [( 4.51732+7.0),  1.44236, ( 2.67606+7.0)]!H
    p=p+3; coordinates(p:p+2) = [( 5.39808+7.0),  2.08574, ( 3.73668+7.0)]!H
    p=p+3; coordinates(p:p+2) = [( 5.55434+7.0),  2.62557, ( 0.92644+7.0)]!H
    p=p+3; coordinates(p:p+2) = [( 4.66419+7.0),  3.33488, (-0.23372+7.0)]!H
    p=p+3; coordinates(p:p+2) = [(-1.45378+7.0),  3.50547, ( 3.98412+7.0)]!H
    p=p+3; coordinates(p:p+2) = [(-0.77806+7.0),  4.38400, ( 2.89675+7.0)]!H      

    p=p+3; coordinates(p:p+2) = [( 1.80482+7.0), (-1.73505+7.0), -0.85955]!H xy
    p=p+3; coordinates(p:p+2) = [( 2.42260+7.0), (-0.35826+7.0), -0.70808]!H
    p=p+3; coordinates(p:p+2) = [(-0.04077+7.0), (-0.15230+7.0), -0.04336]!H
    p=p+3; coordinates(p:p+2) = [(-0.01324+7.0), ( 0.97891+7.0),  0.98093]!H
    p=p+3; coordinates(p:p+2) = [( 4.24637+7.0), ( 1.27773+7.0), -0.86698]!H
    p=p+3; coordinates(p:p+2) = [( 3.83542+7.0), ( 0.83590+7.0), -2.20574]!H
    p=p+3; coordinates(p:p+2) = [(-1.33710+7.0), ( 4.11792+7.0), -0.55783]!H
    p=p+3; coordinates(p:p+2) = [(-0.15858+7.0), ( 4.05926+7.0), -1.54358]!H
    p=p+3; coordinates(p:p+2) = [( 3.70433+7.0), ( 3.38958+7.0),  2.75945]!H
    p=p+3; coordinates(p:p+2) = [( 2.87922+7.0), ( 3.94545+7.0),  1.55536]!H
    p=p+3; coordinates(p:p+2) = [( 0.48484+7.0), ( 1.78382+7.0), -0.82579]!H
    p=p+3; coordinates(p:p+2) = [( 1.83919+7.0), ( 1.86647+7.0), -1.67688]!H
    p=p+3; coordinates(p:p+2) = [( 1.36128+7.0), ( 7.31484+7.0), -0.25489]!H
    p=p+3; coordinates(p:p+2) = [( 0.50437+7.0), ( 5.95037+7.0), -0.40769]!H
    p=p+3; coordinates(p:p+2) = [( 2.15403+7.0), ( 0.94247+7.0),  2.28060]!H
    p=p+3; coordinates(p:p+2) = [( 2.89100+7.0), (-0.27132+7.0),  1.62148]!H
    p=p+3; coordinates(p:p+2) = [( 1.30408+7.0), ( 2.61876+7.0),  2.42595]!H
    p=p+3; coordinates(p:p+2) = [(-0.16326+7.0), ( 2.47088+7.0),  2.68855]!H
    p=p+3; coordinates(p:p+2) = [( 2.63809+7.0), ( 5.08587+7.0), -0.37526]!H
    p=p+3; coordinates(p:p+2) = [( 2.24243+7.0), ( 3.64913+7.0), -0.71816]!H
    p=p+3; coordinates(p:p+2) = [(-2.24224+7.0), ( 3.03902+7.0),  1.19737]!H
    p=p+3; coordinates(p:p+2) = [(-1.71070+7.0), ( 2.17173+7.0),  0.06858]!H
    p=p+3; coordinates(p:p+2) = [( 0.72546+7.0), ( 5.83903+7.0),  1.80321]!H
    p=p+3; coordinates(p:p+2) = [( 1.40607+7.0), ( 5.02224+7.0),  3.03626]!H
    p=p+3; coordinates(p:p+2) = [( 4.51732+7.0), ( 1.44236+7.0),  2.67606]!H
    p=p+3; coordinates(p:p+2) = [( 5.39808+7.0), ( 2.08574+7.0),  3.73668]!H
    p=p+3; coordinates(p:p+2) = [( 5.55434+7.0), ( 2.62557+7.0),  0.92644]!H
    p=p+3; coordinates(p:p+2) = [( 4.66419+7.0), ( 3.33488+7.0), -0.23372]!H
    p=p+3; coordinates(p:p+2) = [(-1.45378+7.0), ( 3.50547+7.0),  3.98412]!H
    p=p+3; coordinates(p:p+2) = [(-0.77806+7.0), ( 4.38400+7.0),  2.89675]!H

    p=p+3; coordinates(p:p+2) = [( 1.80482+7.0), (-1.73505+7.0), (-0.85955+7.0)]!H xyz
    p=p+3; coordinates(p:p+2) = [( 2.42260+7.0), (-0.35826+7.0), (-0.70808+7.0)]!H
    p=p+3; coordinates(p:p+2) = [(-0.04077+7.0), (-0.15230+7.0), (-0.04336+7.0)]!H
    p=p+3; coordinates(p:p+2) = [(-0.01324+7.0), ( 0.97891+7.0), ( 0.98093+7.0)]!H
    p=p+3; coordinates(p:p+2) = [( 4.24637+7.0), ( 1.27773+7.0), (-0.86698+7.0)]!H
    p=p+3; coordinates(p:p+2) = [( 3.83542+7.0), ( 0.83590+7.0), (-2.20574+7.0)]!H
    p=p+3; coordinates(p:p+2) = [(-1.33710+7.0), ( 4.11792+7.0), (-0.55783+7.0)]!H
    p=p+3; coordinates(p:p+2) = [(-0.15858+7.0), ( 4.05926+7.0), (-1.54358+7.0)]!H
    p=p+3; coordinates(p:p+2) = [( 3.70433+7.0), ( 3.38958+7.0), ( 2.75945+7.0)]!H
    p=p+3; coordinates(p:p+2) = [( 2.87922+7.0), ( 3.94545+7.0), ( 1.55536+7.0)]!H
    p=p+3; coordinates(p:p+2) = [( 0.48484+7.0), ( 1.78382+7.0), (-0.82579+7.0)]!H
    p=p+3; coordinates(p:p+2) = [( 1.83919+7.0), ( 1.86647+7.0), (-1.67688+7.0)]!H
    p=p+3; coordinates(p:p+2) = [( 1.36128+7.0), ( 7.31484+7.0), (-0.25489+7.0)]!H
    p=p+3; coordinates(p:p+2) = [( 0.50437+7.0), ( 5.95037+7.0), (-0.40769+7.0)]!H
    p=p+3; coordinates(p:p+2) = [( 2.15403+7.0), ( 0.94247+7.0), ( 2.28060+7.0)]!H
    p=p+3; coordinates(p:p+2) = [( 2.89100+7.0), (-0.27132+7.0), ( 1.62148+7.0)]!H
    p=p+3; coordinates(p:p+2) = [( 1.30408+7.0), ( 2.61876+7.0), ( 2.42595+7.0)]!H
    p=p+3; coordinates(p:p+2) = [(-0.16326+7.0), ( 2.47088+7.0), ( 2.68855+7.0)]!H
    p=p+3; coordinates(p:p+2) = [( 2.63809+7.0), ( 5.08587+7.0), (-0.37526+7.0)]!H
    p=p+3; coordinates(p:p+2) = [( 2.24243+7.0), ( 3.64913+7.0), (-0.71816+7.0)]!H
    p=p+3; coordinates(p:p+2) = [(-2.24224+7.0), ( 3.03902+7.0), ( 1.19737+7.0)]!H
    p=p+3; coordinates(p:p+2) = [(-1.71070+7.0), ( 2.17173+7.0), ( 0.06858+7.0)]!H
    p=p+3; coordinates(p:p+2) = [( 0.72546+7.0), ( 5.83903+7.0), ( 1.80321+7.0)]!H
    p=p+3; coordinates(p:p+2) = [( 1.40607+7.0), ( 5.02224+7.0), ( 3.03626+7.0)]!H
    p=p+3; coordinates(p:p+2) = [( 4.51732+7.0), ( 1.44236+7.0), ( 2.67606+7.0)]!H
    p=p+3; coordinates(p:p+2) = [( 5.39808+7.0), ( 2.08574+7.0), ( 3.73668+7.0)]!H
    p=p+3; coordinates(p:p+2) = [( 5.55434+7.0), ( 2.62557+7.0), ( 0.92644+7.0)]!H
    p=p+3; coordinates(p:p+2) = [( 4.66419+7.0), ( 3.33488+7.0), (-0.23372+7.0)]!H
    p=p+3; coordinates(p:p+2) = [(-1.45378+7.0), ( 3.50547+7.0), ( 3.98412+7.0)]!H
    p=p+3; coordinates(p:p+2) = [(-0.77806+7.0), ( 4.38400+7.0), ( 2.89675+7.0)]!H

    
    p=p+3; coordinates(p:p+2) = [ 1.66322, -0.92324, -0.35186]!O 0
    p=p+3; coordinates(p:p+2) = [-0.45776,  0.78681,  0.11944]!O
    p=p+3; coordinates(p:p+2) = [ 3.53951,  0.77639, -1.32046]!O
    p=p+3; coordinates(p:p+2) = [-0.57027,  4.66525, -0.91234]!O
    p=p+3; coordinates(p:p+2) = [ 2.82467,  3.80518,  2.50783]!O
    p=p+3; coordinates(p:p+2) = [ 1.05912,  2.41662, -1.33542]!O
    p=p+3; coordinates(p:p+2) = [ 1.24000,  6.41845,  0.03330]!O
    p=p+3; coordinates(p:p+2) = [ 2.98116,  0.48514,  2.18258]!O
    p=p+3; coordinates(p:p+2) = [ 0.60944,  1.93689,  2.41997]!O
    p=p+3; coordinates(p:p+2) = [ 2.96096,  4.17215, -0.24367]!O
    p=p+3; coordinates(p:p+2) = [-2.15967,  3.04991,  0.22393]!O
    p=p+3; coordinates(p:p+2) = [ 0.56317,  5.47605,  2.72575]!O
    p=p+3; coordinates(p:p+2) = [ 5.08862,  2.23134,  2.85097]!O
    p=p+3; coordinates(p:p+2) = [ 5.20638,  2.53755, -0.01077]!O
    p=p+3; coordinates(p:p+2) = [-1.42229,  3.63362,  3.02443]!O

    p=p+3; coordinates(p:p+2) = [( 1.66322+7.0), -0.92324, -0.35186]!O x
    p=p+3; coordinates(p:p+2) = [(-0.45776+7.0),  0.78681,  0.11944]!O
    p=p+3; coordinates(p:p+2) = [( 3.53951+7.0),  0.77639, -1.32046]!O
    p=p+3; coordinates(p:p+2) = [(-0.57027+7.0),  4.66525, -0.91234]!O
    p=p+3; coordinates(p:p+2) = [( 2.82467+7.0),  3.80518,  2.50783]!O
    p=p+3; coordinates(p:p+2) = [( 1.05912+7.0),  2.41662, -1.33542]!O
    p=p+3; coordinates(p:p+2) = [( 1.24000+7.0),  6.41845,  0.03330]!O
    p=p+3; coordinates(p:p+2) = [( 2.98116+7.0),  0.48514,  2.18258]!O
    p=p+3; coordinates(p:p+2) = [( 0.60944+7.0),  1.93689,  2.41997]!O
    p=p+3; coordinates(p:p+2) = [( 2.96096+7.0),  4.17215, -0.24367]!O
    p=p+3; coordinates(p:p+2) = [(-2.15967+7.0),  3.04991,  0.22393]!O
    p=p+3; coordinates(p:p+2) = [( 0.56317+7.0),  5.47605,  2.72575]!O
    p=p+3; coordinates(p:p+2) = [( 5.08862+7.0),  2.23134,  2.85097]!O
    p=p+3; coordinates(p:p+2) = [( 5.20638+7.0),  2.53755, -0.01077]!O
    p=p+3; coordinates(p:p+2) = [(-1.42229+7.0),  3.63362,  3.02443]!O


    p=p+3; coordinates(p:p+2) = [ 1.66322, (-0.92324+7.0), -0.35186]!O y
    p=p+3; coordinates(p:p+2) = [-0.45776, ( 0.78681+7.0),  0.11944]!O
    p=p+3; coordinates(p:p+2) = [ 3.53951, ( 0.77639+7.0), -1.32046]!O
    p=p+3; coordinates(p:p+2) = [-0.57027, ( 4.66525+7.0), -0.91234]!O
    p=p+3; coordinates(p:p+2) = [ 2.82467, ( 3.80518+7.0),  2.50783]!O
    p=p+3; coordinates(p:p+2) = [ 1.05912, ( 2.41662+7.0), -1.33542]!O
    p=p+3; coordinates(p:p+2) = [ 1.24000, ( 6.41845+7.0),  0.03330]!O
    p=p+3; coordinates(p:p+2) = [ 2.98116, ( 0.48514+7.0),  2.18258]!O
    p=p+3; coordinates(p:p+2) = [ 0.60944, ( 1.93689+7.0),  2.41997]!O
    p=p+3; coordinates(p:p+2) = [ 2.96096, ( 4.17215+7.0), -0.24367]!O
    p=p+3; coordinates(p:p+2) = [-2.15967, ( 3.04991+7.0),  0.22393]!O
    p=p+3; coordinates(p:p+2) = [ 0.56317, ( 5.47605+7.0),  2.72575]!O
    p=p+3; coordinates(p:p+2) = [ 5.08862, ( 2.23134+7.0),  2.85097]!O
    p=p+3; coordinates(p:p+2) = [ 5.20638, ( 2.53755+7.0), -0.01077]!O
    p=p+3; coordinates(p:p+2) = [-1.42229, ( 3.63362+7.0),  3.02443]!O

    p=p+3; coordinates(p:p+2) = [ 1.66322, -0.92324, (-0.35186+7.0)]!O z
    p=p+3; coordinates(p:p+2) = [-0.45776,  0.78681, ( 0.11944+7.0)]!O
    p=p+3; coordinates(p:p+2) = [ 3.53951,  0.77639, (-1.32046+7.0)]!O
    p=p+3; coordinates(p:p+2) = [-0.57027,  4.66525, (-0.91234+7.0)]!O
    p=p+3; coordinates(p:p+2) = [ 2.82467,  3.80518, ( 2.50783+7.0)]!O
    p=p+3; coordinates(p:p+2) = [ 1.05912,  2.41662, (-1.33542+7.0)]!O
    p=p+3; coordinates(p:p+2) = [ 1.24000,  6.41845, ( 0.03330+7.0)]!O
    p=p+3; coordinates(p:p+2) = [ 2.98116,  0.48514, ( 2.18258+7.0)]!O
    p=p+3; coordinates(p:p+2) = [ 0.60944,  1.93689, ( 2.41997+7.0)]!O
    p=p+3; coordinates(p:p+2) = [ 2.96096,  4.17215, (-0.24367+7.0)]!O
    p=p+3; coordinates(p:p+2) = [-2.15967,  3.04991, ( 0.22393+7.0)]!O
    p=p+3; coordinates(p:p+2) = [ 0.56317,  5.47605, ( 2.72575+7.0)]!O
    p=p+3; coordinates(p:p+2) = [ 5.08862,  2.23134, ( 2.85097+7.0)]!O
    p=p+3; coordinates(p:p+2) = [ 5.20638,  2.53755, (-0.01077+7.0)]!O
    p=p+3; coordinates(p:p+2) = [-1.42229,  3.63362, ( 3.02443+7.0)]!O

    p=p+3; coordinates(p:p+2) = [ 1.66322, (-0.92324+7.0), (-0.35186+7.0)]!O yz
    p=p+3; coordinates(p:p+2) = [-0.45776, ( 0.78681+7.0), ( 0.11944+7.0)]!O
    p=p+3; coordinates(p:p+2) = [ 3.53951, ( 0.77639+7.0), (-1.32046+7.0)]!O
    p=p+3; coordinates(p:p+2) = [-0.57027, ( 4.66525+7.0), (-0.91234+7.0)]!O
    p=p+3; coordinates(p:p+2) = [ 2.82467, ( 3.80518+7.0), ( 2.50783+7.0)]!O
    p=p+3; coordinates(p:p+2) = [ 1.05912, ( 2.41662+7.0), (-1.33542+7.0)]!O
    p=p+3; coordinates(p:p+2) = [ 1.24000, ( 6.41845+7.0), ( 0.03330+7.0)]!O
    p=p+3; coordinates(p:p+2) = [ 2.98116, ( 0.48514+7.0), ( 2.18258+7.0)]!O
    p=p+3; coordinates(p:p+2) = [ 0.60944, ( 1.93689+7.0), ( 2.41997+7.0)]!O
    p=p+3; coordinates(p:p+2) = [ 2.96096, ( 4.17215+7.0), (-0.24367+7.0)]!O
    p=p+3; coordinates(p:p+2) = [-2.15967, ( 3.04991+7.0), ( 0.22393+7.0)]!O
    p=p+3; coordinates(p:p+2) = [ 0.56317, ( 5.47605+7.0), ( 2.72575+7.0)]!O
    p=p+3; coordinates(p:p+2) = [ 5.08862, ( 2.23134+7.0), ( 2.85097+7.0)]!O
    p=p+3; coordinates(p:p+2) = [ 5.20638, ( 2.53755+7.0), (-0.01077+7.0)]!O
    p=p+3; coordinates(p:p+2) = [-1.42229, ( 3.63362+7.0), ( 3.02443+7.0)]!O

    p=p+3; coordinates(p:p+2) = [( 1.66322+7.0), -0.92324, (-0.35186+7.0)]!O xz
    p=p+3; coordinates(p:p+2) = [(-0.45776+7.0),  0.78681, ( 0.11944+7.0)]!O
    p=p+3; coordinates(p:p+2) = [( 3.53951+7.0),  0.77639, (-1.32046+7.0)]!O
    p=p+3; coordinates(p:p+2) = [(-0.57027+7.0),  4.66525, (-0.91234+7.0)]!O
    p=p+3; coordinates(p:p+2) = [( 2.82467+7.0),  3.80518, ( 2.50783+7.0)]!O
    p=p+3; coordinates(p:p+2) = [( 1.05912+7.0),  2.41662, (-1.33542+7.0)]!O
    p=p+3; coordinates(p:p+2) = [( 1.24000+7.0),  6.41845, ( 0.03330+7.0)]!O
    p=p+3; coordinates(p:p+2) = [( 2.98116+7.0),  0.48514, ( 2.18258+7.0)]!O
    p=p+3; coordinates(p:p+2) = [( 0.60944+7.0),  1.93689, ( 2.41997+7.0)]!O
    p=p+3; coordinates(p:p+2) = [( 2.96096+7.0),  4.17215, (-0.24367+7.0)]!O
    p=p+3; coordinates(p:p+2) = [(-2.15967+7.0),  3.04991, ( 0.22393+7.0)]!O
    p=p+3; coordinates(p:p+2) = [( 0.56317+7.0),  5.47605, ( 2.72575+7.0)]!O
    p=p+3; coordinates(p:p+2) = [( 5.08862+7.0),  2.23134, ( 2.85097+7.0)]!O
    p=p+3; coordinates(p:p+2) = [( 5.20638+7.0),  2.53755, (-0.01077+7.0)]!O
    p=p+3; coordinates(p:p+2) = [(-1.42229+7.0),  3.63362, ( 3.02443+7.0)]!O
    
    p=p+3; coordinates(p:p+2) = [( 1.66322+7.0), (-0.92324+7.0), -0.35186]!O xy
    p=p+3; coordinates(p:p+2) = [(-0.45776+7.0), ( 0.78681+7.0),  0.11944]!O
    p=p+3; coordinates(p:p+2) = [( 3.53951+7.0), ( 0.77639+7.0), -1.32046]!O
    p=p+3; coordinates(p:p+2) = [(-0.57027+7.0), ( 4.66525+7.0), -0.91234]!O
    p=p+3; coordinates(p:p+2) = [( 2.82467+7.0), ( 3.80518+7.0),  2.50783]!O
    p=p+3; coordinates(p:p+2) = [( 1.05912+7.0), ( 2.41662+7.0), -1.33542]!O
    p=p+3; coordinates(p:p+2) = [( 1.24000+7.0), ( 6.41845+7.0),  0.03330]!O
    p=p+3; coordinates(p:p+2) = [( 2.98116+7.0), ( 0.48514+7.0),  2.18258]!O
    p=p+3; coordinates(p:p+2) = [( 0.60944+7.0), ( 1.93689+7.0),  2.41997]!O
    p=p+3; coordinates(p:p+2) = [( 2.96096+7.0), ( 4.17215+7.0), -0.24367]!O
    p=p+3; coordinates(p:p+2) = [(-2.15967+7.0), ( 3.04991+7.0),  0.22393]!O
    p=p+3; coordinates(p:p+2) = [( 0.56317+7.0), ( 5.47605+7.0),  2.72575]!O
    p=p+3; coordinates(p:p+2) = [( 5.08862+7.0), ( 2.23134+7.0),  2.85097]!O
    p=p+3; coordinates(p:p+2) = [( 5.20638+7.0), ( 2.53755+7.0), -0.01077]!O
    p=p+3; coordinates(p:p+2) = [(-1.42229+7.0), ( 3.63362+7.0),  3.02443]!O

    p=p+3; coordinates(p:p+2) = [( 1.66322+7.0), (-0.92324+7.0), (-0.35186+7.0)]!O xyz
    p=p+3; coordinates(p:p+2) = [(-0.45776+7.0), ( 0.78681+7.0), ( 0.11944+7.0)]!O
    p=p+3; coordinates(p:p+2) = [( 3.53951+7.0), ( 0.77639+7.0), (-1.32046+7.0)]!O
    p=p+3; coordinates(p:p+2) = [(-0.57027+7.0), ( 4.66525+7.0), (-0.91234+7.0)]!O
    p=p+3; coordinates(p:p+2) = [( 2.82467+7.0), ( 3.80518+7.0), ( 2.50783+7.0)]!O
    p=p+3; coordinates(p:p+2) = [( 1.05912+7.0), ( 2.41662+7.0), (-1.33542+7.0)]!O
    p=p+3; coordinates(p:p+2) = [( 1.24000+7.0), ( 6.41845+7.0), ( 0.03330+7.0)]!O
    p=p+3; coordinates(p:p+2) = [( 2.98116+7.0), ( 0.48514+7.0), ( 2.18258+7.0)]!O
    p=p+3; coordinates(p:p+2) = [( 0.60944+7.0), ( 1.93689+7.0), ( 2.41997+7.0)]!O
    p=p+3; coordinates(p:p+2) = [( 2.96096+7.0), ( 4.17215+7.0), (-0.24367+7.0)]!O
    p=p+3; coordinates(p:p+2) = [(-2.15967+7.0), ( 3.04991+7.0), ( 0.22393+7.0)]!O
    p=p+3; coordinates(p:p+2) = [( 0.56317+7.0), ( 5.47605+7.0), ( 2.72575+7.0)]!O
    p=p+3; coordinates(p:p+2) = [( 5.08862+7.0), ( 2.23134+7.0), ( 2.85097+7.0)]!O
    p=p+3; coordinates(p:p+2) = [( 5.20638+7.0), ( 2.53755+7.0), (-0.01077+7.0)]!O
    p=p+3; coordinates(p:p+2) = [(-1.42229+7.0), ( 3.63362+7.0), ( 3.02443+7.0)]!O


    print*, "things in coordinates should be",n_atoms*3, "and it is: ", p+2
    
    cell(1) = 35.000
    cell(2) = 35.000
    cell(3) = 35.000

    forces(:) = 0.0

    ! Call the scme function.
    do i=1,1
       call scme_calculate(n_atoms, coordinates, cell, forces, u_tot)
    end do

  end subroutine 
  
end module
