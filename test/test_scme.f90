! Copyright (C)  2015-2016  SCMEdev
! Licenced under LGPLv3. See LICENCE for details.


#include "mifu.h"

module test_scme
  use mifu_asserts

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
    do i=1,2000
       call scme_calculate(n_atoms, coordinates, cell, forces, u_tot)
    end do

  end subroutine test_scme_perf

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
           coordinates(p:p+2) = [ -1.0907, -0.8252, -0.9284 ]
    p=p+3; coordinates(p:p+2) = [ -1.4736, -1.3180, -1.6898 ]
    p=p+3; coordinates(p:p+2) = [ -0.8686,  0.0696, -1.2669 ]
    p=p+3; coordinates(p:p+2) = [ -0.1530,  2.5613,  0.8758 ]
    p=p+3; coordinates(p:p+2) = [ -0.7162,  3.3599,  0.8471 ]
    p=p+3; coordinates(p:p+2) = [ -0.8391,  1.8572,  1.0804 ]
    p=p+3; coordinates(p:p+2) = [ -0.6756,  1.5259, -2.2016 ]
    p=p+3; coordinates(p:p+2) = [ -0.3571,  2.1947, -1.5612 ]
    p=p+3; coordinates(p:p+2) = [ -1.4623,  2.0190, -2.5554 ]
    p=p+3; coordinates(p:p+2) = [ -2.1394,  0.5531,  1.1850 ]
    p=p+3; coordinates(p:p+2) = [ -1.8516, -0.1730,  0.5943 ]
    p=p+3; coordinates(p:p+2) = [ -2.6503,  0.1754,  1.8833 ]
    p=p+3; coordinates(p:p+2) = [ -1.9288,  4.7626,  0.8354 ]
    p=p+3; coordinates(p:p+2) = [ -1.6640,  4.9992, -0.1046 ]
    p=p+3; coordinates(p:p+2) = [ -2.1483,  5.5160,  1.3337 ]
    p=p+3; coordinates(p:p+2) = [ -3.6595,  2.6702, -0.1104 ]
    p=p+3; coordinates(p:p+2) = [ -3.1531,  1.9088,  0.2640 ]
    p=p+3; coordinates(p:p+2) = [ -3.2814,  3.4999,  0.2831 ]
    p=p+3; coordinates(p:p+2) = [  1.1160,  4.1296, -1.3211 ]
    p=p+3; coordinates(p:p+2) = [  1.5809,  4.9571, -1.0815 ]
    p=p+3; coordinates(p:p+2) = [  0.9649,  3.6761, -0.4809 ]
    p=p+3; coordinates(p:p+2) = [ -2.9721,  3.2005, -2.6291 ]
    p=p+3; coordinates(p:p+2) = [ -3.3448,  2.9950, -1.7281 ]
    p=p+3; coordinates(p:p+2) = [ -3.7567,  3.1501, -3.2142 ]
    p=p+3; coordinates(p:p+2) = [ -1.2805,  5.2116, -1.8113 ]
    p=p+3; coordinates(p:p+2) = [ -0.4045,  4.8032, -1.9164 ]
    p=p+3; coordinates(p:p+2) = [ -1.8429,  4.6655, -2.4245 ]
    print*, "things in coordinates should be 81 and it is: ", p
    cell(1) = 35.000
    cell(2) = 35.000
    cell(3) = 35.000

    forces(:) = 0.0

    ! Call the scme function.
    do i=1,200
       call scme_calculate(n_atoms, coordinates, cell, forces, u_tot)
    end do

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
           coordinates(p:p+2) = [ 0.6770, -1.0002,  1.3579 ]
    p=p+3; coordinates(p:p+2) = [ 0.3032, -1.8399,  1.0010 ]
    p=p+3; coordinates(p:p+2) = [ 0.9844, -0.4706,  0.6187 ]
    p=p+3; coordinates(p:p+2) = [-0.8010,  1.3294,  1.8051 ]
    p=p+3; coordinates(p:p+2) = [-0.8076,  0.3271,  1.6706 ]
    p=p+3; coordinates(p:p+2) = [-0.1854,  1.4969,  2.5873 ]
    p=p+3; coordinates(p:p+2) = [ 2.3540,  0.1624, -0.5255 ]
    p=p+3; coordinates(p:p+2) = [ 2.9254,  0.2162,  0.3102 ]
    p=p+3; coordinates(p:p+2) = [ 2.7652, -0.5340, -1.0151 ]
    p=p+3; coordinates(p:p+2) = [-0.2063,  4.8690,  1.2810 ]
    p=p+3; coordinates(p:p+2) = [-1.1037,  4.6825,  1.5252 ]
    p=p+3; coordinates(p:p+2) = [-0.0140,  4.3164,  0.5037 ]
    p=p+3; coordinates(p:p+2) = [ 2.8023,  2.7500,  4.2547 ]
    p=p+3; coordinates(p:p+2) = [ 3.1583,  1.8759,  4.3915 ]
    p=p+3; coordinates(p:p+2) = [ 2.8880,  2.9337,  3.3133 ]
    p=p+3; coordinates(p:p+2) = [ 0.7766,  2.5126, -0.2978 ]
    p=p+3; coordinates(p:p+2) = [ 0.0939,  2.0644,  0.2636 ]
    p=p+3; coordinates(p:p+2) = [ 1.3632,  1.8876, -0.7496 ]
    p=p+3; coordinates(p:p+2) = [ 2.3618,  5.5176,  2.1466 ]
    p=p+3; coordinates(p:p+2) = [ 2.4703,  6.4241,  2.3900 ]
    p=p+3; coordinates(p:p+2) = [ 1.4436,  5.4346,  1.7948 ]
    p=p+3; coordinates(p:p+2) = [ 1.4045, -1.0680,  4.0285 ]
    p=p+3; coordinates(p:p+2) = [ 0.9692, -0.1763,  4.1704 ]
    p=p+3; coordinates(p:p+2) = [ 1.3241, -1.2150,  3.0437 ]
    p=p+3; coordinates(p:p+2) = [ 0.2572,  1.5550,  4.5338 ]
    p=p+3; coordinates(p:p+2) = [ 0.9583,  2.2119,  4.5997 ]
    p=p+3; coordinates(p:p+2) = [-0.5835,  2.0534,  4.6742 ]
    p=p+3; coordinates(p:p+2) = [ 3.2170,  2.9856,  1.4426 ]
    p=p+3; coordinates(p:p+2) = [ 3.2209,  3.9213,  1.8040 ]
    p=p+3; coordinates(p:p+2) = [ 2.3920,  2.9388,  0.8785 ]
    p=p+3; coordinates(p:p+2) = [-2.5831,  3.8013,  1.8517 ]
    p=p+3; coordinates(p:p+2) = [-2.7007,  3.9437,  2.8047 ]
    p=p+3; coordinates(p:p+2) = [-2.2112,  2.8837,  1.7868 ]
    p=p+3; coordinates(p:p+2) = [ 0.3907,  4.5914,  4.2507 ]
    p=p+3; coordinates(p:p+2) = [ 0.3022,  4.6365,  3.2729 ]
    p=p+3; coordinates(p:p+2) = [ 1.2973,  4.2601,  4.3560 ]
    p=p+3; coordinates(p:p+2) = [ 3.9521, -0.0128,  4.5110 ]
    p=p+3; coordinates(p:p+2) = [ 3.1401, -0.5928,  4.2823 ]
    p=p+3; coordinates(p:p+2) = [ 4.3636, -0.3274,  5.3284 ]
    p=p+3; coordinates(p:p+2) = [ 4.2349,  0.4765,  1.5088 ]
    p=p+3; coordinates(p:p+2) = [ 4.2995,  0.0572,  2.3993 ]
    p=p+3; coordinates(p:p+2) = [ 4.0242,  1.4243,  1.6858 ]
    p=p+3; coordinates(p:p+2) = [-2.1422,  3.4141,  4.6134 ]
    p=p+3; coordinates(p:p+2) = [-2.5260,  3.5298,  5.4806 ]
    p=p+3; coordinates(p:p+2) = [-1.2304,  3.8854,  4.6447 ]
    print*, "things in coordinates should be",15*3*3, "and it is: ", p+2
    
    cell(1) = 35.000
    cell(2) = 35.000
    cell(3) = 35.000

    forces(:) = 0.0

    ! Call the scme function.
    do i=1,50
       call scme_calculate(n_atoms, coordinates, cell, forces, u_tot)
    end do

  end subroutine 

  subroutine test_scme_cluster_32_perf()
    implicit none
    ! ----------------------------------------
    integer, parameter :: n_atoms = 96 !32 water * 3
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
           coordinates(p:p+2) = [4.714724, 3.794822, 7.227688]
    p=p+3; coordinates(p:p+2) = [5.894676, 4.770402, 3.542073]
    p=p+3; coordinates(p:p+2) = [5.666188, 0.612496, 8.714184]
    p=p+3; coordinates(p:p+2) = [5.026691, 6.274405, 8.828260]
    p=p+3; coordinates(p:p+2) = [3.297651, 7.994921, 7.560544]
    p=p+3; coordinates(p:p+2) = [4.496765, 8.225602, 3.243301]
    p=p+3; coordinates(p:p+2) = [4.014890, 5.744731, 5.563480]
    p=p+3; coordinates(p:p+2) = [1.079311, 6.710587, 7.139180]
    p=p+3; coordinates(p:p+2) = [9.089848, 2.498606, 2.806837]
    p=p+3; coordinates(p:p+2) = [6.285261, 1.842312, 6.286177]
    p=p+3; coordinates(p:p+2) = [1.919755, 3.488047, 4.981780]
    p=p+3; coordinates(p:p+2) = [3.122451, 5.398157, 1.826458]
    p=p+3; coordinates(p:p+2) = [5.918245, 3.795972, 0.757720]
    p=p+3; coordinates(p:p+2) = [8.047289, 1.995078, 9.008911]
    p=p+3; coordinates(p:p+2) = [0.009837, 2.234008, 6.170529]
    p=p+3; coordinates(p:p+2) = [3.827995, 2.751079, 3.193511]
    p=p+3; coordinates(p:p+2) = [7.470885, 4.330890, 5.629196]
    p=p+3; coordinates(p:p+2) = [6.683945, 6.825986, 1.966956]
    p=p+3; coordinates(p:p+2) = [8.162331, 7.135795, 5.146089]
    p=p+3; coordinates(p:p+2) = [5.915112, 7.910117, 6.287469]
    p=p+3; coordinates(p:p+2) = [0.286200, 0.160794, 4.320814]
    p=p+3; coordinates(p:p+2) = [8.609000, 8.720642, 1.588505]
    p=p+3; coordinates(p:p+2) = [2.058608, 7.757951, 0.857307]
    p=p+3; coordinates(p:p+2) = [3.329270, 1.111310, 1.069920]
    p=p+3; coordinates(p:p+2) = [2.897618, 0.379327, 5.028460]
    p=p+3; coordinates(p:p+2) = [2.028507, 7.105495, 4.037806]
    p=p+3; coordinates(p:p+2) = [8.003827, 5.150162, 8.220378]
    p=p+3; coordinates(p:p+2) = [2.583109, 2.046129, 7.618021]
    p=p+3; coordinates(p:p+2) = [6.386365, 1.542067, 3.182596]
    p=p+3; coordinates(p:p+2) = [0.014544, 5.232013, 3.043986]
    p=p+3; coordinates(p:p+2) = [8.068713, 8.207080, 8.095860]
    p=p+3; coordinates(p:p+2) = [1.259448, 3.966680, 9.125192]
    p=p+3; coordinates(p:p+2) = [1.651000, 1.904133, 7.297683]
    p=p+3; coordinates(p:p+2) = [7.489300, 2.584646, 0.373518]
    p=p+3; coordinates(p:p+2) = [5.025503, 0.867717, 0.297158]
    p=p+3; coordinates(p:p+2) = [3.184744, 4.680107, 2.445076]
    p=p+3; coordinates(p:p+2) = [2.330882, 6.942920, 1.234863]
    p=p+3; coordinates(p:p+2) = [4.807957, 8.807140, 2.516223]
    p=p+3; coordinates(p:p+2) = [7.985039, 6.117336, 8.217075]
    p=p+3; coordinates(p:p+2) = [2.905626, 8.443855, 8.342717]
    p=p+3; coordinates(p:p+2) = [6.511206, 4.540879, 4.235224]
    p=p+3; coordinates(p:p+2) = [3.179374, 6.079763, 5.367636]
    p=p+3; coordinates(p:p+2) = [1.559462, 2.958662, 5.703067]
    p=p+3; coordinates(p:p+2) = [6.633376, 7.629552, 5.628650]
    p=p+3; coordinates(p:p+2) = [6.828553, 4.238397, 0.699934]
    p=p+3; coordinates(p:p+2) = [6.237931, 0.963391, 5.973284]
    p=p+3; coordinates(p:p+2) = [5.216592, 7.232110, 6.300030]
    p=p+3; coordinates(p:p+2) = [5.658115, 3.747455, 1.690021]
    p=p+3; coordinates(p:p+2) = [3.447445, 2.264113, 3.993541]
    p=p+3; coordinates(p:p+2) = [0.092429, 4.250129, 3.121253]
    p=p+3; coordinates(p:p+2) = [4.528862, 6.986506, 8.257740]
    p=p+3; coordinates(p:p+2) = [2.576832, 7.443634, 7.377702]
    p=p+3; coordinates(p:p+2) = [1.253327, 0.067623, 4.493787]
    p=p+3; coordinates(p:p+2) = [8.132418, 2.516618, 2.655838]
    p=p+3; coordinates(p:p+2) = [0.210942, 5.384717, 2.093258]
    p=p+3; coordinates(p:p+2) = [2.913229, 0.250497, 1.189487]
    p=p+3; coordinates(p:p+2) = [6.185676, 1.688229, 7.250712]
    p=p+3; coordinates(p:p+2) = [2.902332, 1.466218, 0.286110]
    p=p+3; coordinates(p:p+2) = [1.245411, 3.326594, 0.639775]
    p=p+3; coordinates(p:p+2) = [2.572067, 2.925108, 8.069143]
    p=p+3; coordinates(p:p+2) = [1.322880, 6.470741, 3.860597]
    p=p+3; coordinates(p:p+2) = [5.264514, 9.069863, 8.103554]
    p=p+3; coordinates(p:p+2) = [8.150170, 9.056487, 8.570771]
    p=p+3; coordinates(p:p+2) = [7.491957, 4.821872, 6.509134]
    p=p+3; coordinates(p:p+2) = [4.934731, 3.720091, 8.183964]
    p=p+3; coordinates(p:p+2) = [8.908281, 4.875424, 8.159178]
    p=p+3; coordinates(p:p+2) = [8.552147, 1.681953, 0.668343]
    p=p+3; coordinates(p:p+2) = [5.430863, 5.755537, 8.159240]
    p=p+3; coordinates(p:p+2) = [8.809227, 7.767207, 4.935412]
    p=p+3; coordinates(p:p+2) = [5.119268, 5.022251, 4.049013]
    p=p+3; coordinates(p:p+2) = [4.050740, 5.258455, 6.435706]
    p=p+3; coordinates(p:p+2) = [6.905047, 3.555702, 5.933401]
    p=p+3; coordinates(p:p+2) = [3.688555, 2.147422, 2.492974]
    p=p+3; coordinates(p:p+2) = [5.559858, 1.997668, 3.219644]
    p=p+3; coordinates(p:p+2) = [0.530272, 6.507162, 6.392404]
    p=p+3; coordinates(p:p+2) = [3.896900, 5.464608, 1.289558]
    p=p+3; coordinates(p:p+2) = [2.624854, 7.025622, 3.336256]
    p=p+3; coordinates(p:p+2) = [5.102760, 3.053809, 6.774712]
    p=p+3; coordinates(p:p+2) = [8.567045, 1.930168, 6.827999]
    p=p+3; coordinates(p:p+2) = [3.074774, 8.793845, 4.429651]
    p=p+3; coordinates(p:p+2) = [8.086846, 6.450836, 4.475729]
    p=p+3; coordinates(p:p+2) = [6.474042, 6.007170, 2.529696]
    p=p+3; coordinates(p:p+2) = [0.056309, 0.759693, 5.001359]
    p=p+3; coordinates(p:p+2) = [6.321956, 0.640013, 3.540305]
    p=p+3; coordinates(p:p+2) = [8.766794, 0.052826, 2.350809]
    p=p+3; coordinates(p:p+2) = [7.871889, 8.140986, 1.965417]
    p=p+3; coordinates(p:p+2) = [1.103365, 7.813345, 1.087900]
    p=p+3; coordinates(p:p+2) = [3.076808, 9.113733, 5.896554]
    p=p+3; coordinates(p:p+2) = [2.712701, 3.930271, 5.259852]
    p=p+3; coordinates(p:p+2) = [0.518010, 7.354247, 7.538435]
    p=p+3; coordinates(p:p+2) = [8.718131, 3.095734, 6.028403]
    p=p+3; coordinates(p:p+2) = [7.293086, 8.224417, 7.484842]
    p=p+3; coordinates(p:p+2) = [6.237757, 6.802342, 1.137021]
    p=p+3; coordinates(p:p+2) = [5.216726, 7.616771, 3.304717]
    p=p+3; coordinates(p:p+2) = [0.054559, 1.765512, 3.379760]
    p=p+3; coordinates(p:p+2) = [1.603831, 4.789926, 0.369489]

    print*, "things in coordinates should be",32*3*3, "and it is: ", p+2
    
    cell(1) = 35.000
    cell(2) = 35.000
    cell(3) = 35.000

    forces(:) = 0.0

    ! Call the scme function.
    do i=1,100
       call scme_calculate(n_atoms, coordinates, cell, forces, u_tot)
    end do

  end subroutine 

  subroutine test_scme_cluster_64_perf()
    implicit none
    ! ----------------------------------------
    integer, parameter :: n_atoms = 64*3 !32 water * 3
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
           coordinates(p:p+2) = [ -6.31504583 ,   -6.96366310,    -2.7262795]
    p=p+3; coordinates(p:p+2) = [ -8.45501137 ,   -4.28040600,    4.57241678]
    p=p+3; coordinates(p:p+2) = [ -3.34254265 ,   -7.53874588,    -4.8175392]
    p=p+3; coordinates(p:p+2) = [ -6.00126457 ,   2.73719716 ,   4.28764105 ]
    p=p+3; coordinates(p:p+2) = [ -3.32503390 ,   4.33411789 ,   -7.89316702]
    p=p+3; coordinates(p:p+2) = [ -10.86732292,    2.57142377,    6.23482895]
    p=p+3; coordinates(p:p+2) = [ -4.04160786 ,   -6.00639248,    2.82591510]
    p=p+3; coordinates(p:p+2) = [ 1.37160087  ,  -0.38430426 ,   2.60069180 ]
    p=p+3; coordinates(p:p+2) = [ 0.16538815  ,  3.71141911  ,  1.34374118  ]
    p=p+3; coordinates(p:p+2) = [ -4.60142469 ,   -3.04783869,    3.58867168]
    p=p+3; coordinates(p:p+2) = [ -8.47128773 ,   5.27224493 ,   -5.48457670]
    p=p+3; coordinates(p:p+2) = [ -2.40646720 ,   1.20477355 ,   -1.77750456]
    p=p+3; coordinates(p:p+2) = [ 2.99062514  ,  0.08082870  ,  -5.12260771 ]
    p=p+3; coordinates(p:p+2) = [ -2.17296195 ,   2.63227940 ,   0.51625723 ]
    p=p+3; coordinates(p:p+2) = [ -0.16663843 ,   6.50960827 ,   -3.78184891]
    p=p+3; coordinates(p:p+2) = [ -8.49094486 ,   7.69357872 ,   8.23789406 ]
    p=p+3; coordinates(p:p+2) = [ -4.60740566 ,   -2.26502395,    -2.9574639]
    p=p+3; coordinates(p:p+2) = [ -5.66186476 ,   -9.40780067,    -3.6680538]
    p=p+3; coordinates(p:p+2) = [ -8.64172840 ,   -6.90924168,    -1.0080360]
    p=p+3; coordinates(p:p+2) = [ -6.07779980 ,   -4.89796257,    5.49510241]
    p=p+3; coordinates(p:p+2) = [ 2.34769464  ,  -3.36739612 ,   -6.01602077]
    p=p+3; coordinates(p:p+2) = [ -3.76762342 ,   -7.55055523,    0.48845601]
    p=p+3; coordinates(p:p+2) = [ -3.43009281 ,   1.64490891 ,   4.79866838 ]
    p=p+3; coordinates(p:p+2) = [ -0.46168166 ,   -6.12163162,    1.09510243]
    p=p+3; coordinates(p:p+2) = [ -6.50911379 ,   -0.79509622,    -8.8170967]
    p=p+3; coordinates(p:p+2) = [ 7.49520111  ,  -4.91360664 ,   -1.69814694]
    p=p+3; coordinates(p:p+2) = [ 3.90259004  ,  -2.11219835 ,   2.64090061 ]
    p=p+3; coordinates(p:p+2) = [ -0.68540132 ,   9.46042824 ,   2.36137581 ]
    p=p+3; coordinates(p:p+2) = [ -0.18799661 ,   7.78617859 ,   -6.63894939]
    p=p+3; coordinates(p:p+2) = [ 4.54612160  ,  7.51634836  ,  0.64674801  ]
    p=p+3; coordinates(p:p+2) = [ 3.73654199  ,  0.48572090  ,  -1.15618587 ]
    p=p+3; coordinates(p:p+2) = [ -2.21734905 ,   8.33113861 ,   0.26419133 ]
    p=p+3; coordinates(p:p+2) = [ 4.05333185  ,  -1.98103845 ,   -0.02185570]
    p=p+3; coordinates(p:p+2) = [ 1.90192807  ,  -0.73786831 ,   -2.80972743]
    p=p+3; coordinates(p:p+2) = [ -4.14224911 ,   0.26880389 ,   -3.47064471]
    p=p+3; coordinates(p:p+2) = [ -0.72952807 ,   -1.66858363,    6.68440771]
    p=p+3; coordinates(p:p+2) = [ -2.27017903 ,   -6.87003994,    -2.0870184]
    p=p+3; coordinates(p:p+2) = [ 2.75970435  ,  -7.42528820 ,   -8.08861637]
    p=p+3; coordinates(p:p+2) = [ -2.25251484 ,   -2.58540106,    4.43371773]
    p=p+3; coordinates(p:p+2) = [ 6.79196739  ,  1.08504128  ,  6.35011005  ]
    p=p+3; coordinates(p:p+2) = [ 5.49075556  ,  5.64833307  ,  2.51119447  ]
    p=p+3; coordinates(p:p+2) = [ 2.94602489  ,  1.66261744  ,  -11.18212891]
    p=p+3; coordinates(p:p+2) = [ -0.73105466 ,   3.30720711 ,   9.47586536 ]
    p=p+3; coordinates(p:p+2) = [ -2.29456806 ,   2.04934812 ,   -5.04749346]
    p=p+3; coordinates(p:p+2) = [ 1.78723109  ,  3.34286904  ,  -1.94128072 ]
    p=p+3; coordinates(p:p+2) = [ 1.19450939  ,  7.50113201  ,  3.06086278  ]
    p=p+3; coordinates(p:p+2) = [ 0.25740579  ,  3.77765775  ,  4.20728016  ]
    p=p+3; coordinates(p:p+2) = [ 2.08474922  ,  7.29129648  ,  10.03464317 ]
    p=p+3; coordinates(p:p+2) = [ 6.13780832  ,  -3.29400706 ,   -4.73900652]
    p=p+3; coordinates(p:p+2) = [ 1.31295585  ,  -4.30244303 ,   0.10245909 ]
    p=p+3; coordinates(p:p+2) = [ 5.76587772  ,  -10.96664143,    1.63199544]
    p=p+3; coordinates(p:p+2) = [ 5.30121088  ,  -1.12653327 ,   6.33600140 ]
    p=p+3; coordinates(p:p+2) = [ 6.39001274  ,  4.81748056  ,  -6.36224413 ]
    p=p+3; coordinates(p:p+2) = [ 7.21732378  ,  -0.74047625 ,   0.21252051 ]
    p=p+3; coordinates(p:p+2) = [ 5.90417480  ,  3.87129521  ,  0.18074951  ]
    p=p+3; coordinates(p:p+2) = [ -2.98690438 ,   0.55685866 ,   2.40691853 ]
    p=p+3; coordinates(p:p+2) = [ -0.31108531 ,   -1.99076760,    -1.7842397]
    p=p+3; coordinates(p:p+2) = [ 7.31231642  ,  -3.75379443 ,   0.87101680 ]
    p=p+3; coordinates(p:p+2) = [ 1.21826243  ,  -0.85190046 ,   5.17803812 ]
    p=p+3; coordinates(p:p+2) = [ -1.06496561 ,   5.10914516 ,   6.14282465 ]
    p=p+3; coordinates(p:p+2) = [ -2.22679615 ,   8.84171581 ,   -2.69012761]
    p=p+3; coordinates(p:p+2) = [ 4.01270247  ,  2.12188435  ,  -3.72990346 ]
    p=p+3; coordinates(p:p+2) = [ 1.69374335  ,  4.52970505  ,  -4.40602541 ]
    p=p+3; coordinates(p:p+2) = [ 10.43322468 ,   -0.98995972,    0.03126749]
    p=p+3; coordinates(p:p+2) = [ -5.79312801 ,   -6.26504898,    -2.3759796]
    p=p+3; coordinates(p:p+2) = [ -8.41998100 ,   -3.59314179,    3.90820169]
    p=p+3; coordinates(p:p+2) = [ -3.37715197 ,   -8.40490437,    -4.4418573]
    p=p+3; coordinates(p:p+2) = [ -5.07266521 ,   2.91264296 ,   4.07020903 ]
    p=p+3; coordinates(p:p+2) = [ -2.59695578 ,   4.66874695 ,   -7.44677067]
    p=p+3; coordinates(p:p+2) = [ -10.30360985,    1.91717947,    6.68361759]
    p=p+3; coordinates(p:p+2) = [ -3.95338082 ,   -5.06667852,    2.94704604]
    p=p+3; coordinates(p:p+2) = [ 1.23877645  ,  -0.39831364 ,   3.55441070 ]
    p=p+3; coordinates(p:p+2) = [ 0.13583852  ,  3.45060563  ,  2.28717899  ]
    p=p+3; coordinates(p:p+2) = [ -5.00832415 ,   -3.60508657,    4.28174829]
    p=p+3; coordinates(p:p+2) = [ -8.37912846 ,   6.16458035 ,   -5.07643795]
    p=p+3; coordinates(p:p+2) = [ -2.49952078 ,   0.55142862 ,   -1.11293638]
    p=p+3; coordinates(p:p+2) = [ 3.43661666  ,  0.89750516  ,  -4.74682856 ]
    p=p+3; coordinates(p:p+2) = [ -2.32637024 ,   2.28041053 ,   -0.40385184]
    p=p+3; coordinates(p:p+2) = [ -0.62454700 ,   7.07884359 ,   -4.35806751]
    p=p+3; coordinates(p:p+2) = [ -9.03492355 ,   8.28426933 ,   7.74227190 ]
    p=p+3; coordinates(p:p+2) = [ -3.72664261 ,   -2.64869761,    -2.9965117]
    p=p+3; coordinates(p:p+2) = [ -5.46084690 ,   -9.33606815,    -4.5957007]
    p=p+3; coordinates(p:p+2) = [ -8.22439480 ,   -6.27162790,    -0.4155678]
    p=p+3; coordinates(p:p+2) = [ -6.96759605 ,   -4.64862728,    5.15916824]
    p=p+3; coordinates(p:p+2) = [ 2.04887176  ,  -2.44794869 ,   -6.21551323]
    p=p+3; coordinates(p:p+2) = [ -3.15312243 ,   -7.25187683,    -0.1666311]
    p=p+3; coordinates(p:p+2) = [ -4.03584433 ,   1.57588899 ,   5.53983307 ]
    p=p+3; coordinates(p:p+2) = [ 0.10668758  ,  -5.51180315 ,   0.71371299 ]
    p=p+3; coordinates(p:p+2) = [ -7.34802055 ,   -1.32536137,    -9.0757808]
    p=p+3; coordinates(p:p+2) = [ 7.40926695  ,  -4.66315508 ,   -0.78496742]
    p=p+3; coordinates(p:p+2) = [ 3.17503476  ,  -1.48133922 ,   2.79959226 ]
    p=p+3; coordinates(p:p+2) = [ -0.05739218 ,   8.75185108 ,   2.42537689 ]
    p=p+3; coordinates(p:p+2) = [ 0.24250983  ,  7.51276350  ,  -7.39614630 ]
    p=p+3; coordinates(p:p+2) = [ 5.39381075  ,  8.05878639  ,  0.61367309  ]
    p=p+3; coordinates(p:p+2) = [ 3.20211983  ,  1.06745338  ,  -0.67476559 ]
    p=p+3; coordinates(p:p+2) = [ -3.12351322 ,   8.30922318 ,   0.39998388 ]
    p=p+3; coordinates(p:p+2) = [ 4.00769091  ,  -1.10421813 ,   -0.39931673]
    p=p+3; coordinates(p:p+2) = [ 2.45116353  ,  -0.14531421 ,   -2.29062319]
    p=p+3; coordinates(p:p+2) = [ -4.90046406 ,   0.75393385 ,   -3.82656622]
    p=p+3; coordinates(p:p+2) = [ -0.91604143 ,   -1.57507992,    7.63071060]
    p=p+3; coordinates(p:p+2) = [ -2.97441816 ,   -6.76624346,    -2.7185390]
    p=p+3; coordinates(p:p+2) = [ 3.18323731  ,  -7.34467363 ,   -7.22243023]
    p=p+3; coordinates(p:p+2) = [ -1.48194635 ,   -2.76315689,    3.86267686]
    p=p+3; coordinates(p:p+2) = [ 6.49984837  ,  0.16109729  ,  6.31747103  ]
    p=p+3; coordinates(p:p+2) = [ 4.73815966  ,  5.30496264  ,  2.96629930  ]
    p=p+3; coordinates(p:p+2) = [ 2.38359404  ,  2.36589193  ,  -11.38939476]
    p=p+3; coordinates(p:p+2) = [ -0.46115825 ,   3.61055779 ,   8.62746811 ]
    p=p+3; coordinates(p:p+2) = [ -2.26730871 ,   1.84387732 ,   -5.96944904]
    p=p+3; coordinates(p:p+2) = [ 1.90705061  ,  3.89254904  ,  -1.13334179 ]
    p=p+3; coordinates(p:p+2) = [ 0.77951699  ,  6.69345808  ,  2.85593605  ]
    p=p+3; coordinates(p:p+2) = [ -0.25663465 ,   4.22898912 ,   4.86591101 ]
    p=p+3; coordinates(p:p+2) = [ 1.32062173  ,  7.13611937  ,  9.52436733  ]
    p=p+3; coordinates(p:p+2) = [ 6.34290266  ,  -3.81580305 ,   -5.48040247]
    p=p+3; coordinates(p:p+2) = [ 1.94986438  ,  -4.48227501 ,   0.80672282 ]
    p=p+3; coordinates(p:p+2) = [ 4.84435844  ,  -11.00938797,    1.62436080]
    p=p+3; coordinates(p:p+2) = [ 4.53303432  ,  -0.59858799 ,   6.61539459 ]
    p=p+3; coordinates(p:p+2) = [ 6.31594324  ,  4.06403399  ,  -6.94876432 ]
    p=p+3; coordinates(p:p+2) = [ 6.57920837  ,  -0.37624958 ,   -0.44805548]
    p=p+3; coordinates(p:p+2) = [ 5.21987009  ,  4.09243488  ,  -0.45769536 ]
    p=p+3; coordinates(p:p+2) = [ -2.52341366 ,   1.29267693 ,   2.00923681 ]
    p=p+3; coordinates(p:p+2) = [ -0.97196662 ,   -1.71122682,    -1.1059416]
    p=p+3; coordinates(p:p+2) = [ 7.06028128  ,  -2.97903299 ,   0.28558564 ]
    p=p+3; coordinates(p:p+2) = [ 0.38258818  ,  -1.00604582 ,   5.63397503 ]
    p=p+3; coordinates(p:p+2) = [ -1.87265229 ,   5.07118177 ,   6.69338799 ]
    p=p+3; coordinates(p:p+2) = [ -2.21766090 ,   8.00881004 ,   -2.25659990]
    p=p+3; coordinates(p:p+2) = [ 4.88136101  ,  2.60390139  ,  -3.73689842 ]
    p=p+3; coordinates(p:p+2) = [ 1.04253161  ,  5.26815844  ,  -4.25164032 ]
    p=p+3; coordinates(p:p+2) = [ 10.18956375 ,   -0.70881796,    0.94120270]
    p=p+3; coordinates(p:p+2) = [ -7.14357233 ,   -6.71134567,    -2.3292105]
    p=p+3; coordinates(p:p+2) = [ -9.00186634 ,   -3.98211813,    5.31658316]
    p=p+3; coordinates(p:p+2) = [ -4.16632652 ,   -7.42499590,    -5.2672638]
    p=p+3; coordinates(p:p+2) = [ -6.37168217 ,   2.60768437 ,   3.37017345 ]
    p=p+3; coordinates(p:p+2) = [ -3.58132029 ,   5.06191063 ,   -8.46547031]
    p=p+3; coordinates(p:p+2) = [ -10.81357193,    3.34719133,    6.83979940]
    p=p+3; coordinates(p:p+2) = [ -5.01617479 ,   -6.24317312,    2.96309853]
    p=p+3; coordinates(p:p+2) = [ 0.73124522  ,  -0.91313595 ,   2.24173427 ]
    p=p+3; coordinates(p:p+2) = [ -0.65650594 ,   3.40389252 ,   1.04012930 ]
    p=p+3; coordinates(p:p+2) = [ -5.01135921 ,   -2.13793111,    3.67012882]
    p=p+3; coordinates(p:p+2) = [ -9.21504116 ,   4.93016529 ,   -5.02637339]
    p=p+3; coordinates(p:p+2) = [ -3.14344025 ,   1.09464312 ,   -2.35185099]
    p=p+3; coordinates(p:p+2) = [ 2.50228548  ,  -0.40059093 ,   -4.42516851]
    p=p+3; coordinates(p:p+2) = [ -2.70722580 ,   3.42956281 ,   0.66376978 ]
    p=p+3; coordinates(p:p+2) = [ -0.81140471 ,   6.08111143 ,   -3.13051558]
    p=p+3; coordinates(p:p+2) = [ -7.68518734 ,   8.22941208 ,   8.26448250 ]
    p=p+3; coordinates(p:p+2) = [ -5.14192247 ,   -2.80897188,    -3.5367167]
    p=p+3; coordinates(p:p+2) = [ -5.68517351 ,   -8.49637032,    -3.3568542]
    p=p+3; coordinates(p:p+2) = [ -9.29443645 ,   -6.46964645,    -1.4587142]
    p=p+3; coordinates(p:p+2) = [ -6.02043962 ,   -5.87251663,    5.58129311]
    p=p+3; coordinates(p:p+2) = [ 1.44258010  ,  -3.85883570 ,   -6.05587196]
    p=p+3; coordinates(p:p+2) = [ -3.77164173 ,   -6.90572596,    1.22618222]
    p=p+3; coordinates(p:p+2) = [ -3.34224772 ,   2.55351996 ,   4.67301035 ]
    p=p+3; coordinates(p:p+2) = [ -0.14936532 ,   -7.03945589,    0.97468501]
    p=p+3; coordinates(p:p+2) = [ -6.32388306 ,   -0.16061966,    -9.5393419]
    p=p+3; coordinates(p:p+2) = [ 7.69696808  ,  -4.07921362 ,   -2.07061100]
    p=p+3; coordinates(p:p+2) = [ 3.89378977  ,  -2.14871716 ,   1.62508380 ]
    p=p+3; coordinates(p:p+2) = [ -1.25110734 ,   9.28205872 ,   1.62822235 ]
    p=p+3; coordinates(p:p+2) = [ -0.73031974 ,   8.50955963 ,   -6.82543325]
    p=p+3; coordinates(p:p+2) = [ 4.68234873  ,  7.04517221  ,  1.48108900  ]
    p=p+3; coordinates(p:p+2) = [ 4.22656679  ,  0.95676386  ,  -1.84796762 ]
    p=p+3; coordinates(p:p+2) = [ -1.93985915 ,   7.40440893 ,   0.32687563 ]
    p=p+3; coordinates(p:p+2) = [ 3.65147948  ,  -2.51237321 ,   -0.69508928]
    p=p+3; coordinates(p:p+2) = [ 1.17371511  ,  -0.94052315 ,   -2.24266744]
    p=p+3; coordinates(p:p+2) = [ -4.39947176 ,   -0.62374830,    -3.2033104]
    p=p+3; coordinates(p:p+2) = [ -1.56143308 ,   -1.84061241,    6.26524401]
    p=p+3; coordinates(p:p+2) = [ -1.83610511 ,   -7.74420929,    -2.2727930]
    p=p+3; coordinates(p:p+2) = [ 1.87388873  ,  -7.76346064 ,   -7.93528605]
    p=p+3; coordinates(p:p+2) = [ -3.13204074 ,   -2.77696514,    4.09514523]
    p=p+3; coordinates(p:p+2) = [ 6.47544050  ,  1.62105608  ,  5.58384228  ]
    p=p+3; coordinates(p:p+2) = [ 5.66859770  ,  5.10605860  ,  1.74727118  ]
    p=p+3; coordinates(p:p+2) = [ 2.37403274  ,  0.88985014  ,  -11.07854366]
    p=p+3; coordinates(p:p+2) = [ -1.22229552 ,   2.49213815 ,   9.33124828 ]
    p=p+3; coordinates(p:p+2) = [ -2.83287001 ,   1.40745139 ,   -4.59817171]
    p=p+3; coordinates(p:p+2) = [ 0.92618746  ,  2.91742206  ,  -1.82573462 ]
    p=p+3; coordinates(p:p+2) = [ 2.01779318  ,  7.18451595  ,  3.45678997  ]
    p=p+3; coordinates(p:p+2) = [ 0.80240744  ,  3.18302083  ,  4.75533628  ]
    p=p+3; coordinates(p:p+2) = [ 2.82533002  ,  7.41345215  ,  9.43700218  ]
    p=p+3; coordinates(p:p+2) = [ 6.05682564  ,  -2.42694092 ,   -5.18074989]
    p=p+3; coordinates(p:p+2) = [ 1.73553848  ,  -4.51202154 ,   -0.76760751]
    p=p+3; coordinates(p:p+2) = [ 5.89807510  ,  -10.19888783,    1.12287343]
    p=p+3; coordinates(p:p+2) = [ 5.25727987  ,  -0.92806637 ,   5.40444708 ]
    p=p+3; coordinates(p:p+2) = [ 5.47775412  ,  4.88759518  ,  -5.97809219 ]
    p=p+3; coordinates(p:p+2) = [ 7.13562298  ,  0.08063880  ,  0.73319989  ]
    p=p+3; coordinates(p:p+2) = [ 6.62958002  ,  4.44024515  ,  -0.09421976 ]
    p=p+3; coordinates(p:p+2) = [ -3.19699550 ,   0.77163476 ,   3.30854988 ]
    p=p+3; coordinates(p:p+2) = [ 0.16937001  ,  -2.67924809 ,   -1.37874877]
    p=p+3; coordinates(p:p+2) = [ 7.51060200  ,  -3.42906880 ,   1.75559950 ]
    p=p+3; coordinates(p:p+2) = [ 1.66010559  ,  -0.32027137 ,   5.83107996 ]
    p=p+3; coordinates(p:p+2) = [ -0.69557017 ,   5.98683119 ,   6.06865263 ]
    p=p+3; coordinates(p:p+2) = [ -1.50517571 ,   9.33919430 ,   -2.35033488]
    p=p+3; coordinates(p:p+2) = [ 3.45489883  ,  2.55304313  ,  -3.05991793 ]
    p=p+3; coordinates(p:p+2) = [ 1.78998995  ,  4.13418150  ,  -3.51006246 ]
    p=p+3; coordinates(p:p+2) = [ 9.63492584  ,  -1.43550730 ,   -0.22714230]

    print*, "things in coordinates should be",64*3*3, "and it is: ", p+2
    
    cell(1) = 35.000
    cell(2) = 35.000
    cell(3) = 35.000

    forces(:) = 0.0

    ! Call the scme function.
    do i=1,100
       call scme_calculate(n_atoms, coordinates, cell, forces, u_tot)
    end do

  end subroutine 

end module test_scme

