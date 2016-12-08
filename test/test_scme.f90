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

end module test_scme
