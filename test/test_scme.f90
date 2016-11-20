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
    real*8,dimension(45000)     :: coordinates
    real*8,dimension(3)         :: cell
    real*8,dimension(45000)     :: forces
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
    do i=1,10000
       call scme_calculate(n_atoms, coordinates, cell, forces, u_tot)
    end do

  end subroutine test_scme_perf

end module test_scme
