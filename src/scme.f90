! Copyright (c)  2015-2016  SCMEdev.
! Licenced under the LGPLv3 license. See LICENSE for details.


!> Module for calculating energy and forces of water moecules using the SCME
!! potential that is based on multipole moments, where ach molecule is
!! represented as a multipole expansion up to hexadecapole moment.
!!
!! The SCME potential is described in:
!!
!!      K. T. Wikfeldt, E. R. Batista, F. D. Vila and H. Jonsson
!!      Phys. Chem. Chem. Phys., 2013, 15, 16542
!!
!! Please cite this work if you use the SCME potential in you research.
module scme

  use data_types
  use max_parameters
  use parameters, only: num_cells
  use multipole_parameters, only: d0, q0, o0, h0
  use polariz_parameters, only: dd0, dq0, hp0, qq0
  use molecProperties, only: recoverMolecules, calcCentersOfMass,&
       findPpalAxes, rotatePoles, rotatePolariz, setUnpolPoles, addFields, addDfields
  use calc_lower_order, only: calcEdip_quad
  use calc_higher_order, only: calcEhigh
  use calc_derivs, only: calcDv
  use inducePoles, only: induceDipole, induceQpole
  use forceCM_mod, only: forceCM
  use torqueCM_mod, only: torqueCM
  use atomicForces_mod, only: atomicForces
  use calcEnergy_mod, only: calcEnergy
  use coreInt_mod, only: coreInt
  use dispersion_mod, only: dispersion

  implicit none
  private
  public scme_calculate

contains

  !> The main routine for the SCME potential. Calculates the total energy and
  !! forces on a set of water molecules.
  !!
  !! @param[in] n_atoms : The number of atoms. The number of atoms is assumed
  !!                      to be 3 times the number of water molecules.
  !! @param[in] coords  : The coordinates of the water molecules.
  !!                  coords(l+6*(i-1)) stores the l-th coordinate of the first
  !!                      hydrogen in the i-th molecule
  !!                  coords(l+3+6*(i-1)) stores the l-th coordinate of the second
  !!                      hydrogen in the i-th molecule
  !!                  coords(l+3*(i-1+nHydrogens)) stores the l-th coordinate of the
  !!                      oxygen in the i-th molecule. (nHydrogens is
  !!                      the total number of hydrogen atoms).
  !! @param[in] lattice : The x,y,z dimensions of the rectangular box.
  !! @param[out] fa     : The forces. The array must have space for n_atoms forces.
  !! @param[out] u_tot  : The total energy calculated with the SCME potential.
  subroutine scme_calculate(n_atoms, coords, lattice, fa, u_tot)

    implicit none
    integer, intent(in) :: n_atoms
    real(dp), intent(in) :: coords(n_atoms*3)
    real(dp), intent(in) :: lattice(3)
    real(dp), intent(out) :: fa(n_atoms*3)
    real(dp), intent(out) :: u_tot
    ! ----------------------------------------

    ! Constants and parameters.
    real(dp), parameter :: pi = 3.14159265358979324_dp
    real(dp), parameter :: kk1 = 2.5417709_dp
    real(dp), parameter :: kk2 = 1.88972666351031921149_dp
    real(dp), parameter :: convFactor = 14.39975841_dp / 4.803206799_dp**2
    real(dp), parameter :: rMax = 11.0_dp
    real(dp), parameter :: rMax2 = rMax*rMax
    integer, parameter :: NC = num_cells

    ! Parameter flags for controlling behavior.
    logical*1, parameter :: irigidmolecules = .false.
    logical*1, parameter :: debug = .false.
    logical*1, parameter :: iSlab = .false.
    logical*1, parameter :: addCore = .false.
    logical*1, parameter :: useDMS = .true.

    ! Local flag for controlling behavior.
    logical*1 :: converged

    ! Local variables for energies.
    real(dp) :: uQ, uH, uES, uDisp, uD, uCore

    ! Local arrays for lattice and half lattice.
    real(dp) :: a(3)
    real(dp) :: a2(3)

    ! Center of mass forces and torque.
    real(dp) :: fCM(3,n_atoms/3)
    real(dp) :: fsf(3,n_atoms/3)
    real(dp) :: tau(3,n_atoms/3)

    ! Atomic positions, centers of mass, principal axes.
    real(dp) :: ra(n_atoms*3) 
    real(dp) :: rCM(3,n_atoms/3)
    real(dp) :: x(3,3,n_atoms/3)

    ! Electric fields.
    real(dp) :: eD(3,n_atoms/3)
    real(dp) :: eQ(3,n_atoms/3)
    real(dp) :: eH(3,n_atoms/3)
    real(dp) :: eT(3,n_atoms/3)

    ! Derivatives of E.
    real(dp) :: dEddr(3,3,n_atoms/3)
    real(dp) :: dEqdr(3,3,n_atoms/3)
    real(dp) :: dEhdr(3,3,n_atoms/3)
    real(dp) :: dEtdr(3,3,n_atoms/3)

    ! High order derivatives of the potential.
    real(dp) :: d1v(3,n_atoms/3)
    real(dp) :: d2v(3,3,n_atoms/3)
    real(dp) :: d3v(3,3,3,n_atoms/3)
    real(dp) :: d4v(3,3,3,3,n_atoms/3)
    real(dp) :: d5v(3,3,3,3,3,n_atoms/3)

    ! Work multipoles. They start unpolarized and with the induction
    ! loop we induce dipoles and quadrupoles.
    real(dp) :: dpole0(3,n_atoms/3)
    real(dp) :: qpole0(3,3,n_atoms/3)
    real(dp) :: opole(3,3,3,n_atoms/3)
    real(dp) :: hpole(3,3,3,3,n_atoms/3)
    real(dp) :: dpole(3,n_atoms/3)
    real(dp) :: qpole(3,3,n_atoms/3)

    ! Polarizabilities.
    real(dp) :: dd(3,3,n_atoms/3)
    real(dp) :: dq(3,3,3,n_atoms/3)
    real(dp) :: hp(3,3,3,n_atoms/3)
    real(dp) :: qq(3,3,3,3,n_atoms/3)

    ! Local integers.
    integer :: nM, nO, nH, i, p
    integer :: indO, indH1, indH2

    ! Local arrays for ??? ML
    real(dp) :: uPES(n_atoms*3)
    real(dp) :: qdms(3)
    real(dp) :: dipmom(3)

    ! Input for the potnasa potential.
    real(dp) :: mol(9)
    real(dp) :: grad(9)
    real(dp) :: uPES1
    ! ----------------------------------------


    ! ----------------------------
    ! Set initial Intitial values.
    ! ----------------------------

    ! Total energy.
    u_tot = 0.0_dp

    ! Number of oxygen, hydrogen and moecules.
    nO = n_atoms/3
    nH = nO*2
    nM = nO

    ! Size of the simulation cell.
    a(1) = lattice(1)
    a(2) = lattice(2)
    a(3) = lattice(3)
    a2(1) = lattice(1)/2.0_dp
    a2(2) = lattice(2)/2.0_dp
    a2(3) = lattice(3)/2.0_dp

    ! Recover broken molecules due to periodic boundary conditions.
    call recoverMolecules(coords, ra, nH, nO, a, a2)

    ! Calculate the center of mass for each molecule.
    call calcCentersOfMass(ra, nM, rCM)

    ! Find the rotation matrix x that defines the principal axis.
    call findPpalAxes(ra, nM, x)

    ! Rotate the multipoles d0, dq etc, to allign with the molecules using
    ! the rotation matrix x. The result is stored in the dpole0 etc arrays.
    call rotatePoles(d0, q0, o0, h0, dpole0, qpole0, opole, hpole, nM, x)

    ! call Partridge-Schwenke dipole moment surface routine.
    if (useDMS) then
       do i = 1,nM
          indH1 = 6*(i-1)
          indO  = 3*(i-1+2*nM)
          indH2 = 3+6*(i-1)

          ! Get x, y, z coordinates for H and O.
          do p=1,3
             mol(p) = ra(indO  + p)
             mol(p+3) = ra(indH1  + p)
             mol(p+6) = ra(indH2  + p)
          end do
          call dmsnasa2(mol,qdms)

          ! Calculate dipole moment wrt center of mass.
          dipmom(:) = 0.0_dp
          do p=1,3
             dipmom(p) = dipmom(p) + qdms(1)*(mol(p)-rCM(p,i))
             dipmom(p) = dipmom(p) + qdms(2)*(mol(p+3)-rCM(p,i))
             dipmom(p) = dipmom(p) + qdms(3)*(mol(p+6)-rCM(p,i))
          end do

          ! Set unpolarized dipoles to Partridge-Schwenke dipoles
          ! using conversion constants for eA -> D.
          do p=1,3
             dpole0(p,i) = dipmom(p)*kk1*kk2
          end do
       end do
    end if

    ! NEEDS DOCUMENTATION
    call setUnpolPoles(dpole, qpole, dpole0, qpole0, nM)
    call rotatePolariz(dd0, dq0, qq0, hp0, dd, dq, qq, hp, nM, x)

    call calcEhigh(rCM, opole, hpole, nM, NC, a, a2, uH, eH, dEhdr, rMax2, iSlab)

    ! Here's where the induction loop begins.
    converged = .false.
    do while (.not. converged)

       ! NEEDS DOCUMENTATION
       call calcEdip_quad(rCM, dpole, qpole, nM, NC, a, a2, uD, uQ, eD, dEddr, rMax2, iSlab)

       call addFields(eH, eD, eT, nM)
       call addDfields(dEhdr, dEddr, dEtdr, nM)

       ! Induce dipoles and quadrupoles.
       converged = .true.
       call induceDipole(dpole, dpole0, eT, dEtdr, dd, dq, hp, nM, converged)
       call induceQpole(qpole, qpole0, eT, dEtdr, dq, qq, nM, converged)

    end do

    ! With the polarized multipoles, calculate the derivarives of the
    ! electrostatic potential, up to 5th order.
    call calcDv(rCM, dpole, qpole, opole, hpole, nM, NC, a, a2, d1v, d2v, d3v, d4v, d5v, rMax2, fsf, iSlab)

    ! Compute the force on the center of mass.
    call forceCM(dpole, qpole, opole, hpole, d2v, d3v, d4v, d5v, nM, fsf, fCM)

    ! Compute the torque on the molecule.
    call torqueCM(dpole, qpole, opole, hpole, d1v, d2v, d3v, d4v, nM, tau)

    ! Find 3 forces, one at oxygen and one at each hydrogen such that
    ! the total force and total torque agree with the ones calculated
    ! for the multipoles.
    ! NOTE: This is where 'fa' is calculated.
    call atomicForces(fCM, tau, ra, rCM, nM, fa)

    ! Calculate the energy of interaction between the multipole moments
    ! and the electric field of the other molecules.
    call calcEnergy(dpole0, qpole0, opole, hpole, d1v, d2v, d3v, d4v, nM, u_tot)

    ! Convert the forces and total energy to correct units.
    do i = 1,(3*n_atoms)
       fa(i) = convFactor * fa(i)
    end do
    u_tot = u_tot * convFactor

    ! Store the energy this far for debug printout.
    uES = u_tot

    ! Calculate dispersion forces ??? ML
    call dispersion(ra, fa, uDisp, nM, a, a2)
    u_tot = u_tot + uDisp

    ! Calculate the core contribution to the energy. (only to the energy ??? ML)
    if (addCore) then
       call coreInt(ra, fa, uCore, nM, a, a2)
       u_tot = u_tot + uCore
    end if

    ! Adding intramolecular energy from Partridge-Schwenke PES.
    uPES(:) = 0.0_dp
    if (.not. irigidmolecules) then
       do i=1,nM
          mol(:) = 0.0_dp
          grad(:) = 0.0_dp

          indH1 = 6*(i-1)
          indO  = 3*(i-1+2*nM)
          indH2 = 3+6*(i-1)
          ! Get x, y, z coordinates for H and O.
          do p=1,3
             mol(p) = ra(indO  + p)
             mol(p+3) = ra(indH1  + p)
             mol(p+6) = ra(indH2  + p)
          end do

          call potnasa2(mol,grad,uPES1)
          uPES(i) = uPES1
          u_tot = u_tot + uPES1
          do p=1,3
             fa(indO  + p) = fa(indO + p) - grad(p)
             fa(indH1 + p) = fa(indH1 + p) - grad(p+3)
             fa(indH2 + p) = fa(indH2 + p) - grad(p+6)
          end do
       end do
    end if


! ML: These print-statements makes up half the CPU usage for a
!     dimer input.
!
!ktw    print '(5f16.10)', uTot, uES, uDisp, uCore, sum(uPES)
!    print '(4f16.10)', uTot, uES, uDisp, sum(uPES)
!    print*, size(ra), "is size ra" !JÖ

    return

  end subroutine scme_calculate

end module scme

