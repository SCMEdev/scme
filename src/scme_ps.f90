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

!should be on top of scme_calculate
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

module scme

  use data_types
  !use max_parameters
  use parameters, only: num_cells
  use multipole_parameters, only: d0, q0, o0, h0
  use polariz_parameters, only: dd0, dq0, hp0, qq0
  use molecProperties, only: recoverMolecules, calcCentersOfMass,&
       findPpalAxes, rotatePoles, rotatePolariz, setUnpolPoles, addFields, addDfields, create_rw, calc_cm
  use calc_lower_order, only: calcEdip_quad
  use calc_higher_order, only: calcEhigh
  use calc_derivs, only: calcDv
  use inducePoles, only: induceDipole, induceQpole
  !use forceCM_mod, only: forceCM
  !use torqueCM_mod, only: torqueCM
  
  use atomicForces_mod, only: atomicForces, stillAtomicForces
  
  use calcEnergy_mod, only: calcEnergy
  !use coreInt_mod, only: coreInt
  use sf_disp_tangtoe, only: dispersion, new_dispersion
  use force_torqueCM, only: forceCM, torqueCM
  
 ! the new PS surfaces: 
  use ps_dms, only: vibdms
  use ps_pes, only: vibpes
  use constants, only:A_a0, ea0_Deb, eA_Deb
  use printer_mod, only: printer, printer_h2o_linear, h2o_lin

  implicit none
  private
  public scme_calculate

contains

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
    !real(dp), parameter :: kk1 = ea0_Deb!2.5417709_dp !au2deb
    !real(dp), parameter :: au2deb = kk1 !1.88972612456506198632=A_a0
    !real(dp), parameter :: kk2 = A_a0    !1.88972666351031921149_dp !A2b
    real(dp), parameter :: kk1 = 2.5417709_dp
    real(dp), parameter :: kk2 = 1.88972666351031921149_dp
    real(dp), parameter :: convFactor = 14.39975841_dp / 4.803206799_dp**2 ! e**2[eVÅ]=[eVm / (e*10**10[statcol])**2
    real(dp), parameter :: rMax = 11.0_dp
    real(dp), parameter :: rMax2 = rMax*rMax
    integer, parameter :: NC = num_cells

    ! Parameter flags for controlling behavior.
    logical*1, parameter :: usePS_PES = .true.
    logical*1, parameter :: usePS_DMS = .true.
    logical*1, parameter :: iSlab = .false.
    !logical*1, parameter :: addCore = .false.

    ! Local flag for controlling behavior.
    logical*1, save :: converged

    ! Local variables for energies.
    real(dp), save :: uQ, uH, uES, uDisp, uD, uCore

    ! Local arrays for lattice and half lattice.
    real(dp), save :: a(3)
    real(dp), save :: a2(3)

    ! Center of mass forces and torque.
    real(dp) :: fCM(3,n_atoms/3) ! center of mass force
    real(dp) :: fsf(3,n_atoms/3) 
    real(dp) :: tau(3,n_atoms/3) !center of mass torque

    real(dp) :: ra(n_atoms*3)    ! atomic positioins in stupid format
    real(dp) :: rCM(3,n_atoms/3) ! center of mass positions
    real(dp) :: x(3,3,n_atoms/3) ! rotation matrix

    ! Electric fields.
    real(dp) :: eD(3,n_atoms/3) !
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
    
    ! Multipoles (0 referes to the unpolarized di- and quadrupole)
    real(dp) :: dpole0(3,n_atoms/3)
    real(dp) :: qpole0(3,3,n_atoms/3)
    real(dp) :: opole(3,3,3,n_atoms/3)
    real(dp) :: hpole(3,3,3,3,n_atoms/3)
    real(dp) :: dpole(3,n_atoms/3)
    real(dp) :: qpole(3,3,n_atoms/3)

    ! Polarizabilities (defined in polariz_parameters)
    real(dp) :: dd(3,3,n_atoms/3)
    real(dp) :: dq(3,3,3,n_atoms/3)
    real(dp) :: hp(3,3,3,n_atoms/3)
    real(dp) :: qq(3,3,3,3,n_atoms/3)

    ! Local integers.
    integer, save :: nM, nO, nH, i, p
    integer, save :: indO, indH1, indH2

    ! Local arrays for ??? ML
    !real(dp) :: uPES(n_atoms*3)
    !real(dp), save :: dipmom(3)
    !real(dp), save :: qdms(3)
    

    ! For Partridge-Schwenke surfaces
    !real(dp), save :: mol(9)
    !real(dp), save :: grad(9)
    !real(dp), save :: uPES1
   !new shit: 
    !integer jjj
    real(dp) ps_mol(3,3) !fix this  !!!  ! ps(O,H1,H1 ; x,y,z)
    real(dp) ps_mol_dip(3,3)
    real(dp) ps_grad(9)
    real(dp) ps_pes
    !real(dp) :: qdms(3)
    real(dp) :: dms(3)
    !real(dp), save :: dipmom2(3)
    integer iteration
    type(h2o) :: rw(n_atoms/3)
    real(dp) :: rwCM(3,n_atoms/3)
    integer m
    real(dp) :: uPES(n_atoms/3)
    type(h2o) :: aforces(n_atoms/3)
    real(dp) :: fa_test(n_atoms*3)
    
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
    call recoverMolecules(coords, ra, nH, nO, a, a2,   rw) !JÖ rw
    call create_rw(ra,rw,nM)

!call printer(ra, 'ra')
!call printer(rw, 'rw')

    call calc_cm(rw,rCM,nM)
!call printer(rwCM, 'rwCM')

    ! Calculate the center of mass for each molecule.
    !call calcCentersOfMass(ra, nM, rCM)
!call printer(rCM, 'rCM')



    ! Find the rotation matrix x that defines the principal axis.
    call findPpalAxes(ra, nM, x)

    ! Rotate the multipoles d0, dq etc, to allign with the molecules using
    ! the rotation matrix x. The result is stored in the dpole0 etc arrays.
    call rotatePoles(d0, q0, o0, h0, dpole0, qpole0, opole, hpole, nM, x)

    ! call Partridge-Schwenke dipole moment surface routine.
!    print*, "DIPOLE"
    if (usePS_DMS) then
       do m = 1,nM

          !OHH order in ps_dms: 
          ps_mol_dip(1,:) = rw(m)%o 
          ps_mol_dip(2,:) = rw(m)%h1 
          ps_mol_dip(3,:) = rw(m)%h2
          
          dms = 0
          call vibdms(ps_mol_dip,dms) 
          
          dpole0(:,m) = dms(:)*kk1*kk2!*eA_Deb! dipmom(p)*kk1*kk2!*ea0_Deb*A_a0!     When we fix units this shuld be fixed. 
       end do
    end if
call printer(dpole0, 'dpole0')
    
    !print*, sqrt(sum(dpole0(:,1)**2))
    !print*, sqrt(sum(dpole0(:,1)**2))
    
    ! NEEDS DOCUMENTATION
    call setUnpolPoles(dpole, qpole, dpole0, qpole0, nM)
    call rotatePolariz(dd0, dq0, qq0, hp0, dd, dq, qq, hp, nM, x)

    call calcEhigh(rCM, opole, hpole, nM, NC, a, a2, uH, eH, dEhdr, rMax2, iSlab)
call printer(eH, 'eH')
call printer(uH, 'uH')
call printer(dEhdr, 'dEhdr')
!call printer(opole, 'opole')
!call printer(hpole, 'hpole')


    ! Here's where the induction loop begins.
    converged = .false.
    iteration = 0
    do while (.not. converged)
    iteration = iteration + 1

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
!JÖ    call atomicForces(fCM, tau, ra, rCM, nM, fa)
    call stillAtomicForces(fCM,tau,rw,rCM,nM,aforces)

    ! Calculate the energy of interaction between the multipole moments
    ! and the electric field of the other molecules.
    call calcEnergy(dpole0, qpole0, opole, hpole, d1v, d2v, d3v, d4v, nM, u_tot)
    
    !/////////////////////////////////////////////////////////////////// CHanging units, why?
    ! Convert the forces and total energy to correct units.
    !do i = 1,(3*n_atoms)
    !   fa(i) = convFactor * fa(i)
    !end do
    do m = 1,nM
       aforces(m)%h1 = aforces(m)%h1*convFactor
       aforces(m)%h2 = aforces(m)%h2*convFactor
       aforces(m)%o = aforces(m)%o*convFactor
    enddo
!JÖ    fa(:) = fa(:)*convFactor
    u_tot = u_tot * convFactor

    ! Store the energy this far for debug printout.
    uES = u_tot

    ! Calculate dispersion forces ??? ML
    call new_dispersion(rw, aforces, uDisp, nM, a, a2)
!call printer(uDisp, 'uDisp new')    
!    call dispersion(ra, fa, uDisp, nM, a, a2)
    u_tot = u_tot + uDisp
!call printer(uDisp, 'uDisp old')    


!call printer(fa,'fa')    
!call printer(aforces,'aforces')    
    

    ! Calculate the core contribution to the energy. (only to the energy ??? ML)
    !if (addCore) then
    !   call coreInt(ra, fa, uCore, nM, a, a2)
    !   u_tot = u_tot + uCore
    !end if
    
    ! Adding intramolecular energy from Partridge-Schwenke PES.
    !uPES = 0
    !if (usePS_PES) then
       do m=1,nM
          
          !OHH order in vibpes
          ps_mol(1,:) = rw(m)%o!(indO  + p)
          ps_mol(2,:) = rw(m)%h1!(indH1 + p)
          ps_mol(3,:) = rw(m)%h2!(indH2 + p)

          
          ps_grad = 0
          call vibpes(ps_mol,ps_pes,ps_grad)!,ps_pes) *A2b
          
          u_tot = u_tot + ps_pes
          !uPES(m) = ps_pes


!JÖ          indH1 = 6*(m-1)
!JÖ          indO  = 3*(m-1+2*nM)
!JÖ          indH2 = 3+6*(m-1)
!JÖ          do p=1,3
!JÖ             
!JÖ             !
!JÖ             fa(indO  + p) = fa(indO  + p) - ps_grad(p)  !*h2eV*A2b
!JÖ             fa(indH1 + p) = fa(indH1 + p) - ps_grad(p+3)!*h2eV*A2b
!JÖ             fa(indH2 + p) = fa(indH2 + p) - ps_grad(p+6)!*h2eV*A2b
!JÖ             
!JÖ          end do
          
          aforces(m)%o  = aforces(m)%o   - ps_grad(1:3)  
          aforces(m)%h1 = aforces(m)%h1  - ps_grad(4:6)
          aforces(m)%h2 = aforces(m)%h2  - ps_grad(7:9)
          
       end do
    !end if
    
    !fa_test = fa
!call printer(fa_test,'fa_test')
    !fa=0
    !call h2o_lin(aforces,fa,nM)
    !call h2o_lin(aforces,fa,nM)
!call printer(fa,'fa')
fa_test=0
call printer_h2o_linear(aforces,'printer_h2o_linear(aforces)')
call h2o_lin(aforces,fa,nM)
call printer(fa,'printer( h2o_lin(aforces))')

call printer(aforces,'aforces')
    
    
call printer(u_tot,'u_tot')
!call printer_h2o_linear(aforces,'aforces linear')

!call printer(aforces,'aforces')


    return

  end subroutine scme_calculate

end module scme

!          indH1 = 6*(i-1)
!          indO  = 3*(i-1+2*nM)
!          indH2 = 3+6*(i-1)
!
!          ! Get x, y, z coordinates for H and O.
!          do p=1,3
!             ps_mol_dip(1,p) = ra(indO  + p)
!             ps_mol_dip(2,p) = ra(indH1 + p)
!             ps_mol_dip(3,p) = ra(indH2 + p)
!          end do


! fran vibpes loop
          !indH1 = 6*(m-1)
          !indO  = 3*(m-1+2*nM)
          !indH2 = 3+6*(m-1)
          !do p=1,3
          !   
          !   !
          !   fa(indO  + p) = fa(indO  + p) - ps_grad(p)  !*h2eV*A2b
          !   fa(indH1 + p) = fa(indH1 + p) - ps_grad(p+3)!*h2eV*A2b
          !   fa(indH2 + p) = fa(indH2 + p) - ps_grad(p+6)!*h2eV*A2b
          !   
          !end do
