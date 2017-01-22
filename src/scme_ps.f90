! Copyright (c)  2015-2016  SCMEdev.
! Licenced under the LGPLv3 license. See LICENSE for details.

!_______________________________________________________________________
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

! NB!!!_______________________________________________________________________
! This verison of the code is under development and experimental. 
! It is not equivalent to the code the paper. 
! If you want to use it contact Me at jonatan.ostrom@gmail.com


!UNITS__________________________________________________________________
! The code uses the following units:
! QUANTITY     SYMBOL  (NAME/COMMENT
! charge:           e  (electron charge
! potential:        V  (Volt
! energy:          eV  (electronvolt
! distance:         A  (Angstrom
! dipole:          eA  (
! n-poles:      eA**n  (n=2,3,4 for quadru-, octu- and hexadeca-pole resp.
! mass:             u  (atomic mass unit, probably not important?
! coulomb_k:  14.399_  (27.211_ * 0.529_ = (Eh/eV)*(a0/A) 
!---------
! Conversion factors are from CODATA2014, NIST

module scme

  ! Parameters:
  use data_types, only:A_a0,dp,h2o, num_cells, coulomb_k !, ea0_Deb, eA_Deb,kk1,kk2 
  use multipole_parameters, only: q0, o0, h0 !d0, 
  use polariz_parameters, only: dd0, dq0, hp0, qq0
  
  ! Routines:
  use molecProperties, only: recoverMolecules, rotatePolariz,add_field_gradients, rotate_qoh_poles !rotatePoles,setUnpolPoles,findPpalAxes, calcCentersOfMass,addFields, 
  use calc_lower_order, only: dip_quadField
  use calc_higher_order, only: octu_hexaField
  use calc_derivs, only: calcDv
  use inducePoles, only: induce_dipole, induce_quadrupole
  
  use calcEnergy_mod, only: multipole_energy
  use sf_disp_tangtoe, only:oxygen_dispersion ! dispersion, new_dispersion,
  use force_torqueCM, only: forceCM, torqueCM
  
  use ps_dms, only: vibdms
  use ps_pes, only: vibpes
  use printer_mod, only: printer, h2o_to_linear, xyz_hho_to_linear !printer_h2o_linear, 
  
  use localAxes_mod, only:localAxes2, create_xyz_hho, get_cm,force_and_torque_on_atoms !localAxes, force_torqueOnAtoms,create_rw, calc_cm, create_xyz, 

  implicit none
  private
  public scme_calculate

contains !//////////////////////////////////////////////////////////////

  subroutine scme_calculate(n_atoms, coords, lattice, fa, u_tot)

    implicit none
    integer, intent(in) :: n_atoms
    real(dp), intent(in) :: coords(n_atoms*3)
    real(dp), intent(in) :: lattice(3)
    real(dp), intent(out) :: fa(n_atoms*3)
    real(dp), intent(out) :: u_tot
    ! ----------------------------------------

    ! Constants and parameters.
    real(dp), parameter :: rMax = 11.0_dp
    real(dp), parameter :: rMax2 = rMax*rMax
    integer, parameter :: NC = num_cells

    ! Parameter flags for controlling behavior.
    logical*1, parameter :: iSlab = .false.
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
    
    ! Multipoles (0 referes to the unpolarized di- and quadrupole)
    real(dp) :: dpole(3,n_atoms/3)      ,dpole0(3,n_atoms/3)
    real(dp) :: qpole(3,3,n_atoms/3)    ,qpole0(3,3,n_atoms/3)
    real(dp) :: opole(3,3,3,n_atoms/3)
    real(dp) :: hpole(3,3,3,3,n_atoms/3)

    ! Polarizabilities (defined in polariz_parameters)
    real(dp) :: dd(3,3,n_atoms/3)
    real(dp) :: dq(3,3,3,n_atoms/3)
    real(dp) :: hp(3,3,3,n_atoms/3)
    real(dp) :: qq(3,3,3,3,n_atoms/3)

    ! Local integers.
    integer, save :: nM, nO, nH, i, p
    integer, save :: indO, indH1, indH2
    
    real(dp) ps_grad(3,3)!ps_grad(9)!
    real(dp) ps_pes
    real(dp) :: dms(3)
    integer iteration
    integer m
    real(dp) :: fa_test(n_atoms*3)
    real(dp) :: u_multipole, u_ps
    real(dp) :: xyz_hho(3,3,n_atoms/3) !xyz,hho,nM
    real(dp) :: xa_forces(3,3,n_atoms/3) !xyz,hho,nM
    
    
    !// Routine Starts /////////////////////////////////////////////////
    
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
    call recoverMolecules(coords, ra, nH, nO, a, a2) !JÖ rw
    call create_xyz_hho(ra,xyz_hho,nM)
    
    ! compute centers of mass (cm)
    call get_cm(xyz_hho,rCM,nM)
    
call printer(rCM, 'rCM')    

    !// Multipole Interaction //////////////////////////////////////////
    
    !/ Partridge-Schwenke Dipoles (PSD)
    do m = 1,nM
       dms = 0
       call vibdms(xyz_hho(:,:,m),dms) 
       dpole0(:,m) = dms(:)! *kk1*kk2! *eA_Deb!  (eA comes out)
    end do
    
    !/ Let PSD define the molecular (local) axes by the rotation matrix "x" 
    do m = 1,nM
!       call localAxes(dpole0(:,m),rw(m),x(:,:,m))
       call localAxes2(dpole0(:,m),xyz_hho(:,:,m),x(:,:,m))
    enddo
call printer(x,'x')    

    !/ Rotate the other poles into the local axes coordinate system defined by the dipole
    !call rotatePoles(d0, q0, o0, h0, dpole0, qpole0, opole, hpole, nM, x)
    call rotate_qoh_poles(q0, o0, h0, qpole0, opole, hpole, nM, x)
    
    !/ Save first reference value of dipole and quadrupole before the "induction loop"
    !call setUnpolPoles(dpole, qpole, dpole0, qpole0, nM)
    dPole = dpole0 !JÖ 
    qpole = qpole0 !JÖ
    
    !/ Rotate polarizability tensors into local coordinate systems
    call rotatePolariz(dd0, dq0, qq0, hp0, dd, dq, qq, hp, nM, x) !0=nonrotated
    
    !/ Compute the electric field (F) and gradient (dF) for the static octupole and hexadecapole
    ! why dont we induce those, we have the polarizabilities!?
    call octu_hexaField(rCM, opole, hpole, nM, NC, a, a2, uH, eH, dEhdr, rMax2, iSlab) 
    ! output: uH=scalar energy; eH(3,nM)=field from q,h; dEhdr(3,3,nM)=field gradient

    
    !/ Induce dipole and quadrupole to self consistency
    converged = .false.
    iteration = 0
    do while (.not. converged)
    iteration = iteration + 1

       call dip_quadField(rCM, dpole, qpole, nM, NC, a, a2, uD, uQ, eD, dEddr, rMax2, iSlab)
       ! output: uD,uQ=scalar energies; eD(3,nM)=d+q field; dEddr(3,3,nM)=d+q filed gradient
       
       !call addFields(eH, eD, eT, nM)
       eT = eH + eD !add fields
       call add_field_gradients(dEhdr, dEddr, dEtdr, nM) !dEtdr = dEhdr + dEddr !add field gradients
       

       ! Induce dipoles and quadrupoles.
       converged = .true.

       call induce_dipole(dpole, dpole0, eT, dEtdr, dd, dq, hp, nM, converged)
       call induce_quadrupole(qpole, qpole0, eT, dEtdr, dq, qq, nM, converged)
    end do
    
    !/ Compute filed gradients of the electric fields, to 5th order
    call calcDv(rCM, dpole, qpole, opole, hpole, nM, NC, a, a2, d1v, d2v, d3v, d4v, d5v, rMax2, fsf, iSlab)

    !/ Compute the forces on the centers of mass
    call forceCM(dpole, qpole, opole, hpole, d2v, d3v, d4v, d5v, nM, fsf, fCM)

    !/ Compute the torques on centers of mass
    call torqueCM(dpole, qpole, opole, hpole, d1v, d2v, d3v, d4v, nM, tau)
    
    
    !/ Multipole forces on atoms 
    !/ Compute force on atoms from CM force/torque as "temp. rigid body"
    do m = 1,nM
       xa_forces(:,:,m)=0
       call force_and_torque_on_atoms(tau(:,m),fCM(:,m),xyz_hho(:,:,m), xa_forces(:,:,m), rCM(:,m))
    enddo
    xa_forces=xa_forces*coulomb_k ! Coulomb force constant to get eV/A
    
    
    !/ Mutipole energy
    call multipole_energy(dpole0, qpole0, opole, hpole, d1v, d2v, d3v, d4v, nM, u_multipole)
    u_multipole = u_multipole * coulomb_k ! Coulomb force constant to get eV
    
    
    
    
    
    
    !// Dispersion /////////////////////////////////////////////////////
    !call new_dispersion(rw, aforces, uDisp, nM, a, a2)
    call oxygen_dispersion(xyz_hho, xa_forces, uDisp, nM, a, a2) !uDisp created here
    
    
    !// Partridge-Schwenke PES /////////////////////////////////////////
    u_ps=0
    do m=1,nM
       
       call vibpes(xyz_hho(:,:,m),ps_pes,ps_grad)!,ps_pes) *A2b
       
       u_ps = u_ps + ps_pes
       
       xa_forces(:,:,m) = xa_forces(:,:,m) - ps_grad
       
    end do
    
    !// Output /////////////////////////////////////////////////////////
    u_tot = u_multipole + uDisp + u_ps      !Total system energy (output)
    call xyz_hho_to_linear(xa_forces,fa,nM) !Total forces in fa(nM*9) (output)

    
    
    
    !// Debug output ///////////////////////////////////////////////////!(pipe to file, diff to see change, comment to mute)
    call printer(u_multipole,'u_multipole')
    call printer(uDisp,'uDisp')
    call printer(u_ps,'u_ps')
    call printer(u_tot,'u_tot')
    
    
    !call h2o_to_linear(aforces,fa_test,nM)
    !call printer(fa_test,'aforces linear')
    
    call printer(fa,'xa_forces linear')
    
    
    



    return

  end subroutine scme_calculate

end module scme

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

  !use max_parameters
  !use parameters, only: num_cells
  !use forceCM_mod, only: forceCM
  !use torqueCM_mod, only: torqueCM
  
  !use atomicForces_mod, only: atomicForces, stillAtomicForces
  !use coreInt_mod, only: coreInt
  !use constants, only:A_a0, ea0_Deb, eA_Deb


    !logical*1, parameter :: usePS_PES = .true.
    !logical*1, parameter :: usePS_DMS = .true.
    !logical*1, parameter :: addCore = .false.
    ! Local flag for controlling behavior.


   !new shit: 
    !integer jjj
!    real(dp) ps_mol(3,3) !fix this  !!!  ! ps(O,H1,H1 ; x,y,z)
!    real(dp) ps_mol_dip(3,3)
    !real(dp) :: qdms(3)
    !real(dp), save :: dipmom2(3)
!    type(h2o) :: rw(n_atoms/3)
!    real(dp) :: rwCM(3,n_atoms/3)
!    real(dp) :: uPES(n_atoms/3)
!    type(h2o) :: aforces(n_atoms/3)
!   real(dp) :: rot(3,3,n_atoms/3) ! rotation matrix
!    real(dp) :: xyz(3,n_atoms)

    !call calc_cm(rw,rCM,nM)
    !call create_rw(ra,rw,nM)
!    call create_xyz(ra,xyz,nM)


    ! Local arrays for ??? ML
    !real(dp) :: uPES(n_atoms*3)
    !real(dp), save :: dipmom(3)
    !real(dp), save :: qdms(3)
    
    ! For Partridge-Schwenke surfaces
    !real(dp), save :: mol(9)
    !real(dp), save :: grad(9)
    !real(dp), save :: uPES1


       !OHH order in ps_dms: 
       !ps_mol_dip(:,1) = xyz_hho(:,1,m)!rw(m)%h1 
       !ps_mol_dip(:,2) = xyz_hho(:,2,m)!rw(m)%h2
       !ps_mol_dip(:,3) = xyz_hho(:,3,m)!rw(m)%o 
       !call vibdms(ps_mol_dip,dms) 

       !OHH order in vibpes
       !ps_mol(:,1) = xyz_hho(:,1,m)!rw(m)%h1!(indH1 + p)
       !ps_mol(:,2) = xyz_hho(:,2,m)!rw(m)%h2!(indH2 + p)
       !ps_mol(:,3) = xyz_hho(:,3,m)!rw(m)%o!(indO  + p)
       !call vibpes(ps_mol,ps_pes,ps_grad)!,ps_pes) *A2b


       !uPES(m) = ps_pes
    
    
!JÖ       indH1 = 6*(m-1)
!JÖ       indO  = 3*(m-1+2*nM)
!JÖ       indH2 = 3+6*(m-1)
!JÖ       do p=1,3
!JÖ          
!JÖ          !
!JÖ          fa(indO  + p) = fa(indO  + p) - ps_grad(p)  !*h2eV*A2b
!JÖ          fa(indH1 + p) = fa(indH1 + p) - ps_grad(p+3)!*h2eV*A2b
!JÖ          fa(indH2 + p) = fa(indH2 + p) - ps_grad(p+6)!*h2eV*A2b
!JÖ          
!JÖ       end do


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
          

!       aforces(m)%h1 = 0
!       aforces(m)%h2 = 0
!       aforces(m)%o  = 0
!       call force_torqueOnAtoms(tau(:,m),fCM(:,m),rw(m), aforces(m), rCM(:,m))
!       !/ Convert from stat-units to eV & eV/A
!       aforces(m)%h1 = aforces(m)%h1*coulomb_k
!       aforces(m)%h2 = aforces(m)%h2*coulomb_k
!       aforces(m)%o  = aforces(m)%o* coulomb_k
       
          
       !xa_forces(:,1,m) = xa_forces(:,1,m)   - ps_grad(1:3)
       !xa_forces(:,2,m) = xa_forces(:,2,m)   - ps_grad(4:6)
       !xa_forces(:,3,m) = xa_forces(:,3,m)   - ps_grad(7:9)
       
       !aforces(m)%h1 = aforces(m)%h1  - ps_grad(1:3)
       !aforces(m)%h2 = aforces(m)%h2  - ps_grad(4:6)
       !aforces(m)%o  = aforces(m)%o   - ps_grad(7:9)  
       
       !aforces(m)%h1 = aforces(m)%h1  - ps_grad(:,1)!1:3)
       !aforces(m)%h2 = aforces(m)%h2  - ps_grad(:,2)!4:6)
       !aforces(m)%o  = aforces(m)%o   - ps_grad(:,3)!7:9)  
       
       
       
       
       
