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
  use multipole_parameters, only: q0, o0, h0 !d0, 
  use polariz_parameters, only: dd0, dq0, hp0, qq0
  use molecProperties, only: recoverMolecules, &
       rotatePolariz,  addDfields, &
       rotate_qoh_poles !rotatePoles,setUnpolPoles,findPpalAxes, calcCentersOfMass,addFields, 
  use calc_lower_order, only: dip_quadField
  use calc_higher_order, only: octu_hexaField
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
  use localAxes_mod, only:create_rw, calc_cm, force_torqueOnAtoms, localAxes

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
    real(dp) :: rot(3,3,n_atoms/3) ! rotation matrix
    real(dp) :: u_multipole, u_ps
    ! ----------------------------
    ! Set initial Intitial values.
    ! ----------------------------
    ! Total energy.

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
    
    ! compute centers of mass (cm)
    call calc_cm(rw,rCM,nM)

    !// Multipole Expansion Part ///////////////////////////////////////
    
    !/ Partridge-Schwenke Dipoles (PSD)
    do m = 1,nM
    
       !OHH order in ps_dms: 
       ps_mol_dip(1,:) = rw(m)%o 
       ps_mol_dip(2,:) = rw(m)%h1 
       ps_mol_dip(3,:) = rw(m)%h2
       
       dms = 0
       call vibdms(ps_mol_dip,dms) 
       
       dpole0(:,m) = dms(:)*kk1*kk2!*eA_Deb! dipmom(p)*kk1*kk2!*ea0_Deb*A_a0!     When we fix units this shuld be fixed. 
    end do
    
    !/ Let PSD define the molecular (local) axes by the rotation matrix "x" 
    do m = 1,nM
       call localAxes(dpole0(:,m),rw(m),x(:,:,m))
    enddo
    
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
    call printer(uH*convFactor,'uH')
    
    !/ Induce dipole and quadrupole to self consistency
    converged = .false.
    iteration = 0
    do while (.not. converged)
    iteration = iteration + 1

       ! NEEDS DOCUMENTATION
       call dip_quadField(rCM, dpole, qpole, nM, NC, a, a2, uD, uQ, eD, dEddr, rMax2, iSlab)
       ! output: uD,uQ=scalar energies; eD(3,nM)=d+q field; dEddr(3,3,nM)=d+q filed gradient
       
       !call addFields(eH, eD, eT, nM)
       eT = eH + eD !add fields
       call addDfields(dEhdr, dEddr, dEtdr, nM) !dEtdr = dEhdr + dEddr !add field gradients
       

       ! Induce dipoles and quadrupoles.
       converged = .true.

       call induceDipole(dpole, dpole0, eT, dEtdr, dd, dq, hp, nM, converged)
       call induceQpole(qpole, qpole0, eT, dEtdr, dq, qq, nM, converged)
    end do
    
    !/ Compute filed gradients of the electric fields, to 5th order
    call calcDv(rCM, dpole, qpole, opole, hpole, nM, NC, a, a2, d1v, d2v, d3v, d4v, d5v, rMax2, fsf, iSlab)

    !/ Compute the forces on the centers of mass
    call forceCM(dpole, qpole, opole, hpole, d2v, d3v, d4v, d5v, nM, fsf, fCM)

    !/ Compute the torques on centers of mass
    call torqueCM(dpole, qpole, opole, hpole, d1v, d2v, d3v, d4v, nM, tau)
    !// End Multipoles /////////////////////////////////////////////////
    
    
    !// Forces on atoms ////////////////////////////////////////////////
    !/ Compute force atoms from CM force/torque as "temp. rigid body"
    do m = 1,nM
       aforces(m)%h1 = 0
       aforces(m)%h2 = 0
       aforces(m)%o  = 0
       call force_torqueOnAtoms(tau(:,m),fCM(:,m),rw(m), aforces(m), rCM(:,m))
       !/ Convert from stat-units to eV & eV/A
       aforces(m)%h1 = aforces(m)%h1*convFactor
       aforces(m)%h2 = aforces(m)%h2*convFactor
       aforces(m)%o  = aforces(m)%o*convFactor
    enddo
!call printer(aforces,'aforces')    
    
    
    !/ Calculate the energy of electric field interactions
    call calcEnergy(dpole0, qpole0, opole, hpole, d1v, d2v, d3v, d4v, nM, u_multipole)
    !/ Convert from stat-units to eV & eV/A
    u_multipole = u_multipole * convFactor
    
    
    !// Dispersion /////////////////////////////////////////////////////
    call new_dispersion(rw, aforces, uDisp, nM, a, a2)
    
    
    !// Partridge-Schwenke PES /////////////////////////////////////////
    u_ps=0
    do m=1,nM
       
       !OHH order in vibpes
       ps_mol(1,:) = rw(m)%o!(indO  + p)
       ps_mol(2,:) = rw(m)%h1!(indH1 + p)
       ps_mol(3,:) = rw(m)%h2!(indH2 + p)
       
       ps_grad = 0;ps_pes = 0
       call vibpes(ps_mol,ps_pes,ps_grad)!,ps_pes) *A2b
       
       u_ps = u_ps + ps_pes
       
       aforces(m)%o  = aforces(m)%o   - ps_grad(1:3)  
       aforces(m)%h1 = aforces(m)%h1  - ps_grad(4:6)
       aforces(m)%h2 = aforces(m)%h2  - ps_grad(7:9)
       
    end do
    
    u_tot = u_tot + uDisp + u_multipole + u_ps
    
    call printer(u_multipole,'u_multipole')
    call printer(uDisp,'uDisp')
    call printer(u_ps,'u_ps')
    call printer(u_tot,'u_tot')
    
    call h2o_lin(aforces,fa,nM)
    call printer(aforces,'aforces')



    return

  end subroutine scme_calculate

end module scme

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
