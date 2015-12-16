!-----------------------------------------------------------------------
!                  New potential based on multipole moments.
!     Each molecule is represented as a multipole expansion up to
!     hexadecapole moment.
!-----------------------------------------------------------------------
!     An arrays with positions enters as argument ra().
!     It is assumed that the array of positions ra() has a length equal
!     to 9 times the number of molecules and that its structure is the
!     following:
!           ra(l+6*(i-1)) stores the l-th coordinate of the first
!                         hydrogen in the i-th molecule
!           ra(l+3+6*(i-1)) stores the l-th coordinate of the second
!                         hydrogen in the i-th molecule
!           ra(l+3*(i-1+nHydrogens)) stores the l-th coordinate of the
!                         oxygen in the i-th molecule. (nHydrogens is
!                         the total number of hydrogen atoms)
!     The routine will return an array fa() with the forces on each atom.
!     The array fa() has the same structure as ra().
!     The total potential energy of the configuration is also calculated
!     and returned in 'uTot'.
!
!-----------------------------------------------------------------------

!23456789012345678901234567890123456789012345678901234567890123456789012
!        10        20        30        40        50        60        70

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
  
  !ktw  subroutine scme_calculate(nAtms, raOri, lattice, itagl, fa, uTot, virial)
  subroutine scme_calculate(nAtms, raOri, lattice, fa, uTot)
    ! nAtms: number of atoms, integer
    ! raOri: atom coordinates
    ! lattice: unit cell. lattice = [ax, ay, az]. Only supports rectangular unit cell(?)
    ! fa: forces
    ! uTot: total energy
    ! virial: not implemented
    
    implicit none
    
    !timing
    integer*8 ti, tf, irtc
    real(dp) t1, t2, t3, t4, t5, t6, t7, rMax, rMax2
    real(dp) d1, d2
    
    real(dp) pi, raOri(maxCoo), a(3), a2(3), fa(maxCoo)
    real(dp), dimension(3) :: lattice
    real(dp) fCM(3,maxCoo/3), fsf(3,maxCoo/3), tau(3,maxCoo/3)
    real(dp) uD, uQ, uH, uPol, uDisp, uRho, uCore, uES
    
    !     atomic position, centers of mass, principal axes.
    real(dp) ra(maxCoo), rCM(3,maxCoo/3), x(3,3,maxCoo/3)
    
    !     Electric fields
    real(dp) eD(3,maxCoo/3), eQ(3,maxCoo/3), eH(3,maxCoo/3)
    real(dp) eT(3,maxCoo/3)
    !     Derivatives of E
    real(dp) dEddr(3,3,maxCoo/3), dEqdr(3,3,maxCoo/3)
    real(dp) dEhdr(3,3,maxCoo/3), dEtdr(3,3,maxCoo/3)
    
    !     High order derivatives of the potential
    real(dp) d1v(3,maxCoo/3), d2v(3,3,maxCoo/3), d3v(3,3,3,maxCoo/3)
    real(dp) d4v(3,3,3,3,maxCoo/3), d5v(3,3,3,3,3,maxCoo/3)
    
    real(dp) uTot, faB(maxCoo), u, virial, convFactor
    real(dp) uCov, uAng, Amp, rMin, lambda, gam, rOHmax
    integer nM, nO, nH, nAtms, i, ii, j, k, l, NC, nst, Np
    
    !     Unpolarized multipoles
    !real(dp) d0(3), q0(3,3), o0(3,3,3), h0(3,3,3,3) ! in pol_params
    
    !     Work multipoles. They start unpolarized and with the induction
    !     loop we induce dipoles and quadrupoles.
    real(dp) dpole(3,maxCoo/3), qpole(3,3,maxCoo/3)
    real(dp) opole(3,3,3,maxCoo/3), hpole(3,3,3,3,maxCoo/3)    
    real(dp) dpole0(3,maxCoo/3), qpole0(3,3,maxCoo/3)
    
    !     Polarizabilities
    !real(dp) dd0(3,3), dq0(3,3,3), hp0(3,3,3), qq0(3,3,3,3) !defined in polariz_parameters
    real(dp) dd(3,3,maxCoo/3), dq(3,3,3,maxCoo/3), hp(3,3,3,maxCoo/3)
    real(dp) qq(3,3,3,3,maxCoo/3)
    
    !     Move along tetrahedral axis
    real(dp) m(3,3), mI(3,3), ft(3), dft(3), df(3), ftddr, ddd
    
    !     Study bulk modulus
    real(dp) raOld(maxCoo), scale, dr(3)
    integer iO, iH1, iH2
    !     Check derivatives
    integer at, coord, iMol, iComp
    real(dp) r1(3), r2(3), ro(3)
    real(dp) dx, di(3), fexa, fnum, uOld, r0(3)
    
!ktw    logical*1 firstVisit, converged, debug, iSlab
    logical*1 converged, debug, iSlab
    logical*1 addCore
    integer indO, indH1, indH2
    real(dp) ROH1, ROH2, uPES(maxCoo), psinp(maxCoo,3)
    real(dp) qdms(3), dipmom(3)
    real(dp) mol(9), grad(9), uPES1
!ktw    logical irigidmolecules
    logical*1 irigidmolecules, useDMS
    
    
    ! Debug
    integer     p, q, r, s
    real(dp)      kk1, kk2
    
    irigidmolecules = .false.
!    irigidmolecules = .true.
    useDMS = .true.
!    useDMS = .false.
    debug = .false.
    pi = 3.14159265358979312d0
    !call readPoles(d0, q0, o0, h0) ! now in module multipole_parameters
    !call readPolariz(dd0, dq0, hp0, qq0) ! now in module polariz_parameters
    NC = num_cells       ! Number of cell in each direction (= 0)
    rMax = 11.d0
    rMax2 = rMax*rMax
    !     Convert from Debye^2/angstrom^3 to eV
    convFactor = 14.39975841d0 / 4.803206799d0**2
    iSlab = .FALSE.
    !         iSlab = .TRUE.
!    addCore = .TRUE.
    addCore = .FALSE.
    kk1 = 2.5417709D0
    kk2 = 1.88972666351031921149D0

    !     Size of the simulation cell
    a(1) = lattice(1)
    a(2) = lattice(2)
    a(3) = lattice(3)
    a2(1) = lattice(1)/2
    a2(2) = lattice(2)/2
    a2(3) = lattice(3)/2
    nO = nAtms/3          ! Number of hydrogens
    nH = nO*2             ! Number of oxygens
    nM = nO                ! Number of molecules
    
    uTot = 0.0d0
    
    
    !     Recover broken molecules due to PBC
    call recoverMolecules(raOri, ra, nH, nO, a, a2)
    
    call calcCentersOfMass(ra, nM, rCM)
    
    call findPpalAxes(ra, nM, x)    ! defines rotation matrix(?) x
    
    call rotatePoles(d0, q0, o0, h0, dpole0, qpole0, opole, hpole, nM, x) !defines dpole0... using x
    
    ! call Partridge-Schwenke dipole moment surface routine
    if (useDMS) then
       do i = 1,nM
          indH1 = 6*(i-1)
          indO  = 3*(i-1+2*nM)
          indH2 = 3+6*(i-1)
          !KTW get x, y, z coordinates for H and O
          do p=1,3
             mol(p) = ra(indO  + p)
             mol(p+3) = ra(indH1  + p)
             mol(p+6) = ra(indH2  + p)
          end do
          call dmsnasa2(mol,qdms)
          
          !KTW calculate dipole moment wrt center of mass:
          dipmom(:) = 0.0d0
          do p=1,3
             dipmom(p) = dipmom(p) + qdms(1)*(mol(p)-rCM(p,i))
             dipmom(p) = dipmom(p) + qdms(2)*(mol(p+3)-rCM(p,i))
             dipmom(p) = dipmom(p) + qdms(3)*(mol(p+6)-rCM(p,i))
          end do
       
          !KTW set unpolarized dipoles to Partridge-Schwenke dipoles
          !KTW using conversion constants for eA -> D
          do p=1,3
             dpole0(p,i) = dipmom(p)*kk1*kk2
          end do
       end do
    end if

    call setUnpolPoles(dpole, qpole, dpole0, qpole0, nM)
    call rotatePolariz(dd0, dq0, qq0, hp0, dd, dq, qq, hp, nM, x)
    
    call calcEhigh(rCM, opole, hpole, nM, NC, a, a2, uH, eH, dEhdr, rMax2, iSlab)
    
    !     Here's where the induction loop begins
    converged = .false.
    
    do while (.not. converged)
       call calcEdip_quad(rCM, dpole, qpole, nM, NC, a, a2, uD, uQ, eD, dEddr, rMax2, iSlab)
       
       call addFields(eH, eD, eT, nM)
       call addDfields(dEhdr, dEddr, dEtdr, nM)
       
       !     Induce dipoles and quadrupoles
       converged = .true.
       call induceDipole(dpole, dpole0, eT, dEtdr, dd, dq, hp, nM, converged)
       call induceQpole(qpole, qpole0, eT, dEtdr, dq, qq, nM, converged)
       
    end do                    !End of the induction loop
    
!    d1 = sqrt(dpole0(1,1)**2+dpole0(2,1)**2+dpole0(3,1)**2)
!    print *, 'dip0= ', d1
!    d1 = sqrt(dpole(1,1)**2+dpole(2,1)**2+dpole(3,1)**2)
!    print *, 'dip= ', d1
    
    !     With the polarized multipoles, calculate the derivarives of the
    !     electrostatic potential, up to 5th order.
    
    call calcDv(rCM, dpole, qpole, opole, hpole, nM, NC, a, a2, d1v, d2v, d3v, d4v, d5v, rMax2, fsf, iSlab)
    
    !     Compute the force on the center of mass
    call forceCM(dpole, qpole, opole, hpole, d2v, d3v, d4v, d5v, nM, fsf, fCM)
    
    !     Compute the torque on the molecule
    call torqueCM(dpole, qpole, opole, hpole, d1v, d2v, d3v, d4v, nM, tau)
    
    !     Find 3 forces, one at oxygen and one at each hydrogen such that
    !     the total force and total torque agree with the ones calculated
    !     for the multipoles.    
    call atomicForces(fCM, tau, ra, rCM, nM, fa) ! fa calculated here
    
    !     Calculate the energy of interaction between the multipole moments
    !     and the electric field of the other molecules.
    call calcEnergy(dpole0, qpole0, opole, hpole, d1v, d2v, d3v, d4v, nM, uTot)
    
    !PRB - Don't use this one ---------------------------------------------
    !      call calcEnergy(dpole, qpole, opole, hpole, d1v, d2v, d3v, d4v,
    !    $      nM, uTot)
    !PRB - ----------------------------------------------------------------
    
    !PRB - Difference between energy of polarized and unpolarized molecules?
    !      uPol = 0.d0
    
    !      call calcEnergyI(dpole, dpole0, qpole, qpole0, opole,
    !     $     hpole, d1v, d2v, d3v, d4v, nM, uPol)
    !      print *, uTot*convFactor, uPol*convFactor, (uPol+uTot)*Convfactor
    !      stop
    
    !      print '(2f15.7,$)', uTot*convFactor, uPol*convFactor
    
    !      uTot = uTot - uPol
    
    !PRB - ----------------------------------------------------------------
    
    virial = 0.0d0
    do i = 1, 3*(nH+nO)
       fa(i) = convFactor * fa(i)
       !         virial = virial - ra(i)*fa(i)
    end do
    uTot = uTot * convFactor
    
    uES = uTot
    
    
    ! 112  call dispersion2(ra, fa, uDisp, nM, a, a2, NC, rMax2, iSlab)
    call dispersion(ra, fa, uDisp, nM, a, a2)
    uTot = uTot + uDisp
    
    if (addCore) then
       call coreInt(ra, fa, uCore, nM, a, a2)
       uTot = uTot + uCore
    end if
    
    !     KTW adding intramolecular energy from Partridge-Schwenke PES
    uPES(:) = 0.0d0
    if (.not. irigidmolecules) then
       do i=1,nM
          mol(:) = 0.0d0
          grad(:) = 0.0d0
          
          indH1 = 6*(i-1)
          indO  = 3*(i-1+2*nM)
          indH2 = 3+6*(i-1)
          !     KTW get x, y, z coordinates for H and O
          do p=1,3
             mol(p) = ra(indO  + p)
             mol(p+3) = ra(indH1  + p)
             mol(p+6) = ra(indH2  + p)
          end do
          
          call potnasa2(mol,grad,uPES1)
          uPES(i) = uPES1
          uTot = uTot + uPES1
          do p=1,3
             fa(indO  + p) = fa(indO + p) - grad(p)
             fa(indH1 + p) = fa(indH1 + p) - grad(p+3)
             fa(indH2 + p) = fa(indH2 + p) - grad(p+6)
          end do
          
!          write(*,*) 'forces from scme.f90'
!          write(*,*) 'O', fa(indO+1),fa(indO+2),fa(indO+3)
!          write(*,*) 'H', fa(indH1+1),fa(indH1+2),fa(indH1+3)
!          write(*,*) 'H', fa(indH2+1),fa(indH2+2),fa(indH2+3)

          
       end do
    end if
    
    
!ktw    print '(5f16.10)', uTot, uES, uDisp, uCore, sum(uPES)
    print '(4f16.10)', uTot, uES, uDisp, sum(uPES)

    
    return
    
  end subroutine scme_calculate
  
end module scme

