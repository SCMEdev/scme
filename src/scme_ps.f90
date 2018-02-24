
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
! It is not equivalent to the code in the paper. 
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
  use data_types, only:A_a0,dp, num_cells, coulomb_k !, ea0_Deb, eA_Deb,kk1,kk2    !,h2o
  use multipole_parameters, only: d0, q0, o0, h0 !d0, 
  use polariz_parameters, only: dd0, dq0, hp0, qq0
  
  ! Routines:
  use molecProperties, only: recoverMolecules, rotatePolariz,add_field_gradients, recoverMolecules_new, rotatePoles!rotate_qoh_poles,,setUnpolPoles,findPpalAxes, calcCentersOfMass,addFields, 
  use calc_lower_order, only: dip_quadField
  use calc_higher_order, only: octu_hexaField
  use calc_derivs, only: calcDv
  use inducePoles, only: induce_dipole, induce_quadrupole
  
  use calcEnergy_mod, only: multipole_energy
  use sf_disp_tangtoe, only:oxygen_dispersion ! dispersion, new_dispersion,
  use force_torqueCM, only: forceCM, torqueCM
  
  use ps_dms, only: vibdms
  use ps_pes, only: vibpes
  use printer_mod, only: printer, xyz_hho_to_linear,str,printo,printa !printer_h2o_linear, !, h2o_to_linear
  
  use localAxes_mod, only:dipoleAxes,plusAxes, bisectorAxes, &
                          get_cm,force_and_torque_on_atoms,create_xyz_hho_new 
                          !localAxes, force_torqueOnAtoms,create_rw, calc_cm, create_xyz, !, create_xyz_hho
  
  use qpole, only: get_quadrupoles, expansion_coordinates
  
  use opole, only: get_octupoles
  
  use compressed_utils,only: compress, expand
  use compressed_tensors, only:lin_df, lin_polydf, apple1_df, vector_powers, dfdu_erf, polyinner1, polyinner2, polyinner_matrix, get_stone_field, dfdu, system_stone_field
  use compressed_arrays
  use compressed_tests, only:polycompress_p
  !use detrace_apple, only: detrace_a, ff
  !use comressed_arrays
  
  implicit none
  private
  public scme_calculate

contains !//////////////////////////////////////////////////////////////

  subroutine scme_calculate(n_atoms &
                            ,coords &
                            ,lattice &
                            ,hho_fa &
                            ,u_tot &
                            ,USE_PS_PES &
                            ,USE_FULL_RANK &
                            ,USE_OO_REP &
                            ,USE_ALL_REP &
                            ,USE_VAR_QUAD &
                            ,USE_VAR_OCT &
                            ) 

    implicit none
    integer, intent(in) :: n_atoms
    !real(dp), intent(in) :: coords(n_atoms*3)
    real(dp), intent(in) :: coords(3,n_atoms)
    real(dp), intent(in) :: lattice(3)
    !j real(dp), intent(out) :: fa(n_atoms*3)
    real(dp), intent(out) :: hho_fa(3,n_atoms) !j Better. and: trying to get into quip
    real(dp), intent(out) :: u_tot
    !real(dp), intent(out), optional :: dip_perm(3,n_atoms/3), dip_ind(3,n_atoms/3)
    ! ----------------------------------------

    ! Constants and parameters.
    real(dp), parameter :: rMax = 11.0_dp
    real(dp), parameter :: rMax2 = rMax*rMax
    integer, parameter :: NC = num_cells

    ! Parameter flags for controlling behavior.
    logical*1, parameter :: iSlab = .false.
    logical*1, save :: converged

    
    

    ! Local variables for energies.
    real(dp), save :: uQ, uH, uDisp, uD!, uES, uCore

    ! Local arrays for lattice and half lattice.
    real(dp), save :: a(3)
    real(dp), save :: a2(3)

    ! Center of mass forces and torque.
    real(dp) :: fCM(3,n_atoms/3) ! center of mass force
    real(dp) :: fsf(3,n_atoms/3) 
    real(dp) :: tau(3,n_atoms/3) !center of mass torque

    !real(dp) :: ra(n_atoms*3)    ! atomic positioins in stupid format
    real(dp) :: ra(3,n_atoms)    ! atomic positioins in ...
    real(dp) :: rCM(3,n_atoms/3) ! center of mass positions
    real(dp) :: x(3,3,n_atoms/3) ! rotation matrix
    

    ! Electric fields.
    real(dp) :: eD(3,n_atoms/3) !
    !real(dp) :: eQ(3,n_atoms/3)
    real(dp) :: eH(3,n_atoms/3)
    real(dp) :: eT(3,n_atoms/3)

    ! Derivatives of E.
    real(dp) :: dEddr(3,3,n_atoms/3)
    !real(dp) :: dEqdr(3,3,n_atoms/3)
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
    integer, save :: nM, nO, nH!, i, p
    !integer, save :: indO, indH1, indH2
    
    real(dp) ps_grad(3,3)!ps_grad(9)!
    real(dp) ps_pes
    real(dp) :: dms(3)
    integer iteration
    integer m
    !real(dp) :: fa_test(n_atoms*3)
    real(dp) :: u_multipole, u_ps
    real(dp) :: xyz_hho(3,3,n_atoms/3) !xyz,hho,nM
    real(dp) :: xa_forces(3,3,n_atoms/3) !xyz,hho,nM
    
    
    integer ia
    
    
    real(dp) a_oo, b_oo, u_rep, rr_ox(3), r_ox
    integer ox1, ox2
    
    integer, parameter :: xyz = 3, hho = 3, rra = 3 !to keep track of indices
    integer combi, m1, m2
    real(dp) rr_oo(xyz), rr_oh(xyz,4), rr_hh(xyz,4), r_oo, r_oh, r_hh, aa_hh, aa_oh, aa_oo, A_oh, A_hh
    
    integer :: s !transposing
    !integer :: rank
    
    real(dp) :: cec(xyz,hho,n_atoms/3),cer2(hho,n_atoms/3), rCE(xyz,n_atoms/3)
    
    logical, intent(in), optional:: USE_PS_PES , USE_FULL_RANK , USE_OO_REP   , USE_ALL_REP, USE_VAR_QUAD, USE_VAR_OCT
    logical ::                      PES , FULL , OO_REP   , ALL_REP, VAR_QUAD, VAR_OCT
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Compressed tensors
    integer,parameter :: nx=4,kx=4, px=2, nkx=nx+kx, kpolx=2, pkx=pos_(kx+1), pnx=pos_(nx+1)
    integer nn, p1,p2
    
    real(dp), dimension(pos_(nx+1),n_atoms/3):: qn_scme!qn_perm, qn_tot, 
    
    real(dp), dimension(pos_(kx+1),n_atoms/3) :: phi_scme, phi_comp, f34!phi_perm, phi_tot, 
    real(dp), dimension(pos_(px+1),pos_(px+1),n_atoms/3) :: polz
    
    !real(dp), dimension(pos_(nkx+2)) :: rrr, df
    !real(dp) r2, rr(3), sss(nkx+1)
    
    !integer n1,n2,k1,k2
    
    !real(dp) tol, u_mult1,u_mult2, u_perm1,u_perm2
    !real(dp), dimension(pos_(kpolx+1),n_atoms/3) :: dqn,qn_pol
    !real(dp), dimension(pos_(kpolx+1),n_atoms/3) :: phi_pol, phi_pol2, dphi
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !Default optional arguments
    PES=.true.
    FULL=.true.
    OO_REP=.false.
    ALL_REP=.false.
    VAR_QUAD=.false.
    VAR_OCT=.false.
    
    
    if(present(USE_PS_PES))     PES          = USE_PS_PES
    if(present(USE_FULL_RANK))  FULL         = USE_FULL_RANK
    if(present(USE_OO_REP))     OO_REP       = USE_OO_REP
    if(present(USE_ALL_REP))    ALL_REP      = USE_ALL_REP
    if(present(USE_VAR_QUAD))   VAR_QUAD     = USE_VAR_QUAD
    if(present(USE_VAR_OCT))    VAR_OCT      = USE_VAR_OCT
    
    
    
    
    
    
    
    
    s=2
    
    
    
    
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
    
    
#include "debug.h"

    ! Recover broken molecules due to periodic boundary conditions.
    !call recoverMolecules(coords, ra, nH, nO, a, a2) !JÖ rw
    !call create_xyz_hho(ra,xyz_hho,nM)
    
    call recoverMolecules_new(coords, ra, nM, a, a2) !JÖ rw
    call create_xyz_hho_new(ra,xyz_hho,nM)
    
!tprint(coords, 'coords',s)    
!tprint(ra, 'ra',s)    
!tprint(xyz_hho, 'xyz_hho',s)    
    
    ! compute centers of mass (cm)
    call get_cm(xyz_hho,rCM,nM)
    call expansion_coordinates(xyz_hho,rCE,cec,cer2,nM,hho)
    
    
    
tprint(rCE, 'rCE',s)    
    
tprint(rCM, 'rCM',s)    

    !// Multipole Interaction //////////////////////////////////////////
    
    
    !dip_perm=dpole0 !/output
    
    ! Define local axes according to some rule
    do m = 1,nM
!       call localAxes(dpole0(:,m),rw(m),x(:,:,m))
       !call dipoleAxes(dpole0(:,m),xyz_hho(:,:,m),x(:,:,m))
       call bisectorAxes(xyz_hho(:,:,m),x(:,:,m))
       !call plusAxes(xyz_hho(:,:,m),x(:,:,m))
    enddo
tprint(x,'x',s)    

    !/ Rotate the other poles into the local axes coordinate system defined by the dipole
    
    call rotatePoles(d0, q0, o0, h0, dpole0, qpole0, opole, hpole, nM, x)
    !call rotate_qoh_poles(q0, o0, h0, qpole0, opole, hpole, nM, x)
    !opole_orig = opole
    
    !/ Partridge-Schwenke Dipoles (PSD)
    !do m = 1,nM
    !   dms = 0
    !   call vibdms(xyz_hho(:,:,m),dms) 
    !   dpole0(:,m) = dms(:)! *kk1*kk2! *eA_Deb!  (eA comes out)
    !end do
    
    
    
    !if(VAR_QUAD) 
    !if(VAR_QUAD)call get_quadrupoles(cec,cer2,rCE,qpole0,nM) !overwriting quadrupoles
    !if(VAR_OCT)call get_octupoles(cec,cer2,opole,nM) !overwriting octupoles
    
    
    !/ Save first reference value of dipole and quadrupole before the "induction loop"
    !call setUnpolPoles(dpole, qpole, dpole0, qpole0, nM)
    dPole = dpole0 !JÖ 
    qpole = qpole0 !JÖ
    
    
    qn_scme=0
    call convert_xpoles_to_qn_scme !/ NEW qn_scme
    
    
    
    
    
    
    
    
    !/ Rotate polarizability tensors into local coordinate systems
    call rotatePolariz(dd0, dq0, qq0, hp0, dd, dq, qq, hp, nM, x) !0=nonrotated
    
    
    
    !print*, "rCM",rCM, "opole", opole, "hpole",hpole, "nM",nM, "NC",NC, "a",a, "a2",a2, "uH",uH, "eH",eH, "dEhdr",dEhdr, "rMax2",rMax2, "iSlab",iSlab,"FULL",FULL
    
    
    
    !/ Compute the electric field (F) and gradient (dF) for the permanent octupole and hexadecapole
    ! why dont we induce those, we have the polarizabilities!?
    call octu_hexaField(rCM, opole, hpole, nM, NC, a, a2, uH, eH, dEhdr, rMax2, iSlab,FULL) 
    !! output: uH=scalar energy; eH(3,nM)=field from q,h; dEhdr(3,3,nM)=field gradient
    !call dip_quadField(rCM, dpole, qpole, nM, NC, a, a2, uD, uQ, eD, dEddr, rMax2, iSlab)
    
    !stop
    
    !call printa(eH)
    !call printer(dEhdr,"m",1)
    !stop
    
    f34=0
    print*, qn_scme
    
    call system_stone_field(3,4,1,2,6,rCM,qn_scme,f34)
    
    call printa(f34,t="f34")
    
    stop
    !/ Induce dipole and quadrupole to self consistency
    !converged = .false.
    !iteration = 0
    !do while (.not. converged)
    !iteration = iteration + 1
    !
    !   call dip_quadField(rCM, dpole, qpole, nM, NC, a, a2, uD, uQ, eD, dEddr, rMax2, iSlab)
    !   ! output: uD,uQ=scalar energies; eD(3,nM)=d+q field; dEddr(3,3,nM)=d+q filed gradient
    !   
    !   
    !   
    !   
    !   !call get_stone_field(
    !   !call addFields(eH, eD, eT, nM)
    !   eT = eH + eD !add fields
    !   call add_field_gradients(dEhdr, dEddr, dEtdr, nM) !dEtdr = dEhdr + dEddr !add field gradients
    !   
    !
    !   ! Induce dipoles and quadrupoles.
    !   converged = .true.
    !
    !   call induce_dipole(dpole, dpole0, eT, dEtdr, dd, dq, hp, nM, converged)
    !   call induce_quadrupole(qpole, qpole0, eT, dEtdr, dq, qq, nM, converged)
    !end do
    
    
    
    !dip_ind=dpole
    
    !/ Compute filed gradients of the electric fields, to 5th order
    !call system_stone_field(1,nx,1,kx,nM,rCM,qn_perm,phi_perm)
    call system_stone_field(1,nx,1,kx,nM,rCM,qn_scme,phi_comp)
    
    
    call calcDv(rCM, dpole, qpole, opole, hpole, nM, NC, a, a2, d1v, d2v, d3v, d4v, d5v, rMax2, fsf, iSlab,FULL)
    
    phi_scme=0
    call convert_d1234v_to_phi_scme
    
    print*, "phi COMPRESSED:"
    call printa(phi_comp) 
    
    
    print*, "phi SCME:"
    call printa(phi_scme)
    
    !call printa(phi_perm,2,5,0)
    call printa(phi_scme(2:,:)/phi_comp(2:,:),t="diff")
    
    stop
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    !stop"hej"
!tprint(d1v, 'd1v',s)
!tprint(d2v, 'd2v',s)
!tprint(d3v, 'd3v',s)
!tprint(d4v, 'd4v',s)
!tprint(d5v, 'd5v',s)
    !/ Compute the forces on the centers of mass
    call forceCM(dpole, qpole, opole, hpole, d2v, d3v, d4v, d5v, nM, fsf, fCM)

    !/ Compute the torques on centers of mass
    call torqueCM(dpole, qpole, opole, hpole, d1v, d2v, d3v, d4v, nM, tau)
    
tprint(fCM, 'fCM',s)
tprint(tau, 'tauu',s)
    
    !/ Multipole forces on atoms 
    !/ Compute force on atoms from CM force/torque as "temp. rigid body"
    do m = 1,nM
       xa_forces(:,:,m)=0
       call force_and_torque_on_atoms(tau(:,m),fCM(:,m),xyz_hho(:,:,m), xa_forces(:,:,m), rCM(:,m))
    enddo
    xa_forces=xa_forces*coulomb_k ! Coulomb force constant to get eV/A
    
    
    !/ Mutipole energy
    !/ This should in principle be the induced poles in the field of the static poles
    !/ What is computed is the static poles in the field of the induced poles. 
    !/ But since T is symmetric, this is probably equivalent. 
    !/ 
    !/ But for the forces, they are computed above with the static+induced poles in the static+induced field gradients. Should be static+induced poles in only static field gradients?? 
    !/ Because U_ind = Q_ind F_stat, U_stat = U_stat F_stat, so U_tot = Q_tot F_stat. 
    !/ solution: calcDv(dpole0,qpole0....
    !/           torqueCM(dpole, qpole, ... 
    !/           forceCM(dpole, qpole, ... unchanged syntax      (i.e. total pole in static field)
    !/           multipole_energy(dpole, qpole....            (i.e. total pole in static field)
    !m=4
    !print*, "hpole:"
    !call printo(hpole(:,:,:,:,m))
    !print*, "d4v:"
    !call printo(d4v(:,:,:,:,m))
    
    
    !u_mult1=0
    !u_mult2=0
    !u_perm1 =0
    !u_perm2 =0
    !do m=1,nM
    !    u_mult1 = u_mult1 + sum( qn_scme(:pkx,m)* phi_tot(:pkx,m)*gg_(:pkx) )
    !    u_mult2 = u_mult2 + sum(  qn_tot(:pkx,m)*phi_perm(:pkx,m)*gg_(:pkx) )
    !    u_perm1 = u_perm1 + sum( qn_perm(:,m)*phi_perm(:,m)*gg_(:pkx) )
    !    u_perm2 = u_perm2 + sum( qn_perm(:,m)*phi_perm(:,m)*gg_(:pkx) ) 
    !enddo
    
    
    !print*, 'u_mult1    ', u_mult1* coulomb_k
    !print*, 'u_mult2    ', u_mult2* coulomb_k
    !
    !print*, 'u_perm1    ', u_perm1* coulomb_k
    !print*, 'u_perm2    ', u_perm2* coulomb_k
    
    
    
    call multipole_energy(dpole0, qpole0, opole, hpole, d1v, d2v, d3v, d4v, nM, u_multipole)
    u_multipole = u_multipole * coulomb_k ! Coulomb force constant to get eV
    
    
    
    print*, 'u_multipole', u_multipole
    
    
    
    !// Dispersion /////////////////////////////////////////////////////
    !call new_dispersion(rw, aforces, uDisp, nM, a, a2)
    call oxygen_dispersion(xyz_hho, xa_forces, uDisp, nM, a, a2) !uDisp created here
    
    
    if(PES)then
        !// Partridge-Schwenke PES /////////////////////////////////////////
        u_ps=0
        do m=1,nM
           
           call vibpes(xyz_hho(:,:,m),ps_pes,ps_grad)!,ps_pes) *A2b
           
           u_ps = u_ps + ps_pes
           
           xa_forces(:,:,m) = xa_forces(:,:,m) - ps_grad
           
        end do
    endif!(USE_PS_PES)
    
    !// Output /////////////////////////////////////////////////////////
    u_tot = u_multipole + uDisp + u_ps      !Total system energy (output)
    !call xyz_hho_to_linear(xa_forces,fa,nM) !Total forces in fa(nM*9) (output)
    
    do m = 1,nM
      do ia = 1,3 !o,h,h
        hho_fa(:,(m-1)*3+ia) = xa_forces(:,ia,m)
      enddo
    enddo

    
    !// Add repulsion

!        a_oo = 8.23041822482709224512_dp !fitted to full rank interactions
!        b_oo = 3.48317298891623838841_dp

!        a_oo = 8.28581239473008496178_dp !fitted to reduced rank interactions
!        b_oo = 3.49709147777022710812_dp

    if(OO_REP)then
        
        if(FULL)then
          a_oo = 2872.579_dp!<full
          b_oo = 3.387391_dp!<full       !3.384538_dp
        else
          !stop"put in the fitting parameters for full rank OO repulsion"
          a_oo = 2872.579
          b_oo = 3.387391
        endif
        
        u_rep = 0
        do ox1 = 1,nM-1
           do ox2 = ox1+1,nM
               rr_ox = xyz_hho(:,3,ox1) - xyz_hho(:,3,ox2)
               r_ox = dsqrt( sum(rr_ox**2) )
               u_rep = u_rep + a_oo*dexp (- b_oo*r_ox)
           enddo
        enddo
        u_tot = u_tot + u_rep      !Repulsion energy
    elseif(ALL_REP)then
        
        if(FULL)then
          A_hh =  18.87599_dp !<full
          A_oh =  3658.516_dp !<full
          A_oo =  2342.881_dp !<full
          aa_hh =  3.71965_dp !<full
          aa_oh =  6.59638_dp !<full
          aa_oo = 3.348177_dp !<full
        else
          if(VAR_QUAD)then
            A_hh  = 3.245391_dp ! 
            A_oh  = 2873.137_dp ! 
            A_oo  =  2304.67_dp ! 
            aa_hh = 2.832741_dp ! 
            aa_oh = 6.104925_dp ! 
            aa_oo = 3.338482_dp ! 
          else
            A_hh =  13.87612_dp! 7.851853_dp!  !<redu
            A_oh =  2700.478_dp! 326.0703_dp!  !<redu
            A_oo =  2191.782_dp! 1381.591_dp!  !<redu
            aa_hh = 3.595182_dp!   2.8051_dp!  !<redu 
            aa_oh = 6.272737_dp! 4.854904_dp!  !<redu
            aa_oo = 3.322673_dp! 3.192549_dp!  !<redu
          endif
        endif
 













 
 
        
        !yy  = A_oo*dexp( -      aa_oo*d_oo) 
        !yy = yy + A_oh*dexp( -  aa_oh*d_oh1)
        !yy = yy + A_oh*dexp( -  aa_oh*d_oh2)
        !yy = yy + A_oh*dexp( -  aa_oh*d_h1o)
        !yy = yy + A_oh*dexp( -  aa_oh*d_h2o)
        !yy = yy + A_hh*dexp( -  aa_hh*d_h1h1)
        !yy = yy + A_hh*dexp( -  aa_hh*d_h1h2)
        !yy = yy + A_hh*dexp( -  aa_hh*d_h2h1)
        !yy = yy + A_hh*dexp( -  aa_hh*d_h2h2)        
        
        u_rep = 0
        do m1 = 1,nM-1
           do m2 = m1+1,nM
               rr_oo(:  ) = xyz_hho(:,3,m1) - xyz_hho(:,3,m2)
               rr_oh(:,1) = xyz_hho(:,3,m1) - xyz_hho(:,1,m2)
               rr_oh(:,2) = xyz_hho(:,3,m1) - xyz_hho(:,2,m2)
               rr_oh(:,3) = xyz_hho(:,1,m1) - xyz_hho(:,3,m2)
               rr_oh(:,4) = xyz_hho(:,2,m1) - xyz_hho(:,3,m2)
               rr_hh(:,1) = xyz_hho(:,1,m1) - xyz_hho(:,1,m2)
               rr_hh(:,2) = xyz_hho(:,1,m1) - xyz_hho(:,2,m2)
               rr_hh(:,3) = xyz_hho(:,2,m1) - xyz_hho(:,1,m2)
               rr_hh(:,4) = xyz_hho(:,2,m1) - xyz_hho(:,2,m2)
               
               r_oo = dsqrt(sum(rr_oo(:)**2))
               
               u_rep = u_rep + A_oo*dexp (- aa_oo*r_oo)
               
               do combi = 1,4
                   r_oh = dsqrt(sum(rr_oh(:,combi)**2))
                   r_hh = dsqrt(sum(rr_hh(:,combi)**2))
                   
                   u_rep = u_rep + A_oh*dexp (- aa_oh*r_oh)
                   u_rep = u_rep + A_hh*dexp (- aa_hh*r_hh)
                   
               enddo
               
           enddo
        enddo
        u_tot = u_tot + u_rep      !Repulsion energy
     endif
        
        
        
        
        
        
        
        
    
    
    !// Debug output ///////////////////////////////////////////////////!(pipe to file, diff to see change, comment to mute)
tprint(u_multipole,'u_multipole',s)
tprint(uDisp,'uDisp',s)
tprint(u_ps,'u_ps',s)
tprint(u_tot,'u_tot',s)
    
    
    !call h2o_to_linear(aforces,fa_test,nM)
    !call printer(fa_test,'aforces linear')
    
    !call printer(fa,'xa_forces linear')
!tprint(xa_forces,'xa_forces',s) 
!tprint(hho_fa,'hho_fa',s) 

    !s=0
!call printer(1,'1',s)    
    
    



    contains !/////////////////////////////////////////////////////////////////////////////////
     
    subroutine convert_d1234v_to_phi_scme
        do m = 1,nM !convert field to compressed form
            nn=1
            p1=pos_(nn)+1
            p2=pos_(nn+1)
            !print*,p1,p2
            phi_scme(p1:p2,m) = d1v(:,m) 
            
            nn=2
            p1=pos_(nn)+1
            p2=pos_(nn+1)
            !print*,p1,p2
            phi_scme(p1:p2,m) = compress(reshape(d2v(:,:,m),[3**nn]),nn)
            
            nn=3
            p1=pos_(nn)+1
            p2=pos_(nn+1)
            phi_scme(p1:p2,m) = compress(reshape(d3v(:,:,:,m),[3**nn]),nn)
            
            nn=4
            p1=pos_(nn)+1
            p2=pos_(nn+1)
            phi_scme(p1:p2,m) = compress(reshape(d4v(:,:,:,:,m),[3**nn]),nn)
            
            !polycompress_p(p11f,p12f,p22f,alp)
            call polycompress_p(dd(:,:,m),dq(:,:,:,m),qq(:,:,:,:,m),polz(:,:,m))
            
        enddo
    end
    
    subroutine convert_xpoles_to_qn_scme
    do m = 1,nM !convert multipoles to compressed form
        nn=1
        p1=pos_(nn)+1
        p2=pos_(nn+1)
        !print*,p1,p2
        qn_scme(p1:p2,m) = dpole(:,m) 
        
        nn=2
        p1=pos_(nn)+1
        p2=pos_(nn+1)
        !print*,p1,p2
        qn_scme(p1:p2,m) = compress(reshape(qpole(:,:,m),[3**nn]),nn)!1d0/3d0*
        
        nn=3
        p1=pos_(nn)+1
        p2=pos_(nn+1)
        qn_scme(p1:p2,m) = compress(reshape(opole(:,:,:,m),[3**nn]),nn)!1d0/15d0*
        
        nn=4
        p1=pos_(nn)+1
        p2=pos_(nn+1)
        qn_scme(p1:p2,m) = compress(reshape(hpole(:,:,:,:,m),[3**nn]),nn)!1d0/105d0*
        
        !polycompress_p(p11f,p12f,p22f,alp)
        call polycompress_p(dd(:,:,m),dq(:,:,:,m),qq(:,:,:,:,m),polz(:,:,m))
        
    enddo
    end

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
       
       
       
       
       
