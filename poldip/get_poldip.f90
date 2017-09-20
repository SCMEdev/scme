program jkdls
use printer_mod, only:printer
!use localAxes_mod,only: norm, norm_square
use data_types, only: dp, pi
use qpole,only:expansion_coordinates
use multipole_parameters, only: o0
use localAxes_mod,only: cross, norm_square
use opole
use qpole, only:matrix_trace
use polariz_parameters, only: dd0
implicit none

call main()

contains !/////////////////

subroutine main()
    integer, parameter :: xyz=3,hho=3,hholl=5
    
    real(dp) mass(hho), ang, rad,z_oh,x_h,r_h 
    real(dp) coords(xyz,hholl)
    real(dp) exp_cent(xyz), cec(xyz,hho),cer2(hho)
    real(dp),dimension(3) :: rol1, rol2, roh1, roh2, coords_o, coords_h1, coords_h2, coords_l1, coords_l2, yvec,zvec
    real(dp) square_dist5(hholl), cec5(xyz,hholl), yscale, zscale
    real(dp) alpha(xyz,xyz), iso(hholl)
    integer i
    real(dp) :: iso_alp(hho), isotropy, a(xyz),b(xyz), a2,b2,res_a2_b2, l1(xyz),l2(xyz)
    
    ! for octupole!
    real(dp) oct(xyz,xyz,xyz), charges(hholl)
    
    mass = [1d0,1d0,16d0]
    
    
    ang = 104.3475_dp
    rad = ang/180d0*pi
    r_h = 0.958649_dp
    
    z_oh = cos(rad*0.5d0)*r_h
    x_h = sin(rad*0.5d0)*r_h
    
    coords_h1 = [-x_h,0d0,0d0]
    coords_h2 = [ x_h,0d0,0d0]
    coords_o  = [0d0,0d0, z_oh ]
    
    roh1 = coords_h1 - coords_o
    roh2 = coords_h2 - coords_o
    
    !Lone-pairs 
    yscale = 0.5d0
    zscale = 0.2d0
    yvec = cross(roh1,roh2) * yscale
    zvec = -(roh1+roh2) * zscale
    rol1 = yvec + zvec 
    rol2 = -yvec + zvec 
    
    coords_l1 = rol1 + coords_o
    coords_l2 = rol2 + coords_o
    
  call printer(rol1,'rol1',2)
  call printer(rol2,'rol2',2)
    
    
    coords(:,1) = coords_h1
    coords(:,2) = coords_h2
    coords(:,3) = coords_o
    coords(:,4) = coords_l1
    coords(:,5) = coords_l2
    
    
    call expansion_coordinates(coords(:,1:hho),exp_cent,cec,cer2,nm=1,nsites=3)!nM,nsites
    !cec5(:,1:hho) = cec
    !cec5(:,4) = coords_l1 - exp_cent
    !cec5(:,5) = coords_l2 - exp_cent
    !square_dist5(1:hho) = cer2
    !square_dist5(4) = sum(cec5(:,4)**2)
    !square_dist5(5) = sum(cec5(:,5)**2)
    
    call expansion_coordinates(coords,exp_cent,cec5,square_dist5,nm=1,nsites=5)!nM,nsites
    
    !!! DIP-DIP POLARIZABILITY
    ! 5-SITE ! testing
    
call printer(coords,'coords',2)    
call printer(cec5,'cec5',2)    
call printer(square_dist5,'square_dist5',1)    
    
    iso(1) = 5d0
    iso(2) = 5d0
    iso(3) = 0d0
    iso(4) = 5d0
    iso(5) = 5d0
    
    alpha = 0
    isotropy = 0
    call get_alpha_from_isos(hholl,cec5,square_dist5,iso,isotropy, alpha)
    
call printer(alpha,'alpha traceless anisotropy',2)
    
    print*, 'trace', matrix_trace(alpha)
    do i = 1,xyz
      alpha(i,i) = alpha(i,i) + 9.6697586_dp !the isotropic dipole polarizability (assumed geometry independent here)
 
    enddo
    
call printer(alpha,'alpha',2)
    
    
    
    ! 3-SITE !
    ! 3-site (atomic) traceless anisotropy contribution
    call printer(dd0,'dd0',2)
    
    call get_3_iso_alp(cec,cer2,dd0, iso_alp,isotropy)
  call printer(iso_alp,'iso_alp',2)
  call printer(isotropy,'isotropy',2)  
  
    call get_alpha_from_isos(hho,cec,cer2,iso_alp,isotropy,alpha)
  call printer(alpha,'alpha for found charges + isotropy',2)
    
    ! OCTAPOLE TESTING!
    
    !New coordinates
    
    a = coords_h1 - coords_o
    b = coords_h2 - coords_o
    
    !Lone-pairs 
    yscale = 1d0
    zscale = 0.4d0
    
    a2 = norm_square(a)
    b2 = norm_square(b)
    res_a2_b2 = 1d0/(a2+b2)
    yvec = yscale * cross(a,b) * res_a2_b2
    zvec = - zscale * (roh1+roh2) * res_a2_b2
    l1 = yvec + zvec 
    l2 = -yvec + zvec 
    
    coords_l1 = l1 + coords_o
    coords_l2 = l2 + coords_o
    
  call printer(l1,'l1',2)
  call printer(l2,'l2',2)
    
    
    coords(:,1) = coords_h1
    coords(:,2) = coords_h2
    coords(:,3) = coords_o
    coords(:,4) = coords_l1
    coords(:,5) = coords_l2
    
    call expansion_coordinates(coords,exp_cent,cec5,square_dist5,nm=1,nsites=5)!nM,nsites
    
    !Assign charges
    charges(1) = 2d0
    charges(2) = 2d0
    charges(3) = 3d0
    charges(4) = 4d0
    charges(5) = 4d0
    
    call octupole_tensor(cec5,square_dist5,charges,oct)
    
  call printer(oct,'oct 5q',2)
    
    call octupole_tensor(cec,cer2,charges(1:3),oct)
    
  call printer(oct,'oct 3q',2)
    
    
    
    
end subroutine !//////////






subroutine get_alpha_from_isos(nsites,cec,cer2,iso_alp,isotropy,alpha)
    integer,parameter    :: xyz=3 !,hho=3
    integer,intent(in) :: nsites
    real(dp),intent(in)  :: cec(xyz,nsites),cer2(nsites), iso_alp(nsites),isotropy
    real(dp),intent(out) :: alpha(xyz,xyz)
    integer :: al,be,at,i
    
    alpha = 0
    do at = 1,nsites
      do al = 1,xyz !half matrix + symmetrization is better here
        do be = 1,xyz
        
          
          alpha(al,be) = alpha(al,be) + iso_alp(at) * 3*cec(al,at)*cec(be,at) 
          
          if (al==be) then 
            alpha(al,be) = alpha(al,be) - iso_alp(at)*cer2(at)
          endif
          
        enddo
      enddo
    enddo    
    alpha = alpha*0.5
    
call printer(alpha,'traceless alpha for found charges',2)
    
    do i = 1,xyz
      alpha(i,i) = alpha(i,i) + isotropy
    enddo
    
end subroutine !//////

SUBROUTINE get_3_iso_alp(cec,cer2,dd0, alp,isotropy)
    integer,parameter    :: xyz=3,hho=3, rra=3
    real(dp),intent(in)  :: cec(xyz,hho),cer2(hho),dd0(xyz,xyz)
    real(dp),intent(out) :: alp(hho),isotropy
    !internal
    real(dp) hx2, dd1(xyz,xyz), trace
    integer ii
    

    trace = matrix_trace(dd0)
    isotropy = trace/3d0
    
    dd1 = dd0 
    do ii = 1,xyz !detrace
      dd1(ii,ii) = dd1(ii,ii) - isotropy !1.43290917218_dp   
    enddo
    
    call printer(dd1,'traceless dd0',2)
    
    !compute charges
    hx2=cec(1,1)**2
    alp(1) = ( dd1(1,1) - dd1(2,2) ) / (3*hx2)
    alp(2)=alp(1)
    alp(3)= -2*( dd1(1,1) - alp(1)*(3*hx2-cer2(1)) )/cer2(3)
    
END SUBROUTINE 


end program


!    
!
!call printer(exp_cent,'exp_cent',2)    
!call printer(cec,'cec',2)    
!call printer(cer2,'cer2',2)    
!    
!    
!    
!    OMxxz = o0(1,1,3)
!    OMyyz = o0(2,2,3)
!    
!    z_h = cec(3,1)
!    z_o = cec(3,3)
!    !r_h2 = cer2(2)
!    !x_h as above
!    
!    
!    !q_h = (OMyyz - OMxxz) / (5 * x_h**2 * z_h**2)
!    !
!    !q_o = 2 * (q_h * r_h**2 * z_h - OMyyz) / z_o**3     
!    
!    q_h = (OMxxz - OMyyz) / (5 * cec(1,2)**2 * cec(3,2))
!    
!    q_o = -2d0 * (q_h * cer2(2) * cec(3,2) + OMyyz) / cec(3,3)**3     
!    
!    
!    print*, 'z_h:', z_h, 'z_o:', z_o
!    print*, 'OMxxz:',OMxxz, 'OMyyz:',  OMyyz
!    print*, 'q_o:', q_o, 'q_h:', q_h
!    
!    charges = [q_h,q_h,q_o]
!    
!    call octupole_tensor(cec,cer2,charges,octupole)
!    
!    call printer(octupole,'oct',2)
!    
!    call printer(o0,'o0',2)
!    
!    call printer(octupole-o0,'diff',2)
!    
!    octupole = 0
!    call get_octupoles(cec,cer2,octupole,1)
!    
!    call printer(octupole,'oct',2)
!    
!    call printer(octupole-o0,'diff',2)    
!    
    

!    do al = 1,xyz
!      do be = 1,xyz
!        do at = 1,hholl
!          alpha(al,be) = alpha(al,be) + iso(at)* cec5(al,at) * cec5(be,at) 
!          !!! Detrace !not sure if this is a good Idea. Maybe just use the tracefull. 
!          ! The tracefull needs lone-pairs to get a yy-element, not the traceless. 
!          ! But the physicallity and the sensibility in how the polarizabilities transform w.r.t. 
!          !  the internal coordinates might suffer from using the trace-less trick with no lone-pairs. 
!          !  That is true also for the multipoles above dipole!
!          ! The question is fi they should all have the same lone-pair sites... SHould work, since we can vary the charges. 
!          !
!          ! Both the isotropic and anisotropic contributions need, however, be geonetry dipendent, at some point. 
!          ! with no tracelessness the isotropic part is also allowed to vary, and it would increase
!          ! with increasing molecular size, which could be sensible...?
!          !
!          ! The symmetric stretch mode probably doesn't alter the anisotropy much, but rather the isotropic contribution. 
!          ! Meanwhile the bend-mode probably alters the anisotropic part the most. 
!          ! maybe a sensible approximation is to make the anisotropic part depend only on h1h2, and the isotropic only on oh1+oh2. 
!          ! Anyway, with isotropic sites formalism, this separation cant be made. 
!          !    How to make the site-polarizabilities vary with geometry might however be affected in this way. 
!          !    If we have an isotropic part then this one can be altered only with h1+h2
!          !    while the site-plzties vary only with the bend
!          ! The alternative is to have 4 or 5 site polarizabilities and then make them vary with the geometry to 
!          !    to alter both the isotropic and aniso part. Autmatically giving a larger polarizability with a larger
!          !    molecule may be a natural benefit. 
!          !
!          ! for the lone pais site we may want to have a "reciprocal cross product rule" 
!          !  so the sites get closer not only when the moecule gets straight, but also
!          !  when the hydrogens get further appart, since O and HO dont have lone-pairs. 
!          !  this somehow resembles electronic rearrangement / state change / internal conversion
!          if (al==be) then 
!            alpha(al,be) = alpha(al,be) - iso(at)*square_dist5(at)/3d0
!          endif
!          
!        enddo
!      enddo
!    enddo
