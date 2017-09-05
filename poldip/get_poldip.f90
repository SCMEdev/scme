program jkdls
use printer_mod, only:printer
!use localAxes_mod,only: norm, norm_square
use data_types, only: dp, pi
use multipole_parameters, only: o0
use qpole,only:expansion_coordinates
use localAxes_mod,only: cross
use opole
use qpole, only:matrix_trace
implicit none

call main()

contains !/////////////////

subroutine main()
    integer, parameter :: xyz=3,hho=3,hholl=5
    
    real(dp) mass(hho), ang, rad,z_oh,x_h,r_h !,z_o z_h
!    real(dp) pos_h1(xyz),pos_h2(xyz),pos_o(xyz)
    real(dp) coords(xyz,hholl)
    
    real(dp) exp_cent(xyz), cec(xyz,hho),cer2(hho)
    
    real(dp) z_h,z_o
    
    real(dp) OMxxz,OMyyz, q_h, q_o, charges(hho)
    
    real(dp) octupole(xyz,xyz,xyz)
    
    real(dp),dimension(3) :: rol1, rol2, roh1, roh2, coords_o, coords_h1, coords_h2, coords_l1, coords_l2, yvec,zvec
    
    real(dp) square_dist5(hholl), cec5(xyz,hholl), yscale, zscale
    
    real(dp) iso_h1, iso_h2, iso_l1, iso_l2
    
    real(dp) alpha(xyz,xyz), iso(hholl)
    
    integer al, be, at, i
    
    
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
    
    
    call expansion_coordinates(coords(:,1:hho),exp_cent,cec,cer2,1)
    cec5(:,1:hho) = cec
    cec5(:,4) = coords_l1 - exp_cent
    cec5(:,5) = coords_l2 - exp_cent
    square_dist5(1:hho) = cer2
    square_dist5(4) = sum(cec5(:,4)**2)
    square_dist5(5) = sum(cec5(:,5)**2)
    
    
    
    
call printer(coords,'coords',2)    
call printer(cec5,'cec5',2)    
call printer(square_dist5,'square_dist5',1)    
    
    iso(1) = 5d0
    iso(2) = 5d0
    iso(3) = 0d0
    iso(4) = 5d0
    iso(5) = 5d0
    
    alpha=0
    do al = 1,xyz
      do be = 1,xyz
        do at = 1,hholl
          alpha(al,be) = alpha(al,be) + iso(at)* cec5(al,at) * cec5(be,at) 
          !!! Detrace !not sure if this is a good Idea. Maybe just use the tracefull. 
          ! The tracefull needs lone-pairs to get a yy-element, not the traceless. 
          ! But the physicallity and the sensibility in how the polarizabilities transform w.r.t. 
          !  the internal coordinates might suffer from using the trace-less trick with no lone-pairs. 
          !  That is true also for the multipoles above dipole!
          ! The question is fi they should all have the same lone-pair sites... SHould work, since we can vary the charges. 
          if (al==be) then 
            alpha(al,be) = alpha(al,be) - iso(at)*square_dist5(at)/3d0
          endif
          
        enddo
      enddo
    enddo
    
    !detrace
    !alpha(3,3) = -(alpha(1,1)+alpha(2,2))
call printer(alpha,'alpha traceless anisotropy',2)
    print*, 'trace', matrix_trace(alpha)
    do i = 1,xyz
      alpha(i,i) = alpha(i,i) + 9.5d0 
    enddo
    
call printer(alpha,'alpha',2)
    
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
    
end subroutine


end program
