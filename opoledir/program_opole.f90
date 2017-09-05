program jkdls
use printer_mod, only:printer
!use localAxes_mod,only: norm, norm_square
use data_types, only: dp, pi
use multipole_parameters, only: o0
use qpole,only:expansion_coordinates
use opole
implicit none

call main()

contains !/////////////////

subroutine main()
    integer, parameter :: xyz=3,hho=3
    
    real(dp) mass(hho), ang, rad,z_oh,x_h,r_h !,z_o z_h
    real(dp) pos_h1(xyz),pos_h2(xyz),pos_o(xyz)
    real(dp) molecule(xyz,hho)
    
    real(dp) exp_cent(xyz), cec(xyz,hho),cer2(hho)
    
    real(dp) z_h,z_o
    
    real(dp) OMxxz,OMyyz, q_h, q_o, charges(hho)
    
    real(dp) octupole(xyz,xyz,xyz)
    
    mass = [1d0,1d0,16d0]
    
    
    ang = 104.3475_dp
    rad = ang/180d0*pi
    r_h = 0.958649_dp
    
    z_oh = cos(rad*0.5d0)*r_h
    x_h = sin(rad*0.5d0)*r_h
    
    pos_h1 = [-x_h,0d0,0d0]
    pos_h2 = [ x_h,0d0,0d0]
    pos_o  = [0d0,0d0, z_oh ]
    
    
    
    molecule(:,1) = pos_h1
    molecule(:,2) = pos_h2
    molecule(:,3) = pos_o 

    
call printer(molecule,'molecule',2)    
    
    call expansion_coordinates(molecule,exp_cent,cec,cer2,1)
    

call printer(exp_cent,'exp_cent',2)    
call printer(cec,'cec',2)    
call printer(cer2,'cer2',2)    
    
    
    
    OMxxz = o0(1,1,3)
    OMyyz = o0(2,2,3)
    
    z_h = cec(3,1)
    z_o = cec(3,3)
    !r_h2 = cer2(2)
    !x_h as above
    
    
    !q_h = (OMyyz - OMxxz) / (5 * x_h**2 * z_h**2)
    !
    !q_o = 2 * (q_h * r_h**2 * z_h - OMyyz) / z_o**3     
    
    q_h = (OMxxz - OMyyz) / (5 * cec(1,2)**2 * cec(3,2))
    
    q_o = -2d0 * (q_h * cer2(2) * cec(3,2) + OMyyz) / cec(3,3)**3     
    
    
    print*, 'z_h:', z_h, 'z_o:', z_o
    print*, 'OMxxz:',OMxxz, 'OMyyz:',  OMyyz
    print*, 'q_o:', q_o, 'q_h:', q_h
    
    charges = [q_h,q_h,q_o]
    
    call octupole_tensor(cec,cer2,charges,octupole)
    
    call printer(octupole,'oct',2)
    
    call printer(o0,'o0',2)
    
    call printer(octupole-o0,'diff',2)
    
    octupole = 0
    call get_octupoles(cec,cer2,octupole,1)
    
    call printer(octupole,'oct',2)
    
    call printer(octupole-o0,'diff',2)    
end subroutine
end program
