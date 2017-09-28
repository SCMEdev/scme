program jkdls
use printer_mod, only:printer, str, printo
!use localAxes_mod,only: norm, norm_square
use data_types, only: dp, pi, au_to_debye, a0_A, A_a0, eA_to_debye
use qpole,only:expansion_coordinates
use multipole_parameters, only: o0, q0,d0
use localAxes_mod,only: cross, norm_square
use opole
use qpole, only:matrix_trace
use polariz_parameters, only: dd0, dq0, qq0
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
    real(dp) oct(xyz,xyz,xyz), charges(hholl), levi(xyz,xyz,xyz), three(xyz,xyz,xyz), entries(xyz)
    !real(dp) ddelta(3,3,3,3)
    
    real(dp) hede(3,3,3,3), tempmat(3,3)
    real(dp) gauss_quad(3),trace, rtemp !,quad(3,3)
    integer itemp
    character(5) ctemp
    
    real(dp) out6(6) ,output(15)
    
    mass = [1d0,1d0,16d0]
    
    
    ang = 104.28!<- Batista! 104.3475_dp!<-PS!
    rad = ang/180d0*pi
    r_h = 0.9590!<-Batista! 0.958649_dp!<-PS!
    
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
  call printer(cec,'cec',2)
    !cec5(:,1:hho) = cec
    !cec5(:,4) = coords_l1 - exp_cent
    !cec5(:,5) = coords_l2 - exp_cent
    !square_dist5(1:hho) = cer2
    !square_dist5(4) = sum(cec5(:,4)**2)
    !square_dist5(5) = sum(cec5(:,5)**2)
    
    call expansion_coordinates(coords,exp_cent,cec5,square_dist5,nm=1,nsites=5)!nM,nsites
    
    
    !!! DIP-DIP POLARIZABILITY ///////////////////////////////////////////////////
if(.false.)then 
    ! 5-SITE ! testing
    
  !call printer(coords,'coords',2)    
  !call printer(cec5,'cec5',2)    
  !call printer(square_dist5,'square_dist5',1)    
    
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
    
  call printer(oct,'oct 5q',1)
    
    call octupole_tensor(cec,cer2,charges(1:3),oct)
    
  call printer(oct,'oct 3q',1)
  
  call printer(dq0*eA_to_debye,'data dq0',1)
  call printer(q0*eA_to_debye,'data q0',1)
  call printer(d0*eA_to_debye,'data d0',1)
    
    print*,au_to_debye, eA_to_debye, A_a0
    
  call printer(o0*eA_to_debye,'data o0',1)
    call levicivita(levi)
  call printer(levi,'levi',1)
  
    call octa_trace_entries(three)
  call printer(three,'three',1)

!    call octa_iso_entries(three)
!  call printer(three,'three',1)
    
    
endif
    !Traceless!
    gauss_quad(1) =  1.7423_dp!1.7407 !<-PS-geom
    gauss_quad(2) = -1.6539_dp!-1.6508!<-PS-geom
    gauss_quad(3) = -0.0884_dp!-0.0899!<-PS-geom
    gauss_quad = gauss_quad*3/2!/au_to_debye*A_a0
  call printer (gauss_quad,'gauss_quad traceless mp2',1)
    
    
    !tracefull!!
    !gauss_quad(1) = -4.4517
    !gauss_quad(2) = -7.8432
    !gauss_quad(3) = -6.2824
    !trace = sum(gauss_quad)/3d0
    !print*, 'trace', trace
    !gauss_quad = gauss_quad - trace
    !
    !gauss_quad = gauss_quad*3/2/au_to_debye*A_a0
    !
    !call printer (gauss_quad,'gauss_quad detraced tracefull mp2',1)

    call gaussian_octa(entries)
  call printer(entries,'gauss oct mp2',1) !Why dont I need the 5/2 scale factor here???? because it is in the function:) 
    
    call detrace_hexadecap(out6,output,.true.)
    
    call bat_hexa()
    
    
!    call get_deldel(ddelta)
!  call printer(ddelta,'ddelta',1)
    
    !call hexad_moments(hede)
    
if(.false.)then
    !call printer(dd0,'dd0',1)
    !call printer(dq0,'dq0',1)
    !call printer(dq0(:,1,:),'dq1',1)
    !call printer(dq0(:,2,:),'dq2',1)
    !call printer(dq0(:,3,:),'dq3',1)
    
    !do i = 1,3
    !  !write(ctemp,'(I1)') i
    !  call printer(dq0(i,:,:),'quadrupole induced from '//str(i),1)
    !enddo
    !call printer(qq0,'qq0',1)
    
    
    !call print_row([1.2d0,1.3d0,1.4d0])
    
    tempmat=dd0
    tempmat(1,3) = 1d0
    tempmat(1,2) = 2d0
    
    
    call print_tens(tempmat,2)
    
endif 
    
    !call printo(qq0,[2,4,1,3])
    !print*, "hej"
    !call printo(qq0(:,:,:,1),[2,3,1])
    !print*, "hej"
    !
    !call printo(qq0(3,:,:,1),[1,2])
    !          
    !call printo([1.0d0/3,2.0d0/3,3.0d0/4])
    !
    !call detrace_hede(out6,output,.true.)
    
    !call printer(output,'printer(output)',1)
    
    !call printer(out6,'printer(out6)',1)
    !call printer(qq0,'qq0',1)
    !call print_tens(dd0,2)
end subroutine !//////////

subroutine print_row(vec)
real(dp),intent(in) :: vec(3) 
write(*,'(3f12.5)') vec
endsubroutine

subroutine print_tens(obj,order)
real(dp) :: obj(3,3)
integer order
integer i
do i = 1,3
  if(order==1)then 
    call print_row(obj(i,:))
  elseif(order==2)then
    call print_row(obj(:,i))
  else
    print*, "Not valid order"
  endif
enddo
endsubroutine




!subroutine print_tens3(obj,order)
!real(dp) :: obj(3,3,3)
!integer order(3)
!integer i, io
!do j = 1,3
!  print*, 'outermost'
!  do i = 1,3 !innermost
!    do io = 1,2
!      if(order==[1,2,3])then 
!        call print_row(obj(i,:,j))
!      elseif(order==[2,1,3])then
!        call print_row(obj(:,i,j))
!      elseif(order==[1,2,3])then 
!        call print_row(obj(i,j,:))
!      else
!        print*, "Not valid order"
!      endif
!  enddo
!enddo
!endsubroutine

!function i2s(integ,int_len,left) result(ch)
!integer,intent(in) :: integ, int_len
!integer, intent(in), optional :: left
!character(int_len) ch
!character(5) str_len
!write(str_len,'(I5)') int_len
!write(ch,'(I'//trim(str_len)//')') integ
!if(present(left))then
!  if (left==1) ch = adjustl(ch)
!endif
!endfunction



subroutine levicivita(levi)
    integer,parameter    :: xyz=3
    real(dp), intent(out) :: levi(xyz,xyz,xyz)
    
    levi = 0
    levi(1,2,3) = 1d0
    levi(3,1,2) = 1d0
    levi(2,3,1) = 1d0
    
    levi(1,3,2) = -1d0
    levi(2,1,3) = -1d0
    levi(3,2,1) = -1d0
    
end subroutine
    
subroutine octa_iso_entries(three)
    integer,parameter    :: xyz=3
    real(dp), intent(out) :: three(xyz,xyz,xyz)
    integer al,be
    real(dp) ent
    
    three = 0
    
    do al = 1,xyz
      do be = 1,xyz
        if (be==1) ent=2
        if (be==2) ent=3
        if (be==3) ent=-5
        three(al,be,be) = ent
        if (al /= be) then
          three(be,al,be) = ent
          three(be,be,al) = ent
        endif
      enddo
    enddo
    
end subroutine

subroutine octa_trace_entries(three)
    integer,parameter    :: xyz=3
    real(dp), intent(out) :: three(xyz,xyz,xyz)
    integer al,be,ga
    
    three = 0
    
    do al = 1,xyz
      do be = 1,xyz
        do ga = 1,xyz
          three(al,be,ga) = 10d0/9d0*delta(al,be) + 1d0/9d0*delta(al,ga) + 0.1d0/9d0*delta(ga,be)
        enddo
      enddo
      
    enddo
    
end subroutine

subroutine gaussian_octa(entries)
    integer,parameter    :: xyz=3
    real(dp), intent(out) :: entries(xyz)
    real(dp) a_trace
    !XXX=              0.0000  YYY=              0.0000  ZZZ=              0.0370  XYY=             -0.0000
    !XXY=             -0.0000  XXZ=             -1.1002  XZZ=             -0.0000  YZZ=             -0.0000
    !YYZ=              0.2205  XYZ=             -0.0000

    entries(1) =  0.0370_dp!<-batistaGeom GussMP2!  0.0463!zzz
    entries(2) = -1.1002_dp!<-batistaGeom GussMP2! -1.0939!xxz
    entries(3) =  0.2205_dp!<-batistaGeom GussMP2!  0.2214!yyz
    a_trace = sum(entries)/5d0
    
    entries(1) = entries(1) - a_trace*3
    entries(2) = entries(2) - a_trace
    entries(3) = entries(3) - a_trace
    
    entries = entries*5d0/2d0
    
    
end subroutine    

subroutine get_deldel(ddelta)
    integer,parameter    :: xyz=3 !,hho=3
    real(dp),intent(out) :: ddelta(3,3,3,3)
    integer al,be,ga,de
    do al = 1,3
      do be = 1,3
        do ga = 1,3
          do de = 1,3
            ddelta(al,be,ga,de) = 10d0/9d0*delta(al,be)*delta(ga,de) &
                                + 1d0/9d0*delta(al,de)*delta(ga,be) &
                                + 0.10d0/9d0*delta(al,ga)*delta(de,be)
          enddo
        enddo
      enddo
    enddo
end subroutine    
    


function delta(al,be) result(res)
    integer,parameter    :: xyz=3 !,hho=3
    integer,intent(in) :: al,be
    real(dp) :: res
    
    if(al==be)then
      res = 1d0
    else 
      res = 0d0
    endif
    
end function


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

subroutine detrace_hexadecap(out6, output15, do_print)
    implicit real(dp) (p-z)
    real(dp) :: trace3(3)
    real(dp), intent(out) :: output15(15), out6(6)
    character(10) :: txt15(15)
    integer i
    logical do_print
    
    
    xxxx=             -6.4684_dp; yyyy=             -7.6549_dp; zzzz=             -7.8338_dp; xxxy=             -0.0000_dp;
    xxxz=             -0.0000_dp; yyyx=             -0.0000_dp; yyyz=              0.0000_dp; zzzx=              0.0000_dp;
    zzzy=              0.0000_dp; xxyy=             -2.6905_dp; xxzz=             -1.9449_dp; yyzz=             -2.6498_dp;
    xxyz=              0.0000_dp; yyxz=              0.0000_dp; zzxy=              0.0000_dp;
    
    !Traces
    ! xx: xx,yy,zz
    ! xz: xx,yy,zz
    ! yz: zz,xx,yy
    ! yy: yy,xx,zz
    ! yx: yy,zz,xx
    ! zz: zz,xx,yy
    
    ! Traces AACC, AA = xx,yy,zz //////////////////////////////////////////////////////
    sum_xx = xxxx+xxyy+xxzz
    sum_yy = xxyy+yyyy+yyzz
    sum_zz = xxzz+yyzz+zzzz
    
    ra = sum_xx/35d0
    rb = sum_yy/35d0
    rc = sum_zz/35d0
    
    xxxx_out = xxxx + 3*(-9*ra + rb + rc)
    yyyy_out = yyyy + 3*(ra - 9*rb + rc)
    zzzz_out = zzzz + 3*(ra + rb - 9*rc)
    
    xxyy_out = xxyy - 4*ra - 4*rb + rc
    xxzz_out = xxzz - 4*ra + rb - 4*rc
    yyzz_out = yyzz + ra - 4*rb - 4*rc
    
    ! Traces ABCC, AB = xy,xz,yz //////////////////////////////////////////////////////////
    sum_xy = xxxy+yyyx+zzxy
    sum_xz = xxxz+yyxz+zzzx
    sum_yz = xxyz+yyyz+zzzy
    
    ra_xy = sum_xy/7_dp
    ra_xz = sum_xz/7_dp
    ra_yz = sum_yz/7_dp
    
    xxxy_out = xxxy - ra_xy*3
    yyyx_out = yyyx - ra_xy*3
    zzxy_out = zzxy - ra_xy
    
    xxxz_out = xxxz - ra_xz*3
    yyxz_out = yyxz - ra_xz
    zzzx_out = zzzx - ra_xz*3
    
    xxyz_out = xxyz - ra_yz
    yyyz_out = yyyz - ra_yz*3
    zzzy_out = zzzy - ra_yz*3
    
    
    output15(1)  = xxxx_out
    output15(2)  = xxxz_out
    output15(3)  = zzzy_out
    output15(4)  = xxyz_out
    output15(5)  = yyyy_out
    output15(6)  = yyyx_out
    output15(7)  = xxyy_out
    output15(8)  = yyxz_out
    output15(9)  = zzzz_out
    output15(10) = yyyz_out
    output15(11) = xxzz_out
    output15(12) = zzxy_out
    output15(13) = xxxy_out
    output15(14) = zzzx_out
    output15(15) = yyzz_out
    output15 = output15*9_dp/2_dp
    
    
    if(do_print)then
        
        out6(1) = xxxx_out
        out6(2) = yyyy_out
        out6(3) = zzzz_out
        out6(4) = yyzz_out
        out6(5) = xxzz_out
        out6(6) = xxyy_out
        out6 = out6*9_dp/2_dp
        call printer(out6,'out6',1)
        
        txt15(1)  = 'xxxx_out'
        txt15(2)  = 'xxxz_out'
        txt15(3)  = 'zzzy_out'
        txt15(4)  = 'xxyz_out'
        txt15(5)  = 'yyyy_out'
        txt15(6)  = 'yyyx_out'
        txt15(7)  = 'xxyy_out'
        txt15(8)  = 'yyxz_out'
        txt15(9)  = 'zzzz_out'
        txt15(10) = 'yyyz_out'
        txt15(11) = 'xxzz_out'
        txt15(12) = 'zzxy_out'
        txt15(13) = 'xxxy_out'
        txt15(14) = 'zzzx_out'
        txt15(15) = 'yyzz_out'
        
        do i = 1,15
          print'(f20.15,a,a)',output15(i),'   ',txt15(i)
        enddo
        
        trace3(1) = yyyy_out+xxyy_out+yyzz_out 
        trace3(2) = zzzz_out+xxzz_out+yyzz_out 
        trace3(3) = xxzz_out+xxxx_out+xxyy_out 
        trace3(4) = xxxy_out+yyyx_out+zzxy_out
        trace3(5) = xxxz_out+yyxz_out+zzzx_out
        trace3(6) = xxyz_out+yyyz_out+zzzy_out
        
        do i = 1,6
          print'(a,e15.5)', 'trace:',trace3(i)
        enddo
        
        do i = 1,3
           print*, str(out6(i) - out6(i+3)), ' gauss diff'
        enddo
        
    endif
    
    !! HHmmmm seems like there are only three unique components in reality!?
end subroutine

subroutine bat_hexa()
    implicit real(dp) (x,y,z)
    real(dp) bat6(6)
    zzzz_bat = -1.3637_dp
    xxzz_bat =  1.6324_dp
    yyzz_bat = -0.2687_dp
    xxxx_bat = -0.3575_dp
    xxyy_bat = -1.2749_dp
    yyyy_bat =  1.5436_dp
    
    bat6(1) = xxxx_bat
    bat6(2) = yyyy_bat
    bat6(3) = zzzz_bat
    bat6(4) = yyzz_bat
    bat6(5) = xxzz_bat
    bat6(6) = xxyy_bat
    
    call printer(bat6,'bat6',1)
    
    print*, str(zzzz_bat - xxyy_bat), ' batista diff'
    print*, str(yyyy_bat - xxzz_bat), ' batista diff'
    print*, str(xxxx_bat - yyzz_bat), ' batista diff'
    
end subroutine

subroutine hexad_moments(entries)
    implicit real(dp) (r-z)
    real(dp) entries(6)
    xxxx = -6.5402
    yyyy = -7.7054
    zzzz = -7.8992
    xxyy = -2.7107
    xxzz = -1.9781
    yyzz = -2.6673
    
    entries(1) = xxxx
    entries(2) = yyyy
    entries(3) = zzzz
    entries(4) = xxyy
    entries(5) = xxzz
    entries(6) = yyzz
    
    print*, 'xxxx+xxyy+xxzz',xxxx+xxyy+xxzz
    print*, 'yyyy+xxyy+yyzz',yyyy+xxyy+yyzz
    print*, 'zzzz+yyzz+xxzz',zzzz+yyzz+xxzz
    
    !rescaling while printing
    call printer(entries*7/2,'gaussian hede tracefull ',1)
    
    trace = (xxxx+xxyy+xxzz)/3d0 
    !traceless rescaled
    entries = (entries - trace)*7/2
    
    call printer(entries,'gaussian hede detraced',1)
    
    print*, 'xxxx+xxyy+xxzz',entries(1)+entries(4)+entries(5)
    print*, 'yyyy+xxyy+yyzz',entries(2)+entries(4)+entries(6)
    print*, 'zzzz+yyzz+xxzz',entries(3)+entries(5)+entries(6)
    
    
    
end subroutine


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
