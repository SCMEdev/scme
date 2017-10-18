MODULE detrace_mod
use printer_mod, only: printer,printo, str
use data_types, only:dp
implicit none

!interface detrace
!    module procedure detrace_mono, detrace_dip, detrace_quad, detrace_octa, detrace_hexadeca
!end interface


CONTAINS

subroutine test()
    real(dp) gauss_quad(3), entries(10), entries15(15)
    character(3) char3
    character(4) char4
    integer i
    
    
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
    
    open(55,file="tracefull_octa.txt")
    do i = 1,10
        read(55,*) char3, entries(i)  !
    enddo
    close(55)
    
    
    print*, "octa:"
    call detrace_octa(entries,.true.)                                                                      !          < !!!!!!!!!!!!!!!!!!!!!!!!!!
    !=  3.0000_dp
    != -2.8900_dp
    != -1.1002_dp
    != -7.1200_dp
    != -4.3400_dp
    != -4.0340_dp
    !=  6.3200_dp
    !=  1.2205_dp
    != -2.9800_dp
    !=  7.0370_dp
    
    
    
    
    
    open(56,file="tracefull_hexa.txt")
    do i = 1,15
        read(56,*) char4, entries15(i)  
    enddo
    close(56)
    
    
    call detrace_hexadeca(entries15,.true.)
    
    call bat_hexa()
    
    

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

subroutine detrace(entries)
    real(dp), intent(inout) :: entries(:)
    integer le
    le = size(entries)
    if(le==15)then
      call detrace_hexadeca(entries,.false.)
    elseif(le==10)then
      call detrace_octa(entries,.false.)
    elseif(le==6)then
      call detrace_quad(entries)
    elseif(le==3)then
      print*, "dipole unchanged"
    elseif(le==1)then
    else
      stop"only 15,10,6,3,1 entries"
    endif
    
end subroutine

!subroutine detrace_dip(entries)
!    real(dp), intent(inout) :: entries(3)
!end subroutine
!
!subroutine detrace_mono(entries)
!    real(dp), intent(inout) :: entries
!end subroutine

subroutine detrace_quad(entries)
    implicit real(dp) (x,y,z)
!    logical, intent(in) :: do_print
    real(dp), intent(inout) :: entries(6)
    real(dp) aa
    
    xx = entries(1) 
    !xy = entries(2) 
    !xz = entries(3) 
    yy = entries(4) 
    !yz = entries(5) 
    zz = entries(6) 
    
    aa = (xx+yy+zz)/3
    
    xx_det = xx - aa
    yy_det = yy - aa
    zz_det = zz - aa
    
    entries(1) = xx_det
    entries(4) = yy_det
    entries(6) = zz_det
    
end subroutine


    
    
subroutine detrace_octa(entries,do_print)
    implicit real(dp) (x,y,z)
    logical, intent(in) :: do_print
    real(dp), intent(inout) :: entries(10)
    real(dp) aa_x,aa_y,aa_z,trace_x,trace_y,trace_z
    character(10) :: txt10(10)
    integer i
    
    ! Applequist tricorn ordering
    xxx = entries(1); xzz = entries(6)   
    xxy = entries(2); yyy = entries(7)   
    xxz = entries(3); yyz = entries(8)   
    xyy = entries(4); yzz = entries(9)   
    xyz = entries(5); zzz = entries(10)  
    
    
    trace_z = xxz+yyz+zzz
    trace_y = xxy+yyy+yzz
    trace_x = xxx+xyy+xzz
    
    aa_x = trace_x/5_dp !aa here just means trace/5
    aa_y = trace_y/5_dp
    aa_z = trace_z/5_dp 
    
    xxx_det = xxx - aa_x*3
    xyy_det = xyy - aa_x
    xzz_det = xzz - aa_x
    
    xxy_det = xxy - aa_y
    yyy_det = yyy - aa_y*3
    yzz_det = yzz - aa_y
    
    xxz_det = xxz - aa_z
    yyz_det = yyz - aa_z
    zzz_det = zzz - aa_z*3
    
    xyz_det = xyz !no trace condition 
    
    entries(1)  = xxx_det
    entries(2)  = xxy_det
    entries(3)  = xxz_det
    entries(4)  = xyy_det
    entries(5)  = xyz_det
    entries(6)  = xzz_det
    entries(7)  = yyy_det
    entries(8)  = yyz_det
    entries(9)  = yzz_det
    entries(10) = zzz_det
    
    !entries = entries*5d0/2d0 !5/2 in Stone definition, not in Gaussian Program
    
    if(do_print)then
        
        txt10(1)  = 'xxx_det'
        txt10(2)  = 'xxy_det'
        txt10(3)  = 'xxz_det'
        txt10(4)  = 'xyy_det'
        txt10(5)  = 'xyz_det'
        txt10(6)  = 'xzz_det'
        txt10(7)  = 'yyy_det'
        txt10(8)  = 'yyz_det'
        txt10(9)  = 'yzz_det'
        txt10(10) = 'zzz_det'
        
        do i = 1,10
          print'(f20.15,a,a)',entries(i),'   ',txt10(i)
        enddo      
        
        print*, 'octa trace', xxz_det+yyz_det+zzz_det,'prev:',trace_z
        print*, 'octa trace', xxy_det+yyy_det+yzz_det,'prev:',trace_y
        print*, 'octa trace', xxx_det+xyy_det+xzz_det,'prev:',trace_x
        
        print*, 'C_2v zzz', zzz_det*5d0/2d0
        print*, 'C_2v xxz', xxz_det*5d0/2d0
        print*, 'C_2v yyz', yyz_det*5d0/2d0
        
    endif
    
end subroutine    


subroutine detrace_hexadeca(entries, do_print)
    implicit real(dp) (x,y,z,t)
    real(dp) :: aa,bb,cc,aa_xy,aa_xz,aa_yz
    real(dp) :: trace6(6), out6(6)
    real(dp), intent(inout) :: entries(15)
    character(10) :: txt15(15)
    integer i
    logical do_print
    
    ! Applequist tricorn ordering
    xxxx = entries(1); xxzz = entries(6) ; yyyy = entries(11)
    xxxy = entries(2); xyyy = entries(7) ; yyyz = entries(12)
    xxxz = entries(3); xyyz = entries(8) ; yyzz = entries(13)
    xxyy = entries(4); xyzz = entries(9) ; yzzz = entries(14)
    xxyz = entries(5); xzzz = entries(10); zzzz = entries(15)
    
    !Traces
    ! xx: xx,yy,zz
    ! xz: xx,yy,zz
    ! yz: zz,xx,yy
    ! yy: yy,xx,zz
    ! yx: yy,zz,xx
    ! zz: zz,xx,yy
    
    ! Traces AACC, AA = xx,yy,zz //////////////////////////////////////////////////////
    trace_xx = xxxx+xxyy+xxzz
    trace_yy = xxyy+yyyy+yyzz
    trace_zz = xxzz+yyzz+zzzz
    
    aa = trace_xx/35_dp
    bb = trace_yy/35_dp
    cc = trace_zz/35_dp
    
    xxxx_det = xxxx + 3*(-9*aa + bb + cc)
    yyyy_det = yyyy + 3*(aa - 9*bb + cc)
    zzzz_det = zzzz + 3*(aa + bb - 9*cc)
    
    xxyy_det = xxyy - 4*aa - 4*bb + cc
    xxzz_det = xxzz - 4*aa + bb - 4*cc
    yyzz_det = yyzz + aa - 4*bb - 4*cc
    
    ! Traces ABCC, AB = xy,xz,yz //////////////////////////////////////////////////////////
    trace_xy = xxxy+xyyy+xyzz
    trace_xz = xxxz+xyyz+xzzz
    trace_yz = xxyz+yyyz+yzzz
    
    aa_xy = trace_xy/7_dp !a is just trace/7
    aa_xz = trace_xz/7_dp
    aa_yz = trace_yz/7_dp
    
    xxxy_det = xxxy - aa_xy*3
    xyyy_det = xyyy - aa_xy*3
    xyzz_det = xyzz - aa_xy
    
    xxxz_det = xxxz - aa_xz*3
    xyyz_det = xyyz - aa_xz
    xzzz_det = xzzz - aa_xz*3
    
    xxyz_det = xxyz - aa_yz
    yyyz_det = yyyz - aa_yz*3
    yzzz_det = yzzz - aa_yz*3
    
    
    entries(1)  = xxxx_det
    entries(2)  = xxxy_det
    entries(3)  = xxxz_det
    entries(4)  = xxyy_det
    entries(5)  = xxyz_det
    entries(6)  = xxzz_det
    entries(7)  = xyyy_det
    entries(8)  = xyyz_det
    entries(9)  = xyzz_det
    entries(10) = xzzz_det
    entries(11) = yyyy_det
    entries(12) = yyyz_det
    entries(13) = yyzz_det
    entries(14) = yzzz_det
    entries(15) = zzzz_det
    !entries = entries*9_dp/2_dp !definition in stone and gaussian differ. This is also wrong, the fraction (2n-1)!!/n! = 7*5*3/(4*3*2) = 35/8
    
    
    if(do_print)then
        print*, "Hexadecapole routine printing: (*9_dp/2_dp)"
        
        
        txt15(1)  = 'xxxx_det'
        txt15(2)  = 'xxxz_det'
        txt15(3)  = 'yzzz_det'
        txt15(4)  = 'xxyz_det'
        txt15(5)  = 'yyyy_det'
        txt15(6)  = 'xyyy_det'
        txt15(7)  = 'xxyy_det'
        txt15(8)  = 'xyyz_det'
        txt15(9)  = 'zzzz_det'
        txt15(10) = 'yyyz_det'
        txt15(11) = 'xxzz_det'
        txt15(12) = 'xyzz_det'
        txt15(13) = 'xxxy_det'
        txt15(14) = 'xzzz_det'
        txt15(15) = 'yyzz_det'
        
        do i = 1,15
          print'(f20.15,a,a)',entries(i)*9_dp/2_dp,'   ',txt15(i)
        enddo
        
        trace6(1) = yyyy_det+xxyy_det+yyzz_det 
        trace6(2) = zzzz_det+xxzz_det+yyzz_det 
        trace6(3) = xxzz_det+xxxx_det+xxyy_det 
        trace6(4) = xxxy_det+xyyy_det+xyzz_det
        trace6(5) = xxxz_det+xyyz_det+xzzz_det
        trace6(6) = xxyz_det+yyyz_det+yzzz_det
        
        do i = 1,6
          print'(a,e15.5)', 'trace:',trace6(i)*9_dp/2_dp
        enddo
        
        print*, "batista comparison:"
        out6(1) = xxxx_det
        out6(2) = yyyy_det
        out6(3) = zzzz_det
        out6(4) = yyzz_det
        out6(5) = xxzz_det
        out6(6) = xxyy_det
        out6 = out6*9_dp/2_dp
        
        call printer(out6*9_dp/2_dp,'out6',1)        
       
        do i = 1,3
           print*, str((out6(i) - out6(i+3))*9_dp/2_dp), ' gauss diff'
        enddo
        
    endif
    
    !! HHmmmm seems like there are only three unique components in reality!?
end subroutine

ENDMODULE

