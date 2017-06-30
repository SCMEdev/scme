program test_g
use localAxes_mod
real*8 r, gg(5)
r=0
do while (r<10.0)
 r=r+0.01
 call gs(r,gg)
 write(*,'(6f10.6)') r,gg
enddo
end program

!program test
!use localAxes_mod
!use data_types
!use printer_mod
!use ps_dms
!real(dp) a(3),b(3),c(3), d
!type(h2o) w
!real(dp) x(3,3)
!w%h1 = [ 1.80482, -1.73505, -0.85955]!H[,,]
!w%h2 = [ 2.42260, -0.35826, -0.70808]!H[,,]
!w%o  = [ 1.66322, -0.92324, -0.35186]!O[,,]
!b = [0.0_dp, 1.0_dp, 1.0_dp]
!c = [4.0_dp, 2.0_dp, 2.0_dp]
!a = cross(b,c)
!d = dot(b,c)
!
!call printer( a,'a')
!call printer( d,'d')
!end program

