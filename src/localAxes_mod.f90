module localAxes_mod
use data_types
use printer_mod
implicit none


!private
!public

contains

!subroutine forceOnAtoms(fCM,w,f,rCM)
!
!end subroutine

subroutine gs(s,g)
 real(dp), intent(in) :: s
 real(dp), intent(out) :: g(5)
 real(dp) eee, pi, prefac
 integer k
 eee = erf(s)
 pi = 3.1415926535_dp
 prefac = 2.0_dp/dsqrt(pi)*dexp(-s**2) 
 
 k=0
 g(1) = eee - prefac * ( 2.0_dp**k * s**(2*k+1) ) / dble(facfac(2*k+1))
 
 do k = 1,4
    g(k+1) = g(k) - prefac * ( 2.0_dp**k * s**(2*k+1) ) / dble(facfac(2*k+1))
 enddo
 
end subroutine


recursive function facfac(n) result(res)
 integer, intent(in) :: n
 integer :: res, nn
   res = n
   nn = n-2
   do while (nn>1)
      res = res*nn
      nn = nn-2
   enddo
end function



subroutine force_torqueOnAtoms00(tau,fCM,w,f,rCM)
 real(dp), intent(in) :: tau(3), rCM(3), fCM(3)
 type(h2o), intent(in) :: w
 type(h2o), intent(inout) :: f
 real(dp), dimension(3) :: txr1, txr2, txro, r1, r2, ro
 real(dp) :: I1, I2, Io, Itot, t2
 real(dp), parameter :: mh = 1.0_dp, mo = 16.0_dp, Mm = 2*mh+mo !are these in the right units???
 integer ii
 
   t2 = norm_2(tau)

   r1 = w%h1 - rCM
   r2 = w%h2 - rCM
   ro = w%o  - rCM
   
   
   txr1 = cross(tau,r1)
   txr2 = cross(tau,r2)
   txro = cross(tau,ro)
   
   I1 = mh/t2*norm_2(txr1)
   I2 = mh/t2*norm_2(txr2)
   Io = mo/t2*norm_2(txro)
   
   Itot = I1+I2+Io
   
   f%h1 = f%h1 + mh/Itot*txr1
   f%h2 = f%h2 + mh/Itot*txr2
   f%o  = f%o  + mo/Itot*txro
   
   ! End of torque
   f%h1 = f%h1 + fCM*mh/Mm
   f%h2 = f%h2 + fCM*mh/Mm
   f%o  = f%o  + fCM*mo/Mm
   
end subroutine

subroutine force_torqueOnAtoms(tau,fCM,w,f,rCM)
 integer, parameter :: aim = 3
 real(dp), intent(in) :: tau(3), rCM(3), fCM(3)
 type(h2o), intent(in) :: w
 type(h2o), intent(inout) :: f
 real(dp) :: txr(3,aim), r(3,aim), ft(3,aim), ff(3,aim)
 real(dp) :: I(aim), Itot, t2
 real(dp), parameter :: m(3) = [1.0_dp, 1.0_dp, 16.0_dp], Mm = sum(m) !are these in the right units???
 integer a
 
   t2 = norm_2(tau) !length square of torque

   r(:,1) = w%h1 - rCM !extract oh distances
   r(:,2) = w%h2 - rCM
   r(:,3) = w%o  - rCM
!call printer(r,'r')   

   do a = 1,aim
      txr(:,a) = cross(tau,r(:,a))
      I(a) = m(a)/t2*norm_2(txr(:,a))
   enddo
!call printer(txr,'txr')   
   
   Itot = sum(I)
   
   do a = 1,aim
      ft(:,a) = m(a)/Itot*txr(:,a) ! torque on atoms
      ff(:,a) = fCM*m(a)/Mm        ! force distributed on atoms
   enddo
!call printer(ft,'ft')   
!call printer(ff,'ff')   
   ! end of torque
   
   
   
   f%h1 = f%h1 + ft(:,1) + ff(:,1)
   f%h2 = f%h2 + ft(:,2) + ff(:,2)
   f%o  = f%o  + ft(:,3) + ff(:,3)
!call printer(f,'f')
   
end subroutine





subroutine localAxes(d,w,rot)
 real(dp),  intent(in)  :: d(3)
 type(h2o), intent(in)  :: w
 real(dp),  intent(out) :: rot(3,3)
 real(dp),dimension(3) :: x,y,z, n, oh1, oh2
 
   oh1 = w%h1 - w%o    ! O--H1 and O--H2 vectors
   oh2 = w%h2 - w%o

   n = cross(oh2,oh1)  ! normal vector to plane of molecule 
   n = normalize(n)
   
   z = - normalize(d)  !z-axis in oposite dipole direction (in the plane in PS vibdip)
   x = cross(n,z)      ! x vector (in plane)
   
   y = n               ! If dipole is in plane (it is in PS vibdip)
   !y = cross(x,z)     ! If dipole is not in plane
   
   rot(:,1) = x
   rot(:,2) = y
   rot(:,3) = z
   
end subroutine

recursive function cross(b,c) result(a)
 real(dp), intent(in) :: b(3), c(3)
 real(dp)             :: a(3) 
   a(1) = b(2)*c(3) - b(3)*c(2)
   a(2) = b(3)*c(1) - b(1)*c(3)
   a(3) = b(1)*c(2) - b(2)*c(1)
end function

recursive function dot(b,c) result(a)
 real(dp), intent(in) :: b(3), c(3)
 real(dp)             :: a
   a = b(1)*c(1) + b(2)*c(2) + b(3)*c(3)
end function

recursive function normalize(b) result(res)
 real(dp), intent(in) :: b(3)
 real(dp)             :: res(3)
   res = b*( 1.0_dp/dsqrt(sum(b**2)) )
end function

recursive function norm_2(d) result(d2)
 real(dp), intent(in) :: d(3)
 real(dp) :: d2
   d2 = d(1)*d(1) + d(2)*d(2) + d(3)*d(3)
end function

end module

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

