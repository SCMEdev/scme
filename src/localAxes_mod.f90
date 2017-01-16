module localAxes_mod
use data_types
implicit none


!private
!public

contains

!subroutine forceOnAtoms(fCM,w,f,rCM)
!
!end subroutine

subroutine torqueOnAtoms(tau,w,f,rCM)
 real(dp), intent(in) :: tau(3), rCM(3)
 type(h2o), intent(in) :: w
 type(h2o), intent(inout) :: f
 real(dp), dimension(3) :: txr1, txr2, txro, r1, r2, ro
 real(dp) :: I1, I2, Io, Itot, t2
 real(dp), parameter :: mh = 1.0_dp, mo = 16.0_dp !are these in the right units???
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
   
   f%h1 = mh/Itot*txr1
   f%h2 = mh/Itot*txr2
   f%o  = mo/Itot*txro
   
end subroutine



subroutine test_xyz()
type(xyz) r(2)
r(1)%x = [1.2_dp, 3.4_dp,5.6_dp]
r(2)%x = [4.2_dp, 6.4_dp,1.6_dp]
print*, r(1)%x
print*, r(2)%x
end subroutine

!type(hho) :: w(nM), aforces(nM)
! call torqueOnAtoms2(tau(:,m), rCM(:,m), w(m)%a, fatom(m)%a

subroutine torqueOnAtoms2(tau,rCM,rhho,f)
 integer, parameter :: a_m = 3
 type(xyz), intent(in) :: tau, rCM
 type(xyz), intent(in) :: rhho(a_m)
 type(xyz), intent(inout) :: f(a_m)
 type(xyz) :: r(a_m), txr(a_m)
 real(dp) :: I(a_m), Itot, t2
 !real(dp), dimension(3) :: txr1, txr2, txro, r1, r2, ro
 !real(dp) :: I1, I2, Io, Itot, t2
 real(dp), parameter :: m(a_m) = [1.0_dp, 1.0_dp, 16.0_dp] !are these in the right units???
 integer ii, a
 
   t2 = norm_2(tau%x)
   
   do a = 1,a_m !atoms hho
      r(a)%x = rhho(a)%x - rCM%x
      txr(a)%x = cross(tau%x,r(a)%x)
      I(a) = m(a)/t2*norm_2(txr(a)%x)
   enddo
   
   Itot = sum(I)
   
   do a = 1,a_m
      f(a)%x = m(a)/Itot*txr(a)%x
   enddo
   
end subroutine

subroutine torqueOnAtoms3(tau,rCM,mol,f)
 integer, parameter :: aim = 3
 real(dp), intent(in) :: tau(3), rCM(3)
 real(dp), intent(in) :: mol(3,aim)
 real(dp), intent(inout) :: f(3,aim)
 real(dp) :: r(3), txr(3,aim)
 real(dp) :: I(aim), Itot, t2
 real(dp), parameter :: m(aim) = [1.0_dp, 1.0_dp, 16.0_dp] !are these in the right units???
 integer ii, a
 
   t2 = norm_2(tau)
   
   do a = 1,aim !atoms hho
      r = mol(:,a) - rCM
      txr(:,a) = cross(tau,r)
      I(a) = m(a)/t2*norm_2(txr(:,a))
   enddo
   
   Itot = sum(I)
   
   do a = 1,aim
      f(:,a) = m(a)/Itot*txr(:,a)
   enddo
   
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

