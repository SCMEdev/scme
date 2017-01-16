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
