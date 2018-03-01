module localAxes_mod
use data_types
!use printer_mod
implicit none


!private
!public

contains !/////////////////////////////////////////////////////////////



subroutine create_xyz_hho_new(ra,xyz_hho,nM)
    real(dp), intent(in)   :: ra(3,3*nM) !j  (xyz,hho.nM) ...
    real(dp), intent(out) :: xyz_hho(3,3,nM)!xyz,hho,nM
    integer, intent(in)    :: nM
    integer io,ih1, ih2, m, j!, posi
    real(dp), parameter :: oh_max_A = 1.2_dp !Angstrom
    
    do m = 1,nM !molecules
       ih1 = (m-1)*3+1
       ih2 = (m-1)*3+2
       io = (m-1)*3+3
       
       do j = 1,3 !coords
         
         if( abs(ra(j,io)-ra(j,ih1)) > oh_max_A )write(6,'(a,I3)') 'Long O--H1 in molecule',m
         if( abs(ra(j,io)-ra(j,ih2)) > oh_max_A )write(6,'(a,I3)') 'Long O--H2 in molecule',m
         
         xyz_hho(j,1,m) = ra(j,ih1) !H1
         xyz_hho(j,2,m) = ra(j,ih2) !H2
         xyz_hho(j,3,m) = ra(j,io )  !O
         
       enddo 
     enddo 
   end subroutine create_xyz_hho_new

     
!         rw(m)%r1(1) = ra(iH1+1)
!         rw(m)%r1(2) = ra(iH1+2)
!         rw(m)%r1(3) = ra(iH1+3)
!         rw(m)%r2(1) = ra(iH2+1)
!         rw(m)%r2(2) = ra(iH2+2)
!         rw(m)%r2(3) = ra(iH2+3)
!         rw(m)%ro(1) = ra(io+1)
!         rw(m)%ro(2) = ra(io+2)
!         rw(m)%ro(3) = ra(io+3)

  !----------------------------------------------------------------------+
  !     Routine to calculate the center of mass of each molecule         |
  !----------------------------------------------------------------------+

  subroutine get_cm(rw, rcm, nM) !w=water molecules
  ! Calculates the centers of mass
    implicit none 
    real(dp),  intent(in)   :: rw(3,3,nM)!xyz,hho,nM
    real(dp),  intent(out)  :: rcm(3,nM)!xyz,nM
    integer,   intent(in)   :: nM
    integer m
  
     do m = 1, nM !1 and 16 are weight of H and O
        rcm(:,m) = (1.0_dp*rw(:,1,m) + 1.0_dp*rw(:,2,m) + 16.0_dp * rw(:,3,m)) / 18.0_dp
     end do
     
  end subroutine get_cm


!subroutine forceOnAtoms(fCM,w,f,rCM)
!
!end subroutine

! The damping functions g1,g2... of Anthony Stone
subroutine gs(s,g)
 real(dp), intent(in) :: s
 real(dp), intent(out) :: g(5)
 real(dp) eee, prefac!pi, 
 integer k
 real(dp), parameter :: two_sqrt_pi = 2.0_dp/dsqrt(3.1415926535_dp)
 
 prefac = two_sqrt_pi*dexp(-s**2) 
 eee = erf(s)
 
 !k=0
 !g(1) = eee - prefac * ( 2.0_dp**k * s**(2*k+1) ) / dble(facfac(2*k+1))
 g(1) = eee - prefac * s 
 
 do k = 1,4
    g(k+1) = g(k) - prefac * ( 2.0_dp**k * s**(2*k+1) ) / dble(facfac(2*k+1))
 enddo
 
end subroutine


! This is the double factorial
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


subroutine force_and_torque_on_atoms(tau,fCM,w,f,rCM)
 integer, parameter :: aim = 3
 real(dp), intent(in) :: tau(3), rCM(3), fCM(3)
 !type(h2o), intent(in) :: w
 real(dp), intent(in)  :: w(3,aim) !xyz,hho
 real(dp), intent(out) :: f(3,aim) !xyz,hho
 !type(h2o), intent(inout) :: f
 real(dp) :: txr(3,aim), r(3,aim), ft(3,aim), ff(3,aim)
 real(dp) :: I(aim), Itot, t2
 real(dp), parameter :: m(3) = [1.0_dp, 1.0_dp, 16.0_dp], Mm = m(1)+m(2)+m(3)![j sum(m) didn't work with flang] !are these in the right units??? does i matter? dont think so
 integer a
 
   t2 = norm_square(tau) !length square of torque

   r(:,1) = w(:,1)-rCM   !w%h1 - rCM !extract oh distances
   r(:,2) = w(:,2)-rCM   !w%h2 - rCM
   r(:,3) = w(:,3)-rCM   !w%o  - rCM
!call printer(r,'r')   

   do a = 1,aim
      txr(:,a) = cross(tau,r(:,a))
      I(a) = m(a)/t2*norm_square(txr(:,a))
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
   
   f = ft+ff
   
   !f%h1 = f%h1   + ft(:,1) + ff(:,1)
   !f%h2 = f%h2   + ft(:,2) + ff(:,2)
   !f%o  = f%o    + ft(:,3) + ff(:,3)
!call printer(f,'f')
   
end subroutine








subroutine bisectorAxes(w,rot)
 real(dp), intent(in)  :: w(3,3) !xyz,hho !generalize with xyz,aim (aim=number of atoms in molecule)
 real(dp),  intent(out) :: rot(3,3)
 real(dp) r1(3),r2(3), n(3), bis(3), x(3), y(3), z(3), frac12

   r1 = w(:,1) - w(:,3)!w%h1 - w%o    ! O--H1 and O--H2 vectors
   r2 = w(:,2) - w(:,3)!w%h2 - w%o
   
   frac12 = norm(r1) / norm(r2)
   
   bis = r1 + r2*frac12 ! bisector between OH-vectors
   bis = normalize(bis)
   z = -bis
   
   n = cross(r2,r1)  ! normal vector to plane of molecule 
   n = normalize(n)
   y = n
   
   x = cross(y,z)
   
   rot(:,1) = x
   rot(:,2) = y
   rot(:,3) = z
end subroutine

subroutine dipoleAxes(d,w,rot)
 real(dp),  intent(in)  :: d(3)
 real(dp), intent(in)  :: w(3,3) !xyz,hho !generalize with xyz,aim (aim=number of atoms in molecule)
 real(dp),  intent(out) :: rot(3,3)
 real(dp),dimension(3) :: x,y,z, n, oh1, oh2
 
   oh1 = w(:,1) - w(:,3)!w%h1 - w%o    ! O--H1 and O--H2 vectors
   oh2 = w(:,2) - w(:,3)!w%h2 - w%o

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


subroutine plusAxes(w,rot)
 real(dp), intent(in)  :: w(3,3) !xyz,hho !generalize with xyz,aim (aim=number of atoms in molecule)
 real(dp),  intent(out) :: rot(3,3)
 real(dp) r1(3),r2(3), x(3), y(3), z(3)!n(3), bis(3), 

   r1 = w(:,1) - w(:,3)!w%h1 - w%o    ! O--H1 and O--H2 vectors
   r2 = w(:,2) - w(:,3)!w%h2 - w%o
   
   
   z = - normalize(r2+r1)  !z-axis in oposite dipole direction (in the plane in PS vibdip)
   y = normalize( cross(r2,r1) )
   x = cross(y,z)
   
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

recursive function norm_square(d) result(d2)
 real(dp), intent(in) :: d(3)
 real(dp) :: d2
   d2 = d(1)*d(1) + d(2)*d(2) + d(3)*d(3)
end function

recursive function norm(d) result(d1)
 real(dp), intent(in) :: d(3)
 real(dp) :: d1
   d1 = dsqrt(sum(d**2))
end function



end module

