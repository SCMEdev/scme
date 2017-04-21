module localAxes_mod
use data_types
use printer_mod
implicit none


!private
!public

contains !/////////////////////////////////////////////////////////////

  subroutine create_rw(ra,rw,nM)
    real(dp), intent(in)   :: ra(:)
    type(h2o), intent(out) :: rw(:)
    integer, intent(in) :: nM
    integer iH2, iH1, iO, m, j
    
    do m = 1,nM !molecules
       iH2 = 2 * m
       iH1 = iH2 - 1
       iO  = 2*nM + m !JÖ

       iH1 = 3 * (iH1 - 1)
       iH2 = 3 * (iH2 - 1)
       iO  = 3 * (iO - 1)
       
       do j = 1,3 !coords
         rw(m)%h1(j) = ra(iH1+j) !H1
         rw(m)%h2(j) = ra(iH2+j) !H2
         rw(m)%o(j) = ra(io+j)  !O
       enddo 
     enddo 
   end subroutine create_rw
  
  
  subroutine create_xyz(ra,xyz,nM)
    real(dp), intent(in)   :: ra(3*nM)
    real(dp), intent(out) :: xyz(3,nM)
    integer, intent(in)    :: nM
    integer iH2, iH1, iO, m, j, posi
    real(dp), parameter :: oh_max_A = 2.0_dp !Angstrom
    
    do m = 1,nM !molecules
       iH2 = 2 * m
       iH1 = iH2 - 1
       iO  = 2*nM + m !JÖ

       iH1 = 3 * (iH1 - 1)
       iH2 = 3 * (iH2 - 1)
       iO  = 3 * (iO - 1)
       
       posi = 3*(m-1)
       do j = 1,3 !coords
         
         if( abs(ra(io+j)-ra(iH1+j)) > oh_max_A )write(6,'(a,I3)') 'Long O--H1 in molecule',m
         if( abs(ra(io+j)-ra(iH2+j)) > oh_max_A )write(6,'(a,I3)') 'Long O--H2 in molecule',m
         
         xyz(j,posi+1) = ra(iH1+j) !H1
         xyz(j,posi+2) = ra(iH2+j) !H2
         xyz(j,posi+3) = ra(io+j)  !O
         
       enddo 
     enddo 
   end subroutine create_xyz

  subroutine create_xyz_hho(ra,xyz_hho,nM)
    real(dp), intent(in)   :: ra(3*nM)
    real(dp), intent(out) :: xyz_hho(3,3,nM)!xyz,hho,nM
    integer, intent(in)    :: nM
    integer iH2, iH1, iO, m, j, posi
    real(dp), parameter :: oh_max_A = 2.0_dp !Angstrom
    
    do m = 1,nM !molecules
       iH2 = 2 * m
       iH1 = iH2 - 1
       iO  = 2*nM + m !JÖ

       iH1 = 3 * (iH1 - 1)
       iH2 = 3 * (iH2 - 1)
       iO  = 3 * (iO - 1)
       
       do j = 1,3 !coords
         
         if( abs(ra(io+j)-ra(iH1+j)) > oh_max_A )write(6,'(a,I3)') 'Long O--H1 in molecule',m
         if( abs(ra(io+j)-ra(iH2+j)) > oh_max_A )write(6,'(a,I3)') 'Long O--H2 in molecule',m
         
         xyz_hho(j,1,m) = ra(iH1+j) !H1
         xyz_hho(j,2,m) = ra(iH2+j) !H2
         xyz_hho(j,3,m) = ra(io+j)  !O
         
       enddo 
     enddo 
   end subroutine create_xyz_hho

subroutine create_xyz_hho_new(ra,xyz_hho,nM)
    real(dp), intent(in)   :: ra(3,3*nM) !j  (xyz,hho.nM) ...
    real(dp), intent(out) :: xyz_hho(3,3,nM)!xyz,hho,nM
    integer, intent(in)    :: nM
    integer io,ih1, ih2, m, j, posi
    real(dp), parameter :: oh_max_A = 2.0_dp !Angstrom
    
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
  subroutine calc_cm(rw, wcm, nM) !w=water molecules
  ! Calculates the centers of mass for h2o type
  ! do not use
    implicit none 
    type(h2o), 	intent(in) 	:: rw(:)
    real(dp), 	intent(out) :: wcm(:,:)
    integer, 	intent(in) 	:: nM
!internal:    
    integer m
  
     do m = 1, nM
        wcm(:,m) = (rw(m)%h1 + rw(m)%h2 + 16.0_dp * rw(m)%o) / 18.0_dp
     end do
     
  end subroutine calc_cm

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

subroutine gs(s,g)
 real(dp), intent(in) :: s
 real(dp), intent(out) :: g(5)
 real(dp) eee, pi, prefac
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


recursive function facfac(n) result(res)
! This is the double factorial
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
 real(dp), parameter :: m(3) = [1.0_dp, 1.0_dp, 16.0_dp], Mm = sum(m) !are these in the right units??? does i matter? dont think so
 integer a
 
   t2 = norm_2(tau) !length square of torque

   r(:,1) = w(:,1)-rCM   !w%h1 - rCM !extract oh distances
   r(:,2) = w(:,2)-rCM   !w%h2 - rCM
   r(:,3) = w(:,3)-rCM   !w%o  - rCM
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
   
   f = ft+ff
   
   !f%h1 = f%h1   + ft(:,1) + ff(:,1)
   !f%h2 = f%h2   + ft(:,2) + ff(:,2)
   !f%o  = f%o    + ft(:,3) + ff(:,3)
!call printer(f,'f')
   
end subroutine





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


subroutine localAxes2(d,w,rot)
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

