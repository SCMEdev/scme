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



! Prints h2o type in fa order. 
! It takes type(h2o) and outputs a vector: [H1x H1y H1z H2x H2y H2z H1x H1y H1z H2x H2y H2z ... Ox Oy Oz ...]
subroutine printer_h2o_linear(a,text)
  character(*) text
  type(h2o) :: a(:)
  integer l, m, xyz
   call printer_text(text)
   l = size(a)
   
   do m = 1,l
     do xyz = 1,3
       write(un,'(a,I3,a,3'//dblf//')') '   Mol:',m,' H1:',a(m)%h1(xyz)
     enddo
     do xyz = 1,3
       write(un,'(a,I3,a,3'//dblf//')') '   Mol:',m,' H2:',a(m)%h2(xyz)
     enddo
   enddo
   do m = 1,l
     do xyz = 1,3
       write(un,'(a,I3,a,3'//dblf//')') '   Mol:',m,' O :',a(m)%o(xyz)
     enddo
   enddo

end subroutine 
 
subroutine h2o_to_linear(a,alin,nM)
  integer, intent(in) :: nM
  type(h2o), intent(in) :: a(nM)
  real(dp) , intent(out) :: alin(nM*9)
  integer l, m, xyz
  alin=0
   do m = 1,nM
     do xyz = 1,3
       alin((2*m-2)*3 + xyz) = a(m)%h1(xyz)
       alin((2*m-1)*3 + xyz) = a(m)%h2(xyz)
       alin((2*nM+m-1)*3 + xyz)  = a(m)%o(xyz)
     enddo
   enddo
end subroutine 
 
 
 
  subroutine new_dispersion(rw, fa, uDisp, nM, a, a2)
    
    implicit none
    
    type(h2o), intent(in)  :: rw(:)
    type(h2o), intent(out) :: fa(:)
    real(dp),  intent(out) :: uDisp
    integer,   intent(in)  :: nM
    real(dp),  intent(in)  :: a(3), a2(3) 

!JÖ internal    
    real(dp) t1, t2, df, r
    real(dp) r2, r6, r7, r8, r9, r10, r11
    real(dp) f6, df6, f8, df8, f10, df10
    
    integer n, m, iOn, iOm, i
    real(dp) dr(3), sc
    real(dp), parameter :: C6 = 46.4430d0 * 0.597527378d0
    real(dp), parameter :: C8  = 1141.7000d0 * 0.167324732d0
    real(dp), parameter :: C10 = 33441.0000d0 * 0.046855703d0

    
    uDisp = 0.0_dp
    do n = 1, nM-1
       do m = n+1, nM
          dr = rw(m)%o - rw(n)%o !Om-On
          do i = 1, 3
             if      (dr(i) .gt. a2(i))  then; dr(i) = dr(i) - a(i)
             else if (dr(i) .lt. -a2(i)) then; dr(i) = dr(i) + a(i)
             end if
          end do
          
          r2 = sum(dr**2)
          r = sqrt(r2)
          call tang_toennies_disp(r, f6, df6, f8, df8, f10, df10)
          
          r6 = r**6 !probably better with r2**3
          r7 = r6 * r
          r8 = r7 * r
          r9 = r8 * r
          r10 = r9 * r
          r11 = r10 * r
          
          uDisp = uDisp - C6/r6*f6 - C8/r8*f8 - C10/r10*f10
          
          do i = 1, 3
             df = -C6 * (6.0_dp * f6 / r7 - df6 / r6)
             df = df - C8 * (8.0_dp * f8 / r9 - df8 / r8)
             df = df - C10 * (10.0_dp * f10 / r11 - df10 / r10)
             df = df * dr(i) / r
             !fa(iOn+i) = fa(iOn+i) - df
             !fa(iOm+i) = fa(iOm+i) + df
             fa(n)%o(i) = fa(n)%o(i) - df
             fa(m)%o(i) = fa(m)%o(i) + df
             
          end do
       end do
    end do
    
    return
    
  end subroutine 



! has the h2o type, so not used
!  subroutine create_rw(ra,rw,nM)
!    real(dp), intent(in)   :: ra(:)
!    type(h2o), intent(out) :: rw(:)
!    integer, intent(in) :: nM
!    integer iH2, iH1, iO, m, j
!    
!    do m = 1,nM !molecules
!       iH2 = 2 * m
!       iH1 = iH2 - 1
!       iO  = 2*nM + m !JÖ
!
!       iH1 = 3 * (iH1 - 1)
!       iH2 = 3 * (iH2 - 1)
!       iO  = 3 * (iO - 1)
!       
!       do j = 1,3 !coords
!         rw(m)%h1(j) = ra(iH1+j) !H1
!         rw(m)%h2(j) = ra(iH2+j) !H2
!         rw(m)%o(j) = ra(io+j)  !O
!       enddo 
!     enddo 
!   end subroutine create_rw
  
! This one has to be wrong since ra(3*nM) has too few entries  
!  subroutine create_xyz(ra,xyz,nM)
!    real(dp), intent(in)   :: ra(3*nM)
!    real(dp), intent(out) :: xyz(3,nM)
!    integer, intent(in)    :: nM
!    integer iH2, iH1, iO, m, j, posi
!    real(dp), parameter :: oh_max_A = 2.0_dp !Angstrom
!    
!    do m = 1,nM !molecules
!       iH2 = 2 * m
!       iH1 = iH2 - 1
!       iO  = 2*nM + m !JÖ
!
!       iH1 = 3 * (iH1 - 1)
!       iH2 = 3 * (iH2 - 1)
!       iO  = 3 * (iO - 1)
!       
!       posi = 3*(m-1)
!       do j = 1,3 !coords
!         
!         if( abs(ra(io+j)-ra(iH1+j)) > oh_max_A )write(6,'(a,I3)') 'Long O--H1 in molecule',m
!         if( abs(ra(io+j)-ra(iH2+j)) > oh_max_A )write(6,'(a,I3)') 'Long O--H2 in molecule',m
!         
!         xyz(j,posi+1) = ra(iH1+j) !H1
!         xyz(j,posi+2) = ra(iH2+j) !H2
!         xyz(j,posi+3) = ra(io+j)  !O
!         
!       enddo 
!     enddo 
!   end subroutine create_xyz

! needs to be wring since ra(3*nM) is too small
!  subroutine create_xyz_hho(ra,xyz_hho,nM)
!    real(dp), intent(in)   :: ra(3*nM)
!    real(dp), intent(out) :: xyz_hho(3,3,nM)!xyz,hho,nM
!    integer, intent(in)    :: nM
!    integer iH2, iH1, iO, m, j, posi
!    real(dp), parameter :: oh_max_A = 2.0_dp !Angstrom
!    
!    do m = 1,nM !molecules
!       iH2 = 2 * m
!       iH1 = iH2 - 1
!       iO  = 2*nM + m !JÖ
!
!       iH1 = 3 * (iH1 - 1)
!       iH2 = 3 * (iH2 - 1)
!       iO  = 3 * (iO - 1)
!       
!       do j = 1,3 !coords
!         
!         if( abs(ra(io+j)-ra(iH1+j)) > oh_max_A )write(6,'(a,I3)') 'Long O--H1 in molecule',m
!         if( abs(ra(io+j)-ra(iH2+j)) > oh_max_A )write(6,'(a,I3)') 'Long O--H2 in molecule',m
!         
!         xyz_hho(j,1,m) = ra(iH1+j) !H1
!         xyz_hho(j,2,m) = ra(iH2+j) !H2
!         xyz_hho(j,3,m) = ra(io+j)  !O
!         
!       enddo 
!     enddo 
!   end subroutine create_xyz_hho


!  subroutine calc_cm(rw, wcm, nM) !w=water molecules
!  ! Calculates the centers of mass for h2o type
!  ! do not use
!    implicit none 
!    type(h2o), 	intent(in) 	:: rw(:)
!    real(dp), 	intent(out) :: wcm(:,:)
!    integer, 	intent(in) 	:: nM
!!internal:    
!    integer m
!  
!     do m = 1, nM
!        wcm(:,m) = (rw(m)%h1 + rw(m)%h2 + 16.0_dp * rw(m)%o) / 18.0_dp
!     end do
!     
!  end subroutine calc_cm



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

  subroutine calcCentersOfMass(ra, nM, rCM)

!JÖ    implicit real(dp) (a-h,o-z)
    implicit none !JÖ
    !JÖ state intents and assumed shape:
    real(dp), 	intent(in) 	:: ra(:)
    integer, 	intent(in) 	:: nM
    real(dp), 	intent(out) :: rCM(:,:)
    
    integer iH1, iH2, iO, i, j !JÖ nM, 
!JÖ    real(dp) rCM(3,maxCoo/3), ra(maxCoo)

    do i = 1, nM
       iH2 = 2 * i
       iH1 = iH2 - 1
       iO  = 2*nM + i !JÖ

       iH1 = 3 * (iH1 - 1)
       iH2 = 3 * (iH2 - 1)
       iO  = 3 * (iO - 1)
       do j = 1, 3
          rCM(j,i) = (ra(iH1+j) + ra(iH2+j) + 16.0_dp * ra(iO+j)) / 18.0_dp
       end do
    end do

    return

  end subroutine calcCentersOfMass

  !----------------------------------------------------------------------+
  !     Routine to calculate the principal axes of each molecule         |
  !----------------------------------------------------------------------+
  subroutine findPpalAxesOLD(ra, nM, x)

!JÖ    implicit real(dp) (a-h,o-z)
    implicit none
    
    real(dp), 	intent(in) 	:: ra(:)
    integer, 	intent(in) 	:: nM
    real(dp), 	intent(out) :: x(:,:,:) 

    integer 	:: iH1, iH2, iO, i, j !JÖ nM, 
    real(dp) 	:: r11, r21
!jö    real(dp) x(3,3,maxCoo/3), ra(maxCoo), r11, r21


    do i = 1, nM
       iH2 = 2 * i
       iH1 = iH2 - 1
       iO  = 2 * nM + i

       iH1 = 3 * (iH1 - 1)
       iH2 = 3 * (iH2 - 1)
       iO  = 3 * (iO - 1)
       do j = 1, 3 ! j = x,y,z coord of H1,H2,O
          x(j,1,i) = -(ra(iH1+j) + ra(iH2+j) - 2.0_dp * ra(iO+j)) !x(:,1,i) = -{ H1(:) + H2(:) - 2*O(:) } 
          x(j,2,i) = ra(iH2+j) - ra(iH1+j)                        !x(:,2,i) = H2(:) - H1(:) 
       end do
       r11 = sqrt(x(1,1,i)*x(1,1,i) + x(2,1,i)*x(2,1,i) + x(3,1,i) * x(3,1,i)) !norm(x(:,1,i)),  ":" is x,y,z
       r21 = sqrt(x(1,2,i)*x(1,2,i) + x(2,2,i)*x(2,2,i) + x(3,2,i) * x(3,2,i)) !norm(x(:,2,i))
       do j = 1, 3
          x(j,1,i) = x(j,1,i) / r11
          x(j,2,i) = x(j,2,i) / r21
       end do
       x(1,3,i) = x(2,1,i) * x(3,2,i) - x(3,1,i) * x(2,2,i)
       x(2,3,i) = x(3,1,i) * x(1,2,i) - x(1,1,i) * x(3,2,i)
       x(3,3,i) = x(1,1,i) * x(2,2,i) - x(2,1,i) * x(1,2,i)
    end do
    return

  end subroutine findPpalAxesOLD

  !----------------------------------------------------------------------+
  !     Routine to calculate the principal axes of each molecule         |
  !----------------------------------------------------------------------+
  subroutine findPpalAxes(ra, nM, x)
  !---------------------------------------------------------------------
  ! This routine finds the principal axis of the molecule but does it 
  ! apply to a water molecule with non-fixed geometry?
  !---------------------------------------------------------------------

!JÖ    implicit real(dp) (a-h,o-z)
!JÖ
!JÖ    integer nM, iH1, iH2, iO, i, j
!JÖ    real(dp) x(3,3,maxCoo/3), ra(maxCoo), r11, r21

!JÖ use instead:
    implicit none
    
    real(dp), 	intent(in) 	:: ra(:)
    integer, 	intent(in) 	:: nM
    real(dp), 	intent(out) :: x(:,:,:) 

    integer 	:: iH1, iH2, iO, i, j !JÖ nM, 
    real(dp) 	:: r11, r21


    ! Debug
    integer*4 p

    do i = 1, nM
       iH2 = 2 * i
       iH1 = iH2 - 1
       iO  = 2 * nM + i

       iH1 = 3 * (iH1 - 1)
       iH2 = 3 * (iH2 - 1)
       iO  = 3 * (iO - 1)
       do j = 1, 3 ! j = x,y,z coord of H1,H2,O
          x(j,3,i) = -(ra(iH1+j) + ra(iH2+j) - 2.0_dp * ra(iO+j)) !x(:,1,i) = -{ H1(:) + H2(:) - 2*O(:) } |
          x(j,1,i) = ra(iH2+j) - ra(iH1+j)                        !x(:,2,i) = H2(:) - H1(:)               | These find two vectors in the plane of the molecule
       end do
       r11 = sqrt(x(1,3,i)*x(1,3,i) + x(2,3,i)*x(2,3,i) + x(3,3,i) * x(3,3,i)) !norm(x(:,1,i)),  ":" is x,y,z
       r21 = sqrt(x(1,1,i)*x(1,1,i) + x(2,1,i)*x(2,1,i) + x(3,1,i) * x(3,1,i)) !norm(x(:,2,i))

       do j = 1, 3
          x(j,3,i) = x(j,3,i) / r11 !normalize?
          x(j,1,i) = x(j,1,i) / r21
       end do
       x(1,2,i) = x(2,3,i) * x(3,1,i) - x(3,3,i) * x(2,1,i) !cross product ? find middle vector orthogonal to plane
       x(2,2,i) = x(3,3,i) * x(1,1,i) - x(1,3,i) * x(3,1,i) !do we need all three of them? Yes!
       x(3,2,i) = x(1,3,i) * x(2,1,i) - x(2,3,i) * x(1,1,i)
    end do

    ! Debug
    ! Print the principal axes matrix for each molecule
    !     do i=1,nM
    !       write(10,FMT='(A,I4)') ' Axes for Molecule: ', i
    !       write(10,FMT='(3F10.5)') (x(1,p,i),p=1,3)
    !       write(10,FMT='(3F10.5)') (x(2,p,i),p=1,3)
    !       write(10,FMT='(3F10.5)') (x(3,p,i),p=1,3)
    !     end do

    return

  end subroutine findPpalAxes


  !-----------------------------------------------------------------------
  !        subroutine applyPBC(r, sep, j)
  !
  !                ! not used? I commented it out
  !                ! to use: uncomment, add lattice to argument list and
  !                ! change ax, ay, az to lattice(1), lattice(2), lattice(3)
  !
  !                implicit real(dp) (a-h,o-z)
  !                real(dp) r(3), a(3), sep(3)
  !                integer j
  !
  !                !     Size of the simulation cell
  !                a(1) = ax/2.0d0
  !                a(2) = ay/2.0d0
  !                a(3) = az/2.0d0
  !
  !                do i = 1, 3
  !                        !         if (j*sep(i) .gt. a(i)) then
  !                        if (r(i) .gt. a(i)) then
  !                                r(i) = r(i) - 2.0d0 * a(i)
  !                                !         elseif (j*sep(i) .lt. -a(i)) then
  !                        elseif (r(i) .lt. -a(i)) then
  !                                r(i) = r(i) + 2.0d0 * a(i)
  !                        end if
  !                end do
  !
  !                j = 0
  !
  !                return
  !
  !        end subroutine applyPBC

  !-----------------------------------------------------------------------
  pure subroutine SF(r, swFunc)

    implicit none
    real(dp), intent(in) :: r !JÖ
    real(dp), intent(out) :: swFunc !JÖ
    real(dp) x, x2, x3, rSW, rCut, dr     !JÖ r, swFunc, 
    !JÖ , save ::
!JÖ    real(dp), save :: rL1, rL2, rH1, rH2
    real(dp), parameter :: rL1 = 0.0_dp
    real(dp), parameter :: rL2 = 9.0_dp
    real(dp), parameter :: rH1 = 5.0_dp
    real(dp), parameter :: rH2 = 11.0_dp
    
    !      data rL1, rH1, rL2, rH2 / 1.5d0, 2.7d0, 8.d0, 9.d0 /
!JÖ    data rL1, rH1, rL2, rH2 / 0.0_dp, 5.0_dp, 9.0_dp, 11.0_dp /
    !      data rL1, rH1, rL2, rH2 / 0.0_dp, 5.0_dp, 11.0_dp, 13.0_dp /
!JÖ     save


    !     Kroes
    !$$$      rSW = 8.d0
    !$$$      rCut = 9.d0
    !$$$
    !$$$      x = (r - rSW)/(rCut - rSW)
    !$$$
    !$$$      if (r. lt. rSW) then
    !$$$         swFunc = 1.0d0
    !$$$      else if(r .lt. rCut) then
    !$$$         x2 = x * x
    !$$$         x3 = x2 * x
    !$$$         swFunc = 1.d0 + x3 * (-6.d0 * x2 + 15.d0 * x - 10.d0)
    !$$$      else
    !$$$         swFunc = 0.0d0
    !$$$      end if

    if ((r .ge. rH2) .or. (r .le. rL1)) then
       swFunc = 0.0_dp
    else if ((r .ge. rH1) .and. (r .le. rL2)) then
       swFunc = 1.0_dp
    else if (r .lt. rH1) then

       call tang_toenniesN(r, swFunc, 6)

       !$$$         x = 1.d0 - (r - rL1)/(rH1 - rL1)
       !$$$         x2 = x * x
       !$$$         x3 = x2 * x
       !$$$         swFunc = 1.d0 + x3 * (-6.d0 * x2 + 15.d0 * x - 10.d0)
    else
       x = (r - rL2)/(rH2 - rL2)
       x2 = x * x
       x3 = x2 * x
       swFunc = 1.0_dp + x3 * (-6.0_dp * x2 + 15.0_dp * x - 10.0_dp)
    end if

    !      swFunc = 1.0_dp

    return

  end subroutine SF

  !-----------------------------------------------------------------------
  pure subroutine SFdsf(r, swFunc, dSdr)!JÖ pure

    implicit none
    real(dp), intent(in) :: r !JÖ stated intents
    real(dp), intent(out) :: swFunc, dSdr !JÖ
    real(dp) :: x, x2, x3, rSW, rCut, dr
!JÖ pure    real(dp) :: rL1, rL2, rH1, rH2  !JÖ pure: ^, dr
    
    real(dp), parameter :: rL1 = 0.0_dp   !JÖ pure
    real(dp), parameter :: rH1 = 5.0_dp   !JÖ pure
    real(dp), parameter :: rL2 = 9.0_dp   !JÖ pure
    real(dp), parameter :: rH2 = 11.0_dp  !JÖ pure
!, save
!, save
    !      data rL1, rH1, rL2, rH2 / 1.5d0, 2.7d0, 8.d0, 9.d0 /
!JÖ pure    data rL1, rH1, rL2, rH2 / 0.0_dp, 5.0_dp, 9.0_dp, 11.0_dp /
    !      data rL1, rH1, rL2, rH2 / 0.0_dp, 5.0_dp, 11.0_dp, 13.0_dp /
    !JÖ save

    !     Kroes
    !$$$      rSW = 8.d0
    !$$$      rCut = 9.d0


    !$$$      x = (r - rSW)/(rCut - rSW)
    !$$$
    !$$$      if (r. lt. rSW) then
    !$$$         swFunc = 1.0d0
    !$$$         dSdr = 0.0d0
    !$$$      else if(r .lt. rCut) then
    !$$$         x2 = x * x
    !$$$         x3 = x2 * x
    !$$$         swFunc = 1.d0 + x3 * (-6.d0 * x2 + 15.d0 * x - 10.d0)
    !$$$         dSdr = 30.d0 * x2 * (- x2 + 2.d0 * x - 1.d0) / (rCut-rSW)
    !$$$      else
    !$$$         swFunc = 0.0d0
    !$$$         dSdr = 0.0d0
    !$$$      end if


    if ((r .ge. rH2) .or. (r .le. rL1)) then
       swFunc = 0.0_dp
       dSdr = 0.0_dp
    else if ((r .ge. rH1) .and. (r .le. rL2)) then
       swFunc = 1.0_dp
       dSdr = 0.0_dp
    else if (r .lt. rH1) then

       !c         call switchCloseDsf(r, swFunc, dSdr)
       call tang_toenniesNdF(r, swFunc, dSdr, 6)

       !$$$         x = 1.d0 - (r - rL1)/(rH1 - rL1)
       !$$$         dr = - 1.d0 / (rH1 - rL1)
       !$$$         x2 = x * x
       !$$$         x3 = x2 * x
       !$$$         swFunc = 1.d0 + x3 * (-6.d0 * x2 + 15.d0 * x - 10.d0)
       !$$$         dSdr = 30.d0 * x2 * (- x2 + 2.d0 * x - 1.d0) * dr
    else
       x = (r - rL2)/(rH2 - rL2)
       dr = 1.0_dp / (rH2 - rL2)
       x2 = x * x
       x3 = x2 * x
       swFunc = 1.0_dp + x3 * (-6.0_dp * x2 + 15.0_dp * x - 10.0_dp)
       dSdr = 30.0_dp * x2 * (- x2 + 2.0_dp * x - 1.0_dp) * dr
    end if

    !      swFunc = 1.d0
    !      dSdr = 0.d0

    return

  end subroutine SFdsf

  !----------------------------------------------------------------------+
  !     Add all the fields                                               |
  !----------------------------------------------------------------------+
  subroutine addFields(eH, eD, eT, nM)

    implicit none
!    integer i, j, nM
!    real(dp) eH(3,maxCoo/3), eD(3,maxCoo/3)
!    real(dp) eT(3,maxCoo/3)

    integer , intent(in)   :: nM
    real(dp), intent(in)   :: eH(:,:)
    real(dp), intent(in)   :: eD(:,:)
    real(dp), intent(out) :: eT(:,:)
    
    integer i, j
   
    eT = eH + eD

!JÖ    do i = 1, nM
!JÖ       do j = 1, 3
!JÖ          eT(j,i) = eH(j,i) + eD(j,i)
!JÖ       end do
!JÖ    end do
    return

  end subroutine addFields



  !----------------------------------------------------------------------+
  !     Routine to rotate the multipoles to the orientation of the       |
  !     molecule                                                         |
  !----------------------------------------------------------------------+
  subroutine setUnpolPoles(dpole, qpole, dpole0, qpole0, nM)

    implicit none

    
    real(dp), intent(inout) :: dpole(:,:) !JÖ 
    real(dp), intent(inout) :: qpole(:,:,:)
    real(dp), intent(in) :: dpole0(:,:)
    real(dp), intent(in) :: qpole0(:,:,:)
    integer , intent(in) :: nM
    
    integer i, j, k
    
!JÖ    real(dp) dpole(3,maxCoo/3), qpole(3,3,maxCoo/3)
!JÖ    real(dp) dpole0(3,maxCoo/3), qpole0(3,3,maxCoo/3)

    
    dPole = dpole0 !JÖ 
    qpole = qpole0 !JÖ

!JÖ    do i = 1, nM 
!JÖ       do j = 1, 3
!JÖ          dPole(j,i) = dpole0(j,i)
!JÖ          do k = 1, 3
!JÖ             qpole(k,j,i) = qpole0(k,j,i)
!JÖ          end do
!JÖ       end do
!JÖ    end do
    return

  end subroutine setUnpolPoles
  
      !----------------------------------------------------------------------+
  !> Routine to rotate the multipoles to the orientation of the molecules.
  !! @param[in] d0     : ?
  !! @param[in] q0     : ?
  !! @param[in] o0     : ?
  !! @param[in] h0     : ?
  !! @param[out] dpole : The array is assumed to be zeroized.
  !! @param[out] qpole : ?
  !! @param[out] opole : ?
  !! @param[out] hpole : ?
  !! @param[out] nM    : The number of molecules.
  !! @param[out] x     : ?
subroutine rotatePoles(d0, q0, o0, h0, dpole, qpole, opole, hpole, nM, x)
    implicit none
    real(dp), intent(in) :: d0(3) !JÖ assumed shape can be faster
    real(dp), intent(in) :: q0(3,3)
    real(dp), intent(in) :: o0(3,3,3)
    real(dp), intent(in) :: h0(3,3,3,3)
    integer, intent(in) :: nM
    real(dp), intent(in) :: x(:,:,:)

    real(dp), intent(out) :: dpole(:,:)
    real(dp), intent(out) :: qpole(:,:,:)
    real(dp), intent(out) :: opole(:,:,:,:)
    real(dp), intent(out) :: hpole(:,:,:,:,:)



!    real(dp), intent(in) :: d0(3)
!    real(dp), intent(in) :: q0(3,3)
!    real(dp), intent(in) :: o0(3,3,3)
!    real(dp), intent(in) :: h0(3,3,3,3)
!    integer, intent(in) :: nM
!    real(dp), intent(in) :: x(3,3,nM)
    ! ----------------------------------------

!    real(dp), intent(out) :: dpole(3,nM)
!    real(dp), intent(out) :: qpole(3,3,nM)
!    real(dp), intent(out) :: opole(3,3,3,nM)
!    real(dp), intent(out) :: hpole(3,3,3,3,nM)


    ! local variables.
    integer i, j, k, l, ii, jj, kk, ll, m !jö , save ::

    ! Zeroise.
!JÖ    hpole(1:3, 1:3, 1:3, 1:3, 1:nM) = 0.0_dp
!JÖ    opole(1:3, 1:3, 1:3, 1:nM) = 0.0_dp
!JÖ    qpole(1:3, 1:3, 1:nM) = 0.0_dp
!JÖ    dpole(1:3, 1:nM) = 0.0_dp

    hpole = 0 !JÖ this is enough 
    opole = 0
    qpole = 0
    dpole = 0
    ! Do the calculation.
    do m = 1, nM
       do l = 1, 3
          do ll = 1, 3
             do k = 1, 3
                do kk = 1, 3
                   do j = 1, 3
                      do jj = 1, 3
                         do i = 1, 3
                            do ii = 1, 3
                               hpole(i,j,k,l,m) = hpole(i,j,k,l,m) + &
                                    x(i,ii,m) * x(j,jj,m) * x(k,kk,m) &
                                    * x(l,ll,m) * h0(ii,jj,kk,ll)
                            end do
                         end do
                         opole(j,k,l,m) = opole(j,k,l,m) + &
                              x(j,jj,m)* x(k,kk,m)* x(l,ll,m) * o0(jj,kk,ll)
                      end do
                   end do
                   qpole(k,l,m) = qpole(k,l,m) + x(k,kk,m) * x(l,ll,m)* q0(kk,ll)
                end do
             end do
             dpole(l,m) = dpole(l,m) + x(l,ll,m) * d0(ll)
          end do
       end do
    end do

    return

  end subroutine rotatePoles

