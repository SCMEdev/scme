module calc_derivs
  
  use data_types
!  use max_parameters
  !use molecProperties, only: SFdsf
  use sf_disp_tangtoe, only: SFdsf
  
  implicit none
  
  private
  public calcDv
  
contains
  
  !----------------------------------------------------------------------+
  !     Calculate derivatives of the electric field.                     |
  !----------------------------------------------------------------------+
  subroutine calcDv(rCM, dpole, qpole, opole, hpole, nM, NC, a, a2, ed1, ed2, ed3, ed4, ed5, rMax2, fsf, iSlab,full)

    implicit none
    integer nM, NC, NCz
    real(dp), intent(in) ::  rCM(:,:), a(:), a2(:)
    real(dp), intent(in) ::  dpole(:,:), qpole(:,:,:)
    real(dp), intent(in) ::  opole(:,:,:,:), hpole(:,:,:,:,:)
    !, target :: dpole(:,:), qpole(:,:,:)
    !, target :: opole(:,:,:,:), hpole(:,:,:,:,:)
    logical,intent(in) :: full


    
    real(dp), intent(in) ::   rMax2
    logical*1, intent(in) ::  iSlab
    
    real(dp), intent(out) ::  fsf(:,:)
    real(dp), intent(out) ::  ed1(:,:), ed2(:,:,:), ed3(:,:,:,:)
    real(dp), intent(out) ::  ed4(:,:,:,:,:), ed5(:,:,:,:,:,:) 

                  


!JÖ internal: (trying the save option)     
    real(dp), save :: d1d(3), d2d(3,3), d3d(3,3,3)
    real(dp), save :: d4d(3,3,3,3), d5d(3,3,3,3,3)
    real(dp), save :: d1a(3), d2a(3,3), d3a(3,3,3)
    real(dp), save :: d4a(3,3,3,3), d5a(3,3,3,3,3)
    real(dp), save :: re(3), dr(3), r1, r2, swFunc, dSdr
    integer , save :: i, j, k, l, s, n, m, ii, nx, ny, nz
    integer , save :: in2(2), in3(3), in4(4), in5(5)
    

!JÖ    real(dp) :: d(3), q(3,3), o(3,3,3), h(3,3,3,3)
!JÖ    real(dp), pointer :: d(:), q(:,:), o(:,:,:), h(:,:,:,:)

!, save
!, save
!
!, save
!, save
!
!, save
!, save
!
!, save
!, save
! Save is no improvement, probably because this routine is to high up  
! and the save only consumes cache for the lower level/faster routines.
!
!

!--------------------------------------------    
!    implicit none
!    integer nM, NC, NCz
!    real(dp) rCM(3,maxCoo/3), a(3), a2(3)
!    real(dp) dpole(3,maxCoo/3), qpole(3,3,maxCoo/3)
!    real(dp) opole(3,3,3,maxCoo/3), hpole(3,3,3,3,maxCoo/3)
!    
!    real(dp) ed1(3,maxCoo/3), ed2(3,3,maxCoo/3), ed3(3,3,3,maxCoo/3)
!    real(dp) ed4(3,3,3,3,maxCoo/3), ed5(3,3,3,3,3,maxCoo/3), rMax2
!    
!    real(dp) d1d(3), d2d(3,3), d3d(3,3,3)
!    real(dp) d4d(3,3,3,3), d5d(3,3,3,3,3)
!    
!    real(dp) d1a(3), d2a(3,3), d3a(3,3,3)
!    real(dp) d4a(3,3,3,3), d5a(3,3,3,3,3)
!    
!    real(dp) d(3), q(3,3), o(3,3,3), h(3,3,3,3), fsf(3,maxCoo/3)
!    
!    integer i, j, k, l, s, n, m, ii, nx, ny, nz
!    integer in2(2), in3(3), in4(4), in5(5)
!    real(dp) re(3), dr(3), r1, r2, swFunc, dSdr
!    logical*1 iSlab
!--------------------------


    
    !timming
    integer*8 ti, tf, irtc
    real(dp) t1, t2, t3, t4, t5, t6, t7
    !timming
    
    !do n = 1, nM
    !   do i = 1, 3
    !      ed1(i,n) = 0.0_dp
    !      fsf(i,n) = 0.0_dp
    !      do j = i, 3
    !         ed2(i,j,n) = 0.0_dp
    !         do k = j, 3
    !            ed3(i,j,k,n) = 0.0_dp
    !            do l = k, 3
    !               ed4(i,j,k,l,n) = 0.0_dp
    !               do s = l, 3
    !                  ed5(i,j,k,l,s,n) = 0.0_dp
    !               end do
    !            end do
    !         end do
    !      end do
    !   end do
    !end do
    ed1 = 0
    fsf = 0
    ed2 = 0
    ed3 = 0
    ed4 = 0
    ed5 = 0
!  print*, "some tensor in calcDv"
!  print*, 0,ed1(:,2)
    
    !$$$      t1 = 0.0_dp
    !$$$      t2 = 0.0_dp
    !$$$      t3 = 0.0_dp
    !$$$      t4 = 0.0_dp
    !$$$      ti = irtc()
    
    NCz = NC
    if (iSlab) NCz = 0
    
    !!$omp parallel &
    !!$omp default(none) &
    !!$omp private(d1d,d2d,d3d,d4d,d5d,  d1a,d2a,d3a,d4a,d5a, ed1,ed2,ed3,ed4,ed5,  re, dr, r1, r2, &
    !!$omp         swFunc, dSdr, i, j, k, l, s, m, n, ii, nx, ny, nz, in2,in3,in4,in5  ) &   !,d,q,o,h
    !!$omp shared(nM, NC, a, NCz, rCM, a2, rMax2,dpole,qpole,opole,hpole,iSlab,fsf)
! d,q,o,h
!    !$omp default(private) &

    
    !!$omp do
    !JÖ loop over all water pairs:
    do n = 1, nM
       do m = 1, nM
          
          do nx = -NC, NC
             re(1) = a(1) * nx
             do ny = -NC, NC
                re(2) = a(2) * ny
                crescent:do nz = -NCz, NCz
                   re(3) = a(3) * nz
                   
                   !JÖ is same particle in same box, skip calculation:
                   if ( (n.eq.m) .and. (nx.eq.0) .and. (ny.eq.0) .and. (nz.eq.0)) cycle crescent!goto 11
                   
                   do i = 1, 3
                      dr(i) = rCM(i,n) - rCM(i,m)
                      if      (dr(i) .gt.  a2(i)) then; dr(i) = dr(i) - a(i)
                      else if (dr(i) .lt. -a2(i)) then; dr(i) = dr(i) + a(i)
                      end if !JÖ structure
                      dr(i) = dr(i) + re(i)
                   end do
                   
                   r2 = sum(dr**2) !JÖ dr(1)**2 + dr(2)**2 + dr(3)**2
                   
                   !JÖ Also if the distance is larger than the cutoff, skip calculation:
                   if (r2 .gt. rMax2) cycle crescent!JÖ goto 11
                   
                   r1 = dsqrt(r2)
                   call SFdsf(r1, swFunc, dSdr)
                   
                   
                   !do i = 1, 3
                   !   d(i) = dpole(i,m)
                   !end do
                   !                     call dDpole(d, dr, d1d, d2d, d3d, d4d, d5d)
                   !d = dpole(:,m)
                   !call dDpole(d, dr, d1a, d2a, d3a, d4a, d5a)
                   call dDpole(dpole(:,m), dr, d1a, d2a, d3a, d4a, d5a)
                   
                   !d5a = 0
                   
                   !$$$                     print *, n, m
                   !$$$                     print *, rCM(2,n), rCM(2,m), a(2)
                   !$$$
                   !$$$                     print *, d(1), d(2), d(3)
                   !$$$                     print *, dr(1), dr(2), dr(3)
                   !$$$                     print *, d1a(1), d1a(2), d1a(3)
                   !$$$                     stop
                   
                   !do j = 1, 3
                   !   do i = 1, 3
                   !      q(i,j) = qpole(i,j,m)
                   !   end do
                   !end do
                   
                   !q = qpole(:,:,m)
                   !call dQpole(q, dr, d1d, d2d, d3d, d4d, d5d)
                   call dQpole(qpole(:,:,m), dr, d1d, d2d, d3d, d4d, d5d)
                   
                   d1a = d1a + d1d
                   d2a = d2a + d2d
                   d3a = d3a + d3d
                   if(full)then
                     d4a = d4a + d4d
                     d5a = d5a + d5d
                   endif
                   
                   
                   
                   !call addDerivA(d1a, d2a, d3a, d4a, d5a, d1d, d2d, d3d, d4d, d5d)
                   
                   !$$$                     print *, d1a(1), d1a(2), d1a(3)
                   !$$$                     print *, d1d(1), d1d(2), d1d(3)
                   
                   !do k = 1, 3
                   !   do j = 1, 3
                   !      do i = 1, 3
                   !         o(i,j,k) = opole(i,j,k,m)
                   !      end do
                   !   end do
                   !end do
                   !o = opole(:,:,:,m)
                   !call dOpole(o, dr, d1d, d2d, d3d, d4d, d5d)
                   call dOpole(opole(:,:,:,m), dr, d1d, d2d, d3d, d4d, d5d)
                   !call addDerivA(d1a, d2a, d3a, d4a, d5a, d1d, d2d, d3d, d4d, d5d)
                   d1a = d1a + d1d
                   d2a = d2a + d2d
                   if(full)then
                     d3a = d3a + d3d
                     d4a = d4a + d4d
                     d5a = d5a + d5d
                   endif
                   
                   !$$$                     print *, d1a(1), d1a(2), d1a(3)
                   !$$$                     print *, d1d(1), d1d(2), d1d(3)
                   
                   !do l = 1, 3
                   !   do k = 1, 3
                   !      do j = 1, 3
                   !         do i = 1, 3
                   !            h(i,j,k,l) = hpole(i,j,k,l,m)
                   !         end do
                   !      end do
                   !   end do
                   !end do
                   !h = hpole(:,:,:,:,m)
                   !call dHpole(h, dr, d1d, d2d, d3d, d4d, d5d)
                   call dHpole(hpole(:,:,:,:,m), dr, d1d, d2d, d3d, d4d, d5d)
                   !call addDerivA(d1a, d2a, d3a, d4a, d5a, d1d, d2d, d3d, d4d, d5d)
                   d1a = d1a + d1d
                   if(full)then
                     d2a = d2a + d2d
                     d3a = d3a + d3d
                     d4a = d4a + d4d
                     d5a = d5a + d5d
                   endif
                   
                   call addDeriv(ed1, ed2, ed3, ed4, ed5, d1a, d2a, d3a, d4a, d5a, n, swFunc)
                   
                   !$$$                     print *, d1a(1), d1a(2), d1a(3)
                   !$$$                     print *, ed1(1,1), ed1(2,1), ed1(3,1)
                   !$$$                     stop
                   
                   
                   call addSwitchingForce(d1a, d2a, d3a, d4a, n, dSdr, dr, r1, dpole, qpole,opole, hpole, fsf)
                   
!11              end do
                end do crescent
             end do
          end do
       end do
    end do
    !!$omp end do 
    !!$omp end parallel
    
    !$$$      tf = irtc()
    !$$$      t4 = (tf-ti) * 1e-9
    !$$$
    !$$$      print '(A,f15.6)', 'Calculation: ', t4
    
    !$$$      print '(A,f15.6)', '      Dipoles: ', t1
    !$$$      print '(A,f15.6)', '  Quadrupoles: ', t2
    !$$$      print '(A,f15.6)', '    Octopoles: ', t3
    !$$$      print '(A,f15.6)', 'hexadecapoles: ', t4
    
    
    !     Copy all the permutations. (Is this really necessary??)
    !$$$      ti = irtc()
!  print*, 1,ed1(:,2)
    
    !!$omp parallel do &
    !!$omp private
    
    
    !!!!&
    !!!!!$omp default(private) &
    !!!!!$omp shared(ed2, ed3, ed4, ed5)
    
    
    !!$omp do
    do i = 1, 3
       do j = 1, 3
          !in2(1) = i
          !in2(2) = j
          in2 = [i,j]
          call insertIN(in2, 2)
          
          !do n = 1, nM
          !   ed2(i,j,n) = ed2(in2(1), in2(2), n)
          !end do
          ed2(i,j,:) = ed2(in2(1), in2(2), :)
          
          do k = 1, 3
             !do ii = 1, 2
             !   in3(ii) = in2(ii)
             !end do
             !in3(1:2) = in2
             !in3(3) = k
             in3 = [in2, k]
             call insertIN(in3, 3)
             
             !do n = 1, nM
             !   ed3(i,j,k,n) = ed3(in3(1),in3(2),in3(3),n)
             !end do
             ed3(i,j,k,:) = ed3(in3(1),in3(2),in3(3),:)
             
             do l = 1, 3
                !do ii = 1, 3
                !   in4(ii) = in3(ii)
                !end do
                !in4(1:3) = in3
                !in4(4) = l
                in4 = [in3, l]
                call insertIN(in4, 4)
                
                !do n = 1, nM
                !   ed4(i,j,k,l,n) = ed4(in4(1),in4(2),in4(3),in4(4),n)
                !end do
                ed4(i,j,k,l,:) = ed4(in4(1),in4(2),in4(3),in4(4),:)
                
                do m = 1, 3
                   !do ii = 1, 4
                   !   in5(ii) = in4(ii)
                   !end do
                   !in5(1:4) = in4
                   !in5(5) = m
                   in5 = [in4, m]
                   
                   call insertIN(in5, 5)
                   
                   !do n = 1, nM
                   !   ed5(i,j,k,l,m,n) = ed5(in5(1),in5(2),in5(3) ,in5(4),in5(5),n)
                   !end do
                   ed5(i,j,k,l,m,:) = ed5(in5(1),in5(2),in5(3) ,in5(4),in5(5),:)
                   
                end do
             end do
          end do
       end do
    end do
    !!$omp end do
    !!$omp end parallel 
!  print*, 2,ed1(:,2)
    
    !$$$      tf = irtc()
    !$$$      t4 = (tf-ti) * 1e-9
    !$$$      print '(A,f15.6)', 'Permutations: ', t4
    !$$$      stop
    
    return
    
  end subroutine calcDv
  
  
  !-----------------------------------------------------------------------
  pure subroutine insertIN(i, n) !JÖ pure
  !JÖ this SR removes the n:th value from i(:) moves a few intries up one index to cover the gap
  !and then puts the originally removed value back in the new void. This new place is the entry not larger than 
  !the original value, or the first entry in the vector i(:). 
    implicit none
    integer, intent(inout) :: i(*)
    integer, intent(in)    :: n
    integer ::  j, iaux !JÖ k, , n, j, , save
    
    if (n .ge. 2) then
       iaux = i(n)
       j = n-1
       do while (i(j) .gt. iaux .and. j .ge. 1)
          i(j+1) = i(j)
          j = j - 1
       end do
       i(j+1) = iaux
    end if
    
    return
    
  end subroutine insertIN
  
  !-----------------------------------------------------------------------
  pure integer function delta(i, j) !result(delta)
    implicit none
    integer, intent(in) :: i, j
    !integer delta !, intent(out) :: 
    
    if (i.eq.j) then
       delta = 0
    else
       delta = 0
    end if
    
    return
    
  end function delta
  
  !-----------------------------------------------------------------------
  pure subroutine addDeriv(ed1, ed2, ed3, ed4, ed5, d1d, d2d, d3d, d4d, d5d, n, swFunc)
    
    implicit none

!JÖ inout:    
    real(dp), intent(inout) :: ed1(:,:), ed2(:,:,:), ed3(:,:,:,:)
    real(dp), intent(inout) :: ed4(:,:,:,:,:), ed5(:,:,:,:,:,:)
    
    real(dp), intent(in) :: d1d(:), d2d(:,:), d3d(:,:,:)
    real(dp), intent(in) :: d4d(:,:,:,:), d5d(:,:,:,:,:)
    real(dp), intent(in) :: swFunc
    integer, intent(in) :: n
    
!JÖ internal:     
    integer :: i, j, k, l, s !JÖ , save

!    real(dp) ed1(3,maxCoo/3), ed2(3,3,maxCoo/3), ed3(3,3,3,maxCoo/3)
!    real(dp) ed4(3,3,3,3,maxCoo/3), ed5(3,3,3,3,3,maxCoo/3)
!    
!    real(dp) d1d(3), d2d(3,3), d3d(3,3,3)
!    real(dp) d4d(3,3,3,3), d5d(3,3,3,3,3)
!    real(dp) swFunc
!    
!    integer n, i, j, k, l, s

    
    do i = 1, 3
       ed1(i,n) = ed1(i,n) + d1d(i) * swFunc
       do j = i, 3
          ed2(i,j,n) = ed2(i,j,n) + d2d(i,j) * swFunc
          do k = j, 3
             ed3(i,j,k,n) = ed3(i,j,k,n) + d3d(i,j,k) * swFunc
             do l = k, 3
                ed4(i,j,k,l,n) = ed4(i,j,k,l,n) + d4d(i,j,k,l) * swFunc
                do s = l, 3
                   ed5(i,j,k,l,s,n) = ed5(i,j,k,l,s,n) + d5d(i,j,k,l,s ) * swFunc
                end do
             end do
          end do
       end do
    end do
    
    !ed1(:,n)          = ed1(:,n)         + d1d * swFunc
    !ed2(:,:,n)        = ed2(:,:,n)       + d2d * swFunc
    !ed3(:,:,:,n)      = ed3(:,:,:,n)     + d3d * swFunc
    !ed4(:,:,:,:,n)    = ed4(:,:,:,:,n)   + d4d * swFunc
    !ed5(:,:,:,:,:,n)  = ed5(:,:,:,:,:,n) + d5d * swFunc
    
    return
    
  end subroutine addDeriv
  
  !-----------------------------------------------------------------------
  pure subroutine addDerivA(d1a, d2a, d3a, d4a, d5a, d1d, d2d, d3d, d4d, d5d)
    
    implicit none
!JÖ in/Out    
    real(dp), intent(inout) :: d1a(:), d2a(:,:), d3a(:,:,:)
    real(dp), intent(inout) :: d4a(:,:,:,:), d5a(:,:,:,:,:)
    
    real(dp), intent(in) :: d1d(:), d2d(:,:), d3d(:,:,:)
    real(dp), intent(in) :: d4d(:,:,:,:), d5d(:,:,:,:,:)
!JÖ internal    
    integer :: n, i, j, k, l, s !JÖ , save

!    real(dp) d1a(3), d2a(3,3), d3a(3,3,3)
!    real(dp) d4a(3,3,3,3), d5a(3,3,3,3,3)
!    
!    real(dp) d1d(3), d2d(3,3), d3d(3,3,3)
!    real(dp) d4d(3,3,3,3), d5d(3,3,3,3,3)
!    
!    integer n, i, j, k, l, s

    
    do i = 1, 3
       d1a(i) = d1a(i) + d1d(i)
       do j = i, 3
          d2a(i,j) = d2a(i,j) + d2d(i,j)
          do k = j, 3
             d3a(i,j,k) = d3a(i,j,k) + d3d(i,j,k)
             do l = k, 3
                d4a(i,j,k,l) = d4a(i,j,k,l) + d4d(i,j,k,l)
                do s = l, 3
                   d5a(i,j,k,l,s) = d5a(i,j,k,l,s) + d5d(i,j,k,l,s)
                end do
             end do
          end do
       end do
    end do

    !d1a = d1a + d1d !* swFunc
    !d2a = d2a + d2d !* swFunc
    !d3a = d3a + d3d !* swFunc
    !d4a = d4a + d4d !* swFunc
    !d5a = d5a + d5d !* swFunc

    
    return
    
  end subroutine addDerivA
  
  !-----------------------------------------------------------------------
  pure subroutine addSwitchingForce(d1a, d2a, d3a, d4a, n, dSdr, dr, r1, dpole, qpole, opole, hpole, fsf)
    
    implicit none
    
    real(dp), intent(inout) :: d1a(:), d2a(:,:), d3a(:,:,:), d4a(:,:,:,:)
    real(dp), intent(in)    :: dSdr, dr(:), r1
    real(dp), intent(in)    :: dpole(:,:), qpole(:,:,:)
    real(dp), intent(in)    :: opole(:,:,:,:), hpole(:,:,:,:,:)
    real(dp), intent(out)   :: fsf(:,:)
    integer, intent(in) :: n
    
!JÖ internal:    
    real(dp) :: u
!JÖ , save
    integer :: i, ii, j, k, l, in2(2), in3(3), in4(4)
!JÖ , save

!    real(dp) d1a(3), d2a(3,3), d3a(3,3,3), d4a(3,3,3,3)
!    real(dp) dSdr, dr(3), r1, fsf(3,maxCoo/3), u
!    real(dp) dpole(3,maxCoo/3), qpole(3,3,maxCoo/3)
!    real(dp) opole(3,3,3,maxCoo/3), hpole(3,3,3,3,maxCoo/3)
!    
!    integer i, ii, j, k, l, in2(2), in3(3), in4(4), n
    
    !     Copy all the permutations.
    do i = 1, 3
       do j = 1, 3
          in2(1) = i
          in2(2) = j
          !in2 = [i,j]
          call insertIN(in2, 2)
          
          d2a(i,j) = d2a(in2(1), in2(2))
          
          do k = 1, 3
             do ii = 1, 2
                in3(ii) = in2(ii)
             end do
             in3(3) = k
             !in3 = [in2, k]
             call insertIN(in3, 3)
             
             d3a(i,j,k) = d3a(in3(1),in3(2),in3(3))
             
             do l = 1, 3
                do ii = 1, 3
                   in4(ii) = in3(ii)
                end do
                in4(4) = l
                !in4 = [in3, l]
                call insertIN(in4, 4)
                
                d4a(i,j,k,l) = d4a(in4(1),in4(2),in4(3),in4(4))
                
             end do
          end do
       end do
    end do
    
    u = 0.0_dp
    do i = 1, 3
       u = u + d1a(i) * dpole(i,n)
       do j = 1, 3
          u = u + d2a(i,j) * qpole(i,j,n) / 3.0_dp
          do k = 1, 3
             u = u + d3a(i,j,k) * opole(i,j,k,n) / 15.0_dp
             do l = 1, 3
                u = u + d4a(i,j,k,l) * hpole(i,j,k,l,n) / 105.0_dp
             end do
          end do
       end do
    end do

!JÖ my version, its baed, or not, we'll see    u = 0.0_dp !JÖ is this chit correct, check it!
!JÖ my version, its baed, or not, we'll see    do i = 1, 3
!JÖ my version, its baed, or not, we'll see       u = u + d1a(i) * dpole(i,n)
!JÖ my version, its baed, or not, we'll see       do j = 1, 3
!JÖ my version, its baed, or not, we'll see          u = u + d2a(j,i) * qpole(j,i,n) / 3.0_dp
!JÖ my version, its baed, or not, we'll see          do k = 1, 3
!JÖ my version, its baed, or not, we'll see             u = u + d3a(k,j,i) * opole(k,j,i,n) / 15.0_dp
!JÖ my version, its baed, or not, we'll see             do l = 1, 3
!JÖ my version, its baed, or not, we'll see                u = u + d4a(l,k,j,i) * hpole(l,k,j,i,n) / 105.0_dp
!JÖ my version, its baed, or not, we'll see             end do
!JÖ my version, its baed, or not, we'll see          end do
!JÖ my version, its baed, or not, we'll see       end do
!JÖ my version, its baed, or not, we'll see    end do

    
    u = -u * dSdr / r1
    do i = 1, 3
       fsf(i,n) = fsf(i,n) + u * dr(i)
    end do
    
    return
    
  end subroutine addSwitchingForce
  
  !--------------------------- calcDdField ------------------------------+
  !----------------------------------------------------------------------+
  !     Derivatives of the dipole field                                  |
  !----------------------------------------------------------------------+
  pure subroutine dDpole(d, r, d1d, d2d, d3d, d4d, d5d)
    
    implicit none
    real(dp), intent(in)  :: d(:), r(:) !3
    real(dp), intent(out) :: d1d(:), d2d(:,:), d3d(:,:,:)
    real(dp), intent(out) :: d4d(:,:,:,:), d5d(:,:,:,:,:)

!JÖ internal:    
    !integer  :: dij, dik, dil, dis, djk, djl, djs, dkl, dks, dls
    integer  :: i, j, k, l, s
    real(dp) :: r2, r3, r5, r7, r9, r11, r13, rd
    real(dp) :: t1, t2, y1, y2, y3, z1, z2, z3, w1, w2, w3, w4

    !real(dp) :: dt1k, dt2k, dt1l, ddt2kl, dy2l, dt2l, dy3l, dddt2kls
    !real(dp) :: dt1s, ddy2ls, dz1s, ddt2ks, dy2s, ddt2ls, dt2s
    !real(dp) :: ddy3ls, dz2s, dy3s, dz3s

    
    r2 = r(1)**2 + r(2)**2 + r(3)**2 !sum(r**2) !JÖ 
    r3 = dsqrt(r2) * r2
    r5 = r3 * r2
    r7 = r5 * r2
    r9 = r7 * r2
    r11 = r9 * r2
    r13 = r11 * r2
    
    r3  =        1.0_dp / r3
    r5  =       -3.0_dp / r5
    r7  =       15.0_dp / r7
    r9  =     -105.0_dp / r9
    r11 =      945.0_dp / r11
    r13 =   -10395.0_dp / r13
    
    rd = r(1)*d(1) + r(2)*d(2) + r(3)*d(3) !sum(r*d) !JÖ 
    
    do i = 1, 3
       d1d(i) = d(i) * r3 + rd * r(i) * r5
       
       do j = i, 3
          !     2nd derivative
          t1 = d(i)*r(j) + d(j)*r(i) !+ rd * dij
          t2 = rd * r(i)*r(j)
          
          d2d(i,j) = t1 * r5 + t2 * r7
          
          do k = j, 3
             y2   = t1*r(k) + d(k) * r(i)*r(j)
             y3   = t2 * r(k)
             
             d3d(i,j,k) = y2*r7 + y3*r9
             
             do l = k, 3
                !
                z2 = y2*r(l) + d(l) * r(i)*r(j) * r(k)
                z3 = y3*r(l)
                
                d4d(i,j,k,l) =  z2*r9 + z3*r11
                
                do s = l, 3
                   
                   w3 = z2 * r(s) + d(s) * r(i)*r(j)*r(k) *r(l)
                   
                   w4 = z3 * r(s)
                   
                   d5d(i,j,k,l,s) = w3*r11 + w4 *r13
                   
                end do
             end do
          end do
       end do
    end do
    
    
    return
    
  end subroutine dDpole
  
  
  
  !----------------------------------------------------------------------+
  !     Derivatives of the quadrupole field                              |
  !----------------------------------------------------------------------+
  pure subroutine dQpole(q, r, d1q, d2q, d3q, d4q, d5q)
    
    implicit none
    real(dp), intent(in)  :: q(:,:), r(:)
    real(dp), intent(out) :: d1q(:), d2q(:,:), d3q(:,:,:)
    real(dp), intent(out) :: d4q(:,:,:,:), d5q(:,:,:,:,:)

!    real(dp) q(3,3), r(3), d1q(3), d2q(3,3), d3q(3,3,3)
!    real(dp) d4q(3,3,3,3), d5q(3,3,3,3,3)

    
    integer dij, dik, dil, dis, djk, djl, djs, dkl, dks, dls
    integer i, j, k, l, s
    real(dp) :: r2, r3, r5, r7, r9, r11, r13, r15, rrq
    real(dp) :: t1, t2, t3, y1, y2, y3, z1, z2, z3, z4, w1, w2, w3, w4
    real(dp) :: v(3)
           ! :: 
    real(dp) :: dt2k, dt3k, ddt2kl, dy1l, dt2l, ddt3kl, dy2l, dt3l, dy3l
    real(dp) :: dddt3kls, dt2s, ddt2ls, ddy2ls, ddt2ks, dy1s, dz2s
    real(dp) :: ddt3ks, dy2s, ddt3ls, dt3s, ddy3ls, dz3s, dy3s, dz4s
    
    r2 =  r(1)**2 + r(2)**2 + r(3)**2 !sum(r**2) !JÖ
    r3 = dsqrt(r2) * r2
    r5 = r3 * r2
    r7 = r5 * r2
    r9 = r7 * r2
    r11 = r9 * r2
    r13 = r11 * r2
    r15 = r13 * r2
    
    r5  =        1.0_dp / r5
    r7  =       -5.0_dp / r7
    r9  =       35.0_dp / r9
    r11 =     -315.0_dp / r11
    r13 =     3465.0_dp / r13
    r15 =   -45045.0_dp / r15
    
    rrq = 0.0_dp
    do j = 1, 3
       v(j) = 0.0_dp
       do i = 1, 3
          rrq    = rrq + q(i,j) * r(i) * r(j)
          v(j)   = v(j) + q(i,j) * r(i)
       end do
    end do
    
    do i = 1, 3
       d1q(i) = (2.0_dp * v(i)) * r5 + (rrq * r(i)) * r7
       
       do j = i, 3
          !dij = delta(i,j)
          
          !     2nd derivative
          t1 = 2.0_dp * q(i,j)
          t2 = 2.0_dp * (v(i)*r(j) + v(j)*r(i)) !+ rrq * dij
          t3 = rrq * r(i)*r(j)
          
          d2q(i,j) = t1 * r5 + t2 * r7 + t3 * r9
          
          do k = j, 3
             !dik = delta(i,k)
             !djk = delta(j,k)
             
             !     3rd derivative
             !dt2k = 2.0_dp * (v(i)*djk + v(j)*dik + v(k) * dij) + 2.0_dp &
             !     * (q(i,k)*r(j) + q(j,k)*r(i))
             dt2k =  2.0_dp * (q(i,k)*r(j) + q(j,k)*r(i))
             !dt3k = 2.0_dp * v(k)*r(i)*r(j) + rrq * (r(i)*djk + r(j)*dik)
             dt3k = 2.0_dp * v(k)*r(i)*r(j) !+ rrq * (r(i)*djk + r(j)*dik)
             y1 = t1*r(k) + dt2k
             y2 = t2*r(k) + dt3k
             y3 = t3*r(k)
             
             d3q(i,j,k) = y1*r7 + y2*r9 + y3*r11
             
             do l = k, 3
                !dil = delta(i,l)
                !djl = delta(j,l)
                !dkl = delta(k,l)
                
                !     4th derivative
                !ddt2kl = 2.0_dp * (q(i,l)*djk + q(j,l)*dik + q(k,l) &
                !     * dij) + 2.0_dp * (q(i,k)*djl + q(j,k)*dil)
                !dy1l = t1*dkl + ddt2kl
                
                !z1 = dy1l
                
                !dt2l = 2.0_dp *(v(i)*djl + v(j)*dil + v(l) * dij)+ 2.0_dp* (q(i,l)*r(j) + q(j,l)*r(i))
                dt2l = 2.0_dp* (q(i,l)*r(j) + q(j,l)*r(i))
                
                
                !ddt3kl = 2.0_dp * q(k,l)*r(i)*r(j) + 2.0_dp * v(k)* (r(i)*djl + dil*r(j) ) &
                !     + 2.0_dp * v(l) * (r(i)*djk + r(j)*dik) &
                !     + rrq * (dil*djk + djl*dik)
                ddt3kl = 2.0_dp * q(k,l)*r(i)*r(j) 
                
                !dy2l = dt2l*r(k) + t2*dkl + ddt3kl
                dy2l = dt2l*r(k) + ddt3kl
                
                z2 = y1*r(l) + dy2l
                
                
                !dt3l = 2.0_dp * v(l)*r(i)*r(j) + rrq * (r(i)*djl + r(j)*dil)
                dt3l = 2.0_dp * v(l)*r(i)*r(j) 
                
                !dy3l = dt3l*r(k) + t3 * dkl
                dy3l = dt3l*r(k) 
                z3 = y2*r(l) + dy3l
                
                z4 = y3*r(l)
                d4q(i,j,k,l) =  z1*r7 + z2*r9 + z3*r11 + z4*r13
                
                do s = l, 3
                   !dis = delta(i,s)
                   !djs = delta(j,s)
                   !dks = delta(k,s)
                   !dls = delta(l,s)
                   
                   
                   !     5th derivative
                   
                   !dddt3kls = 2.0_dp * q(k,l)*(r(i)*djs + dis*r(j)) &
                   !     + 2.0_dp * q(k,s) * (r(i)*djl + dil*r(j)) &
                   !     + 2.0_dp * v(k) * (dis*djl + dil*djs) &
                   !     + 2.0_dp * q(l,s) * (r(i)*djk + r(j)*dik) &
                   !     + 2.0_dp * v(l) * (dis*djk + djs*dik) &
                   !     + 2.0_dp * v(s) * (dil*djk + djl*dik)
                   
                   !dt2s = 2.0_dp * (v(i)*djs + v(j)*dis + v(s) * dij) + &
                   !     2.0_dp* (q(i,s)*r(j) + q(j,s)*r(i))
                   dt2s = 2.0_dp* (q(i,s)*r(j) + q(j,s)*r(i))
                   
                   !ddt2ls = 2.0_dp * (q(i,s)*djl + q(j,s)*dil + q(l,s) &
                   !     * dij) + 2.0_dp* (q(i,l)*djs + q(j,l)*dis)
                   
                   !ddy2ls = ddt2ls*r(k) + dt2l*dks + dt2s*dkl + &
                   !     dddt3kls
                   
                   !ddt2ks = 2.0_dp * (q(i,s)*djk + q(j,s)*dik + q(k,s) &
                   !     * dij) + 2.0_dp* (q(i,k)*djs + q(j,k)*dis)
                   
                   !dy1s = t1*dks + ddt2ks
                   
                   !dz2s = y1*dls + dy1s*r(l) + ddy2ls
                   !w1 = z1*r(s) + dz2s
                   !----------
                   !ddt3ks = 2.0_dp * q(k,s)*r(i)*r(j) + 2.0_dp * v(k) &
                   !     * (r(i)*djs + dis*r(j) )+ 2.0_dp * v(s) * (r(i) &
                   !     * djk + r(j)*dik)+ rrq * (dis*djk + djs*dik)
                   ddt3ks = 2.0_dp * q(k,s)*r(i)*r(j) 
                   
                   !dy2s = dt2s*r(k) + t2*dks + ddt3ks
                   dy2s = dt2s*r(k) + ddt3ks
                   
                   
                   !ddt3ls = 2.0_dp * (q(l,s)*r(i)*r(j) + v(l)*(dis*r(j)+ r(i)*djs) ) &
                   !     + 2.0_dp * v(s) * (r(i)*djl + r(j)*dil) &
                   !     + rrq * (dis*djl + djs*dil)
                   ddt3ls = 2.0_dp * q(l,s)*r(i)*r(j)
                   
                   !dt3s = 2.0_dp * v(s)*r(i)*r(j) + rrq * (r(i)*djs +r(j)*dis)
                   dt3s = 2.0_dp * v(s)*r(i)*r(j) 
                   !ddy3ls = ddt3ls*r(k) + dt3l*dks + dt3s * dkl
                   ddy3ls = ddt3ls*r(k)
                   
                   !dz3s = dy2s*r(l) + y2*dls + ddy3ls
                   dz3s = dy2s*r(l) + ddy3ls
                   w2 = z2*r(s) + dz3s
                   
                   !------------
                   !                     dt3s =
                   !dy3s = dt3s*r(k) + t3 * dks
                   dy3s = dt3s*r(k) 
                   
                   !dz4s = dy3s*r(l) + y3*dls
                   dz4s = dy3s*r(l)
                   w3 = z3*r(s) + dz4s
                   
                   !-------------
                   w4 = z4*r(s)
                   
                   !d5q(i,j,k,l,s) = w1*r9 + w2*r11 + w3*r13 + w4*r15
                   d5q(i,j,k,l,s) = w2*r11 + w3*r13 + w4*r15
                   
                end do
             end do
          end do
       end do
    end do
    
    
    return
    
  end subroutine dQpole
  
  
  
  
  !----------------------------------------------------------------------+
  !     Derivatives of the octopole field                                |
  !----------------------------------------------------------------------+
  pure subroutine dOpole(o, r, d1o, d2o, d3o, d4o, d5o)
    
    implicit none
    real(dp), intent(in)  :: o(:,:,:), r(:)
    real(dp), intent(out) :: d1o(:), d2o(:,:), d3o(:,:,:)
    real(dp), intent(out) :: d4o(:,:,:,:), d5o(:,:,:,:,:)

!    real(dp) o(3,3,3), r(3), d1o(3), d2o(3,3), d3o(3,3,3)
!    real(dp) d4o(3,3,3,3), d5o(3,3,3,3,3)
    
    integer dij, dik, dil, dis, djk, djl, djs, dkl, dks, dls
    integer i, j, k, l, s
    real(dp) r2, r3, r5, r7, r9, r11, r13, r15, r17
    real(dp) t1, t2, t3, y1, y2, y3, y4, z1, z2, z3, z4, w1, w2, w3, w4,w5
    real(dp) r3o, v(3), g(3,3)
           !
    real(dp) dt1k, dt2k, dt3k
    real(dp) dt1l, ddt2kl, dy2l, dt2l, ddt3kl, dy3l, dt3l, dy4l
    real(dp) dt1s, dddt2kls, ddy2ls, dz1s, ddt2ks, dy2s, ddt2ls, dt2s, &
         dddt3kls, ddy3ls, dz2s, ddt3ks, dy3s, dt3s, ddt3ls, ddy4ls, &
         dz3s, dy4s, dz4s
    
    
    r2 = r(1)**2 + r(2)**2 + r(3)**2 !sum(r**2) !JÖ 
    r3 = dsqrt(r2) * r2
    r5 = r3 * r2
    r7 = r5 * r2
    r9 = r7 * r2
    r11 = r9 * r2
    r13 = r11 * r2
    r15 = r13 * r2
    r17 = r15 * r2
    
    r7  =       1.0_dp / r7
    r9  =      -7.0_dp / r9
    r11 =      63.0_dp / r11
    r13 =    -693.0_dp / r13
    r15 =    9009.0_dp / r15
    r17 = -135135.0_dp / r17
    
    r3o = 0.0_dp
    do k = 1, 3
       v(k) = 0.0_dp
       do j = 1, 3
          g(j,k) = 0.0_dp
          do i = 1, 3
             r3o    = r3o + o(i,j,k) * r(i) * r(j) * r(k)
             v(k)   = v(k) + o(i,j,k) * r(i) * r(j)
             g(j,k) = g(j,k) + o(i,j,k) * r(i)
          end do
       end do
    end do
    
    do i = 1, 3
       d1o(i) = (3.0_dp * v(i)) * r7 + (r3o * r(i)) * r9
       
       do j = i, 3
          !dij = delta(i,j)
          
          !     2nd derivative
          t1 = 6.0_dp * g(i,j)
          !t2 = 3.0_dp * (v(i)*r(j) + v(j)*r(i)) + r3o * dij
          t2 = 3.0_dp * (v(i)*r(j) + v(j)*r(i))
          t3 = r3o * r(i)*r(j)
          
          d2o(i,j) = t1 * r7 + t2 * r9 + t3 * r11
          
          do k = j, 3
             !dik = delta(i,k)
             !djk = delta(j,k)
             
             !     3rd derivative
             dt1k = 6.0_dp * o(i,j,k)
             !dt2k = 3.0_dp * (v(i)*djk + v(j)*dik + v(k) * dij) & 
             !     + 6.0_dp * (g(i,k)*r(j) + g(j,k)*r(i))
             dt2k = 6.0_dp * (g(i,k)*r(j) + g(j,k)*r(i))
             !dt3k = 3.0_dp * v(k)*r(i)*r(j) + r3o * (r(i)*djk + r(j)*dik)
             dt3k = 3.0_dp * v(k)*r(i)*r(j)
             y1 = dt1k
             y2 = t1*r(k) + dt2k
             y3 = t2*r(k) + dt3k
             y4 = t3*r(k)
             
             d3o(i,j,k) = y1*r7 + y2*r9 + y3*r11 + y4*r13
             
             do l = k, 3
                !dil = delta(i,l)
                !djl = delta(j,l)
                !dkl = delta(k,l)
                
                !     4th derivative
                dt1l = 6.0_dp * o(i,j,l)
                !ddt2kl = 6.0_dp * (g(i,l)*djk + g(j,l)*dik + g(k,l) * dij) &
                !       + 6.0_dp * (o(i,k,l)*r(j) + o(j,k,l)*r(i)) &
                !       + 6.0_dp * (g(i,k)*djl + g(j,k)*dil)
                ddt2kl = 6.0_dp * (o(i,k,l)*r(j) + o(j,k,l)*r(i)) 
                !dy2l = dt1l*r(k) + t1*dkl + ddt2kl
                dy2l = dt1l*r(k) + ddt2kl
                z1 = y1*r(l) + dy2l
                
                !-------
                !dt2l = 3.0_dp * (v(i)*djl + v(j)*dil + v(l) * dij) &
                !     + 6.0_dp * (g(i,l)*r(j) + g(j,l)*r(i))
                dt2l = 6.0_dp * (g(i,l)*r(j) + g(j,l)*r(i))
                
                !ddt3kl = 3.0_dp * (2.0_dp * g(k,l)*r(i)*r(j)  &
                !     + v(k) * (dil*r(j) + r(i)*djl) )  &
                !     + r3o * (dil*djk + djl*dik)  &
                !     + 3.0_dp * v(l) * (r(i)*djk + r(j)*dik)
                ddt3kl = 3.0_dp * 2.0_dp * g(k,l)*r(i)*r(j)  
                !dy3l = dt2l*r(k) + t2*dkl + ddt3kl
                dy3l = dt2l*r(k) + ddt3kl
                z2 = y2*r(l) + dy3l
                
                !-------
                !dt3l = 3.0_dp * v(l)*r(i)*r(j) + r3o * (r(i)*djl + r(j)*dil)
                dt3l = 3.0_dp * v(l)*r(i)*r(j) 
                !dy4l = dt3l*r(k) + t3 * dkl
                dy4l = dt3l*r(k)
                z3 = y3*r(l) + dy4l
                
                z4 = y4*r(l)
                d4o(i,j,k,l) =  z1*r9 + z2*r11 + z3*r13 + z4*r15
                
                
                do s = l, 3
                          
                   dy2s = 6*( o(i,j,s)*r(k) + o(i,k,s)*r(j) + o(j,k,s)*r(i)) 
                   
                   dddt3kls = 6*o(k,l,s)*r(i)*r(j) 
                   
                   ddy3ls = 6*r(k)*(o(i,l,s)*r(j) + o(j,l,s)*r(i)) + dddt3kls
                   
                   dz2s = dy2s*r(l) + ddy3ls
                   
                   w2 = z1*r(s) + dz2s
                   dy3s = 6*(r(k)*(g(i,s)*r(j) + g(j,s)*r(i)) + g(k,s)*r(i)*r(j) )
                   
                   w3 = z2*r(s) + dy3s*r(l) + 6.0_dp * g(l,s)*r(i)*r(j) *r(k) 
                   w4 = z3*r(s) + 3.0_dp * v(s)*r(i)*r(j)*r(k)*r(l)
                   w5 = z4*r(s)
                   
                   d5o(i,j,k,l,s) = w1*r9 + w2*r11 + w3*r13 + w4*r15 + w5*r17
                   
                end do
             end do
          end do
       end do
    end do
    
    
    return
    
  end subroutine dOpole
  
  
  
  !----------------------------------------------------------------------+
  !     Derivatives of the hexadecapole field                            |
  !----------------------------------------------------------------------+
  pure subroutine dHpole(h, r, d1h, d2h, d3h, d4h, d5h)
    
    implicit none
    
    real(dp), intent(in)  :: h(:,:,:,:), r(:)
    real(dp), intent(out) :: d1h(:), d2h(:,:), d3h(:,:,:)
    real(dp), intent(out) :: d4h(:,:,:,:), d5h(:,:,:,:,:)
    
!    real(dp) h(3,3,3,3), r(3), d1h(3), d2h(3,3), d3h(3,3,3)
!    real(dp) d4h(3,3,3,3), d5h(3,3,3,3,3)
    
    integer  :: dij, dik, dil, dis, djk, djl, djs, dkl, dks, dls
    integer  :: i, j, k, l, s
    real(dp) :: r2, r3, r5, r7, r9, r11, r13, r15, r17, r19
    real(dp) :: t1, t2, t3, y1, y2, y3, y4, z1, z2, z3, z4, z5, w1, w2, w3,w4, w5
    real(dp) :: r4h, v(3), g(3,3), d(3,3,3)
           ! 
    real(dp) :: dt1k, dt2k, dt3k
    real(dp) :: dy1l, dt1l, ddt2kl, dy2l, dt2l, ddt3kl, dy3l, dt3l, dy4l
    real(dp) :: ddt1ls, dt1s, dddt2kls, ddy2ls, dy1s, dz2s, ddt2ks, dy2s, &
         dt2s, ddt2ls, dddt3kls, ddy3ls, dz3s, ddt3ks, dy3s, dt3s, &
         ddt3ls, ddy4ls, dz4s, dy4s, dz5s
    
    
    r2 = sum(r**2) !JÖ r(1)**2 + r(2)**2 + r(3)**2
    r3 = dsqrt(r2) * r2
    r5 = r3 * r2
    r7 = r5 * r2
    r9 = r7 * r2
    r11 = r9 * r2
    r13 = r11 * r2
    r15 = r13 * r2
    r17 = r15 * r2
    r19 = r17 * r2
    
    r9  =        1.0_dp / r9
    r11 =       -9.0_dp / r11
    r13 =       99.0_dp / r13
    r15 =    -1287.0_dp / r15
    r17 =    19305.0_dp / r17
    r19 =  -328185.0_dp / r19
    
    r4h = 0.0_dp
    do l = 1, 3
       v(l) = 0.0_dp
       do k = 1, 3
          g(l,k) = 0.0_dp
          do j = 1, 3
             d(l,k,j) = 0.0_dp
             do i = 1, 3
                r4h    = r4h + h(i,j,k,l) * r(i) * r(j) * r(k) * r(l)
                v(l)   = v(l) + h(i,j,k,l) * r(i) * r(j) * r(k)
                g(l,k) = g(l,k) + h(i,j,k,l) * r(i) * r(j)
                d(l,k,j) = d(l,k,j) + h(i,j,k,l) * r(i)
             end do
          end do
       end do
    end do
    
    
    do i = 1, 3
       d1h(i) = (4.0_dp * v(i)) * r9 + (r4h * r(i)) * r11
       
       do j = i, 3
          
          !     2nd derivative
          t1 = 12.0_dp * g(i,j)
          t2 = 4.0_dp * (v(i)*r(j) + v(j)*r(i)) 
          t3 = r4h * r(i)*r(j)
          
          d2h(i,j) = t1 * r9 + t2 * r11 + t3 * r13
          
          do k = j, 3
             
             !     3rd derivative
             y1 = 24.0_dp * d(i,j,k)
             y2 = t1*r(k) + 12.0_dp * (g(i,k)*r(j) + g(j,k)*r(i))
             y3 = t2*r(k) + 4.0_dp * v(k)*r(i)*r(j)
             y4 = t3*r(k)
             
             d3h(i,j,k) = y1*r9 + y2*r11 + y3*r13 + y4*r15
             
             do l = k, 3
                
                !     4th derivative
                z1 = 24.0_dp * h(i,j,k,l)
                z2 = y1*r(l) + 24.0_dp * d(i,j,l)*r(k) + 24.0_dp * (d(i,k,l)*r(j) + d(j,k,l)*r(i))
                z3 = y2*r(l) + 12.0_dp*r(k)*(g(i,l)*r(j) + g(j,l)*r(i)) + 4.0_dp * 3.0_dp * g(k,l)*r(i)*r(j) 
                z4 = y3*r(l) + 4.0_dp * v(l)*r(i)*r(j)*r(k)
                z5 = y4*r(l)
                
                d4h(i,j,k,l) =  z1*r9 + z2*r11 + z3*r13 + z4*r15 +z5*r17
                
                do s = l, 3
                   
                   w1 = z1*r(s) + 24*( h(i,j,k,s)*r(l) + h(i,j,l,s)*r(k) + h(i,k,l,s)*r(j) + h(j,k,l,s)*r(i) )
                   w2 = z2*r(s) + 24*r(l)*( d(i,j,s)*r(k) + d(i,k,s)*r(j) + d(j,k,s)*r(i)) &
					            + 24*(r(k)*( d(i,l,s)*r(j) + d(j,l,s)*r(i)) + d(k,l,s)*r(i)*r(j) )
                   w3 = z3*r(s) + 12 * r(l)*( r(k)*(g(i,s)*r(j) + g(j,s)*r(i))  + g(k,s)*r(i)*r(j)) &
                                + 12 * g(l,s)*r(i)*r(j)*r(k)
                   w4 = z4*r(s) + 4 * v(s)*r(i)*r(j)*r(k)*r(l)
                   
                   w5 = z5*r(s)
                   
                   d5h(i,j,k,l,s) = w1*r11 + w2*r13 + w3*r15 + w4*r17+w5*r19
                   
                end do
             end do
          end do
       end do
    end do
    
    
    return
    
  end subroutine dHpole
  
end module calc_derivs

