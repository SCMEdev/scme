module molforce
  
  use data_types
  use mdutil, only: inv6
  
  implicit none
  
  private
  public molForce3
  
contains
  
  !-----------------------------------------------------------------------
  subroutine molForce3(f1, f2, f3, r1, r2, r3, fTot, tt, flag)
    
    implicit none
    
    integer i, j, ii, jj, ll
    real(dp) r1(3), r2(3), r3(3), rc1(3), rc2(3), rc3(3)
    real(dp) f1(3), f2(3), f3(3), tor(3), fidd(3)
    
    real(dp) b(6), c(6), y(6,6)
    real(dp) tt(3), ftot(3), t1(3), t2(3), t3(3)
    integer flag
    real(dp) f1a(3), f2a(3), f3a(3)
    !      data mF / 0.352d0, -0.864d0, -0.36d0, 0.36d0, 0.48d0, -0.8d0,
    !     $     .864d0, 0.152d0, 0.48d0 /
    !      data mB / 0.352d0, 0.36d0, 0.864d0, -0.864d0, 0.48d0, 0.152d0,
    !     $     -0.36d0, -0.8d0, 0.48d0/
    
    real(dp) mF(3,3,3), mB(3,3,3), sq2, oh
    parameter (sq2=0.707106781186547524d0, oh=0.5d0)
    !$$$      data mF / -sq2, oh, -oh, 0.d0, sq2, sq2, sq2, oh, -oh,
    !$$$     $     -oh, -sq2, -oh, oh, -sq2, oh, -sq2, 0.d0, sq2,
    !$$$     $     sq2, sq2, 0.d0, oh, -oh, -sq2, -oh, oh, -sq2 /
    !$$$
    !$$$      data mB / -sq2, 0.d0, sq2, oh, sq2, oh, -oh, sq2, -oh,
    !$$$     $     -oh, oh, -sq2, -sq2, -sq2, 0.d0, -oh, oh, sq2,
    !$$$     $     sq2, oh, -oh, sq2, -oh, oh, 0.d0, -sq2, -sq2 /
    
    data mF / -0.707106781186547524d0, 0.5d0, -0.5d0, 0.d0, 0.707106781186547524d0, &
         0.707106781186547524d0, 0.707106781186547524d0, 0.5d0, -0.5d0,-0.5d0, &
         -0.707106781186547524d0, -0.5d0, 0.5d0, -0.707106781186547524d0, &
         0.5d0, -0.707106781186547524d0, 0.d0, 0.707106781186547524d0, &
         0.707106781186547524d0, 0.707106781186547524d0, 0.d0, 0.5d0, -0.5d0, &
         -0.707106781186547524d0, -0.5d0, 0.5d0, -0.707106781186547524d0 /
    
    data mB / -0.707106781186547524d0, 0.d0, 0.707106781186547524d0, &
         0.5d0, 0.707106781186547524d0, 0.5d0, -0.5d0, 0.707106781186547524d0, &
         -0.5d0,-0.5d0, 0.5d0, -0.707106781186547524d0, -0.707106781186547524d0, &
         -0.707106781186547524d0, 0.d0, -0.5d0, 0.5d0, 0.707106781186547524d0, &
         0.707106781186547524d0, 0.5d0, -0.5d0,0.707106781186547524d0, -0.5d0, 0.5d0, &
         0.d0, -0.707106781186547524d0, -0.707106781186547524d0 /
    
    call inv6(r1, r2, r3, y, flag)
    
    !$$$c Total force
    
    !$$$      if (flag .gt. 0) then
    !$$$      do i = 1, 3
    !$$$         fTot(i) = f1(i) + f2(i) + f3(i) + fidd(i)
    !$$$      end do
    !$$$      call cross (rc1, f1, t1)
    !$$$      call cross (rc2, f2, t2)
    !$$$      call cross (rc3, f3, t3)
    !$$$c Total torque
    !$$$      do i = 1, 3
    !$$$         tt(i) = t1(i) + t2(i) + t3(i) + tor(i)
    !$$$      end do
    !$$$
    !$$$c----- Start Debugging ----
    !$$$c      print *, 'NEW --------'
    !$$$      print *
    !$$$      print '(6f19.14)', tt(1), tt(2), tt(3), ftot(1), ftot(2),
    !$$$     $     ftot(3)
    !$$$      print '(9g15.5)', (f1(ll), ll=1,3), (f2(ll), ll=1,3),
    !$$$     $     (f3(ll), ll=1,3)
    !$$$      end if
    !----- End Debugging ----
    
    !     The determinant was zero, therefore, we rotated the space. We have
    !     to rotate also the force and torque and then bring things back to
    !     the original orientation.
    if (flag .gt. 0) then
       do ii = 1, 3
          f1a(ii) = 0.d0
          f2a(ii) = 0.d0
          do jj = 1, 3
             f1a(ii) = f1a(ii) + mF(ii,jj,flag) * fTot(jj)
             f2a(ii) = f2a(ii) + mF(ii,jj,flag) * tt(jj)
          end do
       end do
       do ii = 1, 3
          ftot(ii) = f1a(ii)
          tt(ii) = f2a(ii)
       end do
    end if
    
    do i = 1, 3
       b(i) = fTot(i)
       b(i+3) = tt(i)
    end do
    
    do i = 1, 6
       c(i) = 0.d0
       do j = 1, 6
          c(i) = c(i) + y(i,j) * b(j)
       end do
    end do
    
    !     New forces...
    
    f1(1) = c(1)
    f1(2) = 0.d0
    f1(3) = c(3)
    
    f2(1) = 0.d0
    f2(2) = c(2)
    f2(3) = c(6)
    
    f3(1) = c(4)
    f3(2) = c(5)
    f3(3) = 0.d0
    
    !     Here we bring the forces to the original orientation of the space.
    if (flag .gt. 0) then
       do ii = 1, 3
          f1a(ii) = 0.d0
          f2a(ii) = 0.d0
          f3a(ii) = 0.d0
          do jj = 1, 3
             f1a(ii) = f1a(ii) + mB(ii,jj,flag) * f1(jj)
             f2a(ii) = f2a(ii) + mB(ii,jj,flag) * f2(jj)
             f3a(ii) = f3a(ii) + mB(ii,jj,flag) * f3(jj)
          end do
       end do
       do ii = 1, 3
          f1(ii) = f1a(ii)
          f2(ii) = f2a(ii)
          f3(ii) = f3a(ii)
       end do
    end if
    
    !----- Start Debugging ----
    !$$$      if (flag .gt. 0) then
    !$$$      call cross (r1, f1, t1)
    !$$$      call cross (r2, f2, t2)
    !$$$      call cross (r3, f3, t3)
    !$$$c Total torque
    !$$$      do i = 1, 3
    !$$$         tt(i) = t1(i) + t2(i) + t3(i)
    !$$$      end do
    !$$$c Total force
    !$$$      do i = 1, 3
    !$$$         ftot(i) = f1(i) + f2(i) + f3(i)
    !$$$      end do
    !$$$      print '(6f19.10)', tt(1), tt(2), tt(3), ftot(1), ftot(2),
    !$$$     $     ftot(3)
    !$$$c      print '(9g15.5)', (f1(ll), ll=1,3), (f2(ll), ll=1,3),
    !$$$c     $     (f3(ll), ll=1,3)
    !$$$      print *, '------------'
    !$$$      end if
    !----- End Debugging ----
    
    return
    
  end subroutine molForce3
  
end module molforce
