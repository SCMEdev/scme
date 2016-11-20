module mdutil
  
  use data_types
  
  implicit none
  
  private
  public cross, dot, inv6
  
contains
  
  !-----------------------------------------------------------------------
  subroutine cross(r1, r2, r3)
    !     Routine to calculate cross products of 3-vectors. r3 = r1 x r2
    !       Used by tangentialVel - PRB
    
    implicit none
    real(dp) r1(3), r2(3), r3(3)
    
    r3(1) = r1(2)*r2(3) - r1(3)*r2(2)
    r3(2) = r1(3)*r2(1) - r1(1)*r2(3)
    r3(3) = r1(1)*r2(2) - r1(2)*r2(1)
    
    return
    
  end subroutine cross
  
  !-----------------------------------------------------------------------
  real(dp) function dot(r1, r2)
    !     Routine to calculate the dot product of 3-vectors. dot = r1.r2
    
    ! used by shaker, rigidatoms_mod, tangentialvel_mod
    
    real(dp) r1(3), r2(3)
    
    dot = r1(1)*r2(1) + r1(2)*r2(2) + r1(3)*r2(3)
    
    return
    
  end function dot
  
  !-----------------------------------------------------------------------
  subroutine inverseSym(x, y)
    !
    !     Subroutine to invert a 3x3 symmetrical matrix. Input:  x
    !                                                    Output: y
    !
    ! Note: not used -PRB
    
    real(dp) x(3,3), y(3,3), det
    integer :: i, j
    
    det = - x(1,3)*x(1,3) * x(2,2) + 2.0 * x(1,2) * x(1,3) * x(2,3)
    det = det - x(1,1) * x(2,3)*x(2,3) - x(1,2)*x(1,2) * x(3,3)
    det = det + x(1,1) * x(2,2) * x(3,3)
    
    if (abs(det) .lt. 1e-13) then
       !         print *, '2D'
       det = x(2,2)*x(3,3) - x(2,3)*x(2,3)
       y(1,1) = 0.0
       y(1,2) = 0.0
       y(1,3) = 0.0
       y(2,1) = 0.0
       y(3,1) = 0.0
       y(2,2) = x(3,3)
       y(2,3) = -x(2,3)
       y(3,2) = -x(2,3)
       y(3,3) = x(2,2)
       !         print '(A, $)', 'Error in inverse. Determinant ='
       !         print *, det, x(1,1)
       !         stop
    else
       !         print *, '3D'
       y(1,1) = -x(2,3) * x(2,3) + x(2,2) * x(3,3)
       y(1,2) =  x(1,3) * x(2,3) - x(1,2) * x(3,3)
       y(1,3) = -x(1,3) * x(2,2) + x(1,2) * x(2,3)
       y(2,2) = -x(1,3) * x(1,3) + x(1,1) * x(3,3)
       y(2,3) =  x(1,2) * x(1,3) - x(1,1) * x(2,3)
       y(3,3) = -x(1,2) * x(1,2) + x(1,1) * x(2,2)
       y(2,1) = y(1,2)
       y(3,1) = y(1,3)
       y(3,2) = y(2,3)
    end if
    
    do i = 1, 3
       do j = 1, 3
          y(i,j) = y(i,j) / det
       end do
    end do
    
  end subroutine inverseSym
  
  !-----------------------------------------------------------------------
  subroutine inverse(x, y)
    !
    !     Subroutine to invert a 3x3 matrix.
    !
    ! Note: not used -PRB
    
    integer :: i,j
    real(dp) x(3,3), y(3,3), det
    
    det = - x(1,3)*x(2,2)*x(3,1) + x(1,2)*x(2,3)*x(3,1)
    det = det + x(1,3)*x(2,1)*x(3,2) - x(1,1)*x(2,3)*x(3,2)
    det = det - x(1,2)*x(2,1)*x(3,3) + x(1,1)*x(2,2)*x(3,3)
    
    if ((det .lt. 1e-13) .and. (det .gt. -1e-13)) then
       print '(A)', ' Error in inverse(x,y). det = 0.'
       stop
    end if
    
    y(1,1) = -x(2,3) * x(3,2) + x(2,2) * x(3,3)
    y(1,2) =  x(1,3) * x(3,2) - x(1,2) * x(3,3)
    y(1,3) = -x(1,3) * x(2,2) + x(1,2) * x(2,3)
    
    y(2,1) =  x(2,3) * x(3,1) - x(2,1) * x(3,3)
    y(2,2) = -x(1,3) * x(3,1) + x(1,1) * x(3,3)
    y(2,3) =  x(2,1) * x(1,3) - x(1,1) * x(2,3)
    
    y(3,1) = -x(3,1) * x(2,2) + x(2,1) * x(3,2)
    y(3,2) =  x(1,2) * x(3,1) - x(1,1) * x(3,2)
    y(3,3) = -x(1,2) * x(2,1) + x(1,1) * x(2,2)
    
    do i = 1, 3
       do j = 1, 3
          y(i,j) = y(i,j) / det
       end do
    end do
    
  end subroutine inverse
  
  !-----------------------------------------------------------------------
  subroutine inv6(r1, r2, r3, y, flag)
    
    ! used my molforce
    
    implicit none
    integer ii, jj, iRot
    real(dp) a, b, c, d, e, f, g, h, i, y(6,6)
    real(dp) r1(3), r2(3), r3(3)
    real(dp) r1a(3,3), r2a(3,3), r3a(3,3), det, detOld
    integer flag
    
    real(dp) mF(3,3,3), sq2, oh
    parameter (sq2=0.707106781186547524d0, oh=0.5d0)
    !$$$      data mF / -sq2, oh, -oh, 0.d0, sq2, sq2, sq2, oh, -oh,
    !$$$     $     -oh, -sq2, -oh, oh, -sq2, oh, -sq2, 0.d0, sq2,
    !$$$     $     sq2, sq2, 0.d0, oh, -oh, -sq2, -oh, oh, -sq2 /
    
    data mF / -0.707106781186547524d0, 0.5d0, -0.5d0, 0.d0, 0.707106781186547524d0, &
         0.707106781186547524d0, 0.707106781186547524d0, 0.5d0, -0.5d0,-0.5d0, &
         -0.707106781186547524d0, -0.5d0, 0.5d0, -0.707106781186547524d0, 0.5d0, &
         -0.707106781186547524d0, 0.d0, 0.707106781186547524d0,0.707106781186547524d0, &
         0.707106781186547524d0, 0.d0, 0.5d0, -0.5d0, -0.707106781186547524d0, -0.5d0, &
         0.5d0, -0.707106781186547524d0 /
    
    
    !     Rot1 = Rx[45] Ry[-135]
    !     Rot2 = Ry[45] Rz[-135]
    !     Rot2 = Rz[45] Rx[-135]
    
    a = r1(1)
    b = r1(2)
    c = r1(3)
    
    d = r2(1)
    e = r2(2)
    f = r2(3)
    
    g = r3(1)
    h = r3(2)
    i = r3(3)
    
    det = -b*c*d + c*d*e + a*b*f - b*d*f + b*c*g - c*e*g - a*f*h + &
         d*f*h - a*b*i + 2*b*d*i - d*e*i - b*g*i + e*g*i + a*h*i - d*h*i
    
    flag = 0
    if ((det .lt. 0.1) .and. (det .gt. -0.1)) then
       !         print '(f12.4,$)', det
       
       do iRot = 1, 3
          do ii = 1, 3
             r1a(ii,iRot) = 0.d0
             r2a(ii,iRot) = 0.d0
             r3a(ii,iRot) = 0.d0
             do jj = 1, 3
                r1a(ii,iRot) = r1a(ii,iRot) + mF(ii,jj,iRot) * r1(jj)
                r2a(ii,iRot) = r2a(ii,iRot) + mF(ii,jj,iRot) * r2(jj)
                r3a(ii,iRot) = r3a(ii,iRot) + mF(ii,jj,iRot) * r3(jj)
             end do
          end do
       end do
       
       detOld = 0.d0
       do iRot = 1, 3
          a = r1a(1,iRot)
          b = r1a(2,iRot)
          c = r1a(3,iRot)
          
          d = r2a(1,iRot)
          e = r2a(2,iRot)
          f = r2a(3,iRot)
          
          g = r3a(1,iRot)
          h = r3a(2,iRot)
          i = r3a(3,iRot)
          
          det = -b*c*d + c*d*e + a*b*f - b*d*f + b*c*g - c*e*g - a*f*h &
               +d*f*h - a*b*i + 2*b*d*i - d*e*i - b*g*i + e*g*i + a*h*i -d*h*i
          
          !$$$         print '(f12.4, $)', det
          if (abs(det) .gt. abs(detOld)) then
             detOld = det
             flag = iRot
          end if
       end do
       
       a = r1a(1,flag)
       b = r1a(2,flag)
       c = r1a(3,flag)
       
       d = r2a(1,flag)
       e = r2a(2,flag)
       f = r2a(3,flag)
       
       g = r3a(1,flag)
       h = r3a(2,flag)
       i = r3a(3,flag)
       det = detOld
       
       !         print '(f12.4, i6, $)', det, flag
    end if
    
    if ((det .lt. 1.e-10) .and. (det .gt. -1.e-10)) then
       print *, 'Error in inv5(). Det = 0.'
       stop
    end if
    
    y(1,1) = -a*f*h + d*f*h + b*d*i - d*e*i - b*g*i + e*g*i + a*h*i - d*h*i
    y(1,2) = a*f*g - d*f*g - a*d*i + d*d*i
    y(1,3) = -b*d*d + a*d*e + b*d*g - a*e*g
    y(1,4) = -a*d + d*d + a*g - d*g
    y(1,5) = -b*d + d*e + b*g - e*g
    y(1,6) = -a*f + d*f + a*i - d*i
    
    y(2,1) = -b*c*h + c*e*h + b*b*i - b*e*i
    y(2,2) = b*c*g - c*e*g - a*b*i + b*d*i - b*g*i + e*g*i + a*h*i -d*h*i
    y(2,3) = -b*b*d + a*b*e + b*d*h - a*e*h
    y(2,4) = -a*b + b*d + a*h - d*h
    y(2,5) = -b*b + b*e + b*h - e*h
    y(2,6) = -b*c + c*e + b*i - e*i
    
    y(3,1) = -c*f*h + b*f*i + c*h*i - b*i*i
    y(3,2) = c*f*g - c*d*i - f*g*i + d*i*i
    y(3,3) = c*d*e - b*d*f - c*e*g + d*f*h + b*d*i - d*e*i + e*g*i -d*h*i
    y(3,4) = -c*d + c*g + d*i - g*i
    y(3,5) = -b*f + f*h + b*i - h*i
    y(3,6) = -c*f + c*i + f*i - i*i
    
    y(4,1) = -b*c*d + c*d*e + a*b*f - b*d*f + b*c*g - c*e*g - a*b*i+ b*d*i
    y(4,2) = -a*f*g + d*f*g + a*d*i - d*d*i
    y(4,3) = b*d*d - a*d*e - b*d*g + a*e*g
    y(4,4) = a*d - d*d - a*g + d*g
    y(4,5) = b*d - d*e - b*g + e*g
    y(4,6) = a*f - d*f - a*i + d*i
    
    y(5,1) = b*c*h - c*e*h - b*b*i + b*e*i
    y(5,2) = -b*c*d + c*d*e + a*b*f - b*d*f - a*f*h + d*f*h + b*d*i- d*e*i
    y(5,3) = b*b*d - a*b*e - b*d*h + a*e*h
    y(5,4) = a*b - b*d - a*h + d*h
    y(5,5) = b*b - b*e - b*h + e*h
    y(5,6) = b*c - c*e - b*i + e*i
    
    y(6,1) = c*f*h - b*f*i - c*h*i + b*i*i
    y(6,2) = -c*f*g + c*d*i + f*g*i - d*i*i
    y(6,3) = -b*c*d + a*b*f + b*c*g - a*f*h - a*b*i + b*d*i - b*g*i+ a*h*i
    y(6,4) = c*d - c*g - d*i + g*i
    y(6,5) = b*f - f*h - b*i + h*i
    y(6,6) = c*f - c*i - f*i + i*i
    
    do jj = 1, 6
       do ii = 1, 6
          y(ii,jj) = y(ii,jj) / det
       end do
    end do
    
    return
    
  end subroutine inv6
  
end module mdutil
