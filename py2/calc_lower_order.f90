module calc_lower_order
  
  use data_types
  use max_parameters
  use molecProperties, only: SF
  
  implicit none
  
  private
  public calcEdip_quad
  
contains
  
  subroutine calcEdip_quad(rCM, dpole, qpole, nM, NC, a, a2, uD, uQ, eT, dEdr, rMax2, iSlab)
    
    implicit none
    
    real(dp) rCM(3,maxCoo/3), dpole(3,maxCoo/3), qpole(3,3,maxCoo/3), a(3), a2(3)
    real(dp) dEdr(3,3,maxCoo/3), uD, uQ, eT(3,maxCoo/3)
    
    real(dp) re(3), r1, r2, r3, r5, r7, eD(3), eq(3), dr(3), dEdr1(3,3), u, rMax2, swFunc, dSdr
    integer i, j, k, l, nM, NC, NCz, nx, ny, nz
    logical*1 iSlab
    
    NCz = NC
    if (iSlab) NCz = 0
    
    uD = 0.d0
    uQ = 0.d0
    do i = 1, nM
       do j = 1, 3
          eT(j,i) = 0.d0
          do k = 1, 3
             dEdr(k,j,i) = 0.d0
          end do
       end do
       do nx = -NC, NC
          re(1) = a(1) * nx
          do ny = -NC, NC
             re(2) = a(2) * ny
             do nz = -NCz, NCz
                re(3) = a(3) * nz
                
                do j = 1, nM
                   if ( (j.eq.i) .and. (nx.eq.0) .and. (ny.eq.0) .and. (nz.eq.0)) goto 11
                   do k = 1, 3
                      dr(k) = rCM(k,i) - rCM(k,j)
                      if (dr(k) .gt. a2(k)) then
                         dr(k) = dr(k) - a(k)
                      else if (dr(k) .lt. -a2(k)) then
                         dr(k) = dr(k) + a(k)
                      end if
                      dr(k) = dr(k) + re(k)
                   end do
                   r2 = dr(1)**2 + dr(2)**2 + dr(3)**2
                   
                   if (r2 .gt. rMax2) goto 11
                   r1 = sqrt(r2)
                   call SF(r1, swFunc)
                   
                   r3 = r1 * r2
                   r5 = r3 * r2
                   r7 = r5 * r2
                   
                   !     Dipole Field
                   call dField(dr, r2, r3, r5, eD, dpole, u, dEdr1, j)
                   uD = uD + u
                   do k = 1, 3
                      eT(k,i) = eT(k,i) + eD(k) * swFunc
                      do l = 1, 3
                         dEdr(k,l,i) = dEdr(k,l,i) + dEdr1(k,l) * swFunc
                      end do
                   end do
                   
                   !     Quadrupole Field
                   call qField(dr, r2, r5, r7, eq, qpole, u, dEdr1, j)
                   uQ = uQ + u
                   do k = 1, 3
                      eT(k,i) = eT(k,i) + eq(k) * swFunc
                      do l = 1, 3
                         dEdr(k,l,i) = dEdr(k,l,i) + dEdr1(k,l) * swFunc
                      end do
                   end do
                   
11              end do
             end do
          end do
       end do
    end do
    
    return
    
  end subroutine calcEdip_quad
  
  !----------------------------------------------------------------------+
  !     Calculate the dipolar field and its derivative                   |
  !----------------------------------------------------------------------+
  subroutine dField(dr, r2, r3, r5, eD, dpole, u, dEdr, m)
    
    implicit none
    
    real(dp) dpole(3,maxCoo/3), eD(3), dr(3), r2, r3, r5, dEdr(3,3), u
    
    integer i, j, k, m
    real(dp) mDr
    
    mDr = dpole(1,m)*dr(1) + dpole(2,m)*dr(2) + dpole(3,m)*dr(3)
    do i = 1, 3
       eD(i) = (3.d0 * mDr * dr(i) / r2 - dpole(i,m)) / r3
       do j = i, 3
          dEdr(i,j) = (dpole(i,m) * dr(j) + dpole(j,m) * dr(i) - 5.d0 &
               * mDr * dr(i) * dr(j) / r2) * 3.d0 / r5
          if (i.eq.j) then
             dEdr(i,j) = dEdr(i,j) + mDr * 3.d0 / r5
          end if
       end do
    end do
    
    !      u = mDr / r3
    return
    
  end subroutine dField
  
  !----------------------------------------------------------------------+
  !     Calculate the octopolar field and its derivative                 |
  !----------------------------------------------------------------------+
  subroutine qField(dr, r2, r5, r7, eq, qpole, u, dEdr, m)
    
    implicit none
    real(dp) qpole(3,3,maxCoo/3)
    
    integer i, j, m
    real(dp) dr(3), r2, r5, r7, u, dEdr(3,3), eq(3)
    real(dp) v(3), rQr
    
    do i = 1, 3
       v(i) = 0.d0
    end do
    rQr = 0.d0
    do j = 1, 3
       do i = 1, 3
          v(i) = v(i) + qpole(i,j,m) * dr(j)
          rQr = rQr + dr(i) * qpole(i,j,m) * dr(j)
       end do
    end do
    
    do i = 1, 3
       eq(i) = 2.d0 * (2.5d0 * rQr / r2 * dr(i) - v(i)) / r5
       do j = i, 3
          dEdr(i,j) = (-2.d0 * qpole(i,j,m) * r2 + 10.d0 * (v(j) &
               * dr(i) + v(i) * dr(j)) - 35.d0 * rQr * dr(i) * dr(j) / r2) / r7
          if (i.eq.j) then
             dEdr(i,j) = dEdr(i,j) + 5.d0 * rQr / r7
          end if
       end do
    end do
    
    !      u = rQr / r5
    return
    
  end subroutine qField
  
end module calc_lower_order
