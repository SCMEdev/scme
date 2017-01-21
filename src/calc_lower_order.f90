module calc_lower_order
  
  use data_types
!  use max_parameters
  !use molecProperties, only: SF
  use sf_disp_tangtoe, only: SF
  
  implicit none
  
  private
  public dip_quadField
  
contains
  
  subroutine dip_quadField(rCM, dpole, qpole, nM, NC, a, a2, uD, uQ, eT, dEdr, rMax2, iSlab)
    
    implicit none
    
!JÖ    real(dp) rCM(3,maxCoo/3), dpole(3,maxCoo/3), qpole(3,3,maxCoo/3), a(3), a2(3)
!JÖ    real(dp) dEdr(3,3,maxCoo/3), uD, uQ, eT(3,maxCoo/3)

!JÖ here we should change those goto statements to a break/exit statement of named loop

    real(dp), intent(out) :: eT   (:,:) !JÖ 3,maxCoo/3
    real(dp), intent(in)  :: rCM  (:,:)
    real(dp), intent(in)  :: dpole(:,:)
    real(dp), intent(in)  :: qpole(:,:,:) !JÖ 3,3,maxCoo/3
    real(dp), intent(out) :: dEdr (:,:,:)
    real(dp), intent(in)  :: a(:), a2(:) !JÖ 3
    real(dp), intent(out) :: uD, uQ
    real(dp), intent(in)  :: rMax2
    integer,  intent(in)   :: nM, NC 
    logical*1,intent(in) :: iSlab
    
    real(dp), save :: re(3), eD(3), eq(3), dr(3) !JÖ 3
!JÖ    real(dp), save :: uDv(nM),uQv(nM)
    real(dp), save :: dEdr1(3,3) !JÖ 3,3
    real(dp), save :: r1, r2, r3, r5, r7, u, swFunc, dSdr 
    integer , save :: jj, kk, i, j, k, l, NCz, nx, ny, nz
    
!, save
!, save
!, save
!, save
    
    NCz = NC
    if (iSlab) NCz = 0
    
    uD = 0.0_dp
    uQ = 0.0_dp
!    !$omp parallel do shared(rCM, a, a2, NC, uD, uQ, NCz,rMax2, dEdr, nM) &
!    !$omp private(eT, jj, kk, nx, i, j, k, re, ny, nz, dr, r2, r1, swFunc, r3, r5, r7, eD, dpole, u, dEdr1, eq, qpole) &
!    !$omp default(none)
    do i = 1, nM
       do jj = 1, 3
          eT(jj,i) = 0.0_dp
          do kk = 1, 3
             dEdr(kk,jj,i) = 0.0_dp
          end do
       end do
       do nx = -NC, NC
          re(1) = a(1) * nx
          do ny = -NC, NC
             re(2) = a(2) * ny
             do nz = -NCz, NCz
                re(3) = a(3) * nz
                
                bike:do j = 1, nM
                   if ( (j.eq.i) .and. (nx.eq.0) .and. (ny.eq.0) .and. (nz.eq.0)) cycle bike !goto 11
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
                   
                   if (r2 .gt. rMax2) cycle bike !goto 11
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
                   
                !   !for i = 1:3
                !     eT(:,i)     = eT(:,i)     + eD(:)*swFunc
                !       !for j=1:3
                !       dEdr(:,:,i) = dEdr(:,:,i) + dEdr1(:,:)*swFunc
                !   !endfor**2
                   
                   
                   !     Quadrupole Field
                   call qField(dr, r2, r5, r7, eq, qpole, u, dEdr1, j)
                   uQ = uQ + u
                   !do k = 1, 3
                   !   eT(k,i) = eT(k,i) + eq(k) * swFunc
                   !   do l = 1, 3
                   !      dEdr(k,l,i) = dEdr(k,l,i) + dEdr1(k,l) * swFunc
                   !   end do
                   !end do
                      eT(:,i) = eT(:,i) + eq(:) * swFunc
                         dEdr(:,:,i) = dEdr(:,:,i) + dEdr1(:,:) * swFunc
                   
                   
!11              end do
                end do bike
             end do
          end do
       end do
    end do
!    !$omp end parallel do 
    return
    
  end subroutine dip_quadField
  
  !----------------------------------------------------------------------+
  !     Calculate the dipolar field and its derivative                   |
  !----------------------------------------------------------------------+
  pure subroutine dField(dr, r2, r3, r5, eD, dpole, u, dEdr, m)
    
    implicit none
    
!    real(dp) dpole(3,maxCoo/3), eD(3), dr(3), r2, r3, r5, dEdr(3,3), u

    real(dp), intent(in)  :: dpole(:,:) !JÖ 3,maxCoo/3
    real(dp), intent(out) :: dEdr(:,:) !3,3
    real(dp), intent(out) :: eD(:)!3
    real(dp), intent(in)  :: dr(:) !3
    real(dp), intent(in)  :: r2, r3, r5
    real(dp), intent(inout)  :: u
    integer , intent(in)  :: m
    
    !internal: 
    real(dp) :: mDr, temp 
    integer  :: i, j, k
!JÖ    real(dp) mDr
!, save  
!, save      
    mDr = dpole(1,m)*dr(1) + dpole(2,m)*dr(2) + dpole(3,m)*dr(3)
    
    do i = 1, 3
       eD(i) = (3.0_dp * mDr * dr(i) / r2 - dpole(i,m)) / r3
       do j = i, 3
          dEdr(i,j) = (dpole(i,m) * dr(j) + dpole(j,m) * dr(i) - 5.0_dp &
               * mDr * dr(i) * dr(j) / r2) * 3.0_dp / r5
          if (i.eq.j) then
             dEdr(i,j) = dEdr(i,j) + mDr * 3.0_dp / r5
          end if
       end do
    end do
    
    !JÖ>>
    !eD(:) = (3.0_dp * mDr * dr(:) / r2 - dpole(:,m)) / r3
    !do i = 1, 3
    !  do j = i, 3
    !    dEdr(i,j) = (dpole(i,m) * dr(j) + dpole(j,m) * dr(i) - 5.0_dp * mDr * dr(i) * dr(j) / r2) * 3.0_dp / r5
    !  enddo
    !enddo
    
    !temp = mDr * 3.0_dp / r5
    !dEdr(1,1) = dEdr(1,1) + temp ! + mDr * 3.0_dp / r5
    !dEdr(2,2) = dEdr(2,2) + temp ! + mDr * 3.0_dp / r5
    !dEdr(3,3) = dEdr(3,3) + temp ! + mDr * 3.0_dp / r5
    !enddo
    !<<JÖ
    
    !      u = mDr / r3
    return
    
  end subroutine dField
  
  !----------------------------------------------------------------------+
  !     Calculate the octopolar field and its derivative                 |
  !----------------------------------------------------------------------+
  pure subroutine qField(dr, r2, r5, r7, eq, qpole, u, dEdr, m)
    
    implicit none
!JÖ    real(dp) qpole(3,3,maxCoo/3)
    real(dp), intent(in)  :: dr(:)
    real(dp), intent(in)  :: r2, r5, r7
    real(dp), intent(out) :: eq(:) !3
    real(dp), intent(in)  :: qpole(:,:,:), u !JÖ 3,3,maxCoo/3
    real(dp), intent(out) :: dEdr(:,:) !3,3
    integer , intent(in)  :: m

!JÖ internal: 
    integer  :: i, j
    real(dp) :: v(3), rQr, temp 
!, save
!, save
    
!JÖ    integer i, j, m
!JÖ    real(dp) dr(3), r2, r5, r7, u, dEdr(3,3), eq(3)
!JÖ    real(dp) v(3), rQr
    
!JÖ    do i = 1, 3
!JÖ       v(i) = 0.0_dp
!JÖ    end do
    v=0 !JÖ
    
    rQr = 0.0_dp
    do j = 1, 3
       do i = 1, 3
          v(i) = v(i) + qpole(i,j,m) * dr(j)
          rQr = rQr + dr(i) * qpole(i,j,m) * dr(j)
       end do
    end do
    
!oJÖ    eq(:) = 2.0_dp * (2.5_dp * rQr / r2 * dr(:) - v(:)) / r5 !JÖ

    do i = 1, 3
       eq(i) = 2.0_dp * (2.5_dp * rQr / r2 * dr(i) - v(i)) / r5
       do j = i, 3
          dEdr(i,j) = (-2.0_dp * qpole(i,j,m) * r2 + 10.0_dp * (v(j) &
               * dr(i) + v(i) * dr(j)) - 35.0_dp * rQr * dr(i) * dr(j) / r2) / r7
          if (i.eq.j) then
             dEdr(i,j) = dEdr(i,j) + 5.0_dp * rQr / r7
          end if
       end do
    end do

!oJÖ    do i = 1, 3
!oJÖ       do j = i, 3
!oJÖ          dEdr(i,j) = (-2.0_dp * qpole(i,j,m) * r2 + 10.0_dp * (v(j) &
!oJÖ               * dr(i) + v(i) * dr(j)) - 35.0_dp * rQr * dr(i) * dr(j) / r2) / r7
!oJÖ       end do
!oJÖ    end do

    
!    temp = 5.0_dp * rQr / r7
!    dEdr(1,1) = dEdr(1,1) + temp
!    dEdr(2,2) = dEdr(2,2) + temp
!    dEdr(3,3) = dEdr(3,3) + temp
    
    
    !      u = rQr / r5
    return
    
  end subroutine qField
  
end module calc_lower_order
