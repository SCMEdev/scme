module calc_higher_order
  
  use data_types
  use max_parameters
  use molecProperties, only: SF
    use printer_mod, only: printer

  implicit none
  
  private
  public calcEhigh
  
contains
  
  subroutine calcEhigh(rCM, opole, hpole, nM, NC, a, a2, uH, eT, dEdr, rMax2, iSlab)
    
    implicit none
    
!    real(dp) rCM(3,maxCoo/3)
!    real(dp) opole(3,3,3,maxCoo/3), hpole(3,3,3,3,maxCoo/3)
!    
!    integer itagl(maxatoms), nM, i, j, k, l, NC, NCz, nx, ny, nz
!    
!    real(dp) eoT(3,maxCoo/3), ehT(3,maxCoo/3), eT(3,maxCoo/3)
!    real(dp) dEdr(3,3,maxCoo/3)
!    !      real(dp) dEo(3,3,maxCoo/3), dEh(3,3,maxCoo/3)
!    real(dp) dr(3), r1, r2, r3, r5, r7, r9, r11, re(3)
!    real(dp) dEdr1(3,3), eO(3), eH(3), a(3), a2(3), u, uH, rMax2,  swFunc, dSdr
!    logical*1 iSlab

!---------------------------------
    ! in/out
    real(dp), intent(in) ::  rCM(:,:)
    real(dp), intent(in) :: opole(:,:,:,:) 
    real(dp), intent(in) :: hpole(:,:,:,:,:) !3,3,3,3,maxCoo/3
    
    integer,   intent(in) :: nM
    integer,   intent(in) :: NC

    logical*1, intent(in) :: iSlab
    real(dp), intent(in) :: a(3), a2(3)
    real(dp), intent(in) :: rMax2

    real(dp), intent(out)   :: uH
    real(dp), intent(out)   :: dEdr(:,:,:) !3,3,maxCoo/3
    real(dp), intent(out)   ::  eT(:,:)

!----------------------------------------
!    real(dp) ::  rCM(:,:)
!    real(dp) :: opole(:,:,:,:) 
!    real(dp) :: hpole(:,:,:,:,:) !3,3,3,3,maxCoo/3
!    
!    integer :: nM, NC
!
!    logical*1 :: iSlab
!    real(dp) :: a(3), a2(3)
!    real(dp) :: rMax2
!
!    real(dp) :: uH
!    real(dp) :: dEdr(:,:,:) !3,3,maxCoo/3
!    real(dp) ::  eT(:,:)


!JÖ internal: 
    integer itagl(maxatoms), i, j, k, l, NCz, nx, ny, nz
!JÖ    real(dp), dimension(size(eT,1),size(eT,2)) :: eoT, ehT !JÖ automatic array
!    real(dp) 
    
    
    real(dp) dEdr1(3,3)
    real(dp) dr(3), re(3) 
    real(dp) eO(3), eH(3) 
    real(dp) r1, r2, r3, r5, r7, r9, r11
    real(dp) u,  swFunc, dSdr
    

    
    NCz = NC
    if (iSlab) NCz = 0
    
    !calcEhigh(rCM, opole, hpole, nM, NC, a, a2, uH, eT, dEdr, rMax2, iSlab)
    
    !!$omp parallel do &
    !!$omp default(private) &
    !!$omp shared( nM, eoT, ehT, eT, dEdr, NC, NCz, rCM, a, a2, rMax2, uH, opole, hpole, iSlab )
    
     
    do i = 1, nM
       
       !do j = 1, 3
       !   eoT(j,i) = 0.d0
       !   ehT(j,i) = 0.d0
       !   eT(j,i) = 0.d0
       !   do k = 1, 3
       !      dEdr(j,k,i) = 0.d0
       !   end do
       !end do
!JÖ       eoT  = 0
!JÖ       ehT  = 0
       eT   = 0
       dEdr = 0
       
       do nx = -NC, NC
          re(1) = a(1) * nx
          do ny = -NC, NC
             re(2) = a(2) * ny
             do nz = -NCz, NCz
                re(3) = a(3) * nz
                
                bicycle:do j = 1, nM !JÖ named loop instead of "goto 11"
                   if ( (j.eq.i) .and. (nx.eq.0) .and. (ny.eq.0) .and. (nz.eq.0)) cycle bicycle !JÖ goto 11
                   do k = 1, 3
                      dr(k) = rCM(k,i) - rCM(k,j)
                      if      (dr(k) .gt. a2(k))  then; dr(k) = dr(k) - a(k)
                      else if (dr(k) .lt. -a2(k)) then; dr(k) = dr(k) + a(k)
                      end if
                      dr(k) = dr(k) + re(k)
                   end do
                   r2 = dr(1)**2 + dr(2)**2 + dr(3)**2
                   
                   if (r2 .gt. rMax2) cycle bicycle !JÖ goto 11
                   
                   r1 = sqrt(r2)
                   call SF(r1, swFunc)
!call printer(swFunc,'swFunc')                 
                   
                   r3 = r1 * r2
                   r5 = r3 * r2
                   r7 = r5 * r2
                   r9 = r7 * r2
                   r11 = r9 * r2
                   
                   call oField(dr, r2, r7, r9, eO, opole(:,:,:,j), u, dEdr1)!JÖ , j) opole -> opole(:,:,:,j)
                   uH = uH + u
                   !do k = 1, 3
                   !   !                        eoT(k,i) = eoT(k,i) + eO(k) * swFunc
                   !   eT(k,i) = eT(k,i) + eO(k) * swFunc
                   !   do l = 1, 3
                   !      dEdr(k,l,i) = dEdr(k,l,i) + dEdr1(k,l) * swFunc
                   !   end do
                   !end do
                   
                   !do k = 1, 3
                   eT(:,i) = eT(:,i) + eO(:) * swFunc
                   dEdr(:,:,i) = dEdr(:,:,i) + dEdr1(:,:) * swFunc
                   
                   call hField(dr, r2, r9, r11, eH, hpole(:,:,:,:,j), u, dEdr1)!JÖ , j) ! hpole -> hpole(:,:,:,:,j)
                   uH = uH + u
                   !do k = 1, 3
                   !   !                        ehT(k,i) = ehT(k,i) + eH(k) * swFunc
                   !   eT(k,i) = eT(k,i) + eH(k) * swFunc
                   !   do l = 1, 3
                   !      dEdr(k,l,i) = dEdr(k,l,i) + dEdr1(k,l) * swFunc
                   !   end do
                   !end do
                   
                   eT(:,i) = eT(:,i) + eH(:) * swFunc !JÖ
                   dEdr(:,:,i) = dEdr(:,:,i) + dEdr1(:,:) * swFunc
                   
                   
                end do bicycle !JÖ 
!11              end do
             end do
          end do
       end do
    end do
    !!$omp end parallel do
    return
    
  end subroutine calcEhigh
  
  !----------------------------------------------------------------------+
  !     Calculate the octopolar field and its derivative                 |
  !----------------------------------------------------------------------+
  pure subroutine oField(dr, r2, r7, r9, eO, opole, u, dEdr)!JÖ, m)
    
    implicit none
    
!JÖ    real(dp) opole(3,3,3,maxCoo/3)
!JÖ    
!JÖ    !c      real(dp) opole(3,3,3)
!JÖ    real(dp) dr(3), r2, r7, r9, eO(3), u, dEdr(3,3)
!JÖ    
!JÖ    integer i, j, k, m
!JÖ    real(dp) rrrO, v(3), g(3,3), d(3,3)
 
!JÖ  in/out: 
    real(dp), intent(in)  :: dr(:)
    real(dp), intent(in)  :: r2, r7, r9
    real(dp), intent(out) :: eO(:) 
    real(dp), intent(in)  :: opole(:,:,:)!JÖ ,:)
    real(dp), intent(out) :: u
    real(dp), intent(out) :: dEdr(:,:)
!JÖ    integer , intent(in)  :: m

!JÖ internal:    
    integer i, j, k
    real(dp) rrrO, v(3), g(3,3), d(3,3)
   

    rrrO = 0 !JÖ 0.d0
    v    = 0 !JÖ v(3)
    g    = 0 !JÖ g(3)


!JÖ    do i = 1, 3
!JÖ       v(i) = 0.d0
!JÖ       do j = 1, 3
!JÖ          g(i,j) = 0.d0
!JÖ       end do
!JÖ    end do
    
    do i = 1, 3
       do j = 1, 3
          do k = 1, 3
!JÖ             g(i,j) = g(i,j) + opole(i,j,k,m) * dr(k)
!JÖ             v(i)   = v(i)   + opole(i,j,k,m) * dr(j) * dr(k)
!JÖ             rrrO   = rrrO   + opole(i,j,k,m) * dr(i) * dr(j) * dr(k)
             g(i,j) = g(i,j) + opole(i,j,k) * dr(k)
             v(i)   = v(i)   + opole(i,j,k) * dr(j) * dr(k)
             rrrO   = rrrO   + opole(i,j,k) * dr(i) * dr(j) * dr(k)
          end do
       end do
    end do
    
    do i = 1, 3
       eO(i) = (7.0_dp * rrrO/r2 * dr(i) - 3.0_dp * v(i)) / r7
       
       do j = i, 3
          dEdr(i,j) = (-6.0_dp * g(i,j) * r2 + 21.0_dp * (v(i) * dr(j) + &
               v(j) * dr(i)) - 63.0_dp * rrrO * dr(i) * dr(j) / r2) / r9
          if (i .eq. j) then
             dEdr(i,j) = dEdr(i,j) + 7.0_dp * rrrO / r9
          end if
       end do
    end do
    
    u = rrrO / (6.0_dp * r7)
    return
    
  end subroutine oField
  
  !----------------------------------------------------------------------+
  !     Calculate the hexadecapolar field and its derivative             |
  !----------------------------------------------------------------------+
  pure subroutine hField(dr, r2, r9, r11, eH, hpole, u, dEdr)!JÖ, m)
    
    implicit none
    
!    integer i, j, k, l, m
!    real(dp) rrrrH, v(3), g(3,3), d(3,3)
!    
!    real(dp) dr(3), r2, r9, r11, eH(3), hpole(3,3,3,3,maxCoo/3), u, dEdr(3,3)

!JÖ in/out    
    real(dp), intent(in)  :: dr(:)
    real(dp), intent(in)  :: r2, r9, r11
    real(dp), intent(out) :: eH(:)
    real(dp), intent(in)  :: hpole(:,:,:,:) !JÖ ,:)
    real(dp), intent(out) :: u
    real(dp), intent(out) :: dEdr(:,:)
!JÖ    integer , intent(in)  :: m
!JÖ internal: 
    real(dp) rrrrH, v(3), g(3,3), d(3,3)
    integer i, j, k, l
   
    rrrrH  = 0 !JÖ 0.d0
    v      = 0 ! array(3)
    g      = 0 ! array(3)
!JÖ    do i = 1, 3
!JÖ       v(i) = 0.d0
!JÖ       do j = 1, 3
!JÖ          g(j,i) = 0.d0
!JÖ       end do
!JÖ    end do
    
    do i = 1, 3
       do j = 1, 3
          do k = 1, 3
             do l = 1, 3
                !                  g(i,j) = g(i,j) + hpole(i,j,k,l) * dr(k) * dr(l)
                !                  v(l) = v(l) + hpole(i,j,k,l) * dr(i) * dr(j) * dr(l)
                
!JÖ                g(l,k) = g(l,k) + hpole(i,j,k,l,m) * dr(i) * dr(j)
!JÖ                v(i) = v(i) + hpole(i,j,k,l,m) * dr(j) * dr(k) * dr(l)
!JÖ                
!JÖ                rrrrH = rrrrH + hpole(i,j,k,l,m) * dr(i) * dr(j) * dr(k) * dr(l)
                g(l,k) = g(l,k) + hpole(i,j,k,l) * dr(i) * dr(j)
                v(i)   = v(i) +   hpole(i,j,k,l) * dr(j) * dr(k) * dr(l)
                rrrrH  = rrrrH +  hpole(i,j,k,l) * dr(i) * dr(j) * dr(k) * dr(l)
             end do
          end do
       end do
    end do
    
    do i = 1, 3
       eH(i) = (9.0_dp * rrrrH / r2 * dr(i) - 4.0_dp * v(i)) / r9
       
       do j = i, 3
          dEdr(i,j) = (-12.0_dp * g(i,j) + 36.0_dp / r2 * (dr(i) &
               * v(j) + dr(j) * v(i)) - 99.0_dp * rrrrH * dr(i) &
               * dr(j) / r2 / r2) / r9
          if (i .eq. j) then
             dEdr(i,j) = dEdr(i,j) + 9.0_dp * rrrrH / r11
          end if
       end do
    end do
    
    u = rrrrH / (24.0_dp * r9)
    return
    
  end subroutine hField
  
end module calc_higher_order
