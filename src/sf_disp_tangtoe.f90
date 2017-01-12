module sf_disp_tangtoe
! Module w/ swiching func. (SF,SFdsf,swFunc), dipsersion and 3 versions of tang toennies. 
  
  use data_types
  use parameters, only: tt_b !remember to theck if it is correct that only dispersion has a different b!?
  
  implicit none
  
  private
  public SF, SFdsf, dispersion !Tang_ToenniesN, Tang_ToenniesNdF
  
contains

  subroutine dispersion(ra, fa, uDisp, nM, a, a2)
    
    implicit none
    
    real(dp), intent(in) :: ra(:)
    real(dp), intent(out) :: fa(:)
    real(dp), intent(out) :: uDisp
    integer, intent(in) :: nM
    real(dp), intent(in) :: a(3), a2(3)

!JÖ internal    
    real(dp) t1, t2, df, r
    real(dp) r2, r6, r7, r8, r9, r10, r11
    real(dp) f6, df6, f8, df8, f10, df10
    
    integer n, m, iOn, iOm, i
    real(dp) dr(3), sc
    real(dp), parameter :: C6 = 46.4430d0 * 0.597527378d0
    real(dp), parameter :: C8  = 1141.7000d0 * 0.167324732d0
    real(dp), parameter :: C10 = 33441.0000d0 * 0.046855703d0

!------------------
!    integer nM
!    real(dp) t1, t2, df, r, a(3), a2(3), uDisp
!    real(dp) r2, r6, r7, r8, r9, r10, r11
!    real(dp) f6, df6, f8, df8, f10, df10
!    real(dp) ra(maxCoo), fa(maxCoo), C6, C8, C10
!    
!    integer n, m, iOn, iOm, i
!    real(dp) dr(3), sc
!    parameter (C6 = 46.4430d0 * 0.597527378d0, C8  = 1141.7000d0 * 0.167324732d0, &
!         C10 = 33441.0000d0 * 0.046855703d0)
!------------------
    
    !     Dispersion coefficients from Watts-Coker
    !      Set I
    !      parameter (C6 = 37.2484d0*1.15, C8  = 224.48d0*1.15,
    !     $     C10 = 1560.92d0*1.5 * 1.15)
    
    !      parameter (C6 = 37.2484d0, C8  = 224.48d0,
    !     $     C10 = 1560.92d0*1.5)
    
    uDisp = 0.0_dp
    do n = 1, nM-1
       !      do n = 1, 4
       iOn = 3 * (2*nM + n - 1)
       do m = n+1, nM
          !         do m = 5, 8
          iOm = 3 * (2*nM + m - 1)
          
          do i = 1, 3
             dr(i) = ra(iOm+i) - ra(iOn+i)
             if (dr(i) .gt. a2(i)) then
                dr(i) = dr(i) - a(i)
             else if (dr(i) .lt. -a2(i)) then
                dr(i) = dr(i) + a(i)
             end if
          end do
          
          r2 = dr(1)**2 + dr(2)**2 + dr(3)**2
          r = sqrt(r2)
          call tang_toennies_disp(r, f6, df6, f8, df8, f10, df10)
          
          r6 = r**6
          r7 = r6 * r
          r8 = r7 * r
          r9 = r8 * r
          r10 = r9 * r
          r11 = r10 * r
          
          uDisp = uDisp - C6 / r6 * f6 - C8 / r8 * f8 - C10 / r10* f10
          
          do i = 1, 3
             df = -C6 * (6.0_dp * f6 / r7 - df6 / r6)
             df = df - C8 * (8.0_dp * f8 / r9 - df8 / r8)
             df = df - C10 * (10.0_dp * f10 / r11 - df10 / r10)
             df = df * dr(i) / r
             fa(iOn+i) = fa(iOn+i) - df
             fa(iOm+i) = fa(iOm+i) + df
          end do
       end do
    end do
    
    return
    
  end subroutine dispersion
  
  !-----------------------------------------------------------------------
  subroutine tang_toennies_disp(r, f6, df6, f8, df8, f10, df10)
    
    implicit none
    real(dp) b, r, f6, df6, f8, df8, f10, df10
    real(dp) ff6, ff8, ff10, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11
    real(dp) x12, x13, x14
    
    integer k, fact
    real(dp) x, t
    !      parameter (b = 2.48d0)
    parameter (b = 4.4_dp)
    
    x = b * r
    
    f6 = 1.0_dp
    df6 = 0.0_dp
    
    t = exp(-x)
    do k = 0, 6
       f6 = f6 - t
       df6 = df6 - t * b * (-1.0_dp + k / x)
       t = t * x / (k+1)
    end do
    
    f8 = f6
    df8 = df6
    do k = 7, 8
       f8 = f8 - t
       df8 = df8 - t * b * (-1.0_dp + k / x)
       t = t * x / (k+1)
    end do
    
    f10 = f8
    df10 = df8
    do k = 9, 10
       f10 = f10 - t
       df10 = df10 - t * b * (-1.0_dp + k / x)
       t = t * x / (k+1)
    end do
    
    return
    
  end subroutine tang_toennies_disp


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

  
  !-----------------------------------------------------------------------
  pure subroutine Tang_ToenniesN(r, f, n)
    
!JÖ    implicit none
!JÖ    real(dp) r, f, b
!JÖ    integer k, n
!JÖ    real(dp) x, t

    implicit none
    real(dp), intent(in)  :: r
    real(dp), intent(out) :: f
    integer , intent(in)  :: n
    
    integer  :: k
    real(dp) ::  x, t, b
    !JÖ , save
    !JÖ , save
    b = tt_b
    
    x = b * r
    
    f = 1.d0
    
    t = exp(-x)
    do k = 0, n
       f = f - t
       t = t * x / (k+1)
    end do
    
    f = sqrt(f)
    
    !      f = 1.d0
    !      df = 0.d0
    return
    
  end subroutine Tang_ToenniesN
  
  !-----------------------------------------------------------------------
  pure subroutine Tang_ToenniesNdF(r, f, df, n)
    
!JÖ    implicit none
!JÖ    real(dp) r, f, df, b
!JÖ    integer k, n
!JÖ    real(dp) x, t

    implicit none
    real(dp), intent(in)  :: r
    real(dp), intent(out) :: f
    real(dp), intent(out) :: df
    integer , intent(in)  :: n

!internal:    
    integer  :: k
    real(dp) :: x, t, b
    !JÖ , save
    !JÖ , save
    b = tt_b
    
    x = b * r
    
    f = 1.d0
    df = 0.d0
    
    t = exp(-x)
    do k = 0, n
       f = f - t
       df = df - t * b * (-1.d0 + k / x)
       t = t * x / (k+1)
    end do
    
    f = sqrt(f)
    df = df / (2.d0 * f)
    
    !      f = 1.d0
    !      df = 0.d0
    return
    
  end subroutine Tang_ToenniesNdF
  
end module 
