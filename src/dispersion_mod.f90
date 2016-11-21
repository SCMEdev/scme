module dispersion_mod
  
  use data_types
  use max_parameters
  
  implicit none
  
  private
  public dispersion
  
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
          call Tang_Toennies(r, f6, df6, f8, df8, f10, df10)
          
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
  subroutine Tang_Toennies(r, f6, df6, f8, df8, f10, df10)
    
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
    
  end subroutine Tang_Toennies
  
end module dispersion_mod
