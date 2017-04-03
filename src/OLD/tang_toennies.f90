module tang_toennies
  
  use data_types
  use parameters, only: tt_b
  
  implicit none
  
  private
  public Tang_ToenniesN, Tang_ToenniesNdF
  
contains
  
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
  
end module tang_toennies
