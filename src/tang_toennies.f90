module tang_toennies
  
  use data_types
  use parameters, only: tt_b
  
  implicit none
  
  private
  public Tang_ToenniesN, Tang_ToenniesNdF
  
contains
  
  !-----------------------------------------------------------------------
  subroutine Tang_ToenniesN(r, f, n)
    
    implicit none
    real(dp) r, f, b
    integer k, n
    real(dp) x, t
    
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
  subroutine Tang_ToenniesNdF(r, f, df, n)
    
    implicit none
    real(dp) r, f, df, b
    integer k, n
    real(dp) x, t
    
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
