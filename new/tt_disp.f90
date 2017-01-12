subroutine Tang_Toennies(r,fsave)! f6, df6, f8, df8, f10, df10)
    
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    real(dp) b, r,  x, t 
    integer k, fact
    !      parameter (b = 2.48d0)
    parameter (b = 4.4_dp)
    !
    ! NEW shit:
    integer krit(3), ik
    real(dp) fsave(2,3), f,df
    parameter ( krit = [6,8,10] )
    
    x = b * r
    
    f = 1.0_dp
    df = 0.0_dp
    
    t = exp(-x)
    
    ik = 1
    
    do k = 0, 10
       f = f - t
       df = df - t * b * (-1.0_dp + k / x)
       t = t * x / (k+1)
       if (k==krit(ik)) fsave(:,ik)=[f,df]; ik=ik+1
    end do
    
  end subroutine Tang_Toennies
  
!  subroutine Tang_Toennies(r,fsave)! f6, df6, f8, df8, f10, df10)
!    
!    implicit none
!    integer, parameter :: dp = kind(1.0d0)
!    real(dp) b, r !, f6, df6, f8, df8, f10, df10
!    !real(dp) ff6, ff8, ff10, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11
!    !real(dp) x12, x13, x14
!    
!
!    
!    integer k, fact
!    real(dp) x, t
!    !      parameter (b = 2.48d0)
!    parameter (b = 4.4_dp)
!    !
!    integer krit(3), ik , from, till!, start(3)
!    real(dp) fsave(2,3), f,df
!    parameter ( krit = [6,8,10] )!, start = [0,7,9] )
!    
!    x = b * r
!    
!    f = 1.0_dp
!    df = 0.0_dp
!    
!    t = dexp(-x)
!    
!    from=0
!    do ik=1,3
!       till=krit(ik)
!       do k = from, till
!          f = f - t
!          df = df - t * b * (-1.0_dp + k / x)
!          t = t * x / (k+1)
!       end do
!       from=till+1
!       fsave(:,ik)=[f,df]
!    enddo
!
!  end subroutine Tang_Toennies
  
!
Hej Houri och Grace jag jobbar gärna imorgon. 
Idag fick jag veta att man jobbar själv mellan 19 och 21, 
men detta känner jag mig inte bekväm med. 
Det har dragit ut på tiden att hitta lägenhet men jag tror jag har något bra på gång nu. 
Om jag slipper jobba själv 19-21 så är jag väldigt positiv :)
