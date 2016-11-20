module rho
  
  use data_types
  use max_parameters
  use parameters, only: rho_a1, rho_a2, rho_a3, rho_a4, rho_a5, rho_a6, rho_n
  
  implicit none
  
  private
  public calcRho, calcAmp, dDens
  
contains
  
  subroutine calcRho(rho, ra, nM, a, a2, rMax2)
    
    implicit none
    
    real(dp) ra(maxCoo), a(3), a2(3)
    integer nM, i, j, k, iOi, iOj
    real(dp) rho(maxCoo/3), dr(3), dr1, dr2, f, rMax2
    
    do i = 1, nM
       rho(i) = 0.d0
    end do
    
    do i = 1, nM-1
       iOi = 3 * (2*nM + i - 1)
       do j = i+1, nM
          iOj = 3 * (2*nM + j - 1)
          
          do k = 1, 3
             dr(k) = ra(iOj+k) - ra(iOi+k)
             if (dr(k) .gt. a2(k)) then
                dr(k) = dr(k) - a(k)
             else if (dr(k) .lt. -a2(k)) then
                dr(k) = dr(k) + a(k)
             end if
             dr(k) = dr(k)
          end do
          
          dr2 = dr(1)**2 + dr(2)**2 + dr(3)**2
          !            if (dr2 .gt. rMax2) cycle
          dr1 = sqrt(dr2)
          
          call dens(dr1, f)
          
          rho(i) = rho(i) + f
          rho(j) = rho(j) + f
          !            end if
          !            if (i.eq.1) print '(3f20.5)', dr1, 1/dr1**3, rho(i)
       end do
    end do
    
    return
    
  end subroutine calcRho
  
  !-----------------------------------------------------------------------
  subroutine calcAmp(rho, Amp, dAmp, nM)
    
    implicit none
    
    integer nM, i, j, k
    real(dp) rho(maxCoo/3), dr(3), dr1, dr2, f, df, h, dh, rMax2
    real(dp) dAmp(maxCoo/3), Amp(maxCoo/3), rhoAvg
    
    
    do i = 1, nM
       
       !         rhoAvg = rhoAvg + rho(i)
       call hh(rho(i), Amp(i), dAmp(i))
       
       !         print *, rho(i), Amp(i), dAmp(i)
       !         stop
       
    end do
    
    
    !      print *, rhoAvg / nM
    !      stop
    
    return
    
  end subroutine calcAmp
  
  !-----------------------------------------------------------------------
  !     Density
  !
  subroutine dens(r, f)
    
    implicit none
    real(dp) r, f, A, alpha, beta
    parameter (A=2.5d5, alpha=1.5d0, beta=3.d0)
    
    f = A * exp(-r/alpha) / r**beta
    
    !      f = 53000.
    
    return
    
  end subroutine dens
  
  !-----------------------------------------------------------------------
  !     Derivative of the density
  !
  subroutine dDens(r, f)
    
    implicit none
    real(dp) r, f, A, alpha, beta
    parameter (A=2.5d5, alpha=1.5d0, beta=3.d0)
    
    f =- A * exp(-r/alpha) / r**(beta+1) * (beta + r/alpha)
    
    !      f = 0.
    return
    
  end subroutine dDens
  
  !-----------------------------------------------------------------------
  subroutine hh(r, Amp, dAmp)
    
    implicit none
    
    integer n, i
    real(dp) r, rp, Amp, dAmp, a(6), p1, p2
    
    
    a(1) = rho_a1
    a(2) = rho_a2
    a(3) = rho_a3
    a(4)= rho_a4
    a(5)= rho_a5
    a(6)= rho_a6
    n = rho_n
    
    !      if (r .lt. 1465.35283) then
    if (r .lt. 1600.d0) then
       Amp = 0.d0
       dAmp = 0.d0
       !      else if (r .gt. 8275.18577201142) then
    else if (r .gt. 8000.) then
       !     else if (r .gt. 10000.) then
       ! Commented by Fer
       !        Amp = 0.024d0/2d0
       ! New values found for reparametrized potential
       Amp = 0.1750D0/2.0D0
       !        Amp = 0.0610D0/2.0D0
       dAmp = 0.d0
    else
       Amp = a(1)
       dAmp = 0.d0
       rp = 1.d0
       
       do i = 1, n
          dAmp = dAmp + dble(i) * a(i+1) * rp
          rp = rp * r
          Amp = Amp + a(i+1) * rp
       end do
    end if
    
    !      print *, r, Amp, dAmp
    !      stop
    
    return
    
  end subroutine hh
  
end module rho
