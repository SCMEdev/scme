module coreInt_mod
  
  use data_types
  use max_parameters
  use parameters, only: coreInt_c1, coreInt_c2, coreInt_c3, coreInt_c4, coreInt_c5_r
  
  use rho, only: calcRho, calcAmp, dDens
  
  implicit none
  
  private
  public coreInt
  
contains
  
  
  subroutine coreInt(ra, fa, uCore, nM, a, a2)
    
    implicit none
    real(dp), intent(out) :: ra(:), fa(:)
    real(dp), intent(out) :: uCore
    integer,  intent(in) :: nM
    real(dp), intent(in) :: a(3), a2(3)
    
!JÖ internal:
    integer ai(10), kk, jj
    real(dp) t1, t2, df, r, avg
    real(dp) c1, c2, c3, c4, c5,t11, t12, t13
    real(dp),dimension(size(ra)/3) :: Amp, dAmp, rho
    real(dp) f, rMax2
    real(dp) Ri1(3), Ri(3), Ri2(3), Rj1(3), Rj(3), Rj2(3)
    real(dp) vRi1j1(3), vRi2j2(3), vRij(3)
    real(dp) Rij
    real(dp) S
    real(dp) uR_An
    real(dp) an_1, an_2, an_3, an_4
    real(dp) df1, df2
    real(dp) c5_r
    integer*4 Pt_i, Pt_i1, Pt_i2
    integer*4 Pt_j, Pt_j1, Pt_j2
    
    integer i, j, k, n, m, iOn, iOm, iOi, iOj, p
    real(dp) dr(3)

!--------------
!    integer nM, ai(10), kk, jj
!    real(dp) t1, t2, df, uCore, r, a(3), a2(3), avg
!    real(dp) ra(maxCoo), fa(maxCoo), c1, c2, c3, c4, c5,t11, t12, t13
!    real(dp) Amp(maxCoo/3), dAmp(maxCoo/3), rho(maxCoo/3), f, rMax2
!    real(dp) Ri1(3), Ri(3), Ri2(3), Rj1(3), Rj(3), Rj2(3)
!    real(dp) vRi1j1(3), vRi2j2(3), vRij(3)
!    real(dp) Rij
!    real(dp) S
!    real(dp) uR_An
!    real(dp) an_1, an_2, an_3, an_4
!    real(dp) df1, df2
!    real(dp) c5_r
!    integer*4 Pt_i, Pt_i1, Pt_i2
!    integer*4 Pt_j, Pt_j1, Pt_j2
!    
!    integer i, j, k, n, m, iOn, iOm, iOi, iOj, p
!    real(dp) dr(3)
!--------------
    
    c1 = coreInt_c1
    c2 = coreInt_c2
    c3 = coreInt_c3
    c4 = coreInt_c4
    c5_r = coreInt_c5_r
    
    
    call calcRho(rho, ra, nM, a, a2, rMax2)
    call calcAmp(rho, Amp, dAmp, nM)
    
    ! Debug
    !     write(10,FMT='(A,I3,A,F16.10)')
    !    &      'Rho(',int(nM/2),') = ', rho(int(nM/2))
    !     do i=1,nM
    !       write(10,FMT='(A,I3,A,F16.10)')
    !    &        'Rho(',i,') = ', rho(i)
    !     end do
    !     write(10,FMT='(A,I3,A,F16.10)')
    !    &      'Amp(',int(nM/2),') = ', Amp(int(nM/2))
    !     if ( nM .le. 10 ) then
    !       do i=1,nM
    !         write(10,FMT='(I3,3X,2F16.10)')  i, rho(i), Amp(i)
    !       end do
    !     else
    !       write(10,FMT='(I3,3X,2F16.10)') int(nM/2), rho(int(nM/2)),
    !    &                                          Amp(int(nM/2))
    !     end if
    
    uCore = 0.0_dp
    
    
    do n = 1, nM-1
       
       ! Get the index of the first O atom
       iOn = 3 * (2*nM + n - 1)
       
       do m = n+1, nM
          
          ! Get the index of the second O atom
          iOm = 3 * (2*nM + m - 1)
          
          ! Adjust the O-O distance for the PBC's
          do i = 1, 3
             dr(i) = ra(iOm+i) - ra(iOn+i)
             if     (dr(i) .gt. a2(i)) then
                dr(i) = dr(i) - a(i)
             elseif (dr(i) .lt. -a2(i)) then
                dr(i) = dr(i) + a(i)
             end if
          end do
          
          r = dsqrt(dr(1)**2 + dr(2)**2 + dr(3)**2)
          
          t11 = dexp(c3*r)
          t12 = dexp(c3*r/c4)
          t13 = r**c1
          
          ! Debug
          if ( c5_r .ge. 0.0_dp ) then
             c5 = c5_r
          else
             c5 = (Amp(n)+Amp(m))
          end if
          
          t1 = c2*t13*(t11 + c5*t12)
          !         t1 = c2*t13*t11
          
          uCore = uCore + t1
          
          
          t2 = c2*t13*t12
          
          df = (c2*(c1/r+c3)*t13*t11 + c2*c5*(c1/r+c3/c4)*t13*t12) / r
          !         df = (c2*(c1/r+c3)*t13*t11) / r
          
          do i = 1, 3
             fa(iOn+i) = fa(iOn+i) + df * dr(i)
             fa(iOm+i) = fa(iOm+i) - df * dr(i)
          end do
          
          !     Derivative of the embedding part (Amplitude)
          !
          ! Debug
          if (c5_r .lt. 0.0_dp ) then
             do j = 1, nM
                if (j .eq. n) goto 11
                iOj = 3 * (2*nM + j - 1)
                
                do k = 1, 3
                   dr(k) = ra(iOj+k) - ra(iOn+k)
                   if (dr(k) .gt. a2(k)) then
                      dr(k) = dr(k) - a(k)
                   else if (dr(k) .lt. -a2(k)) then
                      dr(k) = dr(k) + a(k)
                   end if
                end do
                
                r = sqrt(dr(1)**2 + dr(2)**2 + dr(3)**2)
                
                call dDens(r, f)
                
                df = dAmp(n) * f * t2 / r
                do k = 1, 3
                   fa(iOn+k) = fa(iOn+k) + df * dr(k)
                   fa(iOj+k) = fa(iOj+k) - df * dr(k)
                end do
11           end do
             
             do j = 1, nM
                if (j .eq. m) goto 12
                iOj = 3 * (2*nM + j - 1)
                
                do k = 1, 3
                   dr(k) = ra(iOj+k) - ra(iOm+k)
                   if (dr(k) .gt. a2(k)) then
                      dr(k) = dr(k) - a(k)
                   else if (dr(k) .lt. -a2(k)) then
                      dr(k) = dr(k) + a(k)
                   end if
                end do
                
                r = sqrt(dr(1)**2 + dr(2)**2 + dr(3)**2)
                call dDens(r, f)
                
                df = dAmp(m) * f * t2 / r
                do k = 1, 3
                   fa(iOm+k) = fa(iOm+k) + df * dr(k)
                   fa(iOj+k) = fa(iOj+k) - df * dr(k)
                end do
12           end do
          end if
          
       end do
    end do
    
  end subroutine coreInt
  
end module coreInt_mod
