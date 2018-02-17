module force_torqueCM
  
  use data_types
  !use max_parameters
  
  implicit none
  
  private
  public forceCM, torqueCM
  
contains

  subroutine forceCM(dpole, qpole, opole, hpole, d2v, d3v, d4v, d5v,nM, fsf, fCM)
    
    implicit none
    
    !     Work multipoles. They start unpolarized and with the induction
    !     loop we induce dipoles and quadrupoles.
    real(dp), intent(in) :: dpole(:,:), qpole(:,:,:)
    real(dp), intent(in) :: opole(:,:,:,:), hpole(:,:,:,:,:)
    
    !     High order derivatives of the potential
    real(dp), intent(in) :: d2v(:,:,:), d3v(:,:,:,:)
    real(dp), intent(in) :: d4v(:,:,:,:,:), d5v(:,:,:,:,:,:)
    
    integer,  intent(in) ::  nM
    real(dp), intent(out) ::  fsf(:,:), fCM(:,:)

!JÖ internal:
    integer n, i, j, k, l, s
    real(dp) f

!--------------------------------------
!    real(dp) dpole(3,maxCoo/3), qpole(3,3,maxCoo/3)
!    real(dp) opole(3,3,3,maxCoo/3), hpole(3,3,3,3,maxCoo/3)
!    
!    !     High order derivatives of the potential
!    real(dp) d2v(3,3,maxCoo/3), d3v(3,3,3,maxCoo/3)
!    real(dp) d4v(3,3,3,3,maxCoo/3), d5v(3,3,3,3,3,maxCoo/3)
!    
!    integer nM
!    real(dp) fCM(3,maxCoo/3), fsf(3,maxCoo/3)
!    
!    integer n, i, j, k, l, s
!    real(dp) f
    
    
    do n = 1, nM
       
       do i = 1, 3
          f = 0.0_dp
          do j = 1, 3
             !     Force on the dipole
             f = f - d2v(j,i,n)*dpole(j,n)
             do k = 1, 3
                !     Force on the quadrupole
                f = f - d3v(k,j,i,n)*qpole(k,j,n) / 3.0_dp
                do l = 1, 3
                   !     Force on the octopole
                   f = f - d4v(l,k,j,i,n)*opole(l,k,j,n) / 15.0_dp
                   do s = 1, 3
                      !     Force on the hexadecapole
                      f = f - d5v(s,l,k,j,i,n)*hpole(s,l,k,j,n) / 105.0_dp
                   end do
                end do
             end do
          end do
          fCM(i,n) = f + fsf(i,n)
       end do
    end do
    
    return
    
  end subroutine forceCM

  
  subroutine torqueCM(dpole, qpole, opole, hpole, d1v, d2v, d3v, d4v, nM, tau)
    
    implicit none
    
    !     Work multipoles. They start unpolarized and with the induction
    !     loop we induce dipoles and quadrupoles.
    real(dp), intent(in) :: dpole(:,:), qpole(:,:,:)
    real(dp), intent(in) :: opole(:,:,:,:), hpole(:,:,:,:,:)
    
    !     High order derivatives of the potential
    real(dp), intent(in) :: d1v(:,:), d2v(:,:,:)
    real(dp), intent(in) :: d3v(:,:,:,:), d4v(:,:,:,:,:)
    
    integer, intent(in)  :: nM
    real(dp), intent(out)  :: tau(:,:) 

!JÖ internal
    integer n, i, j, k!, l, s
    real(dp) t(3)!, asy
!-------------------------------
!    real(dp) dpole(3,maxCoo/3), qpole(3,3,maxCoo/3)
!    real(dp) opole(3,3,3,maxCoo/3), hpole(3,3,3,3,maxCoo/3)
!    
!    !     High order derivatives of the potential
!    real(dp) d1v(3,maxCoo/3), d2v(3,3,maxCoo/3), d3v(3,3,3,maxCoo/3)
!    real(dp) d4v(3,3,3,3,maxCoo/3)
!    
!    integer nM
!    real(dp) tau(3,maxCoo/3), asy
!    
!    integer n, i, j, k, l, s
!    real(dp) t(3)


    
    do n = 1, nM
       
       !     Torque on the dipole PxE
       t(1) = -(dpole(2,n)*d1v(3,n) - dpole(3,n)*d1v(2,n))
       t(2) = -(dpole(3,n)*d1v(1,n) - dpole(1,n)*d1v(3,n))
       t(3) = -(dpole(1,n)*d1v(2,n) - dpole(2,n)*d1v(1,n))
       
       !     Torque on the quadrupoles
       do i = 1, 3
          t(1) = t(1) - 2.0_dp / 3.0_dp * &
               (qpole(2,i,n)*d2v(3,i,n) - qpole(3,i,n)*d2v(2,i,n))
          t(2) = t(2) - 2.0_dp / 3.0_dp * &
               (qpole(3,i,n)*d2v(1,i,n) - qpole(1,i,n)*d2v(3,i,n))
          t(3) = t(3) - 2.0_dp / 3.0_dp * &
               (qpole(1,i,n)*d2v(2,i,n) - qpole(2,i,n)*d2v(1,i,n))
       end do
       
       !     Torque on the octopoles
       do j = 1, 3
          do i = 1, 3
             t(1) = t(1) - 3.0_dp / 15.0_dp * &
                  (opole(2,i,j,n)*d3v(3,i,j,n) - opole(3,i,j,n)*d3v(2,i,j,n))
             t(2) = t(2) - 3.0_dp / 15.0_dp * &
                  (opole(3,i,j,n)*d3v(1,i,j,n) - opole(1,i,j,n)*d3v(3,i,j,n))
             t(3) = t(3) - 3.0_dp / 15.0_dp * &
                  (opole(1,i,j,n)*d3v(2,i,j,n) - opole(2,i,j,n)*d3v(1,i,j,n))
          end do
       end do
       
       !     Torque on the hexadecapoles
       do k = 1, 3
          do j = 1, 3
             do i = 1, 3
                t(1) = t(1) - 4.0_dp / 105.0_dp * &
                     (hpole(2,i,j,k,n)*d4v(3,i,j,k,n) -hpole(3,i,j,k,n)*d4v(2,i,j,k,n))
                t(2) = t(2) - 4.0_dp / 105.0_dp * &
                     (hpole(3,i,j,k,n)*d4v(1,i,j,k,n) -hpole(1,i,j,k,n)*d4v(3,i,j,k,n))
                t(3) = t(3) - 4.0_dp / 105.0_dp * &
                     (hpole(1,i,j,k,n)*d4v(2,i,j,k,n) -hpole(2,i,j,k,n)*d4v(1,i,j,k,n))
             end do
          end do
       end do
       
       do i = 1, 3
          tau(i,n) = t(i)
       end do
    end do
    
    return
    
  end subroutine torqueCM
  
end module 
