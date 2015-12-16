module torqueCM_mod
  
  use data_types
  use max_parameters
  
  implicit none
  
  private
  public torqueCM
  
contains
  
  subroutine torqueCM(dpole, qpole, opole, hpole, d1v, d2v, d3v, d4v, nM, tau)
    
    implicit none
    
    !     Work multipoles. They start unpolarized and with the induction
    !     loop we induce dipoles and quadrupoles.
    real(dp) dpole(3,maxCoo/3), qpole(3,3,maxCoo/3)
    real(dp) opole(3,3,3,maxCoo/3), hpole(3,3,3,3,maxCoo/3)
    
    !     High order derivatives of the potential
    real(dp) d1v(3,maxCoo/3), d2v(3,3,maxCoo/3), d3v(3,3,3,maxCoo/3)
    real(dp) d4v(3,3,3,3,maxCoo/3)
    
    integer nM
    real(dp) tau(3,maxCoo/3), asy
    
    integer n, i, j, k, l, s
    real(dp) t(3)
    
    do n = 1, nM
       
       !     Torque on the dipole PxE
       t(1) = -(dpole(2,n)*d1v(3,n) - dpole(3,n)*d1v(2,n))
       t(2) = -(dpole(3,n)*d1v(1,n) - dpole(1,n)*d1v(3,n))
       t(3) = -(dpole(1,n)*d1v(2,n) - dpole(2,n)*d1v(1,n))
       
       !     Torque on the quadrupoles
       do i = 1, 3
          t(1) = t(1) - 2.d0 / 3.d0 * &
               (qpole(2,i,n)*d2v(3,i,n) - qpole(3,i,n)*d2v(2,i,n))
          t(2) = t(2) - 2.d0 / 3.d0 * &
               (qpole(3,i,n)*d2v(1,i,n) - qpole(1,i,n)*d2v(3,i,n))
          t(3) = t(3) - 2.d0 / 3.d0 * &
               (qpole(1,i,n)*d2v(2,i,n) - qpole(2,i,n)*d2v(1,i,n))
       end do
       
       !     Torque on the octopoles
       do j = 1, 3
          do i = 1, 3
             t(1) = t(1) - 3.d0 / 15.d0 * &
                  (opole(2,i,j,n)*d3v(3,i,j,n) - opole(3,i,j,n)*d3v(2,i,j,n))
             t(2) = t(2) - 3.d0 / 15.d0 * &
                  (opole(3,i,j,n)*d3v(1,i,j,n) - opole(1,i,j,n)*d3v(3,i,j,n))
             t(3) = t(3) - 3.d0 / 15.d0 * &
                  (opole(1,i,j,n)*d3v(2,i,j,n) - opole(2,i,j,n)*d3v(1,i,j,n))
          end do
       end do
       
       !     Torque on the hexadecapoles
       do k = 1, 3
          do j = 1, 3
             do i = 1, 3
                t(1) = t(1) - 4.d0 / 105.d0 * &
                     (hpole(2,i,j,k,n)*d4v(3,i,j,k,n) -hpole(3,i,j,k,n)*d4v(2,i,j,k,n))
                t(2) = t(2) - 4.d0 / 105.d0 * &
                     (hpole(3,i,j,k,n)*d4v(1,i,j,k,n) -hpole(1,i,j,k,n)*d4v(3,i,j,k,n))
                t(3) = t(3) - 4.d0 / 105.d0 * &
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
  
end module torqueCM_mod
