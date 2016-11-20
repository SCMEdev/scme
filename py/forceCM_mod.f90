module forceCM_mod
  
  use data_types
  use max_parameters
  
  implicit none
  
  private
  public forceCM
  
contains
  
  subroutine forceCM(dpole, qpole, opole, hpole, d2v, d3v, d4v, d5v,nM, fsf, fCM)
    
    implicit none
    
    !     Work multipoles. They start unpolarized and with the induction
    !     loop we induce dipoles and quadrupoles.
    real(dp) dpole(3,maxCoo/3), qpole(3,3,maxCoo/3)
    real(dp) opole(3,3,3,maxCoo/3), hpole(3,3,3,3,maxCoo/3)
    
    !     High order derivatives of the potential
    real(dp) d2v(3,3,maxCoo/3), d3v(3,3,3,maxCoo/3)
    real(dp) d4v(3,3,3,3,maxCoo/3), d5v(3,3,3,3,3,maxCoo/3)
    
    integer nM
    real(dp) fCM(3,maxCoo/3), fsf(3,maxCoo/3)
    
    integer n, i, j, k, l, s
    real(dp) f
    
    
    do n = 1, nM
       
       do i = 1, 3
          f = 0.d0
          do j = 1, 3
             !     Force on the dipole
             f = f - d2v(j,i,n)*dpole(j,n)
             do k = 1, 3
                !     Force on the quadrupole
                f = f - d3v(k,j,i,n)*qpole(k,j,n) / 3.d0
                do l = 1, 3
                   !     Force on the octopole
                   f = f - d4v(l,k,j,i,n)*opole(l,k,j,n) / 15.d0
                   do s = 1, 3
                      !     Force on the hexadecapole
                      f = f - d5v(s,l,k,j,i,n)*hpole(s,l,k,j,n) / 105.d0
                   end do
                end do
             end do
          end do
          fCM(i,n) = f + fsf(i,n)
       end do
    end do
    
    return
    
  end subroutine forceCM
  
end module forceCM_mod
