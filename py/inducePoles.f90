module inducePoles
  
  use data_types
  use max_parameters
  
  implicit none
  
  private
  public induceDipole, induceQpole
  
contains
  
  !----------------------------------------------------------------------+
  !     Induce dipole moment                                             |
  !----------------------------------------------------------------------+
  subroutine induceDipole(dpole, dpole0, eT, dEtdr, dd, dq, hp, nM, converged)
    
    implicit none
    integer nM, i, j, k, m
    real(dp) dpole(3,maxCoo/3), dpole0(3,maxCoo/3)
    real(dp) dd(3,3,maxCoo/3), dq(3,3,3,maxCoo/3), hp(3,3,3,maxCoo/3)
    real(dp) eT(3,maxCoo/3), dEtdr(3,3,maxCoo/3), daux
    logical*1 converged
    
    do m = 1, nM
       do i = 1, 3
          daux = dpole(i,m)
          dpole(i,m) = dpole0(i,m)
          do j = 1, 3
             dpole(i,m) = dpole(i,m) + dd(i,j,m) * eT(j,m)
          end do
          do j = 1, 3
             do k = 1, 3
                dpole(i,m) = dpole(i,m) + dq(i,j,k,m) * dEtdr(j,k,m) / 3.d0 + &
                     hp(i,j,k,m) * eT(j,m) * eT(k,m) / 2.d0
             end do
          end do

          if (abs(daux - dpole(i,m)) .gt. 1.e-7) converged = .false.
       end do
    end do
    
    return
    
  end subroutine induceDipole
  
  !----------------------------------------------------------------------+
  !     Induce quadrupole moment                                         |
  !----------------------------------------------------------------------+
  subroutine induceQpole(qpole, qpole0, eT, dEtdr, dq, qq, nM, converged)
    
    implicit none
    
    integer nM, i, j, k, l, m
    real(dp) qpole(3,3,maxCoo/3), dq(3,3,3,maxCoo/3)
    real(dp) qq(3,3,3,3,maxCoo/3), qaux
    real(dp) qpole0(3,3,maxCoo/3), eT(3,maxCoo/3), dEtdr(3,3,maxCoo/3)
    logical*1 converged
    
    do m = 1, nM
       
       do j = 1, 3
          do i = 1, 3
             qaux = qpole(i,j,m)
             qpole(i,j,m) = qpole0(i,j,m)
             
             do k = 1, 3
                qPole(i,j,m) = qPole(i,j,m) + dq(k,i,j,m) * eT(k,m)
             end do
             do k = 1, 3
                do l = 1, 3
                   qPole(i,j,m) = qPole(i,j,m) + qq(i,j,k,l,m) * dEtdr(k,l,m)
                end do
             end do
             if (abs(qaux - qPole(i,j,m)) .gt. 1.e-7) converged = .false.
          end do
       end do
       
    end do
    
    return
    
  end subroutine induceQpole
  
end module inducePoles
