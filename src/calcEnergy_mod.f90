module calcEnergy_mod
  
  use data_types
!  use max_parameters
  
  implicit none
  
  private
  public calcEnergy
  
contains
  
  subroutine calcEnergy(dpole, qpole, opole, hpole, d1v, d2v, d3v,d4v, nM, uTot)
    
    implicit none
    
    !     Work multipoles. They start unpolarized and with the induction
    !     loop we induce dipoles and quadrupoles.
    real(dp), intent(in) :: dpole(:,:), qpole(:,:,:)
    real(dp), intent(in) :: opole(:,:,:,:), hpole(:,:,:,:,:)
    
    !     High order derivatives of the potential
    real(dp), intent(in) :: d1v(:,:), d2v(:,:,:), d3v(:,:,:,:)
    real(dp), intent(in) :: d4v(:,:,:,:,:)
    integer,  intent(in) :: nM
    real(dp), intent(out) :: uTot

    real(dp) ud, uq, uo, uh, du
    
    
    integer n, i, j, k, l, s
    real(dp) u
    
    !PRB -- ud, uq, uo, uh in eV and the conversion factor -------------
    
    real(dp) udev, uqev, uoev, uhev!JÖ, convFactor

!    real(dp) dpole(3,maxCoo/3), qpole(3,3,maxCoo/3)
!    real(dp) opole(3,3,3,maxCoo/3), hpole(3,3,3,3,maxCoo/3)
!    
!    !     High order derivatives of the potential
!    real(dp) d1v(3,maxCoo/3), d2v(3,3,maxCoo/3), d3v(3,3,3,maxCoo/3)
!    real(dp) d4v(3,3,3,3,maxCoo/3), uTot, ud, uq, uo, uh, du
!    
!    integer nM
!    
!    integer n, i, j, k, l, s
!    real(dp) u
!    
!    !PRB -- ud, uq, uo, uh in eV and the conversion factor -------------
!    
!    real(dp) udev, uqev, uoev, uhev, convFactor


    !JÖconvFactor = 14.39975841d0 / 4.803206799d0**2
    !PRB ---------------------------------------------------------------
    
    uTot = 0
    ud   = 0
    uq   = 0
    uo   = 0
    uh   = 0
    
    do n = 1, nM
       
       do i = 1, 3
          !     Energy of the dipole
          du = d1v(i,n) * dpole(i,n)
          uTot = uTot + du
          ud = ud + du
          do j = 1, 3
             !     Energy of the quadrupole
             du = d2v(j,i,n)*qpole(j,i,n) / 3.0_dp
             uTot = uTot + du
             uq = uq + du
             do k = 1, 3
                !     Energy of the octopole
                du = d3v(k,j,i,n)*opole(k,j,i,n) / 15.0_dp
                uTot = uTot + du
                uo = uo + du
                
                do l = 1, 3
                   !     Energy of the hexadecapole
                   du = d4v(l,k,j,i,n)*hpole(l,k,j,i,n) / 105.0_dp
                   uTot = uTot + du
                   uh = uh + du
                end do
             end do
          end do
          
       end do
    end do
    uTot = uTot / 2.0_dp

!JÖ seems unnessecary:    
!JÖ    ud = ud / 2.0_dp
!JÖ    uq = uq / 2.0_dp
!JÖ    uo = uo / 2.0_dp
!JÖ    uh = uh / 2.0_dp
!JÖ    
!JÖ    udev = convFactor * ud
!JÖ    uqev = convFactor * uq
!JÖ    uoev = convFactor * uo
!JÖ    uhev = convFactor * uh
    
!    open(999, file="ESenergies.dat", ACCESS='APPEND')
!    write(999, '(4f16.10)') udev, uqev, uoev, uhev
    
    return
    
  end subroutine calcEnergy
  
  !-----------------------------------------------------------------------
  subroutine calcEnergyI(dpole, dpole0, qpole, qpole0, opole, hpole,d1v, d2v, d3v, d4v, nM, uTot)
    
    implicit none
    
    !     Work multipoles. They start unpolarized and with the induction
    !     loop we induce dipoles and quadrupoles.
    real(dp), intent(in) :: dpole(:,:), qpole(:,:,:)
    real(dp), intent(in) :: dpole0(:,:), qpole0(:,:,:)
    real(dp), intent(in) :: opole(:,:,:,:), hpole(:,:,:,:,:)
    
    !     High order derivatives of the potential
    real(dp), intent(in) :: d1v(:,:), d2v(:,:,:), d3v(:,:,:,:)
    real(dp), intent(in) :: d4v(:,:,:,:,:)
    real(dp), intent(out) :: uTot
    
    integer,  intent(in) :: nM

!JÖ internal    
    integer n, i, j, k, l, s
    real(dp) u
!------------------------------
!    real(dp) dpole(3,maxCoo/3), qpole(3,3,maxCoo/3)
!    real(dp) dpole0(3,maxCoo/3), qpole0(3,3,maxCoo/3)
!    real(dp) opole(3,3,3,maxCoo/3), hpole(3,3,3,3,maxCoo/3)
!    
!    !     High order derivatives of the potential
!    real(dp) d1v(3,maxCoo/3), d2v(3,3,maxCoo/3), d3v(3,3,3,maxCoo/3)
!    real(dp) d4v(3,3,3,3,maxCoo/3), uTot
!    
!    integer nM
!    
!    integer n, i, j, k, l, s
!    real(dp) u


    
    uTot = 0.0_dp
    do n = 1, nM
       
       do i = 1, 3
          !     Energy of the dipole
          uTot = uTot + d1v(i,n) * (dpole(i,n) - dpole0(i,n))
          do j = 1, 3
             !     Energy of the quadrupole
             uTot = uTot + d2v(j,i,n)*(qpole(j,i,n)-qpole0(j,i,n)) / 3.0_dp
          end do
          
       end do
    end do
    uTot = uTot / 2.0_dp
    
    return
    
  end subroutine calcEnergyI
  
  !-----------------------------------------------------------------------
  subroutine polarizationEnergy1(hp, dpole, dpole0, qpole, qpole0,d1v, d2v, nM, uPol)
    
    implicit none
    
    !     Work multipoles. They start unpolarized and with the induction
    !     loop we induce dipoles and quadrupoles.
    real(dp), intent(in) :: dpole(:,:), qpole(:,:,:)
    real(dp), intent(in) :: dpole0(:,:), qpole0(:,:,:)
    
    !     High order derivatives of the potential
    real(dp), intent(in)  :: d1v(:,:), d2v(:,:,:)
    real(dp), intent(in)  :: hp(:,:,:,:)
    real(dp), intent(out) :: uPol
    integer,  intent(in)  ::  nM

!JÖ internal:
    integer i, j, k, n
    real(dp) u
    
    uPol = 0.0_dp
    do n = 1, nM
       !     Dipole polarization
       u = -0.5d0 * ( (dpole(1,n)-dpole0(1,n)) * d1v(1,n) +  &
            (dpole(2,n)-dpole0(2,n)) * d1v(2,n) + &
            (dpole(3,n)-dpole0(3,n)) * d1v(3,n) )
       uPol = uPol + u
       
    end do
    
    return
    
  end subroutine polarizationEnergy1
  
  !-----------------------------------------------------------------------
  subroutine polarizationEnergy(dd, dq, qq, hp, d1v, d2v, nM, uPol)
    
    implicit none
    
    !     Work multipoles. They start unpolarized and with the induction
    !     loop we induce dipoles and quadrupoles.
    
    !     High order derivatives of the potential
    
    real(dp), intent(in)  :: dd(:,:,:), dq(:,:,:,:)
    real(dp), intent(in)  :: hp(:,:,:,:), qq(:,:,:,:,:)
    real(dp), intent(in)  :: d1v(:,:), d2v(:,:,:)
    integer,  intent(in)  :: nM
    real(dp), intent(out) ::  uPol
    
    integer i, j, k, l, n
    real(dp) u
!---------------    
    
!    real(dp) d1v(3,maxCoo/3), d2v(3,3,maxCoo/3), uPol
!    real(dp) dd(3,3,maxCoo/3), dq(3,3,3,maxCoo/3), hp(3,3,3,maxCoo/3)
!    real(dp) qq(3,3,3,3,maxCoo/3)
!    
!    integer i, j, k, l, nM
!    
!    integer n
!    real(dp) u
    
    uPol = 0.0_dp
    do n = 1, nM
       
       !     Dipole-dipole polarization
       do i = 1, 3
          do j = 1, 3
             uPol = uPol + 0.5d0 * dd(i,j,n) *d1v(i,n)*d1v(j,n)
             
             !     Dipole-quadrupole polarization
             do k = 1, 3
                uPol = uPol + dq(i,j,k,n) *d1v(i,n)*d2v(j,k,n) / 3.0_dp
                
                !     Quadrupole-quadrupole polarization
                do l = 1, 3
                   uPol = uPol + qq(i,j,k,l,n) *d2v(i,j,n)*d2v(k,l,n)/ 6.0_dp
                end do
                
                !     first hyperpolarization
                !                  uPol = uPol - hp(i,j,k,n) * d1v(i,n)*d1v(j,n)*d1v(k,n)
                !     $                 / 6.0_dp
             end do
          end do
       end do
    end do
    
    return
    
  end subroutine polarizationEnergy
  
end module calcEnergy_mod
