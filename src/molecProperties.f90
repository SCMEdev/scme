!----------------------------------------------------------------------+
!     Routine to fix molecules broken due to the periodic boundary     |
!     conditions                                                       |
!----------------------------------------------------------------------+

!JÖ changed all .d0 and .0d0 to .0_dp
! and other syntax-oriented changes. 

module molecProperties

  use data_types
!  use max_parameters
  !use tang_toennies, only: Tang_ToenniesN, Tang_ToenniesNdF
  use sf_disp_tangtoe, only: SF, SFdsf
implicit none
  private
  public recoverMolecules, add_field_gradients, rotatePolariz, rotate_qoh_poles , recoverMolecules_new
        !, SF, SFdsf, create_rw, calc_cm, addFields, setUnpolPoles, calcCentersOfMass, findPpalAxes,

contains

  subroutine recoverMolecules(raOri, ra, nH, nO, a, a2)
    !-------------------------------------------------------------------
    ! This routine takes the coords(3*n_atoms) and puts out ar(3*natoms)
    ! (both in the HH HH HH HH HH... O O O O... order). Oxygen position
    ! is copied but in case an H is on the wrong side of the box, it is
    ! moved to correct position beside the Oxygen. 
    !    / Description by Jonatan Öström in Dec 2016
    !-------------------------------------------------------------------
    
    !JÖ reformulate to implicit none and state intention and assumed shape
    implicit none
!JÖ    implicit real(dp) (a-h,o-z)
    
    real(dp), 	intent(in)	:: raOri(:), a(:), a2(:)
    integer, 	intent(in) 	:: nH, nO
    real(dp), 	intent(out) :: ra(:)
    !type(h2o),  intent(out) :: rw(:)

	! Internal
    integer  :: i,l,n, index
    real(dp) :: dist

!JÖ    real(dp) raOri(maxCoo), a(3), a2(3), ra(maxCoo), dist
!JÖ    integer nH,nO
    !HH HH HH HH O O O O is the order expected from raOri
	!JÖ change 1 = 0,nO-1 and n = 0,2 to make reindexing more straightforward in loop 
    do i = 0, nO-1 !JÖ oxygen nr i
       do l = 1, 3 !JÖ the x,y,z directions of atom coordinates and box size a(l)
          do n = 0, 1 !JÖ which hydrogen 
             index = l + 3*(n) + 6*(i)
             dist = raOri(index) - raOri(l + 3*(i+nH)) ! H-position - O-positon
             if (dist .gt. a2(l)) then                 ! a2 = a/2 = boxdimension/2
                ra(index) = raOri(index) - a(l)        
             elseif(dist .lt. -a2(l)) then
                ra(index) = raOri(index) + a(l)
             else
                ra(index) = raOri(index)
             end if
          end do
          ra(l + 3*(i+nH)) = raOri(l + 3*(i+nH)) !oxygen is copied over regardless
       end do
    end do
    !return

  end subroutine recoverMolecules
  
  subroutine recoverMolecules_new(coord, ra, nM, a, a2)
    implicit none
    real(dp), 	intent(in)	:: coord(:,:) !j hho
    real(dp), 	intent(out) :: ra(:,:)    !j hho
    integer, 	intent(in) 	:: nM
    real(dp), 	intent(in)	:: a(:), a2(:)

	! Internal
    integer  :: m,q,h, ih, io
    real(dp) :: dist

    do m = 1, nM !j molecule
       do q = 1, 3 !j x,y,z 
          do h = 1, 2 !j hydrogen 
             ih = (m-1)*3+h
             io = (m-1)*3+3
             dist = coord(q,ih) - coord(q,io)          ! H-position - O-positon
             if (dist .gt. a2(q)) then                 ! a2 = a/2 = boxdimension/2
                ra(q,ih) = coord(q,ih) - a(q)
             elseif(dist .lt. -a2(q)) then
                ra(q,ih) = coord(q,ih) + a(q)
             else
                ra(q,ih) = coord(q,ih)
             end if
          end do
          ra(q,io) = coord(q,io) !oxygen is copied regardless
       end do
    end do
    !return

  end subroutine !





  subroutine rotate_qoh_poles(q0, o0, h0, qpole, opole, hpole, nM, x) !JÖ d0,dpole
    implicit none
    !jö nodip! real(dp), intent(in) :: d0(3) 
    real(dp), intent(in) :: q0(3,3), o0(3,3,3), h0(3,3,3,3), x(3,3,3)
    integer, intent(in) :: nM

    !jö nodip! real(dp), intent(out) :: dpole(:,:)
    real(dp), intent(out) :: qpole(:,:,:), opole(:,:,:,:), hpole(:,:,:,:,:)

! Internal:
    integer i, j, k, l, ii, jj, kk, ll, m 

    ! Zeroise.
    hpole = 0
    opole = 0
    qpole = 0
!    dpole = 0
    
    do m = 1, nM
       do l = 1, 3
          do ll = 1, 3
             do k = 1, 3
                do kk = 1, 3
                   do j = 1, 3
                      do jj = 1, 3
                         do i = 1, 3
                            do ii = 1, 3
                               hpole(i,j,k,l,m) = hpole(i,j,k,l,m) + &
                                    x(i,ii,m) * x(j,jj,m) * x(k,kk,m) &
                                    * x(l,ll,m) * h0(ii,jj,kk,ll)
                            end do
                         end do
                         opole(j,k,l,m) = opole(j,k,l,m) + &
                              x(j,jj,m)* x(k,kk,m)* x(l,ll,m) * o0(jj,kk,ll)
                      end do
                   end do
                   qpole(k,l,m) = qpole(k,l,m) + x(k,kk,m) * x(l,ll,m)* q0(kk,ll)
                end do
             end do
             !jö no dip !dpole(l,m) = dpole(l,m) + x(l,ll,m) * d0(ll)
          end do
       end do
    end do

    return

  end subroutine rotate_qoh_poles


  !----------------------------------------------------------------------+
  !> Routine to rotate ??? to the orientation of the molecule.
  !! @param[in] dd0 : ?
  !! @param[in] dq0 : ?
  !! @param[in] qq0 : ?
  !! @param[in] hp0 : ?
  !! @param[out] dd : The array is assumed to be zeroized.
  !! @param[out] dq : ?
  !! @param[out] qq : ?
  !! @param[out] hp : ?
  !! @param[out] nM : The number of molecules.
  !! @param[out] x  : ?
  subroutine rotatePolariz(dd0, dq0, qq0, hp0, dd, dq, qq, hp, nM, x)

    implicit none
!    real(dp), intent(in) :: hp0(3,3,3)
!    real(dp), intent(in) :: dq0(3,3,3)
!    real(dp), intent(in) :: dd0(3,3)
!    real(dp), intent(in) :: qq0(3,3,3,3)
!    real(dp), intent(out) :: hp(3,3,3,nM)
!    real(dp), intent(out) :: dq(3,3,3,nM)
!    real(dp), intent(out) :: dd(3,3,nM)
!    real(dp), intent(out) :: qq(3,3,3,3,nM)
!    integer, intent(in) :: nM
!    real(dp), intent(in) :: x(3,3,nM)
    ! ----------------------------------------

    real(dp), intent(in) :: hp0(:,:,:)
    real(dp), intent(in) :: dq0(:,:,:)
    real(dp), intent(in) :: dd0(:,:)
    real(dp), intent(in) :: qq0(:,:,:,:)
    real(dp), intent(out) :: hp(:,:,:,:)
    real(dp), intent(out) :: dq(:,:,:,:)
    real(dp), intent(out) :: dd(:,:,:)
    real(dp), intent(out) :: qq(:,:,:,:,:)
    integer , intent(in) :: nM
    real(dp), intent(in) :: x(:,:,:)
    ! ----------------------------------------




    ! Local variables.
    integer i, j, k, l, ii, jj, kk, ll, m

    ! Zeroise.
!JÖ    dd(1:3, 1:3, 1:nM) = 0.0_dp
!JÖ    dq(1:3, 1:3, 1:3, 1:nM) = 0.0_dp
!JÖ    qq(1:3, 1:3, 1:3, 1:3, 1:nM) = 0.0_dp
!JÖ    hp(1:3, 1:3, 1:3, 1:nM) = 0.0_dp

    dd = 0 !JÖ
    dq = 0 !JÖ
    qq = 0 !JÖ
    hp = 0 !JÖ


    ! Do the calculation.
    do m = 1, nM
       do k = 1, 3
          do kk = 1, 3
             do j = 1, 3
                do jj = 1, 3
                   do i = 1, 3
                      do ii = 1, 3
                         do l = 1, 3
                            do ll = 1, 3
                               qq(l,i,j,k,m) = qq(l,i,j,k,m) + &
                                    x(l,ll,m) *  x(i,ii,m)*x(j,jj,m)&
                                    * x(k,kk,m) * qq0(ll,ii,jj,kk)
                            end do
                         end do
                         dq(i,j,k,m) = dq(i,j,k,m) + x(i,ii,m) * x(j,jj,m) &
                              * x(k,kk,m) * dq0(ii,jj,kk)
                         hp(i,j,k,m) = hp(i,j,k,m) + x(i,ii,m) * x(j,jj,m) &
                              * x(k,kk,m) * hp0(ii,jj,kk)
                      end do
                   end do
                   dd(j,k,m) = dd(j,k,m) + x(j,jj,m) * x(k,kk,m)* dd0(jj,kk)
                end do
             end do
          end do
       end do
    end do

    return

  end subroutine rotatePolariz


  !----------------------------------------------------------------------+
  !     Add the derivative of all the fields                             |
  !----------------------------------------------------------------------+
  subroutine add_field_gradients(dEhdr, dEddr, dEtdr, nM)

    implicit none
    integer, intent(in) :: nM !JÖ
!    real(dp) dEhdr(3,3,maxCoo/3)
!    real(dp) dEtdr(3,3,maxCoo/3), dEddr(3,3,maxCoo/3)
    real(dp), intent(out) :: dEtdr(3,3,nM)
    real(dp), intent(in)  :: dEhdr(3,3,nM)
    real(dp), intent(in)  :: dEddr(3,3,nM)

    integer i, j, k !JÖ, nM
    
    !JÖ:
    dEtdr = dEhdr + dEddr !half filled, above and on diagonal
    do j = 2, 3
       do k = 1, j-1
          dEtdr(j,k,:) = dEtdr(k,j,:) !fill rest => symmetric
       end do
    end do

  end subroutine 


end module 

!from addDfields
!    !do i = 1, nM
!    !end do
!    
!    
!    return

!JÖ    do i = 1, nM
!JÖ       do j = 1, 3
!JÖ          do k = 1, 3
!JÖ             dEtdr(k,j,i) = dEhdr(k,j,i) + dEddr(k,j,i)
!JÖ             !               dEtdr(j,k,i) = dEhdr(k,j,i)
!JÖ          end do
!JÖ       end do
!JÖ       do j = 2, 3
!JÖ          do k = 1, j-1
!JÖ             dEtdr(j,k,i) = dEtdr(k,j,i)
!JÖ          end do
!JÖ       end do
!JÖ    end do


