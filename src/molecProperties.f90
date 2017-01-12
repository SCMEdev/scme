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

  private
  public recoverMolecules, calcCentersOfMass, findPpalAxes, rotatePoles,&
       rotatePolariz, setUnpolPoles, addFields, addDfields !, SF, SFdsf

contains

  subroutine recoverMolecules(raOri, ra, nH, nO, a, a2)
    !-------------------------------------------------------------------
    ! Thos routine takes the coords(3*n_atoms) and puts out ar(3*natoms)
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

  !----------------------------------------------------------------------+
  !     Routine to calculate the center of mass of each molecule         |
  !----------------------------------------------------------------------+
  subroutine calcCentersOfMass(ra, nM, rCM)

!JÖ    implicit real(dp) (a-h,o-z)
    implicit none !JÖ
    !JÖ state intents and assumed shape:
    real(dp), 	intent(in) 	:: ra(:)
    integer, 	intent(in) 	:: nM
    real(dp), 	intent(out) :: rCM(:,:)
    
    integer iH1, iH2, iO, i, j !JÖ nM, 
!JÖ    real(dp) rCM(3,maxCoo/3), ra(maxCoo)

    do i = 1, nM
       iH2 = 2 * i
       iH1 = iH2 - 1
       iO  = 2*nM + i !JÖ

       iH1 = 3 * (iH1 - 1)
       iH2 = 3 * (iH2 - 1)
       iO  = 3 * (iO - 1)
       do j = 1, 3
          rCM(j,i) = (ra(iH1+j) + ra(iH2+j) + 16.0_dp * ra(iO+j)) / 18.0_dp
       end do
    end do

    return

  end subroutine calcCentersOfMass

  !----------------------------------------------------------------------+
  !     Routine to calculate the principal axes of each molecule         |
  !----------------------------------------------------------------------+
  subroutine findPpalAxesOLD(ra, nM, x)

!JÖ    implicit real(dp) (a-h,o-z)
    implicit none
    
    real(dp), 	intent(in) 	:: ra(:)
    integer, 	intent(in) 	:: nM
    real(dp), 	intent(out) :: x(:,:,:) 

    integer 	:: iH1, iH2, iO, i, j !JÖ nM, 
    real(dp) 	:: r11, r21
!jö    real(dp) x(3,3,maxCoo/3), ra(maxCoo), r11, r21


    do i = 1, nM
       iH2 = 2 * i
       iH1 = iH2 - 1
       iO  = 2 * nM + i

       iH1 = 3 * (iH1 - 1)
       iH2 = 3 * (iH2 - 1)
       iO  = 3 * (iO - 1)
       do j = 1, 3 ! j = x,y,z coord of H1,H2,O
          x(j,1,i) = -(ra(iH1+j) + ra(iH2+j) - 2.0_dp * ra(iO+j)) !x(:,1,i) = -{ H1(:) + H2(:) - 2*O(:) } 
          x(j,2,i) = ra(iH2+j) - ra(iH1+j)                        !x(:,2,i) = H2(:) - H1(:) 
       end do
       r11 = sqrt(x(1,1,i)*x(1,1,i) + x(2,1,i)*x(2,1,i) + x(3,1,i) * x(3,1,i)) !norm(x(:,1,i)),  ":" is x,y,z
       r21 = sqrt(x(1,2,i)*x(1,2,i) + x(2,2,i)*x(2,2,i) + x(3,2,i) * x(3,2,i)) !norm(x(:,2,i))
       do j = 1, 3
          x(j,1,i) = x(j,1,i) / r11
          x(j,2,i) = x(j,2,i) / r21
       end do
       x(1,3,i) = x(2,1,i) * x(3,2,i) - x(3,1,i) * x(2,2,i)
       x(2,3,i) = x(3,1,i) * x(1,2,i) - x(1,1,i) * x(3,2,i)
       x(3,3,i) = x(1,1,i) * x(2,2,i) - x(2,1,i) * x(1,2,i)
    end do
    return

  end subroutine findPpalAxesOLD

  !----------------------------------------------------------------------+
  !     Routine to calculate the principal axes of each molecule         |
  !----------------------------------------------------------------------+
  subroutine findPpalAxes(ra, nM, x)
  !---------------------------------------------------------------------
  ! This routine finds the principal axis of the molecule but does it 
  ! apply to a water molecule with non-fixed geometry?
  !---------------------------------------------------------------------

!JÖ    implicit real(dp) (a-h,o-z)
!JÖ
!JÖ    integer nM, iH1, iH2, iO, i, j
!JÖ    real(dp) x(3,3,maxCoo/3), ra(maxCoo), r11, r21

!JÖ use instead:
    implicit none
    
    real(dp), 	intent(in) 	:: ra(:)
    integer, 	intent(in) 	:: nM
    real(dp), 	intent(out) :: x(:,:,:) 

    integer 	:: iH1, iH2, iO, i, j !JÖ nM, 
    real(dp) 	:: r11, r21


    ! Debug
    integer*4 p

    do i = 1, nM
       iH2 = 2 * i
       iH1 = iH2 - 1
       iO  = 2 * nM + i

       iH1 = 3 * (iH1 - 1)
       iH2 = 3 * (iH2 - 1)
       iO  = 3 * (iO - 1)
       do j = 1, 3 ! j = x,y,z coord of H1,H2,O
          x(j,3,i) = -(ra(iH1+j) + ra(iH2+j) - 2.0_dp * ra(iO+j)) !x(:,1,i) = -{ H1(:) + H2(:) - 2*O(:) } |
          x(j,1,i) = ra(iH2+j) - ra(iH1+j)                        !x(:,2,i) = H2(:) - H1(:)               | These find two vectors in the plane of the molecule
       end do
       r11 = sqrt(x(1,3,i)*x(1,3,i) + x(2,3,i)*x(2,3,i) + x(3,3,i) * x(3,3,i)) !norm(x(:,1,i)),  ":" is x,y,z
       r21 = sqrt(x(1,1,i)*x(1,1,i) + x(2,1,i)*x(2,1,i) + x(3,1,i) * x(3,1,i)) !norm(x(:,2,i))

       do j = 1, 3
          x(j,3,i) = x(j,3,i) / r11 !normalize?
          x(j,1,i) = x(j,1,i) / r21
       end do
       x(1,2,i) = x(2,3,i) * x(3,1,i) - x(3,3,i) * x(2,1,i) !cross product ? find middle vector orthogonal to plane
       x(2,2,i) = x(3,3,i) * x(1,1,i) - x(1,3,i) * x(3,1,i) !do we need all three of them? Yes!
       x(3,2,i) = x(1,3,i) * x(2,1,i) - x(2,3,i) * x(1,1,i)
    end do

    ! Debug
    ! Print the principal axes matrix for each molecule
    !     do i=1,nM
    !       write(10,FMT='(A,I4)') ' Axes for Molecule: ', i
    !       write(10,FMT='(3F10.5)') (x(1,p,i),p=1,3)
    !       write(10,FMT='(3F10.5)') (x(2,p,i),p=1,3)
    !       write(10,FMT='(3F10.5)') (x(3,p,i),p=1,3)
    !     end do

    return

  end subroutine findPpalAxes


  !----------------------------------------------------------------------+
  !> Routine to rotate the multipoles to the orientation of the molecules.
  !! @param[in] d0     : ?
  !! @param[in] q0     : ?
  !! @param[in] o0     : ?
  !! @param[in] h0     : ?
  !! @param[out] dpole : The array is assumed to be zeroized.
  !! @param[out] qpole : ?
  !! @param[out] opole : ?
  !! @param[out] hpole : ?
  !! @param[out] nM    : The number of molecules.
  !! @param[out] x     : ?
  subroutine rotatePoles(d0, q0, o0, h0, dpole, qpole, opole, hpole, nM, x)
    implicit none
    real(dp), intent(in) :: d0(:) !JÖ assumed shape can be faster
    real(dp), intent(in) :: q0(:,:)
    real(dp), intent(in) :: o0(:,:,:)
    real(dp), intent(in) :: h0(:,:,:,:)
    integer, intent(in) :: nM
    real(dp), intent(in) :: x(:,:,:)

    real(dp), intent(out) :: dpole(:,:)
    real(dp), intent(out) :: qpole(:,:,:)
    real(dp), intent(out) :: opole(:,:,:,:)
    real(dp), intent(out) :: hpole(:,:,:,:,:)



!    real(dp), intent(in) :: d0(3)
!    real(dp), intent(in) :: q0(3,3)
!    real(dp), intent(in) :: o0(3,3,3)
!    real(dp), intent(in) :: h0(3,3,3,3)
!    integer, intent(in) :: nM
!    real(dp), intent(in) :: x(3,3,nM)
    ! ----------------------------------------

!    real(dp), intent(out) :: dpole(3,nM)
!    real(dp), intent(out) :: qpole(3,3,nM)
!    real(dp), intent(out) :: opole(3,3,3,nM)
!    real(dp), intent(out) :: hpole(3,3,3,3,nM)


    ! local variables.
    integer i, j, k, l, ii, jj, kk, ll, m !jö , save ::

    ! Zeroise.
!JÖ    hpole(1:3, 1:3, 1:3, 1:3, 1:nM) = 0.0_dp
!JÖ    opole(1:3, 1:3, 1:3, 1:nM) = 0.0_dp
!JÖ    qpole(1:3, 1:3, 1:nM) = 0.0_dp
!JÖ    dpole(1:3, 1:nM) = 0.0_dp

    hpole = 0 !JÖ this is enough 
    opole = 0
    qpole = 0
    dpole = 0
    ! Do the calculation.
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
             dpole(l,m) = dpole(l,m) + x(l,ll,m) * d0(ll)
          end do
       end do
    end do

    return

  end subroutine rotatePoles

  !----------------------------------------------------------------------+
  !     Routine to rotate the multipoles to the orientation of the       |
  !     molecule                                                         |
  !----------------------------------------------------------------------+
  subroutine setUnpolPoles(dpole, qpole, dpole0, qpole0, nM)

    implicit none

    
    real(dp), intent(inout) :: dpole(:,:) !JÖ 
    real(dp), intent(inout) :: qpole(:,:,:)
    real(dp), intent(in) :: dpole0(:,:)
    real(dp), intent(in) :: qpole0(:,:,:)
    integer , intent(in) :: nM
    
    integer i, j, k
    
!JÖ    real(dp) dpole(3,maxCoo/3), qpole(3,3,maxCoo/3)
!JÖ    real(dp) dpole0(3,maxCoo/3), qpole0(3,3,maxCoo/3)

    
    dPole = dpole0 !JÖ 
    qpole = qpole0 !JÖ

!JÖ    do i = 1, nM 
!JÖ       do j = 1, 3
!JÖ          dPole(j,i) = dpole0(j,i)
!JÖ          do k = 1, 3
!JÖ             qpole(k,j,i) = qpole0(k,j,i)
!JÖ          end do
!JÖ       end do
!JÖ    end do
    return

  end subroutine setUnpolPoles

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
  !     Add all the fields                                               |
  !----------------------------------------------------------------------+
  subroutine addFields(eH, eD, eT, nM)

    implicit none
!    integer i, j, nM
!    real(dp) eH(3,maxCoo/3), eD(3,maxCoo/3)
!    real(dp) eT(3,maxCoo/3)

    integer , intent(in)   :: nM
    real(dp), intent(in)   :: eH(:,:)
    real(dp), intent(in)   :: eD(:,:)
    real(dp), intent(out) :: eT(:,:)
    
    integer i, j
   
    eT = eH + eD

!JÖ    do i = 1, nM
!JÖ       do j = 1, 3
!JÖ          eT(j,i) = eH(j,i) + eD(j,i)
!JÖ       end do
!JÖ    end do
    return

  end subroutine addFields

  !----------------------------------------------------------------------+
  !     Add the derivative of all the fields                             |
  !----------------------------------------------------------------------+
  subroutine addDfields(dEhdr, dEddr, dEtdr, nM)

    implicit none
    integer, intent(in) :: nM !JÖ
!    real(dp) dEhdr(3,3,maxCoo/3)
!    real(dp) dEtdr(3,3,maxCoo/3), dEddr(3,3,maxCoo/3)
    real(dp), intent(out) :: dEtdr(:,:,:)
    real(dp), intent(in)  :: dEhdr(:,:,:)
    real(dp), intent(in)  :: dEddr(:,:,:)

    integer i, j, k !JÖ, nM
    
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
    
    !JÖ:
    dEtdr = dEhdr + dEddr
    !do i = 1, nM
    do j = 2, 3
       do k = 1, j-1
          dEtdr(j,k,:) = dEtdr(k,j,:)
       end do
    end do
    !end do
    
    
    return

  end subroutine addDfields

  !-----------------------------------------------------------------------
  !        subroutine applyPBC(r, sep, j)
  !
  !                ! not used? I commented it out
  !                ! to use: uncomment, add lattice to argument list and
  !                ! change ax, ay, az to lattice(1), lattice(2), lattice(3)
  !
  !                implicit real(dp) (a-h,o-z)
  !                real(dp) r(3), a(3), sep(3)
  !                integer j
  !
  !                !     Size of the simulation cell
  !                a(1) = ax/2.0d0
  !                a(2) = ay/2.0d0
  !                a(3) = az/2.0d0
  !
  !                do i = 1, 3
  !                        !         if (j*sep(i) .gt. a(i)) then
  !                        if (r(i) .gt. a(i)) then
  !                                r(i) = r(i) - 2.0d0 * a(i)
  !                                !         elseif (j*sep(i) .lt. -a(i)) then
  !                        elseif (r(i) .lt. -a(i)) then
  !                                r(i) = r(i) + 2.0d0 * a(i)
  !                        end if
  !                end do
  !
  !                j = 0
  !
  !                return
  !
  !        end subroutine applyPBC

  !-----------------------------------------------------------------------
!  pure subroutine SF(r, swFunc)
!
!    implicit none
!    real(dp), intent(in) :: r !JÖ
!    real(dp), intent(out) :: swFunc !JÖ
!    real(dp) x, x2, x3, rSW, rCut, dr     !JÖ r, swFunc, 
!    !JÖ , save ::
!!JÖ    real(dp), save :: rL1, rL2, rH1, rH2
!    real(dp), parameter :: rL1 = 0.0_dp
!    real(dp), parameter :: rL2 = 9.0_dp
!    real(dp), parameter :: rH1 = 5.0_dp
!    real(dp), parameter :: rH2 = 11.0_dp
!    
!    !      data rL1, rH1, rL2, rH2 / 1.5d0, 2.7d0, 8.d0, 9.d0 /
!!JÖ    data rL1, rH1, rL2, rH2 / 0.0_dp, 5.0_dp, 9.0_dp, 11.0_dp /
!    !      data rL1, rH1, rL2, rH2 / 0.0_dp, 5.0_dp, 11.0_dp, 13.0_dp /
!!JÖ     save
!
!
!    !     Kroes
!    !$$$      rSW = 8.d0
!    !$$$      rCut = 9.d0
!    !$$$
!    !$$$      x = (r - rSW)/(rCut - rSW)
!    !$$$
!    !$$$      if (r. lt. rSW) then
!    !$$$         swFunc = 1.0d0
!    !$$$      else if(r .lt. rCut) then
!    !$$$         x2 = x * x
!    !$$$         x3 = x2 * x
!    !$$$         swFunc = 1.d0 + x3 * (-6.d0 * x2 + 15.d0 * x - 10.d0)
!    !$$$      else
!    !$$$         swFunc = 0.0d0
!    !$$$      end if
!
!    if ((r .ge. rH2) .or. (r .le. rL1)) then
!       swFunc = 0.0_dp
!    else if ((r .ge. rH1) .and. (r .le. rL2)) then
!       swFunc = 1.0_dp
!    else if (r .lt. rH1) then
!
!       call tang_toenniesN(r, swFunc, 6)
!
!       !$$$         x = 1.d0 - (r - rL1)/(rH1 - rL1)
!       !$$$         x2 = x * x
!       !$$$         x3 = x2 * x
!       !$$$         swFunc = 1.d0 + x3 * (-6.d0 * x2 + 15.d0 * x - 10.d0)
!    else
!       x = (r - rL2)/(rH2 - rL2)
!       x2 = x * x
!       x3 = x2 * x
!       swFunc = 1.0_dp + x3 * (-6.0_dp * x2 + 15.0_dp * x - 10.0_dp)
!    end if
!
!    !      swFunc = 1.0_dp
!
!    return
!
!  end subroutine SF
!
!  !-----------------------------------------------------------------------
!  pure subroutine SFdsf(r, swFunc, dSdr)!JÖ pure
!
!    implicit none
!    real(dp), intent(in) :: r !JÖ stated intents
!    real(dp), intent(out) :: swFunc, dSdr !JÖ
!    real(dp) :: x, x2, x3, rSW, rCut, dr
!!JÖ pure    real(dp) :: rL1, rL2, rH1, rH2  !JÖ pure: ^, dr
!    
!    real(dp), parameter :: rL1 = 0.0_dp   !JÖ pure
!    real(dp), parameter :: rH1 = 5.0_dp   !JÖ pure
!    real(dp), parameter :: rL2 = 9.0_dp   !JÖ pure
!    real(dp), parameter :: rH2 = 11.0_dp  !JÖ pure
!!, save
!!, save
!    !      data rL1, rH1, rL2, rH2 / 1.5d0, 2.7d0, 8.d0, 9.d0 /
!!JÖ pure    data rL1, rH1, rL2, rH2 / 0.0_dp, 5.0_dp, 9.0_dp, 11.0_dp /
!    !      data rL1, rH1, rL2, rH2 / 0.0_dp, 5.0_dp, 11.0_dp, 13.0_dp /
!    !JÖ save
!
!    !     Kroes
!    !$$$      rSW = 8.d0
!    !$$$      rCut = 9.d0
!
!
!    !$$$      x = (r - rSW)/(rCut - rSW)
!    !$$$
!    !$$$      if (r. lt. rSW) then
!    !$$$         swFunc = 1.0d0
!    !$$$         dSdr = 0.0d0
!    !$$$      else if(r .lt. rCut) then
!    !$$$         x2 = x * x
!    !$$$         x3 = x2 * x
!    !$$$         swFunc = 1.d0 + x3 * (-6.d0 * x2 + 15.d0 * x - 10.d0)
!    !$$$         dSdr = 30.d0 * x2 * (- x2 + 2.d0 * x - 1.d0) / (rCut-rSW)
!    !$$$      else
!    !$$$         swFunc = 0.0d0
!    !$$$         dSdr = 0.0d0
!    !$$$      end if
!
!
!    if ((r .ge. rH2) .or. (r .le. rL1)) then
!       swFunc = 0.0_dp
!       dSdr = 0.0_dp
!    else if ((r .ge. rH1) .and. (r .le. rL2)) then
!       swFunc = 1.0_dp
!       dSdr = 0.0_dp
!    else if (r .lt. rH1) then
!
!       !c         call switchCloseDsf(r, swFunc, dSdr)
!       call tang_toenniesNdF(r, swFunc, dSdr, 6)
!
!       !$$$         x = 1.d0 - (r - rL1)/(rH1 - rL1)
!       !$$$         dr = - 1.d0 / (rH1 - rL1)
!       !$$$         x2 = x * x
!       !$$$         x3 = x2 * x
!       !$$$         swFunc = 1.d0 + x3 * (-6.d0 * x2 + 15.d0 * x - 10.d0)
!       !$$$         dSdr = 30.d0 * x2 * (- x2 + 2.d0 * x - 1.d0) * dr
!    else
!       x = (r - rL2)/(rH2 - rL2)
!       dr = 1.0_dp / (rH2 - rL2)
!       x2 = x * x
!       x3 = x2 * x
!       swFunc = 1.0_dp + x3 * (-6.0_dp * x2 + 15.0_dp * x - 10.0_dp)
!       dSdr = 30.0_dp * x2 * (- x2 + 2.0_dp * x - 1.0_dp) * dr
!    end if
!
!    !      swFunc = 1.d0
!    !      dSdr = 0.d0
!
!    return
!
!  end subroutine SFdsf
!

end module molecProperties
