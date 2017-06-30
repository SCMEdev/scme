subroutine test_xyz()
type(xyz) r(2)
r(1)%x = [1.2_dp, 3.4_dp,5.6_dp]
r(2)%x = [4.2_dp, 6.4_dp,1.6_dp]
print*, r(1)%x
print*, r(2)%x
end subroutine

!type(hho) :: w(nM), aforces(nM)
! call torqueOnAtoms2(tau(:,m), rCM(:,m), w(m)%a, fatom(m)%a

subroutine torqueOnAtoms2(tau,rCM,rhho,f)
 integer, parameter :: a_m = 3
 type(xyz), intent(in) :: tau, rCM
 type(xyz), intent(in) :: rhho(a_m)
 type(xyz), intent(inout) :: f(a_m)
 type(xyz) :: r(a_m), txr(a_m)
 real(dp) :: I(a_m), Itot, t2
 !real(dp), dimension(3) :: txr1, txr2, txro, r1, r2, ro
 !real(dp) :: I1, I2, Io, Itot, t2
 real(dp), parameter :: m(a_m) = [1.0_dp, 1.0_dp, 16.0_dp] !are these in the right units???
 integer ii, a
 
   t2 = norm_2(tau%x)
   
   do a = 1,a_m !atoms hho
      r(a)%x = rhho(a)%x - rCM%x
      txr(a)%x = cross(tau%x,r(a)%x)
      I(a) = m(a)/t2*norm_2(txr(a)%x)
   enddo
   
   Itot = sum(I)
   
   do a = 1,a_m
      f(a)%x = m(a)/Itot*txr(a)%x
   enddo
   
end subroutine

subroutine torqueOnAtoms3(tau,rCM,mol,f)
 integer, parameter :: aim = 3
 real(dp), intent(in) :: tau(3), rCM(3)
 real(dp), intent(in) :: mol(3,aim)
 real(dp), intent(inout) :: f(3,aim)
 real(dp) :: r(3), txr(3,aim)
 real(dp) :: I(aim), Itot, t2
 real(dp), parameter :: m(aim) = [1.0_dp, 1.0_dp, 16.0_dp] !are these in the right units???
 integer ii, a
 
   t2 = norm_2(tau)
   
   do a = 1,aim !atoms hho
      r = mol(:,a) - rCM
      txr(:,a) = cross(tau,r)
      I(a) = m(a)/t2*norm_2(txr(:,a))
   enddo
   
   Itot = sum(I)
   
   do a = 1,aim
      f(:,a) = m(a)/Itot*txr(:,a)
   enddo
   
end subroutine

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
  pure subroutine SF(r, swFunc)

    implicit none
    real(dp), intent(in) :: r !JÖ
    real(dp), intent(out) :: swFunc !JÖ
    real(dp) x, x2, x3, rSW, rCut, dr     !JÖ r, swFunc, 
    !JÖ , save ::
!JÖ    real(dp), save :: rL1, rL2, rH1, rH2
    real(dp), parameter :: rL1 = 0.0_dp
    real(dp), parameter :: rL2 = 9.0_dp
    real(dp), parameter :: rH1 = 5.0_dp
    real(dp), parameter :: rH2 = 11.0_dp
    
    !      data rL1, rH1, rL2, rH2 / 1.5d0, 2.7d0, 8.d0, 9.d0 /
!JÖ    data rL1, rH1, rL2, rH2 / 0.0_dp, 5.0_dp, 9.0_dp, 11.0_dp /
    !      data rL1, rH1, rL2, rH2 / 0.0_dp, 5.0_dp, 11.0_dp, 13.0_dp /
!JÖ     save


    !     Kroes
    !$$$      rSW = 8.d0
    !$$$      rCut = 9.d0
    !$$$
    !$$$      x = (r - rSW)/(rCut - rSW)
    !$$$
    !$$$      if (r. lt. rSW) then
    !$$$         swFunc = 1.0d0
    !$$$      else if(r .lt. rCut) then
    !$$$         x2 = x * x
    !$$$         x3 = x2 * x
    !$$$         swFunc = 1.d0 + x3 * (-6.d0 * x2 + 15.d0 * x - 10.d0)
    !$$$      else
    !$$$         swFunc = 0.0d0
    !$$$      end if

    if ((r .ge. rH2) .or. (r .le. rL1)) then
       swFunc = 0.0_dp
    else if ((r .ge. rH1) .and. (r .le. rL2)) then
       swFunc = 1.0_dp
    else if (r .lt. rH1) then

       call tang_toenniesN(r, swFunc, 6)

       !$$$         x = 1.d0 - (r - rL1)/(rH1 - rL1)
       !$$$         x2 = x * x
       !$$$         x3 = x2 * x
       !$$$         swFunc = 1.d0 + x3 * (-6.d0 * x2 + 15.d0 * x - 10.d0)
    else
       x = (r - rL2)/(rH2 - rL2)
       x2 = x * x
       x3 = x2 * x
       swFunc = 1.0_dp + x3 * (-6.0_dp * x2 + 15.0_dp * x - 10.0_dp)
    end if

    !      swFunc = 1.0_dp

    return

  end subroutine SF

  !-----------------------------------------------------------------------
  pure subroutine SFdsf(r, swFunc, dSdr)!JÖ pure

    implicit none
    real(dp), intent(in) :: r !JÖ stated intents
    real(dp), intent(out) :: swFunc, dSdr !JÖ
    real(dp) :: x, x2, x3, rSW, rCut, dr
!JÖ pure    real(dp) :: rL1, rL2, rH1, rH2  !JÖ pure: ^, dr
    
    real(dp), parameter :: rL1 = 0.0_dp   !JÖ pure
    real(dp), parameter :: rH1 = 5.0_dp   !JÖ pure
    real(dp), parameter :: rL2 = 9.0_dp   !JÖ pure
    real(dp), parameter :: rH2 = 11.0_dp  !JÖ pure
!, save
!, save
    !      data rL1, rH1, rL2, rH2 / 1.5d0, 2.7d0, 8.d0, 9.d0 /
!JÖ pure    data rL1, rH1, rL2, rH2 / 0.0_dp, 5.0_dp, 9.0_dp, 11.0_dp /
    !      data rL1, rH1, rL2, rH2 / 0.0_dp, 5.0_dp, 11.0_dp, 13.0_dp /
    !JÖ save

    !     Kroes
    !$$$      rSW = 8.d0
    !$$$      rCut = 9.d0


    !$$$      x = (r - rSW)/(rCut - rSW)
    !$$$
    !$$$      if (r. lt. rSW) then
    !$$$         swFunc = 1.0d0
    !$$$         dSdr = 0.0d0
    !$$$      else if(r .lt. rCut) then
    !$$$         x2 = x * x
    !$$$         x3 = x2 * x
    !$$$         swFunc = 1.d0 + x3 * (-6.d0 * x2 + 15.d0 * x - 10.d0)
    !$$$         dSdr = 30.d0 * x2 * (- x2 + 2.d0 * x - 1.d0) / (rCut-rSW)
    !$$$      else
    !$$$         swFunc = 0.0d0
    !$$$         dSdr = 0.0d0
    !$$$      end if


    if ((r .ge. rH2) .or. (r .le. rL1)) then
       swFunc = 0.0_dp
       dSdr = 0.0_dp
    else if ((r .ge. rH1) .and. (r .le. rL2)) then
       swFunc = 1.0_dp
       dSdr = 0.0_dp
    else if (r .lt. rH1) then

       !c         call switchCloseDsf(r, swFunc, dSdr)
       call tang_toenniesNdF(r, swFunc, dSdr, 6)

       !$$$         x = 1.d0 - (r - rL1)/(rH1 - rL1)
       !$$$         dr = - 1.d0 / (rH1 - rL1)
       !$$$         x2 = x * x
       !$$$         x3 = x2 * x
       !$$$         swFunc = 1.d0 + x3 * (-6.d0 * x2 + 15.d0 * x - 10.d0)
       !$$$         dSdr = 30.d0 * x2 * (- x2 + 2.d0 * x - 1.d0) * dr
    else
       x = (r - rL2)/(rH2 - rL2)
       dr = 1.0_dp / (rH2 - rL2)
       x2 = x * x
       x3 = x2 * x
       swFunc = 1.0_dp + x3 * (-6.0_dp * x2 + 15.0_dp * x - 10.0_dp)
       dSdr = 30.0_dp * x2 * (- x2 + 2.0_dp * x - 1.0_dp) * dr
    end if

    !      swFunc = 1.d0
    !      dSdr = 0.d0

    return

  end subroutine SFdsf

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
    real(dp), intent(in) :: d0(3) !JÖ assumed shape can be faster
    real(dp), intent(in) :: q0(3,3)
    real(dp), intent(in) :: o0(3,3,3)
    real(dp), intent(in) :: h0(3,3,3,3)
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

