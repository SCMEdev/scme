module ps_dms
!use constants, only:a0_A
use data_types, only:a0_A, dp

implicit none

integer, parameter, dimension(84) :: idxD1 = [&
       1, 1, 1, 2, 1, 1, 1, 2, 2, 3, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4,&
       1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 5, 1, 1, 1, 1, 1,&
       1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5, 5, 6, 1, 1, 1, 1,&
       1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5,&
       5, 6, 6, 7]

integer, parameter, dimension(84) :: idxD2 = [&
       1, 1, 2, 1, 1, 2, 3, 1, 2, 1, 1, 2, 3, 4, 1, 2, 3, 1, 2, 1,&
       1, 2, 3, 4, 5, 1, 2, 3, 4, 1, 2, 3, 1, 2, 1, 1, 2, 3, 4, 5,&
       6, 1, 2, 3, 4, 5, 1, 2, 3, 4, 1, 2, 3, 1, 2, 1, 1, 2, 3, 4,&
       5, 6, 7, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 1, 2, 3, 4, 1, 2,&
       3, 1, 2, 1]


integer, parameter, dimension(84) :: idxD3 = [&
       1, 2, 1, 1, 3, 2, 1, 2, 1, 1, 4, 3, 2, 1, 3, 2, 1, 2, 1, 1,&
       5, 4, 3, 2, 1, 4, 3, 2, 1, 3, 2, 1, 2, 1, 1, 6, 5, 4, 3, 2,&
       1, 5, 4, 3, 2, 1, 4, 3, 2, 1, 3, 2, 1, 2, 1, 1, 7, 6, 5, 4,&
       3, 2, 1, 6, 5, 4, 3, 2, 1, 5, 4, 3, 2, 1, 4, 3, 2, 1, 3, 2,&
       1, 2, 1, 1]

real(dp), parameter, dimension(84) :: coefD = [&
      -2.1689686086730d-03, 1.4910379754728d-02, 5.3546078430060d-02,&
      -7.4055995388666d-02,-3.7764333017616d-03, 1.4089887256484d-01,&
      -6.2584207687264d-02,-1.1260393113022d-01,-5.7824159269319d-02,&
       1.4360743650655d-02,-1.5469680141070d-02,-1.3036350092795d-02,&
       2.7515837781556d-02, 1.4098478875076d-01,-2.7663168397781d-02,&
      -5.2378176254797d-03,-1.0237198381792d-02, 8.9571999265473d-02,&
       7.2920263098603d-03,-2.6873260551686d-01, 2.0220870325864d-02,&
      -7.0764766270927d-02, 1.2140640273760d-01, 2.0978491966341d-02,&
      -1.9443840512668d-01, 4.0826835370618d-02,-4.5365190474650d-02,&
       6.2779900072132d-02,-1.3194351021000d-01,-1.4673032718563d-01,&
       1.1894031277247d-01,-6.4952851564679d-03, 8.8503610374493d-02,&
       1.4899437409291d-01, 1.3962841511565d-01,-2.6459446720450d-02,&
      -5.0128914532773d-02, 1.8329676428116d-01,-1.5559089125095d-01,&
      -4.0176879767592d-02, 3.6192059996636d-01, 1.0202887240343d-01,&
       1.9318668580051d-01,-4.3435977107932d-01,-4.2080828803311d-02,&
       1.9144626027273d-01,-1.7851138969948d-01, 1.0524533875070d-01,&
      -1.7954071602185d-02, 5.2022455612120d-02,-2.8891891146828d-01,&
      -4.7452036576319d-02,-1.0939400546289d-01, 3.5916564473568d-01,&
      -2.0162789820172d-01,-3.5838629543696d-01, 5.6706523551202d-03,&
       1.3849337488211d-01,-4.1733982195604d-01, 4.1641570764241d-01,&
      -1.2243429796296d-01, 4.7141730971228d-02,-1.8224510249551d-01,&
      -1.8880981556620d-01,-3.1992359561800d-01,-1.8567550546587d-01,&
       6.1850530431280d-01,-6.1142756235141d-02,-1.6996135584933d-01,&
       5.4252879499871d-01, 6.6128603899427d-01, 1.2107016404639d-02,&
      -1.9633639729189d-01, 2.7652059420824d-03,-2.2684111109778d-01,&
      -4.7924491598635d-01, 2.4287790137314d-01,-1.4296023329441d-01,&
       8.9664665907006d-02,-1.4003228575602d-01,-1.3321543452254d-01,&
      -1.8340983193745d-01, 2.3426707273520d-01, 1.5141050914514d-01]


!real*8, parameter :: b2A     = 0.52917721067000000000d0 !nist value bohr in ångström 
!                               
!real*8, parameter :: A2b     = 1.88972612545782819321d0 !Ångström to Bohr
!real*8, parameter :: au2deb  = 2.541746230211d0 !a.u. to debye

real(dp), parameter :: reoh = 0.958649d0
real(dp), parameter :: b1D =  1.0d0
real(dp), parameter :: a   =  0.2999d0
real(dp), parameter :: b   = -0.6932d0
real(dp), parameter :: c0  =  1.0099d0
real(dp), parameter :: c1  = -0.1801d0
real(dp), parameter :: c2  =  0.0892d0
real(dp), parameter :: costhe = -0.24780227221366464506d0

private
public vibdms

contains !//////////////////////////////////////////////////////////////

   subroutine vibdms(rr,dd)
    real(dp), intent(in)  :: rr(3,3)!x(3,3) x(xyz,OHH)
    real(dp), intent(out) :: dd(3)
    real(dp) :: efac, xpow(3,7) !,x1,x2,x3 !7 should be enough
    real(dp) :: p1, p2, pl1, pl2,x(3)
    real(dp) :: v1(3), v2(3), vhh(3), r1, r2, rhh, costh, pc0
    integer inI, inJ, inK,j
    
    
    !for derivative:
    real(dp) :: dp1dr1  
    real(dp) :: dp1dr2  
    real(dp) :: dp1dcabc
    real(dp) :: dp2dr1  
    real(dp) :: dp2dr2  
    real(dp) :: dp2dcabc
    
    integer inI_1, inJ_1, inK_1
    
    real(dp) :: dpc0dr1  
    real(dp) :: dpc0dr2  
    real(dp) :: dpc0dcabc
    real(dp) :: defacdr1 
    real(dp) :: defacdr2 
    
    
    real(dp) :: f1q1r13
    real(dp) :: f1q1r23
    real(dp) :: f2q1r23
    real(dp) :: f2q1r13
    real(dp) :: f1q2r13
    real(dp) :: f1q2r23
    real(dp) :: f2q2r23
    real(dp) :: f2q2r13
    
    real(dp) :: GRADQ(3,3,3)
    real(dp) :: GRADH1(3,3)
    real(dp) :: GRADH2(3,3)
    real(dp) :: GRADO (3,3)
    
      
      v1  = ( rr(:,1) - rr(:,3) )!*A2b !rr(i,2,:) - rr(i,1,:) !H1-O
      v2  = ( rr(:,2) - rr(:,3) )!*A2b !rr(i,3,:) - rr(i,1,:) !H2-O
      v12 = ( rr(:,1) - rr(:,2) )!H1-H2
      !r1  = dsqrt(sum( v1*v1 ))
      !r2  = dsqrt(sum( v2*v2 ))
      !costh = sum(v1*v2) / (r1*r2) !=cos(HOH) 
      
      
      r1  = dsqrt( v1(1)*v1(1) + v1(2)*v1(2) + v1(3)*v1(3) )
      r2  = dsqrt( v2(1)*v2(1) + v2(2)*v2(2) + v2(3)*v2(3) )
      costh = ( v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3) ) /(r1*r2)
      
      r12 = dsqrt( v12(1)*v12(1) + v12(2)*v12(2) + v12(3)*v12(3) )
      
      efac = dexp(-b1D*( (r1-reoh)**2 + (r2-reoh)**2 ) )
      
      x = [ (r1-reoh)/reoh , (r2-reoh)/reoh , (costh - costhe) ]
      
      xpow(:,1) = 1d0
      do j = 2,7 ! powers of the x-vector
          xpow(:,j) = xpow(:,j-1)*x(:)
      enddo
      
      ! Calculate the dipole moment
      
      p1=0d0
      p2=0d0
      
      do j=2,84 !(size_t j = 1; j < 84; ++j) {
          inI = idxD1(j)
          inJ = idxD2(j)
          inK = idxD3(j)
          
          p1 = p1 + coefD(j)*xpow(1,inI)*xpow(2,inJ)*xpow(3,inK)
          p2 = p2 + coefD(j)*xpow(1,inJ)*xpow(2,inI)*xpow(3,inK)
          
          inI_1 = inI - 1
          inJ_1 = inJ - 1
          inK_1 = inK - 1
          
          
          !  derivative:
		  dp1dr1   = dp1dr1   +  coefD(j)*(inI_1)*xpow(1,inI_1)*xpow(2,inJ)*  xpow(3,inK)
		  dp1dr2   = dp1dr2   +  coefD(j)*(inJ_1)*xpow(1,inI)  *xpow(2,inJ_1)*xpow(3,inK)
		  dp1dcabc = dp1dcabc +  coefD(j)*(inK_1)*xpow(1,inI)  *xpow(2,inJ)*  xpow(3,inK_1)
		  dp2dr1   = dp2dr1   +  coefD(j)*(inJ_1)*xpow(1,inJ_1)*xpow(2,inI)*  xpow(3,inK)
		  dp2dr2   = dp2dr2   +  coefD(j)*(inI_1)*xpow(1,inJ)*  xpow(2,inI_1)*xpow(3,inK)
		  dp2dcabc = dp2dcabc +  coefD(j)*(inK_1)*xpow(1,inJ)*  xpow(2,inI)*  xpow(3,inK_1)
      enddo
      
      xx = a0_A
      xx2 = xx**2
      
      pl1 = costh;
      pl2 = 0.5d0*(3d0*pl1*pl1 - 1d0)
      pc0 = a*( r1**b + r2**b )*(c0 + pl1*c1 + pl2*c2)
      
      p1 = coefD(1) + p1*efac + pc0*xx !; // q^H1 in TTM2-F
      p2 = coefD(1) + p2*efac + pc0*xx !; // q^H2 paper
          
      !Dipole moment:         
      dd = ( p1*v1 + p2*v2 )!/a0_A
      
      
      !Partial charges:
      !dd(1) = -(p1 + p2);  
      !dd(2) = p1;          
      !dd(3) = p2;          
      
      
      ! derivative: 
      
      dp1dr1 = dp1dr1 * xx/reoh  
      dp1dr2 = dp1dr2 * xx/reoh 
      dp2dr1 = dp2dr1 * xx/reoh 
      dp2dr2 = dp2dr2 * xx/reoh 
      
      dpc0dr1   = a*b*r1**(b-1)*(c0 + pl1*c1 + pl2*c2)*xx2
      dpc0dr2   = a*b*r2**(b-1)*(c0 + pl1*c1 + pl2*c2)*xx2
      dpc0dcabc = a*( r1**b + r2**b )*(c1 + 0.5d0*(6.0d0*pl1)*c2)*xx
      
      defacdr1 = -2.d0*b1D*(r1 - reoh)*efac*xx
      defacdr2 = -2.d0*b1D*(r2 - reoh)*efac*xx
      
      dp1dr1   = dp1dr1   * efac + p1*defacdr1 + dpc0dr1
      dp1dr2   = dp1dr2   * efac + p1*defacdr2 + dpc0dr2
      dp1dcabc = dp1dcabc * efac + dpc0dcabc
      dp2dr1   = dp2dr1   * efac + p2*defacdr1 + dpc0dr1
      dp2dr2   = dp2dr2   * efac + p2*defacdr2 + dpc0dr2
      dp2dcabc = dp2dcabc * efac + dpc0dcabc
      
      dp1dr1 = dp1dr1 /xx
      dp1dr2 = dp1dr2 /xx
      dp2dr1 = dp2dr1 /xx
      dp2dr2 = dp2dr2 /xx
      
      ! Above derivatives w.r.t internal coordinates
      
      
      !dd = p1*v1 + p2*v2 
      !dd = -(qO+q2)*v1 + -(qO+q2)*v2
      
      
      
      ! derivative w.r.t, to cartesian coordinates below ? 
      ! gotta derive this!
 !------------------------------     
      ! 
      f1q1r13 = (dp1dr1 - (dp1dcabc*costh/r1))/r1;
      f1q1r23 = dp1dcabc/(r1*r2);
      f2q1r23 = (dp1dr2 - (dp1dcabc*costh/r2))/r2;
      f2q1r13 = dp1dcabc/(r2*r1);
      f1q2r13 = (dp2dr1 - (dp2dcabc*costh/r1))/r1;
      f1q2r23 = dp2dcabc/(r1*r2);
      f2q2r23 = (dp2dr2 - (dp2dcabc*costh/r2))/r2;
      f2q2r13 = dp2dcabc/(r2*r1);     
      
      ! first index is atom w.r.t. to which the derivative is
      ! second index is the charge being differentiated
      
      
      
      ! H1: GRADQ(., 1, .)
      !gradient of charge h1(second index) wrt displacement of h1(first index)
      GRADQ(1, 1, :) = f1q1r13*v1(:) + f1q1r23*v2(:)
 
      !gradient of charge h1 wrt displacement of h2
      GRADQ(2, 1, :) = f2q1r13*v1(:) + f2q1r23*v2(:)
 
      !gradient of charge h1 wrt displacement of O
      GRADQ(3, 1, :) = -(GRADQ(1, 1, :) + GRADQ(2, 1, :))
      
      
      !H2: GRADQ(., 2, .)
      !gradient of charge h2 wrt displacement of h1
      GRADQ(1, 2, :) = f1q2r13*v1(:) + f1q2r23*v2(:)
 
      !gradient of charge h2 wrt displacement of h2
      GRADQ(2, 2, :) = f2q2r13*v1(:) + f2q2r23*v2(:)
 
      !gradient of charge h2 wrt displacement of O
      GRADQ(3, 2, :) = -(GRADQ(1, 2, :) + GRADQ(2, 2, :))
 
      
      !O: GRADQ(., 3, .)
      !gradient of charge O wrt displacement of h1
      GRADQ(1, 3, :) = -(GRADQ(1, 1, :) + GRADQ(1, 2, :))
 
      !gradient of charge O wrt displacement of h2
      GRADQ(2, 3, :) = -(GRADQ(2, 1, :) + GRADQ(2, 2, :))
 
      !gradient of charge O wrt displacement of O
      GRADQ(3, 3, :) = -(GRADQ(3, 1, :) + GRADQ(3, 2, :))
      
      !do i = 1, 3
      !   do j = 1, 3
      !      do k = 1, 3
      !         dq3(  1 + (k-1) + 3*(j-1) + 9*(i-1)  ) = GRADQ(i,j,k) 
      !      enddo
      !   enddo
      !enddo
      
      
      
!      ! H1: GRADQ(., 1, .)
!      !gradient of charge h1(second index) wrt displacement of h1(first index)
!      GRADH1(1, :) = f1q1r13*v1(:) + f1q1r23*v2(:)
! 
!      !gradient of charge h1 wrt displacement of h2
!      GRADH1(2, :) = f2q1r13*v1(:) + f2q1r23*v2(:)
! 
!      !gradient of charge h1 wrt displacement of O
!      GRADH1(3, :) = -(GRADH1(1, :) + GRADH1(2, :))
!      
!      
!      !H2: GRADQ(., 2, .)
!      !gradient of charge h2 wrt displacement of h1
!      GRADH2(1, :) = f1q2r13*v1(:) + f1q2r23*v2(:)
! 
!      !gradient of charge h2 wrt displacement of h2
!      GRADH2(2, :) = f2q2r13*v1(:) + f2q2r23*v2(:)
! 
!      !gradient of charge h2 wrt displacement of O
!      GRADH2(3, :) = -(GRADH2(1, :) + GRADH2(2, :))
! 
!      
!      !O: GRADQ(., 3, .)
!      !gradient of charge O wrt displacement of h1
!   !    GRADO(1, :) = -(GRADH1(1, :) + GRADH2(1, :))
!   !  
!   !    !gradient of charge O wrt displacement of h2
!   !    GRADO(2, :) = -(GRADH1(2, :) + GRADH2(2, :))
!   !    
!   !    !gradient of charge O wrt displacement of O
!   !    GRADO(3, :) = -(GRADH1(3, :) + GRADH2(3, :))
!
!       GRADO(:, :) = -(GRADH1(:, :) + GRADH2(:, :))
!       
!       
!      
!      real(dp) :: DH1_H1(3) 
!      real(dp) :: DH2_H1(3) 
!      real(dp) :: DO_H1(3)  
!      real(dp) :: DH1_H2(3) 
!      real(dp) :: DH2_H2(3) 
!      real(dp) :: DO_H2(3)  
!      real(dp) :: DH1_O(3) 
!      real(dp) :: DH2_O(3) 
!      real(dp) :: DO_O(3) 
!      
!      
!       
!       !w.r.t H1:
!      Dq1_H1(:) = f1q1r13*v1(:) + f1q1r23*v2(:)
!      Dq2_H1(:) = f1q2r13*v1(:) + f1q2r23*v2(:)
!      DqO_H1(:)  = -( Dq1_H1(:) + Dq2_H1(:) )
!      
!      !w.r.t H2:
!      Dq1_H2(:) = f2q1r13*v1(:) + f2q1r23*v2(:)
!      Dq2_H2(:) = f2q2r13*v1(:) + f2q2r23*v2(:)
!      DqO_H2(:)  = -( Dq1_H2(:) + Dq2_H2(:) )
!      
!      !w.r.t O:
!      Dq1_O(:) = -( Dq1_H1(:) + Dq1_H2(:) ) 
!      Dq2_O(:) = -( Dq2_H1(:) + Dq2_H2(:) )
!      DqO_O (:) = -( Dq1_O(:) + Dq2_O(:) )
!      
!      Ddx_H1 = 
!      
   end subroutine
end module

      !Dipole moment:         
      !dd = ( p1*v1 + p2*v2 )!/a0_A
      
      
      !Partial charges:
      !q1 = p1;          
      !q2 = p2; 
      !qO = -(p1 + p2);  
      
!      dd = p1*v1 + p2*v2 
!      dd = -(qO+q2)*v1 + -(qO+q2)*v2
      
      
      
      
      
      !r1  = dsqrt( v1(1)*v1(1) + v1(2)*v1(2) + v1(3)*v1(3) )
      !r2  = dsqrt( v2(1)*v2(1) + v2(2)*v2(2) + v2(3)*v2(3) )
      !costh = ( v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3) ) /(r1*r2)
      
      !x1 = (r1 - reoh)/reoh
      !x2 = (r2 - reoh)/reoh
      !x3 = costh - costhe
      !x = [x1,x2,x3]
      
