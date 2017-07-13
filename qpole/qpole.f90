PROGRAM qpole
use printer_mod, only:printer
use localAxes_mod,only: norm, norm_square
use data_types, only: dp, pi
use multipole_parameters, only: q0
implicit none
call main()
contains !///

SUBROUTINE main()
    integer,parameter :: xyz=3,dxyz=3,hho=3, rra=3
    real(dp) wg(xyz,hho), wi(rra)
    real(dp) dp1(xyz,hho),dp2(xyz,hho),dp3(xyz,hho),p(hho)
    real(dp) quadrupole(xyz,xyz)
    !real(dp) dQ_h1(xyz,xyz,xyz),dQ_h2(xyz,xyz,xyz),dQ_o(xyz,xyz,xyz)
    real(dp) ch1(xyz),ch2(xyz),co(xyz)
    real(dp) rMEC(xyz)
    real(dp) ang,rad,rh, oz, hx
    real(dp) :: cec(xyz,hho),cer2(hho), m(hho)
    
    real(dp) dQ(xyz,xyz,dxyz,hho)
    real(dp) dQ1(xyz,xyz,dxyz),dQ2(xyz,xyz,dxyz),dQo(xyz,xyz,dxyz)
    integer a,c

    
    m = [1d0,1d0,16d0]
    
    
    ang = 104.3475_dp
    rad = ang/180d0*pi
    rh = 0.958649_dp
    
    oz = cos(rad*0.5d0)*rh
    hx = sin(rad*0.5d0)*rh
    
    ch1= [-hx,0d0,0d0]
    ch2= [ hx,0d0,0d0]
    co = [0d0,0d0, oz ]
    
    wg(:,1) = ch1
    wg(:,2) = ch2
    wg(:,3) = co
    
    call printer(wg,'wg',2)
    
    call expansion_coordinates(wg,rMEC,cec,cer2)
    call quad_charges(cec,cer2, p, dp1, dp2, dp3)
    call quad_tensor(cec,cer2,p,quadrupole)
    
    call dQ_atomic(cec,cer2,p,quadrupole,m,dQ1,1)
    call dQ_atomic(cec,cer2,p,quadrupole,m,dQ2,2)
    call dQ_atomic(cec,cer2,p,quadrupole,m,dQo,3)
    
    dQ(:,:,:,1)=dq1
    dQ(:,:,:,2)=dq2
    dQ(:,:,:,3)=dqo
    
    
    do a = 1,hho
    do c = 1,xyz
    print*, 'dQ traces', matrix_trace(dQ(:,:,c,a))
    enddo
    enddo
    
call printer(cec,'cec',2)
call printer(cer2,'cer2',2)
    
call printer(p,'p',2)
call printer(dQ,'dQ',2)
call printer(quadrupole,'quadrupole',2)
print*, 'trace', matrix_trace(quadrupole)
    
call printer(q0,'q0',2)
print*, 'trace', matrix_trace(q0)
    
    
END SUBROUTINE main !////////////////


SUBROUTINE quad_charges(cec,cer2, p, dp1, dp2,dp3)
    integer,parameter    :: xyz=3,hho=3, rra=3
    real(dp),intent(in)  :: cec(xyz,hho),cer2(hho)
    real(dp),intent(out) :: p(hho), dp1(xyz,hho),dp2(xyz,hho),dp3(xyz,hho)
    real(dp) hx2
    
    !hx2=cec(1,1)**2
    !p(1) = ( q0(1,1) - q0(2,2) ) / (3*hx2)
    !p(2)=p(1)
    !p(3)= -2*( q0(1,1) - p(1)*(3*hx2-cer2(1)) )/cer2(3)
    
    p(1)=0.617152455762941d0
    p(2)=p(1)
    p(3)=-3.268376558140222d0
    
    dp1 = 0
    dp2 = 0
    dp3 = 0
END SUBROUTINE quad_charges!////////////



SUBROUTINE quad_tensor(cec,cer2,p,q)  !wg = global coordinates
    integer,parameter :: xyz=3, xy=2,hho=3, rra=3, dxyz=3
    real(dp),intent(in)  :: cec(xyz,hho),cer2(hho)
    real(dp),intent(in)  :: p(hho)
    real(dp),intent(out) :: q(xyz,xyz)
    integer a,c
    real(dp)  hx2
    
    
    q(:,:) = 0
    
    do a = 1,hho !add upp the quadrupole compnents from each partial charge
      !on diagonal (traceless)
      q(1,1) = q(1,1) + p(a)*(3*cec(1,a)**2 - cer2(a)) 
      q(2,2) = q(2,2) + p(a)*(3*cec(2,a)**2 - cer2(a))
      !off diagonal (symmetric)
      q(1,2) = q(1,2) + p(a)*cec(1,a)*cec(2,a)
      q(1,3) = q(1,3) + p(a)*cec(1,a)*cec(3,a)
      q(2,3) = q(2,3) + p(a)*cec(2,a)*cec(3,a)
    enddo
    !scale with 1/2 according to definition (See A. J. Stone)
    q(1,1) = q(1,1)*0.5d0
    q(2,2) = q(2,2)*0.5d0
    q(3,3) = - ( q(1,1) + q(2,2) ) !traceless
    !q(3,3) = q(3,3)*0.5d0
    
    !scale with 3/2 according to definition
    q(1,2) = q(1,2)*1.5d0
    q(1,3) = q(1,3)*1.5d0
    q(2,3) = q(2,3)*1.5d0
    
    ! Symmetric matrix
    q(2,1) = q(1,2)
    q(3,1) = q(1,3)
    q(3,2) = q(2,3)
    
    
END SUBROUTINE



SUBROUTINE dQ_atomic(cec,cer2,p,q,m,dQ,a)
    integer,parameter :: xyz=3,hho=3, dxyz=3, xy=2
    real(dp),intent(in)  :: cec(xyz,hho),cer2(hho), p(hho), q(xyz,xyz), m(hho)
    real(dp),intent(out) ::dQ(xyz,xyz,dxyz)
    !real(dp) dQ_h1(xyz,xyz,dxyz), dQ_h2(xyz,xyz,dxyz), dQ_o(xyz,xyz,dxyz)
    integer, intent(in) :: a
    integer c,d,od
    real(dp) mM(hho), mM1(hho)
    
    integer k, c1,c2 !,cs(xyz,dxyz)
        real(dp) ceca(xyz), pa,prea
        
        
        
        !cm scaling
        
        
        prea= 1-m(a)/18d0
        ceca(:)=cec(:,a)
        pa=p(a)
        
        
        
        ! on-diagonal
        do d = 1,xyz
            do c = 1,xyz
                if(d==c)then 
                  dQ(d,d,c) = 2*prea*pa*ceca(c) !diagonal w.r.t same coordinate
                else 
                  dQ(d,d,c) = -prea*pa*ceca(c)  !diagonal w.r.t other coordinate
                endif
            enddo
        enddo
        
        !off-diagonal:        
        dQ(1,2,1) = 1.5*prea*pa*ceca(2)
        dQ(1,2,2) = 1.5*prea*pa*ceca(1)
        dQ(1,2,3) = 0
        
        dQ(2,1,:) = dQ(1,2,:)  !symmetrize
        
        
        dQ(1,3,1) = 1.5*prea*pa*ceca(3)
        dQ(1,3,3) = 1.5*prea*pa*ceca(1)
        dQ(1,3,2) = 0
        
        dQ(3,1,:) = dQ(1,3,:)  !symmetrize
        
        
        dQ(2,3,2) = 1.5*prea*pa*ceca(3)
        dQ(2,3,3) = 1.5*prea*pa*ceca(2)
        dQ(2,3,1) = 0
        
        dQ(3,2,:) = dQ(2,3,:)  !symmetrize
        
    !enddo

END SUBROUTINE


SUBROUTINE expansion_coordinates(wg,rMEC,cec,cer2)
    integer,parameter    :: xyz=3,hho=3
    real(dp),intent(in)  :: wg(xyz,hho)
    real(dp),intent(out) :: rMEC(xyz), cec(xyz,hho),cer2(hho) 
    integer a 
    rMEC(:) = wg(:,1)+wg(:,2)+16*wg(:,3)
    rMEC(:) = rMEC(:)/18d0
    !rMEC(:) = wg(:,3)
    
    do a = 1,hho
      cec(:,a) = wg(:,a) - rMEC(:) !center of expansion coordinates
      cer2(a) = norm_square(cec(:,a)) !their square norms
    enddo
    
END SUBROUTINE






FUNCTION matrix_trace(mat) result(tr)
    real(dp) :: tr
    real(dp),intent(in) :: mat(:,:)
    integer c,sm
    sm =size(mat,1)
    if(sm/=size(mat,2))stop"not square matrix in trace checker"
    tr=0
    do c = 1,sm
      tr = tr + mat(c,c)
    enddo
END FUNCTION !///



END PROGRAM

    !cec(:,1) = wg(:,1) - rMEC(:) !h1
    !cec(:,2) = wg(:,2) - rMEC(:) !h2
    !cec(:,3) = wg(:,3) - rMEC(:) !o ![0d0,0d0,0.06532330d0]![0d0,0d0,0.106532330d0]! this worked out


!SUBROUTINE get_internal(wg,wi)
!    integer,parameter :: xyz=3,hho=3, rra=3
!    real(dp),intent(in)  :: wg(xyz,hho)
!    real(dp),intent(out) :: wi(rra)
!    real(dp) r1,r2,cosa,v1(xyz),v2(xyz)
!    v1 = wg(:,1)-wg(:,3)
!    v2 = wg(:,2)-wg(:,3)
!    r1 = norm(v1)
!    r2 = norm(v2)
!    cosa = sum(v1 * v2)/(r1*r2)
!    
!    wi = [r1,r2,cosa]
!END SUBROUTINE!/////////////


!SUBROUTINE get_rotation_matrix()
!    implicit none 
!    
!    
!END SUBROUTINE!/////////

!    ch1= [-0.5d0,0d0   ,0d0   ]
!    ch2= [0.5d0 ,0d0   ,0d0   ]
!    co = [0d0   ,0d0   ,0.85d0]
    
    !ch1= [-0.757d0, 0d0, 0d0]
    !ch2= [ 0.757d0, 0d0, 0d0]
    !co = [ 0d0  , 0d0, 0.57d0]


!off diagonal tests:
!        do c = 1,dxyz
!            
!            c1 = c1v(c)
!            c2 = c2v(c)
!            k  = kv(c)
!            
!            dQ(c1,c2,c2) = 1.5*mM1(a)*p(a)*cec(c1,a)
!            dQ(c1,c2,c1) = 1.5*mM1(a)*p(a)*cec(c2,a)
!            dQ(c1,c2,k) = 0
!            
!            
!            dQ(c2,c1,:) = dQ(c1,c2,:,a) !symmetrize
!            
!            
!        enddo

        !dQ(1,2,1,a) = 1.5*mM1(a)*p(a)*cec(2,a)
        !dQ(1,2,2,a) = 1.5*mM1(a)*p(a)*cec(1,a)
        !dQ(1,2,3,a) = 0
        !
        !dQ(2,1,:,a) = dQ(1,2,:,a)  !symmetrize
        !
        !
        !dQ(1,3,1,a) = 1.5*mM1(a)*p(a)*cec(3,a)
        !dQ(1,3,3,a) = 1.5*mM1(a)*p(a)*cec(1,a)
        !dQ(1,3,2,a) = 0
        !
        !dQ(3,1,:,a) = dQ(1,3,:,a)  !symmetrize
        !
        !
        !dQ(2,3,2,a) = 1.5*mM1(a)*p(a)*cec(3,a)
        !dQ(2,3,3,a) = 1.5*mM1(a)*p(a)*cec(2,a)
        !dQ(2,3,1,a) = 0
        !
        !dQ(3,2,:,a) = dQ(2,3,:,a)  !symmetrize
        
        
        !k=3
        !do c1 = 1,2
        !  do c2 = c1+1,3
        !    
        !    dQ(c1,c2,c2,a) = 1.5*mM1(a)*p(a)*cec(c1,a)
        !    dQ(c1,c2,c1,a) = 1.5*mM1(a)*p(a)*cec(c2,a)
        !    dQ(c1,c2,k,a) = 0
        !    
        !    
        !    dQ(c2,c1,:,a) = dQ(c1,c2,:,a) !symmetrize
        !    k=k-1
        !    
        !    
        !  enddo
        !enddo !equivalent to noloop
        
        

        
        
        
        
        !do c1 = 1,2
        !  do c2 = c1+1,3
        !    do k = 1,3
        !      
        !      if(k==c1)then
        !      k1=c1
        !      k2=c2
        !      kk=1
        !      elseif(k==c2)then
        !      k1=c2
        !      k2=c1
        !      kk=1
        !      else
        !      k1=k
        !      k2=1 !arbitrary
        !      kk=0
        !      endif
        !      
        !      dQ(c1,c2,k1,a) = kk*1.5*mM1(a)*p(a)*cec(k2,a)
        !      dQ(c2,c1,k1,a) = dQ(c1,c2,k1,a) !symmetrize
        !    enddo
        !  enddo
        !enddo

!! Derivatives of the quadrupole matrix w.r.t global coordinates of atom a in a center of mass expansion
!SUBROUTINE dQ_atomic__const_charge(cec,cer2,p,q,m,dQ,a)
!    integer,parameter :: xyz=3,hho=3, dxyz=3, xy=2
!    real(dp),intent(in)  :: cec(xyz,hho),cer2(hho), p(hho), q(xyz,xyz), m(hho)
!    real(dp),intent(out) ::dQ(xyz,xyz,dxyz)
!    !real(dp) dQ_h1(xyz,xyz,dxyz), dQ_h2(xyz,xyz,dxyz), dQ_o(xyz,xyz,dxyz)
!    integer, intent(in) :: a
!    integer c,d,od
!    real(dp) mM(hho), mM1(hho)
!    
!    integer k, c1,c2 !,cs(xyz,dxyz)
!        
!        
!        
!        
!        !cm scaling
!        mM(:) = m(:)/18d0
!        mM1(:) = 1-mM(:)
!        
!        ! on-diagonal
!        do d = 1,xyz
!            do c = 1,xyz
!                if(d==c)then 
!                  dQ(d,d,c) = 2*mM1(a)*p(a)*cec(c,a) !diagonal w.r.t same coordinate
!                else 
!                  dQ(d,d,c) = -mM1(a)*p(a)*cec(c,a)  !diagonal w.r.t other coordinate
!                endif
!            enddo
!        enddo
!        
!        !off-diagonal:        
!        dQ(1,2,1) = 1.5*mM1(a)*p(a)*cec(2,a)
!        dQ(1,2,2) = 1.5*mM1(a)*p(a)*cec(1,a)
!        dQ(1,2,3) = 0
!        
!        dQ(2,1,:) = dQ(1,2,:)  !symmetrize
!        
!        
!        dQ(1,3,1) = 1.5*mM1(a)*p(a)*cec(3,a)
!        dQ(1,3,3) = 1.5*mM1(a)*p(a)*cec(1,a)
!        dQ(1,3,2) = 0
!        
!        dQ(3,1,:) = dQ(1,3,:)  !symmetrize
!        
!        
!        dQ(2,3,2) = 1.5*mM1(a)*p(a)*cec(3,a)
!        dQ(2,3,3) = 1.5*mM1(a)*p(a)*cec(2,a)
!        dQ(2,3,1) = 0
!        
!        dQ(3,2,:) = dQ(2,3,:)  !symmetrize
!        
!    !enddo
!
!END SUBROUTINE
