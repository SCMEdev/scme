module compressed_tensors 

use printer_mod, only: str, printer, printo, printa
use compressed_utils, bad=>main!, bad=>main2!,only: test_apple_g
use compressed_arrays!, p_=>pos_, l_=>len_!, only: mm_, gg_, pos_, len_, binc_, tmm_ 
implicit none

integer, parameter :: dp = kind(0d0)

    

contains !///////////////////////////////////////////////

subroutine main !Called by 'generic_program.f90'
end subroutine


    

subroutine polarize_stone(alp,phi,dq)
    integer, parameter :: nx=2,mx=2
    real(dp), intent(in) ::  alp(pos_(nx+1),pos_(mx+1))
    !real(dp) ::  alp2(pos_(nx+1),pos_(mx+1))
    real(dp), intent(in) ::  phi(pos_(mx+1))
    real(dp), intent(out) ::  dq(pos_(nx+1))
    
    integer p1,p2,p3,p4!n1,n2,m1,m2
    
    p1=pos_(1)+1
    p2=pos_(1+1)
    p3=pos_(2)+1
    p4=pos_(2+1)
    
    !alp2=alp
    
    !alp2(p1:p2,p3:p4) = alp2(p1:p2,p3:p4)/3d0
    
    dq(p1:p2) = matmul( alp(p1:p2,p1:p4), [phi(p1:p2),phi(p3:p4)/3d0]*gg_(p1:p4) ) 
    dq(p3:p4) = matmul( alp(p3:p4,p1:p4),  phi(p1:p4)*gg_(p1:p4) ) 
    
    
            
end 


function polyinner1(narr,dfarr,nn1,nn2,mm1,mm2) result(marr)!m is rank didfarrerence
    integer :: nn1,nn2,mm1,mm2
    real(dp) :: narr(pos_(nn2+1)), dfarr(pos_(nn2+mm2+1))!narr(:), dfarr(:)
    real(dp) :: marr(pos_(mm2+1))
    integer nn, mm
    integer n1,n2,m1,m2,f1,f2
    
    marr = 0
    do mm = mm1, mm2
        m1 = pos_(mm)+1
        m2 = pos_(mm+1)
        do nn = nn1, nn2
            
            n1 = pos_(nn)+1
            n2 = pos_(nn+1)
            
            f1 = pos_(nn+mm)+1
            f2 = pos_(nn+mm+1)
            
            marr(m1:m2) = marr(m1:m2) + inner( nn+mm, nn, dfarr(f1:f2), narr(n1:n2) )
            
            !print*, mm,nn
        enddo
    enddo
    
end

function get_stone_field(narr,dfarr,nn1,nn2,mm1,mm2) result(marr)!m is rank didfarrerence
    integer,intent(in) :: nn1,nn2,mm1,mm2
    real(dp),intent(in) :: narr(:), dfarr(:)!narr(pos_(nn2+1)), dfarr(pos_(nn2+mm2+1))!
    real(dp) :: marr(pos_(mm2+1))
    real(dp) :: pref
    integer nn, mm
    integer n1,n2,m1,m2,f1,f2
    
    if(pos_(nn2+1)>size(narr))stop"get_stone_field: MOMENT array rank smaller than required multipole-order"
    if(pos_(nn2+mm2+1)>size(dfarr))stop"get_stone_field: gradient-array rank smaller than required interaction-order"
    
    !integer scaling !I am not sure why the scaling factor of 1/(2n-1)!! is needed for the potential to comply with scme. Maybe because we define the moments in some way or because 
    
    marr = 0
    do nn = nn1, nn2
        n1 = pos_(nn)+1
        n2 = pos_(nn+1)
        pref=(-1)**(nn)/dble(intff(2*nn-1,1))
        do mm = mm1, mm2
            m1 = pos_(mm)+1
            m2 = pos_(mm+1)
            
            f1 = pos_(nn+mm)+1
            f2 = pos_(nn+mm+1)
            
            marr(m1:m2) = marr(m1:m2) + pref*inner( nn+mm, nn, dfarr(f1:f2), narr(n1:n2) )
            
            !print*, mm,nn
        enddo
    enddo
    
end

!subroutine system_stone_field(na,nb,nx,ka,kb,kx,nM,rCM,qn,fk)
    !integer, intent(in) :: na,nb,nx,ka,kb,kx,nM
    !real(dp), intent(in) :: rCM(3,nM), qn(pos_(nx+1),nM)
    !real(dp), intent(out) :: fk(pos_(kx+1),nM)
subroutine system_stone_field(na,nb,ka,kb,rCM,qn,fk)
    integer, intent(in) :: na,nb,ka,kb
    real(dp), intent(in) :: rCM(:,:), qn(:,:)
    real(dp), intent(out) :: fk(:,:)
    !internal:
    integer m1,m2, nkb, k2!k1
    real(dp) :: rr(3),r2, sss(nb+kb+1)
    real(dp),dimension(pos_(nb+kb+1)) :: rrr, df
    integer nM
    
    nM=size(qn,2)
    nkb=nb+kb
    fk=0
    
    do m1 = 1,nM
        second:do m2 = 1,nM
            
            if(m1==m2)cycle second
            
            rr = rCM(:,m1)-rCM(:,m2)
            r2 = sum(rr**2)
            
            
            
            call vector_powers(nkb,rr,rrr)!(k,r,rr) 
            !call dfdu_erf(1.6_dp,r2,nkx,sss)!(a,u,nmax,ders) 
            call dfdu(r2,nkb,sss) 
            call lin_polydf(nkb,rrr,sss,df)!(nmax,rrr,sss,df)
            
            !k1=1
            !k1=pos_(0)+1
            k2=pos_(kb+1)
            
            fk(:k2,m1) = fk(:k2,m1) + get_stone_field(qn(:,m2),df,na,nb,ka,kb)
            
            
            
        enddo second
    enddo
end



function polyinner2(narr,dfarr,nn1,nn2,mm1,mm2) result(marr)!m is rank didfarrerence
    integer,  intent(in) :: nn1,nn2,mm1,mm2
    real(dp), intent(in) :: narr(pos_(nn2+1)), dfarr(pos_(nn2+mm2+1))!narr(:), dfarr(:)
    real(dp) :: marr(pos_(mm2+1))
    integer nn, mm
    integer n1,m1,f1 !remove unneeded
    integer nle,mle
    integer mi,ni
    marr = 0
    do mm = mm1, mm2 ! m-rank gradient 
        m1 = pos_(mm)
        mle = len_(mm)
        
        do nn = nn1, nn2 ! sum m-rank gradients from all n-order multipoles
            
            n1 = pos_(nn)
            nle = len_(nn)
            
            f1 = pos_(nn+mm)
            
            do mi = 1, mle ! do the inner product of the central tensor with the multipole
                do ni = 1,nle
                        !!!!!!!!!!!!! instead of plus, use interval indexing !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    marr(m1+mi) = marr(m1+mi) + dfarr(f1+mm_(mi,ni)) * narr(n1+ni) * gg_(n1+ni)
                    
                enddo
            enddo
            !print'(a,*(f10.5))','>> marr', marr
        enddo
    enddo
    
end

function polyinner_matrix(narr,dfarr,nn1,nn2,mm1,mm2) result(marr)!m is rank didfarrerence
    integer,  intent(in) :: nn1,nn2,mm1,mm2
    real(dp), intent(in) :: narr(pos_(nn2+1)), dfarr(pos_(nn2+mm2+1))!narr(:), dfarr(:)
    real(dp) :: marr(pos_(mm2+1))
    integer nn, mm
    integer n1,m1,f1
    integer nle,mle
    integer mi,ni
    real(dp) :: matrix(pos_(mm2+1),pos_(nn2+1))
    
    
    marr = 0
    do mm = mm1, mm2 ! m-rank gradient 
        m1 = pos_(mm)
        mle = len_(mm)
        
        do nn = nn1, nn2 ! sum m-rank gradients from all n-order multipoles
            
            !!alt1 !this thing is for some reason 25% slower with no and full optimization 
            !n1 = pos_(nn)+1
            !n2 = pos_(nn+1)
            !do mi = 1,mle
            !    matrix(m1+mi,n1:n2) = dfarr(f1+mm_(mi,1:nle))
            !enddo

            
            n1 = pos_(nn)
            nle = len_(nn)
            
            f1 = pos_(nn+mm)
            
            do ni = 1, nle !changing this loop to interval indexing is 25% slower with no and full optimization 
            do mi = 1, mle
                matrix(m1+mi,n1+ni) = dfarr(f1+mm_(ni,mi))
            enddo
            enddo
            
        enddo
    enddo
    
    marr = matmul( matrix, narr*gg_(1:pos_(nn2+1)) )
    
    
    
end



function inner(kq, kr, vq, vr) result(vqr) !assumes kq > kr
    integer,  intent(in) :: kq, kr
    real(dp), intent(in) :: vq(len_(kq)), vr(len_(kr))
    real(dp) vqr(len_(kq-kr))
    integer gplace, i, j, maxi, maxj
    
    maxi = len_(kq - kr)
    maxj = len_(kr)
    
    gplace = pos_(kr)
    
    vqr = 0
    do i = 1, maxi
        do j = 1,maxj
            vqr(i) = vqr(i) + vq(mm_(i,j)) * vr(j) * gg_(gplace+j)
        enddo
    enddo
   
end




subroutine lin_polydf(nmax,rrr,sss,df)
    integer n1,n2,n3,n, k1,k2,k3,k
    integer b1,b12,b123, i, in_2k, nmax
    real(dp) su,sss(0:nmax)
    real(dp) df(pos_(nmax+1)), rrr(pos_(nmax+1))
    
    
    n1=0;n2=0;n3=0
    do i = 1, pos_(nmax+1)
        if(i>1)call polynextpow(n1,n2,n3)
        n=n1+n2+n3
        su = 0 
        do k1 = 0,n1/2
            b1 = brakk(n1,k1)
            do k2 = 0,n2/2
                b12 = b1*brakk(n2,k2)
                do k3 = 0,n3/2
                    b123 = b12*brakk(n3,k3)
                    k=k1+k2+k3
                    
                    in_2k = polyfinder(n1-2*k1,n2-2*k2,n3-2*k3)
                    
                    su = su + sss(n-k) * b123 * rrr(in_2k)
                enddo
            enddo
        enddo
        df(i)=su
    enddo    
end


subroutine test_df_matrix
end    
    
    


subroutine create_df_matrix(mmax,grad,nmax,rrr,r2,a_mat,df_matrix)
    
    integer, intent(in) :: mmax, grad, nmax
    real(dp), intent(in) :: rrr(pos_(mmax+nmax+grad+1)), r2
    real(dp), intent(in) :: a_mat(0:mmax,0:nmax) !damping factors matrix
    real(dp), intent(out) :: df_matrix(pos_(mmax+grad+1),pos_(nmax+1))
    
    integer mm, mmg, mgi, mg0,  nn, n0, ni
    real(dp) aa, sss(0:mmax+nmax+grad)
    real(dp) df_temp(len_(mmax+nmax+grad)) !non-poly tensor
    
    df_matrix(1:pos_(grad+1),:) = 0 !must zeroize the first rectangle when using the grad parameter. 
    
    do mm = 0,mmax ! no grad here 
        mmg = mm+grad
        do nn = 0,nmax
            
            aa = a_mat(mm,nn)  ! no grad here because mm defines what potential to use grad only says how much to differentiate it
            call dfdu_erf(aa,r2,nn+mmg,sss) 
            
            call lin_df(nn+mmg,rrr,sss,df_temp)
            
            mg0=pos_(mmg)
            n0=pos_(nn)
            
            do mgi = 1,len_(mmg) ! reshape into large matrix
                do ni = 1, len_(nn)
                    
                    
                    df_matrix(mg0+mgi,n0+ni) = df_temp(mm_(mgi,ni))
                enddo
            enddo
            
        enddo
    enddo
    
    
    
end

subroutine lin_df(nn,rrr,sss,df)
    integer, intent(in) :: nn
    real(dp), intent(out) :: df(len_(nn))
    real(dp), intent(in) ::  rrr(pos_(nn+1)), sss(0:nn)
    
    integer n1,n2,n3, k1,k2,k3,k
    integer b1,b12,b123, i, in_2k
    real(dp) su
    
    
    
    n1=nn;n2=0;n3=0
    do i = 1, len_(nn) !create elements of df
        if(i>1)call nextpow(n1,n2,n3)
        su = 0 
        do k1 = 0,n1/2
            b1 = brakk(n1,k1)
            do k2 = 0,n2/2
                b12 = b1*brakk(n2,k2)
                do k3 = 0,n3/2
                    b123 = b12*brakk(n3,k3)
                    k=k1+k2+k3
                    
                    in_2k = polyfinder(n1-2*k1,n2-2*k2,n3-2*k3)
                    
                    su = su + sss(nn-k) * b123 * rrr(in_2k)
                enddo
            enddo
        enddo
        df(i)=su
    enddo    
end


subroutine dfdu(u,nmax,ders) 
    !Derivatives of the Coulomb potential w.r.t. r^2!
    integer,  intent(in) :: nmax
    real(dp),intent(out) :: ders(0:nmax)
    real(dp), intent(in) :: u !u=r^2
    integer n
    
    ders(0) = 1_dp/sqrt(u)
    do n = 1,nmax
        ders(n) = (-1)**n * intff(2*n-1,1) / u**n * ders(0)
    enddo
    
end

subroutine dfdu_exp(a,u,nmax,ders) 
    !Derivatives of the (1-exp(-r/a))/r potential w.r.t. r^2!
    integer,  intent(in) :: nmax
    real(dp),intent(out) :: ders(0:nmax)
    real(dp), intent(in) :: u, a !u=r^2
    real(dp) r, exp_ra
    !a = 0.7d0 !adjustable parameter
    r = sqrt(u)
    exp_ra = exp(r/a)
    
    
    
    if(nmax.ge.0)ders(0) = 2**0 * (1 - 1d0/exp_ra)/r
    if(nmax.ge.1)ders(1) = 2**1 * (a - a*exp_ra + r)/(2*a*exp_ra*u**1.5d0)
    if(nmax.ge.2)ders(2) = 2**2 * (3*a**2*(-1 + exp_ra) - 3*a*r - u)/(4*a**2*exp_ra*u**2.5d0)
    if(nmax.ge.3)ders(3) = 2**3 * (-15*a**3*(-1 + exp_ra) + 15*a**2*r + 6*a*u + u**1.5d0)/(8*a**3*exp_ra*u**3.5d0)
    if(nmax.ge.4)ders(4) = 2**4 * (105*a**4*(-1 + exp_ra) - 105*a**3*r - 45*a**2*u - 10*a*u**1.5d0 - u**2)/(16*a**4*exp_ra*u**4.5d0)
    if(nmax.ge.5)ders(5) = 2**5 * (-945*a**5*(-1 + exp_ra) + 945*a**4*r + 420*a**3*u + 105*a**2*u**1.5d0 + 15*a*u**2 + u**2.5d0)/(32*a**5*exp_ra*u**5.5d0)
    if(nmax.ge.6)ders(6) = 2**6 * (10395*a**6*(-1 + exp_ra) - 10395*a**5*r - 4725*a**4*u - 1260*a**3*u**1.5d0 - 210*a**2*u**2 - 21*a*u**2.5d0 - u**3)/(64*a**6*exp_ra*u**6.5d0)
    if(nmax.ge.7)ders(7) = 2**7 * (-135135*a**7*(-1 + exp_ra) + 135135*a**6*r + 62370*a**5*u + 17325*a**4*u**1.5d0 + 3150*a**3*u**2 + 378*a**2*u**2.5d0 + 28*a*u**3 + u**3.5d0)/(128*a**7*exp_ra*u**7.5d0)
    if(nmax.ge.8)ders(8) = 2**8 * (2027025*a**8*(-1 + exp_ra) - 2027025*a**7*r - 945945*a**6*u - 270270*a**5*u**1.5d0 - 51975*a**4*u**2 - 6930*a**3*u**2.5d0 - 630*a**2*u**3 - 36*a*u**3.5d0 - u**4)/(256*a**8*exp_ra*u**8.5d0)
    if(nmax.ge.9)ders(9) = 2**9 * (-34459425*a**9*(-1 + exp_ra) + 34459425*a**8*r + 16216200*a**7*u + 4729725*a**6*u**1.5d0 + 945945*a**5*u**2 + 135135*a**4*u**2.5d0 + 13860*a**3*u**3 + 990*a**2*u**3.5d0 + 45*a*u**4 + u**4.5d0)/(512*a**9*exp_ra*u**9.5d0)
    if(nmax.gt.9)stop"rank not implemented in dfdu_exp"
    
    
end

subroutine dfdu_erf(a,u,nmax,ders) 
    !Derivatives of the (1-exp(-r/a))/r potential w.r.t. r^2!
    integer,  intent(in) :: nmax
    real(dp),intent(out) :: ders(0:nmax)
    real(dp), intent(in) :: u, a !u=r^2
    real(dp) r, exp_ra, erf_ra, exp_ua2
    real(dp), parameter :: sqpi = sqrt(acos(-1d0)) !sqrt(pi)
    r = sqrt(u)
    exp_ra = exp(r/a)
    erf_ra = erf(r/a)
    exp_ua2 = exp(u/a**2)
    
    if(nmax.ge.0)ders(0) = 2**0 *( erf_ra/r  )
    if(nmax.ge.1)ders(1) = 2**1 *( 1/(a*exp_ua2*sqpi*u) - erf_ra/(2*u**1.5d0)  )
    if(nmax.ge.2)ders(2) = 2**2 *( -(3*a**2 + 2*u)/(2*a**3*exp_ua2*sqpi*u**2) + (3*erf_ra)/(4*u**2.5d0)  )
    if(nmax.ge.3)ders(3) = 2**3 *( (15*a**4 + 10*a**2*u + 4*u**2)/(4*a**5*exp_ua2*sqpi*u**3) - (15*erf_ra)/(8*u**3.5d0)  )
    if(nmax.ge.4)ders(4) = 2**4 *( -(105*a**6 + 70*a**4*u + 28*a**2*u**2 + 8*u**3)/(8*a**7*exp_ua2*sqpi*u**4) + (105*erf_ra)/(16*u**4.5d0)  )
    if(nmax.ge.5)ders(5) = 2**5 *( (945*a**8 + 630*a**6*u + 252*a**4*u**2 + 72*a**2*u**3 + 16*u**4)/(16*a**9*exp_ua2*sqpi*u**5) - (945*erf_ra)/(32*u**5.5d0)  )
    if(nmax.ge.6)ders(6) = 2**6 *( -(10395*a**10 + 6930*a**8*u + 2772*a**6*u**2 + 792*a**4*u**3 + 176*a**2*u**4 + 32*u**5)/(32*a**11*exp_ua2*sqpi*u**6) + (10395*erf_ra)/(64*u**6.5d0)  )
    if(nmax.ge.7)ders(7) = 2**7 *( (135135*a**12 + 90090*a**10*u + 36036*a**8*u**2 + 10296*a**6*u**3 + 2288*a**4*u**4 + 416*a**2*u**5 + 64*u**6)/(64*a**13*exp_ua2*sqpi*u**7) - (135135*erf_ra)/(128*u**7.5d0)  )
    if(nmax.ge.8)ders(8) = 2**8 *( ((-2*r*(2027025*a**14 + 1351350*a**12*u + 540540*a**10*u**2 + 154440*a**8*u**3 + 34320*a**6*u**4 + 6240*a**4*u**5 + 960*a**2*u**6 + 128*u**7))/(a**15*exp_ua2*sqpi) + 2027025*erf_ra)/(256*u**8.5d0)  )
    if(nmax.ge.9)ders(9) = 2**9 *( ((2*r*(34459425*a**16 + 22972950*a**14*u + 9189180*a**12*u**2 + 2625480*a**10*u**3 + 583440*a**8*u**4 + 106080*a**6*u**5 + 16320*a**4*u**6 + 2176*a**2*u**7 + 256*u**8))/(a**17*exp_ua2*sqpi)- 34459425*erf_ra)/(512*u**9.5d0)  )
    if(nmax.gt.9)stop"rank not implemented in dfdu_erf"
    
    
    
end






subroutine apple1_df(nmax,rrr,rpow,df)
    integer , intent(in) :: nmax
    real(dp),intent(out) :: df(pos_(nmax+1))
    real(dp), intent(in) :: rrr(pos_(nmax+1))
    real(dp), intent(in) :: rpow(2*nmax+1)
    integer n, n1, n2
    
    do n = 0,nmax
        n1 = pos_(n)+1
        n2 = pos_(n+1)
        df(n1:n2) = (-1)**n * rpow(2*n+1) * intff(2*n-1,1) * detracer(rrr(n1:n2),n)
        !print'(a,*(f10.5))','>> marr', df
    enddo
    
end
        
    
    

function apple_potgrad(qq,nn,kk,rpows,dtrrr) result(potgrad)
    integer, intent(in)  :: nn, kk
    real(dp), intent(in) :: qq(len_(nn)), dtrrr(len_(nn+kk)), rpows(:)
    real(dp)             ::  potgrad(len_(kk))
    
    potgrad = (-1)**(kk+1) * intff(2*(nn+kk)-1,2*nn-1) * rpows(nn+kk+1) * inner(nn+kk,nn,dtrrr,qq) !Kom ihåg att detracern är the shiiiiiit ohc inte har (2*n-1)!! i sig
    
end    





function polydet(invec,nmax) result(outvec)
    ! Possibly efficient implementation of Applequist's detracer for polytensors
    ! It requires the trace-matrix tmm. 
    ! The input polytensor AA is used as the temporary array holding the traces
    ! and thus one kan skip the l=0 loop that only copies the original tensor into the temporary array. 
    ! The detraced tensor is written into NEW and in the end AA=new. 
    ! One could probably use pointers to avoid the last copy, or is the compiler fixing that?
    !
    !Equation reference to (Jon Applequist 1989, Traceless cartesian...)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    integer, intent(in) :: nmax
    real(dp), intent(in)  :: invec(pos_(nmax+1))
    real(dp) :: outvec(pos_(nmax+1)), tempvec(pos_(nmax+1)), su
    
    integer g1,g2, cols, l_m
    integer l1,l2,l3,l32, p_m, ti
    
    integer l, i
    integer m, mm(3), br3, br32, br, n1,n2,n3, n2_1, n, p_n
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !Copy input vector, to not destroy it and to use intent(in)
    tempvec = invec
    
    
    outvec(1:4) = tempvec(1:4) !the 0 & 1 ranks
    
    do n = 2,nmax ! loop over ranks
        n2_1 = 2*n-1
        p_n = pos_(n) !polytensor rank position
        
        ! Eq. 2.4:
        do l = 1,n/2 !loop over traces. Starts on 1 in this implementation!
            m = n-2*l       !rank of trace (the non-deltas)
            l_m = len_(m)    !lenght of trace vector
            p_m = pos_(m)    !position in trace-polytensor
            
            cols = len_(l)  !first cols of tm-row that sample A
            g1 = pos_(l)+1  !
            g2 = pos_(l+1)  ! section of apple-g polytensor
            
            do i = 1, l_m !over rows of tm
                tempvec(p_m+i) = dot_product(tempvec(p_n + tmm_(i,1:cols)), gg_(g1:g2)) ! 
            enddo
        enddo
        
        ! Eq. 5.6:
        n1=n;n2=0;n3=0
        do i = 1, len_(n) ! single-loop over the tensor entries n1,n2,n3, or i
            if(i>1)call nextpow(n1,n2,n3)
            
            su = 0
            do l3 = 0, n3/2 !the tripple sum in Eq 5.6 
                br3 = brakk(n3,l3) 
                mm(3) = n3-2*l3
                
                do l2 = 0, n2/2
                    br32 = br3 * brakk(n2,l2)
                    l32 = l3+l2
                    mm(2) = n2-2*l2
                    
                    do l1 = 0, n1/2
                        br = br32 * brakk(n1,l1)
                        l = l32+l1
                        mm(1) = n1-2*l1
                        
                        ti = polyfind(mm)!finder(nn-2*ll)
                        su = su + (-1)**l / dble(intff(n2_1,n2_1-2*l)) * br *  tempvec(ti)
                    enddo
                enddo
            enddo
            
            outvec(pos_(n)+i) = su
            
        enddo
        
    enddo
end


function detracer(AA,n) result (newvec)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    integer, intent(in) :: n
    real(dp), intent(in)  :: AA(len_(n))
    real(dp) :: newvec(len_(n)), tvec(1000)
    real(dp) su
    
    integer g1,g2, cols, rows
    integer l1,l2,l3,l32, t1, ti
    
    integer l, i
    integer m, mm(3), br3, br32, br, n1,n2,n3, n21
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    n21 = 2*n-1
    tvec=0
    do l = 0,n/2        !loop over traces
        m = n-2*l       !rank of trace (the non-deltas)
        rows = len_(m)  !lenght of trace vector
        t1 = pos_(m)    !position in trace-polytensor
        
        cols = len_(l)  !first cols of tm-row that sample A
        
        g1 = pos_(l)+1  !
        g2 = pos_(l+1)  ! section of apple-g polytensor
        
        do i = 1, rows !over rows of tm
            tvec(t1+i) = dot_product(AA(tmm_(i,1:cols)), gg_(g1:g2)) !tm-row * apple-g
        enddo
    enddo

    !print'(a,*(f7.3))', "full trace array",tvec(1:pos_(n+1)) !float måste ha en rutin !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    n1 = n; n2 = 0; n3 = 0
    do i = 1, len_(n)
        if(i>1)call nextpow(n1,n2,n3)
        
        su = 0
        do l3 = 0, n3/2
            br3 = brakk(n3,l3)
            mm(3) = n3-2*l3
            
            do l2 = 0, n2/2
                br32 = br3 * brakk(n2,l2)
                l32 = l3+l2
                mm(2) = n2-2*l2
                
                do l1 = 0, n1/2
                    br = br32 * brakk(n1,l1)
                    l = l32+l1
                    mm(1) = n1-2*l1
                    
                    
                    ti = polyfind(mm)!finder(nn-2*ll)
                    su = su + (-1)**l / dble(intff(n21,n21-2*l)) * br *  tvec(ti)
                enddo
            enddo
        enddo
        
        newvec(i) = su
        
    enddo
    !newvec = newvec/dble(intff(2*n-1,1))
    
    

end



subroutine print_trace_keys
    integer ll(3), nn(3),nn_2(3), n_2, l, it, numbers(100), key(3), ind, nrank, nlen
    !integer i
    logical proceed
    
    ! To investigate all the trace subtractions associated with one specific tensor entry
    nrank = 10
    
    ! initial nkey & nkey/2
    nn = [nrank,0,0]
    
    nlen = sumfac(nrank+1)
    do ind = 1, nlen ! Go through all tensor entries
        
        if (ind>1) call nextpown(nn)
        nn_2 = nn/2
        n_2 = sum(nn_2)
        
        print*, ">>> Entry=["//str(nn,2)//"] index="//str(ind)//' finder(nn) = '//str(finder(nn))
        print*, "     Half=["//str(nn_2,2)//"]"
        
        do l = 1, n_2 ! Increment the l-rank (the number of deltas)
            call init_pow_nl(nn_2,l,ll)
            
            print*, "    l="//str(l)//":"
            
            
            proceed = .true.
            it=0
            do while (proceed) ! Go through all compatible ll-keys
                it=it+1
                key = nn - 2*ll
                !print*, '       ll=['//str(ll,2)//']  finder(ll) = '//str(finder(ll))
                print*, '      key=['//str(key,2)//']  finder(key) = '//str(finder(key))
                call next_pow_nl(nn_2,ll,proceed)
            enddo    
            
            !print*, "    # entries = "//str(it)
            numbers(l)=it
        enddo
        print*, "<<< #entries = "//str(numbers(1:l-1),3)
        print*, ""
    enddo
    !print*, "choices "//str( choose(n_2-l,l)/fac(ll(1))/fac(ll(2))/fac(nn_2(3)) )
    print*, "ABOVE: test_nl ----------------------------------------------------------------------"
end subroutine





    
!function symouter(k1,k2,v1,v2) result(vout)
!    integer, intent(in) :: k1,k2
!    real(dp), intent(in) :: v1(:), v2(:)
!    real(dp) vout( (k1+k2+1)*(k1+k2+2)/2 )
!    integer i, j
!    
!    vout=0
!    do i = 1, sumfac(k1+1)
!    do j = 1, sumfac(k2+1)
!        
!        vout(matr(i,j)) = vout(matr(i,j)) + v1(i)*v2(j)*hhh(i,j,k1,k2)
!        enddo
!        enddo
!    
!end



function symouter(k1,k2,v1,v2) result(vout)
    ! Symmetrizing outer product, correct scaling is applied with the built in h-function = (gi*gj*cho)/gij 
    ! It requires the index-matrix "matr" and the Applequist g-vectors in "gg"
    integer, intent(in) :: k1,k2
    real(dp), intent(in) :: v1(:), v2(:)
    real(dp) vout( (k1+k2+1)*(k1+k2+2)/2 )
    integer i, j
    integer gi,gj,gij, k12, pi, pj, pij, mij, cho
    
    
    k12 = k1+k2
    cho = binc_(k12,k1)
    
    
    pi  = pos_(k1)
    pj  = pos_(k2)
    pij = pos_(k12)
    
    vout=0
    
    do i = 1, len_(k1)
        gi  = gg_(pi + i)
        do j = 1, len_(k2)
            
            mij = mm_(j,i)!faster order, symmetric anyway
            
            gj  = gg_(pj + j)
            gij = gg_(pij + mij)
            
            vout(mij) = vout(mij) + v1(i)*v2(j) * (gi*gj*cho)/gij
            
        enddo
    enddo
    
end



function potgrad(qq,nn,kk,rpows,rrr) 
    real(dp), intent(in) :: qq((nn+1)*(nn+2)/2),rpows(:),rrr(:)
    integer, intent(in) :: nn, kk
    real(dp) potgrad((kk+1)*(kk+2)/2)
    integer i!,j
    integer k1, k2, n1, n2, kni2, sig, k_i, n_i
    
    potgrad = 0
    do i = 0, min(nn,kk) !0!
        
        kni2 = (kk+nn-i)*2
        sig = kk+i+1
        k_i = kk-i
        n_i = nn-i
        
        k1 = pos_(k_i)+1
        k2 = pos_(k_i+1)
        n1 = pos_(n_i)+1
        n2 = pos_(n_i+1)
        
        potgrad(:) = potgrad(:) &
                    + (-1)**sig * intfac(nn,n_i) * intff(kni2-1,2*nn-1) * rpows(kni2+1) & !why +nn in first term???
                    * symouter(k_i, i, rrr(k1:k2),   inner(nn,n_i, qq, rrr(n1:n2))   )
    enddo
    
end
  !integer, parameter :: rpos(8) = [1,      4,     10,     20,     35,     56,     84,    120] ! remove this, use pos_
    
  !print*, 'rrr(0) i potgrad',rrr(0)
    
  !real(dp) temp(1)
  !print'((a,I3,a,2I3,a,2i3,a,2i3))', ' i:',i,',   kni:',kni,kni2, ',   ki:',ki1,ki2,',   ni:',ni1,ni2
  
  !temp = inner(nq,nq-i, qq, rrr(ni1:ni2))
  !print*, 'pref',(-1)**kni * intfac(nq,nq-i) * intff(kni2-1,2*nq-1)
  !print*, 'rpow',kni2+1
  !print*, 'inner', inner(nq,nq-i, qq, rrr(ni1:ni2)), sum(qq*rrr(ni1:ni2))
  !print*, 'symouter 1 ',symouter(kk-i, i, rrr(ki1:ki2),  temp) 
  !print*, 'symouter 2 ',temp(1)*rrr(ki1:ki2)
  !print*, 'temp', temp
  !print*, 'rrr 1 ',rrr(ni1:ni2)
  !print*, 'rrr 2 ',rrr(ki1:ki2)
  !print*, 'rrr 3 ',symouter (1,1,rrr(ni1:ni2),rrr(ni1:ni2))
  
  !print'(*(f10.3))', potgrad
   



pure subroutine vector_powers(k,r,rr) 
    ! Computes successive outer products of 3-vector r(:)
    ! and stores in tricorn polytensor rr(:)
    ! r(3)   position vector
    ! k      highest rank/power
    ! rr     output tricorn-ordered polytensor
    integer, intent(in)   :: k
    real(dp), intent(in)  :: r(3)
    real(dp), intent(out) :: rr(pos_(k+1))
    integer pl,cl,px,cx,i,pz,cz, py, cy
    
    rr(1) = 1
    rr(2:4) = r(:)
    
    do i = 2,k
       ! p=previous, c=current, l=length
       ! x, y, z refer to the position of the x-only, y-only, z-only rows. 
       px = pos_(i-1)+1 !simfacfac(i-1) !+1 because interval indexing, no index addition
       pl = len_(i-1) !sumfac(i)
       cl = len_(i) ! sumfac(i+1)
       
       cx = px+pl
       pz = cx-1
       cz = pz+cl
       cy = cx+pl
       py = cx-i
       
       rr(cx:cy-1) = rr(px:pz) *r(1)
       rr(cy:cz-1) = rr(py:pz) *r(2)
       rr(cz)      = rr(pz)    *r(3)
    enddo
end subroutine





end module
