module compressed_tensors 

use printer_mod, only: str, printer, printo
use compressed_utils, bad=>main!, bad=>main2!,only: test_apple_g
use compressed_arrays!, p_=>pos_, l_=>len_!, only: mm_, gg_, pos_, len_, binc_, tmm_ 
implicit none

integer, parameter :: dp = kind(0d0)

    

contains !///////////////////////////////////////////////

subroutine main !Called by 'generic_program.f90'
end subroutine




subroutine dfdu(u,ders,nmax) 
    integer,  intent(in) :: nmax
    real(dp),intent(out) :: ders(0:nmax)
    real(dp), intent(in) :: u !u=r^2
    integer n
    
    ders(0) = 1_dp/sqrt(u)
    do n = 1,nmax
        ders(n) = (-1)**n * intff(2*n-1,1)/2**n / u**n * ders(0)
    enddo
    
end


subroutine dejun_df(nmax,rrr,df)
    integer n1,n2,n3,n, k1,k2,k3,k
    integer b1,b12,b123, i, in_2k, nmax
    real(dp) su,r2, fff(0:nmax)
    real(dp) df(pos_(nmax+1)), rrr(pos_(nmax+1))

    
    r2=sum(rrr(2:4)**2)
    call dfdu(r2,fff,nmax)
    
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
                    
                    su = su + 2**(n-k) * fff(n-k) * b123 * rrr(in_2k)
                enddo
            enddo
        enddo
        df(i)=su
    enddo    
end


function apple_potgrad(qq,nn,kk,rpows,dtrrr) result(potgrad)
    integer, intent(in)  :: nn, kk
    real(dp), intent(in) :: qq(len_(nn)), dtrrr(len_(nn+kk)), rpows(:)
    real(dp)             ::  potgrad(len_(kk))
    integer i!,j
    !integer k1, k2, n1, n2, kni2, sig, k_i, n_i
    
    
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




function inner(kq, kr, vq, vr) result(vqr) !assumes kq > kr
   integer, intent(in) :: kq, kr
   real(dp), intent(in) :: vq(:), vr(:)
   real(dp) vqr((kq-kr+1)*(kq-kr+2)/2)
   integer kqr, gplace, i, j, maxi, maxj
   
   !if(kr
   
   kqr = kq - kr
   maxi = (kqr+1)*(kqr+2)/2
   maxj = (kr+1)*(kr+2)/2
   
   gplace = pos_(kr)
   
   vqr(:) = 0
   do i = 1, maxi
     do j = 1,maxj
       vqr(i) = vqr(i) + vq(mm_(i,j)) * vr(j) * gg_(gplace+j)
       enddo
       enddo
   
end

    
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
