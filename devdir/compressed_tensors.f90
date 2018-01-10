module compressed_tensors 

use printer_mod, only: str, printer, printo
use compressed_utils, bad=>main!, bad=>main2!,only: test_apple_g

use detrace_apple, bad=>main

use calc_derivs, only:calcDv

use compressed_arrays!, p_=>pos_, l_=>len_!, only: mm_, gg_, pos_, len_, binc_, tmm_ 
implicit none

integer, parameter :: dp = kind(0d0)

    

contains !///////////////////////////////////////////////

subroutine main !Called by 'generic_program.f90'
end subroutine




subroutine get_traces
    integer k,n, i, tvecl, nt, g1,g2
    real(dp) :: testvec(21),tvec(1000)
    n = 5
    
    !testvec = [ 1d0,3d0,5d0,6d0,7d0,9d0,8d0,5d0,4d0,2d0 ]
    !testvec = [ 1d0,3d0,5d0,6d0,7d0,9d0,8d0,5d0,4d0,2d0, 4d0, 2d0, 6d0, 7d0, 1d0 ]
    testvec = [ 1d0,3d0,5d0,6d0,7d0,9d0,8d0,5d0,4d0,2d0, 4d0, 2d0, 6d0, 7d0, 1d0, 2d0, 4d0, 3d0, 1d0, 5d0, 6d0 ]
    
    do k = 1,n/2
        tvecl = sumfac(n-2*k+1)
        nt = len_(k)
        g1 = pos_(k)+1
        g2 = pos_(k+1)
        print'(a,*(I4))', "g1,g2:", g1,g2
        print'(a,*(I4))', "gg=", gg_(g1:g2)
        
        do i = 1, tvecl
            tvec(i) = sum(testvec(tmm_(i,1:nt)*gg_(g1:g2)))
            print'(a,*(I3))',"trace index", tmm_(i,1:nt)
        enddo
        print'(a,*(f7.3))', "trace array", tvec(1:tvecl)    
    enddo
    
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
   



subroutine vector_powers(k,r,rr) 
    ! Computes successive outer products of 3-vector r(:)
    ! and stores in tricorn polytensor rr(:)
    ! r(3)   position vector
    ! k      highest rank/power
    ! rr     output tricorn-ordered polytensor
    integer, intent(in)   :: k
    real(dp), intent(in)  :: r(3)
    real(dp), intent(out) :: rr(sumfacfac(k+1)-1)
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
