module compressed_tensors

use printer_mod, only: str, printer, printo
use compressed_utils, bad=>main!,only: test_apple_g

use detrace_apple, bad=>main

use calc_derivs, only:calcDv

use compressed_arrays!, only: matr, gg, pos00, len00, matchoo, tmatr 
implicit none

integer, parameter :: dp = kind(0d0)

!integer i,j
!integer, parameter :: mat(0:k,0:k) = reshape( [((1, i = 0,k),j=0,k)],shape(mat))
    

contains !///////////////////////////////////////////////





subroutine main
    !call testing
    !call test_sumfac
    !call test_rpow
    !call test_nextpow(7)
    !call compow2(7)
    !call subdiv_pow(5,5)
    
    
    call test_old_field
    call test_potgrad
    !call test_inner
    
    call test_mp_pot
    
    call printo(tmatr,3,1)
    call printo(matr,3,1)
    !call test_intfac_ff
    
    !print*, pos00+1
    !call test_hhh(3,3)
    !print*, ""
    !call subdiv_pow(4,3)
    !print*, ""
    !call subdiv_pow(3,4)
    !print*,"fac 5 / 2 3", fac(5)/( fac(3)*fac(2))
    !call test_factorial(5)
    
    !call test_choose
    !call test_matri(7)
    
    !print*, ""
    !call print_product_index_matrix (6)
    
    !call test_matr(7)
    !call printer(matr, 'matr',1)
    !call test_apple_g
    
    !call test_rrpow
    !print*,""
    !call test_rrr
    
    !call test_next_rev_key2n(4)
    !integer i, j
    !do i = 1,10
    !  do j = 1,10
    !    write(*,'(I6,a)', advance="no") choose(i,j),", "
    !    enddo
    !  print*, ""
    !  enddo
    
    !call h_testing
    
    !call pascal_matrix(10)
    !do i00 = 0,8
    !  print*, sumfac(i00+1),len00(i00)
    !enddo
    
    
    
end subroutine

! Testing /////////////////////////////////////////////////////




subroutine h_testing
integer i, k , n(3)
k=5
n = 0
n(1) = k
do i = 1, sumfac(k+1)
  if(i>1) n = nextpov(n)
  print'(3I3,a,3I4,a,I4)',n,",  ", fac(n(1)), fac(n(2)), fac(n(3)),",  ", fac(n(1))*fac(n(2))*fac(n(3))
  enddo


end



subroutine subdiv_pow(ii,kii)
integer ii, kii
integer i, j, it
integer a1,b1,c1, a2,b2,c2, aa,bb,cc
integer sfii, sfkii
a1 = ii
aa = ii + kii

b1 = 0
c1 = 0
bb = 0
cc = 0

sfii = sumfac(ii+1)
sfkii = sumfac(kii+1)

it = 0
do i = 1, sfii
   if (i>1)call nextpow(a1,b1,c1)
   b2 = 0
   c2 = 0
   a2 = kii
   do j = 1, sfkii
      if (j>1)call nextpow(a2,b2,c2)
      it = it+1
      aa = a1+a2
      bb = b1+b2
      cc = c1+c2
      !print '(3(I3,2I2),7I4)', a1,b1,c1, a2,b2,c2, aa,bb,cc, &
      !      finder([aa,bb,cc]),&
      !      matr(i,j), &
      !      i + j + (vv(i)-1)*(vv(j)-1) - 1, &
      !      i + j + vv1(i)*vv1(j) - 1, &
      !      hh(a1,b1,c1, a2,b2,c2), &
      !      i, j !, it
      
      !write(*,'(I3)', advance="no") hh(a1,b1,c1, a2,b2,c2)
      !write(*,'(*(I3))') hh(a1,b1,c1, a2,b2,c2), &
      !                 ( apple_g([a1,b1,c1])*apple_g([a2,b2,c2])*choose(ii+kii,ii) )/apple_g([aa,bb,cc]), &
      !                 ( gg(pos00(ii)+i)*gg(pos00(kii)+j)*choose(ii+kii,ii) ) / gg(pos00(ii+kii) + matr(i,j)), & 
      !                 00,&
      !                 apple_g([a1,b1,c1]), apple_g([a2,b2,c2]), choose(ii+kii,ii), apple_g([aa,bb,cc]), &
      !                 00, &
      !                 gg(pos00(ii)+i), gg(pos00(kii)+j), choose(ii+kii,ii), gg(pos00(ii+kii) + matr(i,j)), & 
      !                 00
      
      !write(*,'(I6)', advance="no") ( apple_g([a1,b1,c1])*apple_g([a2,b2,c2])*choose(ii+kii,ii) )/apple_g([aa,bb,cc])
      !write(*,'(I6)', advance="no") ( gg(pos00(ii)+i)*gg(pos00(kii)+j)*choose(ii+kii,ii) ) / gg(pos00(ii+kii) + matr(i,j))
      
      
      !write(*,'(I6)', advance="no") ( gg(pos00(ii)+i)*gg(pos00(kii)+j)*choose(ii+kii,ii) ) / gg(pos00(ii+kii) + matr(i,j)) &
      !                              - ( apple_g([a1,b1,c1])*apple_g([a2,b2,c2])*choose(ii+kii,ii) )/apple_g([aa,bb,cc])
      
      write(*,'(I6)', advance="no") hh(a1,b1,c1, a2,b2,c2) - hhh(i,j,ii,kii)
      enddo
   print*,""
   enddo
    
   !do here!
!print'(*(I3))', vv(1:sumfac(ii+kii+1))
!print'(*(I3))', vv1(1:sumfac(ii+kii+1))
!print*, sumfac(ii+kii+1)
end


function hh(a1,b1,c1, a2,b2,c2) !Make a matrix of this sheeeet
    integer a1,b1,c1, a2,b2,c2
    integer aa,bb,cc, hh
    aa = a1+a2
    bb = b1+b2
    cc = c1+c2
    
    hh = ( fac(aa)*fac(bb)*fac(cc) ) / ( fac(a1)*fac(a2)*fac(b1)*fac(b2)*fac(c1)*fac(c2) )

end

subroutine test_hhh(k1, k2)
    integer i1, i2, k1,k2
    integer gp1, gp2, gp12, h, ch
    gp1 = pos00(k1)
    gp2 = pos00(k2)
    gp12 = pos00(k1+k2)
    ch = choose(k1+k2,k1)
    do i1 = 1, sumfac(k1+1)
      do i2 = 1, sumfac(k2+1)
        h = ( gg(gp1+i1)*gg(gp2+i2)*ch ) / gg(gp12 + matr(i1,i2))
        write(*,'((I4))', advance="no") h-hhh(i1,i2,k1,k2)
        enddo
      print*,""
      enddo
    print*, "above, TEST_HHH ------------------------------------------------------------------"
end
function hhh(i,j,ki,kj) !Make a matrix of this sheeeet
    integer, intent(in) :: i,j, ki,kj
    integer kij, hhh, cc, gi,gj,gij
    kij = ki+kj
    cc = matchoo(kij,ki)
    gi = gg(pos00(ki)+i)
    gj = gg(pos00(kj)+j) 
    gij = gg(pos00(kij) + matr(i,j))
    hhh = (gi*gj*cc) / gij
    
end

subroutine test_inner_symouter
    real(dp) :: v1(3), v2(6), v3(10), v4(15), v5(21)
    real(dp) :: f1(3), f2(3,3), f3(3,3,3), f4(3,3,3,3), f5(3,3,3,3,3), of1(3), of2(3,3), of3(3,3,3), of4(3,3,3,3)
    integer i,j,k, l, m
    
    !v1 = [ 1d0, 2d0, 3d0]
    !v2 = [ 1d0, 2d0, 3d0, 4d0, 5d0, 6d0]
    !v3 = [ 1d0, 2d0, 3d0, 4d0, 5d0, 6d0, 7d0, 8d0, 9d0, 10d0]
    !v4 = [ 1d0, 2d0, 3d0, 4d0, 5d0, 6d0, 7d0, 8d0, 9d0, 10d0, 11d0, 12d0, 13d0, 14d0, 15d0]
    !v5 = [ 1d0, 2d0, 3d0, 4d0, 5d0, 6d0, 7d0, 8d0, 9d0, 10d0, 11d0, 12d0, 13d0, 14d0, 15d0, 16d0 ,17d0 ,18d0 ,19d0 ,20d0 ,21d0 ]
    call random_number(v1)
    call random_number(v2)
    call random_number(v3)
    call random_number(v4)
    call random_number(v5)
    
    
    f1 = v1
    f2 = reshape(expand(v2,2),shape(f2))
    f3 = reshape(expand(v3,3),shape(f3))
    f4 = reshape(expand(v4,4),shape(f4))
    f5 = reshape(expand(v5,5),shape(f5))
    
    !call printer(f2,'f2',2)
    !call printer(f3,'f3',2)
    
    of1=0
    of2=0
    of3=0
    do i = 1,3
      do j = 1,3
        do k = 1,3
          
          of1(i) = of1(i) + f2(k,j)*f3(k,j,i)
          
          do l = 1, 3
            of2(j,i) = of2(j,i) + f2(l,k)*f4(l,k,j,i)
            do m = 1, 3
              of3(k,j,i) = of3(k,j,i) + f2(m,l)*f5(m,l,k,j,i)
            enddo
          enddo
        enddo
      enddo
    enddo
    
    
    
    !of1=0
    of2=0
    !of3=0
    of4=0
    do i = 1,3
      do j = 1,3
        do k = 1,3
          
          !of3(i,j,k) = of3(i,j,k) + f1(i)*f2(k,j) + f1(j)*f2(k,i) + f1(k)*f2(i,j)
          
          do l = 1, 3
            !of4(l,k,j,i) = of4(l,k,j,i) + f2(l,k)*f2(j,i)*6 !+ f2(l,j)*f2(k,i) + f2(l,i)*f2(k,j) 
            of4(l,k,j,i) = of4(l,k,j,i) + of1(l)*of3(k,j,i) + of1(k)*of3(l,j,i) + of1(j)*of3(l,k,i) + of1(i)*of3(k,l,j) 
          !  do m = 1, 3
          !    of3(k,j,i) = of3(k,j,i) + f2(m,l)*f5(m,l,k,j,i)
          !  enddo
          enddo
        enddo
      enddo
    enddo
    
    !call printer(of4,
    
    print*,""
    print'(*(f10.4))', compress(reshape(of4,[3**4]),4)
    print'(*(f10.4))', symouter(1,3,inner(3,2,v3,v2),inner(5,2,v5,v2))
    !print*,""
    !print'(*(f10.4))', compress(reshape(of4,[3**4]),4)
    !print'(*(f10.4))', symouter(1,3,v1,v3)
    
    print*, "above, INNER+OUTER ------------------------------------------------------------------"

    
end

subroutine test_inner
    real(dp) :: v1(3), v2(6), v3(10), v4(15), v5(21)
    real(dp) :: f1(3), f2(3,3), f3(3,3,3), f4(3,3,3,3), f5(3,3,3,3,3), of1(3), of2(3,3), of3(3,3,3)
    integer i,j,k, l, m
    
    v1 = [ 1d0, 2d0, 3d0]
    v2 = [ 1d0, 2d0, 3d0, 4d0, 5d0, 6d0]
    v3 = [ 1d0, 2d0, 3d0, 4d0, 5d0, 6d0, 7d0, 8d0, 9d0, 10d0]
    v4 = [ 1d0, 2d0, 3d0, 4d0, 5d0, 6d0, 7d0, 8d0, 9d0, 10d0, 11d0, 12d0, 13d0, 14d0, 15d0]
    v5 = [ 1d0, 2d0, 3d0, 4d0, 5d0, 6d0, 7d0, 8d0, 9d0, 10d0, 11d0, 12d0, 13d0, 14d0, 15d0, 16d0 ,17d0 ,18d0 ,19d0 ,20d0 ,21d0 ]
    
    
    f1 = v1
    f2 = reshape(expand(v2,2),shape(f2))
    f3 = reshape(expand(v3,3),shape(f3))
    f4 = reshape(expand(v4,4),shape(f4))
    f5 = reshape(expand(v5,5),shape(f5))
    
    !call printer(f2,'f2',2)
    !call printer(f3,'f3',2)
    
    of1=0
    of2=0
    of3=0
    do i = 1,3
      do j = 1,3
        do k = 1,3
          
          of1(i) = of1(i) + f2(k,j)*f3(k,j,i)
          
          do l = 1, 3
            of2(j,i) = of2(j,i) + f2(l,k)*f4(l,k,j,i)
            do m = 1, 3
              of3(k,j,i) = of3(k,j,i) + f2(m,l)*f5(m,l,k,j,i)
            enddo
          enddo
        enddo
      enddo
    enddo
    
    
    
    print*,"loop-comparison 3,2-rank"
    print'(*(f10.4))', of1
    print'(*(f10.4))', inner(3,2,v3,v2)
    
    print*,"loop-comparison 4,2-rank"
    print'(*(f10.4))', compress(reshape(of2,[3**2]),2)
    print'(*(f10.4))', inner(4,2,v4,v2)
    
    print*,"loop-comparison 5,2-rank"
    print'(*(f10.4))', compress(reshape(of3,[3**3]),3)
    print'(*(f10.4))', inner(5,2,v5,v2)
    
    
    print*,"scalar=3 5,0-rank"
    print'(*(f10.4))', v5
    print'(*(f10.4))', inner(5,0,v5,[1d0])
    
    print*,"scalar=3 2,0-rank"
    print'(*(f10.4))', v2
    print'(*(f10.4))', inner(2,0,v2,[1d0])
    
    print*," 2,2-rank"
    print'(*(f10.4))', v2
    print'(*(f10.4))', inner(2,2,v2,v2)
    
    print*," 1,1-rank"
    print'(*(f10.4))', v1
    print'(*(f10.4))', inner(1,1,v1,v1)
    
    print*, "above, INNER ------------------------------------------------------------------"

end


function inner(kq, kr, vq, vr) result(vqr) !assumes kq > kr
   integer, intent(in) :: kq, kr
   real(dp), intent(in) :: vq(:), vr(:)
   real(dp) vqr((kq-kr+1)*(kq-kr+2)/2)
   integer kqr, gplace, i, j, maxi, maxj
   
   !if(kr
   
   kqr = kq - kr
   maxi = (kqr+1)*(kqr+2)/2
   maxj = (kr+1)*(kr+2)/2
   
   gplace = pos00(kr)
   
   vqr(:) = 0
   do i = 1, maxi
     do j = 1,maxj
       vqr(i) = vqr(i) + vq(matr(i,j)) * vr(j) * gg(gplace+j)
       enddo
       enddo
   
end

subroutine test_symouter
    real(dp) :: v1(3), v2(6), v3(10), v4(15), v5(21)
    real(dp) :: f1(3), f2(3,3), f3(3,3,3), f4(3,3,3,3), f5(3,3,3,3,3), of1(3), of2(3,3), of3(3,3,3), of4(3,3,3,3)
    integer i,j,k, l
    
    v1 = [ 1d0, 2d0, 3d0]
    v2 = [ 1d0, 2d0, 3d0, 4d0, 5d0, 6d0]
    v3 = [ 1d0, 2d0, 3d0, 4d0, 5d0, 6d0, 7d0, 8d0, 9d0, 10d0]
    v4 = [ 1d0, 2d0, 3d0, 4d0, 5d0, 6d0, 7d0, 8d0, 9d0, 10d0, 11d0, 12d0, 13d0, 14d0, 15d0]
    v5 = [ 1d0, 2d0, 3d0, 4d0, 5d0, 6d0, 7d0, 8d0, 9d0, 10d0, 11d0, 12d0, 13d0, 14d0, 15d0, 16d0 ,17d0 ,18d0 ,19d0 ,20d0 ,21d0 ]
    
    
    f1 = v1
    f2 = reshape(expand(v2,2),shape(f2))
    f3 = reshape(expand(v3,3),shape(f3))
    f4 = reshape(expand(v4,4),shape(f4))
    f5 = reshape(expand(v5,5),shape(f5))
    
    !call printer(f2,'f2',2)
    !call printer(f3,'f3',2)
    
    of1=0
    of2=0
    of3=0
    of4=0
    
    do i = 1,3
      do j = 1,3
        do k = 1,3
          
          of3(i,j,k) = of3(i,j,k) + f1(i)*f2(k,j) + f1(j)*f2(k,i) + f1(k)*f2(i,j)
          
          do l = 1, 3
            !of4(l,k,j,i) = of4(l,k,j,i) + f2(l,k)*f2(j,i)*6 !+ f2(l,j)*f2(k,i) + f2(l,i)*f2(k,j) 
            of4(l,k,j,i) = of4(l,k,j,i) + f1(l)*f3(k,j,i) + f1(k)*f3(l,j,i) + f1(j)*f3(l,k,i) + f1(i)*f3(k,l,j) 
          !  do m = 1, 3
          !    of3(k,j,i) = of3(k,j,i) + f2(m,l)*f5(m,l,k,j,i)
          !  enddo
          enddo
        enddo
      enddo
    enddo
    
    !call printer(of4,
    
    
    print*,"3-rank comparison with loops:"
    print'(*(f10.4))', compress(reshape(of3,[3**3]),3)
    print'(*(f10.4))', symouter(1,2,v1,v2)
    print*,"4-rank comparison with loops:"
    print'(*(f10.4))', compress(reshape(of4,[3**4]),4)
    print'(*(f10.4))', symouter(1,3,v1,v3)
    
    print*, "different order:"
    print'(*(f10.4))', symouter(1,3,v1,v3)
    print'(*(f10.4))', symouter(3,1,v3,v1)
    
    print*, "scalar=3 in pos 2:"
    print'(*(f10.4))', v3
    print'(*(f10.4))', symouter(3,0,v3,[3d0])
    
    print*, "scalar=3 in pos 2:"
    print'(*(f10.4))', v1
    print'(*(f10.4))', symouter(1,0,v1,[3d0])
    
    print*, "scalar=3 in pos1:"
    print'(*(f10.4))', v3
    print'(*(f10.4))', symouter(0,3,[3d0],v3)
    
    print*, "scalar=3 in pos1:"
    print'(*(f10.4))', v1
    print'(*(f10.4))', symouter(0,1,[3d0],v1)
    
    print*, "above, OUTER ------------------------------------------------------------------"
    
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
    
    logical pri
    
    k12 = k1+k2
    cho = matchoo(k12,k1)
    
    pri=.false.!.true.
    if(pri)print*, "choose", cho
    
    pi  = pos00(k1)
    pj  = pos00(k2)
    pij = pos00(k12)
    if(pri)print*, "pi, pj, pij", pi, pj, pij
    
    vout=0
    
    do i = 1, len00(k1)!sumfac(k1+1)
    do j = 1, len00(k2)!sumfac(k2+1)
        
        mij = matr(i,j)
        
        gi  = gg(pi + i)
        gj  = gg(pj + j)
        gij = gg(pij + mij)
        
        !if(pri)print*,  pi + i, pj + j, pij + mij
        !if(pri)print*,  'gi, gj, gij', gi, gj, gij
        
        vout(mij) = vout(mij) + v1(i)*v2(j) * (gi*gj*cho)/gij
        
        if(pri)print*,  'i,j , v1(i), v2(j), h',i,j, v1(i), v2(j), (gi*gj*cho)/gij
        
        enddo
        enddo
    
end


subroutine test_old_field
    integer, parameter :: nm = 2
    real(dp) dpole(3,nm), qpole(3,3,nm), opole(3,3,3,nm), hpole(3,3,3,3,nm) 
    real(dp) d1v(3,nm), d2v(3,3,nm), d3v(3,3,3,nm), d4v(3,3,3,3,nm), d5v(3,3,3,3,3,nm)
    real(dp) a(3), a2(3), rCM(3,nm), fsf(3,nm), rMax2, rMax
    integer NC
    logical*1 iSlab
    logical FULL
    
    real(dp) quad(6), octa(10)!, hexa(15)
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    rMax = 100.1d0
    rMax2 = rMax**2
    rCM(:,1) = [0d0,0d0,0d0]
    rCM(:,2) = [3.4231, 2.74389, 1.54739]
    nc = 1
    a=40d0
    a2=a**2
    Full = .true. 
    iSlab = .false. 
    
    
    dpole=0
    qpole=0
    opole=0
    hpole=0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call random_seed(put=[2,234,1,5,435,4,5,42,3,43,432,4,3,5,23,345,34543])
    
    call random_number(quad)
    print'(a,*(g10.3))', 'quad:',quad
    !qpole(:,:,1) = reshape(expand(opdetr(quad,2),2),shape=[3,3])

    call random_number(octa)
    !opole(:,:,:,1) = reshape(expand(opdetr(octa,3),3),shape=[3,3,3])
    
    !call random_number(hexa)
    !hpole(:,:,:,:,1) = reshape(expand(opdetr(hexa,4),4),shape=[3,3,3,3])
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    dpole=0
    qpole=0
    opole=0
    hpole=0
    dpole(:,1) = [2.345d0, -0.453245d0,0.6564256d0]
    
    call calcDv(rCM, dpole, qpole, opole, hpole, nM, NC, a, a2, d1v, d2v, d3v, d4v, d5v, rMax2, fsf, iSlab,FULL)
    
    print'(a,*(g30.15))','d-1',d1v(:,2)
    print'(a,*(g30.15))','d-2',opdetr(compress( reshape(      d2v(:,:,2),shape=[3**2]),2),2)
    print'(a,*(g30.15))','d-3',opdetr(compress( reshape(    d3v(:,:,:,2),shape=[3**3]),3),3)
    print'(a,*(g30.15))','d-4',opdetr(compress( reshape(  d4v(:,:,:,:,2),shape=[3**4]),4),4)
    print'(a,*(g30.15))','d-5',opdetr(compress( reshape(d5v(:,:,:,:,:,2),shape=[3**5]),5),5)
    
    
    

    print*,'d1v'
    print'(1(g15.3))',d2v(:,:,2) !opdetr(compress( reshape(d1v(:,:,1),shape=[3**2]),2),2) !d2v(:,:,1) !
   
   print*, 'ABOVE: ---------------------test old field --------------------------'
     
end

subroutine test_potgrad
    integer, parameter :: kmax=5, nmax = 3
    integer i
    real(dp) rr(3) 
    real(dp) :: rrr(0:pos00(kmax+1)), quad(6), octa(10), rnorm, rsqe
    real(dp) :: rinvv(2*(kmax+nmax)+1) !, rinvv1(2*(kmax+nmax)+1), rinvv2(2*(kmax+nmax)+1)
    real(dp) :: dd(3), dr, d2v(3,3)
    integer be, ga
    
    rr = [3.4231, 2.74389, 1.54739]
    
    dd = [2.345d0, -0.453245d0,0.6564256d0]
    call  vector_powers(kmax,rr,rrr)
    !print*, rrr(0), rrr(1)
    !print'(a,*(e10.3))',"rrr:",rrr
    
    call random_seed(put=[2,234,1,5,435,4,5,42,3,43,432,4,3,5,23,345,34543])
    
    call random_number(quad)
    call random_number(octa)
    
    print'(a,*(g10.3))','quad:',quad
    print'(a,*(g10.3))','octa:',octa
    
    
    rsqe  = sum(rr**2)!dsqrt(rsq)
    rnorm = dsqrt(rsqe)
    rinvv(1) = 1d0/rnorm
    rinvv(2) = 1d0/rsqe
    
    do i = 3, 2*(kmax+nmax)+1
      rinvv(i) = rinvv(i-2)*rinvv(2)
    enddo
    
    dr = rnorm
    
    ! first dipole potential gradient
    print'(a,*(g30.15))', "d-1 exp", dd/dr**3 - (3*sum(rr*dd)/dr**5) * rr
    print'(a,*(g30.15))', 'd-1 new', potgrad(dd,1,1,rinvv,rrr)
    
    
    ! second dipole potential gradient
    d2v=0
    do be = 1, 3
      do ga = 1, 3 
        d2v(be,ga) =  - 3d0/dr**5 * ( dd(be)*rr(ga)  + dd(ga)*rr(be) + sum(dd*rr)*del(be,ga) ) & ! 
                      + & 
                      3d0*5d0*sum(dd*rr)*rr(be)*rr(ga) / dr**7
        enddo
        enddo
    
    print'(a,*(g30.15))', "d-2 exp", compress(reshape(d2v, shape=[3**2]),2)
    print'(a,*(g30.15))', 'd-2 new', opdetr(potgrad(dd,1,2,rinvv,rrr),2)
    
    
    
    print*, 'above, TEST POTGRAD-----------------------------------------------'
end

subroutine test_mp_pot
    integer, parameter :: nm = 2
    real(dp) dpole(3,nm), qpole(3,3,nm), opole(3,3,3,nm), hpole(3,3,3,3,nm) 
    real(dp) d1v(3,nm), d2v(3,3,nm), d3v(3,3,3,nm), d4v(3,3,3,3,nm), d5v(3,3,3,3,3,nm)
    real(dp) a(3), a2(3), rCM(3,nm), fsf(3,nm), rMax2, rMax
    integer NC
    logical*1 iSlab
    logical FULL
    
    integer, parameter :: kmax=6, nmax = 5
    integer i
    real(dp) rr(3) 
    real(dp) :: rrr(0:pos00(kmax+1)), rnorm, rsqe
    real(dp) :: rinvv(2*(kmax+nmax)+1) !, rinvv1(2*(kmax+nmax)+1), rinvv2(2*(kmax+nmax)+1)
    real(dp) :: dd(3), dr, cq(6), co(10), ch(15)
    !real(dp) :: quad(6), octa(10),
    
    integer n,k
    integer p1,p2,q1,q2
    real(dp) :: phi_old(pos00(6)),phi(pos00(6)), qq(pos00(5))
    
    
    print*, size(phi), sumfacfac(6), 3+6+10+21+15, size(qq), sumfacfac(5)
    
    rr = [3.4231, 2.74389, 1.54739]
    
    
    
    dd = [2.345d0, -0.453245d0,0.6564256d0]
    
    cq = [0.32534, 0.4352345, 1.5324, 1.2543, 1.35435, -1.57964]
    
    co = [0.4352345, 1.5324, 1.2543, 1.35435, -1.57964,0.32534, 0.4352345, 1.5324, 1.2543, 1.35435]
    co = opdetr(co,3)
    
    ch = [2.341,3.52345,3.2465,8.978,6.4356,7.77745,6.43563,7.73094589,3.421,3.4526,2.4564257,9.893543,3.464236,8.979,5.3452]
    ch = opdetr(ch,4)
    
    !dd = 0
    !cq = 0
    !co = 0
    !ch = 0
    
    
    
    qq = 0
    qq = [dd,cq,co,ch]
    
    !print*,"q-thigs"
    !n=1
    !q1 = pos00(n)+1
    !q2 = pos00(n+1)
    !print'(*(I3))',n, q1,q2
    !
    !qq(q1:q2) = dd
    !
    !n=2
    !q1 = pos00(n)+1
    !q2 = pos00(n+1)
    !print'(*(I3))',n, q1,q2
    !
    !qq(q1:q2) = cq
    
    
    print*
    print'(*(f10.4))', qq
    print*, 'size(qq)',size(qq), 3+6+10+15
    
    call  vector_powers(kmax,rr,rrr)
    
    rsqe  = sum(rr**2)!dsqrt(rsq)
    rnorm = dsqrt(rsqe)
    rinvv(1) = 1d0/rnorm
    rinvv(2) = 1d0/rsqe
    
    do i = 3, 2*(kmax+nmax)+1
      rinvv(i) = rinvv(i-2)*rinvv(2)
    enddo
    
    dr = rnorm
    
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    rMax = 100.1d0
    rMax2 = rMax**2
    rCM(:,1) = [0d0,0d0,0d0]
    rCM(:,2) = rr
    nc = 1
    a=40d0
    a2=a**2
    Full = .true. 
    iSlab = .false. 
    
    
    dpole=0
    qpole=0
    opole=0
    hpole=0
    
    
    hpole(:,:,:,:,1) = reshape(expand(ch,4),shape=[3,3,3,3])
    opole(:,:,:,1) = reshape(expand(co,3),shape=[3,3,3])
    qpole(:,:,1) = reshape(expand(cq,2),shape=[3,3])
    dpole(:,1) = dd(:)
    
    call calcDv(rCM, dpole, qpole, opole, hpole, nM, NC, a, a2, d1v, d2v, d3v, d4v, d5v, rMax2, fsf, iSlab,FULL)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
    k=1
    p1 = pos00(k)+1
    p2 = pos00(k+1)
    print'(6(I3))',k, p1,p2
    phi_old(p1:p2) = d1v(:,2)
    
    k=2
    p1 = pos00(k)+1
    p2 = pos00(k+1)
    print'(6(I3))',k, p1,p2
    phi_old(p1:p2) = opdetr(compress( reshape(      d2v(:,:,2),shape=[3**2]),2),2)
    
    k=3
    p1 = pos00(k)+1
    p2 = pos00(k+1)
    print'(6(I3))',k, p1,p2
    phi_old(p1:p2) = opdetr(compress( reshape(    d3v(:,:,:,2),shape=[3**3]),3),3)
    
    k=4
    p1 = pos00(k)+1
    p2 = pos00(k+1)
    print'(6(I3))',k, p1,p2
    phi_old(p1:p2) = opdetr(compress( reshape(  d4v(:,:,:,:,2),shape=[3**4]),4),4)
    
    k=5
    p1 = pos00(k)+1
    p2 = pos00(k+1)
    print'(6(I3))',k, p1,p2
    phi_old(p1:p2) = opdetr(compress( reshape(d5v(:,:,:,:,:,2),shape=[3**5]),5),5)
    
    
            
    
    phi=0
    do n = 1, 4
        q1 = pos00(n)+1
        q2 = pos00(n+1)
            do k = 1,5
                p1 = pos00(k)+1
                p2 = pos00(k+1)
                print'(6I3)',n,k,q1,q2, p1,p2
                phi(p1:p2) = phi(p1:p2) + opdetr(potgrad(qq(q1:q2),n,k,rinvv,rrr),k)
            enddo
    enddo
    
    
    
    
    print'(a,*(g30.15))', 'phi    ',phi
    print'(a,*(g30.15))', 'phi_old',phi_old
    print'(a,*(g30.15))', 'phi+old',phi_old+phi
    
    
    
    
    if(.false.)then
    
    !hexadeca
    print'(a,*(g30.15))', 'd-1 old', d1v(:,2)
    print'(a,*(g30.15))', 'd-1 new', potgrad(ch,4,1,rinvv,rrr)
    
    print'(a,*(g30.15))', 'd-2 old', opdetr(compress( reshape(      d2v(:,:,2),shape=[3**2]),2),2)
    print'(a,*(g30.15))', 'd-2 new', opdetr(potgrad(ch,4,2,rinvv,rrr),2)
    
    print'(a,*(g30.15))', 'd-3 old', opdetr(compress( reshape(    d3v(:,:,:,2),shape=[3**3]),3),3)
    print'(a,*(g30.15))', 'd-3 new', opdetr(potgrad(ch,4,3,rinvv,rrr),3)
    
    print'(a,*(g30.15))', 'd-4 old', opdetr(compress( reshape(  d4v(:,:,:,:,2),shape=[3**4]),4),4)
    print'(a,*(g30.15))', 'd-4 new', opdetr(potgrad(ch,4,4,rinvv,rrr),4)
    
    print'(a,*(g30.15))', 'd-5 old', opdetr(compress( reshape(d5v(:,:,:,:,:,2),shape=[3**5]),5),5)
    print'(a,*(g30.15))', 'd-5 new', opdetr(potgrad(ch,4,5,rinvv,rrr),5)
    
    !octu
    print'(a,*(g30.15))', 'd-1 old', d1v(:,2)
    print'(a,*(g30.15))', 'd-1 new', potgrad(co,3,1,rinvv,rrr)
    
    print'(a,*(g30.15))', 'd-2 old', opdetr(compress( reshape(      d2v(:,:,2),shape=[3**2]),2),2)
    print'(a,*(g30.15))', 'd-2 new', opdetr(potgrad(co,3,2,rinvv,rrr),2)
    
    print'(a,*(g30.15))', 'd-3 old', opdetr(compress( reshape(    d3v(:,:,:,2),shape=[3**3]),3),3)
    print'(a,*(g30.15))', 'd-3 new', opdetr(potgrad(co,3,3,rinvv,rrr),3)
    
    print'(a,*(g30.15))', 'd-4 old', opdetr(compress( reshape(  d4v(:,:,:,:,2),shape=[3**4]),4),4)
    print'(a,*(g30.15))', 'd-4 new', opdetr(potgrad(co,3,4,rinvv,rrr),4)
    
    print'(a,*(g30.15))', 'd-5 old', opdetr(compress( reshape(d5v(:,:,:,:,:,2),shape=[3**5]),5),5)
    print'(a,*(g30.15))', 'd-5 new', opdetr(potgrad(co,3,5,rinvv,rrr),5)
    
    !quadru
    print'(a,*(g30.15))', 'd-2 old', opdetr(compress( reshape(      d2v(:,:,2),shape=[3**2]),2),2)
    print'(a,*(g30.15))', 'd-2 new', opdetr(potgrad(dd,1,2,rinvv,rrr),2)
    
    print'(a,*(g30.15))', 'd-3 old', opdetr(compress( reshape(    d3v(:,:,:,2),shape=[3**3]),3),3)
    print'(a,*(g30.15))', 'd-3 new', opdetr(potgrad(dd,1,3,rinvv,rrr),3)
    
    print'(a,*(g30.15))', 'd-4 old', opdetr(compress( reshape(  d4v(:,:,:,:,2),shape=[3**4]),4),4)
    print'(a,*(g30.15))', 'd-4 new', opdetr(potgrad(dd,1,4,rinvv,rrr),4)
    
    print'(a,*(g30.15))', 'd-5 old', opdetr(compress( reshape(d5v(:,:,:,:,:,2),shape=[3**5]),5),5)
    print'(a,*(g30.15))', 'd-5 new', opdetr(potgrad(dd,1,5,rinvv,rrr),5)
    
    !dip
    print'(a,*(g30.15))', 'd-1 old', d1v(:,2)
    print'(a,*(g30.15))', 'd-1 new', potgrad(cq,2,1,rinvv,rrr)
    
    print'(a,*(g30.15))', 'd-2 old', opdetr(compress( reshape(      d2v(:,:,2),shape=[3**2]),2),2)
    print'(a,*(g30.15))', 'd-2 new', opdetr(potgrad(cq,2,2,rinvv,rrr),2)
    
    print'(a,*(g30.15))', 'd-3 old', opdetr(compress( reshape(    d3v(:,:,:,2),shape=[3**3]),3),3)
    print'(a,*(g30.15))', 'd-3 new', opdetr(potgrad(cq,2,3,rinvv,rrr),3)
    
    print'(a,*(g30.15))', 'd-4 old', opdetr(compress( reshape(  d4v(:,:,:,:,2),shape=[3**4]),4),4)
    print'(a,*(g30.15))', 'd-4 new', opdetr(potgrad(cq,2,4,rinvv,rrr),4)
    
    print'(a,*(g30.15))', 'd-5 old', opdetr(compress( reshape(d5v(:,:,:,:,:,2),shape=[3**5]),5),5)
    print'(a,*(g30.15))', 'd-5 new', opdetr(potgrad(cq,2,5,rinvv,rrr),5)
    
    endif
    

    
   print*, 'ABOVE: dip: grad vs old --------------------------'
     
end



function potgrad(qq,nn,kk,rpows,rrr) 
    real(dp), intent(in) :: qq((nn+1)*(nn+2)/2),rpows(:),rrr(0:)
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
        
        k1 = pos00(k_i)+1
        k2 = pos00(k_i+1)
        n1 = pos00(n_i)+1
        n2 = pos00(n_i+1)
        
        potgrad(:) = potgrad(:) &
                    + (-1)**sig * intfac(nn,n_i) * intff(kni2-1,2*nn-1) * rpows(kni2+1) & !why +nn in first term???
                    * symouter(k_i, i, rrr(k1:k2),   inner(nn,n_i, qq, rrr(n1:n2))   )
    enddo
    
end
  !integer, parameter :: rpos(8) = [1,      4,     10,     20,     35,     56,     84,    120] ! remove this, use pos00
    
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
   




subroutine test_matr(nn)
    integer i, nn, maxi
    if(nn>7)stop"rank cant be larger than 7 in matr intdex matrix"
    maxi = sumfac(nn)
    do i = 1, maxi
      print'(28I3)', matr(i,1:maxi)
    enddo
end subroutine

!subroutine test_matri(nn)
!    integer i, nn, maxi
!    if(nn>7)stop"rank cant be larger than 7 in matr intdex matrix"
!    maxi = sumfac(nn)
!    do i = 1, maxi
!      print'(28I3)', matri(i,1:maxi) - matr(i,1:maxi)
!    enddo
!end subroutine





subroutine test_rrr
    integer, parameter :: k=5
    real(dp) :: rrr(sumfacfac(k+1)-1), rr(3), rrr3(10), rrr2(6), rrr5(sumfac(6)), rrr32(sumfac(6))
    integer i, p1, p2, p3, p4
    call random_seed(put=[2,234,1,5,435,4,5,42,3,43,432,4,3,5,23,345,34543])
    call random_number(rr)
    rr = [1d0,2d0,3d0]
    rr = [0.3810985945,0.5295087287,0.852367145402366]
    print*, size(rrr)
    
    call vector_powers(k,rr,rrr)
    !call rrpow(rr,k,rrr)
    !print'(*(f12.4))', rrr  
    call printer(rrr,'rrr',1)
    
    
    i = 3
    p1 = pos00(i-1)+1
    p2 = pos00(i)
    p3 = pos00(i)+1
    p4 = pos00(i+1)
    
    rrr2=rrr(p1:p2)
    rrr3=rrr(p3:p4)
    
    i=5
    p3 = pos00(i)+1
    p4 = pos00(i+1)
    
    rrr5=rrr(p3:p4)
    
    rrr32 = symouter(3,2,rrr3,rrr2)
    
    
    call printer(rrr5,"5th",1)
    call printer(rrr32,"symouter(3,2)",1)
    call printer(rrr5-rrr32,"diff",1)
    
    print*, "above, TEST RRR ------------------------------------------------------------------"

end


subroutine vector_powers(k,r,rr) 
    ! Do not use this routine, it does not produce th outer products, since it lacks the scaling parameter hhh
    ! Generates k powers of r(3) in tricorn polytensor
    ! r(3)   position vector
    ! k      highest rank/power
    ! rr     output tricorn-ordered polytensor
    integer, intent(in)   :: k
    real(dp), intent(in)  :: r(3)
    real(dp), intent(out) :: rr(0:sumfacfac(k+1)-1)
    integer pl,cl,px,cx,i,pz,cz, py, cy
    
    rr(0) = 1
    rr(1:3) = r(:)
    
    do i = 2,k
       ! p=previous, c=current, l=length
       ! x, y, z refer to the position of the x-only, y-only, z-only rows. 
       px = pos00(i-1)+1 !simfacfac(i-1)
       pl = len00(i-1) !sumfac(i)
       cl = len00(i) ! sumfac(i+1)
       
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
