module compressed_tensors

use printer_mod, only: str, printer, printo
use compressed_utils, bad=>main!,only: test_apple_g
implicit none

integer, parameter :: dp = kind(0d0)

integer j11,j22, i00, matcol, matrow
integer, parameter :: matk = 7, matsize = (matk+1)*(matk+2)/2
integer, parameter :: sumfacv(matk+1) = [(i00*(i00+1)/2, i00 = 1,matk+1)]
integer, parameter :: vv(matsize) = [((j11, j22=1,j11),j11 = 1,matk+1)]
integer, parameter :: vv1(matsize) = vv-1

integer, parameter :: gg(285) = &
[  1,   1,   1, &
   1,   2,   2,   1,   2,   1, &
   1,   3,   3,   3,   6,   3,   1,   3,   3,   1, &
   1,   4,   4,   6,  12,   6,   4,  12,  12,   4,   1,   4,   6,   4,   1, &
   1,   5,   5,  10,  20,  10,  10,  30,  30,  10,   5,  20,  30,  20,   5,   1,   5,  10,  10,   5,   1, &
   1,   6,   6,  15,  30,  15,  20,  60,  60,  20,  15,  60,  90,  &
       60,  15,   6,  30,  60,  60,  30,   6,   1,   6,  15,  20,  15,   6,   1, &
   1,   7,   7,  21,  42,  21,  35, 105, 105,  35,  35, 140, 210, 140,  35,  21, 105, 210, 210, 105,  21, &
        7,  42, 105, 140, 105,  42,   7,   1,   7,  21,  35,  35,  21,   7,   1, &
   1,   8,   8,  28,  56,  28,  56, 168, 168,  56,  70, 280, 420, 280,  70,  56, 280, 560, 560, 280,  56,  &
       28, 168, 420, 560, 420, 168,  28,   8,  56, 168, 280, 280, 168,  56,   8,   1,   8,  28,  56,  70,  &
       56,  28,   8,   1, &
   1,    9,    9,   36,   72,   36,   84,  252,  252,   84,  126,  504,  756,  504,  126,  126,  630, 1260, &
      1260,  630,  126,   84,  504, 1260, 1680, 1260,  504,   84,   36,  252,  756, 1260, 1260,  756,  252, &
        36,    9,   72,  252,  504,  630,  504,  252,   72,    9,    1,    9,   36,   84,  126,  126,   84,   36,    9,    1, &
   1,   10,   10,   45,   90,   45,  120,  360,  360,  120,  210,  840, 1260,  840,  210,  252, 1260, 2520, 2520, 1260,  252, &
       210, 1260, 3150, 4200, 3150, 1260,  210,  120,  840, 2520, 4200, 4200, 2520,  840,  120,   45,  360, 1260, 2520, 3150, &
      2520, 1260,  360,   45,   10,   90,  360,  840, 1260, 1260,  840,  360,   90,   10,    1,   10,   45,  120,  210,  252, & 
       210,  120,   45,   10,    1 &
]

integer, parameter :: matchoo(10,10) = reshape( [&
     1,      2,      3,      4,      5,      6,      7,      8,      9,     10,  &
     2,      1,      3,      6,     10,     15,     21,     28,     36,     45,  &
     3,      3,      1,      4,     10,     20,     35,     56,     84,    120,  &
     4,      6,      4,      1,      5,     15,     35,     70,    126,    210,  &
     5,     10,     10,      5,      1,      6,     21,     56,    126,    252,  &
     6,     15,     20,     15,      6,      1,      7,     28,     84,    210,  &
     7,     21,     35,     35,     21,      7,      1,      8,     36,    120,  &
     8,     28,     56,     70,     56,     28,      8,      1,      9,     45,  &
     9,     36,     84,    126,    126,     84,     36,      9,      1,     10,  &
    10,     45,    120,    210,    252,    210,    120,     45,     10,      1   &
 ], shape(matchoo) )
integer, parameter :: gposarr(10) = [1, 4, 10, 20, 35, 56, 84, 120, 165, 220]
integer, parameter :: gpos(10) = gposarr(1:10)-1

integer, parameter :: matr(28,28) = reshape( [&
 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28,&
 2,  4,  5,  7,  8,  9, 11, 12, 13, 14, 16, 17, 18, 19, 20, 22, 23, 24, 25, 26, 27, 29, 30, 31, 32, 33, 34, 35,&
 3,  5,  6,  8,  9, 10, 12, 13, 14, 15, 17, 18, 19, 20, 21, 23, 24, 25, 26, 27, 28, 30, 31, 32, 33, 34, 35, 36,&
 4,  7,  8, 11, 12, 13, 16, 17, 18, 19, 22, 23, 24, 25, 26, 29, 30, 31, 32, 33, 34, 37, 38, 39, 40, 41, 42, 43,&
 5,  8,  9, 12, 13, 14, 17, 18, 19, 20, 23, 24, 25, 26, 27, 30, 31, 32, 33, 34, 35, 38, 39, 40, 41, 42, 43, 44,&
 6,  9, 10, 13, 14, 15, 18, 19, 20, 21, 24, 25, 26, 27, 28, 31, 32, 33, 34, 35, 36, 39, 40, 41, 42, 43, 44, 45,&
 7, 11, 12, 16, 17, 18, 22, 23, 24, 25, 29, 30, 31, 32, 33, 37, 38, 39, 40, 41, 42, 46, 47, 48, 49, 50, 51, 52,&
 8, 12, 13, 17, 18, 19, 23, 24, 25, 26, 30, 31, 32, 33, 34, 38, 39, 40, 41, 42, 43, 47, 48, 49, 50, 51, 52, 53,&
 9, 13, 14, 18, 19, 20, 24, 25, 26, 27, 31, 32, 33, 34, 35, 39, 40, 41, 42, 43, 44, 48, 49, 50, 51, 52, 53, 54,&
10, 14, 15, 19, 20, 21, 25, 26, 27, 28, 32, 33, 34, 35, 36, 40, 41, 42, 43, 44, 45, 49, 50, 51, 52, 53, 54, 55,&
11, 16, 17, 22, 23, 24, 29, 30, 31, 32, 37, 38, 39, 40, 41, 46, 47, 48, 49, 50, 51, 56, 57, 58, 59, 60, 61, 62,&
12, 17, 18, 23, 24, 25, 30, 31, 32, 33, 38, 39, 40, 41, 42, 47, 48, 49, 50, 51, 52, 57, 58, 59, 60, 61, 62, 63,&
13, 18, 19, 24, 25, 26, 31, 32, 33, 34, 39, 40, 41, 42, 43, 48, 49, 50, 51, 52, 53, 58, 59, 60, 61, 62, 63, 64,&
14, 19, 20, 25, 26, 27, 32, 33, 34, 35, 40, 41, 42, 43, 44, 49, 50, 51, 52, 53, 54, 59, 60, 61, 62, 63, 64, 65,&
15, 20, 21, 26, 27, 28, 33, 34, 35, 36, 41, 42, 43, 44, 45, 50, 51, 52, 53, 54, 55, 60, 61, 62, 63, 64, 65, 66,&
16, 22, 23, 29, 30, 31, 37, 38, 39, 40, 46, 47, 48, 49, 50, 56, 57, 58, 59, 60, 61, 67, 68, 69, 70, 71, 72, 73,&
17, 23, 24, 30, 31, 32, 38, 39, 40, 41, 47, 48, 49, 50, 51, 57, 58, 59, 60, 61, 62, 68, 69, 70, 71, 72, 73, 74,&
18, 24, 25, 31, 32, 33, 39, 40, 41, 42, 48, 49, 50, 51, 52, 58, 59, 60, 61, 62, 63, 69, 70, 71, 72, 73, 74, 75,&
19, 25, 26, 32, 33, 34, 40, 41, 42, 43, 49, 50, 51, 52, 53, 59, 60, 61, 62, 63, 64, 70, 71, 72, 73, 74, 75, 76,&
20, 26, 27, 33, 34, 35, 41, 42, 43, 44, 50, 51, 52, 53, 54, 60, 61, 62, 63, 64, 65, 71, 72, 73, 74, 75, 76, 77,&
21, 27, 28, 34, 35, 36, 42, 43, 44, 45, 51, 52, 53, 54, 55, 61, 62, 63, 64, 65, 66, 72, 73, 74, 75, 76, 77, 78,&
22, 29, 30, 37, 38, 39, 46, 47, 48, 49, 56, 57, 58, 59, 60, 67, 68, 69, 70, 71, 72, 79, 80, 81, 82, 83, 84, 85,&
23, 30, 31, 38, 39, 40, 47, 48, 49, 50, 57, 58, 59, 60, 61, 68, 69, 70, 71, 72, 73, 80, 81, 82, 83, 84, 85, 86,&
24, 31, 32, 39, 40, 41, 48, 49, 50, 51, 58, 59, 60, 61, 62, 69, 70, 71, 72, 73, 74, 81, 82, 83, 84, 85, 86, 87,&
25, 32, 33, 40, 41, 42, 49, 50, 51, 52, 59, 60, 61, 62, 63, 70, 71, 72, 73, 74, 75, 82, 83, 84, 85, 86, 87, 88,&
26, 33, 34, 41, 42, 43, 50, 51, 52, 53, 60, 61, 62, 63, 64, 71, 72, 73, 74, 75, 76, 83, 84, 85, 86, 87, 88, 89,&
27, 34, 35, 42, 43, 44, 51, 52, 53, 54, 61, 62, 63, 64, 65, 72, 73, 74, 75, 76, 77, 84, 85, 86, 87, 88, 89, 90,&
28, 35, 36, 43, 44, 45, 52, 53, 54, 55, 62, 63, 64, 65, 66, 73, 74, 75, 76, 77, 78, 85, 86, 87, 88, 89, 90, 91&
 ], shape(matr) )

integer, parameter :: matri(matsize,matsize) = reshape( [ &
                      ((matrow + matcol -1 + vv1(matcol)*vv1(matrow), matcol = 1,matsize), matrow = 1, matsize) &
                      ], shape(matri) )


contains !///////////////////////////////////////////////





subroutine main
    !call testing
    !call test_sumfac
    !call test_rpow
    !call test_nextpow(7)
    !call compow2(7)
    !call subdiv_pow(5,5)
    !call test_contract
    call test_produce
    
    !print*, ""
    !call subdiv_pow(4,3)
    !print*, ""
    !call subdiv_pow(3,4)
    !print*,"fac 5 / 2 3", fac(5)/( fac(3)*fac(2))
    !call test_factorial(5)
    
    !call test_choose
    !call test_matri(7)
    
    !print*, ""
    !call test_contind(6)
    
    !call test_matr(7)
    !call printer(matr, 'matr',1)
    !call test_apple_g
    !call test_rrpow
    !call test_next_rev_key2n(4)
    !integer i, j
    !do i = 1,10
    !  do j = 1,10
    !    write(*,'(I6,a)', advance="no") choose(i,j),", "
    !    enddo
    !  print*, ""
    !  enddo
    
    
end subroutine

! Testing /////////////////////////////////////////////////////

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
      !                 ( gg(gpos(ii)+i)*gg(gpos(kii)+j)*choose(ii+kii,ii) ) / gg(gpos(ii+kii) + matr(i,j)), & 
      !                 00,&
      !                 apple_g([a1,b1,c1]), apple_g([a2,b2,c2]), choose(ii+kii,ii), apple_g([aa,bb,cc]), &
      !                 00, &
      !                 gg(gpos(ii)+i), gg(gpos(kii)+j), choose(ii+kii,ii), gg(gpos(ii+kii) + matr(i,j)), & 
      !                 00
      
      !write(*,'(I6)', advance="no") ( apple_g([a1,b1,c1])*apple_g([a2,b2,c2])*choose(ii+kii,ii) )/apple_g([aa,bb,cc])
      !write(*,'(I6)', advance="no") ( gg(gpos(ii)+i)*gg(gpos(kii)+j)*choose(ii+kii,ii) ) / gg(gpos(ii+kii) + matr(i,j))
      
      
      !write(*,'(I6)', advance="no") ( gg(gpos(ii)+i)*gg(gpos(kii)+j)*choose(ii+kii,ii) ) / gg(gpos(ii+kii) + matr(i,j)) &
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
    
    hh = fac(aa)*fac(bb)*fac(cc) / ( fac(a1)*fac(a2)*fac(b1)*fac(b2)*fac(c1)*fac(c2) )

end

pure function hhh(i,j,ki,kj) !Make a matrix of this sheeeet
    integer, intent(in) :: i,j, ki,kj
    integer kij, hhh, cc, gi,gj,gij
    kij = ki+kj
    cc = matchoo(kij,ki)
    gi = gg(gpos(ki)+i)
    gj = gg(gpos(kj)+j) 
    gij = gg(gpos(kij) + matr(i,j))
    hhh = (gi*gj*cc) / gij
    
end

subroutine test_contract
    real(dp) :: v1(3), v2(6), v3(10), v4(15), v5(21), ov1(3), ov2(6)
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
    
    
    print'(*(f10.4))', of1
    print'(*(f10.4))', contract(3,2,v3,v2)
    print*,""
    print'(*(f10.4))', compress(reshape(of2,[3**2]),2)
    print'(*(f10.4))', contract(4,2,v4,v2)
    print*,""
    print'(*(f10.4))', compress(reshape(of3,[3**3]),3)
    print'(*(f10.4))', contract(5,2,v5,v2)
    
    !print'(*(f10.4))', contract(5,2,v5,v2)
    
    !print'(*(f10.4))', contract(3,2,v3,v2)-of1
    
    
end

function contract(kbig, ksmall, vbig, vsmall) result(vout) !assumes kbig > ksmall
   integer, intent(in) :: kbig, ksmall
   real(dp), intent(in) :: vbig(:), vsmall(:)
   !real(dp), pointer :: p1, p2 !, target
   real(dp) vout((kbig - ksmall+1)*(kbig - ksmall+2)/2)
   integer kout, gposition, i, j, maxi, maxj
   
   kout = kbig - ksmall
   maxi = (kout+1)*(kout+2)/2
   maxj = (ksmall+1)*(ksmall+2)/2
   
   gposition = gpos(ksmall)
   
   vout(:) = 0
   do i = 1, maxi
     do j = 1,maxj
       vout(i) = vout(i) + vbig(matri(i,j)) * vsmall(j) * gg(gposition+j)
       enddo
       enddo
   
end

subroutine test_produce
    real(dp) :: v1(3), v2(6), v3(10), v4(15), v5(21), ov1(3), ov2(6)
    real(dp) :: f1(3), f2(3,3), f3(3,3,3), f4(3,3,3,3), f5(3,3,3,3,3), of1(3), of2(3,3), of3(3,3,3), of4(3,3,3,3)
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
    
    print*,""
    print'(*(f10.4))', compress(reshape(of3,[3**3]),3)
    print'(*(f10.4))', produce(1,2,v1,v2)
    print*,""
    print'(*(f10.4))', compress(reshape(of4,[3**4]),4)
    print'(*(f10.4))', produce(1,3,v1,v3)
    
    
    
end

function produce(k1,k2,v1,v2) result(vout)
integer, intent(in) :: k1,k2
real(dp), intent(in) :: v1(:), v2(:)
real(dp) vout( (k1+k2+1)*(k1+k2+2)/2 )
integer i, j, k, kout, l1, l2, lout

kout = k1+k2

lout = (kout+1)*(kout+2)/2
l1   = (k1  +1)*(k1  +2)/2
l2   = (k2  +1)*(k2  +2)/2

vout=0
do i = 1, l1
do j = 1, l2
  
  vout(matr(i,j)) = vout(matr(i,j)) + v1(i)*v2(j)*hhh(i,j,k1,k2)
  enddo
  enddo
  

end




subroutine test_matr(nn)
    integer i, nn, maxi
    if(nn>7)stop"rank kant be larger than 7 in matr intdex matrix"
    maxi = sumfac(nn)
    do i = 1, maxi
      print'(28I3)', matr(i,1:maxi)
    enddo
end subroutine

subroutine test_matri(nn)
    integer i, nn, maxi
    if(nn>7)stop"rank kant be larger than 7 in matr intdex matrix"
    maxi = sumfac(nn)
    do i = 1, maxi
      print'(28I3)', matri(i,1:maxi) - matr(i,1:maxi)
    enddo
end subroutine


subroutine compow2(k)
  integer a,b,c, k
  integer too,i, it
  
  it=0
  do too = 0,k
   do i = 0,too
     a=k-too
     b=too-i
     c=i
     it=it+1!((too+1)*too)/2+i+1
     print'(3I2,I5)',a,b,c,it
     enddo
     enddo
end



subroutine test_factorial(k)
integer k, i
do i = 1, k
  print*, fac(i)
  enddo

end

function fac(inn) 
    integer inn, i, fac
    fac=1
    do i = 1, inn
      fac = fac*i
      enddo

end

subroutine test_contind(k)
   implicit integer (i)
   integer j1,j2,j3,j4, j,k, ind(10), ii(2), jj(4), v1(sumfac(k+1)),v2(sumfac(k+1)),v3(sumfac(k+1)) !((k+1)*(k+2)/2),
   real*8 rr(2,3)
   
   rr(1,:) = [1d0,2d0,3d0]
   rr(2,:) = [4d0,5d0,6d0]
   jj(:) = [1,2,3,1]
   ii = [1,2]
   print*, rr(ii,jj)
   
   ileng=sumfac(k+1)
   
   v1=0
   v2=0
   
   v1 = [(i,i=1,ileng)]
   v2 = [((j1, j2=1,j1),j1 = 1,k+1)]
   
   print*, ''
   print'(*(I3))', v1
   print'(*(I3))', v2
   
   print*, ""
   do irow = 1, ileng
       v3 = [(irow + icol -1 + (v2(icol)-1)*(v2(irow)-1), icol = 1,ileng)]
     print'(*(I3))',v3
     enddo
   
!   do iv = 1,ileng
!       ijumps = v2(iv) - 1
!       iones = iv - ijumps - 1
!       
!       print'(*(I3))', v1 + ijumps*v2 + iones
!       
!       enddo
   
!   print*, ""
!   do irow = 1, ileng
!     do icol = 1, ileng
!       v3(icol) = irow + icol -1 + (v2(icol)-1)*(v2(irow)-1)
!       enddo
!     print'(*(I3))',v3
!     enddo
   
!   print*, ""
!   do irow = 1, ileng
!     iset = v2(irow)
!     do icol = 1, ileng
!       ijumps = v2(icol)-1
!       iones = icol - ijumps -1
!       v3(icol) = irow + iones + ijumps*iset
!       
!       
!       enddo
!     print'(*(I3))',v3
!     enddo
     
      
   !print'(*(I3))', v1
   !print'(*(I3))', v1+v2
   !print'(*(I3))', v1+v2+1
   !print'(*(I3))', v1+2*v2+1
   !print'(*(I3))', v1+2*v2+2
   
   
   !do i = 1,6
   
   
   
   
   
   
end subroutine

subroutine test_rrpow
    !call testing
    !call test_sumfac
    !call test_rpow

    real(dp) :: rr(7*8-1)
    call rrpow([1d0,2d0,3d0],5,rr)
    call printer(rr,'rr',1)
    
    !call test_next_rev_key2n(4)
end subroutine


subroutine test_rpow
    real(dp) :: r2(6),r3(10),r4(15),r5(21), r(3)
    r(:) = [1d0,2d0,3d0]

    call rpow(r,r2,r3,r4,r5)
call printer(r,'',1)
call printer(r2,'2',1)
call printer(r3,'3',1)
call printer(r4,'4',1)
call printer(r5,'5',1)

endsubroutine










! END TESTING //////////////////////////////////////////////////////////////




!///////////////////////////////////////////////////////////////////////////////////////


subroutine rrpow(r,k,rr) 
    ! Generates k powers of r(3) in tricorn polytensor
    ! r(3)   position vector
    ! k      highest rank/power
    ! rr     output tricorn-ordered polytensor
    integer, intent(in)   :: k
    real(dp), intent(in)  :: r(3)
    real(dp), intent(out) :: rr(sumfacfac(k+1)-1)
    integer pl,cl,px,cx,i,pz,cz, py, cy
    
    rr(1:3) = r(:)
    
    do i = 2,k
       ! p=previous, c=current, l=length
       ! x, y, z refer to the position of the x-only, y-only, z-only rows. 
       px = sumfacfac(i-1)
       pl = sumfac(i)
       cl = sumfac(i+1)
       
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

subroutine rpow(r,r2,r3,r4,r5) 
!subroutine rpow(r,r2,r3,r4,r5) 
! Generates powers of r(3) in tricorn vecotrs
real(dp) :: r(3), r2(6),r3(10),r4(15),r5(21)





!i=2; pl = sumfac(i); cl = sumfac(i+1)
! r2(1:pl)      = r(1:pl)     *r(1)
! r2(pl+1:cl-1) = r(pl-i+1:pl)*r(2)
! r2(cl)        = r(pl)       *r(3)
!
!i=3; pl = sumfac(i); cl = sumfac(i+1)
! r3(1:pl)      = r2(1:pl)     *r(1)
! r3(pl+1:cl-1) = r2(pl-i+1:pl)*r(2)
! r3(cl)        = r2(pl)       *r(3)
!
!
!i=4; pl = sumfac(i); cl = sumfac(i+1)
! r4(1:pl)      = r3(1:pl)     *r(1)
! r4(pl+1:cl-1) = r3(pl-i+1:pl)*r(2)
! r4(cl)        = r3(pl)       *r(3)
!
!i=5; pl = sumfac(i); cl = sumfac(i+1)
! r5(1:pl)      = r4(1:pl)     *r(1)
! r5(pl+1:cl-1) = r4(pl-i+1:pl)*r(2)
! r5(cl)        = r4(pl)       *r(3)

! 2
r2(1:3) = r(1:3)*r(1)
r2(4:5) = r(2:3)*r(2)
r2(6)   = r(3)  *r(3)


! 3
r3(1:6) = r2(1:6)*r(1)
r3(7:9) = r2(4:6)*r(2)
r3(10)  = r2(6)  *r(3)

!4
r4(1:10)  = r3(1:10)*r(1)
r4(11:14) = r3(7:10)*r(2)
r4(15)    = r3(10)  *r(3)

!5
r5(1:15)  = r4(1:15)  *r(1)
r5(16:20) = r4(11:15) *r(2)
r5(21)    = r4(15)    *r(3)

end





!function contract(qq,qkk,rr,rkk,qr,qr_kk,sym)
!integer qkk, rkk, qr_kk
!real(dp) :: qq(sumfac(qkk+1)), rr(sumfac(rkk+1)), qr(sumfac(qr_kk+1))
!integer k, dqkk, drkk, len3, i, j, lenh
!k = qkk+rkk-qr_kk
!dqkk = qkk-k
!drkk = rkk-k
!
!len3 = sumfac(qr_kk)
!lenh = sumfac(k)
!
!perms = choose(max(dqkk,drkk),min(dqkk,drkk))
!
!do ip = 1,perms
!
!do i = 1, len3
!  
!  qr(i) = 0
!  do j = 1, lenh
!    
!    qr(i) = qr(i) + 
!  
!  
!enddo
!end subroutine































!/////////////////////////////////////////////////



end module
