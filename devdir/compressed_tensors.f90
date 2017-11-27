module compressed_tensors

use printer_mod, only: str, printer, printo
use compressed_utils, bad=>main!,only: test_apple_g

use detrace_apple, bad=>main

use calc_derivs, only:calcDv
implicit none

integer, parameter :: dp = kind(0d0)

integer j11,j22, i00, matcol, matrow
integer, parameter :: matk = 7, matsize = (matk+1)*(matk+2)/2
integer, parameter :: sumfacv(matk+1) = [(i00*(i00+1)/2, i00 = 1,matk+1)]
integer, parameter :: vv(matsize) = [((j11, j22=1,j11),j11 = 1,matk+1)]
integer, parameter :: vv1(matsize) = vv-1

integer, parameter :: gg(0:285) = &
[  1, & ! for the 0th index! Fucking shit should use 0-indexing!
   1,   1,   1, &
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
integer, parameter :: gposarr(0:10) = [0, 1, 4, 10, 20, 35, 56, 84, 120, 165, 220]
integer, parameter :: gpos(0:10) = gposarr(0:10)-1

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
    
    
    !call test_old_field
    call test_potgrad
    !call test_inner
    
    call test_dip_pot
    !call test_intfac_ff
    
    
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
    
    hh = ( fac(aa)*fac(bb)*fac(cc) ) / ( fac(a1)*fac(a2)*fac(b1)*fac(b2)*fac(c1)*fac(c2) )

end

subroutine test_hhh(k1, k2)
    integer i1, i2, k1,k2
    integer gp1, gp2, gp12, h, ch
    gp1 = gpos(k1)
    gp2 = gpos(k2)
    gp12 = gpos(k1+k2)
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
    gi = gg(gpos(ki)+i)
    gj = gg(gpos(kj)+j) 
    gij = gg(gpos(kij) + matr(i,j))
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
   
   gplace = gpos(kr)
   
   vqr(:) = 0
   do i = 1, maxi
     do j = 1,maxj
       vqr(i) = vqr(i) + vq(matri(i,j)) * vr(j) * gg(gplace+j)
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
    cho = choose(k12,k1)
    
    pri=.false.!.true.
    if(pri)print*, "choose", cho
    
    pi  = gpos(k1)
    pj  = gpos(k2)
    pij = gpos(k12)
    if(pri)print*, "pi, pj, pij", pi, pj, pij
    
    vout=0
    
    do i = 1, sumfac(k1+1)
    do j = 1, sumfac(k2+1)
        
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
    real(dp) :: rrr(0:gpos(kmax+1)), quad(6), octa(10), rnorm, rsqe
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

subroutine test_dip_pot
    integer, parameter :: nm = 2
    real(dp) dpole(3,nm), qpole(3,3,nm), opole(3,3,3,nm), hpole(3,3,3,3,nm) 
    real(dp) d1v(3,nm), d2v(3,3,nm), d3v(3,3,3,nm), d4v(3,3,3,3,nm), d5v(3,3,3,3,3,nm)
    real(dp) a(3), a2(3), rCM(3,nm), fsf(3,nm), rMax2, rMax
    integer NC
    logical*1 iSlab
    logical FULL
    
    integer, parameter :: kmax=5, nmax = 3
    integer i
    real(dp) rr(3) 
    real(dp) :: rrr(0:gpos(kmax+1)), rnorm, rsqe
    real(dp) :: rinvv(2*(kmax+nmax)+1) !, rinvv1(2*(kmax+nmax)+1), rinvv2(2*(kmax+nmax)+1)
    real(dp) :: dd(3), dr
    !real(dp) :: quad(6), octa(10),
    
    rr = [3.4231, 2.74389, 1.54739]
    
    dd = [2.345d0, -0.453245d0,0.6564256d0]
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
    
    
    dpole=0
    qpole=0
    opole=0
    hpole=0
    dpole(:,1) = dd(:)
    
    call calcDv(rCM, dpole, qpole, opole, hpole, nM, NC, a, a2, d1v, d2v, d3v, d4v, d5v, rMax2, fsf, iSlab,FULL)
    
    print'(a,*(g30.15))', 'd-1 old', d1v(:,2)
    print'(a,*(g30.15))', 'd-1 new', potgrad(dd,1,1,rinvv,rrr)
    
    print'(a,*(g30.15))', 'd-2 old', opdetr(compress( reshape(      d2v(:,:,2),shape=[3**2]),2),2)
    print'(a,*(g30.15))', 'd-2 new', opdetr(potgrad(dd,1,2,rinvv,rrr),2)
    
    print'(a,*(g30.15))', 'd-3 old', opdetr(compress( reshape(    d3v(:,:,:,2),shape=[3**3]),3),3)
    print'(a,*(g30.15))', 'd-3 new', opdetr(potgrad(dd,1,3,rinvv,rrr),3)
    
    print'(a,*(g30.15))', 'd-4 old', opdetr(compress( reshape(  d4v(:,:,:,:,2),shape=[3**4]),4),4)
    print'(a,*(g30.15))', 'd-4 new', opdetr(potgrad(dd,1,4,rinvv,rrr),4)
    
    print'(a,*(g30.15))', 'd-5 old', opdetr(compress( reshape(d5v(:,:,:,:,:,2),shape=[3**5]),5),5)
    print'(a,*(g30.15))', 'd-5 new', opdetr(potgrad(dd,1,5,rinvv,rrr),5)
    
    

    
   print*, 'ABOVE: dip: grad vs old --------------------------'
     
end



function potgrad(qq,nq,kk,rpows,rrr) 
    real(dp), intent(in) :: qq((nq+1)*(nq+2)/2),rpows(:),rrr(0:)
    integer, intent(in) :: nq, kk
    real(dp) potgrad((kk+1)*(kk+2)/2)
    integer i!,j
    integer ki1, ki2, ni1, ni2, kni, kni2
    
    potgrad = 0
    do i = 0, min(nq,kk) !0!
        
        kni = kk+nq-i
        kni2 = kni*2
        
        ki1 = gpos(kk-i)+1
        ki2 = gpos(kk-i+1)
        ni1 = gpos(nq-i)+1
        ni2 = gpos(nq-i+1)
        
        potgrad(:) = potgrad(:) &
                    + (-1)**kni * intfac(nq,nq-i) * intff(kni2-1,2*nq-1) * rpows(kni2+1) &
                    * symouter(kk-i, i, rrr(ki1:ki2),   inner(nq,nq-i, qq, rrr(ni1:ni2))   )
    enddo
    
end
  !integer, parameter :: rpos(8) = [1,      4,     10,     20,     35,     56,     84,    120] ! remove this, use gpos
    
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

subroutine test_matri(nn)
    integer i, nn, maxi
    if(nn>7)stop"rank cant be larger than 7 in matr intdex matrix"
    maxi = sumfac(nn)
    do i = 1, maxi
      print'(28I3)', matri(i,1:maxi) - matr(i,1:maxi)
    enddo
end subroutine





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
    p1 = gpos(i-1)+1
    p2 = gpos(i)
    p3 = gpos(i)+1
    p4 = gpos(i+1)
    
    rrr2=rrr(p1:p2)
    rrr3=rrr(p3:p4)
    
    i=5
    p3 = gpos(i)+1
    p4 = gpos(i+1)
    
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





end module
