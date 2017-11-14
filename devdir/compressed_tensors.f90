module compressed_tensors

use data_types, only: dp
use printer_mod, only: str, printer, printo
implicit none

integer j11,j22
integer, parameter :: k = 7
integer, parameter :: vv((k+1)*(k+2)/2) = [((j11, j22=1,j11),j11 = 1,k+1)]
integer, parameter :: vv1((k+1)*(k+2)/2) = vv-1
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


contains !///////////////////////////////////////////////




subroutine main
    !call testing
    !call test_sumfac
    !call test_rpow
    !call test_nextpow(7)
    !call compow2(7)
    call subdiv_pow(3,2)
    print*,"fac 5 / 2 3", fac(5)/( fac(3)*fac(2))
    !call test_factorial(5)
    
    
    !print*, ""
    !call test_contind(6)
    
    !call test_matr(7)
    !call printer(matr, 'matr',1)
    !call test_apple_g
    !call test_rrpow
    !call test_next_rev_key2n(4)
end subroutine

! Testing /////////////////////////////////////////////////////

subroutine test_matr(nn)
    integer i, nn, maxi
    if(nn>7)stop"rank kant be larger than 7 in matr intdex matrix"
    maxi = sumfac(nn)
    do i = 1, maxi
      print'(28I3)', matr(i,1:maxi)
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

subroutine test_nextpow(k)
  integer a,b,c, k
  integer ind
  
  a=k
  b=0
  c=0
  
  do ind = 1,sumfac(k+1)
     if(ind>1)call nextpow(a,b,c)
     print'(3I2,I5)',a,b,c,ind
     enddo
     
end

subroutine nextpow(a,b,c)
  integer a,b,c

  if (b>0)then
    b=b-1
    c=c+1
  elseif (a>0)then
    a=a-1
    b=c+1
    c=0
  endif
end

subroutine subdiv_pow(ii,kii)
integer ii, kii
integer i, j, it
integer a1,b1,c1, a2,b2,c2, aa,bb,cc

a1 = ii
aa = ii + kii

b1 = 0
c1 = 0
bb = 0
cc = 0

it = 0
do i = 1, sumfac(ii+1)
   if (i>1)call nextpow(a1,b1,c1)
   b2 = 0
   c2 = 0
   a2 = kii
   do j = 1, sumfac(kii+1)
      if (j>1)call nextpow(a2,b2,c2)
      it = it+1
      aa = a1+a2
      bb = b1+b2
      cc = c1+c2
      print '(3(I3,2I2),7I4)', a1,b1,c1, a2,b2,c2, aa,bb,cc, &
            finder([aa,bb,cc]),&
            matr(i,j), &
            i + j + (vv(i)-1)*(vv(j)-1) - 1, &
            i + j + vv1(i)*vv1(j) - 1, &
            hh(a1,b1,c1, a2,b2,c2), &
            i, j !, it
      
      enddo;enddo
    
   !do here!
print'(*(I3))', vv(1:sumfac(ii+kii+1))
print'(*(I3))', vv1(1:sumfac(ii+kii+1))
print*, sumfac(ii+kii+1)
end

function hh(a1,b1,c1, a2,b2,c2) !Make a matrix of this sheeeet
    integer a1,b1,c1, a2,b2,c2
    integer aa,bb,cc, hh
    aa = a1+a2
    bb = b1+b2
    cc = c1+c2
    
    hh = fac(aa)*fac(bb)*fac(cc) / ( fac(a1)*fac(a2)*fac(b1)*fac(b2)*fac(c1)*fac(c2) )

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

subroutine test_apple_g
print*, "xxyy", apple_g([2,2,0])
print*, "xy", apple_g([1,1,0])

print*, "xxyz", apple_g([2,1,1])
print*, apple_g([1,1,0])
print*, apple_g([1,0,1])
print*, apple_g([0,1,1])

print*, "xxxy", apple_g([3,1,0])
end subroutine

function next(key,tric) 
  integer key(:), tric,next(size(key))
  character(3) :: order
  order="can"
  if(order=="can") next = next_can(key,tric)
  if(order=="lex") next = next_lex(key,tric)
endfunction  

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

subroutine test_sumfac
   integer i
   
   do i = 1,10
     print*, sumfac(i) , i, 1
   enddo
endsubroutine

subroutine testing_many
    
    !call test_sorted
    
    call test_next_key2n(5,1)
    
    !call test_tric_prod()
    
    !call test__expand_compress
    
end subroutine

subroutine test__expand_compress()
    integer rank
    real(dp) :: tricorn(10), full(3,3,3)!, linfu(3**3)
    real(dp) :: tricorn4(15), full4(3,3,3,3)!, linfull4(3**4)!, linfu(3**3)
    
    rank = 3
    tricorn = [1d0, 2d0, 3d0, 4d0, 5d0, 6d0, 7d0, 8d0, 9d0, 10d0]
    full = reshape(expand(tricorn,rank),shape(full),order=[3,1,2])
    call printo(full,[2,1,3])
    
    
    
    print '('//str( triclen(rank) )//'f7.3)', compress(reshape(full,[3**rank]),rank)
    
    rank=4
    tricorn4 = [1d0, 2d0, 3d0, 4d0, 5d0, 6d0, 7d0, 8d0, 9d0, 10d0,11d0,12d0,13d0,14d0,15d0]
    full4 = reshape(expand(tricorn4,rank),shape(full4),order=[3,1,2,4])
    call printo(full4,[1,2,3,4])
    !linfull4 = reshape(full,[3**3])
    print '('//str( triclen(rank) )//'f7.3)', compress(reshape(full4,[3**4]),rank)
    
    
end subroutine


subroutine test_sorted()
    integer, allocatable :: key(:)!, key2(:)
    integer rank
    key = [4,5,2,1,3,1,7,1,3,2,6,1]
    rank = size(key)
    print '('//str(rank)//'I2)', sorted(key)
end subroutine


subroutine test_tric_prod()
    integer rank
    rank=1; print '('//str( ((rank+1)*(rank+2))/2 )//'I6)', tric_prods(rank)
    rank=2; print '('//str( ((rank+1)*(rank+2))/2 )//'I6)', tric_prods(rank)
    rank=3; print '('//str( ((rank+1)*(rank+2))/2 )//'I6)', tric_prods(rank)
    rank=4; print '('//str( ((rank+1)*(rank+2))/2 )//'I6)', tric_prods(rank)
    rank=5; print '('//str( ((rank+1)*(rank+2))/2 )//'I6)', tric_prods(rank)
    rank=6; print '('//str( ((rank+1)*(rank+2))/2 )//'I6)', tric_prods(rank)
    rank=7; print '('//str( ((rank+1)*(rank+2))/2 )//'I6)', tric_prods(rank)
    rank=8; print '('//str( ((rank+1)*(rank+2))/2 )//'I6)', tric_prods(rank)
    rank=9; print '('//str( ((rank+1)*(rank+2))/2 )//'I6)', tric_prods(rank)
    print*, ''
    rank=1; print '('//str( ((rank+1)*(rank+2))/2 )//'I6)', sorted(tric_prods(rank))
    rank=2; print '('//str( ((rank+1)*(rank+2))/2 )//'I6)', sorted(tric_prods(rank))
    rank=3; print '('//str( ((rank+1)*(rank+2))/2 )//'I6)', sorted(tric_prods(rank))
    rank=4; print '('//str( ((rank+1)*(rank+2))/2 )//'I6)', sorted(tric_prods(rank))
    rank=5; print '('//str( ((rank+1)*(rank+2))/2 )//'I6)', sorted(tric_prods(rank))
    rank=6; print '('//str( ((rank+1)*(rank+2))/2 )//'I6)', sorted(tric_prods(rank))
    rank=7; print '('//str( ((rank+1)*(rank+2))/2 )//'I6)', sorted(tric_prods(rank))
    rank=8; print '('//str( ((rank+1)*(rank+2))/2 )//'I6)', sorted(tric_prods(rank))
    rank=9; print '('//str( ((rank+1)*(rank+2))/2 )//'I6)', sorted(tric_prods(rank))

end subroutine

subroutine test_next_key2n(rank,tric) !result(ns)
    integer, intent(in) :: rank, tric
    integer trilen, key(rank), nn(3), upper
    !integer :: ns(3, ((rank+1)*(rank+2))/2)
    integer i, j
    print '('//str(4)//'I2)', next_can([3,3,1,1],0)
    print '('//str(4)//'I2)', next_can([3,3,1,1],1)
    
    print '('//str(4)//'I2)', next_can([1,1,3,3],0)
    print '('//str(4)//'I2)', next_can([1,1,3,3],1)
    
    print '('//str( ((rank+1)*(rank+2))/2 )//'I6)', key2n([1,1,2,3])
    
    trilen = ((rank+1)*(rank+2))/2 ! length of tricorn vector given full tensor rank
    
    key = 1
       
    print*, 'n(1:'//str(rank)//') array:'
    nn = key2n(key)
    i = 1
    call pprint
    
    if(tric==1)upper=trilen
    if(tric==0)upper=3**rank
    
    do i = 2,upper!3**rank!trilen
       key = next(key,tric)
       
       nn = key2n(key)
       call pprint
    enddo
    
    
    
    contains !// ////////////////////////
      
      
      
      subroutine pprint
         print*, "key=",(str(key(j))//" ",j=1,rank ),"  nn=",(str(nn(j))//" ",j=1,3 ),&
                 "  row:"//str(i), "  found row:"//str(finder(nn)), "  g:"//str(apple_g(nn))
      end subroutine
      
end



! END TESTING //////////////////////////////////////////////////////////////




!///////////////////////////////////////////////////////////////////////////////////////
pure function triclen(rank)
    integer, intent(in) :: rank
    integer triclen
    triclen = ((rank+1)*(rank+2))/2
end

function compress(linfull,rank) result(tricorn)
    integer rank, key(rank)
    real(dp) :: linfull(:), tricorn( ((rank+1)*(rank+2))/2 )! 3**rank
    integer i, j, ifull, trilen
    trilen = ((rank+1)*(rank+2))/2
    !if(size(tricorn) /= trilen)stop"Wrong tricorn tensor has wrong length in compress()"
    if(size(linfull) /= 3**rank)stop":: length(linear full tensor) /= 3**rank ::"
    
    key = 1
    tricorn(1) = linfull(1)
    do i = 2,trilen
      key = next(key,1)
      
      ifull = 1
      do j = 1,rank
        ifull = ifull + (key(j)-1)*3**(j-1) 
      enddo
      tricorn(i) = linfull(ifull)
    enddo
end 

function expand(tricorn, rank) result(linfull)
    integer rank, key(rank)
    real(dp) tricorn( : ) , linfull(3**rank)
    integer i, tri_ind, nn(3), trilen
    trilen = (rank+1)*(rank+2)/2
    if(size(tricorn) /= trilen)stop"length(tricorn tensor) /= (rank+1)*(rank+2)/2"
    key = 1
    linfull = 0
    linfull(1) = tricorn(1)
    
    do i = 2,3**rank 
      key = next(key,0)
      nn = key2n(key)
      tri_ind = finder(nn)
      linfull(i) = tricorn(tri_ind)
    enddo
    
end 

function apple_g(n) result(g)
   ! Applequist g()-function. 
   ! Is the number of times unique element occurs in a full symmetric tensor = numer of intex permutations giving the same value
   ! this routine will fail at rank 23
   integer g, n(3), rank, ns(3)
   integer*8 j, nfac, fac1, fac2, gl!using long ints for the factorials to be safe
   ns = sorted(n)
   rank = sum(ns)
   nfac = 1
   fac1 = 1
   fac2 = 1
   
   do j = ns(3)+1, rank
     nfac = nfac*j
     enddo
   
   do j = 1, ns(1)
     fac1 = fac1*j
     enddo
   
   do j = 1, ns(2)
     fac2 = fac2*j
     enddo
   
   gl = nfac/(fac1*fac2)
   
   if(gl>huge(g))stop"Wow, too big value in g()"
   
   g=int(gl,kind(g))
end
   
function key2n(key) result(n)
   integer, intent(in) :: key(:)
   integer i, n(3), rank, ind
   rank = size(key)
   n=0
   ! get number of x,y,z indices (or rather 1,2,3 indices)
   do i = 1,rank
     ind = key(i)
     n(ind) = n(ind) + 1
   enddo
end function

function finder(n) result(row)
   ! Given nx,ny,nz, returns corresponding row in tricorn vector. 
   integer n(3), row, rank, i
   
   rank = sum(n)
   row = 1 + n(3)
   do i = 1,rank - n(1) !addition factorial
     row = row + i
   enddo
end


pure function sumfac(u)
    !third row in pascal matrix (for tricorn lengths)
    integer sumfac
    integer, intent(in) :: u
    sumfac = u*(u+1)/2 ! 2=2!
endfunction

pure function sumfacfac(u)
    !fourth row in pascal matrix ( for tricorn positions in polytensor)
    integer sumfacfac
    integer, intent(in) :: u
    sumfacfac = u*(u+1)*(u+2)/6 ! 6=3!
endfunction



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






function next_lex(key,tric) result(next)
   ! Function returning the next combination of indices (key) given the previous key, IN LEXICOGRAPHICAL ORDER. 
   ! If tric=1 it gives the next tricorn key, if tric=0 it gives the next full-tensor key. 
   integer, intent(in) :: key(:)
   integer, intent(in), optional :: tric
   integer :: next(size(key))
   integer rank, i, j, vali
   
   rank = size(key)
   next = key
   
   outer:do i = rank, 1, -1
     vali=next(i)
     if(vali < 3)then
       next(i) = vali+1
       do j = i+1,rank
         next(j) = 1 + tric*vali !if doing tricorn increment then add this factor. tric is 0 or one
       enddo
       exit outer
     endif
   enddo outer
end function


function next_can(key,tric) result(next)
   ! Function returning the next combination of indices (key) given the previous key. IN CANONCIAL = ANTI-LEXICOGRAPHICAL = FORTRAN-LIKE, ORDERING. 
   ! If tric=1 it gives the next tricorn key, if tric=0 it gives the next full-tensor key. 
   integer, intent(in) :: key(:)
   integer, intent(in), optional :: tric
   integer :: next(size(key))
   integer rank, i, j, vali
   
   rank = size(key)
   next = key
   
   outer:do i = 1,rank!, 1, -1
     vali=next(i)
     if(vali < 3)then
       next(i) = vali+1
       do j = i-1,1, -1
         next(j) = 1 + tric*vali !if doing tricorn increment then add this factor. tric is 0 or one
       enddo
       exit outer
     endif
   enddo outer
end function



pure function sorted(key) 
   ! given an unsorted key, returns surted key. 
   integer, intent(in) :: key(:)
   integer :: sorted(size(key))
   integer :: saveit, minimum, i, j, rank, minposi
   sorted = key
   rank = size(sorted)
   
   do i = 1,rank-1
      saveit = sorted(i)
      minimum = saveit
      minposi = i 
      do j = i+1, rank
        if(sorted(j).le.minimum)then
          minimum = sorted(j)
          minposi = j
        endif
      enddo
      sorted(i)=minimum
      sorted(minposi)=saveit
   enddo
end function



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


function tric_prods(rank) result(prods)
   integer rank, trilen, key(rank)
   integer :: prods( ((rank+1)*(rank+2))/2 ), prod
   integer i, j
   
   trilen = ((rank+1)*(rank+2))/2 ! length of tricorn vector given full tensor rank
   key = 1
   prods(1) = 1
   
   do i = 2,trilen
   !print*, 2
      
      key = next(key,1)
      
      prod=1
      do j = 1,rank
        prod = prod*key(j)
      enddo
      
      prods(i) =  prod
   enddo
end function

end module
