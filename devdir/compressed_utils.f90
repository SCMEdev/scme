module compressed_utils

use printer_mod, only: str, printer, printo
implicit none

integer, parameter :: dp = kind(0d0)
private dp
!public test_apple_g

contains !//////////////////////////////////////////////////
subroutine main
integer i
    !call testing
    !call test_sumfac
    !call test_rpow
    !call test_nextpow(7)
    !call compow2(7)
    !call subdiv_pow(4,4)
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
    call make_gs
    print'(*(a))', (str(sumfacfac(i))//", ", i = 1,10)
    !call test_rrpow
    !call test_next_rev_key2n(4)
end subroutine


subroutine make_gs
integer i, j, n(3), it
it = 0
do i = 1,10
  n = 0
  n(1)= i
  do j = 1, sumfac(i+1)
    !if(j>1)call nextpov(n)
    if(j>1)call nextpow(n(1),n(2),n(3))
    
    !print'(4(I4))', n, apple_g(n)
    write(*,'((I5,a))',advance="no") apple_g(n),','
    it = it+1
    enddo
  print*,'&'
  enddo
print*, "entries", it
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

subroutine test_nextpo_w_v

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

subroutine nextpov(nn)
  integer a,b,c, nn(3)
  a = nn(1)
  b = nn(2)
  c = nn(3)

  if (b>0)then
    b=b-1
    c=c+1
  elseif (a>0)then
    a=a-1
    b=c+1
    c=0
  endif
  
  nn(1) = a
  nn(2) = b
  nn(3) = c

  
end

   
subroutine test_sumfac
   integer i
   
   do i = 1,10
     print*, sumfac(i) , i, 1
   enddo
endsubroutine


pure function sumfac(u)
    !third row in pascal matrix (for tricorn lengths)
    integer sumfac
    integer, intent(in) :: u
    sumfac = u*(u+1)/2 ! 2=2!
endfunction



subroutine test_choose
print*,'choose(6,2)', choose(6,2)
print*,'choose(5,2)', choose(5,2)
print*,'choose(6,3)', choose(6,3)
print*,'choose(2,4)', choose(4,2)

end
!function choose(M,k)
!    implicit none
!    integer, intent(in) :: M, k
!    integer :: choose
!    integer :: p,  i
!    p=min(k,M-k)
!    choose=1
!    do i=1,p
!      choose=(choose*(M-i+1))/i 
!    enddo
!    
!
!end 
function choose(m1,m2)
    implicit none
    integer, intent(in) :: m1, m2
    integer :: choose, MM, kk
    integer :: imax,  i
    MM = max(m1,m2)
    kk = min(m1,m2)
    imax = min(kk,MM-kk)
    choose=1
    do i=1,imax
      choose=(choose*(MM-i+1))/i 
    enddo
    

end 


subroutine test_apple_g
print*, "xxyy", apple_g([2,2,0])
print*, "xy", apple_g([1,1,0])

print*, "xxyz", apple_g([2,1,1])
print*, apple_g([1,1,0])
print*, apple_g([1,0,1])
print*, apple_g([0,1,1])

print*, "xxxy", apple_g([3,1,0])
end subroutine

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



subroutine test_sorted()
    integer, allocatable :: key(:)!, key2(:)
    integer rank
    key = [4,5,2,1,3,1,7,1,3,2,6,1]
    rank = size(key)
    print '('//str(rank)//'I2)', sorted(key)
end subroutine


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

function next(key,tric) 
  integer key(:), tric,next(size(key))
  character(3) :: order
  order="can"
  if(order=="can") next = next_can(key,tric)
  if(order=="lex") next = next_lex(key,tric)
endfunction  

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

pure function sumfacfac(u)
    !fourth row in pascal matrix ( for tricorn positions in polytensor)
    integer sumfacfac
    integer, intent(in) :: u
    sumfacfac = u*(u+1)*(u+2)/6 ! 6=3!
endfunction




endmodule
