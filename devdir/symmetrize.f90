module symmetrize

use data_types, only: dp
use printer_mod, only: str, printer, printo
implicit none

contains !///////////////////////////////////////////////

subroutine main
    call testing
end subroutine

! Testing /////////////////////////////////////////////////////
subroutine testing
    
    call test_sorted
    
    call test_next_key2n(3)
    
    !call test_tric_prod()
    
    call test__expand_compress
    
end subroutine

subroutine test__expand_compress()
    integer rank
    real(dp) :: tricorn(10), full(3,3,3)!, linfu(3**3)
    
    rank = 3
    tricorn = [1d0, 2d0, 3d0, 4d0, 5d0, 6d0, 7d0, 8d0, 9d0, 10d0]
    full = reshape(expand(tricorn,rank),shape(full))
    call printo(full,[2,1,3])
    
    
    print '('//str( triclen(rank) )//'f7.3)', compress(reshape(full,[27]),rank)
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

subroutine test_next_key2n(rank) !result(ns)
    integer, intent(in) :: rank
    integer trilen, key(rank), nn(3)
    !integer :: ns(3, ((rank+1)*(rank+2))/2)
    integer i
    print '('//str(4)//'I2)', next([1,1,3,3],0)
    print '('//str(4)//'I2)', next([1,1,3,3],1)
    
    print '('//str( ((rank+1)*(rank+2))/2 )//'I6)', key2n([1,1,2,3])
    
    trilen = ((rank+1)*(rank+2))/2 ! length of tricorn vector given full tensor rank
    
    key = 1
       
    print*, 'n(1:3) array:'
    nn = key2n(key)
    print*, nn, "row:"//str(1), "  found row:"//str(finder(nn)), "  g:"//str(apple_g(nn))
    do i = 2,trilen
       key = next(key,1)
       
       nn = key2n(key)
       print*, nn, "row:"//str(i), "  found row:"//str(finder(nn)), "  g:"//str(apple_g(nn))
    enddo
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
    real(dp) :: linfull(3**rank), tricorn( ((rank+1)*(rank+2))/2 ) 
    integer i, j, ifull, trilen
    trilen = ((rank+1)*(rank+2))/2
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
    real(dp) tricorn( (rank+1)*(rank+2)/2 ) , linfull(3**rank)
    integer i, tri_ind, nn(3)
    
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
   integer i, n(3), rank
   rank = size(key)
   n=0
   ! get number of x,y,z indices (or rather 1,2,3 indices)
   do i = 1,rank
     if(key(i) == 1)then
       n(1)=n(1)+1
     elseif(key(i) == 2)then
       n(2)=n(2)+1
     elseif(key(i) == 3)then
       n(3)=n(3)+1
     else
       stop"error in key2n, index not 1, 2 or 3"
     endif
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

function next(key,tric)
   ! Function returning the next combination of indices (key) given the previous key. 
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
