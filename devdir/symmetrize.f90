module symmetrize

use data_types, only: dp
use printer_mod, only: str
implicit none
contains

subroutine main
integer, allocatable :: key(:)!, key2(:)
integer rank

key = [4,5,2,1,3,1,7,1,3,2,6,1]
rank = size(key)

print '('//str(rank)//'I2)', sorted(key)

!key2 = next(key)

print '('//str(rank)//'I2)', next([1,1,3,3],0)
print '('//str(rank)//'I2)', next([1,1,3,3],1)

!rank=1; print '('//str( ((rank+1)*(rank+2))/2 )//'I6)', tric_prods(rank)
!rank=2; print '('//str( ((rank+1)*(rank+2))/2 )//'I6)', tric_prods(rank)
!rank=3; print '('//str( ((rank+1)*(rank+2))/2 )//'I6)', tric_prods(rank)
!rank=4; print '('//str( ((rank+1)*(rank+2))/2 )//'I6)', tric_prods(rank)
!rank=5; print '('//str( ((rank+1)*(rank+2))/2 )//'I6)', tric_prods(rank)
!rank=6; print '('//str( ((rank+1)*(rank+2))/2 )//'I6)', tric_prods(rank)
!rank=7; print '('//str( ((rank+1)*(rank+2))/2 )//'I6)', tric_prods(rank)
!rank=8; print '('//str( ((rank+1)*(rank+2))/2 )//'I6)', tric_prods(rank)
!rank=9; print '('//str( ((rank+1)*(rank+2))/2 )//'I6)', tric_prods(rank)
!print*, ''
!rank=1; print '('//str( ((rank+1)*(rank+2))/2 )//'I6)', sorted(tric_prods(rank))
!rank=2; print '('//str( ((rank+1)*(rank+2))/2 )//'I6)', sorted(tric_prods(rank))
!rank=3; print '('//str( ((rank+1)*(rank+2))/2 )//'I6)', sorted(tric_prods(rank))
!rank=4; print '('//str( ((rank+1)*(rank+2))/2 )//'I6)', sorted(tric_prods(rank))
!rank=5; print '('//str( ((rank+1)*(rank+2))/2 )//'I6)', sorted(tric_prods(rank))
!rank=6; print '('//str( ((rank+1)*(rank+2))/2 )//'I6)', sorted(tric_prods(rank))
!rank=7; print '('//str( ((rank+1)*(rank+2))/2 )//'I6)', sorted(tric_prods(rank))
!rank=8; print '('//str( ((rank+1)*(rank+2))/2 )//'I6)', sorted(tric_prods(rank))
!rank=9; print '('//str( ((rank+1)*(rank+2))/2 )//'I6)', sorted(tric_prods(rank))

rank=4; print '('//str( ((rank+1)*(rank+2))/2 )//'I6)', key2n([1,1,2,3])

call tricorn_ns(3)
end subroutine

!subroutine symmetrize6(tricorn)
!real(dp) :: tricorn(6)
!end subroutine




subroutine tricorn_ns(rank) !result(ns)
   integer rank, trilen, key(rank), nn(3)
   !integer :: ns(3, ((rank+1)*(rank+2))/2)
   integer i
   
   trilen = ((rank+1)*(rank+2))/2 ! length of tricorn vector given full tensor rank
   
   key = 1
      
   print*, 'n(1:3) array:'
   nn = key2n(key)
   print*, nn, "row:",1, "found row:",finder(nn), "g:", apple_g(nn)
   do i = 2,trilen
      key = next(key,1)
      
      nn = key2n(key)
      print*, nn, "row:",i, "found row:",finder(nn), "g:", apple_g(nn)
   enddo
end

function apple_g(n) result(g)
   ! Apple quist g()-function. 
   ! Is the number of times unique element occurs in a full symmetric tensor
   integer g, n(3), rank, f,f1,f2,f3,ns(3), i
   
   rank = sum(n)
   ns = sorted(n)
   f=1
   f1=1
   f2=1
   f3=1
   do i = 1,ns(1)
     f1 = f1*i
     enddo
   
   do i = 1,ns(2)
     f2 = f2*i
     enddo
   
   do i = 1,ns(3)
     f3 = f3*i
     enddo
   
   do i = 1,rank
     f = f*i
     enddo
   
   g=f/(f1*f2*f3)
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

end module
































!/////////////////////////////////////////////////

#ifdef OUT
subroutine expand(tricorn,fullt,rank)
integer rank, trilen, key(rank)
real(dp) :: tricorn( ((rank+1)*(rank+2))/2 ), linfull(3**rank)
integer :: prods( ((rank+1)*(rank+2))/2 ), prod
integer i, j, posit, key(rank)

trilen = ((rank+1)*(rank+2))/2 ! length of tricorn vector given full tensor rank

key = 1
prods(1) = 1

do i = 2,trilen
   
   key = next(key)
   
   prod=1
   do j = 1,rank
     prod = prod*key(j)
   enddo
   
   prods(i) =  prod
enddo

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
#endif
