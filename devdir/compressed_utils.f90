module compressed_utils

use printer_mod, only: str, printer, printo
implicit none

integer, parameter :: dp = kind(0d0)
private dp
!public test_apple_g

contains !//////////////////////////////////////////////////

function del(a,b)
    integer a, b, del
    del=0
    if (a==b) del=1
end


subroutine print_trace_index_matrix(k)
    integer k
    integer j1,j2, icol, irow, ileng,v2(sumfac(k+1)), v4(sumfac(k+1)),v3(sumfac(k+1)), tleng !((k+1)*(k+2)/2),

    ileng=sumfac(k+1)
    tleng = sumfac((k/2+1))
    v2 = [((j1, j2=1,j1),j1 = 1,k+1)]-1 !for index matrix
    v3 = 2*v2![((j1, j2=1,j1,2),j1 = 1,k+1,2)] !for trace matrix
    v4 = v2*v2 !for trace matrix
    print*,"______"
    print'(a,*(I3))',"v2: ", v2
    print'(a,*(I3))',"v3: ", v3
    print'(a,*(I3))',"v4: ", v4
    
    print*
    do irow = 1, ileng
        write(*,'(*(I3,a),a)'),((icol-1)*2 + v4(icol) + irow + (v3(icol))*(v2(irow)),',', icol = 1,tleng)
    enddo
end subroutine

subroutine print_product_index_matrix(k)
    !integer, parameter :: k = 7
    integer k
    integer j1,j2, icol, irow, ileng,v2(sumfac(k+1))
    
    ileng=sumfac(k+1)
    
    v2 = [((j1, j2=1,j1),j1 = 1,k+1)]-1 !for index matrix
    
    print*,"______"
    print'(a,*(I3))',"v2: ", v2
    
    
    print*
    do irow = 1, ileng
        write(*,'(*(I2,a),a)'),(irow + icol -1 + (v2(icol))*(v2(irow)),',', icol = 1,ileng)
    enddo
    
    !Print tight
    !print*, sumfac(k+1), "cols/rows"
    !do irow = 1, sumfac(k+1)
    !    !v3 = [(irow + icol -1 + (v2(icol)-1)*(v2(irow)-1), icol = 1,ileng)]
    !    !print'(*(I4))',v3
    !    !write(*,'(*(a))',advance="no"),(str(irow + icol -1 + (v2(icol)-1)*(v2(irow)-1))//',', icol = 1,ileng)
    !    write(*,'(*(a))'),(str(irow + icol -1 + (v2(icol)-1)*(v2(irow)-1))//',', icol = 1,ileng)
    !    !print'(*(I4))',v3
    !    
    !    enddo
    
end subroutine

subroutine print_choose_matrix
    integer, parameter :: k = 19
    integer i,j
    !integer , parameter :: k = 10
    !integer, parameter :: init(0:k,0:k) = 
    integer mat(0:k,0:k), temp
    
    ! Put 1 on edges and diagonal
    mat(:,0)=1
    mat(0,:)=1
    do i = 0,k
        mat(i,i) = 1
    enddo
    
    ! Fill lower triangle with sum rule, symmetrize
    do i = 2, k
        do j = 1,i-1
            temp = mat(i-1,j-1) + mat(i-1, j)
            mat(i,j) = temp
            mat(j,i) = temp
        enddo
    enddo
    
    ! Print for copying
    write(*,'(a)') "integer, parameter :: choose_matrix(0:"//str(k)//",0:"//str(k)//") = reshape( [&"
    call printo(mat,0,1)
    write(*,'(a)') " ], shape(choose_matrix) )"
end subroutine


subroutine print_apple_g
    integer i, j, n(3), it, gsiz, linelen, gval, glen
    character(1) zero
    print*,"____"
    print*,">>>>  g-arrays:"
    it = 0
    zero = "0"
    do i = 0,10
      glen = sumfac(i+1)
      if(i>9)zero=""
      write(*,'(a)',advance="no") "integer, parameter :: g"//trim(zero)//str(i)//"("//str(glen)//") = ["
      if (i>6) write(*,'(a)') "& "
      
      linelen = 0
      
      n = 0
      n(1)= i
      do j = 1, glen
        if(j>1)call nextpown(n)
        
        gval = apple_g(n)
        gsiz = ceiling(log10(dble(gval+1))) 
        linelen = linelen + gsiz+1
        write(*,'(I'//str(gsiz)//')',advance="no") gval
        it = it+1
        if (j<glen)then
          write(*,'(a)',advance="no") ','
          if(linelen>110)then
            print'(a)', '& '
            linelen = 0
          endif
        else 
          write(*,'(a)') "]"
        endif
      enddo
      !print*!,'&'
      enddo
    print*, "entries", it
end subroutine



subroutine main
!integer i
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
    !call make_gs
    !print'(*(a))', (str(sumfacfac(i))//", ", i = 1,10)
    !call test_rrpow
    !call test_next_rev_key2n(4)
    
    !call test_next_key2n(7,1)
    
    !call print_choose_matrix
    
    !call test_nextpown
    !call print_apple_g
    
    !call print_product_index_matrix(6)
    !call print_trace_index_matrix(6)
    !call print_traces(6)
    
    !call test_choose
    
    
    !call test_brakk
    call main2
    
    
end subroutine

subroutine test_nr_of_trace_elements
    integer k, n, su
    do n = 12,1,-1
      su = 0
      do k= 1,n/2
        su = su + sumfac(n-2*k+1) 
      enddo
      print'(a)',"rank="//str(n)//" self="//str(sumfac(n+1))//" sum="//str(su)
    enddo
end subroutine



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
integer, parameter :: k = 10
integer i,j, mat(0:k,0:k)
print*,'choose(6,2)', choose(6,2)
print*,'choose(5,2)', choose(5,2)
print*,'choose(6,3)', choose(6,3)
print*,'choose(2,4)', choose(4,2)
print*,'choose(4,0)', choose(4,0)
print*,'choose(0,4)', choose(0,4)
print*,'choose(0,2)', choose(0,2)

do i = 0, k
  do j = 0, k
    mat(i,j) = choose(i,j)
  enddo
enddo
call printo(mat,3)
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
pure function choose(m1,m2)
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




subroutine print_traces(rank)
    integer :: mbar(3), kbar(3), sup(3), rank, i, j, k
    
    
    do k = 1, rank/2
        print*
        print*, "<<"//str(k)//">>-fold trace:"
        
        
        do i = 1,min(3000,sumfac(rank-2*k+1))
            !print'(I2,a)',i,"th entry of trace array:"
        
        
        
            
            
            if(i==1)then
                mbar = [rank-2*k,0,0]
            else
                call nextpown(mbar)
            endif
                
            !print*," mbar =["//str(mbar,2)//"], finder(mbar)="//str(finder(mbar))//", i="//str(i)
            
            
            do j = 1, sumfac(k+1)
                if(j==1)then
                    kbar = [k, 0, 0]
                else
                    call nextpown(kbar)
                endif 
                sup = kbar*2 + mbar
                !print*, "new:"
                !print'(*(a,3I4))','  kbar:',kbar, ' 2key:', 2*kbar, ' mbar+2key:', sup
                !print'(*(a,I4))','  kbar:',finder(kbar), ' 2key:',finder(2*kbar), ' mbar+2key:',finder(sup)
                !print'(a)', "trace: finder(["//str(sup,2)//"])="//str(finder(sup),3)
                write(*,'(I3)',advance="no") finder(sup)
            enddo
            print*!, "<< finished trace"
        enddo
            
            
        
    enddo
end subroutine
    


function finder(n) result(row)
   ! Given nx,ny,nz, returns corresponding row in tricorn vector. 
   integer n(3), row, n23 !rank, i
   
   n23 = n(2)+n(3)
   row = 1 + n(3) + n23*(n23+1)/2
   
   !do i = 1,rank - n(1) !addition factorial
   !  row = row + i
   !enddo
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


subroutine test_intfac_ff
    integer a, b
    character(3) ach, bch
    call get_command_argument(1, ach)
    call get_command_argument(2, bch)
    
    read(ach,*) a
    read(bch,*) b
    
    print*, intfac(a,b)
    print*, intff(a,b)
    !print*, intfac(5,5)
    !print*, intfac(5,4)
    !print*, intfac(5,3)
    !print*, intfac(7,3)
    !print*, intfac(7,5)
    !
    !print*, intff(5,5)
    !print*, intff(6,4)
    !print*, intff(5,3)
    !print*, intff(7,3)
    !print*, intff(7,5)
    !print*, intff(17,7)
    
    print*, "above, TEST_INTFAC_FF ------------------------------------------------------------------"
    
end
   

function intfac(aa,bb) result(cc)
    ! Factorial in range: intfac(A,B) = A!/B!
    integer aa, bb
    integer i, cc
    if( aa<0 .or. bb<0 )stop"negative number in intfac"
    cc = 1
    do i = bb+1,aa
      cc = cc*i
    enddo
end

function intff(aa,bb) result(cc)
    ! Double factorial in range: intff(A,B) = A!!/B!! if both A and B are odd (or even)
    integer aa, bb
    integer i, cc
    if (mod(aa-bb,2) .ne. 0)stop"number interval not divisible by 2 in intff"
    if( aa<0 .or. bb<0 )stop"negative number in intff"
    cc = 1
    do i = bb+2,aa, 2
      cc = cc*i
    enddo
end

subroutine test_fac(k)
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

function lfac(inn) result(fac)
    integer*8 inn, i, fac
    fac=1
    do i = 1, inn
      fac = fac*i
      enddo

end


subroutine test_brakk
    integer n, l
    do n = 1, 20
    do l = 1, n/2
    print*, n,l, brakk(n,l), brakk_bad(n,l)
    enddo
    enddo
end subroutine

function vecbrakk(nn,ll) result(res)
    integer nn(3), ll(3)
    integer res
    res = brakk(nn(1),ll(1)) * brakk(nn(2),ll(2)) * brakk(nn(3),ll(3)) 
end function


function brakk(n,l) result(p)
    ! This routine calculates the expression n!/(2^l*l!*(n-2l)),
    ! i.e. how many ways to select index pairs (for kroneker-delta) 
    ! from a group of indices without regards to order within the pair or group.
    ! It is done through a numerically stable recursion relation that I derived.
    integer,intent(in) :: n, l
    integer i, l2
    integer p
    
    l2 = 2*l
    
    if(n<l2)stop"n<l2 in brakk"
    if (n>19 .and. l>5)stop"Too big value in 'brakk()'"
    
    p=1
    do i = l2+1, n   !recursion relation
      p = (p*i)/(i-l2)
    enddo
    do i = 3,l2-1, 2 !remaining double factorial
      p = p*i
    enddo
    
end function



subroutine test_nextpown!_apple
  integer n(3), k
  integer ind
  
  k = 4
  n(1) = k
  n(2) = 0
  n(3) = 0
  
  do ind = 1,sumfac(k+1)
     if(ind>1)call nextpown(n)
     print'(3I2,2I5)',n,ind, apple_g(n)
     enddo
     
end

subroutine nextpown(nn)
  integer, intent(inout) :: nn(3) 
  integer a,b,c 
  a=nn(1);b=nn(2);c=nn(3)
  if (b>0)then
    b=b-1
    c=c+1
  elseif (a>0)then
    a=a-1
    b=c+1
    c=0
  endif
  nn = [a,b,c]
end

!subroutine nextpown(n)
!  integer, intent(inout) :: n(3) 
!  integer a,b,c 
!  if (n(2)>0)then
!    n(2)=n(2)-1
!    n(3)=n(3)+1
!  elseif (a>0)then
!    n(1)=n(1)-1
!    n(2)=n(3)+1
!    n(3)=0
!  endif
!end

subroutine test_next_pow_nl!_apple
  integer n(3), k
  integer ind
  
  k = 4
  n(1) = k
  n(2) = 0
  n(3) = 0
  
  do ind = 1,sumfac(k+1)
     if(ind>1)call nextpown(n)
     print'(3I2,2I5)',n,ind, apple_g(n)
     enddo
     
end


                                                                       subroutine main2; call test_nl; end

subroutine test_nl
    integer ll(3), nn(3),nn_2(3), n_2, l, it, numbers(100)
    !integer i
    logical proceed
    nn = [7,3,6]*2
    nn_2 = nn/2
    n_2 = sum(nn_2)
    !l = 4
    do l = 1, n_2
        call init_nl(nn_2,l,ll)
        
        print*
        
        print*, '__nn/2_=_['//str(nn_2,2)//']__l_=_'//str(l)//"__"
        
        proceed = .true.
        it=0
        do while (proceed)
            it=it+1
            print*, 'll = '//str(ll,2)//"  finder(ll) = "//str(finder(ll))
            call next_pow_nl(nn_2,ll,proceed)
        enddo    
        
        print*, "# entries = "//str(it)
        numbers(l)=it
    enddo
    print*, "# entries = "//str(numbers(1:l-1),3)
    
    !print*, "choices "//str( choose(n_2-l,l)/fac(ll(1))/fac(ll(2))/fac(nn_2(3)) )
end subroutine

subroutine init_nl(nn_2, l, ll)
    
    integer, intent(in)  :: nn_2(3), l
    integer, intent(out) :: ll(3)
    
    integer ls, i
    
    ls = l
    
    ll=0
    do i = 1,3 
      do while (ls>0 .and. nn_2(i)>ll(i))
        ll(i) = ll(i)+1
        ls = ls -1
      enddo
    enddo
    
end subroutine




subroutine next_pow_nl(nn_2,ll, proceed)
    integer, intent(in)    :: nn_2(3) 
    integer, intent(inout) :: ll(3) 
    integer a,b,c , ra,rb,rc 
    logical proceed
    
    proceed = .true.
    ra=nn_2(1);rb=nn_2(2);rc=nn_2(3)
    a=ll(1);b=ll(2);c=ll(3)
    
    if (b>0 .and.rc>c)then
        b = b-1
        c = c+1
    elseif (a>0 .and.rb>b)then
        a = a-1
        b = b+1
        do while (c>0 .and.rb>b)
            b = b+1
            c = c-1
        enddo
    elseif (a>0 .and.rc>c)then
        a = a-1
        c = c+1
    else
        proceed = .false.
    endif
    ll = [a,b,c]
end
















! bad but not yet trash !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function brakk_bad(n,l) result(res)
    integer,intent(in) :: n, l
    integer n_2l, a,b,  i
    integer*8 p1,p2
    integer res
    n_2l = n-2*l
    b = max(l,n_2l)
    a = min(l,n_2l)
    
    p1=1
    p2=1
    
    do i = b+1,n
     p1 = p1*i
    enddo
    
    do i = 2, a
     p2 = p2*i
    enddo
    
    res = int(p1/p2)/2**l
end function

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

function nextpov(nn)
  integer a,b,c, nn(3), nextpov(3)
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
  
  nextpov(1) = a
  nextpov(2) = b
  nextpov(3) = c

  
end

endmodule
