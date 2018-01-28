module compressed_utils

use printer_mod, only: str, printer, printo
use compressed_arrays 
implicit none

integer, parameter :: dp = kind(0d0)
private dp
!public test_apple_g

contains !//////////////////////////////////////////////////

subroutine main
call print_trace_lengths    
end subroutine


subroutine print_trace_lengths
    integer i , n, acc, vall, temp(0:100), n05
    
    n = 4
    do n = 0, 15
        acc=0
        n05=n/2
        do i = 1, n05
            vall = len_(n-2*i)
            acc = acc+vall
            temp(i) =vall
            !print*, vall, acc
        enddo
        print*, "n="//str(n,2)//", len="//str(acc,3)//" <= "//str(temp(1:n05))!//" => "//str(sum(temp(1:n05)))
    enddo
    
    call printo([1, 23, 54, 1, 5, 7, 8, 909, 87, 0],0)
    print*, str([1, 23, 54, 1, 5, 7, 8, 909, 87, 0],0)
end subroutine


function del(a,b)
    integer a, b, del
    del=0
    if (a==b) del=1
end


subroutine print_trace_index_matrix(k)
    integer k
    integer j1,j2, icol, irow, ileng,v2(sumfac(k+1)), v4(sumfac(k+1)),v3(sumfac(k+1)), tleng !((k+1)*(k+2)/2),
    integer matrix(sumfac(k+1),sumfac(k/2+1))
    
    ileng = sumfac(k+1)
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
        
        
        !write(*,'(*(I3,a),a)'),((icol-1)*2 + v4(icol) + irow + (v3(icol))*(v2(irow)),',', icol = 1,tleng)
        do icol = 1, tleng
            matrix(irow,icol) = (icol-1)*2 + v4(icol) + irow + (v3(icol))*(v2(irow))
        enddo
    enddo
    print'(a)',"integer, parameter :: tmatrix("//str(ileng)//","//str(tleng)//") = reshape([ &"
    call printo(matrix,0,1)
    print'(a)',"], shape(tmatrix), order=[2,1])"

    
end subroutine

subroutine print_long_index_matrix(k)
    !integer, parameter :: k = 7
    integer k
    integer j1,j2, icol, irow, ileng,v2(sumfac(k+1)), vall, linelen
    
    ileng=len_(k)
    
    v2 = [((j1, j2=1,j1),j1 = 1,k+1)]-1 !for index matrix
    
    print*,"______"
    print'(a,*(I3))',"v2: ", v2
    
    
    print*
    !do irow = 1, ileng
    !    write(*,'(*(I2,a),a)'),(irow + icol -1 + (v2(icol))*(v2(irow)),',', icol = 1,ileng)
    !enddo
    
    linelen = 0
    print'(a)', "integer, parameter :: matr("//str(ileng)//","//str(ileng)//") = reshape( [ &"
    do irow = 1, ileng
        do icol = 1, ileng
            vall = irow + icol -1 + (v2(icol))*(v2(irow))
            linelen = linelen + intsize(vall)+1
            
            if (linelen>110) then 
                write(*,'(a)') "& "
                linelen=0
            endif
            !write(*,'(a)',advance="no") str(vall)
            
            if (irow*icol < ileng**2) then
                write(*,'(a)',advance="no") str(vall)//","
            else
                write(*,'(a)',advance="no") str(vall)
            endif
            
        enddo
    enddo
    print'(a)', ""
    print'(a)', "], shape(matr) )"
    
    
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

subroutine print_square_index_matrix(k)
    !integer, parameter :: k = 7
    integer k
    integer j1,j2, icol, irow, ileng,v2(len_(k)), matrix(len_(k),len_(k))
    
    ileng=len_(k)
    
    v2 = [((j1, j2=1,j1),j1 = 1,k+1)]-1 !for index matrix
    
    print*,"______"
    print'(a,*(I3))',"v2: ", v2
    
    
    print*
    !do irow = 1, ileng
    !    write(*,'(*(I2,a),a)'),(irow + icol -1 + (v2(icol))*(v2(irow)),',', icol = 1,ileng)
    !enddo
    
    do irow = 1, ileng
        do icol = 1, ileng
            matrix(icol,irow) = irow + icol -1 + (v2(icol))*(v2(irow))
            
        enddo
    enddo
    print'(a)', "integer, parameter :: matr("//str(ileng)//","//str(ileng)//") = reshape( [ &"
    call printo(matrix,0,1)
    print'(a)', "], shape(matr) )"
    
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

function intsize(ii) 
    integer ii, intsize
    intsize = ceiling(log10(dble(ii+1)))
end


subroutine print_apple_g(imax)
    integer i, j, n(3), it, gsiz, linelen, gval, glen, imax
    character(1) zero
    character(:), allocatable :: gname
    character(1000) :: allgs
    print*,"____"
    print*,">>>>  g-arrays:"
    it = 0
    zero = "0"
    allgs="["
    
    do i = 0,imax
      glen = sumfac(i+1)
      if(i>9)zero=""
      gname = "g"//trim(zero)//str(i)
      write(*,'(a)',advance="no") "integer, parameter :: "//gname//"("//str(glen)//") = ["
      if (i>6) write(*,'(a)') "& "
      
      if(i<imax)then 
        allgs = trim(allgs)//gname//","
      else
        allgs = trim(allgs)//gname//"]"
      endif
      
      linelen = 0
      
      n = 0
      n(1)= i
      do j = 1, glen
        if(j>1)call nextpown(n)
        
        gval = apple_g(n)
        gsiz = intsize(gval)
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
    !print*, "entries", it
    print*
    gname = "gg"//str(imax)
    print'(a)', "integer, parameter :: "//gname//"("//str(it)//") = "//trim(allgs)
    print'(a)', "integer, parameter :: gg_(pos_(matk+1)) = "//gname//"(1:pos_(matk+1))"
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




pure function sumfac(u)
    !third row in pascal matrix (for tricorn lengths)
    integer sumfac
    integer, intent(in) :: u
    sumfac = u*(u+1)/2 ! 2=2!
endfunction


function hh(n1, n2) !Make a matrix of this sheeeet
    integer n1(3),n2(3), nn(3), hh
    nn = n1+n2
    hh = vfac(nn)/( vfac([n1,n2]) )
    !( fac(nn(1))*fac(nn(2))*fac(nn(3)) ) / ( fac(n1(1))*fac(n1(2))*fac(n1(3)) * fac(n2(1))*fac(n2(2))*fac(n2(3)) )
end



function hh_abc(a1,b1,c1, a2,b2,c2) result(hh) !Make a matrix of this sheeeet
    integer a1,b1,c1, a2,b2,c2
    integer aa,bb,cc, hh
    aa = a1+a2
    bb = b1+b2
    cc = c1+c2
    
    hh = ( fac(aa)*fac(bb)*fac(cc) ) / ( fac(a1)*fac(a2)*fac(b1)*fac(b2)*fac(c1)*fac(c2) )

end

function hhh(i,j,ki,kj) !Make a matrix of this sheeeet
    integer, intent(in) :: i,j, ki,kj
    integer kij, hhh, cc, gi,gj,gij
    kij = ki+kj
    cc = choose(kij,ki)
    gi = gg_(pos_(ki)+i)
    gj = gg_(pos_(kj)+j) 
    gij = gg_(pos_(kij) + mm_(i,j))
    hhh = (gi*gj*cc) / gij
    
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
    print*, "ABOVE: print_traces ---------------------------------------"
end subroutine
    


subroutine test_polyfinder
integer nn(3), i, n, it
it = 0
print*, "key, pos, finder, polyf, it, polyf-it "
do n = 0, 7
    nn=[n,0,0]
    do i = 1, len_(n)
        if (i>1)call nextpown(nn)
        it = it+1
        print'(a,*(I4))'," ["//str(nn)//"]",pos_(n), finder(nn), polyfind(nn), it, polyfind(nn)-it
    enddo
enddo
print*, str(pos_)
end

function finder(n) result(row)
   ! Given nx,ny,nz, returns corresponding row in tricorn _vector_. 
   integer n(3), row, n23 !rank, i
   
   n23 = n(2)+n(3)
   row = 1 + n(3) + n23*(n23+1)/2
   
   !do i = 1,rank - n(1) !addition factorial
   !  row = row + i
   !enddo
end

function polyfind(nn) result(row)
   ! Given nx,ny,nz, returns corresponding row in tricorn _polytensor_. 
   integer nn(3), row, n32, n321 !rank, i
   
   n32 = nn(2)+nn(3)
   row = 1 + nn(3) + len_(n32-1)
   
   n321 = n32+nn(1)
   row = row + pos_(n321)
   
end




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

function fac(inn) !short iteger factorial
    integer inn, i, fac
    fac=1
    do i = 1, inn
      fac = fac*i
      enddo

end

function lfac(inn) result(fac) !long integer factorial
    integer*8 inn, i, fac
    fac=1
    do i = 1, inn
      fac = fac*i
      enddo

end

function vfac(vec) result(prod)
    ! Returns the product of the factorials of the elements of an integer vector
    integer prod, vec(:), i, j, si
    si = size(vec)
    prod = 1
    do i = 1, si
        do j = 2, vec(i)
            prod = prod*j
        enddo
    enddo
end function



function vecbrakk(nn,ll) result(res)
    ! Returns the product of the brakks (see 'brakk' funciton) of the elements of two same-length integer vectors
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


subroutine nextpown(nn)
    ! Updates the argument to the next power vector in lex ordering
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


subroutine init_pow_nl(nn_2, l, ll)
    ! Returns the first l-vector given an n/2-vector
    ! l-vector is the powers of delta; n-vector is the powers of xyz
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


endmodule
