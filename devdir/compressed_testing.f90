module compressed_testing

use compressed_arrays
use compressed_tensors, bad=>main
use compressed_utils, bad=>main

implicit none


contains !//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

subroutine main
call &
!test_polyfinder
test_detracer
!all_tests
!get_traces !-----------
!call &
!test_brakk
!print_traces(5)
!print_trace_keys
!test_intfac_ff
end subroutine


subroutine all_tests
    
    !call get_traces
    call test_sumfac
    !call test_nextpow(7)
    !call test_nextpow_v_wn(10)

    
    call test_subdiv_pow_h(5,5)
    call test_subdiv_pow_h(4,3)
    call test_subdiv_pow_h(3,4)
    
    call test_next_key2n(4,1)
    call test_sorted
    call test_expand_compress
    call test_brakk
    call test_nextpown
    call test_next_pown
    call print_trace_keys
    
    call h_testing
    call test_hhh(3,4)
    
    call test_inner_symouter
    call test_inner
    call test_symouter
    
    call test_old_field
    call test_potgrad
    call test_mp_pot
    
    call printo(tmm_,0)
    call printo(mm_,0,0)
    call printo(mm_,0,1)
    call printo([1,2,3,4,5,6,77,7777,777,9],0,0)
    !call test_intfac_ff!takes two command line arguments
    
    call test_fac(5)
    call test_choose
    
    call print_product_index_matrix (6)
    
    call test_matr(7)
    call test_apple_g
    
    call test_rrr
    
    call test_polyfinder
    
    
    print*, '------------------------------------------------------------'
    print*, 'ABOVE: all_tests -------------------------------------------'
    print*, '------------------------------------------------------------'

end subroutine


! UTILS ____________________________________________________________________________________________________________________________
subroutine test_sumfac
    integer i
    
    do i = 1,10
        print*, sumfac(i) , i, 1
    enddo
    print*, 'ABOVE: test_sumfac -------------------------------------------'
endsubroutine

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
    print*, 'ABOVE: test_choose -------------------------------------------'
end

subroutine test_apple_g
    print*, "xxyy", apple_g([2,2,0])
    print*, "xy", apple_g([1,1,0])
    
    print*, "xxyz", apple_g([2,1,1])
    print*, apple_g([1,1,0])
    print*, apple_g([1,0,1])
    print*, apple_g([0,1,1])
    
    print*, "xxxy", apple_g([3,1,0])
    print*, 'ABOVE: test_apple_g -------------------------------------------'
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
    
    
    
    print*, 'ABOVE: test_next_key2n -------------------------------------------'
    contains !// ////////////////////////
      
      
      
      subroutine pprint
         print*, "key=",(str(key(j))//" ",j=1,rank ),"  nn=",(str(nn(j))//" ",j=1,3 ),&
                 "  row:"//str(i), "  found row:"//str(finder(nn)), "  g:"//str(apple_g(nn))
      end subroutine
      
end

subroutine test_sorted()
    integer, allocatable :: key(:)!, key2(:)
    integer rank
    key = [4,5,2,1,3,1,7,1,3,2,6,1]
    rank = size(key)
    print '('//str(rank)//'I2)', sorted(key)
    print*, 'ABOVE: test_sorted -------------------------------------------'
end subroutine

subroutine test_expand_compress
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
    
    
    print*, 'ABOVE: test_expand_compress -------------------------------------------'
end subroutine

subroutine test_fac(k)
    integer k, i
    do i = 1, k
        print*, fac(i)
    enddo
    
    print*, 'ABOVE: test_fac -------------------------------------------'
end

subroutine test_brakk
    integer n, l
    do n = 1, 10!20
        do l = 1, n/2
            print*, n,l, brakk(n,l), brakk_naive(n,l)
        enddo
    enddo
   print*, 'ABOVE: test_brakk -------------------------------------------'
end subroutine

function brakk_naive(n,l) result(res)
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
        
   print*, 'ABOVE: test_nextpown -------------------------------------------'
end

subroutine test_next_pown!_nl!_apple
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
        
    print*, 'ABOVE: test_next_pown -------------------------------------------'
end


subroutine test_subdiv_pow_h(ii,kii)
    integer ii, kii
    integer i, j, it
    integer n1(3), n2(3), nn(3)
    integer sfii, sfkii
    
    n1=0
    n1(1)=ii
    
    
    nn=0
    nn(1)=ii+kii
    
    
    sfii = sumfac(ii+1)
    sfkii = sumfac(kii+1)
    
    it = 0
    do i = 1, sfii
        if (i>1)call nextpown(n1)
        
        n2=0
        n2(1)=kii
        
        do j = 1, sfkii
            if (j>1)call nextpown(n2)
            it = it+1
            nn = n1+n2
            !*
            write(*,'(I6)', advance="no") hh(n1,n2) - hhh(i,j,ii,kii)
            enddo
        print*,""
    enddo
    !**    
    Print*, "ABOVE: test_subdiv_pow_h------------------------------------------------------------"
end

!*
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
            !                 ( gg(pos_(ii)+i)*gg(pos_(kii)+j)*choose(ii+kii,ii) ) / gg(pos_(ii+kii) + matr(i,j)), & 
            !                 00,&
            !                 apple_g([a1,b1,c1]), apple_g([a2,b2,c2]), choose(ii+kii,ii), apple_g([aa,bb,cc]), &
            !                 00, &
            !                 gg(pos_(ii)+i), gg(pos_(kii)+j), choose(ii+kii,ii), gg(pos_(ii+kii) + matr(i,j)), & 
            !                 00
            
            !write(*,'(I6)', advance="no") ( apple_g([a1,b1,c1])*apple_g([a2,b2,c2])*choose(ii+kii,ii) )/apple_g([aa,bb,cc])
            !write(*,'(I6)', advance="no") ( gg(pos_(ii)+i)*gg(pos_(kii)+j)*choose(ii+kii,ii) ) / gg(pos_(ii+kii) + matr(i,j))
            
            
            !write(*,'(I6)', advance="no") ( gg(pos_(ii)+i)*gg(pos_(kii)+j)*choose(ii+kii,ii) ) / gg(pos_(ii+kii) + matr(i,j)) &
            !                              - ( apple_g([a1,b1,c1])*apple_g([a2,b2,c2])*choose(ii+kii,ii) )/apple_g([aa,bb,cc])
            
            !write(*,'(I6)', advance="no") hh(a1,b1,c1, a2,b2,c2) - hhh(i,j,ii,kii)
            !write(*,'(I6)', advance="no") hh(n1(1),n1(2),n1(3), n2(1),n2(2),n2(3)) - hhh(i,j,ii,kii)

!**
    !do here!
    !print'(*(I3))', vv(1:sumfac(ii+kii+1))
    !print'(*(I3))', vv1(1:sumfac(ii+kii+1))
    !print*, sumfac(ii+kii+1)


! TENSORS __________________________________________________________________________________________________________________________
subroutine test_detracer
    integer, parameter :: n = 5
    real(dp) ::  compvec(len_(n)), AA(len_(n)),newvec(len_(n))
    

    !testvec = [1,0,0,2,0,3]
    
    !testvec = [ 1d0,3d0,5d0,6d0,7d0,9d0,8d0,5d0,4d0,2d0 ]
    !testvec = [ 1d0,3d0,5d0,6d0,7d0,9d0,8d0,5d0,4d0,2d0, 4d0, 2d0, 6d0, 7d0, 1d0 ]
    !testvec = [ 1d0,2d0,3d0,4d0,5d0,6d0,7d0,8d0,9d0,10d0,11d0,12d0,13d0,14d0,15d0]/10d0
    AA = [ 1d0,2d0,3d0,4d0,5d0,6d0,7d0,8d0,9d0,10d0,11d0,12d0,13d0,14d0,15d0,16d0,17d0,18d0,19d0,20d0,21d0]/10d0

    newvec = detracer(AA,n)
    
    
    compvec = opdetr(AA,n)
    print'(a,*(f10.2))', 'testvec', AA
    print'(a,*(f10.2))', 'newvec ', newvec
    print'(a,*(f10.2))', 'compvec', compvec
    print'(a,*(f10.2))', 'frac   ', newvec/compvec
    print*, "ABOVE: test_detracer ---------------------------------------"

end


subroutine h_testing
    integer i, k , n(3)
    k=5
    n = 0
    n(1) = k
    do i = 1, sumfac(k+1)
        if(i>1) call nextpown(n)
        print'(3I3,a,3I4,a,I4)',n,",  ", fac(n(1)), fac(n(2)), fac(n(3)),",  ", fac(n(1))*fac(n(2))*fac(n(3))
    enddo
    
    print*, "ABOVE: h_testing ----------------------------------------------------------------------"
end

subroutine test_hhh(k1, k2)
    integer i1, i2, k1,k2
    integer gp1, gp2, gp12, h, ch
    gp1 = pos_(k1)
    gp2 = pos_(k2)
    gp12 = pos_(k1+k2)
    ch = choose(k1+k2,k1)
    do i1 = 1, sumfac(k1+1)
      do i2 = 1, sumfac(k2+1)
        h = ( gg_(gp1+i1)*gg_(gp2+i2)*ch ) / gg_(gp12 + mm_(i1,i2))
        write(*,'((I4))', advance="no") h-hhh(i1,i2,k1,k2)
        enddo
      print*,""
      enddo
    print*, "above, TEST_HHH ------------------------------------------------------------------"
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
    
    print*, "ABOVE: test_inner ------------------------------------------------------------------"

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
    
    print*, "ABOVE: test_symouter ------------------------------------------------------------------"
    
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
    real(dp) :: rrr(pos_(kmax+1)), quad(6), octa(10), rnorm, rsqe
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
    
    
    
    print*, 'ABOVE: TEST POTGRAD-----------------------------------------------'
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
    real(dp) :: rrr(pos_(kmax+1)), rnorm, rsqe
    real(dp) :: rinvv(2*(kmax+nmax)+1) !, rinvv1(2*(kmax+nmax)+1), rinvv2(2*(kmax+nmax)+1)
    real(dp) :: dd(3), dr, cq(6), co(10), ch(15)
    !real(dp) :: quad(6), octa(10),
    
    integer n,k
    integer p1,p2,q1,q2
    real(dp) :: phi_old(pos_(6)),phi(pos_(6)), qq(pos_(5))
    
    
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
    qq = [0d0,dd,cq,co,ch]
    
    !print*,"q-thigs"
    !n=1
    !q1 = pos_(n)+1
    !q2 = pos_(n+1)
    !print'(*(I3))',n, q1,q2
    !
    !qq(q1:q2) = dd
    !
    !n=2
    !q1 = pos_(n)+1
    !q2 = pos_(n+1)
    !print'(*(I3))',n, q1,q2
    !
    !qq(q1:q2) = cq
    
    
    print*
    print'(*(f10.4))', qq(2:)
    print*, 'size(qq)',size(qq), 1+3+6+10+15
    
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
    
    
    phi_old(1) = 0
    phi(1) = 0
    
    
    k=1
    p1 = pos_(k)+1
    p2 = pos_(k+1)
    print'(6(I3))',k, p1,p2
    phi_old(p1:p2) = d1v(:,2)
    
    k=2
    p1 = pos_(k)+1
    p2 = pos_(k+1)
    print'(6(I3))',k, p1,p2
    phi_old(p1:p2) = opdetr(compress( reshape(      d2v(:,:,2),shape=[3**2]),2),2)
    
    k=3
    p1 = pos_(k)+1
    p2 = pos_(k+1)
    print'(6(I3))',k, p1,p2
    phi_old(p1:p2) = opdetr(compress( reshape(    d3v(:,:,:,2),shape=[3**3]),3),3)
    
    k=4
    p1 = pos_(k)+1
    p2 = pos_(k+1)
    print'(6(I3))',k, p1,p2
    phi_old(p1:p2) = opdetr(compress( reshape(  d4v(:,:,:,:,2),shape=[3**4]),4),4)
    
    k=5
    p1 = pos_(k)+1
    p2 = pos_(k+1)
    print'(6(I3))',k, p1,p2
    phi_old(p1:p2) = opdetr(compress( reshape(d5v(:,:,:,:,:,2),shape=[3**5]),5),5)
    
    
            
    
    phi=0
    do n = 1, 4
        q1 = pos_(n)+1
        q2 = pos_(n+1)
            do k = 1,5
                p1 = pos_(k)+1
                p2 = pos_(k+1)
                print'(6I3)',n,k,q1,q2, p1,p2
                phi(p1:p2) = phi(p1:p2) + opdetr(potgrad(qq(q1:q2),n,k,rinvv,rrr),k)
            enddo
    enddo
    
    
    
    
    print'(a,*(g30.15))', 'phi    ',phi(2:)
    print'(a,*(g30.15))', 'phi_old',phi_old(2:)
    print'(a,*(g30.15))', 'phi+old',phi_old(2:)+phi(2:)
    
    
    
    
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
    

    
   print*, 'ABOVE: test_mp_pot  --------------------------'
     
end


subroutine test_matr(nn)
    integer i, nn, maxi
    if(nn>7)stop"rank cant be larger than 7 in matr intdex matrix"
    maxi = sumfac(nn)
    do i = 1, maxi
      print'(28I3)', mm_(i,1:maxi)
    enddo
   print*, 'ABOVE: test_matr ----------------------------------'
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
    real(dp) :: rrr(pos_(k+1)), rr(3), rrr3(10), rrr2(6), rrr5(len_(5)), rrr32(len_(5))
    real(dp) :: rrr4(len_(4)), rrr31(len_(4)), rrr22(len_(4)), rrr21(len_(3))
    integer i, p1, p2
    call random_seed(put=[2,234,1,5,435,4,5,42,3,43,432,4,3,5,23,345,34543])
    call random_number(rr)
    rr = [1d0,2d0,3d0]
    rr = [0.3810985945,0.5295087287,0.852367145402366]
    print*, size(rrr)
    
    call vector_powers(k,rr,rrr)
    !call rrpow(rr,k,rrr)
    !print'(*(f12.4))', rrr  
    call printer(rrr,'rrr',1)
    
    
    i = 2
    p1 = pos_(i)+1
    p2 = pos_(i+1)
    rrr2=rrr(p1:p2)
    
    i = 3
    p1 = pos_(i)+1
    p2 = pos_(i+1)
    rrr3=rrr(p1:p2)
    
    i = 4
    p1 = pos_(i)+1
    p2 = pos_(i+1)
    rrr4=rrr(p1:p2)
    
    i=5
    p1 = pos_(i)+1
    p2 = pos_(i+1)
    
    rrr5=rrr(p1:p2)
    
    rrr32 = symouter(3,2,rrr3,rrr2)
    rrr22 = symouter(2,2,rrr2,rrr2)
    rrr31 = symouter(3,1,rrr3,rr)
    rrr21 = symouter(2,1,rrr2,rr)
    
    
    call printer(rrr5,"5th",1)
    call printer(rrr32,"symouter(3,2)",1)
    call printer(rrr32/rrr5,"frac32",1)
    call printer(rrr22/rrr4,"frac22",1)
    call printer(rrr31/rrr4,"frac31",1)
    call printer(rrr21/rrr3,"frac21",1)
    
    print*, "above, TEST RRR ------------------------------------------------------------------"

end










end module
