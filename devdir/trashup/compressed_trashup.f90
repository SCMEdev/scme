subroutine detracing_old
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer l, i, nn(3), ll(3)
 !   real(dp) :: testvec(21), newvec(21), compvec(21), tvec(1000)
    real(dp) su
    
        integer, parameter :: n = 5
    !real(dp) :: testvec(21), newvec(21), compvec(21), tvec(1000)
    real(dp) :: testvec(len_(n)), newvec(len_(n)), compvec(len_(n)), tvec(1000)

    
    !first implementation
    integer g1,g2, cols, rows, key(3), nh, nnh(3), tpos, t1, n_2l
    logical proceed
    
    !real(dp) :: testvec(15), newvec(15), compvec(15), tvec(1000)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
    !testvec = [1,0,0,2,0,3]
    !testvec = [ 1d0,3d0,5d0,6d0,7d0,9d0,8d0,5d0,4d0,2d0 ]
    !testvec = [ 1d0,3d0,5d0,6d0,7d0,9d0,8d0,5d0,4d0,2d0, 4d0, 2d0, 6d0, 7d0, 1d0 ]
    !testvec = [ 1d0,2d0,3d0,4d0,5d0,6d0,7d0,8d0,9d0,10d0,11d0,12d0,13d0,14d0,15d0]/10d0
    testvec = [ 1d0,2d0,3d0,4d0,5d0,6d0,7d0,8d0,9d0,10d0,11d0,12d0,13d0,14d0,15d0,16d0,17d0,18d0,19d0,20d0,21d0]/10d0
    
    tvec=0
    
    do l = 0,n/2 !fill in the trace-arrays (polytensor)
        n_2l = n-2*l        !rank of trace
        rows = len_(n_2l)   !lenght of trace (sub-)array
        t1 = pos_(n_2l)     !position in trace-polytensor
        
        cols = len_(l)  !first entries of tm-row that sample A
        
        g1 = pos_(l)+1  !
        g2 = pos_(l+1)  !section of apple-g polytensor
        
        print*, "l = "//str(l)
        print'(a,*(I4))', "g1,g2:", g1,g2
        print'(a,*(I4))', "gg=", gg_(g1:g2)
        
        do i = 1, rows !over rows of trace-index-matrix (tm)
            
            tvec(t1+i) = dot_product(testvec(tmm_(i,1:cols)), gg_(g1:g2)) !tm-row * apple-g
            print'(a,*(I3))',"trace index  ", tmm_(i,1:cols)
            print'(a,*(f7.3))',"trace vallues", testvec(tmm_(i,1:cols))
            
            !print'(a,*(I3))',"apple-g num", gg_(g1:g2)
        enddo
        print'(a,*(f7.3))', "partial trace array", tvec(t1+1:t1+rows)    
        
    enddo
    print'(a,*(f7.3))', "full trace array",tvec(1:pos_(n+1)) !float måste ha en rutin !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    
    
    nn = [n,0,0]
    do i = 1, len_(n)
        if(i>1)call nextpown(nn)
        nnh = nn/2
        nh = sum(nnh)
        
        su = 0
        do l = 0, nh
            call init_pow_nl(nnh,l,ll)
            n_2l = n-2*l
            
            tpos = pos_(n_2l)
            proceed = .true.
            
            do while (proceed)
                key = nn - 2*ll
                t1 = finder(key)
    !            print*, "nn="//str(nn)//", nn/2="//str(nnh)//", ll="//str(ll)//", key="//str(key)//&
    !            ", pos="//str(tpos,3)//" =>idx="//str(t1,3)//" =>val="//str(tvec(tpos+t1))
                print*, t1, tpos+t1
                su = su + (-1)**l * intff(2*(n-l)-1,1) * vecbrakk(nn,ll) * tvec(tpos+t1)
                call next_pow_nl(nnh,ll,proceed)
            enddo
        enddo
        !print*, "su", su
        newvec(i) = su
        
        
    enddo
    
    
    
    
     
    compvec = detrace_a(testvec,n)
    print'(a,*(f10.2))', 'testvec', testvec
    print'(a,*(f10.2))', 'newvec ', newvec
    print'(a,*(f10.2))', 'compvec', compvec
    print'(a,*(f10.2))', 'frac   ', newvec/compvec
    
    
    print*, "ABOVE: detracing_old ---------------------------------------"
end subroutine


subroutine detracing_new
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer l, i, nn(3), ll(3)
    integer, parameter :: n = 5
    !real(dp) :: testvec(21), newvec(21), compvec(21), tvec(1000)
    real(dp) :: testvec(len_(n)), newvec(len_(n)), compvec(len_(n)), AA(len_(n)), tvec(1000)
    !real(dp) :: testvec(15), newvec(15), compvec(15), tvec(1000)
    real(dp) su
    
    !reimplementation
    integer l1,l2,l3,l32, ti, ni
    
    integer m, mi, mm(3), li, mpo, br3, br32, br, n1,n2,n3
    !real(dp) :: tAA(1000)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    
    !testvec = [1,0,0,2,0,3]
    
    !testvec = [ 1d0,3d0,5d0,6d0,7d0,9d0,8d0,5d0,4d0,2d0 ]
    !testvec = [ 1d0,3d0,5d0,6d0,7d0,9d0,8d0,5d0,4d0,2d0, 4d0, 2d0, 6d0, 7d0, 1d0 ]
    !testvec = [ 1d0,2d0,3d0,4d0,5d0,6d0,7d0,8d0,9d0,10d0,11d0,12d0,13d0,14d0,15d0]/10d0
    testvec = [ 1d0,2d0,3d0,4d0,5d0,6d0,7d0,8d0,9d0,10d0,11d0,12d0,13d0,14d0,15d0,16d0,17d0,18d0,19d0,20d0,21d0]/10d0
    
    tvec=0
    AA = testvec
    do l = 0, n/2
        
        m = n-2*l
        mpo = pos_(m)
        
        mm = [m,0,0]
        do mi = 1,len_(m)
            if(mi>1)call nextpown(mm)
            
            su = 0
            ll = [l,0,0]
            do li = 1,len_(l)
                if(li>1)call nextpown(ll)
                
                
                nn = mm + 2*ll
                ni = finder(nn) 
                !print*, "ni=",ni
                
                su = su + fac(l)/vfac(ll) * AA(ni)
            enddo
            tvec(mpo+mi) = su
            !print*,"su=",su
            !print*, "mpo+mi", mpo+mi
        enddo
    enddo
    
    print'(a,*(f7.3))', "full trace array",tvec(1:pos_(n+1)) !float måste ha en rutin !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    n1 = n; n2 = 0; n3 = 0
    do i = 1, len_(n)
        if(i>1)call nextpow(n1,n2,n3)
        
        su = 0
        do l3 = 0, n3/2
            br3 = brakk(n3,l3)
            mm(3) = n3-2*l3
            
            do l2 = 0, n2/2
                br32 = br3 * brakk(n2,l2)
                l32 = l3+l2
                mm(2) = n2-2*l2
                
                do l1 = 0, n1/2
                    br = br32 * brakk(n1,l1)
                    l = l32+l1
                    mm(1) = n1-2*l1
                    
                    
                    ti = polyfind(mm)!finder(nn-2*ll)
                    su = su + (-1)**l * intff(2*(n-l)-1,1) * br *  tvec(ti)
                    !print*, "vecbrakk(nn,ll)", vecbrakk(nn,ll)
                    !print*, "ti=",ti, "tvec(ti)=", tvec(ti)
                    !print*, "doing"
                enddo
            enddo
        enddo
        
        newvec(i) = su
        
    enddo
    
    compvec = detrace_a(testvec,n)
    print'(a,*(f10.2))', 'testvec', testvec
    print'(a,*(f10.2))', 'newvec ', newvec
    print'(a,*(f10.2))', 'compvec', compvec
    print'(a,*(f10.2))', 'frac   ', newvec/compvec
    
    
    print*, "ABOVE: detracing new ---------------------------------------"
end

subroutine test_nextpow_v_wn(k)
    integer a,b,c, k, nn(3) , nn2(3)
    integer ind
    
    a=k
    b=0
    c=0
    nn = [a,b,c]
    nn2 = nn
    
    do ind = 1,sumfac(k+1)
        if(ind>1)then 
            call nextpow(a,b,c)
            nn = nextpov(nn)
            call nextpown(nn2)
        endif
        print*, "a,b,c:"//str([a,b,c])//", nn="//str(nn)//", nn2="//str(nn2)//",    index="//str(ind)
    enddo
     
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


subroutine pascal_matrix(size)
  integer :: size, row1(size), row2(size), i1,i2, i12, pmat(size,size), temp, choom(size,size)
  
  pmat(1,:)=1
  pmat(:,1)=1
  
  do i1 = 2,size
    
    do i2 = 2,size
      pmat(i1,i2) = pmat(i1,i2-1)+pmat(i1-1,i2)
    enddo
    
    
  enddo
  !do i1 = 1, size
  !  print'(*(I6))',pmat(i1,:)
  !enddo
  
  print*
  
  do i1 = 1, size
    print'(*(I6))', (pmat(i1,i2),i2 = 1, size)
  enddo
  
  print*
  
  choom=0
  print*
  do i1 = 1, size
    do i2 = 1, size
      if(i2.ge.i1)then
        temp = pmat(i1,i2-i1+1)
      else
        temp = 0
      endif
      write(*,'(1(I6))',advance="no") temp
    enddo
    print*
  enddo
  
  !choom=0
  print*, "hej"
  do i1 = 1, size
    do i2 = i1, size
      !if(i2.ge.i1)then
        temp = pmat(i1,i2-i1+1)
        choom(i1,i2) = temp
        choom(i2,i1) = temp
    enddo
  enddo
  
  do i1 = 1, size
    print'(*(I6))', (choom(i1,i2),i2 = 1, size)
  enddo

  
  print*
  print'(*(I6))', temp_position(0:)
  print'(*(I6))', pos00(0:)
  
  print*, sumfacfac(5), pos00(5)
  
  print*, "hej" , 0**0, 0**7, 7**0
  
end


subroutine test_rrpow
    !call testing
    !call test_sumfac
    !call test_rpow
    integer, parameter :: k=5
    real(dp) :: rr(sumfacfac(k+1)-1)
    call rrpow([1d0,2d0,3d0],k,rr)
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


subroutine rrpow(r,k,rr) 
    ! Do not use this routine, it does not produce th outer products, since it lacks the scaling parameter hhh
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

function regouter(k1,k2,v1,v2) result(vout) 
! make this into a regular outer product of outer products of the r-vector. 
! on the other hand, this may not really be nessecary since we can allready produce the outer product array with the old routine. 
! and the operation is only valid for operations on the same source vector, so having a general routine may not be very usefull. 
! it can basically only be used for creating the rrr-array, which can be produced in a different way. 
! successive outer product of different vectors are not symmetric and can not be espressed in compressed form. 
    integer, intent(in) :: k1,k2
    real(dp), intent(in) :: v1((k1+1)*(k1+2)/2), v2((k2+1)*(k2+2)/2)
    real(dp) vout( (k1+k2+1)*(k1+k2+2)/2 )
    integer i, k12
    
    
    k12 = k1+k2

    
    
    
    vout=0
    
    do i = 1, sumfac(k12+1)
        
        !mij = matr(i,j)
        
        
        
        !vout(mij) = vout(mij) + v1(i)*v2(j) 
        
        !if(pri)print*,  'i,j , v1(i), v2(j), h',i,j, v1(i), v2(j), (gi*gj*cho)/gij
        
        enddo
    


end

!subroutine vector_powers(k,rr,rrr)
!    integer, intent(in) :: k
!    real(dp), intent(in) :: rr(3)
!    real(dp), intent(out) :: rrr(0:sumfacfac(k+1)-1)
!    integer i, p1, p2, p3, p4
!    rrr = 0
!    rrr(0) = 1
!    rrr(1:3) = rr
!    do i = 2,k
!      ! start/end positions of subvectors in the rrr vector. 
!      ! (Need +1 since they are 0-indexed, because usually i or j are added to them.) 
!      p1 = gpos(i-1)+1
!      p2 = gpos(i)
!      p3 = gpos(i)+1
!      p4 = gpos(i+1)
!      
!      !print'(*(I3))', i, p1, p2, p3, p4
!      
!      rrr(p3:p4) = symouter(i-1,1,rrr(p1:p2),rr)
!      enddo
!
!end

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

