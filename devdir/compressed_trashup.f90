

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

