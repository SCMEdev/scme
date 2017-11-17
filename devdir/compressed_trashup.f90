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

