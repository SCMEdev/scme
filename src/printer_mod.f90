! This module print tensors of rank 1 to 5 with the generic interface "printer of all the s ranks. 
module printer_mod
    !use data_types, only: dp!, h2o
    
    implicit none
    
    integer, parameter :: dp = kind(0d0)
    integer, parameter :: un=6
    character(*), parameter :: intf='(I7)'
    character(*), parameter :: dblf='f24.15' !double format
    
    ! Printer interfaces 
    interface printer
        module procedure p_real, p_int, p_real_vec, p_real_mat, p_real_3d, p_real_4d, p_real_5d!, p_h2o, p_h2os
    end interface
    
    ! Interface for tensor printed in given order
    interface printo !cat_tens
        module procedure print_order1, print_order2, print_order3, print_order4, print_integer_matrix
    end interface
    
    ! Interface for an str() funciton like in python
    interface str
        module procedure i2s, i2s_align, f2s, l2s, iv2s!, sp2s
    end interface
    
    ! Private/Public
    private
    public printer, xyz_hho_to_linear, str, printo !, printer_h2o_linear, h2o_to_linear
    
    
contains !/////////////////////////////////////////////


function i2s(inte) result(ch)
    integer, intent(in) ::  inte!, length
    character(:), allocatable :: ch
    character(50) temp
    write(temp,'(I50)') inte
    ch = trim(adjustl(temp))
end function

function i2s_align(inte, space) result(ch)
    integer, intent(in) ::  inte!, length
    integer, intent(in) :: space
    character(:), allocatable :: ch
    character(space) temp
    write(temp,'(I'//i2s(space)//')') inte
    ch = adjustr(temp)
end function



function iv2s(int_vec,space) result(ch)
    integer, intent(in) ::  int_vec(:)
    integer, intent(in), optional :: space!, length
    integer xspace
    character(:), allocatable :: ch
    character(1000) temp
    xspace = 5
    if(present(space)) xspace = space
    write(temp,'(*(I'//i2s(xspace)//'))') int_vec
    ch = trim(adjustl(temp))
endfunction


function f2s(rea,dec) result(ch)
    real(dp), intent(in) :: rea
    integer, intent(in), optional :: dec 
    character(:), allocatable :: ch
    character(30) temp
    integer decim
    
    if(present(dec))then
        decim=dec
    else
        decim=7
    endif
    
    write(temp,'(f30.'//i2s(decim)//')') rea
    ch = trim(adjustl(temp))
endfunction

!function sp2s(rea,dec) result(ch)
!    real, intent(in) :: rea
!    integer, intent(in), optional :: dec 
!    character(:), allocatable :: ch
!    character(30) temp
!    integer decim
!    
!    if(present(dec))then
!        decim=dec
!    else
!        decim=7
!    endif
!    
!    write(temp,'(f30.'//i2s(decim)//')') rea
!    ch = trim(adjustl(temp))
!endfunction


function l2s(logi) result(ch)
    logical logi
    character(:), allocatable :: ch
    if (logi)then
     ch = "True"
    else
     ch = "False"
    endif
endfunction


!/////////////////////////////// new version
subroutine print_order4(tensor,order)
    real(dp), intent(in)        :: tensor(3,3,3,3)
    integer,intent(in),optional :: order(4)
    
    integer       :: ic(4),ordc(4) !corrected order&index
    integer       :: i,j,k,l
    character(12) :: ctemp(3)
    character(2)  :: b1,b2
    
    if (present(order)) then
       !order = [1,2,3,4]
       ordc = order
       ordc(1) = order(2)
       ordc(2) = order(1)
    else
       ordc = [2,1,3,4]
    endif
    
    do l = 1,3
        ic(ordc(4)) = l
        print*, ""!"4-slice"!//str(l)
        do k = 1,3
            ic(ordc(3)) = k
            print*, ""!"  3-slice:"//str(k)
            do j = 1,3
                ic(ordc(2)) = j
                do i = 1,3
                    ic(ordc(1)) = i
                      
                    write (ctemp(i),'(f12.7)') tensor(ic(1),ic(2),ic(3),ic(4))
                    
                enddo
                
                call brackets(j,b1,b2)
                print*, '   '//b1//ctemp(1)//ctemp(2)//ctemp(3)//b2
            enddo
        enddo
    enddo
    
end subroutine

subroutine print_order3(tensor,order)
    real(dp), intent(in)        :: tensor(3,3,3)
    integer,intent(in),optional :: order(3)
    
    integer       :: ic(3),ordc(3) !corrected order&index
    integer       :: i,j,k
    character(12) :: ctemp(3)
    character(2)  :: b1,b2
    
    if (present(order)) then
       ordc(1) = order(2)
       ordc(2) = order(1)
       ordc(3) = order(3)
    else
       ordc = [2,1,3]
    endif
    
    do k = 1,3
        ic(ordc(3)) = k
        print*
        do j = 1,3
            ic(ordc(2)) = j
            do i = 1,3
                ic(ordc(1)) = i
                write (ctemp(i),'(f12.7)') tensor(ic(1),ic(2),ic(3))
            enddo
            
            call brackets(j,b1,b2)
            print*, '   '//b1//ctemp(1)//ctemp(2)//ctemp(3)//b2
        enddo
    enddo
    
end subroutine


subroutine print_order2(tensor,order)
    real(dp), intent(in)        :: tensor(3,3)
    integer,intent(in),optional :: order(2)
    
    !corrected order&index
    integer       :: ic(2),ordc(2) 
    integer       :: i,j
    
    character(12) :: ctemp(3)
    character(2)  :: b1, b2
    
    if (present(order)) then
       ordc(1) = order(2)
       ordc(2) = order(1)
    else
       ordc = [2,1]
    endif
    
    print*
    do j = 1,3
        ic(ordc(2)) = j
        do i = 1,3
            ic(ordc(1)) = i
            write (ctemp(i),'(f12.7)') tensor(ic(1),ic(2))
        enddo
        ! Printing
        call brackets(j,b1,b2)
        print*, '   '//b1//ctemp(1)//ctemp(2)//ctemp(3)//b2
    enddo
    
end subroutine

subroutine print_order1(tensor)
    real(dp), intent(in) :: tensor(3)
    
    integer       :: j
    character(12) :: ctemp
    character(2)  :: b1,b2
    
    print*
    do j = 1,3
        write (ctemp,'(f12.7)') tensor(j)
        call brackets(j,b1,b2)
        print*, '   '//b1//' '//trim(adjustl(ctemp))//b2
    enddo
end subroutine

subroutine brackets(j,b1,b2)
    integer,intent(in) :: j
    character(2), intent(out) :: b1,b2
    if (j==1)then !Matrix brackets
      b1=" /"; b2=" \" !"
    elseif (j==2)then
      b1=" |"; b2=" |" !"
    else
      b1=" \"; b2=" /" !"
    endif
end subroutine


subroutine print_integer_matrix(matrix,space,copy)
  integer, optional :: space, copy
  integer :: matrix(:,:),imax, jmax, i, j, xspace, maxx(size(matrix,2))
  character(5) :: endsep 
  character(1) :: sep
  
  imax = size(matrix,1)
  jmax = size(matrix,2)
  
  
  xspace = 5
  if(present(space)) xspace = space
  
  if(xspace==0)then
    do j = 1,jmax
      maxx(j) = maxval(matrix(:,j))
      
    enddo
    maxx = ceiling(log10(dble(maxx+1)))
  else
    maxx = xspace
  endif  
  
  
  sep = " "
  endsep = " "
  if(present(copy)) then 
    if (copy>0)then
      sep = ","
      endsep = " & ! "
    endif
  endif
  
    
  
    
  print*, "____"
  print*, ">>>>   Integer Matrix: "//i2s(imax)//" x "//i2s(jmax)
  
  
  
  do i = 1,imax 
    do j = 1, jmax
      write(*,'(I'//i2s(maxx(j))//',a)',advance="no") matrix(i,j),sep 
    enddo
    write(*,'(a)') endsep
  enddo
    
    
    
  
end subroutine 



!////////////////////////////////////////////////////// Special Printers:


subroutine xyz_hho_to_linear(a,alin,nM)
  integer, intent(in) :: nM
  real(dp), intent(in) :: a(3,3,nM)!xyz,hho,nM
  real(dp), intent(out) :: alin(nM*9)
  integer m, xyz
  alin=0
   do m = 1,nM
     do xyz = 1,3
       alin((2*m-2)*3 + xyz)     = a(xyz,1,m) !a(m)%h1(xyz) 
       alin((2*m-1)*3 + xyz)     = a(xyz,2,m) !a(m)%h2(xyz) 
       alin((2*nM+m-1)*3 + xyz)  = a(xyz,3,m) !a(m)%o(xyz)  
     enddo
   enddo
end subroutine 





!//////////////////////////////////////////////////// Tensor Printers for "printer" interface above:
! p_ means printer
! printer_text is at the bottom
! printers for: tensors of rank 1-5, type(h2o) objects, 

!s=0:noprint
!s=1:print
!s=2:print transpose


!subroutine p_h2o(a,text,s)
!  character(*) text
!  type(h2o) :: a
!  integer s
!  if(s==0)return
!   call printer_text(text)
!   write(un,'(a,3'//dblf//')') '   H1:',a%h1(:)
!   write(un,'(a,3'//dblf//')') '   H2:',a%h2(:)
!   write(un,'(a,3'//dblf//')') '   O :',a%o(:)
!end subroutine 
!
!subroutine p_h2os(a,text)
!  character(*) text
!  type(h2o) :: a(:)
!  integer l, m
!   call printer_text(text)
!   l = size(a)
!   do m = 1,l
!     write(un,'(a,I3)') '   H2O Nr:',m
!     write(un,'(a,3'//dblf//')') '   H1:',a(m)%h1(:)
!     write(un,'(a,3'//dblf//')') '   H2:',a(m)%h2(:)
!     write(un,'(a,3'//dblf//')') '   O :',a(m)%o(:)
!   enddo
!end subroutine 
   
     
  
!s=0:noprint,,,,,,, this was stupid!!!
!s=1:print
!s=2:print transpose

!debug = 1 => print
!debug = 0 => noprint




subroutine p_int(a,text,s)!,production)
  character(*) text
  integer a
  integer s
  !logical production
  !if(production)return
   call printer_text_scalar(text)
   write(un,intf) a
end subroutine

subroutine p_real(a,text,s)!,production)
  character(*) text
  real(dp) a
  integer s
  !logical production
  !if(production)return
   call printer_text_scalar(text)
   write(un,'('//dblf//')') a
end subroutine

subroutine p_real_vec(a,text,s)!,production)
  character(*) text
  character(3) :: ch_cols
  character(20):: forma
  real(dp) a(:)
  integer length, i
  integer s
  !logical production
  !if(production)return

  length = size(a)
  call printer_text(text)
  
  if(s==1)then
    do i = 1,length
      write(un,'('//dblf//')') a(i)
    enddo
  elseif(s==2)then
    write(ch_cols,'(I3)') length
    forma = '('//ch_cols//dblf//')'
    write(un,forma) a(:)
  endif
end subroutine


subroutine p_real_mat(a,text,s)!,production)
  character(*) :: text
  character(3) :: ch_cols
  character(20):: forma
  integer rows,cols,i
  real(dp)     :: a(:,:)
  integer s
  !logical production
  !if(production)return

   rows = size(a,1)
   cols = size(a,2)

   write(ch_cols,'(I3)') cols
   forma = '('//ch_cols//dblf//')'

   call printer_text(text)
   if(s==1)then
     do i = 1,rows
       write(un,forma) a(i,:)
     enddo
   elseif(s==2)then
     do i = 1,cols
       write(un,forma) a(:,i)
     enddo
   endif
end subroutine

subroutine p_real_3d(a,text,s)!,production)
  character(*) :: text
  character(3) :: ch_cols
  character(20):: forma
  integer rows,cols,slices,slice,i
  real(dp)     :: a(:,:,:)
  integer s
  !logical production
  !if(production)return

   rows = size(a,1)
   cols = size(a,2)
   slices = size(a,3)

   write(ch_cols,'(I3)') cols
   forma = '('//ch_cols//dblf//')'

   call printer_text(text)
   do slice = 1,slices
      write(un,'(a,I2,a)') "   slice nr. ",slice, ':'
      if(s==1)then
        do i = 1,rows
          write(un,forma) a(i,:,slice)
        enddo
      elseif(s==2)then
        do i = 1,cols
          write(un,forma) a(:,i,slice)
        enddo
      endif
   enddo
end subroutine

subroutine p_real_4d(a,text,s)!,production)
  character(*) :: text
  character(3) :: ch_cols
  character(20):: forma
  integer rows,cols,slices,slice,i, i4s, i4
  real(dp)     :: a(:,:,:,:)
  integer s
  !logical production
  !if(production)return

   rows = size(a,1)
   cols = size(a,2)
   slices = size(a,3)
   i4s = size(a,4)

   write(ch_cols,'(I3)') cols
   forma = '('//ch_cols//dblf//')'

   call printer_text(text)
   do i4 = 1,i4s
      write(un,'(a,I2,a)') "   4D index nr. ",i4, ':'
      do slice = 1,slices
         write(un,'(a,I2,a)') "      Slice nr. ",slice, ':'
         if(s==1)then
         do i = 1,rows
           write(un,forma) a(i,:,slice,i4 )
         enddo
         elseif(s==2)then
         do i = 1,cols
           write(un,forma) a(:,i,slice,i4 )
         enddo
         endif
      enddo
   enddo
end subroutine

subroutine p_real_5d(a,text,s)!,production)
  character(*) :: text
  character(3) :: ch_cols
  character(20):: forma
  integer rows,cols,slices,slice,i, i4s, i4, i5s, i5
  real(dp)     :: a(:,:,:,:,:)
  integer s
  !logical production
  !if(production)return

   rows = size(a,1)
   cols = size(a,2)
   slices = size(a,3)
   i4s = size(a,4)
   i5s = size(a,5)

   write(ch_cols,'(I3)') cols
   forma = '('//ch_cols//dblf//')'

   call printer_text(text)
   do i5 = 1,i5s
      write(un,'(a,I2,a)') "   5D index nr. ",i5, ':'
      do i4 = 1,i4s
         write(un,'(a,I2,a)') "      4D index nr. ",i4, ':'
         do slice = 1,slices
            write(un,'(a,I2,a)') "         Slice nr. ",slice, ':'
            if(s==1)then
            do i = 1,rows
              write(un,forma) a(i,:,slice,i4,i5)
            enddo
            elseif(s==2)then
            do i = 1,cols
              write(un,forma) a(:,i,slice,i4,i5)
            enddo
            endif
         enddo
      enddo
   enddo
end subroutine



! Internal ///////////////////
subroutine printer_text(text)
  character(*) text
   write(un,'(A)') '/'//text//': ------ '
end subroutine

subroutine printer_text_scalar(text)
  character(*) text
   write(un,'(A)',advance='no') '/'//text//':  '
end subroutine


end module


!!//////////////////////////////////////////////////////////
!program test
!use printer_mod
!use data_types
!implicit none
!real(dp) :: a(3,3), b(3,3,3), c(3,3,3,3), d(3,3,3,3,3)
!a(:,1) = [1.0,2.0,3.0]
!a(:,2) = [10.0,2.0,3.0]
!a(:,3) = [100.0,2.0,3.0]
!b(:,:,1) = a
!b(:,:,2) = a+1
!b(:,:,3) = a+2
!c(:,:,:,1) = b
!c(:,:,:,2) = b+1
!c(:,:,:,3) = b+2
!d(:,:,:,:,1) = c
!d(:,:,:,:,2) = c+1
!d(:,:,:,:,3) = c+2
!
!
!
!
!call printer(3.0_dp,'hej')
!call printer(3,'hej')
!call printer([3.0_dp,3.0_dp,3.0_dp],'hej')
!call printer(a,'hej')
!call printer(b,'hej,detta är b')
!call printer(c,'hej,detta är c   c')
!call printer(d,'hej,detta är d   d     d  d  d')
!
!end program
!
