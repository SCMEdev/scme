! This module print tensors of rank 1 to 5 with the generic interface "printer of all the s ranks. 
module printer_mod

  use data_types, only: dp, h2o
  
  implicit none
  integer, parameter :: un=6
  character(*), parameter :: intf='(I7)'
  character(*), parameter :: dblf='f24.15' !double format
!  character(*), parameter :: dblf='(f24.16)'
  
interface printer
 module procedure p_real, p_int, p_real_vec, p_real_mat, p_real_3d, p_real_4d, p_real_5d, p_h2o, p_h2os
end interface

  private
  public printer, printer_h2o_linear


contains !/////////////////////////////////////////////

!////////////////////////////////////////////////////// Special Printers:

! Prints h2o type in fa order. 
! It takes type(h2o) and outputs a vector: [H1x H1y H1z H2x H2y H2z H1x H1y H1z H2x H2y H2z ... Ox Oy Oz ...]
subroutine printer_h2o_linear(a,text)
  character(*) text
  type(h2o) :: a(:)
  integer l, m, xyz
   call printer_text(text)
   l = size(a)
   
   do m = 1,l
     do xyz = 1,3
       write(un,'(a,I3,a,3'//dblf//')') '   Mol:',m,' H1:',a(m)%h1(xyz)
       write(un,'(a,I3,a,3'//dblf//')') '   Mol:',m,' H2:',a(m)%h2(xyz)
     enddo
   enddo
   do m = 1,l
     do xyz = 1,3
       write(un,'(a,I3,a,3'//dblf//')') '   Mol:',m,' O :',a(m)%o(xyz)
     enddo
   enddo

end subroutine 
 


!//////////////////////////////////////////////////// Tensor Printers for "printer" interface above:
! p_ means printer
! printer_text is at the bottom
! printers for: tensors of rank 1-5, type(h2o) objects, 


subroutine p_h2o(a,text)
  character(*) text
  type(h2o) :: a
   call printer_text(text)
   write(un,'(a,3'//dblf//')') '   H1:',a%h1(:)
   write(un,'(a,3'//dblf//')') '   H2:',a%h2(:)
   write(un,'(a,3'//dblf//')') '   O :',a%o(:)
end subroutine 

subroutine p_h2os(a,text)
  character(*) text
  type(h2o) :: a(:)
  integer l, m
   call printer_text(text)
   l = size(a)
   do m = 1,l
     write(un,'(a,I3)') '   H2O Nr:',m
     write(un,'(a,3'//dblf//')') '   H1:',a(m)%h1(:)
     write(un,'(a,3'//dblf//')') '   H2:',a(m)%h2(:)
     write(un,'(a,3'//dblf//')') '   O :',a(m)%o(:)
   enddo
end subroutine 
   
     
  

subroutine p_int(a,text)
  character(*) text
  integer a
   call printer_text(text)
   write(un,intf) a
end subroutine

subroutine p_real(a,text)
  character(*) text
  real(dp) a
   call printer_text(text)
   write(un,'('//dblf//')') a
end subroutine

subroutine p_real_vec(a,text)
  character(*) text
  real(dp) a(:)
  integer length, i
  length = size(a)
  call printer_text(text)
  
  do i = 1,length
    write(un,'('//dblf//')') a(i)
  enddo
  
end subroutine


subroutine p_real_mat(a,text)
  character(*) :: text
  character(3) :: ch_cols
  character(20):: forma
  integer rows,cols,i
  real(dp)     :: a(:,:)

   rows = size(a,1)
   cols = size(a,2)

   write(ch_cols,'(I3)') cols
   forma = '('//ch_cols//dblf//')'

   call printer_text(text)
   do i = 1,rows
     write(un,forma) a(i,:)
   enddo
end subroutine

subroutine p_real_3d(a,text)
  character(*) :: text
  character(3) :: ch_cols, ch_slices
  character(20):: forma
  integer rows,cols,slices,slice,i
  real(dp)     :: a(:,:,:)

   rows = size(a,1)
   cols = size(a,2)
   slices = size(a,3)

   write(ch_cols,'(I3)') cols
   forma = '('//ch_cols//dblf//')'

   call printer_text(text)
   do slice = 1,slices
      write(un,'(a,I2,a)') "   slice nr. ",slice, ':'
      do i = 1,rows
        write(un,forma) a(i,:,slice)
      enddo
   enddo
end subroutine

subroutine p_real_4d(a,text)
  character(*) :: text
  character(3) :: ch_cols
  character(20):: forma
  integer rows,cols,slices,slice,i, i4s, i4
  real(dp)     :: a(:,:,:,:)

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
         do i = 1,rows
           write(un,forma) a(i,:,slice,i4 )
         enddo
      enddo
   enddo
end subroutine

subroutine p_real_5d(a,text)
  character(*) :: text
  character(3) :: ch_cols
  character(20):: forma
  integer rows,cols,slices,slice,i, i4s, i4, i5s, i5
  real(dp)     :: a(:,:,:,:,:)

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
            do i = 1,rows
              write(un,forma) a(i,:,slice,i4,i5)
            enddo
         enddo
      enddo
   enddo
end subroutine



! Internal ///////////////////
subroutine printer_text(text)
  character(*) text
   write(un,'(A)') '/'//text//': ------ '
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
