program tqws
      
      use iso_c_binding
      implicit none
      
      interface 
         subroutine potnasa(xxx,dvdr,vvv) bind(C,name='potnasa2_')
            use iso_c_binding
            real(c_double), dimension(0:8) :: xxx,dvdr
            real(c_double) :: vvv
         end subroutine
      end interface 
      
      integer, parameter :: n = 2
      real(c_double) :: mol(n,0:8), mol2(0:8), dpot(0:8), vv
      integer i, j
      real*8 , parameter :: A2b = 1.889725989d0 !Ångström to Bohr
      real*8 , parameter :: au2deb = 2.541746230211d0 !a.u. to debye
      !(c_double)
      !      x(n,xyz,HHO)
      mol(1,0:2)  = [  0.000d0,  0.000d0,  0.000d0]!*A2b ! O
      mol(1,3:5)  = [  0.018d0, -0.739d0,  0.521d0]!*A2b ! H
      mol(1,6:8)  = [ -0.815d0, -0.673d0, -0.592d0]!*A2b ! H 
      mol(2,0:2)  = [  1.675d0, -2.803d0,  0.911d0]!*A2b ! O
      mol(2,3:5)  = [  1.942d0, -1.964d0,  1.223d0]!*A2b ! H
      mol(2,6:8)  = [  2.319d0, -3.450d0,  1.136d0]!*A2b ! H
      
      do i = 1,n
          dpot=0;vv=0
          call potnasa(mol(i,:),dpot,vv)
          do j=0,8,3
             print*, dpot(j), dpot(j+1), dpot(j+2)
          enddo
          print*,'-------------'
          print*, 'pot:',vv
          print*,'-------------'
      enddo
      
end program 
