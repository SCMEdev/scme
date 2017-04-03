program test_ps_cpp 
      use iso_c_binding
      implicit none
      
      interface 
         subroutine dmsnasa(xxx,ddd) bind(C,name='dmsnasa2_')
            use iso_c_binding
            real(c_double), dimension(9) :: xxx
            real(c_double), dimension(3) :: ddd
         end subroutine
      end interface 
      
      integer, parameter :: n = 2
      real*8 mol(n,0:8), qdms(3)!, A2b, au2deb
      integer i
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
         qdms = 0
         call dmsnasa(mol(i,:),qdms)
         print*, sqrt(sum(qdms**2)) *au2deb !'cpp total dipole     :',
         print*,               qdms *au2deb !'cpp dipole components:',
      enddo
      
end program 
