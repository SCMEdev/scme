program test_ps_cpp
      
      implicit none
      !integer, parameter :: n = 2
      real*8 mol(9), mol2(9), qdms(3)!, A2b, au2deb
      integer i
      real*8 , parameter :: A2b = 1.889725989d0 !Ångström to Bohr
      real*8 , parameter :: au2deb = 2.541746230211d0 !a.u. to debye
      
      !      x(n,xyz,HHO)
      mol(1:3)  = [  0.000,  0.000,  0.000]!*A2b ! O
      mol(4:6)  = [  0.018, -0.739,  0.521]!*A2b ! H
      mol(7:9)  = [ -0.815, -0.673, -0.592]!*A2b ! H 
      
      mol2(1:3) = [  1.675, -2.803,  0.911]!*A2b ! O
      mol2(4:6) = [  1.942, -1.964,  1.223]!*A2b ! H
      mol2(7:9) = [  2.319, -3.450,  1.136]!*A2b ! H

      qdms = 0
      call dmsnasa2(mol,qdms)
      print*, 'cpp total dipole     :',sqrt(sum(qdms**2))/0.393456
      print*, 'cpp dipole components:',qdms/0.393456
      
      qdms = 0
      call dmsnasa2(mol2,qdms)
      print*, 'cpp total dipole     :',sqrt(sum(qdms**2))/0.393456
      print*, 'cpp dipole components:',qdms/0.393456



!      call vibdip(x,d,n)
!      do i = 1,n
!          print*, sum(d(i,:)**2)*au2deb
!          print*, d(i,:)*au2deb
!      enddo
end program 

!x(1,1,:) = 
!x(1,2,:) = 
!x(1,3,:) = 
!
!x(2,1,:) = 
!x(2,2,:) = 
!x(2,3,:) = 
