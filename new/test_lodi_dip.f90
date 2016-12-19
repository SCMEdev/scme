    program testlodi
      implicit none
      integer, parameter :: n = 2
      real*8 x(n,3,3), d(n,2)
      integer i
      real*8 , parameter :: A2b = 1.889725989d0 !Ångström to Bohr
      real*8 , parameter :: au2deb = 2.541746230211d0 !a.u. to debye
      real*8 v1(3), v2(3), r1,r2,cabc,d1,d2
      
      
      
      x(1,1,:) = [  0.000,  0.000,  0.000]*A2b ! O
      x(1,2,:) = [  0.018, -0.739,  0.521]*A2b ! H
      x(1,3,:) = [ -0.815, -0.673, -0.592]*A2b ! H 
      x(2,1,:) = [  1.675, -2.803,  0.911]*A2b ! O
      x(2,2,:) = [  1.942, -1.964,  1.223]*A2b ! H
      x(2,3,:) = [  2.319, -3.450,  1.136]*A2b ! H
      
      ! x is in bohr
      
      do i = 1,2
          v1 = x(i,3,:) - x(i,1,:) !JÖ changed order
          r1 = sqrt(sum( v1**2 ))
          
          v2 = x(i,2,:) - x(i,1,:) !JÖ changed order
          r2 = sqrt(sum( v2**2 ))
          
          cabc = ( v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3) ) / (r1*r2)      
          
          call DIPS(d1,d2,r1,r2,cabc)
          d(i,:) = [d1,d2]
      enddo



      
      do i = 1,n
          print*, sum(d(i,:)**2)*au2deb
          print*, d(i,:)*au2deb
      enddo
    end program 

