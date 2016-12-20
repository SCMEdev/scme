    program testdip
      use ps_dms, only: dmsnasa
      implicit none
      integer, parameter :: n = 2
      real*8 x(n,3,3), d(3)!, A2b, au2deb
      integer i,j
      real*8 , parameter :: A2b = 1.889725989d0 !Ångström to Bohr
      real*8 , parameter :: au2deb = 2.541746230211d0 !a.u. to debye
      !      x(n,xyz,HHO)

      x(1,1,:) = [  0.000d0,  0.000d0,  0.000d0]!*A2b ! O
      x(1,2,:) = [  0.018d0, -0.739d0,  0.521d0]!*A2b ! H
      x(1,3,:) = [ -0.815d0, -0.673d0, -0.592d0]!*A2b ! H 
      x(2,1,:) = [  1.675d0, -2.803d0,  0.911d0]!*A2b ! O
      x(2,2,:) = [  1.942d0, -1.964d0,  1.223d0]!*A2b ! H
      x(2,3,:) = [  2.319d0, -3.450d0,  1.136d0]!*A2b ! H


      do i = 1,n
          call dmsnasa(x(i,:,:),d)
          print*, sqrt(sum(d(:)**2))*au2deb
          print*, d(:)*au2deb
      enddo
    end program 
