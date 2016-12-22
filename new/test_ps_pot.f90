!-----------------------------------------------------------------------
! This is text
! Under this text is source code
!-----------------------------------------------------------------------


!newer:
program ps_test
      use ps_pot, only: vibpes
      implicit none
      integer, parameter :: n = 2
      real*8 pot, x(n,3,3),dpot(9)!n,OHH,xyz
      !real*8 A2b!, pi, deg2rad, conv(3)
      integer i, j
      !pi = 3.1415926535d0
      real*8, parameter :: A2b = 1.0d0/0.529177249d0 !Angstr to Bohr
      
      !conversion to bohr happens here, in the code, bohr should come in too
      
      x(1,2,:) = [  0.018d0, -0.739d0,  0.521d0]!*A2b ! H
      x(1,3,:) = [ -0.815d0, -0.673d0, -0.592d0]!*A2b ! H 
      x(2,2,:) = [  1.942d0, -1.964d0,  1.223d0]!*A2b ! H
      x(2,3,:) = [  2.319d0, -3.450d0,  1.136d0]!*A2b ! H
      x(1,1,:) = [  0.000d0,  0.000d0,  0.000d0]!*A2b ! O
      x(2,1,:) = [  1.675d0, -2.803d0,  0.911d0]!*A2b ! O

      
      do i = 1,n
          
          call vibpot(x(i,:,:),pot,dpot)
          
          
          do j = 1,9,3
              print*, dpot(j), dpot(j+1), dpot(j+2)
          enddo
          print*, '---------------'
          print*, 'pot:', pot
          print*, '---------------'
          
      enddo
      
end program












!program ps_test
!      use ps_pot, only: vibpot
!      implicit none
!      integer, parameter :: n = 6
!      real*8 pot(n), x(n,3,3)!n,OHH,xyz
!      !real*8 A2b!, pi, deg2rad, conv(3)
!      integer p, i
!      !pi = 3.1415926535d0
!      real*8, parameter :: A2b = 1.0d0/0.529177249d0 !Angstr to Bohr
!      !deg2rad = pi/180.d0
!      !conv = [A2au, A2au, deg2rad]
!      !p=1;   conf(p,:) = [0.95d0,1.01d0,106.0d0]*conv
!      !p=p+1; conf(p,:) = [0.93d0,0.93d0,104.0d0]*conv
!      !p=p+1; conf(p,:) = [0.94d0,0.94d0,105.0d0]*conv
!      !p=p+1; conf(p,:) = [1.34d0,1.24d0,108.0d0]*conv
!      
!      !conversion to bohr happens here, in the code, bohr should come in too
!      
!      x(1,1,:) = [-1.502169d+00, -1.913590d-01,  1.434927d+00]*A2b ! O
!      x(1,2,:) = [-2.006698d+00, -4.223270d-01,  2.219847d+00]*A2b ! H
!      x(1,3,:) = [-6.010540d-01, -5.969720d-01,  1.553718d+00]*A2b ! H 
!      x(2,1,:) = [-1.744575d+00, -3.823480d-01, -1.309144d+00]*A2b ! O
!      x(2,2,:) = [-2.516835d+00, -7.667650d-01, -1.733766d+00]*A2b ! H
!      x(2,3,:) = [-1.888941d+00, -4.796530d-01, -3.476240d-01]*A2b ! H
!      x(3,1,:) = [-5.604090d-01,  2.017830d+00, -1.219840d-01]*A2b ! O
!      x(3,2,:) = [-9.898310d-01,  1.592736d+00, -8.774190d-01]*A2b ! H
!      x(3,3,:) = [-9.477200d-01,  1.533567d+00,  6.252280d-01]*A2b ! H
!      x(4,1,:) = [ 9.648030d-01, -1.165765d+00,  1.439987d+00]*A2b ! O
!      x(4,2,:) = [ 1.542224d+00, -3.936920d-01,  1.344373d+00]*A2b ! H
!      x(4,3,:) = [ 9.795570d-01, -1.522041d+00,  5.278330d-01]*A2b ! H
!      x(5,1,:) = [ 9.747050d-01, -1.401503d+00, -1.335970d+00]*A2b ! O
!      x(5,2,:) = [ 1.470709d+00, -5.709330d-01, -1.277710d+00]*A2b ! H
!      x(5,3,:) = [ 6.516100d-02, -1.118951d+00, -1.522886d+00]*A2b ! H
!      x(6,1,:) = [ 2.002280d+00,  1.057824d+00, -1.245020d-01]*A2b ! O
!      x(6,2,:) = [ 2.674716d+00,  1.735342d+00, -2.379950d-01]*A2b ! H
!      x(6,3,:) = [ 1.141637d+00,  1.532266d+00, -1.401210d-01]*A2b ! H      
!      
!      call vibpot(x,pot,n)
!      do i = 1,n
!          print*, 'pot = ',pot(i), "a.u. = ",pot(i)*27.211396132, "eV"
!      enddo
!end program


!    coordinates(1)  =  0.018       ! x H1 on W1
!    coordinates(2)  = -0.739       ! y H1 on W1
!    coordinates(3)  =  0.521       ! z H1 on W1
!    
!    coordinates(4)  = -0.815       ! x H2 on W1
!    coordinates(5)  = -0.673
!    coordinates(6)  = -0.592!< done
!    
!    coordinates(7)  =  1.942       ! x H1 on W2
!    coordinates(8)  = -1.964       ! y H1 on W2
!    coordinates(9)  =  1.223       ! z H1 on W2!< done
!    
!    coordinates(10) =  2.319       ! x H2 on W2
!    coordinates(11) = -3.450
!    coordinates(12) =  1.136!< done
!    
!    coordinates(13) =  0.000       ! x O on W1
!    coordinates(14) =  0.000
!    coordinates(15) =  0.000!< done
!    
!    coordinates(16) =  1.675       ! x O on W2
!    coordinates(17) = -2.803
!    coordinates(18) =  0.911

      !deg2rad = pi/180.d0
      !conv = [A2au, A2au, deg2rad]
      !p=1;   conf(p,:) = [0.95d0,1.01d0,106.0d0]*conv
      !p=p+1; conf(p,:) = [0.93d0,0.93d0,104.0d0]*conv
      !p=p+1; conf(p,:) = [0.94d0,0.94d0,105.0d0]*conv
      !p=p+1; conf(p,:) = [1.34d0,1.24d0,108.0d0]*conv




!    program testdip
!      use ps_dip, only: vibdip
!      implicit none
!      integer, parameter :: n = 6
!      real*8 x(n,3,3), d(n,3), A2b, au2deb
!      integer i
!      !      x(n,xyz,HHO)
!      A2b = 1.889725989d0 !Ångström to Bohr
!      au2deb = 2.541746230211d0 !a.u. to debye
!      
!
!
!      x(1,:,1) = [-1.502169d+00, -1.913590d-01,  1.434927d+00]*A2b ! O
!      x(1,:,2) = [-2.006698d+00, -4.223270d-01,  2.219847d+00]*A2b ! H
!      x(1,:,3) = [-6.010540d-01, -5.969720d-01,  1.553718d+00]*A2b ! H 
!      x(2,:,1) = [-1.744575d+00, -3.823480d-01, -1.309144d+00]*A2b ! O
!      x(2,:,2) = [-2.516835d+00, -7.667650d-01, -1.733766d+00]*A2b ! H
!      x(2,:,3) = [-1.888941d+00, -4.796530d-01, -3.476240d-01]*A2b ! H
!      x(3,:,1) = [-5.604090d-01,  2.017830d+00, -1.219840d-01]*A2b ! O
!      x(3,:,2) = [-9.898310d-01,  1.592736d+00, -8.774190d-01]*A2b ! H
!      x(3,:,3) = [-9.477200d-01,  1.533567d+00,  6.252280d-01]*A2b ! H
!      x(4,:,1) = [ 9.648030d-01, -1.165765d+00,  1.439987d+00]*A2b ! O
!      x(4,:,2) = [ 1.542224d+00, -3.936920d-01,  1.344373d+00]*A2b ! H
!      x(4,:,3) = [ 9.795570d-01, -1.522041d+00,  5.278330d-01]*A2b ! H
!      x(5,:,1) = [ 9.747050d-01, -1.401503d+00, -1.335970d+00]*A2b ! O
!      x(5,:,2) = [ 1.470709d+00, -5.709330d-01, -1.277710d+00]*A2b ! H
!      x(5,:,3) = [ 6.516100d-02, -1.118951d+00, -1.522886d+00]*A2b ! H
!      x(6,:,1) = [ 2.002280d+00,  1.057824d+00, -1.245020d-01]*A2b ! O
!      x(6,:,2) = [ 2.674716d+00,  1.735342d+00, -2.379950d-01]*A2b ! H
!      x(6,:,3) = [ 1.141637d+00,  1.532266d+00, -1.401210d-01]*A2b ! H
!
!
!      
!      call vibdip(x,d,n)
!      
!      do i = 1,n
!          print*, sum(d(i,:)**2)*au2deb
!      enddo
!    end program        
!       
!
       
       
       
!solve this:        
       !c5z(:) = ( f5z*c5z(:) + fbasis*cbasis(:) + fcore*ccore(:) + frest*crest(:) )*cm2au 
       !c5z(:) = [tc5z(1)*2d0,tc5z(2:245)]
       
!      print*, pot
       !phh2    = phh2*0.529177249d0
!       c5z(:) = c5z(:)
!       b1      = 2.0d0
!         roh=roh/0.529177249d0
!         alphaoh=alphaoh*0.529177249d0
!       phh1    = 16.94879431193463d0
!         phh1=phh1*f5z
!         phh1=phh1*4.556335d-6
!       deoh    = 42290.92019288289d0
!         deoh=deoh*f5z
!         deoh=deoh*4.556335d-6
!       reoh    = 0.958649d0
         !do  i=1,245
         !enddo 
!      end if

