!!!!#define call test(A)

program slkdjf

!integer, parameter :: dp = kind(0d0)

#ifdef HEXAD 
use detrace_mod !, only:test
integer i
integer, parameter :: le=15
real(dp) tricorn(le)
do i = 1,le
  read(*,*) tricorn(i)
enddo
call detrace_hexadeca(tricorn,.false.)
do i = 1,le
  print*, tricorn(i)
enddo
#endif

#ifdef OCTA 
use detrace_mod !, only:test
integer i
integer, parameter :: le=10
real(dp) tricorn(le)
do i = 1,le
  read(*,*) tricorn(i)
enddo
call detrace_octa(tricorn,.false.)
do i = 1,le
  print*, tricorn(i)
enddo
#endif

#ifdef GENERAL
use detrace_mod !, only:test
integer i
integer, parameter :: maxl=16
real(dp) tricorn(maxl)
integer le, vari
do i = 1,maxl
  read(*,*,iostat=vari) tricorn(i)
  if(vari<0)then
    le = i-1
    exit
  endif
enddo
call detrace(tricorn(1:le))
do i = 1,le
  print*, tricorn(i)
enddo
#endif

#ifdef APPLE
use detrace_apple, only:xtrace !,only:qpole_main => main
integer, parameter :: dp = kind(0d0)
integer i
!character(4) dum
integer, parameter :: le=15
real(dp) inten(le),outen(le) 
!print*, "hej"
do i = 1,le
  read(*,*) inten(i)
enddo
call xtrace(inten,outen,4)
do i = 1,le
  print*, outen(i) /7d0/5d0/3d0
enddo
#endif



#ifdef YO
use YO
call main
#endif


#ifdef STUPID
print*, size([1.0,2.0]), size([1.0]) !, size(5.0)
#endif



end program
