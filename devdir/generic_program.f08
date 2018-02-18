!!!!#define call test(A)

program slkdjf

!integer, parameter :: dp = kind(0d0)


#ifdef MODU
use MODU

#ifdef SUB
call SUB
#else
call main
#endif

#endif


#ifdef STUPID
print*, size([1.0,2.0]), size([1.0]) !, size(5.0)
#endif



end program
