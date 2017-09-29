!!!!#define call test(A)

program slkdjf
use detrace_mod, only:test
#ifdef RUN_TEST
call test()
#endif
#ifdef RUN_OTHER 
print*, "Other"
#endif
end program
