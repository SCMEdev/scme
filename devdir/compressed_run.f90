module compressed_run

use compressed_tests, bad => main

implicit none


contains !//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

subroutine main
call &
test_polydet
!test_polyfinder
!test_detracer
!all_tests
!get_traces !-----------
!call &
!test_brakk
!print_traces(5)
!print_trace_keys
!test_intfac_ff
end subroutine



end module