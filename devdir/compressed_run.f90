module compressed_run

use compressed_tests, bad => main

implicit none


contains !//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

subroutine main
call &
!test_intfac_ff
test_df

!test_polynextpow_n
!test_mp_pot
!test_polyinner
!print_long_index_matrix(15)
!print_square_index_matrix(15)
!print_trace_index_matrix(15)
!print_trace_ind
!print_apple_g(15)
!test_detracers
!test_detracer_linearity
!test_polydet
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
