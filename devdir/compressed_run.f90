module compressed_run

use compressed_tests, bad => main

implicit none


contains !//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

subroutine main
call &
!test_next_lex_and_can
test_polarize
!test_polarize2
!test_compress_expand_subdivided
!test_next_set2pown(5,1)
!test_stone_field
!call &
!test_old_field
!geom2xyz
!test_intfac_ff
!call &
!test_df
!test_polarize2
!call &
!test_printoa
!test_polynextpow_n
!test_mp_pot
!test_polyinner
!print_long_index_matrix(15)
!call &
!print_square_index_matrix(15)
!print_trace_index_matrix(15)
!print_trace_ind
!print_apple_g(15)
!call &
!test_detracers
!test_detracer_linearity
!test_polydet
!test_polyfinder
!call &
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
