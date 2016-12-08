!!!JÖ  ! Copyright (C)  2015-2016  SCMEdev
!!!JÖ  ! Licenced under LGPLv3. See LICENCE for details.
!!!JÖ  
!!!JÖ  
!!!JÖ  ! Include the macros.
!!!JÖ  #include "mifu.h"
!!!JÖ  
!!!JÖ  ! The main routine for the test suite.
!!!JÖ  program testrunner
!!!JÖ  
!!!JÖ    ! Use the test module(s) here.
!!!JÖ    use test_scme
!!!JÖ  
!!!JÖ    ! Register each test here.
!!!JÖ  
!!!JÖ    ! Basic functionality tests for scme.
!!!JÖ    MIFU_REGISTER_TEST(test_scme_monomer1)
!!!JÖ    MIFU_REGISTER_TEST(test_scme_dimer1)
!!!JÖ  !  MIFU_REGISTER_TEST(test_scme_perf)
!!!JÖ    MIFU_REGISTER_TEST(test_scme_cluster_perf)
!!!JÖ  
!!!JÖ    ! End the tests.
!!!JÖ    MIFU_END()
!!!JÖ  
!!!JÖ  end program testrunner

program sausage
use test_scme
call test_scme_cluster_32_perf()
end program
