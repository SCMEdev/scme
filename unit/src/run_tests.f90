! Copyright (C)  2015-2016  SCMEdev
! Licenced under LGPLv3. See LICENCE for details.


! Include the macros.
#include "mifu.h"

! The main routine for the test suite.
program testrunner

  ! Use the test module(s) here.
  use test_scme

  ! Register each test here.

  ! Basic functionality tests for scme.
  MIFU_REGISTER_TEST(test_scme_monomer1)
  MIFU_REGISTER_TEST(test_scme_dimer1)
  MIFU_REGISTER_TEST(test_scme_perf)
!  MIFU_REGISTER_TEST(test_scme_cluster_perf)

  ! End the tests.
  MIFU_END()

end program testrunner

