Flat profile:

Each sample counts as 0.01 seconds.
  %   cumulative   self              self     total           
 time   seconds   seconds    calls  ms/call  ms/call  name    
 58.59      1.33     1.33     5060     0.26     0.27  __calc_lower_order_MOD_calcedip_quad
 13.66      1.64     0.31                             __tang_toennies_MOD_tang_toenniesn
 10.57      1.88     0.24    59400     0.00     0.03  __scme_MOD_scme_calculate
  8.37      2.07     0.19       20     9.50     9.50  __calc_derivs_MOD_calcdv
  2.64      2.13     0.06     5060     0.01     0.01  __molecproperties_MOD_sfdsf
  1.76      2.17     0.04       20     2.00     2.00  __inducepoles_MOD_inducedipole
  1.32      2.20     0.03  5029200     0.00     0.00  __mdutil_MOD_cross
  0.88      2.22     0.02  2235200     0.00     0.00  __tang_toennies_MOD_tang_toenniesndf
  0.88      2.24     0.02       20     1.00     1.01  __calc_higher_order_MOD_calcehigh
  0.44      2.25     0.01     5060     0.00     0.00  __inducepoles_MOD_induceqpole
  0.44      2.26     0.01       20     0.50     0.50  __molecproperties_MOD_addfields
  0.44      2.27     0.01       20     0.50     0.50  __molecproperties_MOD_rotatepolariz
  0.00      2.27     0.00    19800     0.00     0.00  __mdutil_MOD_inv6
  0.00      2.27     0.00     8800     0.00     0.00  ps::dmsnasa(double const*, double*)
  0.00      2.27     0.00     5060     0.00     0.00  __calcenergy_mod_MOD_calcenergy
  0.00      2.27     0.00     5060     0.00     0.00  __molecproperties_MOD_adddfields
  0.00      2.27     0.00      640     0.00     0.00  __molforce_MOD_molforce3
  0.00      2.27     0.00       20     0.00     0.00  __forcecm_mod_MOD_forcecm
  0.00      2.27     0.00       20     0.00     0.00  __molecproperties_MOD_calccentersofmass
  0.00      2.27     0.00       20     0.00     0.00  __molecproperties_MOD_findppalaxes
  0.00      2.27     0.00       20     0.00     0.00  __molecproperties_MOD_recovermolecules
  0.00      2.27     0.00       20     0.00     0.00  __molecproperties_MOD_rotatepoles
  0.00      2.27     0.00       20     0.00     0.00  __molecproperties_MOD_setunpolpoles
  0.00      2.27     0.00       20     0.00     0.00  __test_scme_MOD_test_scme_monomer1
  0.00      2.27     0.00       20     0.00     0.00  __torquecm_mod_MOD_torquecm

 %         the percentage of the total running time of the
time       program used by this function.

cumulative a running sum of the number of seconds accounted
 seconds   for by this function and those listed above it.

 self      the number of seconds accounted for by this
seconds    function alone.  This is the major sort for this
           listing.

calls      the number of times this function was invoked, if
           this function is profiled, else blank.

 self      the average number of milliseconds spent in this
ms/call    function per call, if this function is profiled,
	   else blank.

 total     the average number of milliseconds spent in this
ms/call    function and its descendents per call, if this
	   function is profiled, else blank.

name       the name of the function.  This is the minor sort
           for this listing. The index shows the location of
	   the function in the gprof listing. If the index is
	   in parenthesis it shows where it would appear in
	   the gprof listing if it were to be printed.

Copyright (C) 2012-2015 Free Software Foundation, Inc.

Copying and distribution of this file, with or without modification,
are permitted in any medium without royalty provided the copyright
notice and this notice are preserved.

		     Call graph (explanation follows)


granularity: each sample hit covers 2 byte(s) for 0.44% of 2.27 seconds

index % time    self  children    called     name
                                  20             __scme_MOD_scme_calculate <cycle 1> [1]
                               59400             __calc_derivs_MOD_calcdv <cycle 1> [5]
[1]     78.0    0.24    1.53   59400+20      __scme_MOD_scme_calculate <cycle 1> [1]
                1.33    0.05    5060/5060        __calc_lower_order_MOD_calcedip_quad [3]
                0.06    0.00    5060/5060        __molecproperties_MOD_sfdsf [6]
                0.04    0.00      20/20          __inducepoles_MOD_inducedipole [7]
                0.02    0.00      20/20          __calc_higher_order_MOD_calcehigh [9]
                0.01    0.00    5060/5060        __inducepoles_MOD_induceqpole [11]
                0.01    0.00      20/20          __molecproperties_MOD_rotatepolariz [13]
                0.01    0.00      20/20          __molecproperties_MOD_addfields [12]
                0.00    0.00    5060/5060        __calcenergy_mod_MOD_calcenergy [26]
                0.00    0.00    5060/5060        __molecproperties_MOD_adddfields [27]
                0.00    0.00      20/20          __molecproperties_MOD_calccentersofmass [30]
                0.00    0.00      20/20          __molecproperties_MOD_findppalaxes [31]
                0.00    0.00      20/20          __molecproperties_MOD_rotatepoles [33]
                0.00    0.00      20/20          __molecproperties_MOD_setunpolpoles [34]
                0.00    0.00      20/20          __forcecm_mod_MOD_forcecm [29]
                0.00    0.00      20/20          __torquecm_mod_MOD_torquecm [36]
                0.00    0.00      20/20          __molecproperties_MOD_recovermolecules [32]
                                  20             __calc_derivs_MOD_calcdv <cycle 1> [5]
                                  20             __scme_MOD_scme_calculate <cycle 1> [1]
-----------------------------------------------
                1.33    0.05    5060/5060        __scme_MOD_scme_calculate <cycle 1> [1]
[3]     60.8    1.33    0.05    5060         __calc_lower_order_MOD_calcedip_quad [3]
                0.03    0.00 5009400/5029200     __mdutil_MOD_cross [8]
                0.02    0.00 2226400/2235200     __tang_toennies_MOD_tang_toenniesndf [10]
-----------------------------------------------
                                                 <spontaneous>
[4]     13.7    0.31    0.00                 __tang_toennies_MOD_tang_toenniesn [4]
-----------------------------------------------
                                  20             __scme_MOD_scme_calculate <cycle 1> [1]
[5]      8.4    0.19    0.00      20         __calc_derivs_MOD_calcdv <cycle 1> [5]
                0.00    0.00   19800/19800       __mdutil_MOD_inv6 [24]
                0.00    0.00    8800/8800        ps::dmsnasa(double const*, double*) [25]
                               59400             __scme_MOD_scme_calculate <cycle 1> [1]
-----------------------------------------------
                0.06    0.00    5060/5060        __scme_MOD_scme_calculate <cycle 1> [1]
[6]      2.6    0.06    0.00    5060         __molecproperties_MOD_sfdsf [6]
-----------------------------------------------
                0.04    0.00      20/20          __scme_MOD_scme_calculate <cycle 1> [1]
[7]      1.8    0.04    0.00      20         __inducepoles_MOD_inducedipole [7]
-----------------------------------------------
                0.00    0.00   19800/5029200     __calc_higher_order_MOD_calcehigh [9]
                0.03    0.00 5009400/5029200     __calc_lower_order_MOD_calcedip_quad [3]
[8]      1.3    0.03    0.00 5029200         __mdutil_MOD_cross [8]
-----------------------------------------------
                0.02    0.00      20/20          __scme_MOD_scme_calculate <cycle 1> [1]
[9]      0.9    0.02    0.00      20         __calc_higher_order_MOD_calcehigh [9]
                0.00    0.00   19800/5029200     __mdutil_MOD_cross [8]
                0.00    0.00    8800/2235200     __tang_toennies_MOD_tang_toenniesndf [10]
-----------------------------------------------
                0.00    0.00    8800/2235200     __calc_higher_order_MOD_calcehigh [9]
                0.02    0.00 2226400/2235200     __calc_lower_order_MOD_calcedip_quad [3]
[10]     0.9    0.02    0.00 2235200         __tang_toennies_MOD_tang_toenniesndf [10]
-----------------------------------------------
                0.01    0.00    5060/5060        __scme_MOD_scme_calculate <cycle 1> [1]
[11]     0.4    0.01    0.00    5060         __inducepoles_MOD_induceqpole [11]
-----------------------------------------------
                0.01    0.00      20/20          __scme_MOD_scme_calculate <cycle 1> [1]
[12]     0.4    0.01    0.00      20         __molecproperties_MOD_addfields [12]
-----------------------------------------------
                0.01    0.00      20/20          __scme_MOD_scme_calculate <cycle 1> [1]
[13]     0.4    0.01    0.00      20         __molecproperties_MOD_rotatepolariz [13]
-----------------------------------------------
                0.00    0.00   19800/19800       __calc_derivs_MOD_calcdv <cycle 1> [5]
[24]     0.0    0.00    0.00   19800         __mdutil_MOD_inv6 [24]
-----------------------------------------------
                0.00    0.00    8800/8800        __calc_derivs_MOD_calcdv <cycle 1> [5]
[25]     0.0    0.00    0.00    8800         ps::dmsnasa(double const*, double*) [25]
-----------------------------------------------
                0.00    0.00    5060/5060        __scme_MOD_scme_calculate <cycle 1> [1]
[26]     0.0    0.00    0.00    5060         __calcenergy_mod_MOD_calcenergy [26]
-----------------------------------------------
                0.00    0.00    5060/5060        __scme_MOD_scme_calculate <cycle 1> [1]
[27]     0.0    0.00    0.00    5060         __molecproperties_MOD_adddfields [27]
-----------------------------------------------
                0.00    0.00     640/640         __atomicforces_mod_MOD_atomicforces [38]
[28]     0.0    0.00    0.00     640         __molforce_MOD_molforce3 [28]
-----------------------------------------------
                0.00    0.00      20/20          __scme_MOD_scme_calculate <cycle 1> [1]
[29]     0.0    0.00    0.00      20         __forcecm_mod_MOD_forcecm [29]
-----------------------------------------------
                0.00    0.00      20/20          __scme_MOD_scme_calculate <cycle 1> [1]
[30]     0.0    0.00    0.00      20         __molecproperties_MOD_calccentersofmass [30]
-----------------------------------------------
                0.00    0.00      20/20          __scme_MOD_scme_calculate <cycle 1> [1]
[31]     0.0    0.00    0.00      20         __molecproperties_MOD_findppalaxes [31]
-----------------------------------------------
                0.00    0.00      20/20          __scme_MOD_scme_calculate <cycle 1> [1]
[32]     0.0    0.00    0.00      20         __molecproperties_MOD_recovermolecules [32]
-----------------------------------------------
                0.00    0.00      20/20          __scme_MOD_scme_calculate <cycle 1> [1]
[33]     0.0    0.00    0.00      20         __molecproperties_MOD_rotatepoles [33]
-----------------------------------------------
                0.00    0.00      20/20          __scme_MOD_scme_calculate <cycle 1> [1]
[34]     0.0    0.00    0.00      20         __molecproperties_MOD_setunpolpoles [34]
-----------------------------------------------
                0.00    0.00      20/20          __test_scme_MOD_test_scme_cluster_32_perf [51]
[35]     0.0    0.00    0.00      20         __test_scme_MOD_test_scme_monomer1 [35]
-----------------------------------------------
                0.00    0.00      20/20          __scme_MOD_scme_calculate <cycle 1> [1]
[36]     0.0    0.00    0.00      20         __torquecm_mod_MOD_torquecm [36]
-----------------------------------------------
                                 640             __atomicforces_mod_MOD_atomicforces [38]
[38]     0.0    0.00    0.00       0+640     __atomicforces_mod_MOD_atomicforces [38]
                0.00    0.00     640/640         __molforce_MOD_molforce3 [28]
                                 640             __atomicforces_mod_MOD_atomicforces [38]
-----------------------------------------------

 This table describes the call tree of the program, and was sorted by
 the total amount of time spent in each function and its children.

 Each entry in this table consists of several lines.  The line with the
 index number at the left hand margin lists the current function.
 The lines above it list the functions that called this function,
 and the lines below it list the functions this one called.
 This line lists:
     index	A unique number given to each element of the table.
		Index numbers are sorted numerically.
		The index number is printed next to every function name so
		it is easier to look up where the function is in the table.

     % time	This is the percentage of the `total' time that was spent
		in this function and its children.  Note that due to
		different viewpoints, functions excluded by options, etc,
		these numbers will NOT add up to 100%.

     self	This is the total amount of time spent in this function.

     children	This is the total amount of time propagated into this
		function by its children.

     called	This is the number of times the function was called.
		If the function called itself recursively, the number
		only includes non-recursive calls, and is followed by
		a `+' and the number of recursive calls.

     name	The name of the current function.  The index number is
		printed after it.  If the function is a member of a
		cycle, the cycle number is printed between the
		function's name and the index number.


 For the function's parents, the fields have the following meanings:

     self	This is the amount of time that was propagated directly
		from the function into this parent.

     children	This is the amount of time that was propagated from
		the function's children into this parent.

     called	This is the number of times this parent called the
		function `/' the total number of times the function
		was called.  Recursive calls to the function are not
		included in the number after the `/'.

     name	This is the name of the parent.  The parent's index
		number is printed after it.  If the parent is a
		member of a cycle, the cycle number is printed between
		the name and the index number.

 If the parents of the function cannot be determined, the word
 `<spontaneous>' is printed in the `name' field, and all the other
 fields are blank.

 For the function's children, the fields have the following meanings:

     self	This is the amount of time that was propagated directly
		from the child into the function.

     children	This is the amount of time that was propagated from the
		child's children to the function.

     called	This is the number of times the function called
		this child `/' the total number of times the child
		was called.  Recursive calls by the child are not
		listed in the number after the `/'.

     name	This is the name of the child.  The child's index
		number is printed after it.  If the child is a
		member of a cycle, the cycle number is printed
		between the name and the index number.

 If there are any cycles (circles) in the call graph, there is an
 entry for the cycle-as-a-whole.  This entry shows who called the
 cycle (as parents) and the members of the cycle (as children.)
 The `+' recursive calls entry shows the number of function calls that
 were internal to the cycle, and the calls entry for each member shows,
 for that member, how many times it was called from other members of
 the cycle.

Copyright (C) 2012-2015 Free Software Foundation, Inc.

Copying and distribution of this file, with or without modification,
are permitted in any medium without royalty provided the copyright
notice and this notice are preserved.

Index by function name

  [25] ps::dmsnasa(double const*, double*) [24] __mdutil_MOD_inv6 [6] __molecproperties_MOD_sfdsf
   [5] __calc_derivs_MOD_calcdv [27] __molecproperties_MOD_adddfields [28] __molforce_MOD_molforce3
   [9] __calc_higher_order_MOD_calcehigh [12] __molecproperties_MOD_addfields [1] __scme_MOD_scme_calculate
   [3] __calc_lower_order_MOD_calcedip_quad [30] __molecproperties_MOD_calccentersofmass [4] __tang_toennies_MOD_tang_toenniesn
  [26] __calcenergy_mod_MOD_calcenergy [31] __molecproperties_MOD_findppalaxes [10] __tang_toennies_MOD_tang_toenniesndf
  [29] __forcecm_mod_MOD_forcecm [32] __molecproperties_MOD_recovermolecules [35] __test_scme_MOD_test_scme_monomer1
   [7] __inducepoles_MOD_inducedipole [13] __molecproperties_MOD_rotatepolariz [36] __torquecm_mod_MOD_torquecm
  [11] __inducepoles_MOD_induceqpole [33] __molecproperties_MOD_rotatepoles (2) <cycle 1>
   [8] __mdutil_MOD_cross     [34] __molecproperties_MOD_setunpolpoles
