Flat profile:

Each sample counts as 0.01 seconds.
  %   cumulative   self              self     total           
 time   seconds   seconds    calls  ms/call  ms/call  name    
 79.00      1.58     1.58     5850     0.27     0.29  __calc_lower_order_MOD_calcedip_quad
  8.00      1.74     0.16    31550     0.01     0.06  __scme_MOD_scme_calculate
  5.00      1.84     0.10   684400     0.00     0.00  __tang_toennies_MOD_tang_toenniesn
  4.00      1.92     0.08       50     1.60     1.60  __calc_derivs_MOD_calcdv
  1.50      1.95     0.03     5850     0.01     0.01  __inducepoles_MOD_induceqpole
  1.50      1.98     0.03     5850     0.01     0.01  __molecproperties_MOD_adddfields
  0.50      1.99     0.01                             ps::dmsnasa(double const*, double*)
  0.50      2.00     0.01                             ps::potnasa(double const*, double*, double*)
  0.00      2.00     0.00  1239000     0.00     0.00  __molecproperties_MOD_sf
  0.00      2.00     0.00    10500     0.00     0.00  __molecproperties_MOD_sfdsf
  0.00      2.00     0.00     5850     0.00     0.00  __inducepoles_MOD_inducedipole
  0.00      2.00     0.00     5850     0.00     0.00  __molecproperties_MOD_addfields
  0.00      2.00     0.00     5800     0.00     0.00  __tang_toennies_MOD_tang_toenniesndf
  0.00      2.00     0.00      750     0.00     0.00  __mdutil_MOD_inv6
  0.00      2.00     0.00      750     0.00     0.00  __molforce_MOD_molforce3
  0.00      2.00     0.00       50     0.00     0.00  __atomicforces_mod_MOD_atomicforces
  0.00      2.00     0.00       50     0.00     0.02  __calc_higher_order_MOD_calcehigh
  0.00      2.00     0.00       50     0.00     0.00  __calcenergy_mod_MOD_calcenergy
  0.00      2.00     0.00       50     0.00     0.00  __dispersion_mod_MOD_dispersion
  0.00      2.00     0.00       50     0.00     0.00  __forcecm_mod_MOD_forcecm
  0.00      2.00     0.00       50     0.00     0.00  __molecproperties_MOD_calccentersofmass
  0.00      2.00     0.00       50     0.00     0.00  __molecproperties_MOD_findppalaxes
  0.00      2.00     0.00       50     0.00     0.00  __molecproperties_MOD_recovermolecules
  0.00      2.00     0.00       50     0.00     0.00  __molecproperties_MOD_rotatepolariz
  0.00      2.00     0.00       50     0.00     0.00  __molecproperties_MOD_rotatepoles
  0.00      2.00     0.00       50     0.00     0.00  __molecproperties_MOD_setunpolpoles
  0.00      2.00     0.00       50     0.00     0.00  __torquecm_mod_MOD_torquecm

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


granularity: each sample hit covers 2 byte(s) for 0.50% of 2.00 seconds

index % time    self  children    called     name
[1]     99.0    0.24    1.74      50+31550   <cycle 1 as a whole> [1]
                0.16    1.74   31550             __scme_MOD_scme_calculate <cycle 1> [3]
                0.08    0.00      50             __calc_derivs_MOD_calcdv <cycle 1> [6]
-----------------------------------------------
                                                 <spontaneous>
[2]     99.0    0.00    1.98                 __test_scme_MOD_test_scme_cluster_15_perf [2]
                0.24    1.74      50/50          __scme_MOD_scme_calculate <cycle 1> [3]
-----------------------------------------------
                               31500             __calc_derivs_MOD_calcdv <cycle 1> [6]
                0.24    1.74      50/50          __test_scme_MOD_test_scme_cluster_15_perf [2]
[3]     95.0    0.16    1.74   31550         __scme_MOD_scme_calculate <cycle 1> [3]
                1.58    0.10    5850/5850        __calc_lower_order_MOD_calcedip_quad [4]
                0.03    0.00    5850/5850        __molecproperties_MOD_adddfields [8]
                0.03    0.00    5850/5850        __inducepoles_MOD_induceqpole [7]
                0.00    0.00      50/50          __calc_higher_order_MOD_calcehigh [11]
                0.00    0.00    5850/5850        __molecproperties_MOD_addfields [24]
                0.00    0.00    5850/5850        __inducepoles_MOD_inducedipole [23]
                0.00    0.00      50/50          __molecproperties_MOD_recovermolecules [34]
                0.00    0.00      50/50          __molecproperties_MOD_calccentersofmass [32]
                0.00    0.00      50/50          __molecproperties_MOD_findppalaxes [33]
                0.00    0.00      50/50          __molecproperties_MOD_rotatepoles [36]
                0.00    0.00      50/50          __molecproperties_MOD_setunpolpoles [37]
                0.00    0.00      50/50          __molecproperties_MOD_rotatepolariz [35]
                0.00    0.00      50/50          __forcecm_mod_MOD_forcecm [31]
                0.00    0.00      50/50          __torquecm_mod_MOD_torquecm [38]
                0.00    0.00      50/50          __atomicforces_mod_MOD_atomicforces [28]
                0.00    0.00      50/50          __calcenergy_mod_MOD_calcenergy [29]
                0.00    0.00      50/50          __dispersion_mod_MOD_dispersion [30]
                                  50             __calc_derivs_MOD_calcdv <cycle 1> [6]
-----------------------------------------------
                1.58    0.10    5850/5850        __scme_MOD_scme_calculate <cycle 1> [3]
[4]     84.0    1.58    0.10    5850         __calc_lower_order_MOD_calcedip_quad [4]
                0.10    0.00  678600/684400      __tang_toennies_MOD_tang_toenniesn [5]
                0.00    0.00 1228500/1239000     __molecproperties_MOD_sf [21]
-----------------------------------------------
                0.00    0.00    5800/684400      __calc_higher_order_MOD_calcehigh [11]
                0.10    0.00  678600/684400      __calc_lower_order_MOD_calcedip_quad [4]
[5]      5.0    0.10    0.00  684400         __tang_toennies_MOD_tang_toenniesn [5]
-----------------------------------------------
                                  50             __scme_MOD_scme_calculate <cycle 1> [3]
[6]      4.0    0.08    0.00      50         __calc_derivs_MOD_calcdv <cycle 1> [6]
                0.00    0.00   10500/10500       __molecproperties_MOD_sfdsf [22]
                0.00    0.00    5800/5800        __tang_toennies_MOD_tang_toenniesndf [25]
                               31500             __scme_MOD_scme_calculate <cycle 1> [3]
-----------------------------------------------
                0.03    0.00    5850/5850        __scme_MOD_scme_calculate <cycle 1> [3]
[7]      1.5    0.03    0.00    5850         __inducepoles_MOD_induceqpole [7]
-----------------------------------------------
                0.03    0.00    5850/5850        __scme_MOD_scme_calculate <cycle 1> [3]
[8]      1.5    0.03    0.00    5850         __molecproperties_MOD_adddfields [8]
-----------------------------------------------
                                                 <spontaneous>
[9]      0.5    0.01    0.00                 ps::dmsnasa(double const*, double*) [9]
-----------------------------------------------
                                                 <spontaneous>
[10]     0.5    0.01    0.00                 ps::potnasa(double const*, double*, double*) [10]
-----------------------------------------------
                0.00    0.00      50/50          __scme_MOD_scme_calculate <cycle 1> [3]
[11]     0.0    0.00    0.00      50         __calc_higher_order_MOD_calcehigh [11]
                0.00    0.00    5800/684400      __tang_toennies_MOD_tang_toenniesn [5]
                0.00    0.00   10500/1239000     __molecproperties_MOD_sf [21]
-----------------------------------------------
                0.00    0.00   10500/1239000     __calc_higher_order_MOD_calcehigh [11]
                0.00    0.00 1228500/1239000     __calc_lower_order_MOD_calcedip_quad [4]
[21]     0.0    0.00    0.00 1239000         __molecproperties_MOD_sf [21]
-----------------------------------------------
                0.00    0.00   10500/10500       __calc_derivs_MOD_calcdv <cycle 1> [6]
[22]     0.0    0.00    0.00   10500         __molecproperties_MOD_sfdsf [22]
-----------------------------------------------
                0.00    0.00    5850/5850        __scme_MOD_scme_calculate <cycle 1> [3]
[23]     0.0    0.00    0.00    5850         __inducepoles_MOD_inducedipole [23]
-----------------------------------------------
                0.00    0.00    5850/5850        __scme_MOD_scme_calculate <cycle 1> [3]
[24]     0.0    0.00    0.00    5850         __molecproperties_MOD_addfields [24]
-----------------------------------------------
                0.00    0.00    5800/5800        __calc_derivs_MOD_calcdv <cycle 1> [6]
[25]     0.0    0.00    0.00    5800         __tang_toennies_MOD_tang_toenniesndf [25]
-----------------------------------------------
                0.00    0.00     750/750         __molforce_MOD_molforce3 [27]
[26]     0.0    0.00    0.00     750         __mdutil_MOD_inv6 [26]
-----------------------------------------------
                0.00    0.00     750/750         __atomicforces_mod_MOD_atomicforces [28]
[27]     0.0    0.00    0.00     750         __molforce_MOD_molforce3 [27]
                0.00    0.00     750/750         __mdutil_MOD_inv6 [26]
-----------------------------------------------
                0.00    0.00      50/50          __scme_MOD_scme_calculate <cycle 1> [3]
[28]     0.0    0.00    0.00      50         __atomicforces_mod_MOD_atomicforces [28]
                0.00    0.00     750/750         __molforce_MOD_molforce3 [27]
-----------------------------------------------
                0.00    0.00      50/50          __scme_MOD_scme_calculate <cycle 1> [3]
[29]     0.0    0.00    0.00      50         __calcenergy_mod_MOD_calcenergy [29]
-----------------------------------------------
                0.00    0.00      50/50          __scme_MOD_scme_calculate <cycle 1> [3]
[30]     0.0    0.00    0.00      50         __dispersion_mod_MOD_dispersion [30]
-----------------------------------------------
                0.00    0.00      50/50          __scme_MOD_scme_calculate <cycle 1> [3]
[31]     0.0    0.00    0.00      50         __forcecm_mod_MOD_forcecm [31]
-----------------------------------------------
                0.00    0.00      50/50          __scme_MOD_scme_calculate <cycle 1> [3]
[32]     0.0    0.00    0.00      50         __molecproperties_MOD_calccentersofmass [32]
-----------------------------------------------
                0.00    0.00      50/50          __scme_MOD_scme_calculate <cycle 1> [3]
[33]     0.0    0.00    0.00      50         __molecproperties_MOD_findppalaxes [33]
-----------------------------------------------
                0.00    0.00      50/50          __scme_MOD_scme_calculate <cycle 1> [3]
[34]     0.0    0.00    0.00      50         __molecproperties_MOD_recovermolecules [34]
-----------------------------------------------
                0.00    0.00      50/50          __scme_MOD_scme_calculate <cycle 1> [3]
[35]     0.0    0.00    0.00      50         __molecproperties_MOD_rotatepolariz [35]
-----------------------------------------------
                0.00    0.00      50/50          __scme_MOD_scme_calculate <cycle 1> [3]
[36]     0.0    0.00    0.00      50         __molecproperties_MOD_rotatepoles [36]
-----------------------------------------------
                0.00    0.00      50/50          __scme_MOD_scme_calculate <cycle 1> [3]
[37]     0.0    0.00    0.00      50         __molecproperties_MOD_setunpolpoles [37]
-----------------------------------------------
                0.00    0.00      50/50          __scme_MOD_scme_calculate <cycle 1> [3]
[38]     0.0    0.00    0.00      50         __torquecm_mod_MOD_torquecm [38]
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

   [9] ps::dmsnasa(double const*, double*) [7] __inducepoles_MOD_induceqpole [21] __molecproperties_MOD_sf
  [10] ps::potnasa(double const*, double*, double*) [26] __mdutil_MOD_inv6 [22] __molecproperties_MOD_sfdsf
  [28] __atomicforces_mod_MOD_atomicforces [8] __molecproperties_MOD_adddfields [27] __molforce_MOD_molforce3
   [6] __calc_derivs_MOD_calcdv [24] __molecproperties_MOD_addfields [3] __scme_MOD_scme_calculate
  [11] __calc_higher_order_MOD_calcehigh [32] __molecproperties_MOD_calccentersofmass [5] __tang_toennies_MOD_tang_toenniesn
   [4] __calc_lower_order_MOD_calcedip_quad [33] __molecproperties_MOD_findppalaxes [25] __tang_toennies_MOD_tang_toenniesndf
  [29] __calcenergy_mod_MOD_calcenergy [34] __molecproperties_MOD_recovermolecules [38] __torquecm_mod_MOD_torquecm
  [30] __dispersion_mod_MOD_dispersion [35] __molecproperties_MOD_rotatepolariz [1] <cycle 1>
  [31] __forcecm_mod_MOD_forcecm [36] __molecproperties_MOD_rotatepoles
  [23] __inducepoles_MOD_inducedipole [37] __molecproperties_MOD_setunpolpoles
