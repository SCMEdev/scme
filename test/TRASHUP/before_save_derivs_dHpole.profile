Flat profile:

Each sample counts as 0.01 seconds.
  %   cumulative   self              self     total           
 time   seconds   seconds    calls  ms/call  ms/call  name    
 47.30      0.35     0.35     5850     0.06     0.07  __calc_lower_order_MOD_calcedip_quad
 18.92      0.49     0.14    31500     0.00     0.02  __scme_MOD_scme_calculate
 12.16      0.58     0.09   690200     0.00     0.00  ps::dmsnasa(double const*, double*)
  9.46      0.65     0.07       50     1.40     1.42  __calc_derivs_MOD_calcdv
  4.05      0.68     0.03     5900     0.01     0.01  __molecproperties_MOD_adddfields
  4.05      0.71     0.03       50     0.60     0.62  __calc_higher_order_MOD_calcehigh
  2.70      0.73     0.02       50     0.40     0.40  __inducepoles_MOD_inducedipole
  1.35      0.74     0.01      100     0.10     0.10  __molecproperties_MOD_rotatepolariz
  0.00      0.74     0.00  1255350     0.00     0.00  __mdutil_MOD_inv6
  0.00      0.74     0.00     5850     0.00     0.00  __calcenergy_mod_MOD_calcenergy
  0.00      0.74     0.00     5850     0.00     0.00  __inducepoles_MOD_induceqpole
  0.00      0.74     0.00      750     0.00     0.00  __molforce_MOD_molforce3
  0.00      0.74     0.00      100     0.00     0.00  __molecproperties_MOD_rotatepoles
  0.00      0.74     0.00       50     0.00     0.00  __forcecm_mod_MOD_forcecm
  0.00      0.74     0.00       50     0.00     0.00  __molecproperties_MOD_findppalaxes
  0.00      0.74     0.00       50     0.00     0.00  __molecproperties_MOD_recovermolecules
  0.00      0.74     0.00       50     0.00     0.00  __test_scme_MOD_test_scme_monomer1
  0.00      0.74     0.00       50     0.00     0.00  __torquecm_mod_MOD_torquecm

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


granularity: each sample hit covers 2 byte(s) for 1.35% of 0.74 seconds

index % time    self  children    called     name
                                  50             __scme_MOD_scme_calculate <cycle 1> [1]
                               31500             __calc_derivs_MOD_calcdv <cycle 1> [5]
[1]     90.4    0.14    0.53   31500+50      __scme_MOD_scme_calculate <cycle 1> [1]
                0.35    0.09    5850/5850        __calc_lower_order_MOD_calcedip_quad [3]
                0.03    0.00      50/50          __calc_higher_order_MOD_calcehigh [6]
                0.03    0.00    5900/5900        __molecproperties_MOD_adddfields [7]
                0.02    0.00      50/50          __inducepoles_MOD_inducedipole [8]
                0.01    0.00     100/100         __molecproperties_MOD_rotatepolariz [9]
                0.00    0.00    5850/5850        __calcenergy_mod_MOD_calcenergy [21]
                0.00    0.00    5850/1255350     __mdutil_MOD_inv6 [20]
                0.00    0.00    5850/5850        __inducepoles_MOD_induceqpole [22]
                0.00    0.00     100/100         __molecproperties_MOD_rotatepoles [24]
                0.00    0.00      50/50          __molecproperties_MOD_findppalaxes [26]
                0.00    0.00      50/50          __forcecm_mod_MOD_forcecm [25]
                0.00    0.00      50/50          __torquecm_mod_MOD_torquecm [29]
                0.00    0.00      50/50          __molecproperties_MOD_recovermolecules [27]
                                  50             __calc_derivs_MOD_calcdv <cycle 1> [5]
                                  50             __scme_MOD_scme_calculate <cycle 1> [1]
-----------------------------------------------
                0.35    0.09    5850/5850        __scme_MOD_scme_calculate <cycle 1> [1]
[3]     59.3    0.35    0.09    5850         __calc_lower_order_MOD_calcedip_quad [3]
                0.09    0.00  678600/690200      ps::dmsnasa(double const*, double*) [4]
                0.00    0.00 1228500/1255350     __mdutil_MOD_inv6 [20]
-----------------------------------------------
                0.00    0.00    5800/690200      __calc_derivs_MOD_calcdv <cycle 1> [5]
                0.00    0.00    5800/690200      __calc_higher_order_MOD_calcehigh [6]
                0.09    0.00  678600/690200      __calc_lower_order_MOD_calcedip_quad [3]
[4]     12.2    0.09    0.00  690200         ps::dmsnasa(double const*, double*) [4]
-----------------------------------------------
                                  50             __scme_MOD_scme_calculate <cycle 1> [1]
[5]      9.6    0.07    0.00      50         __calc_derivs_MOD_calcdv <cycle 1> [5]
                0.00    0.00    5800/690200      ps::dmsnasa(double const*, double*) [4]
                0.00    0.00   10500/1255350     __mdutil_MOD_inv6 [20]
                               31500             __scme_MOD_scme_calculate <cycle 1> [1]
-----------------------------------------------
                0.03    0.00      50/50          __scme_MOD_scme_calculate <cycle 1> [1]
[6]      4.2    0.03    0.00      50         __calc_higher_order_MOD_calcehigh [6]
                0.00    0.00    5800/690200      ps::dmsnasa(double const*, double*) [4]
                0.00    0.00   10500/1255350     __mdutil_MOD_inv6 [20]
-----------------------------------------------
                0.03    0.00    5900/5900        __scme_MOD_scme_calculate <cycle 1> [1]
[7]      4.1    0.03    0.00    5900         __molecproperties_MOD_adddfields [7]
-----------------------------------------------
                0.02    0.00      50/50          __scme_MOD_scme_calculate <cycle 1> [1]
[8]      2.7    0.02    0.00      50         __inducepoles_MOD_inducedipole [8]
-----------------------------------------------
                0.01    0.00     100/100         __scme_MOD_scme_calculate <cycle 1> [1]
[9]      1.4    0.01    0.00     100         __molecproperties_MOD_rotatepolariz [9]
-----------------------------------------------
                0.00    0.00    5850/1255350     __scme_MOD_scme_calculate <cycle 1> [1]
                0.00    0.00   10500/1255350     __calc_derivs_MOD_calcdv <cycle 1> [5]
                0.00    0.00   10500/1255350     __calc_higher_order_MOD_calcehigh [6]
                0.00    0.00 1228500/1255350     __calc_lower_order_MOD_calcedip_quad [3]
[20]     0.0    0.00    0.00 1255350         __mdutil_MOD_inv6 [20]
-----------------------------------------------
                0.00    0.00    5850/5850        __scme_MOD_scme_calculate <cycle 1> [1]
[21]     0.0    0.00    0.00    5850         __calcenergy_mod_MOD_calcenergy [21]
-----------------------------------------------
                0.00    0.00    5850/5850        __scme_MOD_scme_calculate <cycle 1> [1]
[22]     0.0    0.00    0.00    5850         __inducepoles_MOD_induceqpole [22]
-----------------------------------------------
                0.00    0.00     750/750         __atomicforces_mod_MOD_atomicforces [31]
[23]     0.0    0.00    0.00     750         __molforce_MOD_molforce3 [23]
-----------------------------------------------
                0.00    0.00     100/100         __scme_MOD_scme_calculate <cycle 1> [1]
[24]     0.0    0.00    0.00     100         __molecproperties_MOD_rotatepoles [24]
-----------------------------------------------
                0.00    0.00      50/50          __scme_MOD_scme_calculate <cycle 1> [1]
[25]     0.0    0.00    0.00      50         __forcecm_mod_MOD_forcecm [25]
-----------------------------------------------
                0.00    0.00      50/50          __scme_MOD_scme_calculate <cycle 1> [1]
[26]     0.0    0.00    0.00      50         __molecproperties_MOD_findppalaxes [26]
-----------------------------------------------
                0.00    0.00      50/50          __scme_MOD_scme_calculate <cycle 1> [1]
[27]     0.0    0.00    0.00      50         __molecproperties_MOD_recovermolecules [27]
-----------------------------------------------
                0.00    0.00      50/50          __test_scme_MOD_test_scme_cluster_15_perf [50]
[28]     0.0    0.00    0.00      50         __test_scme_MOD_test_scme_monomer1 [28]
-----------------------------------------------
                0.00    0.00      50/50          __scme_MOD_scme_calculate <cycle 1> [1]
[29]     0.0    0.00    0.00      50         __torquecm_mod_MOD_torquecm [29]
-----------------------------------------------
                                 750             __atomicforces_mod_MOD_atomicforces [31]
[31]     0.0    0.00    0.00       0+750     __atomicforces_mod_MOD_atomicforces [31]
                0.00    0.00     750/750         __molforce_MOD_molforce3 [23]
                                 750             __atomicforces_mod_MOD_atomicforces [31]
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

   [4] ps::dmsnasa(double const*, double*) [22] __inducepoles_MOD_induceqpole [23] __molforce_MOD_molforce3
   [5] __calc_derivs_MOD_calcdv [20] __mdutil_MOD_inv6     [1] __scme_MOD_scme_calculate
   [6] __calc_higher_order_MOD_calcehigh [7] __molecproperties_MOD_adddfields [28] __test_scme_MOD_test_scme_monomer1
   [3] __calc_lower_order_MOD_calcedip_quad [26] __molecproperties_MOD_findppalaxes [29] __torquecm_mod_MOD_torquecm
  [21] __calcenergy_mod_MOD_calcenergy [27] __molecproperties_MOD_recovermolecules (2) <cycle 1>
  [25] __forcecm_mod_MOD_forcecm [9] __molecproperties_MOD_rotatepolariz
   [8] __inducepoles_MOD_inducedipole [24] __molecproperties_MOD_rotatepoles
