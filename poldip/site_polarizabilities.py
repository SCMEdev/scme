# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 18:03:32 2017

@author: j
"""

import numpy as np

data_dip_pol = np.array([
9.932447,
9.438219,
9.638610
])*1.88972666351031921149**(-3)
print data_dip_pol

data_tracepart = sum(data_dip_pol)/3.0
print data_tracepart

data_traceless = data_dip_pol -  data_tracepart

print data_traceless

data_dipquad_pol = np.array([
-6.602274,
-2.613041,
-1.016789,
 3.819314,
-2.802525
])*1.88972666351031921149**(-4)
print sum(data_dipquad_pol)