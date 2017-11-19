# -*- coding: utf-8 -*-
"""

st: 25806676
"""

from __future__ import absolute_import, division, print_function
import numpy as np


def FTCSWB(phiOld, c, nt):
    "Advection of profile in phiOld using FTCS using non-dimensional"
    "advection coeffient, c"
    
    nx = len(phiOld)
    
    phi = phiOld.copy()
    "Using the Warming and Beam to modify the FTCS scheme"    
    for j in range(int(nt)):
        for i in range(nx):
            phi[i] = phiOld[i] - 0.5*c*((3 - c)*phiOld[i] - \
                     (4 - 2*c)*phiOld[(i-1)%nx] + (1 - c)*phiOld[(i-2)%nx])
        phiOld = phi.copy()
            
    return phi
    
