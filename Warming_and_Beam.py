# -*- coding: utf-8 -*-
"""

st: 25806676
"""

from __future__ import absolute_import, division, print_function
import numpy as np


def FTCSWB(phiOld, c, nt):
    "Scheme FTCSWB is using Warming and Beam method to establish"
    "the function which is centred in time and forward in space"
    "advection Courant Number, c = u*dt/dx"
    
    #argument testing
    if nt <= 0:
        raise ValueError('Arguement nt to scheme FTCSWB should be positive')
    if not isinstance(phiOld, np.ndarray):
        raise ValueError('Arguement phiOld to scheme FTCSWB should be an numpy array')
    if not (int(nt) == nt):
        raise ValueError('Arguement nt to scheme FTCSWB should be an integer')
    
    nx = len(phiOld)
    
    phi = phiOld.copy()
    
    for j in range(int(nt)):
        for i in range(nx):
            phi[i] = phiOld[i] - 0.5*c*((3 - c)*phiOld[i] - \
                     (4 - 2*c)*phiOld[(i-1)%nx] + (1 - c)*phiOld[(i-2)%nx])
        phiOld = phi.copy()
            
    return phi
    
