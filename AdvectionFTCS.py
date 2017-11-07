# Numerical schemes for simulating advection.py

from __future__ import absolute_import, division, print_function
import numpy as np


def FTCS(phiOld, c, nt):
    "Advection of profile in phiOld using FTCS using non-dimensional"
    "advection coeffient, c"
    
    nx = len(phiOld)
    
    phi = phiold.copy()
    
    for j in xrange(int(nt)):
        for i in xrange(nx):
            phi[i] = phiold[i] - 0.5*c*(phiold[i+1] - phiold[(i-1)%nx])
            
    return phi
    
   
        
        


