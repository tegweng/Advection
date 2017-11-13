# Numerical schemes for simulating advection.py

from __future__ import absolute_import, division, print_function
import numpy as np

"""
id 25806676
"""

def FTCS(phiOld, c, nt):
    "Advection of profile in phiOld using FTCS using non-dimensional"
    "advection coeffient, c"
    
    nx = len(phiOld)
    
    phi = phiOld.copy()
    
    for j in range(int(nt)):
        for i in range(nx-1):
            phi[i] = phiOld[i] - 0.5*c*(phiOld[i+1] - phiOld[(i-1)%nx])
        #assign the calculated results of next time step to phiOld
        phiOld = phi.copy()
            
    return phi
    
   
        
        


