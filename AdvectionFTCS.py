# Numerical schemes for simulating advection.py

from __future__ import absolute_import, division, print_function
import numpy as np

"""
id 25806676
"""

def FTCS(phiOld, c, nt):
"Advection of profile in phiOld using FTCS using non-dimensional"
"advection Courant number, c"
    
    #arguement tesing
    if nt <= 0:
        raise ValueError('Arguement nt to scheme FTCS should be positive')
    if not isinstance(phiOld, np.ndarray):
        raise ValueError('Arguement phiOld to scheme FTCS should be an numpy array')
    if not (int(nt) == nt):
        raise ValueError('Arguement nt to scheme FTCS should be an integer')
    
    nx = len(phiOld)
    
    phi = phiOld.copy()
    
    for j in range(int(nt)):
        for i in range(0,nx):
            phi[i] = phiOld[i] - 0.5*c*(phiOld[(i+1)%nx] - phiOld[(i-1)%nx])
        #assign the calculated results of next time step to phiOld
        phiOld = phi.copy()
            
    return phi
    
   
        
        


