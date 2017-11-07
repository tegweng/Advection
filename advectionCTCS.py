# Numerical schemes for simulating diffusion for outer code diffusion.py

from __future__ import absolute_import, division, print_function
import numpy as np


def CTCS(phiOld, c, nt):
    
    """
    id 25825273
    """
    
    # arguments test
    if nt<=0:
        raise ValueError('Error in CTCS: Argument nt to CTCS should be > 0')
    if not(int(nt) == nt):
        raise ValueError('Error in CTCS:\
                         Argument nt to CTCS should be an integer')
    if not(isinstance(float(c),float) and float(c) > 0):
        raise TypeError('Error in CTCS:\
                        Argument c to CTCS should be a positive float')
    if not(isinstance(phiOld,np.ndarray)):
        raise TypeError('Error in CTCS:\
                        Argument phiOld to CTCS should be an array')
    
    nx = len(phiOld)
    
    # new time-step array for phi
    phiNew = phiOld.copy()
    phiCurrent = phiOld.copy()
    
    # FTCS for computing the first time step
    for j in range(0,nx):
        phiCurrent[j] = phiOld[j] - c*0.5*(phiOld[(j+1)%nx] - phiOld[(j-1)%nx])
    # CTCS for all time steps
    for it in range(nt):
        
        # spatial points
        for j in range(0,nx):
            phiNew[j] = phiOld[j] - c*\
                                   (phiCurrent[(j+1)%nx] - phiCurrent[(j-1)%nx])
        # output to phiOld for the next time-step
        phiOld = phiCurrent.copy()
        phiCurrent = phiNew.copy()
        
    return phiNew









