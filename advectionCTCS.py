# Numerical schemes for simulating diffusion for outer code diffusion.py

from __future__ import absolute_import, division, print_function
import numpy as np

def TV(phi):
    """
    Code for computing Total Variation for a function phi
    (with periodic boundaries).
    """
    if not(isinstance(phi,np.ndarray)):
        raise TypeError('Error in TV: Argument phi to TV should be an array')
        
    nx = len(phi)    
    I=0.
    for i in range(nx):
        I+=abs(phi[(i+1)%nx]-phi[i])
    return I
    
def CTCS(phiOld, c, nt):
    
    """
    Scheme solving the linear advection equation numerically by a CTCS finite
    difference scheme, assuming c as Courant number and nt as the number of
    time steps.
    id 25825273
    """
    
    # arguments test
    if not(int(nt) == nt):
        raise TypeError(\
            'Error in CTCS: Argument nt to CTCS should be an integer')
    if not(isinstance(float(c),float) and float(c) > 0):
        raise TypeError(\
            'Error in CTCS: Argument c to CTCS should be a positive float')
    if not(isinstance(phiOld,np.ndarray)):
        raise TypeError(\
            'Error in CTCS: Argument phiOld to CTCS should be an array')
    
    nx = len(phiOld)
    
    # new time-step array for phi
    phiNew = phiOld.copy()
    phiCurrent = phiOld.copy()
    TotalVariation = np.zeros(nt)
    
    # FTCS for computing the first time step
    for j in range(0,nx):
        phiCurrent[j] = phiOld[j] - c*0.5*(phiOld[(j+1)%nx] - phiOld[(j-1)%nx])
    # total variation for first time step    
    TotalVariation[0]=TV(phiCurrent)
    
    # CTCS for all time steps
    for it in range(1,nt):
        # spatial points
        for j in range(0,nx):
            phiNew[j] = phiOld[j] - c*\
                                   (phiCurrent[(j+1)%nx] - phiCurrent[(j-1)%nx])
        
        # output to phiOld for the next time-step
        phiOld = phiCurrent.copy()
        phiCurrent = phiNew.copy()
        # total variation for time steps 2 to nt
        TotalVariation[it]=TV(phiNew)
        
    return phiNew #TotalVariation
