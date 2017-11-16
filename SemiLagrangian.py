# Student ID: 25818629
# Semi-Lagrangian numerical scheme used to solve the advection equation

from __future__ import absolute_import, division, print_function
import numpy as np

def SemiLag(phiOld, c, nt, u, dt):
    '''
    Advection of profile in phiOld using Semi-Lagrangian Advection in
    one dimension using the Courant number, c.
    '''
    
     #Error catching
    if nt <= 0:
        raise ValueError('Error in FTBS: Argument nt should be > 0')
    if not(int(nt) == nt):
        raise ValueError('Error in FTBS: Argument nt should be an integer')
    if not(isinstance(phiOld,np.ndarray)):
        raise TypeError('Error in FTBS: Argument phiOld should be an array')
    
    #Number of space steps is the length of the original phi
    nx = len(phiOld)   
    
    #New time-step array for phi
    phi = phiOld.copy()
    
    for it in range(nt):
        for j in range(0, nx):  
            k = int(np.floor(j - c))     
            beta = j - k - c
            
            #Cubic Lagrange interpolation
            phiOld[j] = -1/6.*beta*(1 - beta)*(2 - beta)*phiOld[(k-1)%nx] \
                         + 1/2.*(1 + beta)*(1 - beta)*(2 - beta)*phiOld[k] \
                         + 1/2.*(1 + beta)*beta*(2 - beta)*phiOld[(k+1)%nx] \
                         - 1/6.*(1 + beta)*beta*(1 - beta)*phiOld[(k+2)%nx]                         

            #Simple interpolation
            #phiOld[jd] = phiOld[k] + beta*(phiOld[(k+1)%nx] - phiOld[k])
            phi[j] = phiOld[jd]
            
        phiOld = phi.copy()
        
    return phi
