# Student ID: fr818629
# FTBS numerical scheme used to solve the advection equation

from __future__ import absolute_import, division, print_function
import numpy as np

#Define a function to solve advection equation using FTBS
def FTBS(phiOld, c, nt):
    '''
    Advection of profile in phiOld using the Courant number, c.
    '''
    
    #Number of space steps is the length of the orifginal phi
    nx = len(phiOld)
    
    #new time-step array for phi
    phi = phiOld.copy()
    
    #FTBS for all time steps using arrays over time and space
    #Using modulo arithmetic for periodic boundaries
    for it in range(nt):
        for i in range(1, nx-1):
            phi[i] = phiOld[i] - c*(phiOld[i] - phiOld[(i-1)%nx])
        phiOld = phi.copy()
        
    return phi