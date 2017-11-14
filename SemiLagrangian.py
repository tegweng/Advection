# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 14:16:46 2017

@author: fr818629
"""

from __future__ import absolute_import, division, print_function
import numpy as np

def SemiLag(phiOld, c, nt, u, dt):
    '''
    Advection of profile in phiOld using Semi-Lagrangian Advection in
    one dimension using the Courant number, c.
    '''
    
    #Number of space steps is the length of the original phi
    nx = len(phiOld)   
    
    #New time-step array for phi
    phi = phiOld.copy()
    
    for it in xrange(nt):
        for j in xrange(0, nx):
            k = np.floor(j - c)
            beta = j - k - c
            jd = j - u*dt
            
            #Cubic Lagrange interpolation
            phiOld[jd] = -1/6.*beta*(1 - beta)*(2 - beta)*phiOld[(k-1)%nx] \
                         + 1/2.*(1 + beta)*(1 - beta)*(2 - beta)*phiOld[k] \
                         + 1/2.*(1 + beta)*beta*(2 - beta)*phiOld[(k+1)%nx] \
                         - 1/6.*(1 + beta)*beta*(1 - beta)*phiOld[(k+2)%nx]                         

            #Simple interpolation
            #phiOld[jd] = phiOld[k] + beta*(phiOld[(k+1)%nx] - phiOld[k])
            phi[j] = phiOld[jd]
        phiOld = phi.copy()
        
    return phi