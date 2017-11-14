# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 14:16:39 2017

@author: 25803263
"""

# TVS scheme for simulating advection for outer code advection.py

from __future__ import absolute_import, division, print_function
import numpy as np

def TVD(phiOld, c, nt, u):
    
    nx = len(phiOld)
    #don't need to use int(nx) as should already be an integer
    
    #setting up the arrays for the loop later
    phiH = np.zeros(nx)
    phiL = np.zeros(nx)
    VLL = np.zeros(nx)
    r = np.zeros(nx)
    phiTVD = np.zeros(nx)
    phi = np.zeros(nx)
    
    #making sure c and u are floats
    c = float(c)
    u = float(u)
    
    #data checking
    if nt<=0:
        raise ValueError('Argument nt to TVD must be positive and non zero.')     
    if not(isinstance(phi,np.ndarray)):
        raise TypeError('Argument phiOld to TVD must be a numpy array')
        
    #time loop
    for it in xrange(nt):
        #calculating the spatial points at each time step
        #Lax-Wendroff as the high-order flux
       
    
        for i in xrange(0,nx):
            phiH[i] = 0.5 * (1 + c) * phiOld[i] + 0.5 * (1 - c) * phiOld[(i+1)%nx]
    
        #First-order upwinds as the low-order flux
      

            if u < 0:
                phiL[i] = phiOld[i]
            else: 
                phiL[i] = phiOld[(i+1)%nx]
                
        #calculating r
            g = (phiOld[(i+1)%nx] - phiOld[i])
        #need to think about what happens if the denominator is zero
            r[i] = (phiOld[i] - phiOld[(i-1)%nx]) / (phiOld[(i+1)%nx] - phiOld[i])
    
        #calculating Van Leer Limiter    
       
            if g==0:
                VLL[i] = 2
                #as when r -> infinity VLL -> 2
            else:
                VLL[i] = ( r[i] + abs(r[i]) ) / (1 + abs(r[i]) )
    
        #calculating values of phi at j+1/2 etc
        

            phiTVD[i] = VLL[i] * phiH[i] + (1 - VLL[i]) * phiL[i]
        
        #calculating phi at n+1
        
            phi[i] = phiOld[i] - c * (phiTVD[i] - phiTVD[(i-1)%nx])
        #setting phiOld to be phi to return to top of loop
        phiOld = phi.copy()
    
    return phi


try:
    TVD(np.zeros(6), 1, 0, 0.5)
except ValueError:
    pass
else:
    print('Error in TVD, an error should be raised if nt<=0')

try:
    TVD(0,1,4, -5)
except TypeError:
    pass
else:
    print('Error in TVD, an error should be raised if phiOld is not a numpy array')