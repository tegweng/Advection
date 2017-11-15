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
    
    phi = np.zeros(nx)
    
    #making sure c and u are floats
    c = float(c)
    u = float(u)
    
    #data checking
    if nt<=0:
        raise ValueError('Argument nt to TVD must be positive and non zero.')
        
    if not(isinstance(phiOld,np.ndarray)):
        raise TypeError('Argument phiOld to TVD must be a numpy array')
        
    #time loop
    for it in xrange(nt):
        #calculating the spatial points at each time step
        phiTVD = flux(phiOld.copy(), c, u)
        
        for i in xrange(1,nx):
        #calculating phi at n+1
            
            phi[i] = phiOld[i] - c * (phiTVD[i] - phiTVD[i-1])
            
        #setting phiOld to be phi to return to top of loop
        phiOld = phi.copy()
    
    return phi

def flux(phiOld, c, u):
    """
    Function to calculate the values of x_j+1/2 so can use values at i-1
    in the loop for calculating phi
    """
    nx = len(phiOld)
    #don't need to use int(nx) as should already be an integer
    
    #making sure c and u are floats
    c = float(c)
    u = float(u)
    #data checking
    if not(isinstance(phiOld,np.ndarray)):
        raise TypeError('Argument phiOld to TVD must be a numpy array')
    #setting up the arrays for the loop later
    phiH = np.zeros(nx)
    phiL = np.zeros(nx)
    Lim = np.zeros(nx)
    r = np.zeros(nx)
    phiTVD = np.zeros(nx)
    
    #Loop over spatial coordinates
    for i in xrange(0,nx):
        #Lax-Wendroff as the high-order flux
        phiH[i] = 0.5 * (1 + c) * phiOld[i] + 0.5 * (1 - c) * phiOld[(i+1)%nx]
    
        #First-order upwind as the low-order flux
      
        if u >= 0:
            phiL[i] = phiOld[i]
        else: 
            phiL[i] = phiOld[(i+1)%nx]
                
        #calculating r
        g = phiOld[(i+1)%nx] - phiOld[i]
        h = phiOld[i] - phiOld[i-1]
            
        r[i] = h / g
    
        #calculating Limiter function    
        
        #Van Leer limiter function    
        if g == 0:
            if h ==0:
                Lim[i] = 1
                #as when they both tend to zero the answer will tend to one
            else: 
                Lim[i] = 2
                #as when r -> infinity VLL -> 2
        else:
            Lim[i] = ( r[i] + abs(r[i]) ) / (1 + abs(r[i]) )
        """   
        #constant limiter function to recover Lax-Wendroff scheme
        Lim[i] = 1

        #Koren limiter - third order accurate for sufficiently smooth data
        Lim[i] = max(0, min(2 * r[i], (2 + r[i]) / 3, 2))
        
        """
        #calculating values of phi at j+1/2 etc
        
        phiTVD[i] = Lim[i] * phiH[i] + (1 - Lim[i]) * phiL[i]
            
    return phiTVD

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
    
try:
    flux(0,1,4)
except TypeError:
    pass
else:
    print('Error in flux, an error should be raised if phiOld is not a numpy array')