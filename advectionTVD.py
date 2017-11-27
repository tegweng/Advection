# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 14:16:39 2017

@author: 25803263
"""

# TVS scheme for simulating advection for outer code advection.py

from __future__ import absolute_import, division, print_function
import numpy as np

def TVD(phiOld, c, nt, u, limiter):
    """
    Code to use implement the advection scheme TVD (Total Variation Diminishing)
    using the Courant number c and one of three limiter functions:
    Van Leer, a fixed value to recover the Lax-Wendroff scheme, or the
    Koren limiter.
    id: 25803263
    """
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
    for it in range(nt):
        #calculating the spatial points at each time step
        phiTVD = flux(phiOld.copy(), c, u, limiter)
        
        for i in range(1,nx):
        #calculating phi at n+1
            
            phi[i] = phiOld[i] - c * (phiTVD[i] - phiTVD[i-1])
        """
        #checking that the scheme is in fact TVD
        #set TV back to zero at each time step
        if it == 0:
            TVold = np.zeros(len(phi)-1)
        else:
            TVold = TV.copy()
        TV=np.zeros(len(phi)-1)
        for i in xrange(0, len(phi)-2):
            
            TV[i] = abs(phi[i+1]-phi[i])
        if it > 0:
            if np.sum(TVold) >= np.sum(TV):
                print(str(limiter) + ' at time step ' + str(it) + ' - yes!')
            else:
                print(str(limiter) + ' at time step ' + str(it) + ' - no!')
        """
        #setting phiOld to be phi to return to top of loop
        phiOld = phi.copy()
        
    
    return phi

def flux(phiOld, c, u, limiter = "VanLeer"):
    """
    Function to calculate the values of x_j+1/2 so can use values at i-1
    in the loop for calculating phi.
    id: 25803263
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
    for i in range(0,nx):
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
    
        if limiter == "VanLeer":
        
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
        elif limiter == "LaxWend":
            #constant limiter function to recover Lax-Wendroff scheme
            Lim[i] = 1

        elif limiter == "Koren":
            #Koren limiter - third order accurate for sufficiently smooth data
            Lim[i] = max(0, min(2 * r[i], (2 + r[i]) / 3, 2))
        
        
        #calculating values of phi at j+1/2 etc
        
        phiTVD[i] = Lim[i] * phiH[i] + (1 - Lim[i]) * phiL[i]
            
    return phiTVD
