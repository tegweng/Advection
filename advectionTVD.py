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
    #Lax-Wendroff as the high-order flux
    phiH = np.zeros(nx)
    
    for i in xrange(0,nx):
        phiH[i] = 0.5 * (1 + c) * phiOld[i] +0.5 * (1 - c) * phiOld[(i+1)%nx]
    
    #First-order upwinds as the low-order flux
    phiL = np.zeros(nx)
    
    for i in xrange(0, nx):    
        if u < 0:
            phiL[i] = phiOld[i]
        else: 
            phiL[i] = phiOld[(i+1)%nx]
    
    #calculating r
    r = np.zeros(nx) 
    
    for i in xrange(0, nx):
        r[i] = (phiOld[i] - phiOld[(i-1)%nx]) / (phiOld[(i+1)%nx] - phiOld[i])
    
    #calculating Van Leer Limiter    
    VLL = np.zeros(nx)
    
    for i in xrange(0,nx):
        VLL[i] = ( r[i] + abs(r[i]) ) / (1 + abs(r[i]) )
    
    #calculating new value of phi
    phi = np.zeros(0,nx)
    for i in xrange(0, nx):
        phi[i] = VLL[i] * phiH[i] + (1 - VLL[i]) * phiL[i]
        
    #calculating the solution after so many time steps
    
    
    
    return phi
