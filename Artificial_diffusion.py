# -*- coding: utf-8 -*-
"""
Created on Tue Nov 07 14:08:55 2017

@author: 25803263
"""

# BTBS scheme for simulating advection for outer code advection.py

from __future__ import absolute_import, division, print_function
import numpy as np

# The linear algebra package for BTBS (for solving the matrix equation)
import scipy.linalg as la


def BTCS(phi, c, nt):
    "Diffusion of profile in phi using BTBS using the\
    Courant number, c, assuming fixed value boundary conditions"
    
    nx = len(phi)
    
    #array representing BTBS
    M=np.zeros([nx,nx])
    
  
    #setting initial 
    for i in range(0,nx):
        M[i,(i-1)%nx] = -c/2
        M[i,i] = 1
        M[i,(i+1)%nx] = c/2
    
    #BTBS for all time steps
    
    for it in range(int(nt)):

        #as have to do the following for all time steps it is in the for loop
        phi = la.solve(M, phi)
    
    return phi

def Artificial_diffusion2(phiOld, c, nt, dx, dt, k):

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
    d=k*dt/dx**2
    
    # new time-step array for phi
    phiNew = phiOld.copy()
    phiCurrent = phiOld.copy()

    # FTCS for both advection and diffusion to compute the first time step
    for j in range(0,nx):
        phiCurrent[j] = phiOld[j] - c*0.5*(phiOld[(j+1)%nx] - phiOld[(j-1)%nx])\
                        + d*(phiOld[(j+1)%nx] - 2*phiOld[j] + phiOld[(j-1)%nx])
    # CTCS for advection, FTCS for diffusion
    for it in range(nt):

        # spatial points
        for j in range(0,nx):
            phiNew[j] = phiOld[j] - \
                        c*(phiCurrent[(j+1)%nx] - phiCurrent[(j-1)%nx]) + \
                        2*d*(phiOld[(j+1)%nx] - 2*phiOld[j] + phiOld[(j-1)%nx])
        # output to phiOld for the next time-step
        phiOld = phiCurrent.copy()
        phiCurrent = phiNew.copy()

    return phiNew

def Artificial_diffusion4(phiOld, c, nt, dx, dt, k):

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
    d=k*dt/dx**4

    # new time-step array for phi
    phiNew = phiOld.copy()
    phiCurrent = phiOld.copy()

    # FTCS for both advection and diffusion to compute the first time step
    for j in range(0,nx):
        phiCurrent[j] = phiOld[j] - c*0.5*(phiOld[(j+1)%nx] - phiOld[(j-1)%nx])\
        - d*(phiOld[(j+2)%nx] - 4*phiOld[(j+1)%nx] + 6*phiOld[j] -\
                    4*phiOld[(j-1)%nx] + phiOld[(j-2)%nx])
    # CTCS for advection, FTCS for diffusion
    for it in range(nt):

        # spatial points
        for j in range(0,nx):
            phiNew[j] = phiOld[j] - \
                        c*(phiCurrent[(j+1)%nx] - phiCurrent[(j-1)%nx]) - \
                        2*d*(phiOld[(j+2)%nx] - 4*phiOld[(j+1)%nx] + \
                            6*phiOld[j] - 4*phiOld[(j-1)%nx] + phiOld[(j-2)%nx])
        # output to phiOld for the next time-step
        phiOld = phiCurrent.copy()
        phiCurrent = phiNew.copy()

    return phiNew












