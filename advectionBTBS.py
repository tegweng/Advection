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


def BTBS(phi, c, nt):
    "Diffusion of profile in phi using BTBS using the\
    Courant number, c, assuming fixed value boundary conditions"
    
    nx = len(phi)
    
    #array representing BTBS
    M=np.zeros([nx,nx])
    
    #zero gradient boundary conditions for initial time-step
    
    M[0,0] = 0.5
    M[-1,-1] = 0.5

    
    #setting initial 
    for i in xrange(1,nx-1):
        M[i,i-1] = -c
        M[i,i] = 1+c
    
    #BTBS for all time steps
    
    for it in xrange(int(nt)):
        #RHS for zero gradient boundary conditions - have to keep resetting the \
        #boundary to be 0 at each time step
        phi[0] = 0.
        phi[-1] = 0.
        #as have to do the following for all time steps it is in the for loop
        phi = la.solve(M, phi)
    
    return phi