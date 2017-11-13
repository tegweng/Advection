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
   
    #setting initial 
    for i in range(0,nx):
        M[i,(i-1)%nx] = -c
        M[i,i] = 1+c
    
    #BTBS for all time steps
    
    for it in range(int(nt)):
        #as have to do the following for all time steps it is in the for loop
        phi = la.solve(M, phi)
    
    return phi