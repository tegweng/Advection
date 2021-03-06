# -*- coding: utf-8 -*-
"""
Created on Tue Nov 07 14:08:55 2017
Code for implementing artificial diffusion into the linear advection scheme.
@author: 25825273
"""

from __future__ import absolute_import, division, print_function
import numpy as np

def TV(phi):
    """
    Code for computing Total Variation for a function phi
    (with periodic boundaries).
    """
    if not(isinstance(phi,np.ndarray)):
        raise TypeError('Error in TV: Argument phi to TV should be an array')
        
    nx = len(phi)    
    I=0.
    for i in range(nx):
        I+=abs(phi[(i+1)%nx]-phi[i])
    return I

def Artificial_diffusion(phiOld, c, nt, dx, dt, d, orderAD=2):

    """
    Scheme for solving linear advection equation with an added artificial
    diffusion of order (orderAD) which can be 2nd (default),
          DΦ
          -- - k*del^2(Φ) = 0       where D/Dt is the material derivative,
          Dt
    or 4th,
          DΦ
          -- + k*del^4(Φ) = 0
          Dt
    The function returns the array for phi at time nt*dt and the value of the
    diffusivity constant k (accordingly with the order of diffusivity selected)
    id 25825273
    """

    # arguments test
    if not(isinstance(nt,int)):
        raise TypeError('Error in Artificial_diffusion:\
                Argument nt to Artificial_diffusion should be an integer')
    if not(isinstance(float(c),float) and float(c) > 0):
        raise TypeError('Error in Artificial_diffusion:\
                Argument c to Artificial_diffusion should be a positive float')
    if not(isinstance(float(dx),float) and float(dx) > 0):
        raise TypeError('Error in Artificial_diffusion:\
                Argument dx to Artificial_diffusion should be a positive float')
    if not(isinstance(float(dt),float) and float(dt) > 0):
        raise TypeError('Error in Artificial_diffusion:\
                Argument dt to Artificial_diffusion should be a positive float')
    if not(isinstance(float(d),float) and float(d) >= 0):
        raise TypeError('Error in Artificial_diffusion:\
                Argument d to Artificial_diffusion should be a positive float')
    if not(isinstance(phiOld,np.ndarray)):
        raise TypeError('Error in Artificial_diffusion:\
                Argument phiOld to Artificial_diffusion should be an array')

    nx = len(phiOld)
    
    # calculating value of dimensional diffusivity for informative purposes
    k = d*dx**orderAD/dt
    
    # new time-step array for phi
    phiNew = phiOld.copy()
    phiCurrent = phiOld.copy()
    # initialising Total Variation for monotonicity analysis
    TotalVariation = np.zeros(nt)
    
    if (orderAD==2):
        # FTCS for both advection and diffusion to compute the first time step
        for j in range(0,nx):
            phiCurrent[j] = phiOld[j] - c*0.5*(phiOld[(j+1)%nx] - \
                      phiOld[(j-1)%nx]) + d*(phiOld[(j+1)%nx] - 2*phiOld[j]\
                             + phiOld[(j-1)%nx])
        
        # total variation for first time step
        TotalVariation[0]=TV(phiCurrent)
        
        # CTCS for advection, FTCS for diffusion
        for it in range(1,nt):
            # spatial points
            for j in range(0,nx):
                phiNew[j] = phiOld[j] - \
                c*(phiCurrent[(j+1)%nx] - phiCurrent[(j-1)%nx]) + \
                2*d*(phiOld[(j+1)%nx] - 2*phiOld[j] + phiOld[(j-1)%nx])
            
            # total variation for time step 2 to nt
            TotalVariation[it]=TV(phiNew)
            
            # output to phiOld for the next time-step
            phiOld = phiCurrent.copy()
            phiCurrent = phiNew.copy()
    
    if (orderAD==4):
        # FTCS for both advection and diffusion to compute the first time step
        for j in range(0,nx):
            phiCurrent[j] = phiOld[j] - c*0.5*(phiOld[(j+1)%nx]\
                      - phiOld[(j-1)%nx]) - d*(phiOld[(j+2)%nx]\
                               - 4*phiOld[(j+1)%nx] + 6*phiOld[j]\
                               - 4*phiOld[(j-1)%nx] + phiOld[(j-2)%nx])
            
        # total variation for first time step
        TotalVariation[0]=TV(phiCurrent)
        
        # CTCS for advection, FTCS for diffusion
        for it in range(1, nt):
            # spatial points
            for j in range(0,nx):
                phiNew[j] = phiOld[j] - \
                c*(phiCurrent[(j+1)%nx] - phiCurrent[(j-1)%nx]) - \
                2*d*(phiOld[(j+2)%nx] - 4*phiOld[(j+1)%nx] + \
                            6*phiOld[j] - 4*phiOld[(j-1)%nx] + phiOld[(j-2)%nx])
            
            # total variation for time step 2 to nt
            TotalVariation[it]=TV(phiNew)
            
            # output to phiOld for the next time-step
            phiOld = phiCurrent.copy()
            phiCurrent = phiNew.copy()
    
    if not((orderAD==2) or (orderAD==4)):
        raise ValueError(\
            'Error in Artificial_diffusion:\
           Argument orderAD to Artificial_diffusion should be 2 or 4 (integer)')

    return phiNew, k, TotalVariation
