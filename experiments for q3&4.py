# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 23:31:44 2017

@author: 25806676
"""

#Designing the experiment to compare 
#the CTCS, FTBS and Warming and Beam scheme

from __future__ import absolute_import, division, print_function

import matplotlib.pyplot as plt

import numpy as np

# read in all the linear advection schemes, initial conditions and other
# code associated with this application (use with runfile if exec not supported)

execfile("advectionBTBS.py")
execfile("Advection_FTBS.py")
execfile("advectionFTCS.py")
execfile("advectionCTCS.py")
execfile("diagnostics.py")
execfile("initialConditions.py")
execfile("Artificial_diffusion.py")
execfile("advectionTVD.py")
execfile("SemiLagrangian.py")
execfile("advectionTVD.py")
execfile("Warming_and_Beam.py")

# Define a function to compare the schemes' monotonicity by calculating their 
# total variation from the time steps nt=0 until 60.
def ex_1(xmin = 0., xmax = 1., nx = 40, u = 1, n = 60, c = 0.125, \
         dt = 0.003125, squareWaveMin = 0.0, squareWaveMax = 0.5, \
         func = squareWave, name_fig = 'experiment1'):
             
    dx = (xmax - xmin)/nx
    L2WB = np.zeros(n)
    L2CTCS = np.zeros(n)
    L2FTBS = np.zeros(n)
    # spatial points for plotting and for defining initial conditions
    x = np.zeros(nx)
    for j in range(nx):
        x[j] = xmin + j*dx
        
    #initial conditions (Each line is a different initial condition)
    
    phiOld = func(x,squareWaveMin, squareWaveMax)

    
    TV1 = np.zeros(n)
    TV2 = np.zeros(n)
    TV3 = np.zeros(n)
    L2WB = np.zeros(n)
    L2CTCS = np.zeros(n)
    L2FTBS = np.zeros(n)
    nt = np.zeros(n)
    
    # When nt=0, there is no numerical solutions
    nt[0] = 0
    # the value of total variation should be the same as phiOld's
    TV1[0] = total_variation(phiOld)
    TV2[0] = total_variation(phiOld)
    TV3[0] = total_variation(phiOld)
    
    for i in range(1,n):
        nt[i] = 0  + 1*i
        T = nt[i]*dt
        phiAnalytic = func((x - u * T)%1.,squareWaveMin,squareWaveMax)
        
        phiWB = WB(phiOld.copy(), c, nt[i])
        phiCTCS = CTCS(phiOld.copy(), c, nt[i])
        phiFTBS = FTBS(phiOld.copy(), c, nt[i])
        
        TV1[i] = total_variation(phiWB)
        TV2[i] = total_variation(phiCTCS)
        TV3[i] = total_variation(phiFTBS)
        

    font = {'size' : 10}
    plt.rc('font', **font)
    plt.figure(1)
    plt.clf()
    plt.ion()
    plt.plot(x, phiOld, label='Initial', color='black',linestyle='--')
    plt.plot(x, phiAnalytic, label='Analytic', color='black', linewidth=2)
    plt.plot(x, phiCTCS, label = 'CTCS', color = 'red', linewidth=1.5)
    plt.plot(x, phiFTBS, label = 'FTBS', color = 'green', linewidth=1.5)
    plt.plot(x, phiWB, label='WB', color='blue', linewidth=1.5)
    plt.xlabel('$x$')
    plt.ylabel('$phi$')
    plt.axhline(0, linestyle=':', color='black', linewidth=2)
    plt.legend()
    plt.savefig('Plots/' + func.__name__ + '.pdf') 
    
    # plot the value of total variation           
    plt.figure(2)
    plt.clf()
    plt.ion()
    plt.plot(nt, TV1, label='TV of WB', color = 'blue', linewidth=1.5)
    plt.plot(nt, TV2, label='TV of CTCS', color = 'red', linewidth=1.5)
    plt.plot(nt, TV3, label='TV of FTBS', color = 'green', linewidth=1.5)
    plt.xlabel('$nt$')
    plt.ylabel('$Total  Variation$')
    plt.axhline(2, linestyle=':', color='black', linewidth=2)
    plt.legend(bbox_to_anchor=(0.45,1.05))
    plt.savefig('Plots/' + name_fig + 'a_' + \
               func.__name__ + '.pdf')
    
    # Compare the L2 error of CTCS, FTBS and WB from nt = 20
    for i in range(0,n):
        nt[i] = 20  + 1*i
        T = nt[i]*dt
        phiOld = func(x,squareWaveMin, squareWaveMax)
        phiAnalytic = func((x - u * T)%1.,squareWaveMin,squareWaveMax)
        
        phiWB = WB(phiOld.copy(), c, nt[i])
        phiCTCS = CTCS(phiOld.copy(), c, nt[i])
        phiFTBS = FTBS(phiOld.copy(), c, nt[i])
        
        L2WB[i] = L2ErrorNorm(phiWB, phiAnalytic)
        L2CTCS[i] = L2ErrorNorm(phiCTCS, phiAnalytic)
        L2FTBS[i] = L2ErrorNorm(phiFTBS, phiAnalytic)
    
        
    # plot the L2 error               
    plt.figure(3)
    plt.clf()
    plt.ion()
    plt.plot(nt, L2WB, label='L2 error of WB', color = 'blue')
    plt.plot(nt, L2CTCS, label='L2 error of CTCS', color = 'red')
    plt.plot(nt, L2FTBS, label='L2 error of FTBS', color = 'green')
    plt.xlabel('$nt$')
    plt.ylabel('$L2 Error$')
    plt.legend(bbox_to_anchor=(0.45,1.05))
    plt.savefig('Plots/' + name_fig + 'b_' + \
               func.__name__ + '.pdf')
         

# the funtion of calculating total variation
def total_variation(phi):
    
    nx = len(phi)
    tvd = 0.
    
    for i in range(0,nx-1):
        tvd = tvd + np.abs(phi[i+1] - phi[i])
        
    return tvd

# experiment for question 4
def ex_2(xmin = 0., xmax = 1., nx = 40, n = 30, u = 1, T= 0.125, d=0.1012, \
         squareWaveMin = 0.0, squareWaveMax = 0.5, c = 0.125, \
         func = squareWave, name_fig='experiment2', limiter = "LaxWend"):
        
    dx = np.zeros(n)
    dt = np.zeros(n)
    nt = np.zeros(n)
    nx = np.zeros(n)
    
    tr1 = np.zeros(n)
    tr2 = np.zeros(n)
    tr3 = np.zeros(n)
    tr4 = np.zeros(n)
    tr5 = np.zeros(n)
    tr6 = np.zeros(n)
    
    # set nt and nx from 50 but do calculating every 50 steps
    for i in range(0,n):
        nt[i] = 50 + 50*i
        nx[i] = 50 + 50*i
        
        dx[i] = (xmax - xmin)/nx[i]
        dt[i] = T/nt[i]
        x = np.zeros(int(nx[i]))
        for j in range(int(nx[i])):
            x[j] = xmin + j*dx[i]
            
        phiOld = func(x, squareWaveMin, squareWaveMax)
        phiAnalytic = func((x - u*T)%1.,squareWaveMin, squareWaveMax)
        
        phiTVD = TVD(phiOld.copy(), c, int(nt[i]), u, limiter)
        phidiff2 = Artificial_diffusion(phiOld.copy(), c, int(nt[i]), dx[i], dt[i], 0.1012, 2)[0]
        phiSemiLag = SemiLag(phiOld.copy(), c, int(nt[i]), u, dt[i])
        phiWB = WB(phiOld.copy(), c, int(nt[i]))
        phiCTCS = CTCS(phiOld.copy(), c, int(nt[i]))
        phiFTBS = FTBS(phiOld.copy(), c, int(nt[i]))
    
        tr1[i] = trun(phiAnalytic, phiCTCS)
        tr2[i] = trun(phiAnalytic, phiFTBS)
        tr3[i] = trun(phiAnalytic, phiWB)
        tr4[i] = trun(phiAnalytic, phiTVD)
        tr5[i] = trun(phiAnalytic, phiSemiLag)
        tr6[i] = trun(phiAnalytic, phidiff2)
    
    font = {'size' : 12}
    plt.rc('font', **font)
    plt.figure(1)
    plt.clf()
    plt.ion()
    plt.plot(dt, tr1, label='CTCS', color='red')
    plt.plot(dt, tr2, label='FTBS', color='green')
    plt.plot(dt, tr3, label='WB', color='blue')
    plt.plot(dt, tr4, label='TVD', color='pink')
    plt.plot(dt, tr5, label='SemiLag', color='purple')
    plt.plot(dt, tr6, label='diff2', color='black')
    plt.axhline(0, linestyle=':', color='black', linewidth=2)
    plt.xlabel('$dt$')
    plt.ylabel('$truncation error$')
    plt.legend(bbox_to_anchor=(1.0,1.0))
    plt.savefig('Plots/' + func.__name__ + 'aa.pdf') 
        
        
#def calculate the  average truncation error    
def trun(phi1, phi2):
    
    nx = len(phi1)
    t = 0.
    
    for i in range(0,nx):
        t = t + np.abs((phi2[i]-phi1[i]))
     
    T = t/nx
    
    return T    
    

    


    
