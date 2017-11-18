#!/usr/bin/python

# Outer code for setting up the advection problem on a uniform
# grid and calling the function to perform the advection and plot.

from __future__ import absolute_import, division, print_function

import matplotlib.pyplot as plt

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
execfile("Warming and Beam.py")
"""
runfile("advectionBTBS.py")
runfile("Advection_FTBS.py")
runfile("advectionFTCS.py")
runfile("advectionCTCS.py")
runfile("diagnostics.py")
runfile("initialConditions.py")
runfile("Artificial_diffusion.py")
runfile("advectionTVD.py")
runfile("SemiLagrangian.py")
runfile("Warming and Beam.py")
"""

def main(xmin = 0., xmax = 1., nx = 41, T = 0.125, nt = 40, u = 1, d=0.1, \
         squareWaveMin = 0.0, squareWaveMax = 0.5, \
         func = squareWave, name_fig='attempt'):

    """
    Advect an initial function between squareWaveMin and squareWaveMax on a 
    domain between x = xmin and x = xmax split over nx spatial steps 
    with Courant number c and time step dt for nt time steps.
    There are two separate initial conditions which are contained in the file
    initialConditions.py. These can be commented out in the code below as
    required.
    """
        
    # code for fixed c, nt and nx (modification is needed also in the main 
    # arguments) changing dt so that T remains the same, then calculating the
    # courant number c. All simulations are of the same duration.
    
    # default parameters set in the function arguments
    dx = (xmax - xmin)/nx
    dt = T/nt                  # time step imposing correct time length T
    c = dt*u/dx                #calculating c now we have dt
   
    print("Courant number = ", c)
    print("dx = ", dx, " dt = ", dt, " nt = ", nt)
    print("end time = ", nt*dt)
    
    
    # spatial points for plotting and for defining initial conditions
    x = np.zeros(nx)
    for j in range(nx):
        x[j] = xmin + j*dx
    #print('x = ', x)
    
    #initial conditions (Each line is a different initial condition)
    
    phiOld = func(x,squareWaveMin, squareWaveMax)
    
    # analytic solution of the advection equation (in domain [0,1) )
    # using modulo to keep the solution in the domain
   
    phiAnalytic = func(x - u * T,squareWaveMin,squareWaveMax)
    


    # advection using various diffusion schemes (for Artificial_diffusion we are
    # interested only in the first item of the tuple returned)
    phiTVD = TVD(phiOld.copy(), c, nt, u)
    phiArt_diff2 = Artificial_diffusion(phiOld.copy(), c, nt, dx, dt, d, 2)[0]
    phiArt_diff4 = Artificial_diffusion(phiOld.copy(), c, nt, dx, dt, d, 4)[0]
    phiSemiLag = SemiLag(phiOld.copy(), c, nt, u, dt)
    phiFTCSWB = FTCSWB(phiOld.copy(), c, nt)
    phiBTBS = BTBS(phiOld.copy(), c, nt)
    phiCTCS = CTCS(phiOld.copy(), c, nt)
    phiFTCS = FTCS(phiOld.copy(), c, nt)
    phiFTBS = FTBS(phiOld.copy(), c, nt)
    
    print("d_2 is ", k*dt/dx**2)
    print("d_4 is ", k*dt/dx**4)
    
    
    # plot the solutions of the linear finite difference schemes 
    font = {'size' : 10}
    plt.rc('font', **font)
    plt.figure(1)
    plt.clf()
    plt.ion()
    plt.plot(x, phiOld, label='Initial', color='black')
    plt.plot(x, phiAnalytic, label='Analytic', color='black', linestyle='--', \
             linewidth=2)
    plt.plot(x, phiCTCS, label='CTCS', color='blue')
    plt.plot(x, phiFTBS, label='FTBS', color='red')
    plt.plot(x, phiFTCS, label = 'FTCS', color = 'orange')
    plt.plot(x, phiBTBS, label = 'BTBS', color = 'green')
    plt.axhline(0, linestyle=':', color='black')
    plt.xlim([0,1])
    plt.ylim([-1,2])
    plt.legend()
    plt.xlabel('$x$')
    plt.title("dt = {:.5f}, c = {:.3f}".format(dt, c))
    plt.savefig('Plots/' + name_fig + '_' + \
                func.__name__ + '.pdf')
   
