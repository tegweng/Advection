#!/usr/bin/python

# Outer code for setting up the advection problem on a uniform
# grid and calling the function to perform the advection and plot.

from __future__ import absolute_import, division, print_function

import matplotlib.pyplot as plt

# read in all the linear advection schemes, initial conditions and other
# code associated with this application (use with runfile if exec not supported)
"""
execfile("advectionBTBS.py")
execfile("Advection_FTBS.py")
execfile("advectionFTCS.py")
execfile("advectionCTCS.py")
execfile("diagnostics.py")
execfile("initialConditions.py")
execfile("advectionTVD.py")

"""
runfile("advectionBTBS.py")
runfile("Advection_FTBS.py")
runfile("advectionFTCS.py")
runfile("advectionCTCS.py")
runfile("diagnostics.py")
runfile("initialConditions.py")
runfile("advectionTVD.py")

def main(xmin = 0., xmax = 1., nx = 41, T = 0.125, nt = 40, u = 1, \
         squareWaveMin = 0.0, squareWaveMax = 0.5, scheme = BTBS, \
         func = cosine, name_fig='attempt'):
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
    dx = (xmax - xmin)/(nx-1)
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
    #phiOld = squareWave(x, squareWaveMin, squareWaveMax)
    phiOld = func(x,squareWaveMin, squareWaveMax)
    
    # analytic solution of the advection equation (in domain [0,1) )
    #using modulo to keep the solution in the domain
   
    phiAnalytic = func(x - u * T,squareWaveMin,squareWaveMax)
    
    # diffusion using various diffusion schemes
    phiScheme = scheme(phiOld.copy(), c, nt)
    
    # calculate and print out error norms
    L2errScheme = L2ErrorNorm(phiScheme, phiAnalytic)
    print(scheme.__name__ + "L2 error norm = ", L2errScheme)
    
    
    # plot the solutions
    font = {'size' : 10}
    plt.rc('font', **font)
    plt.figure(1)
    plt.clf()
    plt.ion()
    plt.plot(x, phiOld, label='Initial', color='black')
    plt.plot(x, phiAnalytic, label='Analytic', color='black', linestyle='--', \
             linewidth=2)
    plt.plot(x, phiScheme, label=scheme.__name__, color='blue')
    plt.axhline(0, linestyle=':', color='black')
    plt.xlim([0,1])
    plt.ylim([0,2])
    plt.legend()
    plt.xlabel('$x$')
    plt.title("dt = {:.5f}, c = {:.3f}".format(dt, c))
    plt.savefig(name_fig + '(c=' + str(c) +')' + scheme.__name__+ '_' + \
                func.__name__ + '.pdf')
    """
    # plot the errors
    plt.figure(2)
    plt.clf()
    plt.ion()
    
    
    # defining the error vectors that are used in the graph and for evaluating
    # the extremes (m) on the y-axis for plotting
    errorFTCS = phiFTCS - phiAnalytic
    errorBTCS = phiBTCS - phiAnalytic
    m = max(abs(errorFTCS).max(), abs(errorBTCS).max())
    
    plt.plot(x, errorFTCS, label='Error FTCS', color='blue')
    plt.plot(x, errorBTCS, label='Error BTCS', color='red')
    plt.axhline(0, linestyle=':', color='black')
    plt.ylim([-m,m])
    plt.legend()
    plt.xlabel('$x$')
    plt.title("t = {:.2f}, d = {:.2f}".format(nt*dt, d))
    plt.savefig('Plots/' + name_fig + '(t=' + str(int(nt*dt)) + ')_errors.pdf')
    """
    
    return x, phiScheme, phiAnalytic


#def nrms_error_graph(N,d_fixed):
    
    """
    Code for studying the order of convergence of the FTCS and BTCS schemes.
    It must be used with main() working with fixed d and changeable dt,
    preferably commenting the code for plotting the results of main, in order
    to avoid the code to generate 2*N graphs, that is two for each iteration.
    """
    
    # initializing the vector which stores the L2 Error Norm values
    #vector = np.zeros((N,3))
    #for it in range(N):
    #    vector[it,:] = main(nx=21+it*50,d=d_fixed)
    
    # defining the vectors ued for calculating the order of convergence for both
    # FTCS(n) and BTCS(m) schemes by computing the slope between points
    #n = np.zeros(N-1)
    #for i in range(N-1):
    #    n[i] = (np.log(vector[i+1,1]) - np.log(vector[i,1]))/\
    #            (np.log(vector[i+1,0]) - np.log(vector[i,0]))
    
    #m = np.zeros(N-1)
    #for i in range(N-1):
    #    m[i] = (np.log(vector[i+1,2]) - np.log(vector[i,2]))/\
    #            (np.log(vector[i+1,0]) - np.log(vector[i,0]))
    
    # printing out the result by averaging the values stored in n, m
    #print("Order of convergence: FTCS = {:.2f}, BTCS = {:.2f}".format(n.mean()\
    #      ,m.mean()))
    
    # plotting the points onto a log-log graph
    #font = {'size' : 10}
    #plt.rc('font', **font)
    #plt.figure(3)
    #plt.clf()
    #plt.ion()
    #plt.loglog(vector[:,0].transpose(), vector[:,1].transpose(),'bx',\
    #           label='FTCS (mean slope = {:.2f})'.format(n.mean()))
    #plt.loglog(vector[:,0].transpose(), vector[:,2].transpose(),'ro',\
    #           label='BTCS (mean.slope = {:.2f})'.format(m.mean()))
    #plt.loglog(vector[:,0].transpose(),8*(vector[:,0].transpose())**2)
    #plt.legend()
    #plt.xlabel('$\Delta x$')
    #plt.ylabel('$\ell_2$ Norm Error')
    #plt.title('Log-log plot of $\ell_2$ Norm Error with d={:.1f}'.format(\
    #          d_fixed))
    #plt.savefig('Plots/L2errorPlot{}.pdf'.format(int(d_fixed*10)))
