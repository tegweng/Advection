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
execfile("Warming_and_Beam.py")
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
runfile("Warming_and_Beam.py")
"""

def main(xmin = 0., xmax = 1., nx = 41, T = 0.125, nt = 40, u = 1, d_2=0.1012, \
         d_4 = 0.05, squareWaveMin = 0.0, squareWaveMax = 0.5, \
         func = squareWave, name_fig='attempt', limiter = "Vanleer"):
             
    """
    Advect an initial function between squareWaveMin and squareWaveMax on a 
    domain between x = xmin and x = xmax split over nx spatial steps 
    with Courant number c and time step dt for nt time steps.
    There are two separate initial conditions which are contained in the file
    initialConditions.py. These are passed to the function as arguments using 
    "func = ".
    
    The limiter argument refers to the limiter used for the advection scheme 
    TVD.
        
    All students created the code in this function together.
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
   
    phiAnalytic = func((x - u * T)%1.,squareWaveMin,squareWaveMax)
    


    # advection using various diffusion schemes (for Artificial_diffusion we are
    # interested only in the first item of the tuple returned)
    phiTVDVL = TVD(phiOld.copy(), c, nt, u, limiter)
    phiTVDLW = TVD(phiOld.copy(), c, nt, u, limiter = "LaxWend")
    phiTVDK = TVD(phiOld.copy(), c, nt, u, limiter = "Koren")
    phiArt_diff2 = Artificial_diffusion(phiOld.copy(), c, nt, dx, dt, d_2, 2)[0]
    phiArt_diff4 = Artificial_diffusion(phiOld.copy(), c, nt, dx, dt, d_4, 4)[0]
    phiSemiLag = SemiLag(phiOld.copy(), c, nt, u, dt)
    phiWB = WB(phiOld.copy(), c, nt)
    phiBTBS = BTBS(phiOld.copy(), c, nt)
    phiCTCS = CTCS(phiOld.copy(), c, nt)
    phiFTCS = FTCS(phiOld.copy(), c, nt)
    phiFTBS = FTBS(phiOld.copy(), c, nt)

    #values of artificial diffusion coefficients (dimensional)
    print("k_2 is ", d_2*dx**2/dt)
    print("k_4 is ", d_4*dx**4/dt)
    
    
    # plot the solutions of the linear finite difference schemes 
    font = {'size' : 8}
    plt.rc('font', **font)
    plt.figure(1)
    plt.clf()
    plt.ion()
    plt.axhline(0, linestyle=':', color='gray')
    plt.plot(x, phiOld, label='Initial', color='gray')
    plt.plot(x, phiAnalytic, label='Analytic', color='gray', linestyle='--', \
             linewidth=2)
    """
    plt.plot(x, phiCTCS, label='CTCS', color='blue')
    plt.plot(x, phiFTBS, label='FTBS', color='red')
    plt.plot(x, phiFTCS, label = 'FTCS', color = 'orange')
    plt.plot(x, phiBTBS, label = 'BTBS', color = 'green')
    """
    plt.plot(x, phiTVDK, label='TVD', color='blue')
    plt.plot(x, phiSemiLag, label = 'SemiLag', color = 'green')
    plt.plot(x, phiWB, label='WB', color='maroon')
    plt.plot(x, phiArt_diff2, label='Diffusion (2nd)', color='red')
    plt.plot(x, phiArt_diff4, label = 'Diffusion (4th)', color = 'orange')
    
    plt.xlim([0,1])
    plt.ylim([-1,2])
    plt.legend(borderaxespad=0.)
    plt.xlabel('$x$')
    plt.ylabel('$\phi$')
    plt.title("dt = {:.5f}, c = {:.3f}, ($d_2$ = {:.3f}, $d_4 = {:.3f}$)"\
              .format(dt, c, d_2, d_4))
    plt.savefig('Plots/' + name_fig + '_' + \
                func.__name__ + '.pdf')
    

   #calculate and print out error norms
    print("FTBS L2 error norm = ", L2ErrorNorm(phiFTBS, phiAnalytic))
    print("CTCS L2 error norm = ", L2ErrorNorm(phiCTCS, phiAnalytic))
    print("FTCS L2 error norm = ", L2ErrorNorm(phiFTCS, phiAnalytic))
    print("BTBS L2 error norm = ", L2ErrorNorm(phiBTBS, phiAnalytic))
    print("TVD L2 error norm = ", L2ErrorNorm(phiTVDK, phiAnalytic))
    print("Semi-lagrangian L2 error norm = ", L2ErrorNorm(phiSemiLag, phiAnalytic))
    print("Warming and beam L2 error norm = ", L2ErrorNorm(phiWB, phiAnalytic))
    print("Artificial diffusion L2 error norm = ", 
          L2ErrorNorm(phiArt_diff2, phiAnalytic))
    print("Artificial hyperdiffusion L2 error norm = ", 
          L2ErrorNorm(phiArt_diff4, phiAnalytic))

dxs = [0.00024, 0.0024, 0.024, 0.24]
gradFTBS = [0, 0 , 0, 0]
gradCTCS = [0, 0 , 0, 0]
gradFTCS = [0, 0 , 0, 0]
gradBTBS = [0, 0 , 0, 0]
gradTVDK = [0, 0 , 0, 0]
gradSL = [0, 0 , 0, 0]
gradWB = [0, 0 , 0, 0]
gradART2 = [0, 0 , 0, 0]
gradART4 = [0, 0 , 0, 0]



plt.figure(4)
plt.loglog(dxs,order(dxs,gradFTBS), label = 'FTBS', color = 'blue')
plt.loglog(dxs,order(dxs, gradCTCS), label = 'BTBS', color = 'red')
plt.loglog(dxs,order(dxs, gradFTCS), label = 'FTBS', color = 'green')
plt.loglog(dxs,order(dxs, gradCTCS), label = 'BTBS', color = 'yellow')
plt.loglog(dxs,order(dxs, gradTVDK), label = 'TVD', color = 'orange')
plt.loglog(dxs,order(dxs, gradSL), label = 'SL', color = 'brown')
plt.loglog(dxs,order(dxs, gradWB), label = 'WB', color = 'purple')
plt.loglog(dxs,order(dxs, gradART2), label = 'AD 2', color = 'navy')
plt.loglog(dxs,order(dxs, gradART4), label = 'AD 4', color = 'pink')
plt.legend(bbox_to_anchor=(1, 0.5))
plt.xlabel('$\Delta x$')
plt.savefig('plots/order_accuracy.pdf')