"""
author for all the differences between this file and advection.py
id 25825273
"""

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
execfile("Artificial_diffusion.py")
execfile("advectionTVD.py")
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


def main(xmin = 0., xmax = 1., nx = 41, T = 0.125, c=0.125, u = 1, d_2=0.1012, \
         d_4=0.05, squareWaveMin = 0.0, squareWaveMax = 0.5, nt=40, Q1=0, Q2=0,\
         func = squareWave, name_fig='attempt', limiter = "LaxWend", n=40):
    """
    Advect an initial function between squareWaveMin and squareWaveMax on a 
    domain between x = xmin and x = xmax split over nx spatial steps 
    with Courant number c and time step dt for nt time steps.
    There are two separate initial conditions which are contained in the file
    initialConditions.py. These can be commented out in the code below as
    required.
    Passing Q1 or Q2 to main selects which kind of operations are requested:
        Q1= True is for plotting numerical solutions of the schemes
        Q2= True is for plotting total variation of the schemes
        Q1 and Q2 = False is for the function to return l2 norm error values
                    for all of the schemes
    """
    
    
    if (Q1==True or Q2==True):
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
    
    if (Q1==False and Q2==False):
        # default parameters set in the function arguments
        nx=n
        nt=n
        dx = (xmax - xmin)/n
        T = c
        dt = T/n
        
        print("Courant number = ", u*dt/dx)
        print("dx = ", dx, " dt = ", dt, " nt = ", n)
        print("end time = ", nt*dt)
    
    # initialising spatial grid
    x = np.zeros(nx)
    for j in range(nx):
        x[j] = xmin + j*dx
    
    # initialising temporal grid (for monotonicity study)
    t = np.zeros(nt)
    for j in range(nt):
        t[j] = 0. + j*dt
    #print('x = ', x)
    
    #initial conditions (Each line is a different initial condition)
    
    phiOld = func(x,squareWaveMin, squareWaveMax)
    
    # analytic solution of the advection equation (in domain [0,1) )
    # using modulo to keep the solution in the domain
    
    phiAnalytic = func((x - u * T)%1.,squareWaveMin,squareWaveMax)
    
    
    # advection using various diffusion schemes (for Artificial_diffusion we are
    # interested only in the first item of the tuple returned)
    phiTVD = TVD(phiOld.copy(), c, nt, u, limiter)
    phiArt_diff2 = Artificial_diffusion(phiOld.copy(), c, nt, dx, dt, d_2, 2)[0]
    phiArt_diff4 = Artificial_diffusion(phiOld.copy(), c, nt, dx, dt, d_4, 4)[0]
    phiSemiLag = SemiLag(phiOld.copy(), c, nt, u, dt)
    phiWB = WB(phiOld.copy(), c, nt)
    phiBTBS = BTBS(phiOld.copy(), c, nt)
    phiCTCS, PhiCTCS_TV = CTCS(phiOld.copy(), c, nt)
    phiFTCS = FTCS(phiOld.copy(), c, nt)
    phiFTBS, PhiFTBS_TV = FTBS(phiOld.copy(), c, nt)
    
    ArtDiff2_TV=Artificial_diffusion(phiOld.copy(), c, nt, dx, dt, d_2, 2)[2]
    ArtDiff4_TV=Artificial_diffusion(phiOld.copy(), c, nt, dx, dt, d_4, 4)[2]
    
    #values of artificial diffusion coefficients (dimensional)
    #print("k_2 is ", d*dx**2/dt)
    #print("k_4 is ", d*dx**4/dt)
    
    font = {'size' : 15}
    plt.rc('font', **font)
        
    if (Q1==True):
        # plot the solutions of PART II finite difference schemes 
        plt.figure(1)
        plt.clf()
        plt.ion()
        plt.plot(x, phiOld, label='Initial', color='black')
        plt.plot(x, phiAnalytic, label='Analytic', color='black',\
                 linestyle='--', linewidth=2)
        plt.plot(x, phiSemiLag, label = 'Semi Lag', color = 'maroon')
        #plt.plot(x, phiTVD, label='TVD', color='blue')
        
        plt.plot(x, phiCTCS, label='CTCS', color='blue')
        plt.plot(x, phiFTCS, label = 'FTCS', color = 'orange')
        
        #plt.plot(x, phiArt_diff2, label='2nd-order diff.', color='red')
        
        #plt.plot(x, phiArt_diff4, label = '4th-order diff.', color = 'green')
        
        plt.axhline(0, linestyle=':', color='black')
        plt.xlim([0,1])
        plt.ylim([-1,2])
        plt.legend(loc=0)
        plt.xlabel('$x$')
        plt.ylabel('$\phi$')
        plt.title("dt = {:.5f}, c = {:.3f}".format(dt, c))
        plt.savefig('Plots/' + name_fig + '_' + \
                    func.__name__ + '.pdf')
        break
        
    if (Q2==True):
        plt.figure(2)
        plt.clf()
        plt.ion()
        plt.plot(t, PhiFTBS_TV, label='FTBS.', color='orange')
        plt.plot(t, PhiCTCS_TV, label='CTCS.', color='blue')
        plt.plot(t, ArtDiff2_TV, label='2nd-order diff.', color='red')
        plt.plot(t, ArtDiff4_TV, label = '4th-order diff.', color = 'green')
        plt.xlim([0,T])
        plt.title('Total Variation for nt='+str(nt)+' time steps dt = '\
                  + str(dt)+'s')
        plt.xlabel('t [s]')
        plt.ylabel('Total Variation')
        plt.legend(loc='center right')
        plt.savefig('Plots/' + name_fig + '_' + \
                func.__name__ + '.pdf')
        break
    
    if (Q1==False or Q2==False):    
        # calculate error norms
        L2errFTBS = L2ErrorNorm(phiFTBS, phiAnalytic)
        L2errCTCS = L2ErrorNorm(phiCTCS, phiAnalytic)
        L2errSemiLag = L2ErrorNorm(phiSemiLag, phiAnalytic)
        L2errWB = L2ErrorNorm(phiWB, phiAnalytic)
        L2errArt_Diff2 = L2ErrorNorm(phiArt_diff2, phiAnalytic)
        L2errArt_Diff4 = L2ErrorNorm(phiArt_diff4, phiAnalytic)
        L2errTVD = L2ErrorNorm(phiTVD, phiAnalytic)
        return dx, L2errFTBS, L2errCTCS, L2errSemiLag, L2errWB,\
    L2errArt_Diff2, L2errArt_Diff4, L2errTVD

def d_values(xmin = 0., xmax = 1., nx = 41, T = 0.125, c=0.125, u = 1,\
         squareWaveMin = 0.0, squareWaveMax = 0.5, nt=40,\
         func = squareWave, N=1000, dmin=0.0005, orderAD=2, graph=True):
    """
    Function to estimate the best value for non-dimensional diffusion to be
    utilised for general plots in main
    """
    
    dt = T/nt                  # time step imposing correct time length T
    dx = (xmax - xmin)/nx
    c = dt*u/dx                #calculating c now we have dt
   
    print("Courant number = ", c)
    print("dx = ", dx, " dt = ", dt, " nt = ", nt)
    print("end time = ", nt*dt)
    """
    """
    # spatial points for plotting and for defining initial conditions
    x = np.zeros(nx)
    for j in range(nx):
        x[j] = xmin + j*dx
    #print('x = ', x)
    
    #initial conditions (Each line is a different initial condition)
    
    phiOld = func(x,squareWaveMin, squareWaveMax)
    
    # analytic solution of the advection equation (in domain [0,1) )
    # using modulo to keep the solution in the domain
   
    phiAnalytic = func((x - u*T)%1.0,squareWaveMin,squareWaveMax)
    
    diffusion=np.zeros(N)
    error=np.zeros(N)
    for i in range(N):
        d=dmin*i
        diffusion[i]=d
        phiArt_diff = Artificial_diffusion(phiOld.copy(),\
                                           c, nt, dx, dt, d, orderAD)[0]
        error[i]= L2ErrorNorm(phiArt_diff, phiAnalytic)
        
    if graph==True:
        plt.plot(diffusion,error)
        
    return x, diffusion
    

def nrms_error_graph(T = 0.125, c=0.125, u = 1, func = squareWave,\
                     name_fig='attempt', limiter = "LaxWend",\
                     N = 10, d_2 = 0.1012, d_4=0.05):
    
    """
    Code for studying the order of convergence of the schemes used in main(),
    plottin N points on a loglog graph for each scheme under study.
    """
    # number of values returned by main
    length=8
    # initializing the vector which stores the L2 Error Norm values
    vector = np.zeros((N,length))
    for it in range(N):
        vector[it,:] = np.array(main(c = c, d_2=d_2, d_4=d_4, T=0.125,\
              func = func, n = 20 + it*20, Q1 = 0, Q2 = 0))
    
    # defining the matrix used to calculate the order of convergence for all of
    # the schemes by computing the slope between points
    m = np.zeros((N-1,length-1))
    for it in range(length-1):
        for i in range(N-1):
            m[i,it] = (np.log(vector[i+1,it+1]) - np.log(vector[i,it+1]))/\
                    (np.log(vector[i+1,0]) - np.log(vector[i,0]))
    
    #plotting the points onto a log-log graph
    font = {'size' : 10}
    plt.rc('font', **font)
    plt.figure('Error graph')
    plt.clf()
    plt.ion()
    plt.loglog(vector[:,0].transpose(), vector[:,1].transpose(),'x',\
               label='FTBS (slope = {:.2f})'.format(m[:,0].mean()))
    plt.loglog(vector[:,0].transpose(), vector[:,2].transpose(),'d',\
               label='CTCS (slope = {:.2f})'.format(m[:,1].mean()))
    plt.loglog(vector[:,0].transpose(), vector[:,3].transpose(),'o',\
               label='SL (slope = {:.2f})'.format(m[:,2].mean()))
    plt.loglog(vector[:,0].transpose(), vector[:,4].transpose(),'v',\
               label='WB (slope = {:.2f})'.format(m[:,3].mean()))
    plt.loglog(vector[:,0].transpose(), vector[:,5].transpose(),'h',\
               label='AD 2nd (slope = {:.2f})'.format(m[:,4].mean()))
    plt.loglog(vector[:,0].transpose(), vector[:,6].transpose(),'s',\
               label='AD 4th (slope = {:.2f})'.format(m[:,5].mean()))
    plt.loglog(vector[:,0].transpose(), vector[:,7].transpose(),'*',\
               label='TVD ' + limiter +\
               ' (slope = {:.2f})'.format(m[:,6].mean()))
    plt.loglog(vector[:,0].transpose(),20*(vector[:,0].transpose())**2)
    plt.legend(loc='lower left', bbox_to_anchor=(0.56,0.01))
    plt.xlabel('$\Delta x$')
    plt.ylabel('$\ell_2$ Norm Error')
    plt.title('Log-log plot of $\ell_2$ Norm Error with $d_2=${}, $d_4=${}'.
              format(d_2,d_4))
    plt.savefig('Plots/OrderAccuracy.pdf')
