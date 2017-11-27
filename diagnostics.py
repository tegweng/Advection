# Various function for plotting results and for calculating error measures

### Copy out most of this code. Code commented with 3#s (like this) ###
### is here to help you to learn python and need not be copied      ###

### If you are using Python 2.7 rather than Python 3, import various###
### functions from Python 3 such as to use real number division     ###
### rather than integer division. ie 3/2  = 1.5  rather than 3/2 = 1###
from __future__ import absolute_import, division, print_function

### The numpy package for numerical functions and pi                ###
import numpy as np

import math

#import package for graphing
import matplotlib.pyplot as plt


def L2ErrorNorm(phi, phiExact):
    "Calculates the L2 error norm (RMS error) of phi in comparison to"
    "phiExact, ignoring the boundaries"
    
    #remove one of the end points
    phi = phi[1:-1]
    phiExact = phiExact[1:-1]
    
    # calculate the error and the error norms
    phiError = phi - phiExact
    L2 = np.sqrt(sum(phiError**2)/sum(phiExact**2))

    return L2
    
def disprel(c, dx):
    """
    25803263
    Function to calculate the dispersion relation for CTCS and plot 
    the relation with the w = k. const the actual solution. u = 1.
    """
    k = np.arange(0, int(np.pi / dx)+1, 1)
    kdx = k * dx
    ang = np.zeros(len(k))
    cang = np.zeros(len(k))
    alpha = np.zeros(len(k))
    pmode = np.zeros(len(k))
    cmode = np.zeros(len(k))
    tmode = np.zeros(len(k))
    
    for i in xrange(0,int(np.pi / dx)+1):
        alpha[i] = math.asin(c * math.sin(kdx[i]))
         
        if alpha[i] == 0:
            pmode[i] = 1
        else:
            pmode[i] = alpha[i] / (c * kdx[i])
    
        cmode[i] = - pmode[i]
    
        tmode[i] = kdx[i]
        
        ang[i] = float(kdx[i]) * pmode[i]
        cang[i] = -ang[i]
    
    font = {'size' : 8}
    plt.rc('font', **font)
    plt.figure(1)
    plt.clf()
    plt.ion()
    plt.plot(kdx, pmode, label='Physical', color='black')
    plt.plot(kdx, cmode, label='Computational', color='red')
    plt.axhline(0, linestyle=':', color='black')
    plt.xlim([0,np.pi])
    plt.ylim([-1.,1.])
    plt.xticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi], 
               (0, '$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'))
    plt.legend(borderaxespad=0.)
    plt.xlabel('$k\Delta x$')
    plt.ylabel('$u_n/u$')
    plt.title('dx =' + str(dx) +' c =' + str(c) )
    plt.savefig('Plots/dispersionrelation.pdf')
    
    font = {'size' : 8}
    plt.rc('font', **font)
    plt.figure(2)
    plt.clf()
    plt.ion()
    plt.plot(kdx, ang, label='Physical', color='black')
    plt.plot(kdx, tmode, label='True', color='black', linestyle='--', \
             linewidth=2)
    plt.plot(kdx, cang, label='Computational', color='red')
    plt.axhline(0, linestyle=':', color='black')
    plt.xlim([0,np.pi])
    plt.ylim([-1.5,1.5])
    plt.xticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi], 
               (0, '$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'))
    plt.legend(borderaxespad=0.)
    plt.xlabel('$k\Delta x$')
    plt.ylabel('$\omega_n/u$')
    plt.title('dx =' + str(dx) +' c =' + str(c) )
    plt.savefig('Plots/dispersionrelation2.pdf')
    
    return  ang
    
def compare(c1, c2, dx):
    
    """
    Code to compare the dispersion relations depending on the value of c. 
    Note that you need to comment out the figures in the function disprel 
    to get accurrate results.
    
    25803263
    """
    k = np.arange(0, int(np.pi / dx)+1, 1)
    kdx = k * dx
    
    font = {'size' : 10}
    plt.rc('font', **font)
    plt.figure(3)
    plt.clf()
    plt.ion()
    plt.plot(kdx, disprel(c1,dx), label=str(c1), color='black')
    plt.plot(kdx, disprel(c2,dx), label=str(c2), color='red')
    plt.axhline(0, linestyle=':', color='black')
    plt.xlim([0,np.pi])
    plt.ylim([-1.5,1.5])
    plt.xticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi], 
               (0, '$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'))
    plt.legend(borderaxespad=0.)
    plt.xlabel('$k\Delta x$')
    plt.ylabel('$\omega_n/u$')
    plt.title('c2 =' + str(c2) +' c1 =' + str(c1) )
    plt.savefig('Plots/dispersionrelationcompare.pdf')
    
def order(dxs, errors):
    grad = np.zeros(len(errors)-1)
    for i in xrange(0,len(grad)):
        grad[i]=(math.log(abs(errors[i+1]))-math.log(abs(errors[i])))/(math.log(abs(dxs[i+1]))-math.log(abs(dxs[i])))
    print('The accuracy is' + str(grad))

