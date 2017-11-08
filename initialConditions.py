import numpy as np

# Initial conditions function for diffusion

def squareWave(x,alpha,beta):
    """
    id 25825273
    id 25803263
    """
    
    phi = np.zeros_like(x)
    
    # The grid spacing (assumed uniform)
    dx = x[1] - x[0]
    
    # Set phi away from the end points (assume zero at the end points)
    for j in range(0,len(x)-2):
        # edges of the grid box (using west and east notation)
        xw = x[j] - 0.5*dx
        xe = x[j] + 0.5*dx
        
        #integral quantity of phi
        phi[j] = max((min(beta, xe) - max(alpha, xw))/dx, 0)
    phi[len(x)-1]=phi[0]
    return phi


def cosine(x, beta, alpha):
    "A wave as a function of position, x, which is 0 for x > alpha"

    "id 25803263"
    
    
    phi = np.zeros_like(x)
    
    # Set phi away from the end points (assume zero at the end points)
    for j in range(0,len(x)-2):
        if (beta < x[j] < alpha):
            phi[j] = 0.5*(1 - np.cos(4*np.pi*x[j]))
    phi[len(x)-1]=phi[0]
    return phi






def naive(x,alpha,beta):
    phi=np.zeros_like(x)
    dx=x[1]-x[0]
    west_limit=int(alpha/dx)
    east_limit=int(beta/dx)
    for j in range(west_limit+1,east_limit+1):
        phi[j]=1.
    return phi
