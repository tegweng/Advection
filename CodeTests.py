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


# testing Artificial_diffusion              id 25825273
try:
    Artificial_diffusion(np.zeros(6), 0.125, 40, 0.05, 0.05, 0.1, 1)
except ValueError:
    pass
else:
    print('Error in Artificial_diffusion:\
          an error should be raised if orderAD is different from 2 or 4')

try:
    Artificial_diffusion(np.zeros(6), 0.125, 40.5, 0.05, 0.05, 0.1, 2)
except TypeError:
    pass
else:
    print('Error in Artificial_diffusion:\
          an error should be raised if nt is not an integer')

try:
    Artificial_diffusion([0,1,2], 0.125, 40, 0.05, 0.05, 0.1, 2)
except TypeError:
    pass
else:
    print('Error in Artificial_diffusion:\
          an error should be raised if phiOld is not a numpy array')
    
try:
    Artificial_diffusion(np.zeros(6), 0.125, 40, -0.05, 0.05, 0.1, 2)
except TypeError:
    pass
else:
    print('Error in Artificial_diffusion:\
          an error should be raised if dx is less than zero')

try:
    Artificial_diffusion(np.zeros(6), 0.125, 40, 0.05, -0.05, 0.1, 2)
except TypeError:
    pass
else:
    print('Error in Artificial_diffusion:\
          an error should be raised if dt is less than zero')

try:
    Artificial_diffusion(np.zeros(6), 0.125, 40, 0.05, 0.05, -0.1, 2)
except TypeError:
    pass
else:
    print('Error in Artificial_diffusion:\
          an error should be raised if d is less than zero')

try:
    Artificial_diffusion(np.zeros(6), -0.125, 40, 0.05, 0.05, 0.1, 2)
except TypeError:
    pass
else:
    print('Error in Artificial_diffusion:\
          an error should be raised if c is less than zero')

# testing CTCS                              id 25825273
try:
    CTCS(np.zeros(6), 0.125, 40.5)
except TypeError:
    pass
else:
    print('Error in CTCS:\
          an error should be raised if nt is not an integer')

try:
    CTCS(np.zeros(6), 0, 40)
except TypeError:
    pass
else:
    print('Error in CTCS:\
          an error should be raised if c is less than or equal to zero')

try:
    CTCS([0,1,2], 0.125, 40)
except TypeError:
    pass
else:
    print('Error in CTCS:\
          an error should be raised if phiOld is not a numpy array')

#testing TVD and flux functions idL 25803263


try:
    TVD(np.zeros(6), 1, 0, 0.5)
except ValueError:
    pass
else:
    print('Error in TVD, an error should be raised if nt<=0')

try:
    TVD(0,1,4, -5)
except TypeError:
    pass
else:
    print('Error in TVD, an error should be raised if phiOld is not a numpy array')
try:
    flux(0,1,4)
except TypeError:
    pass
else:
    print('Error in flux, an error should be raised if phiOld is not a numpy array')

#testing BTBS id:25803263
  
   
try:
    BTBS(np.zeros(6), 1, 0)
except ValueError:
    pass
else:
    print('Error in BTBS, an error should be raised if nt<=0')

try:
    BTBS(0,1,4)
except TypeError:
    pass
else:
    print('Error in BTBS, an error should be raised if phiOld is not a numpy array')


#Testing FTBS ID: 25818629
try:
    FTBS(np.zeros(8), 1, 0)
except ValueError:
    pass
else:
    print('Error in FTBS, error should be raised if nt <= 0')

try:
    FTBS(np.zeros(8), 1, 0.2)
except ValueError:
    pass
else:
    print('Error in FTBS, error should be raised if nt is not an integer')

try:
    FTBS(0, 1, 2)
except TypeError:
    pass
else:
    print('Error in FTBS, error should be raised if phiOld is not a numpy array')


#Testing Semi-Lagrangian ID: 25818269
try:
    SemiLag(np.zeros(8), 1, 0, 1, 0.5)
except ValueError:
    pass
else:
    print('Error in SemiLag, error should be raised if nt <= 0')
    
try:
    SemiLag(np.zeros(8), 1, 0.2, 1, 0.5)
except ValueError:
    pass
else:
    print('Error in SemiLag, error should be raised if nt is not an integer')

try:
    SemiLag(0, 1, 2, 1, 0.5)
except TypeError:
    pass
else:
    print('Error in SemiLag, error should be raised if phiOld is not a numpy \
          array')


#Testing FTCS id:25806676
try:
    FTCS(np.zeros(6), 1, -1)
except ValueError:
    pass
else:
    print('Error in FTCS: error should be raised if nt <= 0')

try:
    FTCS(6, 6, 6)
except TypeError:
    pass
else:
    print('Error in FTCS:\
          error should be raised if phiOld is not a numpy array')

try:
    FTCS(np.zeros(6), 1, 0.6)
except ValueError:
    pass
else:
    print('Error in FTCS: error should be raised if nt is not an integer')  
    
    
#Testing FTCSWB id:25806676
try:
    FTCSWB(np.zeros(8), 1, -1)
except ValueError:
    pass
else:
    print('Error in FTCSWB: error should be raised if nt <= 0')

try:
    FTCSWB(8, 8, 8)
except TypeError:
    pass
else:
    print('Error in FTCSWB:\
          error should be raised if phiOld is not a numpy array')

try:
    FTCSWB(np.zeros(8), 1, 0.6)
except ValueError:
    pass
else:
    print('Error in FTCSWB: error should be raised if nt is not an integer')







