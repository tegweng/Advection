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


# testing Artificial_diffusion
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

# testing CTCS
try:
    CTCS(np.zeros(6), 0.125, 40, 0.05, 0.05, 0.1, 1)
except TypeError:
    pass
else:
    print('Error in CTCS:\
          an error should be raised if orderAD is different from 2 or 4')

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
















