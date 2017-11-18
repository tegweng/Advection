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


















