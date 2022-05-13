#------------------------------------------------------------------------------
# © 2021-2022. Triad National Security, LLC. All rights reserved.  This
# program was produced under U.S. Government contract 89233218CNA000001
# for Los Alamos National Laboratory (LANL), which is operated by Triad
# National Security, LLC for the U.S.  Department of Energy/National
# Nuclear Security Administration. All rights in the program are
# reserved by Triad National Security, LLC, and the U.S. Department of
# Energy/National Nuclear Security Administration. The Government is
# granted for itself and others acting on its behalf a nonexclusive,
# paid-up, irrevocable worldwide license in this material to reproduce,
# prepare derivative works, distribute copies to the public, perform
# publicly and display publicly, and to permit others to do so.
#------------------------------------------------------------------------------

from singularity_eos import IdealGas
import numpy as np
import math

SMALL = 1e-20 # to avoid dividing by zero

#class IdealGas(object):
#    def __init__(self, gm1, Cv):
#        pass
#
#    def PressureFromDensityInternalEnergy(self, rho, sie, lmbda):
#        return 0.0
#
#    def BulkModulusFromDensityInternalEnergy(self, rho, sie, lmbda):
#        return 0.0
#
#    def Finalize(self):
#        pass
#
#    def nlambda(self):
#        return 1

def PressureSoundSpeedFromDensityEnergyDensity(rho, uu, eos, # inputs
                                               P, cs, Ncells # outputs
                                               ):
    """
    Loop through a flattened list of cells and fill pressure and sound
    speed using a given EOS. 
    """
    # Allocate lambda object, which is used for caching by the tables.
    # Lambdas can do other things, such as take extra parameters as well.
    #
    # Technically this object should be allocated outside this function
    # and passed in, so that it can be persistent. But it's here for
    # example's sake.
    #
    # For parallel execution, you need one lambda vector per thread.
    # Since we're running in serial, we'll just use one.
    nlambda = eos.nlambda() # get number of elements per lambda
    lmbda = np.zeros(nlambda, dtype=np.double)

    # Loop through the cells and use the two function calls
    for i in range(Ncells):
      sie = uu[i] / (rho[i] + SMALL) # convert to specific internal energy
      P[i] = eos.PressureFromDensityInternalEnergy(rho[i], sie, lmbda)
      bmod = eos.BulkModulusFromDensityInternalEnergy(rho[i], sie, lmbda)
      cs[i] = math.sqrt(bmod / (rho[i] + SMALL))
  

# Parameters for ideal gas
gm1 = 0.6
Cv = 2

# We initialize the eos two ways to show how it can be done
# First initialize the EOS the "easy" way:
eos1 = IdealGas(gm1, Cv)


# If you're on device, you need to call
# eos1.GetOnDevice();
# to call the EOS on device
# However, since we're doing everything on host here, we don't make this call

# Make some arrays
N = 50
rho = np.zeros(N, dtype=np.double)
uu  = np.zeros(N, dtype=np.double)
P   = np.zeros(N, dtype=np.double)
cs  = np.zeros(N, dtype=np.double)

# Here we fill the rho and uu arrays with something sensible.
for i in range(N):
    rho[i] = 1 + 0.1 * math.sin((2 * math.pi * i) / N)
    uu[i] = 1e-2 * rho[i]

# Call it!
PressureSoundSpeedFromDensityEnergyDensity(rho, uu, eos1, P, cs, N)

# And let's print out the final value just for fun
print("The final values are:")
print("rho =", rho[N - 1])
print("uu  =", uu[N - 1])
print("P   =", P[N - 1])
print("cs  =", cs[N - 1])

# It's good practice to call Finalize() after you're done using an EOS.
# This usually only does anything if you're on device, but for GPU-data it
# will call the appropriate destructor.
eos1.Finalize()
