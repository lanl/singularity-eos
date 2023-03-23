#!/usr/bin/env python

# © 2021-2023. Triad National Security, LLC. All rights reserved.  This
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

import numpy as np
from units import cgs
try:
  import mpi4py
except ImportError:
  pass
import h5py
import argparse

parser = argparse.ArgumentParser(
  description = "Generate a table of ideal gas EOS to be read by Stellar Collapse readers")
parser.add_argument('-g', '--gamma', default = 1.4, type=float,
                    help = 'adiabatic index')
parser.add_argument('--mu', default = cgs['MP'], type=float,
                    help = 'mean molecular weight of gas particles')
parser.add_argument('--nrho', default = 234, type=int,
                    help = 'Number of points in density')
parser.add_argument('--ntemp', default = 136, type=int,
                    help = 'Number of points in temperature')
parser.add_argument('--nye', default = 50, type=int,
                    help = 'Number of points in electron fraction')
parser.add_argument('--rhomin', default = 1.e-10, type=float,
                    help = 'Minimum density (g/cm^3)')
parser.add_argument('--rhomax', default = 1.e-4, type=float,
                    help = 'Maximum density (g/cm^3)')
parser.add_argument('--tempmin', default = 1.e2, type=float,
                    help = 'Minimum temperature (K)')
parser.add_argument('--tempmax', default = 1.e10, type=float,
                    help = 'Maximum temperature (K)')
parser.add_argument('--yemin', default = 0.0, type=float,
                    help = 'Minimum electron fraction')
parser.add_argument('--yemax', default = 0.55, type=float,
                    help = 'Maximum electron fraction')
parser.add_argument('-o', '--output', default = None, type=str,
                    help = 'Name of output file')
args = parser.parse_args()

# Arguments
filename = args.output
mu = args.mu
gam = args.gamma
rho_min = args.rhomin
rho_max = args.rhomax
n_rho = args.nrho
temp_min = args.tempmin
temp_max = args.tempmax
temp_min_MeV = temp_min*cgs['KBOL']/cgs['MEV']
temp_max_MeV = temp_max*cgs['KBOL']/cgs['MEV']
n_temp = args.ntemp
ye_min = args.yemin
ye_max = args.yemax
n_ye = args.nye
energy_shift = 0.0 # All energies positive for gamma law

# Base grid
lrho = np.linspace(np.log10(rho_min), np.log10(rho_max), n_rho)
ltemp = np.linspace(np.log10(temp_min), np.log10(temp_max), n_temp)
ltemp_MeV = np.linspace(np.log10(temp_min_MeV), np.log10(temp_max_MeV), n_temp)
ye = np.linspace(ye_min, ye_max, n_ye)
rho = np.power(10., lrho)
temp = np.power(10., ltemp)

# 3D grid quantities
YE, T, RHO = np.meshgrid(ye, temp, rho, indexing='ij')
PRS = RHO*cgs['KBOL']*T/mu
UU = PRS/(gam - 1.)
EPS = UU/RHO
EF = RHO*cgs['CL']**2 + UU*gam
CS2 = cgs['CL']**2*gam*PRS/EF
CS = np.sqrt(CS2)
ENT = PRS*RHO**(-gam)
GAMMA = np.ones_like(PRS)
DPDRHOE = (gam - 1.)*EPS
DPDERHO = (gam - 1.)*RHO
DEDT = np.ones_like(PRS)*EPS/T
LP = np.log10(PRS)
LE = np.log10(EPS)
Xa = np.zeros_like(PRS)
Xh = np.zeros_like(PRS)
Xp = YE
Xn = 1. - Xp
Abar = np.ones_like(PRS)
Zbar = np.ones_like(PRS)

# Output
if filename is None:
  print("No filename provided! Exiting...")
  import sys
  sys.exit()
with h5py.File(filename, 'w') as f:
  f.create_dataset('energy_shift', data = np.array([energy_shift]))
  f.create_dataset('pointsrho', data = np.array([n_rho]))
  f.create_dataset('pointstemp', data = np.array([n_temp]))
  f.create_dataset('pointsye', data = np.array([n_ye]))
  f.create_dataset('cs2', data = CS2)
  f.create_dataset('dpdrhoe', data=DPDRHOE)
  f.create_dataset('dpderho', data=DPDERHO)
  f.create_dataset('dedt', data = DEDT*cgs['MEV']/cgs['KBOL'])
  f.create_dataset('entropy', data = ENT)
  f.create_dataset('gamma', data = GAMMA)
  f.create_dataset('logenergy', data = LE)
  f.create_dataset('logpress', data = LP)
  f.create_dataset('logrho', data = lrho)
  f.create_dataset('logtemp', data = ltemp_MeV)
  f.create_dataset('ye', data = ye)
  f.create_dataset("Xa", data = Xa)
  f.create_dataset("Xh", data = Xh)
  f.create_dataset("Xp", data = Xp)
  f.create_dataset("Xn", data = Xn)
  f.create_dataset("Abar", data = Abar)
  f.create_dataset("Zbar", data = Zbar)
