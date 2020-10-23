#!/usr/bin/env python                                                          #
#                                                                              #
# Autor: Milena C.                                                             #
# Date: 7/22/2019                                                              #
# Use: Producing ALP text files to be passed to HEASOFT to produce XSPEC       #
#      ALP spectra with a given input normalization and energy range.          #
#------------------------------------------------------------------------------#

import os
import sys
import argparse
import numpy as np

# Defining the strings for spectral files
bkg = 'bn121225417_LAT-LLE_bkgspectra.bak'

# Defining the energy range, as determined by the LLE background
from astropy.io import fits
hdu = fits.open(bkg)
energy = hdu[2].data # reported in keV
E_MeV = (energy['E_MIN'] + energy['E_MAX'])/2000. # average bin energy in MeV


# Calculating the ALP stuff
sys.path.append('/Users/milena/Desktop/2nd_yr_project/ALPs/ALPs_analysis/') # path to calc_alp_signal script
from scipy import integrate
from calc_alp_signal import ALPSNSignal

alp_sn18 = ALPSNSignal(Mprog = 18.)
ts = np.linspace(0.,25.,len(E_MeV))
dndedt_alp18 = alp_sn18(E_MeV, ts, g10=0.1)
SED18 = integrate.simps(dndedt_alp18, ts)/ts.max()

np.save('bn121225417_SED18.npy', SED18)
np.save('bn121225417_EMeV.npy', E_MeV)
