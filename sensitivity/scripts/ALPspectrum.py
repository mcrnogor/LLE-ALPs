#!/usr/bin/env python                                                          #
#                                                                              #
# Autor: Milena Crnogorcevic                                                   #
# Date created: 7/22/2019                                                      #
# Last edit: 10/20/2020                                                        #
# Use: Producing the ALP spectra for a given energy range, duration, and       #
# progenitor mass.                                                             #
#------------------------------------------------------------------------------#
# Example to run the script:
# python ALPspectrum.py -m 10 -b ../data/bn081024891_LAT-LLE_bkgspectra.bak --duration 25. --ALPcode /Users/milena/Desktop/2nd_yr_project/ALPs/ALPs_analysis/ -show

# importing general packages
import os
import sys
import argparse
import numpy as np
import logging
from astropy.io import fits
from scipy import integrate


__description__ = 'Produce ALP spectra for different progenitor masses in the LLE GRB energy range.'

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__,
                                 formatter_class=formatter)
PARSER.add_argument('-m', '--mass', type=int, required=True, choices=[10, 18],
                    help='input the progenitor mass')
PARSER.add_argument('-b', '--bkg', type=str, required=True,
                    help='input the .bak file of the GRB observation of interest')
PARSER.add_argument('--duration', type=float, default = 25.,
                    help='input the time in seconds over which the ALP signal should be averaged')
PARSER.add_argument('--ALPcode', type=str, required=True,
                    help='input the path to the calc_alp_signal script')
PARSER.add_argument('-s', '--show', type=bool, required=False, default=False,
                    help='show the plot of the spectrum')

def ALP_spectrum(**kwargs):
    # Defining the strings for spectral files
    mass = kwargs['mass']
    bkg = kwargs['bkg']
    t_max = kwargs['duration']
    ALPscript = kwargs['ALPcode']

    # Defining the energy range of interest, as determined by the LLE background
    hdu = fits.open(bkg)
    energy = hdu[2].data # reported in keV
    E_MeV = (energy['E_MIN'] + energy['E_MAX'])/2000. # average bin energy in MeV

    #Calculating the ALP spectrum

    sys.path.append(ALPscript) # path to calc_alp_signal script
    from calc_alp_signal import ALPSNSignal

    alp_sn = ALPSNSignal(Mprog = mass)
    ts = np.linspace(0.,t_max,len(E_MeV))
    dndedt_alp = alp_sn(E_MeV, ts, g10=0.1)
    SED = integrate.simps(dndedt_alp, ts)/ts.max()
    logging.info('Saving the ALP spectrum and the energy range...')
    np.save('SED_%s.npy'%mass, SED)
    np.save('EMeV.npy', E_MeV)

    if kwargs['show']:
        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot(E_MeV, SED)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        plt.xscale('log')
        plt.xlim(5,1000)
        plt.xlabel('Energy, [MeV]', fontsize=16)
        plt.ylabel('dN/(dE $\Delta$t) [$10^{52}$ counts s$^{-1}$ MeV$^{-1}$]', fontsize=16)
        plt.title('%s-solar-mass progenitor, averaged over %s s'%(mass,t_max))
        plt.show()

if __name__ == '__main__':
    args = PARSER.parse_args()
    ALP_spectrum(**args.__dict__)
