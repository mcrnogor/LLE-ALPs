#!/usr/bin/env python                                                          #
#                                                                              #
# Autor: Milena Crnogorcevic                                                   #
# Date created: 7/22/2019                                                      #
# Last edit: 10/20/2020                                                        #
# Use: Fitting the ALP model to the simulated spectra.                         #
#------------------------------------------------------------------------------#
# Example to run the script (adjust the paths to the SED, energy, bkg, fak, and rsp):
# python ALP_fit.py

# Importing the necessary packages
import numpy as np
import argparse
from scipy.interpolate import interp1d

# importing XSPEC
import xspec as xs
xs.Xset.allowPrompting = False #keeps pyxspec from hanging, waiting a response to a prompt
xs.Xset.allowNewAttributes = True


__description__ = 'Fit ALP model for different normalization values for a given response and background.'

N = np.logspace(np.log10(8.4e-8), np.log10(8.4e2), 30) # normalization values, [1e-52 cm^-2]

E_LLE = np.load(args.emev) # LLE energy bins, [MeV]
SED = np.load(args.sed) # the output of Manuel's code, N = 1 cm^-2

f = interp1d(E_LLE*1000, SED/1000, 'cubic', fill_value="extrapolate") # in keV and per keV

# defining Xspec model
def ALP(engs, params, flux):
    energy_binsizes = np.ediff1d(engs)
    # finding the average bin value for energy
    energy = (np.asarray(engs[:-1])+np.asarray(engs[1:]))/2
    flux[:] = f(energy) * energy_binsizes

ALPInfo = ()
xs.AllModels.addPyMod(ALP, ALPInfo, 'add', spectrumDependent=False)


normALP = np.zeros((30,2001))
chi_stat = np.zeros((30,2001))
Fit_st = np.zeros((30,2001))

xs.AllData.clear()   # clear all data, if any.

for j in range(30):
    for i in range(1,2001):
        xs.AllData.clear()
        s=xs.Spectrum("bn150416773_fakeit_%s.fak{%s}" %(j,i))
        s.response = "bn150416773_LAT-LLE_weightedrsp.rsp"
        s.background= ("bn150416773_LAT-LLE_bkgspectra.bak{1}")
        xs.AllModels.clear()
        m=xs.Model("ALP")
        xs.Fit.statMethod = "pgstat"
        #List of value floats [val,delta,min,bot,top,max].
        #m.ALP.norm.values=[N[j], N[j]*0.1, N[j]*0.1, N[j]*0.95, N[j]*1.05, N[j]*1.95]
        m.ALP.norm.frozen = False
        xs.Fit.perform()
        normALP[j,i] = m.ALP.norm.values[0]
        chi_stat[j,i] = xs.Fit.testStatistic
        Fit_st[j,i] = xs.Fit.statistic

    np.save('TS_pgfit_ALP.npy', chi_stat)
    np.save('pgfit_ALP.npy', Fit_st)
    np.save('norm_ALP.npy', normALP)
