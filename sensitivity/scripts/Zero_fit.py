#!/usr/bin/env python                                                          #
#                                                                              #
# Autor: Milena Crnogorcevic                                                   #
# Date created: 7/22/2019                                                      #
# Last edit: 10/20/2020                                                        #
# Use: Fitting the Zero model to the simulated spectra.                         #
#------------------------------------------------------------------------------#
# Example to run the script (adjust the paths to the SED, energy, bkg, fak, and rsp):
# python Zero_fit.py

# Importing the necessary packages
import numpy as np
import argparse
from scipy.interpolate import interp1d

# importing XSPEC
import xspec as xs
xs.Xset.allowPrompting = False #keeps pyxspec from hanging, waiting a response to a prompt
xs.Xset.allowNewAttributes = True

def zero(engs, params, flux):
    energy_binsizes = np.ediff1d(engs)
    flux[:] = energy_binsizes * 1e-150

zeroInfo = ()
xs.AllModels.addPyMod(zero, zeroInfo, 'add',  spectrumDependent=False)

normZERO = np.zeros((30,2001))
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
        m=xs.Model("zero")
        xs.Fit.statMethod = "pgstat"
        xs.Fit.perform()
        chi_stat[j,i] = xs.Fit.testStatistic
        Fit_st[j,i] = xs.Fit.statistic

np.save('TS_pgfit_ZERO.npy', chi_stat)
np.save('pgfit_ZERO.npy', Fit_st)     
