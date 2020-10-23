#!/usr/bin/env python                                                          #
#                                                                              #
# Autor: Milena Crnogorcevic                                                   #
# Date created: 7/22/2019                                                      #
# Last edit: 10/20/2020                                                        #
# Use: Producing the simulated ALP spectra for different normalization values  #
# for a given LLE background and instrument response.                          #
#------------------------------------------------------------------------------#
# Example to run the script:
# python fakeit.py --emev EMeV.npy --sed SED_10.npy -b ../data/bn081024891_LAT-LLE_bkgspectra.bak -r ../data/bn081024891_LAT-LLE_weightedrsp.rsp -show


# Importing the necessary packages
import numpy as np
import argparse
from scipy.interpolate import interp1d

# importing XSPEC
import xspec as xs
xs.Xset.allowPrompting = False #keeps pyxspec from hanging, waiting a response to a prompt
xs.Xset.allowNewAttributes = True


__description__ = 'Produce ALP spectra for different normalization values for a given response and background.'

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__,
                                 formatter_class=formatter)

PARSER.add_argument( '--emev', type=str, required=True,
                    help='input the .npy array pointing to the energy array of the GRB observation of interest')
PARSER.add_argument('--sed', type=str, required=True,
                    help='input the .npy array pointing to the ALP spectrum array for a given mass')
PARSER.add_argument('-b', '--bkg', type=str, required=True,
                    help='input the .bak file of the GRB observation of interest')
PARSER.add_argument('-r', '--rsp', type=str, required=True,
                    help='input the .rsp file of the GRB observation of interest')
PARSER.add_argument('-s', '--show', type=bool, required=False, default=False,
                    help='show an example of the plot of the simulated spectrum')

args = PARSER.parse_args()

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

# Faking the data
def fake_data(**kwargs):
    bkg = kwargs['bkg'] +"{1}"
    rsp = kwargs['rsp']
    
    xs.AllData.clear()
    xs.AllModels.clear()
    m = xs.Model('ALP')

    for i in range(len(N)):
        m.ALP.norm.values = N[i]
        m.ALP.norm.frozen = True
        xs.AllData.clear()
        fs = xs.FakeitSettings(background=bkg,response=rsp, fileName='fake_spec_%s.fak' %i)
        sl = 2000*[fs] # making 2000 simulations
        xs.AllData.fakeit(2000, sl)


    if kwargs['show']:
        import matplotlib.pyplot as plt
        xs.AllData.clear()
        for i in range(1,2001):
            xs.AllData +="fake_spec_20.fak{%i}"%i
        xs.AllModels.clear()
        m=xs.Model("ALP")
        xs.Fit.statMethod = "pgstat"

        xs.Plot.device="/xs"
        xs.Plot.xAxis="MeV"
        xs.Plot.add=True
        xs.Plot.background=True
        xs.Plot.xLog=True
        xs.Plot.yLog=True
        xs.Plot.show()
        xs.Plot("ufspec")

        xerr_uf = np.zeros((2001, len(E_LLE)))
        yerr_uf = np.zeros((2001, len(E_LLE)))
        spec_uf = np.zeros((2001, len(E_LLE)))

        for i in range(1, 2001):
            xerr_uf[i,:] = np.asarray(xs.Plot.xErr(i))
            yerr_uf[i,:] = np.asarray(xs.Plot.yErr(i))
            spec_uf[i,:] = np.asarray(xs.Plot.y(i))
        for i in range(1,2001):
            plt.errorbar(E_LLE, spec_uf[i,:], xerr=xerr_uf[i,:], yerr=yerr_uf[i,:])
        plt.plot(E_LLE, SED*N[20],'r', zorder = 2001, label='Theoretical spectrum')
        plt.xlabel('Energy, [MeV]')
        plt.ylabel('[photons cm$^{-2}$ s$^{-1}$ MeV$^{-1}$] ')
        plt.xscale('log')
        plt.grid()
        plt.legend()
        plt.show()

if __name__ == '__main__':
    args = PARSER.parse_args()
    fake_data(**args.__dict__)
