#!/usr/bin/env python                                                          #
#                                                                              #
# Autor: Milena Crnogorcevic                                                   #
# Date created: 7/22/2019                                                      #
# Last edit: 10/20/2020                                                        #
# Use: Producing the theoretical ALP spectra shown in the paper.               #
#------------------------------------------------------------------------------#
# Example to run the script:
# python ALPpaper_plots.py --ALPcode /Users/milena/Desktop/2nd_yr_project/ALPs/ALPs_analysis/

# importing general packages
import os
import sys
import argparse
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator, ScalarFormatter

plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.unicode'] = True
plt.rcParams['mathtext.rm'] = 'Times New Roman'
plt.rcParams['mathtext.it'] = 'Times New Roman:italic'
plt.rcParams['mathtext.bf'] = 'Times New Roman:bold'

plt.rc('font', family='serif', size=18)
plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16
plt.rcParams['xtick.major.size'] = 5
plt.rcParams['ytick.major.size'] = 5


__description__ = 'Produce the ALP spectra shown in the paper, with arbitrary energy and time binning.'

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__,
                                 formatter_class=formatter)
PARSER.add_argument('--ALPcode', type=str, required=True,
                    help='input the path to the calc_alp_signal script')

args = PARSER.parse_args()

# calculating the ALP spectrum
sys.path.append(args.ALPcode) # path to calc_alp_signal script
from calc_alp_signal import ALPSNSignal

# defining the energy range
EMeV = np.linspace(1.,1000.,500)
ts = np.linspace(0.,25.,500)
ee, tt = np.meshgrid(EMeV,ts, indexing = 'ij')


# mass = 10 solar masses
alp_sn10 = ALPSNSignal(Mprog = 10.)
dndedt_alp10 = alp_sn10(EMeV, ts, g10=0.1)

# mass = 18 solar masses
alp_sn18 = ALPSNSignal(Mprog = 18.)
dndedt_alp18 = alp_sn18(EMeV, ts, g10=0.1)

# plotting the temporal and energy dependence of the signal for 10-solar-mass progenitor
levels = MaxNLocator(nbins=100).tick_values(0.0, 0.0005)
cmap = plt.cm.get_cmap('magma')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

plt.figure(figsize=(7,7))
c = plt.pcolormesh(ee, tt, dndedt_alp10, cmap=cmap, norm=norm)
plt.colorbar(format='%.1e').set_label(label='dN/(dE dt) $(10^{52}$ $\mathrm{MeV}^{-1} \mathrm{s}^{-1})$',size=20,weight='bold')
plt.gca().set_xscale('log')
plt.gca().set_xlabel('Energy (MeV)', fontsize=20)
plt.gca().set_ylabel('Time (s)', fontsize=20)
plt.xlim(10,300)
plt.ylim(0,25)
plt.annotate('10-M$_{\\odot}$ progenitor', xy=(0.5, 0.93), xycoords='axes fraction', size=20, color='white', weight='bold')
plt.show()

# 18-solar-mass progenitor
plt.figure(figsize=(7,7))
c = plt.pcolormesh(ee, tt, dndedt_alp18, cmap=cmap, norm=norm)
plt.colorbar(format='%.1e').set_label(label='dN/(dE dt) $(10^{52}$ $\mathrm{MeV}^{-1} \mathrm{s}^{-1})$',size=20,weight='bold')
plt.gca().set_xscale('log')
plt.gca().set_xlabel('Energy (MeV)', fontsize=20)
plt.gca().set_ylabel('Time (s)', fontsize=20)
plt.xlim(10,300)
plt.ylim(0,25)
plt.annotate('18-M$_{\\odot}$ progenitor', xy=(0.5, 0.93), xycoords='axes fraction', size=20, color='white', weight='bold')
plt.show()

# spectra for 10- and 18-solar-mass progenitors averaged over 25 seconds

SED10 = integrate.simps(dndedt_alp10, tt, axis = 1)/tt.max() # using the average for 10 solar masses
SED18 = integrate.simps(dndedt_alp18, tt, axis = 1)/tt.max() # using the average for 18 solar masses

SED_Emax10 = EMeV[np.where(SED10 == np.max(SED10))]
SED_Emax18 = EMeV[np.where(SED18 == np.max(SED18))]

plt.figure(figsize=(9,6))
plt.plot(EMeV, SED10, color = '#008080',  label = '10-M$_{\\odot}$ progenitor')
plt.plot(EMeV, SED18, color = 'C3', linestyle = '-',alpha=0.5, label = '18-M$_{\\odot}$ progenitor')
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.xscale('log')
plt.axvline(SED_Emax10, color='k',linestyle = ':', alpha=0.5)
plt.xlim(5,800)
plt.text(24, 0.000, 'E$_{peak}$$\sim$70 MeV ', fontsize=18)
plt.xlabel('Energy (MeV)', fontsize=18)
plt.ylabel('dN/(dE $\Delta$t) ($10^{52}$ counts s$^{-1}$ MeV$^{-1}$ cm$^{-2}$)', fontsize=18)
plt.legend(fontsize=18)
plt.tight_layout()
plt.show()
