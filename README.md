# LLE-ALPs
 LLE distance limits for ALPs
 
 [MC Oct 22, 2020]:
 This directory contains the following subdirectories:
 - data:
 Data directory has all the relevant spectral, response, and background files needed for the sensitivity study. .bak files correspond to the LAT-LLE background files, .rsp are response files, and .pha are the GRB spectral files. 
Within the data directory, there are three additional subdirectories "high_background", "medium_background", and "low_background". These subdirectories contain the simulated data (.fak files, numbers within the .fak file name correspond to the normalization from N = np.logspace(8.4e-6, 8.4e2, 30)) for the coarse grid. Fine grid numbering is the same, for the normalizations linearly distributed between the two critical normalizations found with the coarse grid. 
The notebooks within the data directory are jupyter representation of python scripts in /scripts.

1_simulating_spectra.ipynb: using the XSPEC fakeit
3_ALP_fitting.ipynb: fitting the ALP model to the simulated data
4_Zero_fitting.ipynb: fitting the Zero model to the simulated data
5_difference.ipynb: LLR calculation

notebooks 2_simulated_spectra_plots.ipynb and 6_Distance-g_medium.ipynb are playground notebooks, and are just checking different spectra and results.

- scripts:
python scripts for running the sensitivity analysis (script adaptation of the notebooks above)

- figures:
some of the figures that are in the paper [will be updated]
