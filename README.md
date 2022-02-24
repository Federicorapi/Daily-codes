# PhD thesis
Codes I write on the go during the PhD at ENS - Paris

Thz folder contains the following scripts:
* BraggCOMSOL.m - To launch a COMSOL simulation of an optical cavity directly from MatLAB. 
* Bragg_0_1.m - To simulate on MatLAB an optical cavity made by Bragg Mirrors and to follow a given mode at different incidence angles
* Effective_index.m - Starting from a given Bragg mirror, the script otpimises each layer thickness in order to minimise the drift in frequency of a given mode
* Reflectivity_plots.m - Compute reflection and transmission of a given material for different thicknesses
* Scan.py - Show the scan obtained from the µPL manip
* Multidiel_Federico.m Computes reflectivity of a Bragg mirror
* n2r.m - converts refractive indices to reflection coefficients of M-layer structure. It is used inside multidiel_Federico.m
* plot_single_spectrum_from_h5.m - Allows to plot one spectrum obtained from the µPL manip
* plot_spectra.m - Once the map is plotted, allows to select the pixel from which we want to plot the spectrum
* scan_h5.m -  Show the scan obtained from the µPL manip
* polarization.m - Plots polarisation measurements 
