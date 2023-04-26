# Signal Synchronization Strategies and Time Domain SETI with Gaia DR3

This is the GitHub repository for Nilipour et al. (2023).

Retrieving Data
-------
The relevant Gaia data can be downloaded with `queryGaia.py`. In particular, the script can easily retrieve all Gaia stars, all Gaia variable stars, all stars in the Gaia Catalogue of Nearby Stars (GCNS), and all GCNS variable stars.

Run
``` 
python queryGaia.py
python queryGaia.py -n 10000000 -ns 0 -usr {Gaia username} -pw {Gaia password}
```
to get the files
```
../GCNS_var.fits
../Gaia_var.fits
```
with all the variable stars in GCNS and Gaia. A Gaia account is required to download more than 3 million rows at a time. By default, the number of rows retrieved (```-n```) is set to 100,000. Although Gaia DR3 contains nearly 12 million variable stars, the number of those with reported distances and epoch photometry is slightly less than 10 million. 

To retrieve variable and non-variable stars, add the argument ```-v 0```. These will be saved to the files
```
../GCNS.fits
../Gaia.fits
```

And, to include Gaia variable star classifications, add the argument ```-c 1```. These will be saved to the files
```
../GCNS_var_class.fits
../Gaia_var_class.fits
```
The methodology for Gaia variable star classification is given by Eyer et al. (2022).

Run 
```
python queryGaia.py -h
```
for a further description of the arguments.

Sample Light Curve Analysis
--------
The Lomb-Scargle Periodogram python notebook contains a short tutorial on the Astropy implementation of the Lomb-Scargle periodogram for finding periodicity in unevenly sampled data. Included is a sample analysis of a short-timescale variable star light curve from Gaia.

The Eclipsing Binary and RR Lyrae python notebooks give examples of how to use the Lomb-Scargle periodogram to find periodicity in variables that are not necessarily sinusoidal in nature. The eclipsing binary notebook also shows how to fit a double Gaussian function to the folded data.


Gaia Light Curve Analysis
--------
The Ellipsoid Lightcurves python notebooks are used to plot the Gaia light curves of stars that have crossed the SETI Ellipsoid during Gaia's observation period. In addition to the light curves in all three bands, the ellipsoid crossing time and its error are displayed as vertical lines. Various parameters before and after crossing time are also compared, including variable periods, median fluxes, and the phase and amplitude of variability. 

The SN Light XCorr python notebook contains a potentially interesting cross-correlation analysis between supernova light curves and Gaia stellar light curves, based on the possibility that an intelligent civilization may send signals copying bright astrophysical events they observe. Although the potential of this technique is tightly constrained by the sparsity of Gaia photometry data, it may be better applied to data from other telescopes or surveys, such as TESS. 

The primary focus of Nilipour et al. (2023) is a novel variability analysis that searches for changes in the parameter of stellar light curves at the SETI Ellipsoid or Seto crossing times. This analysis is contained in the Lightcurves Parameter Comparison notebook, which goes through the full steps of the analysis, from calculating all candidate targets using both methods (which are then saved to the ```stars.csv``` file) to demonstrating the ranking system used to determine the most interesting light curves.


SETI Ellipsoid Analysis
--------
The SN1987A Ellipsoid python notebooks plot the SETI Ellipsoid, using SN 1987A as the conspicuous event, for subsets of the Gaia data that are retrieved as described above. The notebooks also contain various analyses performed on the stars that are in or on, or those that have recently crossed, the ellipsoid.


Seto Method
---------
The Seto notebook implements the search framework proposed by [Seto (2021)](https://iopscience.iop.org/article/10.3847/1538-4357/ac0c7b) with Gaia DR3. The method is explained in the Seto (2021) paper, as well as in Nilipour et al. (2023). An animation is also provided here.

![](https://github.com/anilipour/Gaia-DR3-Time-Domain-SETI/blob/main/Figures/setoAnimation.gif)


Other Python Files
--------
The rest of the ```.py``` files include various functions that are used extensively in the notebooks, from calculating crossing times to plotting light curves. 


Dependencies
---------
* Astropy
* Scipy
* NumPy
* Astroquery
* Matplotlib 
* Scipy
* [cubehelix](https://github.com/jradavenport/cubehelix)

Citations
---------
1. [Davenport, James R. A. et al. "Searching the SETI Ellipsoid with Gaia." (2022).](https://arxiv.org/abs/2206.04092)
2. [Eyer, L. et al. "Gaia Data Release 3. Summary of the variability processing and analysis." (2022).](https://arxiv.org/abs/2206.06416)
3. [Seto, Naoki. "Search for Galactic Civilizations Using Historical Supernovae." (2021)](https://iopscience.iop.org/article/10.3847/1538-4357/ac0c7b)


[![DOI](https://zenodo.org/badge/503060847.svg)](https://zenodo.org/badge/latestdoi/503060847)
