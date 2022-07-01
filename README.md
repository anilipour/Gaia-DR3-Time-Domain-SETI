# Gaia DR3: SETI Ellipsoid
## Exploring the SETI Ellipsoid with Gaia DR3 Data



Retrieving Data
-------
Four sets of data can be downloaded with `ellipsoid.py`: all Gaia stars, all Gaia variable stars, all stars in the Gaia Catalogue of Nearby Stars (GCNS), and all GCNS variable stars

Run
``` 
python ellipsoid.py
python ellipsoid.py --num 6000000 --ns 0 --u {Gaia username} --p {Gaia password}
```
to get the files
```
../GCNS_var.fits
../Gaia_var.fits
```
with all the variable stars in GCNS and Gaia. A Gaia account is required to download more than 3 million rows at a time. By default, the number of rows retrieved (```--num```) is set to 100,000. Although Gaia DR3 contains nearly 12 million variable stars, the number of those with reported distances and epoch photometry is slightly less than 6 million. 

To retrieve variable and non-variable stars, add the argument ```--v 0```. These will be saved to the files
```
../GCNS.fits
../Gaia.fits
```

And, to include Gaia variable star classifications, add the argument ```--c 1```. These will be saved to the files
```
../GCNS_var_class.fits
../Gaia_var_class.fits
```

In addition, `ellipsoid.py` contains several useful functions for SETI ellipsoid calculations, light curve plotting, and Gaia querying.


Light Curve Analysis
--------
The Lomb-Scargle Periodogram python notebook contains a short tutorial on the Astropy implementation of the Lomb-Scargle periodogram for finding periodicity in unevenly sampled data. Included is an analysis of a short-timescale variable star light curve from Gaia.

The Ellipsoid Lightcurves python notebook is used to plot the Gaia light curves of stars that have crossed the SETI Ellipsoid during Gaia's observation period. In addition to the light curves in all three bands, the ellipsoid crossing time and its error are displayed as vertical lines.

**to add**
- Other light curve analyses (such as double gaussian for eclipsing binary)


SETI Ellipsoid Analysis
--------
The SN1987A Ellipsoid python notebooks plot the SETI Ellipsoid, using SN 1987A as the conspicuous event, for subsets of the Gaia data that are retrieved as described above. The notebooks also contain various analyses performed on the stars that are in or on, or those that have recently crossed, the ellipsoid.


Dependencies
---------
* Astropy
* Scipy
* NumPy
* Astroquery
* Matplotlib 
* Scipy
* [cubehelix](https://github.com/jradavenport/cubehelix)


