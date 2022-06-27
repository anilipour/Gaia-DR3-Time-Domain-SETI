# Gaia-DR3-SETI-Ellipsoid
### Exploring the SETI Ellipsoid with Gaia DR3 Data


Four sets of data can be downloaded with ellipsoid.py: all Gaia stars, all Gaia variable stars, all stars in the Gaia Catalogue of Nearby Stars (GCNS), and all GCNS variable stars

Run
```
python ellipsoid.py
python ellipsoid.py --num 12000000 --ns 0 --u {Gaia username} --p {Gaia password}
```
to get the files
```
../GCNS_var.fits
../Gaia_var.fits
```
with all the variable stars in GCNS and Gaia. A Gaia account is required to download more than 3 million rows at a time. By default, the number of rows retrieved (```num```) is set to 100,000. To retrieve variable and non-variable stars, add the argument ``` --v 0 ```. These will be saved to the files
```
../GCNS.fits
../Gaia.fits
```

