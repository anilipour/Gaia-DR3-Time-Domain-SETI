import numpy as np


from astropy.coordinates import SkyCoord
import astropy.coordinates as coord
from astropy import constants as const
from astropy import units as u
from astropy.time import Time
from astropy.table import Table, QTable
import os
import scipy.stats

from astroquery.gaia import Gaia
import matplotlib.pyplot as plt


def readFile(filename):
    filepath = f'../{filename}.fits'
    stars = Table.read(filepath, format='fits')
    c1 = SkyCoord(ra = stars['ra'],
              dec = stars['dec'],
              distance = stars['dist'],
              frame = 'icrs')

    return c1, stars
    



def onEllipsoid(event, t0,  c1, stars, tol):
    # event is a SkyCoord with RA, Dec, distance of event (i.e. SN 1987A)
    # t0 is the time event was observed on Earth
    # c1 is SkyCoord with RA, Dec, distance of stars
    # stars is astropy Table with necessary info 
    # tol is the tolerance for stars being on the ellipse
    # Returns table of stars currently within (-tol, tol) of being on the ellipse 
    
    t1 = Time.now()
    dt = t1-t0 # time since event was observed on Earth
    
    
    # Ellipsoid geometry
    c = event.distance.to('lyr') / 2 # dist to foci from ellipse center (half the distance from Earth to SN1987A, the two foci)   
    a = (((dt.to('s') * const.c) / 2) + c).to('lyr') # semi-major axis of ellipse
    # first term is the distance along the major axis beyond Earth -- it is half the distance light would travel in the time since the event was observed
     
    
    # Star geometry
    d1 = stars['dist'] # dist to stars
    d2 = c1.separation_3d(event) # dist from all GCNS stars to SN 1987A
    
    
    # Which stars are within some tolerance of being ON the ellipse?
    OYES = np.abs((d1.to('lyr').value + d2.to('lyr').value) - (2 * a.to('lyr').value)) <= tol
    # check if sum of focal point distances minus the major axis diameter is within the interval (-tol, tol)
    
    return c1[OYES], stars[OYES]


def crossEllipsoid(event, t0,  c1, stars, tol, start, end=None):
    # event is a SkyCoord with RA, Dec, distance of event (i.e. SN 1987A)
    # t0 is the time event was observed on Earth
    # c1 is SkyCoord with RA, Dec, distance of stars
    # stars is astropy Table with necessary info 
    # tol is the tolerance for stars crossing the ellipse
    # start is the start time from when we consider stars crossing the ellipse
    # Returns table of stars that have crossed the ellipse since start 
    
    if end==None:
        end = Time.now()
    c = event.distance.to('lyr') / 2
    d1 = stars['dist'] # dist to stars
    d2 = c1.separation_3d(event) # dist from all GCNS stars to SN 1987A
    
    etime = d2.to('lyr') + d1.to('lyr') - (2*c)
    
    gstart = start
    gstart_t0 = (gstart-t0).to('year').value # years between start time and time the event was observed on Earth

    crossings = ((etime.value >= gstart_t0-tol) & # since start (minus tolerance)
      (etime.value < ((end-t0).to('year').value+tol)) # up to (tolerance) years from now
     ) 

    
    return c1[crossings], stars[crossings]


def crossErrorEllipsoid(event, t0,  c1, stars, tol, start=Time('2014-07-25T10:30'), end=Time('2017-05-28T08:44')):
    # event is a SkyCoord with RA, Dec, distance of event (i.e. SN 1987A)
    # t0 is the time event was observed on Earth
    # c1 is SkyCoord with RA, Dec, distance of stars
    # stars is astropy Table with necessary info 
    # tol is the tolerance for stars crossing the ellipse
    # start is the start time from when we consider stars crossing the ellipse
    # end is the end time from whn we consider stars crossing the ellipse
    # The default start and end times are the data collection period of Gaia DR3,
    # as detailed at cosmos.esa.int/web/gaia/dr3
    # Returns table of stars that have crossed the ellipse between start and end,
    # and the error in the crossing time is also between start and end 
    
    c = event.distance.to('lyr') / 2
    d1 = stars['dist'] # dist to stars
    d2 = c1.separation_3d(event) # dist from all GCNS stars to SN 1987A
    
    derror = (((stars['dist84'] - stars['dist']) + (stars['dist'] - stars['dist16']))/2).to('lyr').value
    # distance error to the stars dominates the error in crossing time
    etime = d2.to('lyr') + d1.to('lyr') - (2*c)
    
    gstart_t0 = (start-t0).to('year').value # years between start time and time the event was observed on Earth

    crossings = ((etime.value - derror >= gstart_t0-tol) & # since start (minus tolerance)
      (etime.value + derror < ((end-t0).to('year').value+tol)) # up to (tolerance) years from now
    )

    
    return c1[crossings], stars[crossings]
     

def getLC(lc_table, source_id, band):
    # Returns light curve for a specific band and source id from a light curve
    # table of all sources in all bands
    source_mask = lc_table['source_id'] == source_id
    source_lc = lc_table[source_mask]
    band_mask = source_lc['band'] == band
    band_lc = source_lc[band_mask]

    return band_lc
    
    
def varClass(source_id):
    # Returns the Gaia variable classification for a given source id

    query = f"SELECT best_class_name \
            FROM gaiadr3.vari_classifier_result \
            WHERE gaiadr3.vari_classifier_result.source_id = '{source_id}'"

    job = Gaia.launch_job_async(query, output_format='fits')
    results = job.get_results()
    return results['best_class_name']


def varClasses(source_id_tup):
    # Returns the Gaia variable classifications for a tuple of source ids

    query = f"SELECT best_class_name, source_id \
            FROM gaiadr3.vari_classifier_result \
            WHERE gaiadr3.vari_classifier_result.source_id IN {source_id_tup}"

    job = Gaia.launch_job_async(query, output_format='fits')
    results = job.get_results()
    return results

def xTime(star, c0=None, t0=None):
    # Returns the ellipsoid crossing time of the SkyCoord object(s)
    # if c0 and t0 are None, the ellipsoid defaults to SN 1987A

    if not c0 or not t0:
        # Properties of SN1987A
        t0 = Time({'year': 1987, 'month': 2, 'day': 23}, format='ymdhms')

        c0_radec = SkyCoord.from_name('SN 1987A')

        # Panagia (1999) https://ui.adsabs.harvard.edu/abs/1999IAUS..190..549P/abstract
        d0 = 51.4 * u.kpc
        d0_err = 1.2 * u.kpc

        c0 = SkyCoord(ra=c0_radec.ra, dec=c0_radec.dec, distance=d0)
    
    
    c = c0.distance.to('lyr') / 2
    d1 = star.distance # dist to stars
    d2 = star.separation_3d(c0) # dist from all GCNS stars to SN 1987A
    
    etime = d2.to('lyr') + d1.to('lyr') - (2*c)
    xtime = etime.value*u.year + t0
    return xtime.jd

def medLC(xtime, lc):
    leftY, rightY = [], []

    for key in lc:
        if key < xtime:
            leftY.append(lc[key])
        else:
            rightY.append(lc[key])
    
    medLeft = np.nanmedian(leftY)
    medRight = np.nanmedian(rightY)
    madLeft = scipy.stats.median_absolute_deviation(leftY)
    madRight = scipy.stats.median_absolute_deviation(rightY)
    return medLeft, medRight, madLeft, madRight


def getMedLC(sid, star, lcDict, y='flux', time=None):
    # Calculates the medians before and after crossing time
    
    if time==None:
        xtime = xTime(star)[0]
    else:
        xtime = time
    lcTimes = Time(lcDict[sid]['time'].value + 2455197.5, format='jd')
    
    if y != 'flux':
        lcY = lcDict[sid]['mag']

    else:
        lcY = lcDict[sid]['flux']

    lc = dict(zip(lcTimes.value, lcY.value))
    

    leftMed, rightMed, leftMAD, rightMAD = medLC(xtime, lc)

    return leftMed, rightMed, leftMAD, rightMAD




def compMedLC(sid, star, lcDict, y='flux', time=None):
    leftMed, rightMed, leftMAD, rightMAD, = getMedLC(sid, star, lcDict, y, time)
    
    if rightMed > leftMed + 3*leftMAD or rightMed < leftMed - 3*leftMAD:
        return True
    elif leftMed > rightMed + 3*rightMAD or leftMed < rightMed - 3*rightMAD:      
        return True
    else:
        return False
 


def plotLC(sid, star, star_row, glcDict, bplcDict, rplcDict, y='flux', median=True, normalize=True, time=None):
    # Plots the light curves for a star in the G, BP, and RP bands
    # Needs the source id of the star, the SkyCoord object corresponding to the star,
    # the row of the table corresponding to the star,
    # and three dictionaries containing the source ids of stars as keys and the
    # light curve tables as items

    if time==None:
        xtime = xTime(star)[0]
    else:
        xtime=time
    derror = (((star_row['dist84'] - star_row['dist']) + (star_row['dist'] - star_row['dist16']))/2).to('lyr').value*365.25

    glcTimes = Time(glcDict[sid]['time'].value + 2455197.5, format='jd')
    blcTimes = Time(bplcDict[sid]['time'].value + 2455197.5, format='jd')
    rlcTimes = Time(rplcDict[sid]['time'].value + 2455197.5, format='jd')
    
    if y != 'flux':
        glcY = glcDict[sid]['mag']
        blcY = bplcDict[sid]['mag']
        rlcY = rplcDict[sid]['mag']
        gerr = None
        berr = None
        rerr = None

    else:
        glcY = glcDict[sid]['flux']
        blcY = bplcDict[sid]['flux']
        rlcY = rplcDict[sid]['flux']
        gerr = glcDict[sid]['flux_error']
        berr = bplcDict[sid]['flux_error']
        rerr = rplcDict[sid]['flux_error']

    fig, ax = plt.subplots(3, sharex=True, figsize=[8,6], dpi=150)

    glcMed = np.nanmedian(glcY)
    blcMed = np.nanmedian(blcY)
    rlcMed = np.nanmedian(rlcY)

    glcY /= glcMed
    blcY /= blcMed
    rlcY /= rlcMed

    gerr /= glcMed
    berr /= blcMed
    rerr /= rlcMed
    
    ax[0].errorbar(glcTimes.value, glcY, yerr=gerr, fmt='o', c='green', label='G', ms=3, elinewidth=1)
    ax[1].errorbar(blcTimes.value, blcY, yerr=berr, fmt='o', c='blue', label='BP', ms=3, elinewidth=1)
    ax[2].errorbar(rlcTimes.value, rlcY, yerr=rerr, fmt='o', c='red', label='RP', ms=3, elinewidth=1)


    ax[0].vlines(xtime, ymin=min((glcY-gerr).value), ymax=max((glcY+gerr).value), color='green')
    ax[1].vlines(xtime, ymin=min((blcY-berr).value), ymax=max((blcY+berr).value), color='blue')
    ax[2].vlines(xtime, ymin=min((rlcY-rerr).value), ymax=max((rlcY+rerr).value), color='red')
    
    ax[0].vlines([xtime+derror, xtime-derror],  ymin=min((glcY-gerr).value), ymax=max((glcY+gerr).value), color='green', ls='dashed')
    ax[1].vlines([xtime+derror, xtime-derror],  ymin=min((blcY-berr).value), ymax=max((blcY+berr).value), color='blue', ls='dashed')
    ax[2].vlines([xtime+derror, xtime-derror],  ymin=min((rlcY-rerr).value), ymax=max((rlcY+rerr).value), color='red', ls='dashed')



    ax[0].legend(loc='best')
    ax[1].legend(loc='best')
    ax[2].legend(loc='best')

    ax[2].set_xlabel('BJD')

    if y != 'flux':
        if normalize:
            ax[1].set_ylabel('Normlized Mag')
        else:
            ax[1].set_ylabel('Mag')
    else:
        if normalize:
            ax[1].set_ylabel('Normalized Flux')
        else:
            ax[1].set_ylabel('Flux (e/sec)')

    ax[0].set_title(f'Source {sid}')
    
    if median:
        gLeftMed, gRightMed, gLeftMAD, gRightMAD = getMedLC(sid, star, glcDict, y, time)
        bLeftMed, bRightMed, bLeftMAD, bRightMAD = getMedLC(sid, star, bplcDict, y, time)
        rLeftMed, rRightMed, rLeftMAD, rRightMAD = getMedLC(sid, star, rplcDict, y, time)
        ax[0].hlines(gLeftMed, xmin=min(glcTimes.value), xmax=xtime, color='green', linewidth=0.5)
        ax[0].hlines([gLeftMed+3*gLeftMAD, gLeftMed-3*gLeftMAD], xmin=min(glcTimes.value), xmax=xtime, color='green', linewidth=0.5, ls = 'dashed')
        ax[0].hlines(gRightMed, xmin=xtime, xmax=max(glcTimes.value), color='green', linewidth=0.5)
        ax[0].hlines([gRightMed+3*gRightMAD, gRightMed-3*gRightMAD], xmin=xtime, xmax=max(glcTimes.value), color='green', linewidth=0.5, ls = 'dashed')

        ax[1].hlines(bLeftMed, xmin=min(blcTimes.value), xmax=xtime, color='blue', linewidth=0.5)
        ax[1].hlines([bLeftMed+3*bLeftMAD, bLeftMed-3*bLeftMAD], xmin=min(blcTimes.value), xmax=xtime, color='blue', linewidth=0.5, ls = 'dashed')
        ax[1].hlines(bRightMed, xmin=xtime, xmax=max(blcTimes.value), color='blue', linewidth=0.5)
        ax[1].hlines([bRightMed+3*bRightMAD, bRightMed-3*bRightMAD], xmin=xtime, xmax=max(blcTimes.value), color='blue', linewidth=0.5, ls = 'dashed')

        ax[2].hlines(rLeftMed, xmin=min(rlcTimes.value), xmax=xtime, color='red', linewidth=0.5)
        ax[2].hlines([rLeftMed+3*rLeftMAD, rLeftMed-3*rLeftMAD], xmin=min(rlcTimes.value), xmax=xtime, color='red', linewidth=0.5, ls = 'dashed')
        ax[2].hlines(rRightMed, xmin=xtime, xmax=max(rlcTimes.value), color='red', linewidth=0.5)
        ax[2].hlines([rRightMed+3*rRightMAD, rRightMed-3*rRightMAD], xmin=xtime, xmax=max(rlcTimes.value), color='red', linewidth=0.5, ls = 'dashed')