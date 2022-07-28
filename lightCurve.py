import numpy as np

from astropy.time import Time
import scipy.stats

import ellipsoid 

from astroquery.gaia import Gaia
import matplotlib.pyplot as plt


def retrieveLC(source_id):
    # Retrieves light curve of a single source id from Gaia
    query = f"SELECT vari.* \
            FROM gaiadr3.vari_summary as vari \
            WHERE source_id = '{source_id}'"

    job     = Gaia.launch_job_async(query, output_format='fits')
    results = job.get_results() # results is the row of the source_id in vari_summary
    retrieval_type = 'EPOCH_PHOTOMETRY'          # Options are: 'EPOCH_PHOTOMETRY', 'MCMC_GSPPHOT', 'MCMC_MSC', 'XP_SAMPLED', 'XP_CONTINUOUS', 'RVS', 'ALL'
    data_structure = 'INDIVIDUAL'     # Options are: 'INDIVIDUAL', 'COMBINED', 'RAW'
    data_release   = 'Gaia DR3'     # Options are: 'Gaia DR3' (default), 'Gaia DR2'


    datalink  = Gaia.load_data(ids=results['source_id'], data_release = data_release, retrieval_type=retrieval_type, data_structure = data_structure, verbose = False, format='votable', output_file=None)
    dl_key  = [inp for inp in datalink.keys()]
    product = datalink[dl_key[0]][0]
    product_tb = product.to_table()  # Export to Astropy Table object.

    return product_tb
    
def getLC(lc_table, source_id, band):
    # Returns light curve for a specific band and source id from a light curve
    # table of all sources in all bands
    source_mask = lc_table['source_id'] == source_id
    source_lc = lc_table[source_mask]
    band_mask = source_lc['band'] == band
    band_lc = source_lc[band_mask]

    return band_lc

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
        xtime = ellipsoid.xTime(star)[0]
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
 


def plotLC(sid, star, star_row, glcDict, bplcDict, rplcDict, y='flux', median=True, normalize=True, time=None, cross=True):
    # Plots the light curves for a star in the G, BP, and RP bands
    # Needs the source id of the star, the SkyCoord object corresponding to the star,
    # the row of the table corresponding to the star,
    # and three dictionaries containing the source ids of stars as keys and the
    # light curve tables as items
    
    if cross:
        if time==None:
            xtime = ellipsoid.xTime(star)[0]
        else:
            xtime=time
    
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

    if cross:
        derror = (((star_row['dist84'] - star_row['dist']) + (star_row['dist'] - star_row['dist16']))/2).to('lyr').value*365.25

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

    ax[0].set_title(f'Source {sid} Light Curve')
    
    if median and cross:
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