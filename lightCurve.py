import numpy as np

from astropy.time import Time
from astropy.timeseries import LombScargle, TimeSeries
import scipy.stats
from astropy import units as u

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

def comparePeriod(lc, xtime, min_freq=1, max_freq=10):
    lcTimes = Time(lc['time'].value + 2455197.5, format='jd')
    lcFlux = lc['flux']
    lcFerr = lc['flux_error']

    lcTS = TimeSeries(time = lcTimes, data={'flux' : lcFlux/np.median(lcFlux), 'flux_error' : lcFerr/np.median(lcFlux)})
    
    freqs = np.linspace(min_freq/u.d, max_freq/u.d, 10000)
    power = LombScargle(lcTS.time, lcTS['flux'], lcTS['flux_error']).power(freqs)
    best_freq = freqs[np.argmax(power)]


    right_mask = lcTS.time.value > xtime
    left_mask = lcTS.time.value <= xtime

    lcTS_right = TimeSeries(time = lcTimes[right_mask], data={'flux' : lcFlux[right_mask]/np.median(lcFlux[right_mask]), 'flux_error' : lcFerr[right_mask]/np.median(lcFlux[right_mask])})
    lcTS_left = TimeSeries(time = lcTimes[left_mask], data={'flux' : lcFlux[left_mask]/np.median(lcFlux[left_mask]), 'flux_error' : lcFerr[left_mask]/np.median(lcFlux[left_mask])})

    power_right = LombScargle(lcTS_right.time, lcTS_right['flux'], lcTS_right['flux_error']).power(freqs)
    power_left = LombScargle(lcTS_left.time, lcTS_left['flux'], lcTS_left['flux_error']).power(freqs)

    best_freq_right = freqs[np.argmax(power_right)]
    best_freq_left = freqs[np.argmax(power_left)]

    return best_freq, best_freq_right, best_freq_left

def comparePeriodPlot(lcDict, sid, star_row, star, min_freq=0.05, max_freq=10, save=False, savefolder=None):
    xtime = ellipsoid.xTime(star)[0]
    
    lc = lcDict[str(sid)]
    lcTimes = Time(lc['time'].value + 2455197.5, format='jd')
    lcFlux = lc['flux']
    lcFerr = lc['flux_error']

    lcTS = TimeSeries(time = lcTimes, data={'flux' : lcFlux/np.median(lcFlux), 'flux_error' : lcFerr/np.median(lcFlux)})
    
    freqs = np.linspace(min_freq/u.d, max_freq/u.d, 20000)
    power = LombScargle(lcTS.time, lcTS['flux'], lcTS['flux_error']).power(freqs)
    best_freq = freqs[np.argmax(power)]


    right_mask = lcTS.time.value > xtime
    left_mask = lcTS.time.value <= xtime

    lcTS_right = TimeSeries(time = lcTimes[right_mask], data={'flux' : lcFlux[right_mask]/np.median(lcFlux[right_mask]), 'flux_error' : lcFerr[right_mask]/np.median(lcFlux[right_mask])})
    lcTS_left = TimeSeries(time = lcTimes[left_mask], data={'flux' : lcFlux[left_mask]/np.median(lcFlux[left_mask]), 'flux_error' : lcFerr[left_mask]/np.median(lcFlux[left_mask])})

    power_right = LombScargle(lcTS_right.time, lcTS_right['flux'], lcTS_right['flux_error']).power(freqs)
    power_left = LombScargle(lcTS_left.time, lcTS_left['flux'], lcTS_left['flux_error']).power(freqs)

    best_freq_right = freqs[np.argmax(power_right)]
    best_freq_left = freqs[np.argmax(power_left)]


    # Create subplots grid
    grid = plt.GridSpec(3, 4, wspace=0.2, hspace=0.45)
    fig = plt.figure(figsize=(6, 4), dpi=250)

    full_LSP = fig.add_subplot(grid[0, 0:])
    left_LSP = fig.add_subplot(grid[1, :2])
    right_LSP = fig.add_subplot(grid[1, 2:], sharey=left_LSP)
    lc_plot = fig.add_subplot(grid[2, 0:])

    ## Plot Periodograms
    full_LSP.plot(freqs, power, lw=0.5, color='g', label='Full Periodogram')
    left_LSP.plot(freqs, power_left, color='g', lw=0.5, label='Left Periodogram')
    right_LSP.plot(freqs, power_right, color='g', lw=0.5, label='Right Periodogram')

    # Add a dashed line at the peak frequency
    full_LSP.vlines(best_freq.value, ymin=power[np.where(freqs==best_freq)[0][0]], ymax=1.2, lw=0.5, ls='--', label=f'Peak Freq: {best_freq.value:0.3f}'+' d$^{-1}$')
    left_LSP.vlines(best_freq_left.value, ymin=power_left[np.where(freqs==best_freq_left)[0][0]], ymax=1.2, lw=0.5, ls='--', label=f'Peak Freq: {best_freq_left.value:0.3f}'+' d$^{-1}$')
    right_LSP.vlines(best_freq_right.value, ymin=power_right[np.where(freqs==best_freq_right)[0][0]], ymax=1.2, lw=0.5, ls='--', label=f'Peak Freq: {best_freq_right.value:0.3f}'+' d$^{-1}$')

    # Legend
    full_LSP.legend(prop={'size': 4}, loc=2)
    left_LSP.legend(prop={'size': 4}, loc=2)
    right_LSP.legend(prop={'size': 4}, loc=2)

    
    ## Light Curve Plot

    # Distance (and crossing time) error
    derror = (((star_row['dist84'] - star_row['dist']) + (star_row['dist'] - star_row['dist16']))/2).to('lyr').value*365.25

    # Plot crossing time and errors
    lc_plot.axvline(xtime, 0, 1, color='green', lw=1)
    lc_plot.axvline(xtime+derror, 0, 1, color='green', ls='dashed', lw=1)
    lc_plot.axvline(xtime-derror, 0, 1, color='green', ls='dashed', lw=1)
            
    # Plot light curve
    lc_plot.errorbar(lcTS.time.value, lcTS['flux'], yerr=lcTS['flux_error'], fmt='.', c='green', label='G', ms=3, elinewidth=1)
    
    # Get medians before and after crossing time
    leftMed, rightMed, leftMAD, rightMAD = getMedLC(str(sid), star, lcDict)/np.median(lcFlux)

    # Plot medians
    lc_plot.hlines(leftMed, xmin=min(lcTS.time.value), xmax=xtime, color='green', linewidth=0.5)
    lc_plot.hlines([leftMed+3*leftMAD, leftMed-3*leftMAD], xmin=min(lcTS.time.value), xmax=xtime, color='green', linewidth=0.5, ls = 'dashed')
    lc_plot.hlines(rightMed, xmin=xtime, xmax=max(lcTS.time.value), color='green', linewidth=0.5)
    lc_plot.hlines([rightMed+3*rightMAD, rightMed-3*rightMAD], xmin=xtime, xmax=max(lcTS.time.value), color='green', linewidth=0.5, ls = 'dashed')


    # Set axis labels
    full_LSP.set_xlabel('Frequency', fontsize=6)
    full_LSP.set_ylabel('Power', fontsize=6)

    left_LSP.set_xlabel('Frequency', fontsize=6)
    left_LSP.set_ylabel('Power', fontsize=6)

    right_LSP.set_xlabel('Frequency', fontsize=6)
    right_LSP.set_ylabel('Power', fontsize=6)

    lc_plot.set_xlabel('BJD', fontsize=6)
    lc_plot.set_ylabel('Normalized Flux', fontsize=6)

    full_LSP.xaxis.set_label_coords(0.5, -0.2)
    left_LSP.xaxis.set_label_coords(0.5, -0.2)
    right_LSP.xaxis.set_label_coords(0.5, -0.2)

    full_LSP.tick_params(axis='both', which='major', labelsize=5, pad=2)
    left_LSP.tick_params(axis='both', which='major', labelsize=5, pad=2)
    right_LSP.tick_params(axis='x', which='major', labelsize=5, pad=2)
    lc_plot.tick_params(axis='both', which='major', labelsize=5)
    
    right_LSP.get_yaxis().set_visible(False)

    lc_plot.xaxis.get_offset_text().set_fontsize(5)

    full_LSP.set_ylim((0,1.24))
    left_LSP.set_ylim((0,1.24))

    # Set titles
    full_LSP.set_title(f'Source {sid} Lomb-Scargle Periodograms and Light Curve', fontsize=8)
    lc_plot.set_title('G Band Light Curve', fontsize=6, y=0.93)

    if save:
        plt.savefig(f'{savefolder}/{sid}', transparent=False, facecolor='white')
        plt.close()

    
