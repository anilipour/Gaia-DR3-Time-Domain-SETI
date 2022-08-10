import numpy as np

from astropy.time import Time
from astropy.timeseries import LombScargle, TimeSeries
from astropy.table import vstack
import scipy.stats
from astropy import units as u

import ellipsoid 

from astroquery.gaia import Gaia
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit


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


def chunks(lst, n):
    ""
    "Split an input list into multiple chunks of size =< n"
    ""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def retrieveLCs(source_ids):
    # Retrieves light curve of a list of source ids
    dl_threshold = 5000               # DataLink server threshold. It is not possible to download products for more than 5000 sources in one single call.
    ids_chunks   = list(chunks(source_ids, dl_threshold))
    datalink_all = []


    retrieval_type = 'EPOCH_PHOTOMETRY'          # Options are: 'EPOCH_PHOTOMETRY', 'MCMC_GSPPHOT', 'MCMC_MSC', 'XP_SAMPLED', 'XP_CONTINUOUS', 'RVS', 'ALL'
    data_structure = 'COMBINED'     # Options are: 'INDIVIDUAL', 'COMBINED', 'RAW'
    data_release   = 'Gaia DR3'     # Options are: 'Gaia DR3' (default), 'Gaia DR2'
    dl_key         = f'{retrieval_type}_{data_structure}.xml'


    ii = 0
    for chunk in ids_chunks:
        ii = ii + 1
        print(f'Downloading Chunk #{ii}; N_files = {len(chunk)}')
        datalink  = Gaia.load_data(ids=chunk, data_release = data_release, retrieval_type=retrieval_type, format = 'votable', data_structure = data_structure)
        datalink_all.append(datalink)

    temp = [inp[dl_key][0].to_table() for inp in datalink_all]
    mergedLC = vstack(temp)   

    return mergedLC
    
def getLC(lc_table, source_id, band):
    # Returns light curve for a specific band and source id from a light curve
    # table of all sources in all bands
    source_mask = lc_table['source_id'] == source_id
    source_lc = lc_table[source_mask]
    band_mask = source_lc['band'] == band
    band_lc = source_lc[band_mask]

    return band_lc

def lcBandDict(lc_table):
    glcDict = {}
    bplcDict = {}
    rplcDict = {}
    source_ids = list(dict.fromkeys(lc_table['source_id']))
    for source_id in source_ids:
        glcDict[str(source_id)] = getLC(lc_table, source_id, 'G')
        bplcDict[str(source_id)] = getLC(lc_table, source_id, 'BP')
        rplcDict[str(source_id)] = getLC(lc_table, source_id, 'RP')

    return glcDict, bplcDict, rplcDict

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


def comparePeriodFoldPlot(lcDict, sid, star_row, star, min_freq=0.05, max_freq=10, save=False, savefolder=None, ecl=True, print_param=True):
    xtime = ellipsoid.xTime(star)[0]
    
    lc = lcDict[str(sid)]
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

    # Get folded light curves
    if ecl:
        # set epoch_time to the same for all folds to compare phase
        epoch_time = lcTS_left['time'][np.where(lcTS_left['flux']==np.max(lcTS_left['flux']))[0][0]]
        pLC, lcTS_folded = fitDoubleGaussian(lcTS, best_freq, epoch_time)
        pRight, lcTS_folded_right = fitDoubleGaussian(lcTS_right, best_freq_right, epoch_time)
        pLeft, lcTS_folded_left = fitDoubleGaussian(lcTS_left, best_freq_left, epoch_time)
    else:
        lcTS_folded = lcTS.fold(period=1./best_freq, normalize_phase=True, epoch_phase=0)
        lcTS_folded_right = lcTS_right.fold(period=1./best_freq_right, normalize_phase=True, epoch_phase=0)
        lcTS_folded_left = lcTS_left.fold(period=1./best_freq_left, normalize_phase=True, epoch_phase=0)


    # Create subplots grid
    grid = plt.GridSpec(4, 6, wspace=0.2, hspace=0.45)
    fig = plt.figure(figsize=(6, 6.67), dpi=250)

    full_LSP = fig.add_subplot(grid[0, 0:])
    left_LSP = fig.add_subplot(grid[1, :3])
    right_LSP = fig.add_subplot(grid[1, 3:], sharey=left_LSP)
    lc_plot = fig.add_subplot(grid[2, 0:])
    folded_lcPlot = fig.add_subplot(grid[3, 0:2])
    left_FLP = fig.add_subplot(grid[3, 2:4], sharey=folded_lcPlot)
    right_FLP = fig.add_subplot(grid[3, 4:6], sharey=folded_lcPlot)

    ## Plot Periodograms
    full_LSP.plot(freqs, power, lw=0.5, color='g', label='Full Periodogram')
    left_LSP.plot(freqs, power_left, color='g', lw=0.5, label='Left Periodogram')
    right_LSP.plot(freqs, power_right, color='g', lw=0.5, label='Right Periodogram')

    # Add a dashed line at the peak frequency
    full_LSP.vlines(best_freq.value, ymin=power[np.where(freqs==best_freq)[0][0]], ymax=1.2, lw=0.5, ls='--', label=f'Peak Freq: {best_freq.value:0.3f}'+' d$^{-1}$')
    left_LSP.vlines(best_freq_left.value, ymin=power_left[np.where(freqs==best_freq_left)[0][0]], ymax=1.2, lw=0.5, ls='--', label=f'Peak Freq: {best_freq_left.value:0.3f}'+' d$^{-1}$')
    right_LSP.vlines(best_freq_right.value, ymin=power_right[np.where(freqs==best_freq_right)[0][0]], ymax=1.2, lw=0.5, ls='--', label=f'Peak Freq: {best_freq_right.value:0.3f}'+' d$^{-1}$')

    # Legend
    full_LSP.legend(prop={'size': 4}, loc=0)
    left_LSP.legend(prop={'size': 4}, loc=0)
    right_LSP.legend(prop={'size': 4}, loc=0)

    
    ## Light Curve Plot

    # Distance (and crossing time) error
    derror = (((star_row['dist84'] - star_row['dist']) + (star_row['dist'] - star_row['dist16']))/2).to('lyr').value*365.25

    # Plot crossing time and errors
    lc_plot.axvline(xtime, 0, 1, color='green', lw=1)
    lc_plot.axvline(xtime+derror, 0, 1, color='green', ls='dashed', lw=1)
    lc_plot.axvline(xtime-derror, 0, 1, color='green', ls='dashed', lw=1)
            
    # Plot light curve
    lc_plot.errorbar(lcTS.time.value, lcTS['flux'], yerr=lcTS['flux_error'], fmt='.', c='green', label='G', ms=1.5, elinewidth=1)
    
    # Plot folded light curves
    folded_lcPlot.errorbar(lcTS_folded.time.value, lcTS_folded['flux'], yerr=lcTS_folded['flux_error'], fmt='.', c='green', label='Full Folded', ms=2, elinewidth=1)
    left_FLP.errorbar(lcTS_folded_left.time.value, lcTS_folded_left['flux'], yerr=lcTS_folded_left['flux_error'], fmt='.', c='green', label='Left Folded', ms=2, elinewidth=1)
    right_FLP.errorbar(lcTS_folded_right.time.value, lcTS_folded_right['flux'], yerr=lcTS_folded_right['flux_error'], fmt='.', c='green', label='Right Folded', ms=2, elinewidth=1)

    # If eclipsing binary, plot double Gaussian fit
    if ecl:
        xs = np.linspace(-0.5, 0.5, 100)
        folded_lcPlot.plot(xs, negDoubleGaussian(xs, pLC[0], pLC[1], pLC[2], pLC[3], pLC[4], pLC[5], pLC[6]), c='g')
        left_FLP.plot(xs, negDoubleGaussian(xs, pLeft[0], pLeft[1], pLeft[2], pLeft[3], pLeft[4], pLeft[5], pLeft[6]), c='g')
        right_FLP.plot(xs, negDoubleGaussian(xs, pRight[0], pRight[1], pRight[2], pRight[3], pRight[4], pRight[5], pRight[6]), c='g')



    # Get medians before and after crossing time
    leftMed, rightMed, leftMAD, rightMAD = getMedLC(str(sid), star, lcDict)/np.median(lcFlux)

    # Plot medians
    lc_plot.hlines(leftMed, xmin=min(lcTS.time.value), xmax=xtime, color='green', linewidth=0.5)
    lc_plot.hlines([leftMed+3*leftMAD, leftMed-3*leftMAD], xmin=min(lcTS.time.value), xmax=xtime, color='green', linewidth=0.5, ls = 'dashed')
    lc_plot.hlines(rightMed, xmin=xtime, xmax=max(lcTS.time.value), color='green', linewidth=0.5)
    lc_plot.hlines([rightMed+3*rightMAD, rightMed-3*rightMAD], xmin=xtime, xmax=max(lcTS.time.value), color='green', linewidth=0.5, ls = 'dashed')

    # Legends
    folded_lcPlot.legend(prop={'size': 4}, loc=0)
    left_FLP.legend(prop={'size': 4}, loc=0)
    right_FLP.legend(prop={'size': 4}, loc=0)

    ## Set labels

    # Axis labels
    full_LSP.set_xlabel('Frequency', fontsize=6)
    full_LSP.set_ylabel('Power', fontsize=6)

    left_LSP.set_xlabel('Frequency', fontsize=6)
    left_LSP.set_ylabel('Power', fontsize=6)

    right_LSP.set_xlabel('Frequency', fontsize=6)

    lc_plot.set_xlabel('BJD', fontsize=6)
    lc_plot.set_ylabel('Normalized Flux', fontsize=6)

    folded_lcPlot.set_xlabel('Phase', fontsize=6)
    folded_lcPlot.set_ylabel('Normalized Flux', fontsize=6)

    left_FLP.set_xlabel('Phase', fontsize=6)

    right_FLP.set_xlabel('Phase', fontsize=6)

    full_LSP.xaxis.set_label_coords(0.5, -0.2)
    left_LSP.xaxis.set_label_coords(0.5, -0.2)
    right_LSP.xaxis.set_label_coords(0.5, -0.2)

    full_LSP.tick_params(axis='both', which='major', labelsize=5, pad=2)
    left_LSP.tick_params(axis='both', which='major', labelsize=5, pad=2)
    right_LSP.tick_params(axis='x', which='major', labelsize=5, pad=2)
    lc_plot.tick_params(axis='both', which='major', labelsize=5)
    folded_lcPlot.tick_params(axis='both', which='major', labelsize=5)
    left_FLP.tick_params(axis='both', which='major', labelsize=5)
    right_FLP.tick_params(axis='both', which='major', labelsize=5)

    right_LSP.get_yaxis().set_visible(False)
    right_FLP.get_yaxis().set_visible(False)
    left_FLP.get_yaxis().set_visible(False)

    lc_plot.xaxis.get_offset_text().set_fontsize(5)

    full_LSP.set_ylim((0,1.24))
    left_LSP.set_ylim((0,1.24))

    # Set titles
    full_LSP.set_title(f'Source {sid} Lomb-Scargle Periodograms and Light Curve', fontsize=8)
    lc_plot.set_title('G Band Light Curve', fontsize=6, y=0.95)

    if save:
        plt.savefig(f'{savefolder}/{sid}', transparent=False, facecolor='white')
        plt.close()

    if print_param and ecl:
        if pLeft[0] < pLeft[3]:
            ld1, ld2 = pLeft[2], pLeft[5]
            lp1, lp2 = pLeft[0], pLeft[3]
        else:
            ld1, ld2 = pLeft[5], pLeft[2]
            lp1, lp2 = pLeft[3], pLeft[0]

        if pRight[0] < pRight[3]:
            rd1, rd2 = pRight[2], pRight[5]
            rp1, rp2 = pRight[0], pRight[3]
        else:
            rd1, rd2 = pRight[5], pRight[2]
            rp1, rp2 = pRight[3], pRight[0]
        
        print(f'Left Median: {leftMed:0.3f}, Right Median: {rightMed:0.3f}')
        print(f'Left Freq: {best_freq_left.value:0.3f}, Right Freq: {best_freq_left.value:0.3f}')
        print(f'Left Depth 1: {ld1:0.3f}, Right Depth 1: {rd1:0.3f}')
        print(f'Left Depth 2: {ld2:0.3f}, Right Depth 2: {rd2:0.3f}')
        print(f'Left Phase 1: {lp1:0.3f}, Right Phase 1: {rp1:0.3f}')
        print(f'Left Phase 2: {lp2:0.3f}, Right Phase 2: {rp2:0.3f}')


def negDoubleGaussian(x, mu1, sig1, d1, mu2, sig2, d2, b):
    g1 = d1*np.exp(-np.power(x - mu1, 2.) / (2 * np.power(sig1, 2.)))
    g2 = d2*np.exp(-np.power(x - mu2, 2.) / (2 * np.power(sig2, 2.)))
    ndg = -g1 - g2 + b
    return ndg

def fitDoubleGaussian(lc, frequency, epoch_time, p0=None):
    folded_LC = lc.fold(period=2./frequency, normalize_phase=True, epoch_time=epoch_time)
    mean1, mean2 = folded_LC['time'][np.where(folded_LC['flux']==np.min(folded_LC['flux']))[0][0]].value, None
    while mean2 == None:
        checkLow = folded_LC['time'][np.where(folded_LC['flux']==np.min(folded_LC['flux']))[0][0]].value
        if np.abs(checkLow - mean1) < 0.3:
            folded_LC.remove_row(np.where(folded_LC['flux']==np.min(folded_LC['flux']))[0][0])
        else:
            mean2 = checkLow
    
    folded_LC = lc.fold(period=2./frequency, normalize_phase=True, epoch_time=epoch_time)
    depth1, depth2 = folded_LC['flux'][folded_LC['time'] == mean1].value.data[0], folded_LC['flux'][folded_LC['time'] == mean2].value.data[0]
    
    if p0==None:
        p0 = np.array([mean1, 0.1, depth1, mean2, 0.1, depth2, np.max(folded_LC['flux'].value)])
    # print(p0, folded_LC, lc)
    p1, cov = curve_fit(negDoubleGaussian, folded_LC['time'].value, folded_LC['flux'].value, p0, maxfev=10000)
    return p1, folded_LC