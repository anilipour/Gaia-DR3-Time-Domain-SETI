from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import constants as const
from astropy import units as u

import numpy as np
from scipy.stats import binned_statistic_2d

import matplotlib.pyplot as plt
from astropy.table import hstack, Table




def ringsNow(event, t0, d0_err, stars):
    # event is a SkyCoord object of the reference event (with ra, dec, distance)
    # t0 is a Time object of the Earth observation time of the reference event
    # d0_err is the error in distance to the reference event
    # stars is a SkyCoord object of all the stars to be searched

    d0 = event.distance

    deltaTE = Time.now() - t0 # time difference since arrival of SN signal to Earth
    tauE = (const.c/d0.to('m'))*deltaTE.to('s') # Seto's 'normalized time'
    q = 1 # Seto's mode -- we use the primary mode q = 1
    
    
    thetaPlus = (np.arctan(q) + np.arccos((1+tauE.value)/np.sqrt(1+q**2)))*u.rad
    thetaMinus = (np.arctan(q) - np.arccos((1+tauE.value)/np.sqrt(1+q**2)))*u.rad
    # angles to search

    depthPlus = d0*np.cos(thetaPlus)
    depthMinus = d0*np.cos(thetaMinus)
    # effective search depths for the two angles

    dTheta = (tauE*d0_err/(np.abs(np.cos(thetaPlus) - np.sin(thetaPlus))*d0))*u.rad
    # finite width for rings based on distance error

    separations = stars.separation(event)
    # angular separation between stars and reference event

    maskPlus = np.abs(separations.degree - thetaPlus.to('deg').value) < dTheta.to('deg').value
    maskMinus = np.abs(separations.degree - thetaMinus.to('deg').value) < dTheta.to('deg').value
    ringPlus = stars[maskPlus]
    ringMinus = stars[maskMinus]
    # rings without distance consideration

    distMaskMinus = ringMinus.distance < depthMinus
    distMaskPlus = ringPlus.distance < depthPlus
    
    distRingPlus = ringPlus[distMaskPlus]
    distRingMinus = ringMinus[distMaskMinus]

    return distRingPlus, distRingMinus


def xtimes(event, stars, returnTable=False, star_table=None):
    d0 = event.distance.to('lyr')

    separations = stars.separation(event)

    xtime = d0.to('m')/const.c*(np.cos(separations)+np.sin(separations)-1)
    within = stars.distance.to('lyr').value < d0.value*np.cos(separations)


    # xtimes are the time since Earth observation the stars will cross the observing line
    if returnTable and star_table:
        return stars[within], xtime[within], star_table[within]
    else:
        return stars[within], xtime[within]

def plotXtimes(stars, xtime, t0, event, name, q=1):
    maxTauE = np.sqrt(1+q**2) - 1 # if tauE is greater than this, there are no angle solutions
    d0 = event.distance

    H, xe, ye, bn = binned_statistic_2d(stars.ra.wrap_at(180*u.degree).radian, stars.dec.radian,
                    values = xtime.to('yr').value + t0.decimalyear, 
                    statistic='mean', bins=175) # mean binning for contours
    
    fig = plt.figure(figsize=(8,5), dpi=150)
    ax = fig.add_subplot(111, projection="mollweide")


    ax.scatter(stars.ra.wrap_at(180*u.degree).radian, stars.dec.radian,
                c=(xtime.to('yr').value + t0.decimalyear), cmap=plt.cm.autumn, rasterized=True, s=1)


    maxYear = np.floor(((maxTauE*d0.to('m')/const.c).to('yr')+t0).decimalyear)-2
    contourLevels = np.linspace(Time.now().value.year, maxYear, 5).astype(int)
    cs = plt.contour(xe[1:], ye[1:], H.T, colors='k', linewidths=1, levels=contourLevels) # binned contours
    plt.clabel(cs, fontsize=8)



    ax.set_xticklabels(['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h'])
    ax.tick_params(axis='both', labelsize=6)
    ax.grid(True)
    ax.set_title(f'Gaia Stars - Seto Signalling Schematic for SN{name}', fontsize=10)
    plt.show()


def rings(event, t0, stars, start, end, returnTable=False, star_table=None):
    deltaStart, deltaEnd = start - t0, end - t0

    # get the crossing times of all the stars
    
    if returnTable and star_table:
        c1, xtime, starTable = xtimes(event, stars, returnTable=returnTable, star_table=star_table)
    else:
        c1, xtime = xtimes(event, stars, returnTable=returnTable, star_table=star_table)
    
    

    # get the stars with crossings times within the specified dates 
    inRangeC = (xtime.to('yr').value > deltaStart.to('yr').value) & (xtime.to('yr').value < deltaEnd.to('yr').value)


    if returnTable and star_table:
        starTable = hstack([starTable, Table({'xtime' : xtime+t0})])
        return c1[inRangeC], starTable[inRangeC]
    else:
        return c1[inRangeC]