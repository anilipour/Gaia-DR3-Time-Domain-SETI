from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import constants as const
from astropy import units as u

import numpy as np


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