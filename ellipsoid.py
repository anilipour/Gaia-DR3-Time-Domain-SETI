import numpy as np


from astropy.coordinates import SkyCoord
from astropy import constants as const
from astropy import units as u
from astropy.time import Time
from astropy.table import Table

from astroquery.gaia import Gaia


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

    if c0==None or t0==None:
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

