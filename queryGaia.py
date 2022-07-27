import argparse
from astroquery.gaia import Gaia
import os
from astropy.coordinates import SkyCoord
from astropy.table import QTable
from astropy import units as u

def checkBool(arg):
    if arg == 1:
        return True
    else:
        return False

def queryGaia(var, gcns, num, save, c):
    # Retrieves data from Gaia archive
    # var = variable stars or all stars?
    # gcns = GCNS or all Gaia stars?
    # num = max number of stars to retrieve
    # save = save to a file? 
    # c = include variable classification?
    
    if gcns == True:
        dist_col = 'dist_50'
        cartesian = True
        x_col, y_col, z_col = 'xcoord_50', 'ycoord_50', 'zcoord_50'
        ra_col, dec_col = 'ra', 'dec'
        dist84_col, dist16_col = 'dist_84', 'dist_16'
        gmag_col, bpmag_col, rpmag_col = 'phot_g_mean_mag', 'phot_bp_mean_mag', 'phot_rp_mean_mag'
        
        query_select = f'SELECT TOP {num} gcns.{dist_col}, gcns.{x_col}, gcns.{y_col}, gcns.{z_col}, \
                        gcns.{ra_col}, gcns.{dec_col}, gcns.{dist84_col}, gcns.{dist16_col}, \
                        gcns.{gmag_col}, gcns.{bpmag_col}, gcns.{rpmag_col}, gcns.source_id'
        
        
        if var == True:
            
            if c:
                query_select += ', vclass.best_class_name'
                sf = 'GCNS_var_class.fits'
                query = f"{query_select} \
                        FROM gaiadr3.vari_summary AS var \
                        INNER JOIN external.gaiaedr3_gcns_main_1 AS gcns ON var.source_id=gcns.source_id \
                        INNER JOIN gaiadr3.vari_classifier_result AS vclass ON vclass.source_id = var.source_id"
                        
            else:
                query = f'{query_select} \
                        FROM gaiadr3.vari_summary AS var \
                        JOIN external.gaiaedr3_gcns_main_1 AS gcns ON var.source_id=gcns.source_id'
                sf = 'GCNS_var.fits'
            
        elif var == False:
            query = f'{query_select} \
            FROM external.gaiaedr3_gcns_main_1 AS gcns'
            sf = 'GCNS.fits'
            
    elif gcns == False:
        dist_col = 'r_med_photogeo'
        cartesian = False
        ra_col, dec_col = 'ra', 'dec'
        dist84_col, dist16_col = 'r_hi_photogeo', 'r_lo_photogeo'
        gmag_col, bpmag_col, rpmag_col = 'phot_g_mean_mag', 'phot_bp_mean_mag', 'phot_rp_mean_mag'
        
        query_select = f'SELECT TOP {num} dist.{dist_col}, \
                        source.{ra_col}, source.{dec_col}, dist.{dist84_col}, dist.{dist16_col}, \
                        source.{gmag_col}, source.{bpmag_col}, source.{rpmag_col}, source.source_id'
               

        if var == True:
            if c:
                query_select += ', vclass.best_class_name'
                sf = 'Gaia_var_class.fits'
                query = f"{query_select} \
                        FROM gaiadr3.vari_summary AS var \
                        INNER JOIN gaiadr3.gaia_source AS source ON var.source_id=source.source_id \
                        INNER JOIN external.gaiaedr3_distance AS dist ON var.source_id=dist.source_id \
                        INNER JOIN gaiadr3.vari_classifier_result AS vclass ON vclass.source_id = var.source_id \
                        WHERE source.has_epoch_photometry='true'"
                        
            else:
                query = f"{query_select} \
                        FROM gaiadr3.vari_summary AS var \
                        INNER JOIN gaiadr3.gaia_source AS source ON var.source_id=source.source_id \
                        INNER JOIN external.gaiaedr3_distance AS dist ON var.source_id=dist.source_id \
                        WHERE source.has_epoch_photometry='true'"
                sf = 'Gaia_var.fits'
        
        elif var == False:
            query = f"{query_select} \
            FROM gaiadr3.gaia_source AS source \
            JOIN external.gaiaedr3_distance AS dist ON source.source_id=dist.source_id"
            sf = 'Gaia.fits'
            

    job = Gaia.launch_job_async(query, output_format='fits')
    
    results = job.get_results()

    
    
    c1 = SkyCoord(ra = results[ra_col],
              dec = results[dec_col],
              distance = results[dist_col],
              frame = 'icrs')
    
    if cartesian == True:
        xcoord, ycoord, zcoord = results[x_col], results[y_col], results[z_col]
    else:
        xcoord = c1.transform_to('galactocentric').x.to('pc') + 8122*u.pc
        ycoord = c1.transform_to('galactocentric').y.to('pc')
        zcoord = c1.transform_to('galactocentric').z.to('pc') - 20.8*u.pc
    
    if c:
        stars = QTable([results['source_id'], results[ra_col], results[dec_col], results[dist_col], results[dist84_col], results[dist16_col], xcoord, ycoord, zcoord, results[gmag_col], results[bpmag_col], results[rpmag_col], results['best_class_name']],
                   names=('id', 'ra', 'dec', 'dist', 'dist84', 'dist16', 'x', 'y', 'z', 'g', 'bp', 'rp', 'class'))
    else:
        stars = QTable([results['source_id'], results[ra_col], results[dec_col], results[dist_col], results[dist84_col], results[dist16_col], xcoord, ycoord, zcoord, results[gmag_col], results[bpmag_col], results[rpmag_col]],
                   names=('id', 'ra', 'dec', 'dist', 'dist84', 'dist16', 'x', 'y', 'z', 'g', 'bp', 'rp'))


    savefile = f'../{sf}'
    if save:
        if os.path.exists(savefile):
            os.remove(savefile)    
        stars.write(savefile, format='fits')

        print(f'{len(stars)} stars retrieved and saved to {savefile}')

    else:
        print(f'{len(stars)} stars retrieved')
    return c1, stars

def login(username, password):
    Gaia.login(user=username, password=password)

def logout():
    Gaia.logout()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Select targets on SETI Ellipsoid'
        )
    parser.add_argument(
        '-v', '--variable', type = int, default=1,
        help='Only variable stars (1) or all stars (0)'
        )
    parser.add_argument(
        '-ns', '--GCNS', type = int, default=1,
        help='Use Gaia Catalogue of Nearby Stars (GCNS) (1) or entire Gaia catalogue (0)?'
        )
    parser.add_argument(
        '-n', '--number', type=int, default=100000,
        help='Maximum number of query targets'
        )
    parser.add_argument(
        '-usr', '--username', type=str, default=None,
        help='Gaia login username'
        )
    parser.add_argument(
        '-pw', '--password', type=str, default=None,
        help='Gaia login password'
        )
    parser.add_argument(
        '-c', '--classification', type=int, default=0,
        help='Include variable classification (1) or not (0)?'
        )
    parser.add_argument(
        '-s', '--save', type=int, default=1,
        help='Save results to a fits file (1) or not (0)?'
    )
  

    args = parser.parse_args()

    v, ns, c, s = checkBool(args.variable), checkBool(args.GCNS), checkBool(args.classification), checkBool(args.save)
    num = args.number
    username = args.username
    password = args.password

    if username and password:
        login(username, password)

    
    c1, stars = queryGaia(v, ns, num, s, c)

    logout()

    