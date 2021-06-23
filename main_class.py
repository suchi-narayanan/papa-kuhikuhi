from astroquery.mast import Catalogs, Observations, Tesscut
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astroquery.utils.tap.core import TapPlus as tap
from astropy.stats import mad_std, LombScargle
from astroquery.gaia import Gaia
from utils import *

import os, sys, argparse, astroquery
import numpy as np
import pandas as pd
import astropy.units as u


gaia = tap(url="https://gea.esac.esa.int/tap-server/tap")

class papa:
    'parent class'
    def __init__(self, radius, path=None, ralist=None, declist=None):
        "path: file containing target coordinate"
        
        self.path = path
        self.radius = radius
        self.ra = None
        self.dec = None
        
        ## Input list file ##
        if self.path is not None:
            dat = pd.read_csv(self.path)
            c = SkyCoord(ra=dat['RA'], dec=dat['Dec'], unit=('hourangle', 'deg'), frame='icrs')
            self.ra = c.ra.degree
            self.dec = c.dec.degree
        
        if ralist is not None and declist is not None: #provide your own list of coordinates
            self.ra = ralist
            self.dec = declist
            
        if self.ra is None or self.dec is None:
            raise ValueError('Coordinates have not been loaded in!')
            
            
    def check_inputfile(self):
        "check whether input file is right/not null"
        "print the number of input targets"
        pass


    def Tess_query(self, ra=None, dec=None):

        if (ra is None) and (dec is None):
            ra, dec = self.ra, self.dec
        else:
            ra, dec = [ra], [dec] 

        out_id =[]
        out_Tmag = []
        out_ra =[]
        out_dec =[]
        out_sector =[]

        for i in range(len(ra)):
            print(self.radius, ra[i], dec[i])
            catalogData = Catalogs.query_object("%s %s"%(ra[i], dec[i]), radius = self.radius, catalog = "TIC")
            if not catalogData:
                print('No TIC counterpart')
            else:
                #print( catalogData[:nstar]['ID', 'Tmag', 'Jmag', 'ra', 'dec', 'objType'] )
                ticid = int(catalogData[0]['ID'])
                #print('TIC counterpart = ',ticid)
                Ra = catalogData[0]['ra']
                Dec = catalogData[0]['dec']
                coord = SkyCoord(Ra, Dec, unit="deg")
                sector_table = Tesscut.get_sectors(coord)
                
            if not sector_table:
                print('TESS has not observed this star')
            else:
                #print(sector_table)
                out_id.append(ticid)
                out_Tmag.append(catalogData[0]['Tmag'])
                out_ra.append(Ra)
                out_dec.append(Dec)
                out_sector.append(len(sector_table))
        return out_id, out_Tmag, out_ra, out_dec, out_sector


    def Gaia_query(self, ra=None, dec=None):
        "gaia EDR3"      
        gaia_dat = []

        if (ra is None) and (dec is None):
            ra, dec = self.ra, self.dec
        else:
            ra, dec = [ra], [dec] 

        for i in range(len(ra)):
            "check the unit of radius"
            query = f"SELECT gaia_source.source_id, gaia_source.ra, gaia_source.dec,\
            gaia_source.phot_g_n_obs, gaia_source.phot_g_mean_flux,\
            gaia_source.phot_g_mean_flux_error,\
            gaia_source.bp_rp, gaia_source.g_rp, gaia_source.phot_g_mean_mag,\
            gaia_source.parallax, gaia_source.pmra, gaia_source.pmdec \
            FROM gaiaedr3.gaia_source  WHERE CONTAINS( \
            POINT('ICRS',gaiaedr3.gaia_source.ra,gaiaedr3.gaia_source.dec), \
            CIRCLE('ICRS',{ra[i]},{dec[i]}, {self.radius}))=1\
            AND g_rp > -10\
            ORDER BY ra ASC"
            job = gaia.launch_job_async(query)
            r = job.get_results()

            if len(r) > 0:
                gaia_dat.append(r)
            else:
                gaia_dat.append(0)

        return gaia_dat   

    def ASASN_query(self, ra=None, dec=None):
        '''
        Inputs
        ra : Right Ascension of Target in degrees
        dec : Declination of Target in degrees
        radius : Size of cone search query in arcminutes (default 0.5)
        Outputs
        tab: pandas DataFrame containing ASAS-SN parameters
        
        TODO: Wrap this in object so we pass a Class instance
        '''

        if (ra is None) and (dec is None):
            ra, dec = self.ra, self.dec
        else:
            ra, dec = ra, dec

        out_asassn = []
        for i in range(len(ra)):
            try:
                ra, dec, radius = np.float(ra[i]), np.float(dec[i]), np.float(self.radius*60.)
            except:
                raise ValueError('Input of RA=%s, Dec=%s, Radius=%s not allowed.' %(ra, dec, radius))
            
            baseurl = 'https://asas-sn.osu.edu/variables'
            exturl = '?ra=%s&dec=%s&radius=%s&vmag_min=&vmag_max=&amplitude_min=&amplitude_max=&period_min=&period_max=&lksl_min=&lksl_max=&class_prob_min=&class_prob_max=&parallax_over_err_min=&parallax_over_err_max=&name=&references[]=I&references[]=II&references[]=III&references[]=IV&references[]=V&references[]=VI&sort_by=raj2000&sort_order=asc&show_non_periodic=true&show_without_class=true&asassn_discov_only=false&'\
                        %(str(ra), str(dec), str(radius))
            out_asassn.append(pd.read_html(baseurl+exturl)[0])           
            
        return out_asassn
    
    def combined_query(self):
        '''
        combined gaia-tess-ASASSN query

        '''
        for i in range(len(self.ra)):
            Tess_query(self.ra[i], self.dec[i])
            Gaia_query(self.ra[i], self.dec[i])
            ASASN_query(self.ra[i], self.dec[i])

    
def test():
    """
    unit test for ff.py
    """

    dir='/Users/zjw/Desktop/test/std_info.dat'
    radius = 0.00001
    papa_test = papa(dir, radius)
    print(papa_test.ra[0], papa_test.dec[0])
    re = papa_test.Gaia_query() 
    print(re)

if __name__ == '__main__':
    test()

        





      
