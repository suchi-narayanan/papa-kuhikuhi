from astroquery.mast import Catalogs, Tesscut
import os
import numpy as np
from astroquery.mast import Observations
from astropy.coordinates import SkyCoord
from astropy.io import fits
import sys
from astropy.stats import mad_std, LombScargle
import argparse
from astroquery.gaia import Gaia
import astroquery
from astroquery.utils.tap.core import TapPlus as tap
import pandas as pd
import astropy.units as u
gaia = tap(url="https://gea.esac.esa.int/tap-server/tap")


class papa:
    'parent class'
    def __init__(self, path, radius):
        "path: file containing target coordinate"
        
        self.path = path
        self.radius = radius
        dat = pd.read_csv(self.path)
        c = SkyCoord(ra=dat['RA'], dec=dat['Dec'], unit=('hourangle', 'deg'), frame='icrs')

        self.ra = c.ra.degree
        self.dec = c.dec.degree
        
        

    def check_inputfile(self):
        "check whether input file is right/not null"
        "print the number of input targets"
        pass
    def Tess_query(self):
        #super().__init__()
        #result = Catalogs.query_criteria(catalog="Tic", ID=self.ticid)
        out_id =[]
        out_Tmag = []
        out_ra =[]
        out_dec =[]
        out_sector =[]

        for i in range(len(self.ra)):
            catalogData = Catalogs.query_object(self.ra[i], self.dec[i], radius = self.radius, catalog = "TIC")
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
        return out.id, out_Tmag, out_ra, out_dec, out_sector

    def Gaia_query(self):
        "gaia EDR3"      
        gaia_dat = [] 

        for i in range(len(self.ra)):
            "check the unit of radius"
            query = f"SELECT gaia_source.source_id, gaia_source.ra, gaia_source.dec,\
            gaia_source.phot_g_n_obs, gaia_source.phot_g_mean_flux,\
            gaia_source.phot_g_mean_flux_error,\
            gaia_source.bp_rp, gaia_source.g_rp, gaia_source.phot_g_mean_mag,\
            gaia_source.parallax, gaia_source.pmra, gaia_source.pmdec \
            FROM gaiaedr3.gaia_source  WHERE CONTAINS( \
            POINT('ICRS',gaiaedr3.gaia_source.ra,gaiaedr3.gaia_source.dec), \
            CIRCLE('ICRS',{self.ra[i]},{self.dec[i]}, {self.radius}))=1\
            AND g_rp > -10\
            ORDER BY ra ASC"
            job = gaia.launch_job_async(query)
            r = job.get_results()
            if len(r) > 0:
                gaia_dat.append(r)
            else:
                gaia_dat.append(0)

        return gaia_dat   

    def ASASN_query(self):
        '''
        Inputs
        ra : Right Ascension of Target in degrees
        dec : Declination of Target in degrees
        radius : Size of cone search query in arcminutes (default 0.5)
        Outputs
        tab: pandas DataFrame containing ASAS-SN parameters
        
        TODO: Wrap this in object so we pass a Class instance
        '''
        out_asassn = []
        for i in range(len(self.ra)):
            try:
                ra, dec, radius = np.float(self.ra[i]), np.float(self.dec[i]), np.float(self.radius)
            except:
                raise ValueError('Input of RA=%s, Dec=%s, Radius=%s not allowed.' %(ra, dec, radius))
            
            baseurl = 'https://asas-sn.osu.edu/variables'
            exturl = '?ra=%s&dec=%s&radius=%s&vmag_min=&vmag_max=&amplitude_min=&amplitude_max=&period_min=&period_max=&lksl_min=&lksl_max=&class_prob_min=&class_prob_max=&parallax_over_err_min=&parallax_over_err_max=&name=&references[]=I&references[]=II&references[]=III&references[]=IV&references[]=V&references[]=VI&sort_by=raj2000&sort_order=asc&show_non_periodic=true&show_without_class=true&asassn_discov_only=false&'\
                        %(str(ra), str(dec), str(radius))
            out_asassn.append(pd.read_html(baseurl+exturl)[0])           
            
        return out_asassn
    
    
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

        





      
