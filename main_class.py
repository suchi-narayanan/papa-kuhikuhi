from astroquery.mast import Catalogs, Tesscut
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astroquery.utils.tap.core import TapPlus as tap
import os, sys, argparse
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')


gaia = tap(url="https://gea.esac.esa.int/tap-server/tap")

class papa:
    def __init__(self, radius, path=None, ralist=None, declist=None):
        '''
        
        '''
    

        self.path = path
        self.radius = radius
        self.ra = None
        self.dec = None
        self.merge_cols = ['TIC', 'Lit Name', 'DR2_source_id', 'eDR3_source_id', 'ASAS-SN Name',
       'Tmag', 'Num_sectors', 'eDR3_ra', 'eDR3_dec', 'DR2_ra', 'DR2_dec', 'parallax', 'phot_g_mean_mag',
       'bp_rp', 'ASAS-SN_Period', 'ASAS-SN_Var'] # Merged column ordering
        
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
            ra, dec = ra, dec 

        if type(ra) == np.float64 or type(ra) == float: # list compatibility
            ra, dec = [ra], [dec] 

        out_id =[]
        out_Tmag = []
        out_ra =[]
        out_dec =[]
        out_sector =[]
        out_gaia = []


        for i in range(len(ra)):
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
                out_id.append(ticid)
                out_Tmag.append(catalogData[0]['Tmag'])
                out_ra.append(Ra)
                out_dec.append(Dec)
                out_sector.append(len(sector_table))
                out_gaia.append(catalogData[0]['GAIA'])
                 
        return out_id, out_Tmag, out_ra, out_dec, out_sector, out_gaia


    def Gaia_query(self, ra=None, dec=None):
        "gaia EDR3"      
        gaia_dat = []

        if (ra is None) and (dec is None):
            ra, dec = self.ra, self.dec
        else:
            ra, dec = ra, dec
            
        if type(ra) == np.float64 or type(ra) == float: # list compatibility
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
            gaia_dat.append(r)

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
            
        if type(ra) == np.float64 or type(ra) == float: # list compatibility
            ra, dec = [ra], [dec]  
            
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
        if type(self.ra) == float: # single star
            ra, dec = [self.ra], [self.dec]
        else:
            ra, dec = self.ra, self.dec
        
        for i in range(len(ra)):
            tss = self.Tess_query(ra[i], dec[i])
            tss_df = pd.DataFrame({'TIC': tss[0],
                       'Tmag': tss[1],
                      'DR2_ra': tss[2],
                      'DR2_dec': tss[3],
                      'Num_sectors': tss[4],
                      'DR2_source_id': tss[5]})
            
            
            g = self.Gaia_query(ra[i], dec[i])
            g = g[0][['source_id', 'ra', 'dec', 'phot_g_mean_mag', 'bp_rp', 'parallax',
                 ]].to_pandas()
            g = g.rename(columns={'source_id': 'eDR3_source_id',
                     'ra': 'eDR3_ra',
                     'dec': 'eDR3_dec',
                     })
            
            asas = self.ASASN_query(ra[i], dec[i])[0][['ASAS-SN Name', 'Other Names', 'Period', 'Type']]
            asas = asas.rename(columns={'Other Names': 'Lit Name',
                   'Period': 'ASAS-SN_Period',
                    'Type': 'ASAS-SN_Var'})
            
            merged = pd.concat([g, asas, tss_df], axis=1)[self.merge_cols]
            
            if i > 0:
                stack = pd.concat([stack, merged])
            else:
                stack = merged
            
        return stack

    
def test():
    """
    unit test for ff.py
    """

    dir='test.dat'
    radius = 0.0005
    papa_test = papa(radius, path=dir)
    
    print(papa_test.ra[0], papa_test.dec[0])
    re = papa_test.combined_query()
    print(re)

if __name__ == '__main__':
    test()

        





      
