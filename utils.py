import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import gaia.tap

def query_asas_sn(ra='', dec='', radius=0.5):
    '''
    Inputs
    ra : Right Ascension of Target in degrees
    dec : Declination of Target in degrees
    radius : Size of cone search query in arcminutes (default 0.5)


    Outputs
    tab: pandas DataFrame containing ASAS-SN parameters
    
    TODO: Wrap this in object so we pass a Class instance
    ''' 
    
    try:
        ra, dec, radius = np.float(ra), np.float(dec), np.float(radius)
    except:
        raise ValueError('Input of RA=%s, Dec=%s, Radius=%s not allowed.' %(ra, dec, radius))
    
    baseurl = 'https://asas-sn.osu.edu/variables'
    exturl = '?ra=%s&dec=%s&radius=%s&vmag_min=&vmag_max=&amplitude_min=&amplitude_max=&period_min=&period_max=&lksl_min=&lksl_max=&class_prob_min=&class_prob_max=&parallax_over_err_min=&parallax_over_err_max=&name=&references[]=I&references[]=II&references[]=III&references[]=IV&references[]=V&references[]=VI&sort_by=raj2000&sort_order=asc&show_non_periodic=true&show_without_class=true&asassn_discov_only=false&'\
                %(str(ra), str(dec), str(radius))
    
    return pd.read_html(baseurl+exturl)[0]
    

def obs_cmd_data():
    '''
        Query Gaia for an Observed CMD, return a Matplotlib Figure.

    Outputs
    bp_rp: Gaia DR2 Bp_Rp values for observed CMD
    g_abs: Absolute Gaia G-band magnitude for observed CMD
    
    '''
    query = gaia.tap.query("""SELECT bp_rp_index / 40 AS bp_rp, g_abs_index / 10 AS g_abs, n
                                FROM (
                                SELECT 
                                floor(bp_rp * 40) AS bp_rp_index,
                                floor((phot_g_mean_mag + 5 * log10(parallax) - 10) * 10) AS g_abs_index,
                                count(*) AS n
                                FROM gaiadr2.gaia_source
                                WHERE parallax_over_error > 1
                                AND a_g_val IS NOT NULL
                                AND random_index < 850000
                                GROUP BY bp_rp_index, g_abs_index
                                ) AS subquery """)
    
    bp_rp, g_abs = query['bp_rp'].data, query['g_abs'].data
    
    return bp_rp, g_abs
