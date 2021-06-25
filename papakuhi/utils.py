import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import astroquery
from astroquery.utils.tap.core import TapPlus as tap

def obs_cmd_data():
    '''
        Query Gaia for an Observed CMD, return a Matplotlib Figure.

    Outputs
    bp_rp: Gaia DR2 Bp_Rp values for observed CMD
    g_abs: Absolute Gaia G-band magnitude for observed CMD
    
    '''
    gaia = tap(url="https://gea.esac.esa.int/tap-server/tap")

    query = """SELECT bp_rp_index/40 AS bp_rp, g_abs_index/10 AS g_abs
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
                                ) AS subquery"""
    job = gaia.launch_job_async(query)
    obj = job.get_results()
    bprp = obj['bp_rp']
    g_abs =  obj['g_abs']
    
    return bp_rp, g_abs
