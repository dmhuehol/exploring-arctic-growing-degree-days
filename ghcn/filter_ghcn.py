'''filter_ghcn
Filter Global Historical Climatology Network (GHCN) stations based on 
some characteristics to obtain useful subsets the full dataset.
See the comments within filt_d for input format.
Filter by:
    Region: Implemented. See gddt_region_library for available regions.
    Variable: Implemented. Core elements are: PRCP (precip, tenths mm),
        SNOW (snowfall, mm), SNWD (snow depth, mm), TMAX (max temp, 
            tenths deg C), and TMIN (min temp, tenths deg C)
        www1.ncdc.noaa.gov/pub/data/ghcn/daily/readme.txt for full list
    Period of record: Implemented.

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''
import sys
sys.path.append(
    "/Users/dhueholt/Documents/Github/" 
    + "growing-degree-days-treelines-poles/")

from icecream import ic
import numpy as np

import fun_ghcn as f_g
import gddt_region_library as g_rlib

inv_path = "/Users/dhueholt/Documents/ghcn_data/"
inv_name = "ghcnd-inventory.txt"
filt_d = {
    #  reg_oi: 'global' or object from gddt_region_library
    "reg_oi": g_rlib.Arctic50N(),
    #  tol: latlon distance tolerance for point, 'closest' for closest
    "tol": 0.25,
    #  var_oi: single code as string or list of codes
    "var_oi": ['TMIN', 'TMAX'],
    #  span_type: '>' or '<' for greater than or less than span_por
    "span_type": '>',
    #  por: period of record in years
    "por": 1,
}
save_st = True
out_path = "/Users/dhueholt/Documents/ghcn_data/station_lists/"
out_name = 'auto' # ghcn_reg_var_span_pctvalid.dat

df_inv = f_g.open_inv(inv_path, inv_name)
df_inv_filt_reg = f_g.manage_area_ghcn(
    df_inv, filt_d['reg_oi'], tol=filt_d['tol'])
df_inv_filt_var = f_g.manage_var(df_inv_filt_reg, filt_d['var_oi'])
df_inv_filt_span = f_g.manage_span(
    df_inv_filt_var, filt_d['span_type'], filt_d['por'])

if save_st:
    if out_name == 'auto':
        reg_str = filt_d['reg_oi']['reg_abv']
        if np.size(filt_d['var_oi']) > 1:
            var_str = "-".join(filt_d['var_oi'])
        else:
            var_str = filt_d['var_oi']
        por_str = filt_d['span_type'] + str(filt_d['por'])
        out_name = 'ghcn_' + reg_str + '_' + var_str + '_' + por_str + '.csv'
    ic(out_name)
    np.savetxt(
        out_path + out_name, df_inv_filt_span, delimiter=",", fmt="%s")