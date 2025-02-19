'''filter_ghcn
Filter Global Historical Climatology Network (GHCN) stations to obtain a
subset of interest. 

Filter options are handled by a FilterParams instance. See documentation 
in classes_gddt for details on syntax. Implemented options are:
    Region: Select stations within a region or near a given point . See 
        gddt_region_library for available regions.
    Variable: Select stations that have certain kinds of weather 
        observations. Core elements are: PRCP (precip, tenths mm), SNOW 
        (snowfall, mm), SNWD (snow depth, mm), TMAX (max temp, tenths 
        deg C), and TMIN (min temp, tenths deg C). Full list is at:
        www1.ncdc.noaa.gov/pub/data/ghcn/daily/readme.txt
    Period of record: Select stations with a period of record greater or
        less than some number of years. WARNING: span information in the 
        station list is not always correct! It's best to be liberal with 
        this here (i.e., set to 1) and use check_spancov_ghcn later to 
        constrain period of record from the data entries themselves.
'''
import sys
sys.path.append(
    "/Users/dhueholt/Documents/Github/" 
    + "exploring-arctic-growing-degree-days/")

from icecream import ic
import numpy as np

import classes_gddt as cg
import fun_ghcn as f_g
import gddt_region_library as g_rlib

inv_path = "/Users/dhueholt/Documents/ghcn_data/"
inv_name = "ghcnd-inventory.txt"
fp = cg.FilterParams(
    por=1, reg_oi=g_rlib.Arctic50N(), span_type='>', tol=None, 
    var_oi=['TMIN', 'TMAX'])
save_st = True
out_path = "/Users/dhueholt/Documents/ghcn_data/station_lists/"
#  Default naming scheme: ghcn_reg_var_span_pctvalid.dat
out_name = 'auto'

df_inv = f_g.open_inv(inv_path, inv_name)
df_inv_filt_reg = f_g.manage_area_ghcn(df_inv, fp.reg_oi, tol=fp.tol)
df_inv_filt_var = f_g.manage_var(df_inv_filt_reg, fp.var_oi)
df_inv_filt_span = f_g.manage_span(
    df_inv_filt_var, fp.span_type, fp.por)
if save_st:
    if out_name == 'auto':
        reg_str = fp.reg_oi['reg_abv']
        if np.size(fp.var_oi) > 1:
            var_str = "-".join(fp.var_oi)
        else:
            var_str = fp.var_oi
        por_str = fp.span_type + str(fp.por)
        out_name = 'ghcn_' + reg_str + '_' + var_str + '_' + por_str + '.csv'
    ic(out_name)
    np.savetxt(
        out_path + out_name, df_inv_filt_span, delimiter=",", fmt="%s")