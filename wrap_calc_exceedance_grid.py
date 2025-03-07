''' wrap_calc_exceedance_grid
Script to calculate exceedance for gridded data and save as a netCDF 
file. These files can then be used to make maps of crossover and
no-analog states by applying different thresholds. To calculate 
a gridded file of crossover year from exceedance data, see
wrap_calc_crossover_grid.
'''
import sys

from icecream import ic
import numpy as np
import xarray as xr

import classes_gddt as cg
import fun_process as fproc
import gddt_region_library as g_rlib

dp_gdd = cg.DataParams(
    path='/Users/dhueholt/Documents/gddt_data/gdd/lens2/', 
    tok='*arc*.nc', var='gdd5_sum', flag_raw_ds=True, 
    flag_raw_da=True, flag_time_slice=True, flag_manage_rlz=True, 
    flag_land_mask=False, flag_roi=False)
dp_gdd_base = cg.DataParams(
    path='/Users/dhueholt/Documents/gddt_data/gdd/preindustrial/', 
    tok='*.nc', var='gdd5_sum', flag_raw_ds=False, flag_raw_da=True, 
    flag_time_slice=True, flag_manage_rlz=True, 
    flag_land_mask=True, flag_roi=True)
#  setp_gdd.base_yrs determines period of samples to compare to
setp_gdd = cg.SetParams(
    area_stat='mean', base_yrs=[0, 2000], mask_flag='none', 
    reg_oi=g_rlib.BrooksRange_colonist(), rlz='all', window=None, 
    yrs=[1850, 2100])
setp_gdd_base = cg.SetParams(
    area_stat='mean', base_yrs=[0, 2000], mask_flag='none', 
    reg_oi=setp_gdd.reg_oi, rlz='all', window=None, yrs=[0, 2000])
out_flag = True
base_str = str(setp_gdd.base_yrs[0]) + '-' + str(setp_gdd.base_yrs[1])
period_str = str(setp_gdd.yrs[0]) + '-' + str(setp_gdd.yrs[1])
out_path = '/Users/dhueholt/Documents/gddt_data/LENS2/exceedance/'
out_name = 'gexc_' + period_str + '_' + 'base' + base_str + '.nc'

dict_data = fproc.common_opener(dp=dp_gdd, setp=setp_gdd)
da_rlz = dict_data['manage_rlz'].compute()
dict_base = fproc.common_opener(dp=dp_gdd_base, setp=setp_gdd_base)
da_base_all = dict_base['raw_da'].compute()
da_base_period = da_base_all.sel(
    year=slice(setp_gdd.base_yrs[0], setp_gdd.base_yrs[1]))
np_base_samples_ravel = np.ravel(da_base_period.data)
try:
    np_base_samples = np_base_samples_ravel.reshape(
        (len(da_base_period.realization) * len(da_base_period.year), 
        len(da_base_period.lat), len(da_base_period.lon)))
except AttributeError:
    np_base_samples = np_base_samples_ravel.reshape(
        (1 * len(da_base_period.year), 
        len(da_base_period.lat), len(da_base_period.lon)))
years = da_rlz.year.data
list_exceed = list()
for year_count, year in enumerate(years):
    #  Helps track progress
    ic(year, np.shape(list_exceed))
    da_rlz_year = da_rlz.sel(year=year)
    np_compare_to_base = da_rlz_year.data
    dict_g = fproc.exceed_subceed(np_base_samples, np_compare_to_base)
    np_gexc = np.array(dict_g['gexc'])
    list_exceed.append(np_gexc)
np_exceed_yearlatlon = np.stack(list_exceed, axis=0)
#  For data with realization dimension
try:
    ds_exceed = xr.Dataset(
        data_vars = {
            "exceedance": (
                ("year", "realization", "lat", "lon"), np_exceed_yearlatlon),       
        },
        coords = {
            "year": years,
            "realization": da_rlz.realization.data,
            "lon": da_rlz.lon.data,
            "lat": da_rlz.lat.data
        },
        attrs = {
            "long_name": 'member exceedance',
            "length_of_base_period": len(np_base_samples)
        }
    )
#  For data with no realization dimension
except AttributeError:
    out_name = 'ensmn-forced_' + out_name
    ds_exceed = xr.Dataset(
        data_vars = {
            "exceedance": (("year", "lat", "lon"), np_exceed_yearlatlon),       
        },
        coords = {
            "year": years,
            "lon": da_rlz.lon.data,
            "lat": da_rlz.lat.data
        },
        attrs = {
            "long_name": 'unknown',
            "length_of_base_period": len(np_base_samples)
        }
    )
ic(ds_exceed)
if out_flag:
    ic(out_path + out_name)
    ds_exceed.to_netcdf(out_path + out_name)