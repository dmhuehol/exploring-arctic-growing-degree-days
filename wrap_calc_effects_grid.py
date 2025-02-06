'''wrap_calc_effects_grid
Script to calculate an effect size statistic for gridded data and save
as a netCDF file. Supports the following methods as inputs to setp_gdd.effect:
    'gexc': Samples exceeded by each realization (Gexc in robustness)

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University 
'''
import sys

from icecream import ic
import numpy as np
import xarray as xr

import classes_gddt as cg
import fun_process as fproc
import gddt_region_library as g_rlib

dp_gdd = cg.DataParams(
    path='/Users/dhueholt/Documents/gddt_data/gdd/nomask/', 
    tok='*arc*.nc', var='gdd5_sum', flag_raw_ds=False, flag_raw_da=True, 
    flag_time_slice=True, flag_manage_rlz=True, 
    flag_land_mask=False, flag_roi=False)
dp_gdd_base = cg.DataParams(
    path='/Users/dhueholt/Documents/gddt_data/gdd/preindustrial/', 
    tok='*.nc', var='gdd5_sum', flag_raw_ds=False, flag_raw_da=True, 
    flag_time_slice=True, flag_manage_rlz=True, 
    flag_land_mask=True, flag_roi=True)
setp_gdd = cg.SetParams(
    area_stat=None, base_yrs=[0, 2000],  mask_flag='none', effect='gexc', 
    reg_oi='global', rlz='all', window=None, yrs=[1850, 2100])
setp_gdd_base = cg.SetParams(
    area_stat='mean', base_yrs=[0, 2000], beat=4000, mask_flag='none', 
    effect='gexc',
    reg_oi=setp_gdd.reg_oi, rlz='all', 
    window=None, yrs=[0, 2000])
out_flag = True
base_str = str(setp_gdd.base_yrs[0]) + '-' + str(setp_gdd.base_yrs[1])
per_str = str(setp_gdd.yrs[0]) + '-' + str(setp_gdd.yrs[1])
out_path = '/Users/dhueholt/Documents/gddt_data/LENS2/effects/'
out_name = 'gexc_' + per_str + '_' + 'base' + base_str + '.nc'

d_data = fproc.common_opener(dp=dp_gdd, setp=setp_gdd)
da_all = d_data['raw_da'].compute()
d_base = fproc.common_opener(dp=dp_gdd_base, setp=setp_gdd_base)
da_base = d_base['raw_da'].compute()
da_base_period = da_base.sel(year=slice(setp_gdd.base_yrs[0], setp_gdd.base_yrs[1]))
np_bp = np.ravel(da_base_period.data)
try:
    np_x1 = np_bp.reshape(
        (len(da_base_period.realization) * len(da_base_period.year), 
        len(da_base_period.lat), len(da_base_period.lon)))
except AttributeError:
    np_x1 = np_bp.reshape(
        (1 * len(da_base_period.year), 
        len(da_base_period.lat), len(da_base_period.lon)))
da_rlz = d_data['time_slice'].compute()
years = da_rlz.year.data
l_gexc_yrs = list()
for yrc, yr in enumerate(years):
    ic(yr)
    da_rlz_yr = da_rlz.sel(year=yr)
    np_x2 = da_rlz_yr.data
    d_gexc = fproc.gexc(np_x1, np_x2)
    np_gexc = np.array(d_gexc['above'])
    l_gexc_yrs.append(np_gexc)
    ic(np.shape(l_gexc_yrs))
np_gexc_yrll = np.stack(l_gexc_yrs, axis=0)
ds_gexc = xr.Dataset(
    data_vars = {
        "gexc": (("year", "realization", "lat", "lon"), np_gexc_yrll),
    },
    coords = {
        "year": years,
        "realization": da_rlz.realization.data,
        "lon": da_rlz.lon.data,
        "lat": da_rlz.lat.data,
    },
    attrs = {
        "long_name": 'gexc'
    }
)
ic(ds_gexc)
if out_flag:
    ic(out_path + out_name)
    ds_gexc.to_netcdf(out_path + out_name)