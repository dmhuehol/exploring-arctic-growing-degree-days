'''wrap_ens_guide
Wraps the prep_guide and guide in fun_process to return ensemble member 
indices corresponding to various properties for input region and time.

'guide_by': Choose property to guide by over period of interest.
Valid inputs are:
    'mean': Mean
    'stdev': Standard deviation
    'max': Absolute maximum
    'min': Absolute minimum
    'trend': Linear trend

TO BE ADDED
properties: 
stats:

Written by Daniel Hueholt
Graduate Research Assistant
'''
import sys

from icecream import ic
import numpy as np
from scipy import stats
import xarray as xr

import fun_process as fproc
import gddt_region_library as g_rlib

d_d = {
    "p": '/Users/dhueholt/Documents/gddt_data/gdd/nomask/',
    "tok": '*arc*.nc',
    "var": 'gdd5_sum'
}
set_d = { #'all', True to plot all + mean as bonus rlz
    "yrs": [1930, 1949],
    "mask_flag": None,
    "mask": '/Users/dhueholt/Documents/gddt_data/mask/cesm_atm_mask.nc',
    "reg_oi": g_rlib.Tasiilaq(),
    "area_stat": 'mean',
    "qoi": 0.25, # Quantile value for guide function
    "guide_by": 'trend', # Choose property to guide by
}
d_d = d_d

ds_in = xr.open_mfdataset(
    d_d['p'] + d_d['tok'], concat_dim='realization', combine='nested', 
    chunks={'time': 10000}, coords='minimal')
da_in = ds_in[d_d['var']]
da_in = da_in.squeeze()
time_wheel = slice(set_d["yrs"][0], set_d["yrs"][1])
da_toi = da_in.sel(year=time_wheel)
da_mask = fproc.apply_mask(da_toi, set_d)
da_roi, loc_str, _ = fproc.manage_area(da_mask, set_d)

da_guide = fproc.prep_guide(da_roi, set_d)
ind_dict = fproc.guide(da_guide, set_d)
ic(ind_dict)
