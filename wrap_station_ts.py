'''wrap_station_ts
Script to wrap timeseries plotting functions in fun_plots for
station data such as Danish Meteorological Institute station data. 

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
########################################################################
###############################################################################
'''
from glob import glob
import sys

from icecream import ic
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
xr.set_options(keep_attrs=True)

import fun_process as fproc
import fun_plots as fpl
import gddt_region_library as g_rlib

d_d = {
    "p": '/Users/dhueholt/Documents/gddt_data/gdd/dmi/',
    "tok": '*34360_4360*.nc', # 34216, 34360
    "var": 'gdd5_sum'
}
set_d = {
    "yrs": [2011, 2020],
    "bad_thr": 36,
}
d_d = d_d
files = glob(d_d['p'] + d_d['tok'])
# ex = fproc.pieces(files[0])
# ic(ex)
ds_in = xr.open_mfdataset(
    d_d['p'] + d_d['tok'], concat_dim='station', combine='nested', 
    coords='minimal')
da_in = ds_in[d_d['var']]
time_wheel = slice(set_d["yrs"][0], set_d["yrs"][1])
da_toi = da_in.sel(year=time_wheel)
da_squeeze = da_toi.squeeze()

da_ndd = ds_in['ndd'].compute().squeeze().sel(year=time_wheel)
np_ndd = da_ndd.data
da_ndd_bad_yr = da_ndd.year.data[np_ndd > set_d['bad_thr']]
# da_nobad = da_squeeze.drop_sel(year=da_ndd_bad_yr)
da_nobad = da_squeeze.where(np_ndd < set_d['bad_thr'])  

plot_this = da_nobad
# plot_this = da_nobad / da_nobad.max()
# plot_this = (da_nobad - da_nobad.mean(dim='year')) / da_nobad.std(dim='year')
x_data = da_nobad.year

plot_d = {
    # "plot_each": False, # 'station?' Plot each member separately
    # "type": 'spaghetti', # 'station?'
    "o_path": '/Users/dhueholt/Documents/gddt_fig/20240513_era5coh/',
    "o_name": '2_' + 'ts_34360-4360',
        # + ex['d_id'] + d_d['var'] + '_' + ex['loc'] + '_' 
        # + str(set_d['yrs'][0]) + str(set_d['yrs'][1]),
    "title": 'Tasiilaq blended annual GDD ' + str(set_d['yrs'][0]) 
        + '-' + str(set_d['yrs'][1]),
    # "title": ex['d_long'] + ' ' + ex['loc'] + ' ' + ex["var_long"] 
        # + ' ' + str(set_d['yrs'][0]) + '-' + str(set_d['yrs'][1]),
    "y_label": d_d['var'],
    "x_label": 'years',
    "x_lim": set_d['yrs'],
    "y_lim": None,
    "col": '#ff80ed',
    # "r_col": '#a7cb6d',
    "leg_flag": False,
    "mn_flag": True,
    "dpi": 400
}
plt.figure()
plt.rcParams.update({'font.family': 'Catamaran'})
plt.rcParams.update({'font.weight': 'light'}) #normal, bold, heavy, light, ultrabold, ultralight
plt.rcParams.update({'font.size': 12})
fpl.plot_timeseries(plot_this, x_data, plot_d)