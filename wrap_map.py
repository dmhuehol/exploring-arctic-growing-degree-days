'''wrap_map
Script to wrap map plotting functions in fun_plots. '''

import sys

import cartopy
from icecream import ic
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

import fun_process as fproc
import fun_plots as fpl

d_d = {
    "p": '/Volumes/Polycrystal/Data/LENS2/daily_TREFHT/refined/',
    "tok": '*antar*.nc',
    "var": 'TREFHT'
}
set_d = { #'all', True to plot all + mean as bonus rlz
    'yrs': [1850, 2100],
    'rlz': 'mean',
    'plot_all': False
}
plot_d = {
    "o_path": '/Users/dhueholt/Documents/gddt_fig/20240509_rlzanomera/',
    "o_name": 'LENS2_daily2mtemp_Antarctic_' + str(set_d['yrs'][0]) 
        + str(set_d['yrs'][1]-1) + 'rlz',
    "figsize": (6, 3), #(5,5) (10,4)
    "proj": 'Antarctic',
    "cmap": 'turbo',
    "cb_bool": True,
    "cb_vals": [-25, 25],
    "cb_ticks": [-25, -15, -5, 0, 5, 15, 25],
    "cb_extent": 'neither',
    "cb_label": '\u00b0C',
    "title": 'LENS2: Daily 2m temperature mn ' + str(set_d['yrs'][0]) 
        + '-' + str(set_d['yrs'][1]),
    "dpi": 400
}
########################################################################
ds_in = xr.open_mfdataset(
    d_d['p'] + d_d['tok'], concat_dim='realization', combine='nested', 
    chunks={'time': 10000}, coords='minimal')
da_in = ds_in[d_d['var']]
da_in = da_in.squeeze()
da_rlz = da_in
if 'realization' in da_in.dims:
    da_rlz = fproc.manage_rlz(da_in, set_d)
da_time = da_rlz.mean(dim='time')
plot_this = da_time

plt.figure()
plt.rcParams.update({'font.family': 'Catamaran'})
plt.rcParams.update({'font.weight': 'light'}) #normal, bold, heavy, light, ultrabold, ultralight
plt.rcParams.update({'font.size': 12})
fpl.plot_globe(plot_this, plot_d)