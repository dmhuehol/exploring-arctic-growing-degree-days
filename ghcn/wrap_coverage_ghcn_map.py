########################################################################
###############################################################################
''' wrap_coverage_map
Plot map of stations shaded by coverage.

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''
from glob import glob
import sys
sys.path.append(
    "/Users/dhueholt/Documents/Github/" 
    + "growing-degree-days-treelines-poles/")

import cmocean as cmo
from icecream import ic
import matplotlib.cm as mcm
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
xr.set_options(keep_attrs=True)
import pandas as pd

import fun_ghcn as fg
import fun_plots as fpl

yr = 30
sc_d = {
    "sc_f": "~/Documents/ghcn_data/spancov_lists/spancov_Arctic50N_gt1.csv",
    "span_thr": 365 * yr + np.floor(yr / 4),
    "cov_thr": 90,
    #  cov_type: '>' for greater than, '<' for less than, '=' for equal to
    "cov_type": '>',
}
df_sc = fg.get_spancov(sc_d)
n_st = ic(len(df_sc))
plot_d = {
    "plot_all": False,
    "o_path": '/Users/dhueholt/Documents/gddt_fig/20240618_slpCompAndMoreCompAndGHCN/',
    "o_name": 'GHCN_above50N_cov' + str(sc_d['cov_thr']) + '_span' + str(yr),
    "o_prefix": '6_cmap_',
    "figsize": (7, 6), #(6,3) (5,5) (10,4)
    "proj": 'Arctic', #'Arctic', 'Antarctic'
    "cmap": cmo.cm.algae,
    "cb_bool": True,
    #  "cb_vals": [min, max] or 'auto'
    "cb_vals": [sc_d['cov_thr'], 100],
    "cb_ticks": None,
    "cb_extent": 'neither',
    "cb_label": '% cov',
    "title": 'GHCN stations above 50N with >' + str(yr) + 'yr record' + ' n=' + str(n_st),
    "title_size": 14,
    "dpi": 400
}

df_to_plot = df_sc['coverage']
df_lat = df_sc['lat']
df_lon = df_sc['lon']
plt.figure()
plt.rcParams.update({'font.family': 'Catamaran'})
plt.rcParams.update({'font.weight': 'light'})
plt.rcParams.update({'font.size': 10})
plot_d['o_name'] = plot_d['o_prefix'] + plot_d['o_name']
fpl.plot_globe_ng(df_to_plot, df_lat, df_lon, plot_d)
