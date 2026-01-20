"""wrap_trend_hist
Plot histograms of trends from GHCN stations and CESM. Count data beyond
threshold and
"""
import sys
sys.path.append(
    "/Users/danielhueholt/Documents/GitHub/" \
    + "exploring-arctic-growing-degree-days/")

from icecream import ic
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pyfonts import load_google_font, set_default_font
import seaborn as sn

import classes_gddt as cg
import fun_process as fproc
import gddt_region_library as g_rlib

data_path = "/Users/danielhueholt/Data/gddt_data/trends/"
data_fn = "hist_GHCN_abv50N_cov97_span10_trend18732022_withLENS2.csv"
dp_model = cg.DataParams(
    path='/Users/danielhueholt/Data/gddt_data/gdd/lens2/', tok='*arc*.nc',
    var='gdd5_sum', flag_raw_ds=True, flag_raw_da=True, flag_time_slice=True,
    flag_manage_rlz=True, flag_land_mask=False, flag_roi=True)
setp_model = cg.SetParams(
    area_stat='pass', reg_oi=g_rlib.Arctic50N(), rlz='all', yrs=[1873, 2022])
msg_dist_bool = True
threshold = -50
threshold_type = '<'
threshold_stat = 'prob'
x_max_mag = 35 # 160 for 10-year, 35 for 30-year
ppar = cg.PlotParams(
    bw=5, color_comp='#6f8c31', color_r='#4e318c', dpi='.pdf', figsize=(7, 7),
    font='Catamaran', leg_dict={"bool": True, "frameon": False, "size": 22},
    o_name=data_fn.replace('.csv', ''),
    o_path='/Users/danielhueholt/Documents/Figures/arc-gdd_fig/20260116_cont/',
    o_prefix='', title='Histogram of all trends', y_lim=[0, 0.007])
use_font = load_google_font(ppar.font, danger_not_verify_ssl=True)
set_default_font(use_font)

fn_pcs = data_fn.split('_')
df_trends = pd.read_csv(data_path + data_fn)
unique_intervals = sorted([*{*df_trends["interval"]}])
all_keys = df_trends.keys()
lens2_keys = [k for k in all_keys if "LENS2" in k]

np_all_trends = df_trends["slope"].to_numpy()
np_all_trends_lens2 = df_trends[lens2_keys].to_numpy()
np_ravel_trends_lens2 = np.ravel(np_all_trends_lens2)

fig, ax = plt.subplots(figsize=ppar.figsize)
ax.spines[['right', 'top']].set_visible(False)
sn.histplot(
    x=np_all_trends, binwidth=ppar.bw, color=ppar.color_comp,
    stat='probability', label='GHCN', alpha=0.5)
sn_hist_lens2 = sn.histplot(
    x=np_ravel_trends_lens2, binwidth=ppar.bw, color=ppar.color_r,
    stat='probability', label='LENS2', alpha=0.5)
if ppar.leg_dict["bool"]:
    plt.legend(
        fontsize=ppar.leg_dict["size"], frameon=ppar.leg_dict["frameon"])
if threshold is not None:
    if threshold > 0:
        ppar.o_prefix += 'high-extreme_'
        plt.xlim([threshold, x_max_mag])
    else:
        ppar.o_prefix += 'low-extreme_'
        plt.xlim([-1 * x_max_mag, threshold])
    #  Report number/fraction of values past threshold to shell
    past_threshold_ghcn = fproc.count_threshold(
        np_all_trends, threshold, type=threshold_type, stat='fraction')
    past_threshold_lens2 = fproc.count_threshold(
        np_ravel_trends_lens2, threshold, type=threshold_type, stat='fraction')
    ic(past_threshold_ghcn, past_threshold_lens2)
if ppar.y_lim != 'auto':
    plt.ylim(ppar.y_lim)
if ppar.o_name is None:
    raise TypeError("Must specify output filename (o_name is None)")
if isinstance(ppar.dpi, str):
   if 'pdf' in ppar.dpi:
       ppar.o_name += '.pdf'
plt.savefig(ppar.o_path + ppar.o_prefix + ppar.o_name)
