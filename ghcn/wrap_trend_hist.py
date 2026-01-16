"""wrap_trend_hist
Plot histograms of trends from GHCN stations and CESM.
"""
import sys
sys.path.append(
    "/Users/danielhueholt/Documents/GitHub/" \
    + "exploring-arctic-growing-degree-days/")

from icecream import ic
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sn

import classes_gddt as cg
import fun_calc_var as fcv
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
fn_pcs = data_fn.split('_')
msg_dist_bool = True
ppar = cg.PlotParams(
    bw=5, figsize=(7, 7),
    o_name=data_fn.replace('.csv', ''),
    o_path='/Users/danielhueholt/Documents/Figures/arc-gdd_fig/20260115_cont/',
    title='Histogram of all trends')

df_trends = pd.read_csv(data_path + data_fn)
unique_intervals = sorted([*{*df_trends["interval"]}])
all_keys = df_trends.keys()
lens2_keys = [k for k in all_keys if "LENS2" in k]

np_all_trends = df_trends["slope"].to_numpy()
np_all_trends_lens2 = df_trends[lens2_keys].to_numpy()
np_ravel_trends_lens2 = np.ravel(np_all_trends_lens2)

fig, ax = plt.subplots(figsize=ppar.figsize)
plt.rcParams.update({"font.family": 'Open Sans'})
plt.rcParams.update({"font.weight": 'light'})
plt.rcParams.update({"font.size": 10})
sn.histplot(
    x=np_all_trends, binwidth=ppar.bw, color='#8eb8ad', stat='probability',
    label='GHCN')
sn.histplot(
    x=np_ravel_trends_lens2, binwidth=ppar.bw, color='#69838d',
    stat='probability', label='LENS2')
plt.legend()
if ppar.o_name is None:
    raise TypeError("Must specify output filename (o_name is None)")
plt.savefig(ppar.o_path + ppar.o_name)
