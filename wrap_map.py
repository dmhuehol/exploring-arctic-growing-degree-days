''' wrap_map
Open LENS2 data, manage realizations, select slice, and plot map using 
fun_plots. Written to handle processes manually to facilitate 
exploratory analysis of an arbitrary variable.
'''
import sys

from icecream import ic
import matplotlib.cm as mcm
import matplotlib.pyplot as plt

import classes_gddt as cg
import fun_process as fproc
import fun_plots as fpl

dp = cg.DataParams(
    path='/Users/dhueholt/Documents/gddt_data/LENS2/monthly_SST/AMJJAS/',
    tok='*.nc', var='SST', flag_raw_ds=True, flag_raw_da=True, 
    flag_time_slice=True, flag_manage_rlz=True, flag_land_mask=True, 
    flag_roi=False)
setp = cg.SetParams(mask_flag='ocean', yrs=[2000, 2010], rlz='mean')
ppar = cg.PlotParams(
    cb_bool=True, cb_extent='neither', cb_label='deg C', cb_vals=[-2, 35],
    cmap=mcm.turbo, dpi=400, figsize=(6, 3), o_bool=True,
    o_path='/Users/dhueholt/Documents/gddt_fig/20250217_wrapmap/',
    o_name='LENS2_amjjas-sst_global_' + str(setp.yrs[0]) 
        + str(setp.yrs[1]-1) + 'rlzmn', proj='EqualEarth180',)
d_open = fproc.common_opener(dp, setp)

#  Define data to plot. Requires hand-tuning for variable of interest.
da_rlz_mn = d_open['land_mask'].compute()
da_timerlz_mn = da_rlz_mn.mean(dim='time')
da_timerlz_mn_celsius = da_timerlz_mn - 273.15
plot_this = da_timerlz_mn_celsius.squeeze()

plt.rcParams.update({'font.family': 'Catamaran'})
#  'font.weight': normal, bold, heavy, light, ultrabold, ultralight
plt.rcParams.update({'font.weight': 'light'})
plt.rcParams.update({'font.size': 12})
fpl.plot_globe(plot_this, ppar)