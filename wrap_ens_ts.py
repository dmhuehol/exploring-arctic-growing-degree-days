'''wrap_ens_ts
Script to wrap ensemble timeseries plotting functions in fun_plots. 

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''
import sys

from icecream import ic
import matplotlib.pyplot as plt
import xarray as xr

import classes_gddt as cg
import fun_process as fproc
import fun_plots as fpl
import gddt_region_library as g_rlib
xr.set_options(keep_attrs=True)

plot_stdev = False
dp_gdd = cg.DataParams(
    path='/Users/dhueholt/Documents/gddt_data/gdd/preindustrial/', 
    tok='*arc*.nc', var='gdd5_sum', flag_raw_ds=True, 
    flag_raw_da=True, flag_time_slice=False, flag_manage_rlz=True, 
    flag_land_mask=False, flag_roi=True)
setp_gdd = cg.SetParams(
    area_stat='mean', mask_flag='none', reg_oi=g_rlib.BelangerIslandCell(), rlz='all', 
    yrs=[0, 2000])
ppar = cg.PlotParams(
    # color='#ff80ed', color_r='#ffa6f2', dpi=400, label='LENS2', 
    color='#ffa6f2', color_r='#ffa6f2', dpi=400, label='Preindustrial', 
    mn_bool=True, o_bool=True, plot_all=False,
    o_path='/Users/dhueholt/Documents/gddt_fig/20241029_gddts/',
    o_prefix='',
    ts_type = 'spaghetti', x_label='year', x_lim=[0, 2000], 
    y_label='GDD', y_lim=[0, 2021.69])#y_lim='auto')

d_data = fproc.common_opener(dp=dp_gdd, setp=setp_gdd)
da_roi = d_data['roi']
if plot_stdev:
    da_roi = da_roi.std(dim='realization')
plot_this = da_roi
x_data = da_roi.year
name_dict = fproc.namer(da_roi, setp_gdd)
ppar.title = name_dict['a_title_norlz']
ppar.o_name = name_dict['a_name_norlz']

plt.rcParams.update({'font.family': 'Catamaran'})
#  'font.weight': normal, bold, heavy, light, ultrabold, ultralight
plt.rcParams.update({'font.weight': 'light'})
plt.rcParams.update({'font.size': 12})
if ppar.plot_all:
    for rlzc, rlz in enumerate(plot_this.realization.data):
        act_this = plot_this.isel(realization=rlz)
        act_mem = da_roi.attrs['members'][rlzc]
        org_name = ppar.o_name
        org_title = ppar.title
        ppar.o_name = ppar.o_name + '_rlz' + str(act_mem)
        ppar.title = org_title + ' rlz ' + str(act_mem)
        fpl.plot_timeseries(act_this, x_d=x_data, ppar=ppar) # Only plot that makes sense for single realization
        plt.close()
        ppar.o_name = org_name
        ppar.title = org_title
else:
    # ppar.title = name_dict['a_title'] + ' ' + setp_gdd.reg_oi['reg_str']
    if plot_stdev:
        ppar.title = 'LENS2 GDD standard dev ' + setp_gdd.reg_oi['reg_str']
    else:
        ppar.title = 'LENS2 GDD ' + setp_gdd.reg_oi['reg_str']
    ppar.o_name = ppar.o_prefix + name_dict['a_name']
    if 'realization' not in plot_this.dims:
        fpl.plot_timeseries(plot_this, x_d=x_data, ppar=ppar)
    elif ppar.ts_type == 'spread':
        raise NotImplementedError('Spread not implemented')
        # fpl.plot_timeseries_spread(plot_this, x_data, ppar)
    elif ppar.ts_type == 'spaghetti':
        fpl.plot_timeseries_spaghetti(plot_this, x_data, ppar)
    else:
        ic('Check plot type!')