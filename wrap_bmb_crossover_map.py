'''wrap_bmb_crossover_map
Plot difference in crossover maps between biomass burning ensembles.

This requires separate pre-calculated gridded files of crossover data
with dimensions time, lat, lon (and realizations optional) for both the
CMIP6 forcings and the smoothed biomass burning forcings. To make these
files, run wrap_calc_exceedance_grid and wrap_calc_crossover_grid.
'''
import sys

from icecream import ic
from matplotlib import colormaps as cm
import matplotlib.pyplot as plt

import classes_gddt as cg
import fun_calc_var as fcv
import fun_plots as fpl
import fun_process as fproc

cmn_path = '/Users/danielhueholt/Documents/Figures/arc-gdd_fig/20260116_cont/'
dp_crossover_cmip6 = cg.DataParams(
    path='/Users/danielhueholt/Data/gddt_data/LENS2/exceedance/crossover/',
    tok='forcingcmip6_forcedcrossover_threshold80.nc', var='crossover',
    flag_raw_ds=True, flag_raw_da=True, flag_time_slice=False,
    flag_manage_rlz=True, flag_land_mask=True, flag_roi=True)
dp_crossover_smoothed = cg.DataParams(
    path='/Users/danielhueholt/Data/gddt_data/LENS2/exceedance/crossover/',
    tok='forcingsmoothed_forcedcrossover_threshold80.nc', var='crossover',
    flag_raw_ds=True, flag_raw_da=True, flag_time_slice=False,
    flag_manage_rlz=True, flag_land_mask=True, flag_roi=True)
setp_crossover = cg.SetParams(
    area_stat='pass', base_yrs=[0, 2000],
    mask='/Users/danielhueholt/Data/gddt_data/mask/cesm_atm_mask.nc',
    mask_flag='land', reg_oi='global', rlz='all', yrs=[1850, 2100])
ppar_crossover = cg.PlotParams(
    cb_bool=True, cb_extent='neither', cb_label='crossover year',
    cb_vals=[-20, 20],
    # cmap=fpl.crossover_n(n=10), dpi=800,
    cmap=cm['RdBu'], dpi=800,
    figsize=(5,4), o_bool=True, o_name='', o_path=cmn_path, o_prefix='',
    plot_crossover_dict=dict(
        forced_dict=dict(bool=True, threshold=(80,)),
        member_dict=dict(bool=False)),
    plot_each_member=False, proj='Arctic', quantile=None, title='',
    title_size=10)

crossover_dict_cmip6 = fproc.common_opener(dp=dp_crossover_cmip6, setp=setp_crossover)
da_crossover_cmip6 = crossover_dict_cmip6['land_mask'].compute()
crossover_attrs_cmip6 = crossover_dict_cmip6['raw_ds'].attrs
crossover_dict_smoothed = fproc.common_opener(dp=dp_crossover_smoothed, setp=setp_crossover)
da_crossover_smoothed = crossover_dict_smoothed['land_mask'].compute()
crossover_attrs_smoothed = crossover_dict_smoothed['raw_ds'].attrs
if ppar_crossover.plot_crossover_dict['forced_dict']['bool']:
    msg_crossover = 'Plotting forced crossover at ' \
        + str(crossover_attrs_cmip6['threshold']) + '% threshold'
    ppar_crossover.o_name = 'bmbdifference_forcedcrossover_LENS2_threshold' \
        + str(crossover_attrs_cmip6['threshold'])
    ppar_crossover.title = 'CMIP6 - Smoothed forced crossover \n (ensemble mean >' \
        + str(crossover_attrs_cmip6['threshold']) + '% of Preindustrial samples)'
else:
    raise ValueError(
        'Check inputs! Ensure forced crossover are not None.')
ic(msg_crossover)
name_dict = fproc.namer(crossover_dict_cmip6["raw_da"], setp_crossover)

ppar_crossover.o_name = ppar_crossover.o_prefix + ppar_crossover.o_name
ic(ppar_crossover.o_name, ppar_crossover.title)
plt.rcParams.update({'font.family': 'Catamaran'})
#  'font.weight': normal, bold, heavy, light, ultrabold, ultralight
plt.rcParams.update({'font.weight': 'normal'})
plt.rcParams.update({'font.size': 12})
plot_this = da_crossover_cmip6 - da_crossover_smoothed
fcv.check_stats(plot_this)
fig, ax = fpl.plot_globe(plot_this, ppar_crossover)
