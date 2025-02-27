'''wrap_anom_map
Script to calculate anomalies against base period and plot as map.
'''
import sys

import cftime
from icecream import ic
import matplotlib.cm as mcm
import matplotlib.pyplot as plt
import xarray as xr

import classes_gddt as cg
import fun_process as fproc
import fun_plots as fpl
import gddt_region_library as g_rlib
xr.set_options(keep_attrs=True)

cmn_path = '/Users/dhueholt/Documents/gddt_fig/20250227_continueRemaking/'
cmn_yrs = [2082, 2091]
dp = cg.DataParams(
    path='/Users/dhueholt/Documents/gddt_data/LENS2/monthly_SST/AMJJAS/',
    tok='*.nc', var='SST', flag_raw_ds=True, flag_raw_da=True, 
    flag_time_slice=True, flag_manage_rlz=True, flag_land_mask=False, 
    flag_roi=False)
setp = cg.SetParams(
    dims='realization', mask_flag=None, yrs=cmn_yrs, 
    yrs_rel_to=cmn_yrs, rlz=76, z_flag=True)
dp_base = cg.DataParams(
    path='/Users/dhueholt/Documents/gddt_data/LENS2/monthly_SST/AMJJAS/',
    tok='*.nc', var='SST', flag_raw_ds=True, flag_raw_da=True, 
    flag_time_slice=True, flag_manage_rlz=True, flag_land_mask=False, 
    flag_roi=False)
setp_base = cg.SetParams(
    mask_flag=None, yrs=cmn_yrs, yrs_rel_to=cmn_yrs, rlz='all',
    z_flag=True)
ppar = cg.PlotParams(
    cb_bool=True, cb_extent='neither', cb_label='auto', cb_vals=[-1.5, 1.5],
    cmap=fpl.balance_n(n=18), dpi=400, figsize=(10, 4), 
    o_bool=False, o_path=cmn_path, o_prefix='', plot_each_member=False, 
    proj='EqualEarth180', title_size=9)
ppar_region = cg.PlotParams(
    alpha=0.53, cb_bool=False, cb_ticks=[-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6], 
    cb_vals=[-0.6, 0.6], cmap=mcm.Greys, color='#000000', dpi=400, 
    edgecolors='#000000', marker_size=30, o_bool=True, o_path=cmn_path, 
    o_prefix='', plot_each_member=False, proj='EqualEarth180', set_bad=False, 
    title_size=9)
setp_reg = cg.SetParams(reg_oi=g_rlib.BrooksRange_colonist())

# Set up data
# ic('Before opening')
dict_data = fproc.common_opener(dp, setp)
da_rlzoi = dict_data['manage_rlz']
time_slice_base = slice(
    cftime.DatetimeNoLeap(setp.yrs_rel_to[0], 1, 1, 0, 0, 0, 0),
    cftime.DatetimeNoLeap(setp.yrs_rel_to[1], 12, 31, 14, 24, 0, 0)
)
dict_data_base = fproc.common_opener(dp_base, setp_base)
da_base_rlzoi = dict_data_base['manage_rlz']
da_anom = fproc.calc_anomaly(da_rlzoi, da_base_rlzoi, setp)

plot_this = da_anom.mean(dim='time').squeeze()
ppar.title = 'SST anomalies for ' + str(cmn_yrs[0]) + '-' \
    + str(cmn_yrs[1]) + ', member ' + str(setp.rlz)
ppar.o_name = 'SSTanom_' + str(cmn_yrs[0]) + '-' + str(cmn_yrs[1]) \
    + '_rlz' + str(setp.rlz)
plt.rcParams.update({'font.family': 'Catamaran'})
#  'font.weight': normal, bold, heavy, light, ultrabold, ultralight
plt.rcParams.update({'font.weight': 'light'})
plt.rcParams.update({'font.size': 12})
if ppar_region.o_bool:
    ppar_region.o_name = ppar.o_name
    ppar_region.title = ppar.title
    fig, ax_comp = fpl.plot_globe(plot_this, ppar)
    if len(setp_reg.reg_oi["reg_lats"]) > 1:
        reg_ones = fpl.mask_region(setp_reg.reg_oi)
        fpl.plot_globe(reg_ones, ppar_region, ax=ax_comp)
    else:
        fpl.plot_globe_ng(
            ppar_region.color, setp_reg.reg_oi["reg_lats"], 
            setp_reg.reg_oi["reg_lons"], ppar_region, ax=ax_comp)
else:
    fpl.plot_globe(plot_this, ppar)