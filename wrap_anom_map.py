'''wrap_anom_map
Script to calculate anomalies against base period and 
wrap map plotting functions in fun_plots.

'guide_by': Choose property to guide by over period of interest.
Valid inputs are:
    'period_mean': Mean
    'period_stdev': Standard deviation
    'period_max': Absolute maximum
    'period_min': Absolute minimum
    'period_trend': Linear trend

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University '''
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

cmn_path = '/Users/dhueholt/Documents/gddt_fig/20241208_furtheraes/'
cmn_yrs = [2002, 2011]
dp_cmpst = cg.DataParams(
    path='/Users/dhueholt/Documents/gddt_data/LENS2/monthly_SST/AMJJAS/',
    tok='*.nc', var='SST', flag_raw_ds=True, flag_raw_da=True, 
    flag_time_slice=True, flag_manage_rlz=True, flag_land_mask=False, 
    flag_roi=False)
setp_cmpst = cg.SetParams(
    dims='realization', mask_flag=None, yrs=cmn_yrs, 
    yrs_rel_to=cmn_yrs, rlz=24, z_flag=True)
dp_base = cg.DataParams(
    path='/Users/dhueholt/Documents/gddt_data/LENS2/monthly_SST/AMJJAS/',
    tok='*.nc', var='SST', flag_raw_ds=True, flag_raw_da=True, 
    flag_time_slice=True, flag_manage_rlz=True, flag_land_mask=False, 
    flag_roi=False)
setp_base = cg.SetParams(
    mask_flag=None, yrs=cmn_yrs, yrs_rel_to=cmn_yrs, rlz='all',
    z_flag=True)
ppar_cmpst = cg.PlotParams(
    cb_bool=True, cb_extent='neither', cb_label='auto', cb_vals=[-1.5, 1.5],
    cmap=fpl.balance_n(n=18), dpi=400, figsize=(10, 4), 
    o_bool=False, o_path=cmn_path, o_prefix='', plot_all=False, 
    proj='EqualEarth180', title_size=9)
ppar_reg = cg.PlotParams(
    alpha=0.53, cb_bool=False, cb_ticks=[-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6], 
    cb_vals=[-0.6, 0.6], cmap=mcm.Greys, color='#000000', dpi=400, 
    edgecolors='#000000', marker_size=30, o_bool=True, o_path=cmn_path, 
    o_prefix='', plot_all=False, proj='EqualEarth180', set_bad=False, 
    title_size=9)
setp_reg = cg.SetParams(reg_oi=g_rlib.BrooksRange_colonist())

# Set up data
# ic('Before opening')
d_cmpst = fproc.common_opener(dp_cmpst, setp_cmpst)
da_cmpst = d_cmpst['manage_rlz']
time_slice_base = slice(
    cftime.DatetimeNoLeap(setp_cmpst.yrs_rel_to[0], 1, 1, 0, 0, 0, 0),
    cftime.DatetimeNoLeap(setp_cmpst.yrs_rel_to[1], 12, 31, 14, 24, 0, 0)
)
d_base = fproc.common_opener(dp_base, setp_base)
da_base = d_base['manage_rlz']
da_anom = fproc.calc_anomaly(da_cmpst, da_base, setp_cmpst)

plot_this = da_anom.mean(dim='time').squeeze()
ppar_cmpst.title = 'SST anomalies for ' + str(cmn_yrs[0]) + '-' \
    + str(cmn_yrs[1]) + ', member ' + str(setp_cmpst.rlz)
ppar_cmpst.o_name = 'SSTanom_' + str(cmn_yrs[0]) + '-' + str(cmn_yrs[1]) \
    + '_rlz' + str(setp_cmpst.rlz)
plt.rcParams.update({'font.family': 'Catamaran'})
#  'font.weight': normal, bold, heavy, light, ultrabold, ultralight
plt.rcParams.update({'font.weight': 'light'})
plt.rcParams.update({'font.size': 12})
if ppar_reg.o_bool:
    ppar_reg.o_name = ppar_cmpst.o_name
    ppar_reg.title = ppar_cmpst.title
    fig, ax_comp = fpl.plot_globe(plot_this, ppar_cmpst)
    if len(setp_reg.reg_oi["reg_lats"]) > 1:
        reg_ones = fpl.mask_region(setp_reg.reg_oi)
        fpl.plot_globe(reg_ones, ppar_reg, ax=ax_comp)
    else:
        fpl.plot_globe_ng(
            ppar_reg.color, setp_reg.reg_oi["reg_lats"], 
            setp_reg.reg_oi["reg_lons"], ppar_reg, ax=ax_comp)
else:
    fpl.plot_globe(plot_this, ppar_cmpst)