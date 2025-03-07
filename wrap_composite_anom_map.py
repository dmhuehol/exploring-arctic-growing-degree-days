'''wrap_composite_anom_map
Script to composite and plot anomalies in one field based on properties
of another guiding variable. Optionally, pairs the composite plot with a 
timeseries of the guiding variables in the relevant time period and 
region.

This is configured to make composites of SST based on properties of 
growing degree days (GDD). The script includes options to plot 
composites of sea level pressure (PSL) and ice fraction (ICEFRAC) as 
well. It can be straightforwardly extended to other variables by 
obtaining the appropriate data and modifying the DataParams and 
PlotParams instances.

The GuideParams class is used to determine the properties of the 
composite.
'guide_by': Choose property to guide by over period of interest.
Valid inputs are:
    'mean': Mean
    'stdev': Standard deviation
    'max': Absolute maximum
    'min': Absolute minimum
    'trend': Linear trend
See the documentation in classes_gddt.GuideParams for further details.
'''
import sys

from icecream import ic
import matplotlib.cm as mcm

import classes_gddt as cg
import fun_process as fproc
import fun_plots as fpl
import gddt_region_library as g_rlib
import puppets

#  type: 'coe_hpc' or 'local'
dp_gdd, dp_gdd_alltimes, dp_psl, dp_sst, dp_icefrac, cmn_path = fproc.get_params(
    type='local', cmn_path='')
ip = cg.IntervalParams(span=10, strt_yr=1850, end_yr=2100, type='rolling')
setp_gdd = cg.SetParams(
    area_stat='mean', mask_flag='none', reg_oi=g_rlib.BrooksRange_colonist(), 
    rlz='all', yrs=[1850, 1859])
setp_psl = cg.SetParams(
    dims=['time', 'realization'], mask_flag=None, rlz='all', yrs=[1850, 1859], 
    yrs_rel_to=[1850, 1859], z_flag=True)
setp_sst = cg.SetParams(
    dims=['time', 'realization'], mask_flag=None, rlz='all', yrs=[1850, 1859], 
    yrs_rel_to=[1850, 1859], z_flag=True)
gp_gdd = cg.GuideParams(guide_by='mean', composite_key='min20', qoi=0.25)
ppar_base_composite = cg.PlotParams(
    cb_bool=True, cb_extent='neither', cb_label='auto', 
    cb_ticks=[-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6], cb_vals=[-0.6, 0.6],
    cmap=fpl.balance_n(n=18), dpi=400, figsize=(10, 4), 
    o_bool=False, o_path=cmn_path, o_prefix='', plot_each_member=False, 
    proj='EqualEarth180', title_size=9)
ppar_region = cg.PlotParams( # FRAME
    alpha=0.53, cb_bool=False, cb_ticks=[-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6], 
    cb_vals=[-0.6, 0.6], cmap=mcm.Greys, color='#000000', dpi=400, 
    edgecolors='#000000', marker_size=30, o_bool=True, o_path=cmn_path, 
    o_prefix='', plot_each_member=False, proj='EqualEarth180', set_bad=False, 
    title_size=9)
ppar_super = cg.PlotParams(
    cb_bool=True, cb_extent='neither', cb_label='auto', 
    cb_ticks=[-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6], cb_vals=[-0.6, 0.6],
    cmap=fpl.balance_n(n=18), dpi=400, figsize=(10, 4), o_bool=False,
    o_path=cmn_path, o_prefix='super_', plot_each_member=False, 
    proj='EqualEarth180', title_size=9)
ppar_gdd_ts_gray = cg.PlotParams(
    color='#ff80ed', color_r='#c4c4c4', label='LENS2 members', o_bool=False, 
    lw=1, o_path=cmn_path, x_label='year', x_lim=[1850, 2100], 
    y_label='GDD base 5', y_lim='auto')
ppar_gdd_ts = cg.PlotParams( # FRAME
    color='#ff80ed', color_r='#ffa6f2', dpi=400,
    label='LENS2 members: ' + gp_gdd.composite_key,
    o_bool=True, o_path=cmn_path, x_label='year', x_lim=[1850, 2100], 
    y_label='GDD base 5', y_lim='auto')
super_interval = 30
pair_ts = True
avg_all_composites = True

puppets.paired_composites_ts(
    dp_guide=dp_gdd, setp_guide=setp_gdd, dp_composite=dp_sst, 
    setp_composite=setp_sst, dp_guide_alltimes=dp_gdd_alltimes, ip=ip, 
    gp=gp_gdd, ppar=ppar_base_composite, ppar_region=ppar_region, 
    ppar_super=ppar_super, ppar_guide_ts_all=ppar_gdd_ts_gray, 
    ppar_guide_ts=ppar_gdd_ts, pair_timeseries=pair_ts, 
    super_interval=super_interval, avg_all_composites=avg_all_composites)