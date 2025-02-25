'''wrap_multi_composite_anom_map
Script to composite and plot anomalies in one field based on properties
of another guiding variable. Optionally, pairs the composite plot with a 
timeseries of the guiding varaiables in the relevant time period and 
region. 

This is exactly the same as wrap_composite_anom_map, except it makes
composites in parallel using multiprocessing to allow more figures to be
generated at once. It is configured to make composites of SST based on 
properties of growing degree days (GDD). The script includes options to plot 
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
from multiprocessing import Process
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
    type='coe_hpc', cmn_path='')
ip = cg.IntervalParams(span=10, strt_yr=1850, end_yr=2100, type='rolling')
setp_gdd = cg.SetParams(
    area_stat='mean', mask_flag='none', reg_oi=(
        g_rlib.Tasiilaq(), g_rlib.UpperQuebec(), 
        g_rlib.SiberiaNorthBorealForest(),
        ),
    rlz='all', yrs=[1850, 1859])
setp_psl = cg.SetParams(
    dims=['time', 'realization'], mask_flag=None, rlz='all', yrs=[1850, 1859], 
    yrs_rel_to=[1850, 1859], z_flag=True)
setp_sst = cg.SetParams(
    dims=['time', 'realization'], mask_flag=None, rlz='all', yrs=[1850, 1859], 
    yrs_rel_to=[1850, 1859], z_flag=True)
gp_gdd = cg.GuideParams(guide_by='mean', composite_key=(
    'max20', 'min20',
    ), qoi=0.25)
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
    label='LENS2 members: ', #composite_key added in loop
    o_bool=True, o_path=cmn_path, x_label='year', x_lim=[1850, 2100], 
    y_label='GDD base 5', y_lim='auto')
super_interval = 30
pair_ts = True
avg_all_composites = True

if __name__=='__main__':
    for roi in setp_gdd.reg_oi:
        loop_setp_gdd = cg.SetParams(
            dims=['time', 'realization'], mask_flag=None, rlz='all', 
            reg_oi=roi, yrs=[1850, 1859], yrs_rel_to=[1850, 1859], z_flag=True)
        for coi in gp_gdd.composite_key:
            loop_gp_gdd = cg.GuideParams(
                guide_by='mean', composite_key=coi, qoi=0.25)
            ppar_gdd_ts.label = ppar_gdd_ts.label + loop_gp_gdd.composite_key
            ic(loop_setp_gdd.reg_oi, loop_gp_gdd.composite_key, ppar_gdd_ts.label)
            keyword_args = {
                "dp_guide": dp_gdd,
                "setp_guide": loop_setp_gdd,
                "dp_composite": dp_sst,
                "setp_composite": setp_sst,
                "dp_guide_alltimes": dp_gdd_alltimes,
                "ip": ip,
                "gp": loop_gp_gdd,
                "ppar": ppar_base_composite,
                "ppar_region": ppar_region,
                "ppar_super": ppar_super,
                "ppar_guide_ts_all": ppar_gdd_ts_gray,
                "ppar_guide_ts": ppar_gdd_ts,
                "pair_timeseries": pair_ts,
                "super_interval": super_interval,
                "avg_all_composites": avg_all_composites,
            }
            p = Process(
                target=puppets.paired_composites_ts, kwargs=keyword_args)
            p.start()
            ppar_gdd_ts.label = ppar_gdd_ts.label.replace(
                loop_gp_gdd.composite_key, '')