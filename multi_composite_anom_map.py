'''multi_composite_anom_map
Script to composite and plot anomalies in a separate geophysical field 
based on growing degree days properties. Optionally, pairs the composite
plot with a timeseries of the growing degree days in the relevant time
period and region. Composites are made in parallel using multiprocessing
for efficiency.

'guide_by': Choose property to guide by over period of interest.
Valid inputs are:
    'mean': Mean
    'stdev': Standard deviation
    'max': Absolute maximum
    'min': Absolute minimum
    'trend': Linear trend

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University 
'''
from multiprocessing import Process
import sys

from icecream import ic
import matplotlib.cm as mcm
import matplotlib.pyplot as plt

import classes_gddt as cg
import fun_process as fproc
import fun_plots as fpl
import gddt_region_library as g_rlib
import puppets

#  type: 'coe_hpc' or 'local'
dp_gdd, dp_gdd_roi_alltimes, dp_psl, dp_sst, dp_icefrac, cmn_path = fproc.get_params(
    type='coe_hpc', cmn_path='/home/dhueholt/gddt_fig/frames/TEMP_PSL1/')
ip = cg.IntervalParams(span=1, strt_yr=1850, end_yr=2100, type='rolling')
setp_gdd = cg.SetParams(
    area_stat='mean', mask_flag='none', reg_oi=(
        g_rlib.Tasiilaq(), 
        ), rlz='all', 
    yrs=[1850, 1859])
setp_psl = cg.SetParams(
    mask_flag=None, yrs=[1850, 1859], yrs_rel_to=[1850, 1859], z_flag=True)
setp_sst = cg.SetParams(
    mask_flag=None, yrs=[1850, 1859], yrs_rel_to=[1850, 1859], z_flag=True)
gp_gdd = cg.GuideParams(guide_by='mean', cmpst_key=('max20', 'min20'))
ppar_comp = cg.PlotParams(
    cb_bool=True, cb_extent='neither', cb_label='auto', cb_vals=[-0.8, 0.8],
    cmap=fpl.balance_n(n=16), dpi=400, figsize=(10, 4), 
    o_bool=False, o_path=cmn_path, o_prefix='', plot_all=False, proj='EqualEarth180', 
    title_size=9)
ppar_reg = cg.PlotParams( # FRAME
    alpha=0.53, cb_bool=False, cb_vals=[-0.8, 0.8], cmap=mcm.Greys, color='#000000', 
    dpi=400, edgecolors='#000000', marker_size=30, o_bool=True, o_path=cmn_path, 
    o_prefix='', plot_all=False, proj='EqualEarth180', set_bad=False, title_size=9)
ppar_super = cg.PlotParams(
    cb_bool=True, cb_extent='neither', cb_label='auto', cb_vals=[-0.8, 0.8],
    cmap=fpl.balance_n(n=16), dpi=400, figsize=(10, 4), o_bool=True,
    o_path=cmn_path, o_prefix='super_', plot_all=False, proj='EqualEarth180',
    title_size=9)
ppar_gdd_ts_gray = cg.PlotParams(
    color='#ff80ed', color_r='#c4c4c4', label='LENS2 members', o_bool=False, 
    lw=1, o_path=cmn_path, x_label='year', x_lim=[1850, 2100], 
    y_label='GDD base 5', y_lim='auto')
ppar_gdd_ts = cg.PlotParams( # FRAME
    color='#ff80ed', color_r='#ffa6f2', dpi=400,
    label='LENS2 members: ',
    o_bool=True, o_path=cmn_path, x_label='year', x_lim=[1850, 2100], 
    y_label='GDD base 5', y_lim='auto')

if __name__=='__main__':
    for roi in setp_gdd.reg_oi:
        loop_setp_gdd = cg.SetParams(
            area_stat='mean', mask_flag='none', reg_oi=roi, rlz='all', 
            yrs=[1850, 1859])
        for coi in gp_gdd.cmpst_key:
            loop_gp_gdd = cg.GuideParams(guide_by='mean', cmpst_key=coi)
            ppar_gdd_ts.label = ppar_gdd_ts.label + loop_gp_gdd.cmpst_key
            ic(loop_setp_gdd.reg_oi, loop_gp_gdd.cmpst_key, ppar_gdd_ts.label)
            p = Process(target=puppets.composites, args=(
                dp_gdd, loop_setp_gdd, dp_psl, setp_psl, dp_gdd_roi_alltimes,
                ip, loop_gp_gdd, ppar_comp, ppar_reg, ppar_super, 
                ppar_gdd_ts_gray, ppar_gdd_ts))
            p.start()
            ppar_gdd_ts.label = ppar_gdd_ts.label.replace(
                loop_gp_gdd.cmpst_key, '')