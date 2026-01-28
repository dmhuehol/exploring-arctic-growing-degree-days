'''wrap_crossover_ts
Plot crossover for a point or region as a timeseries:
    1) Forced crossover
        Defined when the ensemble mean exceeds a certain threshold of
        the preindustrial base period; by default 80% of samples,
        although others can be used.
    2) Member crossover
        Defined for each ensemble member when it exceeds a certain
        threshold of the preindustrial base period; by default 90% of
        samples. A threshold of 100% exceedance defines a no-analog
        state.

The calculation itself addresses the question: "How many samples in one
period (np_compare_to_base) exceed the samples from a baseline
(np_base_samples)?" As configured here, this compares samples from the
CESM2 Large Ensemble to the CESM2 Preindustrial control. Beyond a
certain threshold, this denotes "crossover" from a climate state of
internal climate variability to one where the warming trend dominantes.

The underlying statistic is a non-parametric expression of effect size.
The choice of threshold in the definition involves subjectivity. See the
accompanying manuscript Hueholt et al. "Exploring the Influence of
Internal Climate Variability and Forced Change on Arctic Greening" for
discussion of the values used by default here, and
wrap_crossover_sensitivity for a sensitivity analysis.
'''
import sys

from icecream import ic
import matplotlib.pyplot as plt
import numpy as np

import classes_gddt as cg
import fun_plots as fpl
import fun_process as fproc
import gddt_region_library as g_rlib

roll_windows = cg.IntervalParams().create_intvls(
    strt_yr=1850, end_yr=2101, spn=10, type='rolling')
dp_gdd = cg.DataParams(
    path='/Users/danielhueholt/Data/gddt_data/gdd/lens2/',
    tok='*arc*.nc', var='gdd5_sum', flag_raw_ds=True,
    flag_raw_da=True, flag_time_slice=False, flag_manage_rlz=True,
    flag_land_mask=False, flag_roi=True)
dp_gdd_base = cg.DataParams(
    path='/Users/danielhueholt/Data/gddt_data/gdd/preindustrial/',
    tok='*.nc', var='gdd5_sum', flag_raw_ds=False, flag_raw_da=True,
    flag_time_slice=True, flag_manage_rlz=True,
    flag_land_mask=True, flag_roi=True)
#  setp_gdd.base_yrs determines period of samples to compare to
#  setp_gdd.window is the window applied to the plots
setp_gdd = cg.SetParams(
    area_stat='mean', base_yrs=[0, 2000], mask_flag='none',
    reg_oi=g_rlib.BrooksRange_colonist(), rlz='all', window=10,
    yrs=[1850, 2100])
setp_gdd_base = cg.SetParams(
    area_stat='mean', base_yrs=[0, 2000], mask_flag='none',
    reg_oi=setp_gdd.reg_oi, rlz='all', window=None, yrs=[0, 2000])
#  PlotParams settings for FORCED CROSSOVER and RANGE OF MEMBER
#  CROSSOVERS with STORYLINE MEMBERS 24 AND 52 HIGHLIGHTED, (Figure 2a):
ppar_forced_range_storyline2452 = cg.PlotParams(
    dpi=400, leg_bool=False, lw=0.5,
    mn_dict=dict(
        bool=True, color='#515151', label='Ensemble mean', linestyle='-',
        lw=3),
    o_bool=True, o_name='crossover' + '_' + setp_gdd.reg_oi['reg_abv'] \
        + '_base' + str(setp_gdd.base_yrs[0]) + '-' \
        + str(setp_gdd.base_yrs[1]),
    o_path='/Users/danielhueholt/Documents/Figures/arc-gdd_fig/20260120_ts/',
    o_prefix='', plot_as_percent=True,
    plot_crossover_dict=dict(
        forced_dict=dict(
                bool=True, exceed_color=('#515151',),
                exceed_linestyle=('--',),
                exceed_lw=(0.4, 0.4, 0.4, 1.2, 0.4, 0.4),
                threshold=(50, 60, 70, 80, 90, 100),
                years_alpha=(0, 0, 0, 1, 0, 0),
                years_color=('#515151',), years_linestyle=('--',),
                years_lw=(1.3,)),
        member_dict=dict(
                bool=True, exceed_alpha=(0,), exceed_color=('#97215c',),
                exceed_linestyle=('--',), exceed_lw=(0.3,),
                range_dict=dict(
                    alpha=0.6, bool=True, color='#d285ae', edgecolor=None,
                    range=[0.1, 0.9]),
                threshold=(90,), years_alpha=(0,), years_color=('#97215c',),
                years_linestyle=('--',), years_lw=(0.3,))),
    rlz_dict=dict(
        bool=True, color='#cccccc', label='Members', linestyle='-', lw=0.5),
    storyline_dict=dict(
        select=[24, 52], color=['#b66363', '#2098ae'],
        label=['member 24', 'member 52'],
        linestyle=['-', '-'], lw=[1.5, 1.5]),
    title=setp_gdd.reg_oi['reg_str'],
    ts_type='spaghetti', x_label='', x_lim=[1850, 2100],
    y_label='exceed percent of Preindustrial control years', y_lim=[0, 100.05])
#  PlotParams settings for MEMBER CROSSOVER and ONLY STORYLINE 24
#  (Figure 2b):
ppar_storyline24 = cg.PlotParams(
    dpi=400, leg_bool=False, lw=0.5,
    mn_dict=dict(
        bool=False, color='#515151', label='Ensemble mean', linestyle='-',
        lw=3),
    o_bool=True, o_name='crossover' + '_' + setp_gdd.reg_oi['reg_abv'] \
        + '_base' + str(setp_gdd.base_yrs[0]) + '-' \
        + str(setp_gdd.base_yrs[1]),
    o_path='/Users/danielhueholt/Documents/Figures/arc-gdd_fig/20260120_ts/',
    o_prefix='', plot_as_percent=True,
    plot_crossover_dict=dict(
        forced_dict=dict(bool=False,),
        member_dict=dict(
                bool=True, exceed_alpha=(1,), exceed_color=('#b66363',),
                exceed_linestyle=('--',), exceed_lw=(1.2,),
                range_dict=dict(bool=False), threshold=(90,), years_alpha=(1,),
                years_color=('#b66363',), years_linestyle=('--',),
                years_lw=(1.2,))),
    rlz_dict=dict(bool=False,),
    storyline_dict=dict(
        select=[24,], color=['#b66363',], label=['member 24',],
        linestyle=['-',], lw=[2,]),
    title=setp_gdd.reg_oi['reg_str'], ts_type='spaghetti', x_label='',
    x_lim=[1850, 2100],
    y_label='exceed percent of Preindustrial control years', y_lim=[0, 100.05],
    yticks=[30, 60, 90])
#  PlotParams settings for MEMBER CROSSOVER and ONLY STORYLINE 52
#  (Figure 2c):
ppar_storyline52 = cg.PlotParams(
    dpi=400, leg_bool=False, lw=0.5,
    mn_dict=dict(
        bool=False, color='#515151', label='Ensemble mean', linestyle='-',
        lw=3),
    o_bool=True, o_name='crossover' + '_' + setp_gdd.reg_oi['reg_abv'] \
        + '_base' + str(setp_gdd.base_yrs[0]) + '-' \
        + str(setp_gdd.base_yrs[1]),
    o_path='/Users/danielhueholt/Documents/Figures/arc-gdd_fig/20260120_ts/',
    o_prefix='', plot_as_percent=True,
    plot_crossover_dict=dict(
        forced_dict=dict(bool=False,),
        member_dict=dict(
                bool=True, exceed_alpha=(1,), exceed_color=('#2098ae',),
                exceed_linestyle=('--',), exceed_lw=(1.2,),
                range_dict=dict(bool=False), threshold=(90,), years_alpha=(1,),
                years_color=('#2098ae',), years_linestyle=('--',),
                years_lw=(1.2,))),
    rlz_dict=dict(bool=False,),
    storyline_dict=dict(
        select=[52,], color=['#2098ae',], label=['member 52',],
        linestyle=['-',], lw=[2,]),
    title=setp_gdd.reg_oi['reg_str'], ts_type='spaghetti', x_label='',
    x_lim=[1850, 2100],
    y_label='exceed percent of Preindustrial control years', y_lim=[0, 100.05],
    yticks=[30, 60, 90])
#  PlotParams settings for MODIFICATION FOR EXPLORATORY ANALYSIS
#  currently exploring meaning and interpretation of no-analog states
ppar_exploratory = cg.PlotParams(
    dpi=400, leg_bool=False, lw=0.5,
    mn_dict=dict(
        bool=False, color='#515151', label='Ensemble mean', linestyle='-',
        lw=3),
    o_bool=True, o_name='crossover' + '_' + setp_gdd.reg_oi['reg_abv'] \
        + '_base' + str(setp_gdd.base_yrs[0]) + '-' \
        + str(setp_gdd.base_yrs[1]),
    o_path='/Users/danielhueholt/Documents/Figures/arc-gdd_fig/20260120_ts/',
    o_prefix='exploratory_', plot_as_percent=True,
    plot_crossover_dict=dict(
        forced_dict=dict(
                bool=False, exceed_color=('#515151',),
                exceed_linestyle=('--',),
                exceed_lw=(0.4, 0.4, 0.4, 1.2, 0.4, 0.4),
                threshold=(50, 60, 70, 80, 90, 100),
                years_alpha=(0, 0, 0, 1, 0, 0),
                years_color=('#515151',), years_linestyle=('--',),
                years_lw=(1.3,)),
        member_dict=dict(
                bool=True, exceed_alpha=(1,), exceed_color=('#515151',),
                exceed_linestyle=('--',), exceed_lw=(1.2,),
                range_dict=dict(
                    alpha=0.6, bool=True, color='#d285ae', edgecolor=None,
                    range=[0.1, 0.9]),
                threshold=(100,), years_alpha=(0,), years_color=('#97215c',),
                years_linestyle=('--',), years_lw=(0.3,))),
    rlz_dict=dict(
        bool=True, color='#cccccc', label='Members', linestyle='-', lw=0.5),
    # storyline_dict=dict(
    #     select=[24,52], color=['#b66363', '#2098ae'],
    #     label=['member 24', 'member 52'],
    #     linestyle=['-', '-'], lw=[1.5, 1.5]),
    storyline_dict=dict(
        select=[13,], color=['#515151',],
        label=['member 13',],
        linestyle=['-',], lw=[0.5,]),
    title=setp_gdd.reg_oi['reg_str'],
    ts_type='spaghetti', x_label='', x_lim=[1850, 2100],
    y_label='exceed percent of Preindustrial control years', y_lim=[0, 100.05])
#  Choose PlotParams to use here
ppar_to_use = ppar_exploratory

dict_data = fproc.common_opener(dp=dp_gdd, setp=setp_gdd)
da_data_roi = dict_data['roi'].compute()
dict_base = fproc.common_opener(dp=dp_gdd_base, setp=setp_gdd_base)
da_base_roi = dict_base['roi'].compute()
da_base_period = da_base_roi.sel(
    year=slice(setp_gdd.base_yrs[0], setp_gdd.base_yrs[1]))
years_to_plot = np.arange(setp_gdd.yrs[0], setp_gdd.yrs[1] + 1)
np_exceed, l_exceed_rolling_avg = fproc.common_calc_exceed(
    da_data_roi, da_base_period, years_to_plot, ppar_to_use, setp_gdd)

plt.figure()
plt.rcParams.update({'font.family': 'Catamaran'})
plt.rcParams.update({'font.weight': 'light'}) #normal, bold, heavy, light, ultrabold, ultralight
plt.rcParams.update({'font.size': 12})
if setp_gdd.window is None:
    ppar_to_use.o_name = 'nowindow_' + ppar_to_use.o_name
    fpl.plot_exceed_crossover_ts(
        np_exceed, x_d=years_to_plot, ppar=ppar_to_use)
else:
    plot_this = np.array(l_exceed_rolling_avg).T
    years_windows_plot = years_to_plot[:-(setp_gdd.window - 1)]
    ppar_to_use.o_name = str(setp_gdd.window) + 'yrwindow_' \
        + ppar_to_use.o_name
    fpl.plot_exceed_crossover_ts(
        plot_this, x_d=years_windows_plot, ppar=ppar_to_use)
