''' wrap_effects
Plot timeseries of effect size over time for region_library object. Supports
the following methods as inputs to setp_gdd.effect:
    'cliffs': Cliffs delta statistic
    'cliffs_mean': Cliffs delta statistic of ensemble mean
    'gexc': Samples exceeded by each realization (Gexc in robustness)
    'robustness': Robustness from Hueholt et al. 2022

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''
import sys

from icecream import ic
import matplotlib.pyplot as plt
import numpy as np

import classes_gddt as cg
import fun_plots as fpl
import fun_process as fproc
import gddt_region_library as g_rlib

plot_hist_gexc = False
#  Manual animation of existing exceedance histograms
# fpl.images_mp4(
#     '/Users/dhueholt/Documents/gddt_fig/20241028_threshAndStory/', '*hist*', 
#     '/Users/dhueholt/Documents/gddt_fig/20241028_threshAndStory/', 
#     'exc_hist_brooksrangepoint.mp4', fps=6)
roll_wind = cg.IntervalParams().create_intvls(strt_yr=1850, end_yr=2101, spn=10, type='rolling')
dp_gdd = cg.DataParams(
    path='/Users/dhueholt/Documents/gddt_data/gdd/lens2/', 
    tok='*arc*.nc', var='gdd5_sum', flag_raw_ds=True, 
    flag_raw_da=True, flag_time_slice=False, flag_manage_rlz=True, 
    flag_land_mask=False, flag_roi=True)
dp_gdd_base = cg.DataParams(
    path='/Users/dhueholt/Documents/gddt_data/gdd/preindustrial/', 
    tok='*.nc', var='gdd5_sum', flag_raw_ds=False, flag_raw_da=True, 
    flag_time_slice=True, flag_manage_rlz=True, 
    flag_land_mask=True, flag_roi=True)
setp_gdd = cg.SetParams(
    area_stat='mean', base_yrs=[0, 2000], beat=4000, mask_flag='none', 
    effect='gexc',
    reg_oi=g_rlib.BrooksRange_colonist(), rlz='all', 
    window=10, yrs=[1850, 2100])
setp_gdd_base = cg.SetParams(
    area_stat='mean', base_yrs=[0, 2000], beat=4000, mask_flag='none', 
    effect='gexc',
    reg_oi=setp_gdd.reg_oi, rlz='all', 
    window=None, yrs=[0, 2000])
ppar = cg.PlotParams(
    # color='#a50b5e', color_r='#a50b5e', dpi=400, 
    color='#2098ae', color_r='#97215C', dpi=400, 
    forced_crossover_bool=False, label='LENS2', lw=0.5, #0.5, 3
    member_crossover_bool=True,
    mn_bool=False, o_bool=True, plot_all=False,
    o_path='/Users/dhueholt/Documents/gddt_fig/20250220_noanstory/',
    o_name=setp_gdd.effect + '_' + setp_gdd.reg_oi['reg_abv'] \
        + '_base' + str(setp_gdd.base_yrs),
    o_prefix='',
    storyline=None, title=setp_gdd.reg_oi['reg_str'],
    ts_type = 'spaghetti', x_label='', x_lim=[1850, 2100],
    # y_label = '', yticks=[],
    y_label='exceed percent of Preindustrial control years', #'members >80% (4,000) of base period years', 
    y_lim=[0, 100.05])
ppar_hist = cg.PlotParams(
    bw=5, o_path=ppar.o_path, 
    o_name='gexc_hist' + '_' + setp_gdd.reg_oi['reg_abv'], stat='count', 
    x_label='percent of preindustrial control years', x_lim=[-1, 101], 
    y_label='samples')

d_data = fproc.common_opener(dp=dp_gdd, setp=setp_gdd)
da_roi = d_data['roi'].compute()
d_base = fproc.common_opener(dp=dp_gdd_base, setp=setp_gdd_base)
da_base_roi = d_base['roi'].compute()
base_period = da_base_roi.sel(year=slice(setp_gdd.base_yrs[0], setp_gdd.base_yrs[1]))
np_x1 = np.ravel(base_period.data)
list_effect = list()
years = np.arange(setp_gdd.yrs[0], setp_gdd.yrs[1] + 1)
#  Stores Gexc values (number of base samples exceeded by member at a year)
list_gexc = list()
list_gm = list()
for yr in years:
    da_act_yr = da_roi.sel(year=yr)
    np_x2 = da_act_yr.data
    if 'cliffs' in setp_gdd.effect:
        if setp_gdd.effect == 'cliffs_mean':
            np_x2 = da_act_yr.mean(dim='realization').data
        yr_effect = fproc.cliffs_delta(np_x1, np_x2)
        list_effect.append(yr_effect)
    else:
        d_gexc = fproc.gexc(np_x1, np_x2)
        dgm = fproc.gexc(np_x1, [np.mean(np_x2),])
        np_gexc = np.array(d_gexc['above'])
        list_gexc.append(np_gexc)
        list_gm.append(dgm['above'])
        count_robust = fproc.threshold_count(np_gexc, beat=setp_gdd.beat)
        list_effect.append(count_robust)
if 'cliffs' in setp_gdd.effect:
    ppar.lw = 3
    ppar.y_lim = [-1, 1]
    if setp_gdd.effect == 'cliffs_mean':
        ppar.y_label = 'cliffs delta of ensemble mean'
    else:
        ppar.y_label = 'cliffs delta'
    np_effect = np.array(list_effect)
elif setp_gdd.effect == 'gexc':
    ppar.y_lim = [0, 100.05]
    np_effect = np.array(list_gexc) / len(np_x1) * 100
    np_all_gexc = np_effect
elif setp_gdd.effect == 'robustness':
    ppar.y_lim = [0, 100.05]
    ppar.lw = 3
    ppar.o_name = 'beat' + str(setp_gdd.beat) + '_' + ppar.o_name
    np_effect = np.array(list_effect)
    np_all_gexc = np.ravel(np.array(list_gexc))

plt.figure()
plt.rcParams.update({'font.family': 'Catamaran'})
plt.rcParams.update({'font.weight': 'light'}) #normal, bold, heavy, light, ultrabold, ultralight
plt.rcParams.update({'font.size': 12})
if setp_gdd.window is None:
    ppar.o_name = 'nowindow_' + ppar.o_name
    fpl.plot_effect_timeseries(np_effect, x_d=years, ppar=ppar)
else:
    l_effect_rolling_avg = list()
    lera_gm = fproc.moving_average(np.squeeze(np.array(list_gm)).T, setp_gdd.window)
    for rlz in np.arange(0, 100):
        rolling_avg = fproc.moving_average(np_effect[:, rlz], setp_gdd.window)
        l_effect_rolling_avg.append(rolling_avg)
    plot_this = np.array(l_effect_rolling_avg)
    ppar.o_name = str(setp_gdd.window) + 'yrwindow_' + ppar.o_name
    fpl.plot_effect_timeseries(plot_this.T, x_d=years[:-(setp_gdd.window - 1)], lera_gm=np.array(lera_gm).T, ppar=ppar)
if plot_hist_gexc:
    yrs_ky = str(years[0]) + '-' + str(years[-1])
    d_hist = {yrs_ky: np_all_gexc}
    save_o_name = ppar_hist.o_name
    ppar_hist.o_name = yrs_ky + ppar_hist.o_name
    ppar_hist.title = setp_gdd.reg_oi['reg_str'] + ' exceedance hist ' + yrs_ky
    fpl.plot_hist(d_hist, ppar=ppar_hist)
    ppar_hist.o_name = save_o_name