'''wrap_bmb_crossover_hist
Plot difference in crossover distributions between biomass burning
ensembles.

This requires pre-calculated gridded files with subsamples of crossover
data with dimensions time, lat, lon (and realizations optional) for both
the CMIP6 forcings and the smoothed biomass burning forcings. To make
these files, run wrap_calc_exceedance_grid and wrap_calc_crossover_grid.
'''
import sys

from icecream import ic
from matplotlib import colormaps as cm
import matplotlib.pyplot as plt
import numpy as np
from pyfonts import load_google_font, set_default_font
import seaborn as sn

import classes_gddt as cg
import fun_calc_var as fcv
import fun_plots as fpl
import fun_process as fproc

cmn_path = '/Users/danielhueholt/Documents/Figures/arc-gdd_fig/20260116_cont/'
l_cmip6 = list()
l_smoothed = list()
l_compare_in_cmip6 = list()
l_compare_in_smoothed = list()
l_compare_between = list()
n_samples = 25
ppar = cg.PlotParams(
    bw=2, color_comp='#6f8c31', color_r='#4e318c', dpi='.pdf', figsize=(7, 7),
    font='Catamaran',
    leg_dict={"bool": True, "frameon": False, "loc": 'upper left', "size": 18},
    o_name='',
    o_path='/Users/danielhueholt/Documents/Figures/arc-gdd_fig/20260116_cont/',
    o_prefix='', title='Histogram of all trends', x_lim=[-40, 30])
use_font = load_google_font(ppar.font, danger_not_verify_ssl=True)
set_default_font(use_font)

for sample in np.arange(1, n_samples + 1):
    dp_crossover_cmip6 = cg.DataParams(
        path='/Users/danielhueholt/Data/gddt_data/LENS2/exceedance/crossover/',
        tok='rs' + str(sample) + '-cmip6_forcedcrossover_threshold80.nc',
        var='crossover',
        flag_raw_ds=True, flag_raw_da=True, flag_time_slice=False,
        flag_manage_rlz=True, flag_land_mask=True, flag_roi=True)
    dp_crossover_smoothed = cg.DataParams(
        path='/Users/danielhueholt/Data/gddt_data/LENS2/exceedance/crossover/',
        tok='rs' + str(sample) + '-smoothed_forcedcrossover_threshold80.nc',
        var='crossover',
        flag_raw_ds=True, flag_raw_da=True, flag_time_slice=False,
        flag_manage_rlz=True, flag_land_mask=True, flag_roi=True)
    setp_crossover = cg.SetParams(
        area_stat='pass', base_yrs=[0, 2000],
        mask='/Users/danielhueholt/Data/gddt_data/mask/cesm_atm_mask.nc',
        mask_flag='land', reg_oi='global', rlz='all', yrs=[1850, 2100])
    crossover_dict_cmip6 = fproc.common_opener(dp=dp_crossover_cmip6, setp=setp_crossover)
    crossover_dict_smoothed = fproc.common_opener(dp=dp_crossover_smoothed, setp=setp_crossover)
    l_cmip6.append(crossover_dict_cmip6['land_mask'].compute())
    l_smoothed.append(crossover_dict_smoothed['land_mask'].compute())
for sample in range(n_samples):
    comparison_sample = (sample + 1) % n_samples
    l_compare_in_cmip6.append(
        np.ravel(l_cmip6[sample] - l_cmip6[comparison_sample]))
    l_compare_in_smoothed.append(
        np.ravel(l_smoothed[sample] - l_smoothed[comparison_sample]))
    l_compare_between.append(np.ravel(l_cmip6[sample] - l_smoothed[sample]))
ppar.o_name = 'difference-hist.pdf'
ppar.title = 'Difference between subsamples'

plot_cmip6 = np.ravel(l_compare_in_cmip6)
plot_smoothed = np.ravel(l_compare_in_smoothed)
plot_between = np.ravel(l_compare_between)
ic(np.shape(l_compare_in_cmip6), np.shape(plot_cmip6))
fig, ax = plt.subplots(figsize=ppar.figsize)
plt.rcParams.update({'font.family': 'Catamaran'})
#  'font.weight': normal, bold, heavy, light, ultrabold, ultralight
plt.rcParams.update({'font.weight': 'normal'})
plt.rcParams.update({'font.size': 12})
ax.spines[['right', 'top']].set_visible(False)
sn.histplot(
    plot_cmip6, binwidth=ppar.bw, label='within CMIP6', color='#f07024',
    alpha=0.5, stat='probability', element='step')
sn.histplot(
    plot_smoothed, binwidth=ppar.bw, label='within smoothed', color='#1c542c',
    alpha=0.5, stat='probability', element='step')
sn.histplot(
    plot_between, binwidth=ppar.bw, label='CMIP6 - smoothed', color='#99d6e5',
    alpha=0.5, stat='probability', element='step')
if ppar.x_lim != 'auto':
    plt.xlim(ppar.x_lim)
if ppar.leg_dict["bool"]:
    plt.legend(
        fontsize=ppar.leg_dict["size"], frameon=ppar.leg_dict["frameon"],
        loc=ppar.leg_dict["loc"])
out_name_full = ppar.o_path + ppar.o_prefix + ppar.o_name
ic(out_name_full)
plt.savefig(out_name_full)
