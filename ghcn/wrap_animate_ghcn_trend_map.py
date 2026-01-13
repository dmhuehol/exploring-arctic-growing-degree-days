''' wrap_animate_ghcn_trend_map
Plot maps of stations shaded by trend over a time period and animate
resulting images as mp4.
'''
import sys
sys.path.append(
    "/Users/danielhueholt/Documents/GitHub/"
    + "exploring-arctic-growing-degree-days/")

from icecream import ic
import matplotlib.pyplot as plt
import pandas as pd

import classes_gddt as cg
import fun_ghcn as fg
import fun_plots as fpl

data_path = "/Users/danielhueholt/Data/gddt_data/trends/"
data_fn = "hist_GHCN_abv50N_cov97_span10_trend18732022.csv"
fn_pcs = data_fn.split('_')
msg_dist_bool = True
#  The three PlotParams instances correspond to layers of the plot:
#  ppar: base layer, to plot GHCN stations as dots shaded by trend.
#  ppar_nan: to plot GHCN stations where trends are NaN (e.g., percent
#      change when an endpoint is zero)
#  ppar_out_dist: for denoting stations outside LENS2 distribution
ppar = cg.PlotParams(
    anim_d=dict(frame_tok='*.png', mp4_path=None, mp4_name=None),
    anim_flag=False, ax_facecolor='#cccccc', cb_bool=True, cb_extent='neither',
    cb_label='% change', cb_ticks=None, cb_vals=[-20, 20],
    cmap=fpl.gdd_trend(), dpi=400, figsize=(7, 6), frame_flag=True,
    marker_size=8, o_bool=False, o_name='',
    o_path='/Users/danielhueholt/Documents/Figures/arc-gdd_fig/20260113_comparison/',
    o_prefix='', proj='Arctic', title='', title_size=14)
ppar_nan = cg.PlotParams(
    ax_facecolor='#cccccc', cb_bool=False, cb_vals=[-20, 20], color='#ff7b5a',
    dpi=400, figsize=(7, 6), frame_flag=True, marker='*', marker_size=16,
    o_bool=False, o_prefix='emergeonly_',
    o_path='/Users/danielhueholt/Documents/Figures/arc-gdd_fig/20260113_comparison/',
    proj='Arctic')
ppar_out_dist = cg.PlotParams(
    alpha=0.65, anim_d=dict(frame_tok='*.png', mp4_path=None, mp4_name=None),
    anim_flag=True, ax_facecolor='#cccccc', cb_bool=False, cb_vals=[-20, 20],
    color=None, dpi=400, edgecolors='#000000', figsize=(7, 6), frame_flag=True,
    marker_size=30, o_bool=True, o_prefix='',
    o_path='/Users/danielhueholt/Documents/Figures/arc-gdd_fig/20260113_comparison/',
    proj='Arctic')

df_trends = pd.read_csv(data_path + data_fn)
unique_intervals = sorted([*{*df_trends['interval']}])
for loop_count, interval in enumerate(unique_intervals):
    ppar.o_name = fn_pcs[1] + '_' + fn_pcs[2] + '_' + fn_pcs[3] + '_' \
        + fn_pcs[4] + '_trend' + interval
    ppar.title = 'Trends for GHCN stations above 50N ' + interval + ' cov >97'
    df_interval = df_trends[df_trends['interval'] == interval]
    df_interval_out_dist = df_interval[~df_interval["in_nearest_lens2_dist"]]
    if msg_dist_bool:
        df_interval_indist = df_interval[df_interval["in_nearest_lens2_dist"]]
        msg_in_lens2 = fg.report_in_lens2(df_interval_indist, df_interval)
        ic(msg_in_lens2 + interval)
    df_interval_nan = df_interval[df_interval["nan_trend_bool"]]
    ppar_nan.cmap = fpl.rep_color(color=ppar_nan.color, d=df_interval_nan)

    fig = plt.figure(figsize=ppar.figsize)
    plt.rcParams.update({"font.family": 'Catamaran'})
    plt.rcParams.update({"font.weight": 'light'})
    plt.rcParams.update({"font.size": 10})
    ax = fpl.set_proj(fig, ppar)
    #  Base layer: all lat, lon, trends
    ppar.o_name = ppar.o_prefix + ppar.o_name
    ppar.title = ppar.title +  ' n=' + str(len(df_interval.index))
    fpl.plot_globe_ng(
        df_interval["slope"], df_interval["lat"], df_interval["lon"], ppar,
        ax=ax)
    #  Second layer: NaN lat, lon, trends (usually no data)
    ppar_nan.o_name = ppar_nan.o_prefix + ppar.o_name
    ax = fpl.plot_globe_ng(
        df_interval_nan["slope"], df_interval_nan["lat"],
        df_interval_nan["lon"], ppar_nan, ax=ax)
    #  Third layer: in/outside model distribution
    ppar_out_dist.o_name = ppar_out_dist.o_prefix + ppar.o_name
    ppar_out_dist.title = ppar.title
    ax = fpl.plot_globe_ng(
        df_interval_out_dist["slope"], df_interval_out_dist["lat"],
        df_interval_out_dist["lon"], ppar_out_dist, ax=ax)
    plt.close()

por_str = fn_pcs[4].replace('span', '')
if ppar.anim_flag:
    fpl.movie_maker(ppar, unique_intervals)
if ppar_out_dist.anim_flag:
    fpl.movie_maker(ppar_out_dist, unique_intervals)
