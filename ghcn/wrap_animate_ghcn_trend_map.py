''' wrap_animate_ghcn_trend_map
Plot maps of stations shaded by trend over a time period and animate 
resulting images as mp4.

To make an individual still image: 
    1) Obtain a single interval using scp.por and IntervalParams
    2) Set ppar.frame_flag=True and ppar.anim_flag=False

To animate existing images:
    1) Set ppar.frame_flag=False
'''
import collections
from itertools import compress
import sys
sys.path.append(
    "/Users/dhueholt/Documents/Github/" 
    + "exploring-arctic-growing-degree-days/")
import time

from icecream import ic
import matplotlib.pyplot as plt
import polars as pl

import classes_gddt as cg
import fun_calc_var as fcv
import fun_ghcn as fg
import fun_plots as fpl
import fun_process as fproc
import gddt_region_library as g_rlib

d_path = "/Users/dhueholt/Documents/gddt_data/Arctic50N_gt1/"
sc_f = "~/Documents/ghcn_data/spancov_lists/spancov_Arctic50N_gt1.csv"
scp = cg.SpanCovParams(f=sc_f, por=30, cov_thr=97, cov_type='>')
ip = cg.IntervalParams(
    span=scp.por, strt_yr=1873, end_yr=2022, type='noverlap')
transform = 'percent_change'

dp_lens2 = cg.DataParams(
    path='/Users/dhueholt/Documents/gddt_data/gdd/lens2/', tok='*arc*.nc',
    var='gdd5_sum', flag_raw_ds=True, flag_raw_da=True, flag_time_slice=True,
    flag_manage_rlz=True, flag_land_mask=False, flag_roi=True)
setp_lens2 = cg.SetParams(
    area_stat='pass', reg_oi=g_rlib.Arctic50N(), rlz='all')
animate_only = True
check_dist_bool = True
#  The three PlotParams instances correspond to layers of the plot:
#  ppar: base layer, to plot GHCN stations as dots shaded by trend.
#  ppar_nan_trends: to plot GHCN stations where trends are NaN (e.g., 
#      when percent change is plotted but one of the endpoints is zero)
#  ppar_outside_dist: for denoting stations outside LENS2 distribution
ppar = cg.PlotParams(
    anim_d=dict(frame_tok='total*.png', mp4_path=None, mp4_name=None), 
    anim_flag=False, ax_facecolor='#cccccc', cb_bool=True, cb_extent='neither',
    cb_label='% change', cb_ticks=None, cb_vals=[-20, 20], 
    cmap=fpl.gdd_trend(), dpi=400, figsize=(7, 6), frame_flag=False, 
    marker_size=8, o_bool=False, o_name=None, 
    o_path='/Users/dhueholt/Documents/gddt_fig/20250220_refactoringAndNoAnStory/',
    o_prefix='', proj='Arctic', title=None, title_size=14)
ppar_nan_trends = cg.PlotParams(
    ax_facecolor='#cccccc', cb_bool=False, cb_vals=[-20, 20], color='#ff7b5a',
    dpi=400, figsize=(7, 6), frame_flag=False, marker='*', marker_size=16, 
    o_bool=False, o_prefix='emergeonly_', 
    o_path='/Users/dhueholt/Documents/gddt_fig/20250220_refactoringAndNoAnStory/',
    proj='Arctic')
ppar_outside_dist = cg.PlotParams(
    alpha=0.65, anim_d=dict(frame_tok='*.png', mp4_path=None, mp4_name=None), 
    anim_flag=True, ax_facecolor='#cccccc', cb_bool=False, cb_vals=[-20, 20], 
    color=None, dpi=400, edgecolors='#000000', figsize=(7, 6), frame_flag=True, 
    marker_size=30, o_bool=True, o_prefix='', 
    o_path='/Users/dhueholt/Documents/gddt_fig/20250220_refactoringAndNoAnStory/', 
    proj='Arctic')

nan_trend_count = 0
for loop_count, interval in enumerate(ip.intervals):
    #  Need basic name information even if frame_flag is False!
    yr_str = str(interval[0]) + '-' + str(interval[1])
    if check_dist_bool:
        ppar_outside_dist.o_name = 'GHCN_abv50N_cov' + scp.cov_thr_str \
            + '_span' + scp.por_str + '_trend' + yr_str
        ppar_outside_dist.title = 'Trends for GHCN stations above 50N ' \
            + yr_str + ' cov >97 '
    else:
        ppar.o_name = 'GHCN_abv50N_cov' + scp.cov_thr_str + '_span' + scp.por_str \
            + '_trend' + yr_str
        ppar.title = 'Trends for GHCN stations above 50N ' + yr_str + ' cov >97 '
    if ppar.frame_flag is False:
        break
    if check_dist_bool:
        setp_lens2.yrs = interval
        lens2_d = fproc.common_opener(dp=dp_lens2, setp=setp_lens2)
        da_lens2 = lens2_d['roi'].compute()
        years = da_lens2.year
        d_trend = fcv.calc_lin_reg_vec(years, da_lens2)
        da_lens2_grad = d_trend['grad']
    ic(cg.TrackProg(cli=loop_count, iter=ip.intervals).message())
    df_sc = fg.get_spancov(scp)
    t_slc = (pl.datetime(interval[0], 1, 1), pl.datetime(interval[1], 12, 31))
    lf_trend = collections.defaultdict(list)
    lf_nan_trend = collections.defaultdict(list)
    stid = df_sc["st_id"]
    tic = time.time()
    # Calculate trend for each station
    for cc, act_st in enumerate(stid):
        act_name = act_st + '_refined.csv'
        ldf_act = pl.scan_csv(d_path + act_name)
        d_lr = fg.trend_ghcn(ldf_act, scp, t_slc)
        grad = d_lr['grad']
        if grad is None:
            continue
        else:
            ll = ldf_act.select(pl.first(['LATITUDE', 'LONGITUDE'])).collect()
            ll_lat = ll.select('LATITUDE').item()
            ll_lon = ll.select('LONGITUDE').item()
            if check_dist_bool:
                ll_lon360 = ll_lon % 360
                np_lens2_grad_ll = da_lens2_grad.sel(
                    lat=ll_lat, lon=ll_lon360, method='nearest').data
                check_in_lens2_dist = fproc.check_in_dist(
                    np_lens2_grad_ll, grad)
            else:
                check_in_lens2_dist = None
            if d_lr['data'][0] != 0:
                lr_trnd = d_lr['data'][-1] - d_lr['data'][0]
                per_chng =  lr_trnd / d_lr['data'][0] * 100
                lf_trend["grad"].append(grad)
                lf_trend["per_chng"].append(per_chng)
                lf_trend["lat"].append(ll_lat)
                lf_trend["lon"].append(ll_lon)
                lf_trend["in_dist"].append(check_in_lens2_dist)
            else:
                nan_trend_count = nan_trend_count + 1
                lf_nan_trend["grad"].append(100)
                lf_nan_trend["per_chng"].append(100)
                lf_nan_trend["lat"].append(ll_lat)
                lf_nan_trend["lon"].append(ll_lon)
                lf_nan_trend["in_dist"].append(check_in_lens2_dist)
    toc = time.time() - tic; ic(toc)
    #  Allow any additional manual transforms
    if transform == "percent_change":
        plot_grad = lf_trend["per_chng"]
    else:
        plot_grad = lf_trend["grad"]
    plot_lat = lf_trend["lat"]
    plot_lon = lf_trend["lon"]
    ppar_nan_trends.cmap = fpl.rep_color(
        color=ppar_nan_trends.color, d=lf_nan_trend["grad"])
    
    fig = plt.figure(figsize=ppar.figsize)
    plt.rcParams.update({"font.family": 'Catamaran'})
    plt.rcParams.update({"font.weight": 'light'})
    plt.rcParams.update({"font.size": 10})
    ax = fpl.set_proj(fig, ppar)
    ppar.o_name = ppar.o_prefix + ppar.o_name
    ppar.title = ppar.title +  ' n=' + str(len(plot_grad))
    fpl.plot_globe_ng(plot_grad, plot_lat, plot_lon, ppar, ax=ax)
    ppar_nan_trends.o_name = ppar_nan_trends.o_prefix + ppar.o_name
    ax = fpl.plot_globe_ng(
        lf_nan_trend["per_chng"], lf_nan_trend["lat"], lf_nan_trend["lon"], 
        ppar_nan_trends, ax=ax)
    if check_dist_bool:
        outside_dist = [not elem for elem in lf_trend["in_dist"]]
        grad_outside_dist = list(compress(plot_grad, outside_dist))
        lat_outside_dist = list(compress(plot_lat, outside_dist))
        lon_outside_dist = list(compress(plot_lon, outside_dist))
        ppar_outside_dist.o_name = ppar_outside_dist.o_prefix + ppar.o_name
        ppar_outside_dist.title = ppar.title
        fpl.plot_globe_ng(
            grad_outside_dist, lat_outside_dist, lon_outside_dist, 
            ppar_outside_dist, ax=ax)
    plt.close()
nan_trend_msg = 'Number of NaN trends: ' + str(nan_trend_count)
ic(nan_trend_msg)

if ppar.anim_flag:
    fpl.movie_maker(ppar, ip.intervals, scp, yr_str)
if ppar_outside_dist.anim_flag:
    fpl.movie_maker(ppar_outside_dist, ip.intervals, scp, yr_str)