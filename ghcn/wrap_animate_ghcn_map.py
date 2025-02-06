''' wrap_animate_ghcn_map
Plot maps of stations shaded by trend over a time period, and animate
the resulting images as an mp4.

To make an individual still image, select a single interval using 
scp.por and IntervalParams, and set ppar.frame_flag=True and 
anim_flag=False.


########################################################################
###############################################################################
Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''
import collections
from itertools import compress
import math
import sys
sys.path.append(
    "/Users/dhueholt/Documents/Github/" 
    + "growing-degree-days-treelines-poles/")
import time

from icecream import ic
import matplotlib.pyplot as plt
import polars as pl

import classes_gddt as cg
import fun_calc_var as fcv
import fun_ghcn as fg
import fun_plots as fpl
import fun_process as fproc
import gddt_region_library as grlib

d_path = "/Users/dhueholt/Documents/gddt_data/Arctic50N_gt1/"
sc_f = "~/Documents/ghcn_data/spancov_lists/spancov_Arctic50N_gt1.csv"
scp = cg.SpanCovParams(f=sc_f, por=5, cov_thr=97, cov_type='>')
ip = cg.IntervalParams(span=scp.por, strt_yr=1873, end_yr=2022, type='noverlap')
calc = 'percent_change'
in_dist_bool = True
ppar = cg.PlotParams(
    anim_d=dict(frame_tok='total*.png', mp4_path=None, mp4_name=None), 
    anim_flag=False, ax_facecolor='#cccccc', cb_bool=True, cb_extent='neither',
    cb_label='% change', cb_ticks=None, cb_vals=[-20, 20], cmap=fpl.gdd_trend(),
    dpi=400, figsize=(7, 6), frame_flag=True, marker_size=8, o_bool=False, o_name=None, 
    o_path='/Users/dhueholt/Documents/gddt_fig/20240725_inL2DistOptAndComps/',
    o_prefix='', proj='Arctic', title=None, title_size=14)
ppar_nan = cg.PlotParams(
    ax_facecolor='#cccccc', cb_bool=False, cb_vals=[-20, 20], color='#ff7b5a',
    dpi=400, figsize=(7, 6), frame_flag=False, marker='*', marker_size=16, 
    o_bool=False, o_prefix='emergeonly_', o_path='/Users/dhueholt/Documents/gddt_fig/20240725_inL2DistOptAndComps/',
    proj='Arctic')
dp_l2 = cg.DataParams(
    path='/Users/dhueholt/Documents/gddt_data/gdd/nomask/', tok='*arc*.nc',
    var='gdd5_sum', flag_raw_ds=True, flag_raw_da=True, flag_time_slice=True,
    flag_manage_rlz=True, flag_land_mask=False, flag_roi=True)
setp_l2 = cg.SetParams(
    area_stat='pass', reg_oi=grlib.Arctic50N(), rlz='all')
ppar_in_dist = cg.PlotParams(
    alpha=0.65, anim_d=dict(frame_tok='*.png', mp4_path=None, mp4_name=None), 
    anim_flag=True, ax_facecolor='#cccccc', cb_bool=False, cb_vals=[-20, 20], color=None, dpi=400, 
    edgecolors='#000000', figsize=(7, 6), frame_flag=True, marker_size=30, 
    o_bool=True, o_prefix='', 
    o_path='/Users/dhueholt/Documents/gddt_fig/20240725_inL2DistOptAndComps/', proj='Arctic')

zero_count = 0
for atyc, aty in enumerate(ip.intervals):
    #  Need basic name information even if frame_flag is False!
    yr_str = str(aty[0]) + '-' + str(aty[1])
    ppar.o_name = 'GHCN_abv50N_cov' + scp.cov_thr_str + '_span' + scp.por_str \
        + '_trend' + yr_str
    ppar.title = 'Trends for GHCN stations above 50N ' + yr_str + ' cov >97 '
    if ppar.frame_flag is False:
        break
    if in_dist_bool:
        setp_l2.yrs = aty
        l2_d = fproc.common_opener(dp=dp_l2, setp=setp_l2)
        da_l2 = l2_d['roi'].compute()
        years = da_l2.year
        d_trend = fcv.calc_lin_reg_vec(years, da_l2)
        da_l2_grad = d_trend['grad']
    ic(cg.TrackProg(cli=atyc, iter=ip.intervals).message())
    df_sc = fg.get_spancov(scp)
    t_slc = (pl.datetime(aty[0], 1, 1), pl.datetime(aty[1], 12, 31))
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
            if in_dist_bool:
                ll_lon360 = ll_lon % 360
                np_l2_grad_ll = da_l2_grad.sel(
                    lat=ll_lat, lon=ll_lon360, method='nearest').data
                check_in_l2_dist = fproc.check_in_dist(np_l2_grad_ll, grad)
            else:
                check_in_l2_dist = None
            if d_lr['data'][0] != 0:
                lr_trnd = d_lr['data'][-1] - d_lr['data'][0]
                per_chng =  lr_trnd / d_lr['data'][0] * 100
                lf_trend["grad"].append(grad)
                lf_trend["per_chng"].append(per_chng)
                lf_trend["lat"].append(ll_lat)
                lf_trend["lon"].append(ll_lon)
                lf_trend["in_dist"].append(check_in_l2_dist)
            else:
                ic()
                lf_nan_trend["grad"].append(100)
                lf_nan_trend["per_chng"].append(100)
                lf_nan_trend["lat"].append(ll_lat)
                lf_nan_trend["lon"].append(ll_lon)
                lf_nan_trend["in_dist"].append(check_in_l2_dist)
    toc = time.time() - tic; ic(toc)
    #  Allow any additional manual transforms
    plot_grad = lf_trend["per_chng"] if calc == 'percent_change' else lf_trend["grad"]
    plot_lat = lf_trend["lat"]
    plot_lon = lf_trend["lon"]
    ppar_nan.cmap = fpl.rep_color(color=ppar_nan.color, d=lf_nan_trend["grad"])
    
    fig = plt.figure(figsize=ppar.figsize)
    plt.rcParams.update({"font.family": 'Catamaran'})
    plt.rcParams.update({"font.weight": 'light'})
    plt.rcParams.update({"font.size": 10})
    ax = fpl.set_proj(fig, ppar)
    ppar.o_name = ppar.o_prefix + ppar.o_name
    ppar.title = ppar.title +  ' n=' + str(len(plot_grad))
    fpl.plot_globe_ng(plot_grad, plot_lat, plot_lon, ppar, ax=ax)
    ppar_nan.o_name = ppar_nan.o_prefix + ppar.o_name
    ax = fpl.plot_globe_ng(
        lf_nan_trend["per_chng"], lf_nan_trend["lat"], lf_nan_trend["lon"], ppar_nan, ax=ax)
    if in_dist_bool:
        not_in_dist = [not elem for elem in lf_trend["in_dist"]]
        plot_grad_in_dist = list(compress(plot_grad, not_in_dist))
        plot_lat_in_dist = list(compress(plot_lat, not_in_dist))
        plot_lon_in_dist = list(compress(plot_lon, not_in_dist))
        ppar_in_dist.o_name = ppar_in_dist.o_prefix + ppar.o_name
        ppar_in_dist.title = ppar.title
        fpl.plot_globe_ng(plot_grad_in_dist, plot_lat_in_dist, plot_lon_in_dist, ppar_in_dist, ax=ax)
    plt.close()

if ppar.anim_flag:
    fpl.movie_maker(ppar, ip.intervals, scp, yr_str)
if ppar_in_dist.anim_flag:
    fpl.movie_maker(ppar_in_dist, ip.intervals, scp, yr_str)