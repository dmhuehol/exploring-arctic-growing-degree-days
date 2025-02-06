'''plot_ghcn_ts
Plot GHCN data as a timeseries. At some point this may be wrapped into
wrap_station_ts; for now the different file formats and styles make it
easier to start out separately.

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''
import sys
sys.path.append(
    "/Users/dhueholt/Documents/Github/" 
    + "growing-degree-days-treelines-poles/")

from icecream import ic
import matplotlib.pyplot as plt
import numpy as np
import polars as pl

import classes_gddt as cg
import fun_calc_var as fcv
import fun_ghcn as fg
import fun_plots as fp

d_path = "/Users/dhueholt/Documents/gddt_data/Arctic50N_gt1/"
sc_f = "~/Documents/ghcn_data/spancov_lists/spancov_Arctic50N_gt1.csv"
scp = cg.SpanCovParams(f=sc_f, por=30, cov_thr=97, cov_type='>')
intvls = fp.create_intvls(strt_yr=1993, end_yr=2022, spn=scp.por)
o_path = '/Users/dhueholt/Documents/gddt_fig/20240709_trendsOnTs/'

df_sc = fg.get_spancov(scp)
n_st = ic(len(df_sc))
stid = df_sc["st_id"]
for atyc, aty in enumerate(intvls):
    for cc, act_st in enumerate(stid):
        if cc % 100 == 0:
            ic(cg.TrackProg(cli=cc, iter=stid).message())
        #  Manual station code
        # act_st = 'UKE00156852'
        df_act = df_sc.filter(pl.col('st_id') == act_st)
        lat = df_act.select(pl.col('lat')).item()
        lon = df_act.select(pl.col('lon')).item()
        act_ll = str(lat) + 'N,' + str(lon) + 'E'
        act_cov = round(df_act.select(pl.col('coverage')).item(), 2)
        dp = cg.DataParams(path=d_path, tok=act_st + '_refined.csv')
        cmn_name = act_st + '_ts'
        cmn_title = act_st + ' ' + act_ll
        ppar_gdd = cg.PlotParams(
            color='#ff80ed', o_bool=False, o_path=o_path, o_name=cmn_name, 
            title=cmn_title, x_label='year', x_lim='auto', y_label='GDD base 5', 
            y_lim='auto')
        ppar_ndd = cg.PlotParams(
            color='#cccccc', o_bool=True, o_path=o_path, o_name=cmn_name, 
            title=cmn_title, y_label='Missing days', y_lim=[0, 366], 
            yticks=np.arange(0, 365 + 73, 73))
        ppar_lr = cg.PlotParams(
            color='#ce4a4a', label='linear', leg_bool=True, o_bool=False, o_path=o_path, o_name=cmn_name, 
            title=cmn_title, x_label='year', x_lim='auto', y_label='GDD base 5')
        ppar_tsr = cg.PlotParams(
            color='#d07dce', label='Theil-Sen', leg_bool=True, o_bool=False, o_path=o_path, o_name=cmn_name, 
            title=cmn_title, x_label='year', x_lim='auto', y_label='GDD base 5')
        ppar_rr = cg.PlotParams(
            color='#ffb65a', label='RANSAC', leg_bool=True, o_bool=True, o_path=o_path, o_name=cmn_name, 
            title=cmn_title, x_label='year', x_lim='auto', y_label='GDD base 5')
        df_ghcn = pl.read_csv(dp.path + dp.tok)
        df_gdd = df_ghcn.filter(pl.col('GDD_SUM').is_not_nan())
        df_dt = fg.add_dtns(df_gdd)
        t_slc = (pl.datetime(aty[0], 1, 1), pl.datetime(aty[1], 12, 31))
        df_tslc = df_dt.filter(
            pl.col("DATETIME").is_between(t_slc[0], t_slc[1]))
        np_gdd = df_tslc.select(pl.col('GDD_SUM')).to_numpy()
        np_yrs = df_tslc.select(pl.col('DATETIME')).to_series().dt.year().to_numpy()
        if cg.Sentinel(data=np_yrs, match=aty).not_in_span():
            continue
        if act_cov < 100:
            np_ndd = df_tslc.select(pl.col('NDD_SUM')).to_numpy()
        d_lr = fcv.calc_lin_reg(np_yrs, np_gdd)
        lr_data = d_lr['grad'] * np_yrs + d_lr['intcpt']
        d_tsr = fcv.calc_theilsen_reg(np_yrs, np_gdd)
        tsr_data = d_tsr['grad'] * np_yrs + d_tsr['intcpt']
        est = fcv.calc_ransac_reg(np_yrs, np_gdd)
        rr_data = est.predict(np_yrs.reshape(-1, 1))
        
        plt.rcParams.update({'font.family': 'Catamaran'})
        #  'font.weight': normal, bold, heavy, light, ultrabold, ultralight
        plt.rcParams.update({'font.weight': 'light'})
        plt.rcParams.update({'font.size': 10})
        fig, ax1 = plt.subplots()
        fp.plot_timeseries(np_gdd, x_d=np_yrs, ppar=ppar_gdd)
        fp.plot_timeseries(lr_data, x_d=np_yrs, ppar=ppar_lr)
        fp.plot_timeseries(tsr_data, x_d=np_yrs, ppar=ppar_tsr)
        fp.plot_timeseries(rr_data, x_d=np_yrs, ppar=ppar_rr)
        if act_cov < 100:
            ax2 = ax1.twinx()
            fp.plot_timeseries(np_ndd, x_d=np_yrs, ppar=ppar_ndd)
        plt.close()