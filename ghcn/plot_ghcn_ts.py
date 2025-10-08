'''plot_ghcn_ts
Plot GHCN station data as a timeseries.

This supports two main modes of operation:
    1) Plot timeseries for all stations matching span and coverage 
    parameters in a given time period; useful for exploring dataset
    2) Plot the whole record for all station(s) given station ID(s)

Run with default settings, this runs in mode 2) to reproduce the 
subpanels of Figure 5 in the manuscript Hueholt et al. "Exploring the 
Influence of Internal Climate Variability and Forced Change on Arctic 
Greening" which this codebase accompanies.

To run in mode 1), change ip to generate the desired intervals (the 
commented line demonstrates for 1993-2022), change stid to 'auto', and
change plot_require_whole to True.

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''
import sys
sys.path.append(
    "/Users/danielhueholt/Documents/GitHub/" 
    + "exploring-arctic-growing-degree-days/")

from icecream import ic
import matplotlib.pyplot as plt
import numpy as np
import polars as pl

import classes_gddt as cg
import fun_calc_var as fcv
import fun_ghcn as fg
import fun_plots as fp

d_path = "/Users/danielhueholt/Documents/Data/gddt_data/Arctic50N_gt1/"
sc_f = "~/Documents/Data/ghcn_data/spancov_lists/spancov_Arctic50N_gt1.csv"
scp = cg.SpanCovParams(f=sc_f, por=30, cov_thr=97, cov_type='>')
#  As written, ip effectively plots whole record for every station
ip = cg.IntervalParams(span=224, strt_yr=1800, end_yr=2024, type='noverlap')
# ip = cg.IntervalParams(span=30, strt_yr=1993, end_yr=2022, type='noverlap')
#  True/False to plot linear regression line
plot_linearreg = False
#  Require time period of record to match entire interval to plot
plot_require_whole = False
#  List of station IDs, or 'auto' to plot all available; matches Fig. 5
stid = ['UKE00105900', 'CA002403206', 'SWE00139842', 'SWE00138226',
        'RSM00021982', 'RSM00024688']
o_path = '/Users/danielhueholt/Documents/Figures/arc-gdd_fig/20251008_spinup/'

df_sc = fg.get_spancov(scp)
n_st = ic(len(df_sc))
if stid == 'auto':
    stid = df_sc["st_id"]
for atyc, aty in enumerate(ip.intervals):
    for cc, act_st in enumerate(stid):
        if (cc % 100 == 0) & (cc != 0):
            ic(cg.TrackProg(cli=cc, iter=stid).message())
        #  Manual station code
        # act_st = 'UKE00105900'
        df_act = df_sc.filter(pl.col('st_id') == act_st)
        lat = df_act.select(pl.col('lat')).item()
        lon = df_act.select(pl.col('lon')).item()
        act_ll = str(lat) + 'N,' + str(lon) + 'E'
        act_cov = round(df_act.select(pl.col('coverage')).item(), 2)
        dp = cg.DataParams(path=d_path, tok=act_st + '_refined.csv')
        cmn_name = act_st + '_ts'
        cmn_title = act_st + ' ' + act_ll
        df_ghcn = pl.read_csv(dp.path + dp.tok)
        df_gdd = df_ghcn.filter(pl.col('GDD_SUM').is_not_nan())
        df_dt = fg.add_dtns(df_gdd)
        t_slc = (pl.datetime(aty[0], 1, 1), pl.datetime(aty[1], 12, 31))
        df_tslc = df_dt.filter(
            pl.col("DATETIME").is_between(t_slc[0], t_slc[1]))
        np_gdd = df_tslc.select(pl.col('GDD_SUM')).to_numpy()
        np_yrs = df_tslc.select(
            pl.col('DATETIME')).to_series().dt.year().to_numpy()
        if cg.Sentinel(data=np_yrs, match=aty).not_in_span():
            if plot_require_whole:
                continue
        if act_cov < 100:
            np_ndd = df_tslc.select(pl.col('NDD_SUM')).to_numpy()
        d_lr = fcv.calc_lin_reg(np_yrs, np_gdd)
        lr_data = d_lr['grad'] * np_yrs + d_lr['intcpt']
        
        plt.rcParams.update({'font.family': 'Catamaran'})
        #  'font.weight': normal, bold, heavy, light, ultrabold, ultralight
        plt.rcParams.update({'font.weight': 'light'})
        plt.rcParams.update({'font.size': 10})
        fig, ax1 = plt.subplots()
        ppar_gdd = cg.PlotParams(
            color='#6f8c31', o_bool=False, o_path=o_path, o_name=cmn_name, 
            title=cmn_title, x_label='year', x_lim='auto',
            y_label='GDD base 5', y_lim='fix_nonnegative',)
        ppar_ndd = cg.PlotParams(
            color='#cccccc', o_bool=True, o_path=o_path, o_name=cmn_name, 
            title=cmn_title, y_label='Missing days', y_lim=[0, 366], 
            yticks=np.arange(0, 365 + 73, 73))
        #  By default, plot with y-tick labels; for the paper figures, 
        #  I disable these manually in fp.apply_params_ts and add them
        #  back in Keynote for aesthetics.
        fp.plot_timeseries(np_gdd, x_d=np_yrs, ppar=ppar_gdd)
        if plot_linearreg:
            ppar_lr = cg.PlotParams(
                color='#ce4a4a', label='linear', leg_bool=True, o_bool=False, 
                o_path=o_path, o_name=cmn_name, title=cmn_title, 
                x_label='year', x_lim='auto', y_label='GDD base 5')
            fp.plot_timeseries(lr_data, x_d=np_yrs, ppar=ppar_lr)
        if act_cov < 100:
            ax2 = ax1.twinx()
            fp.plot_timeseries(np_ndd, x_d=np_yrs, ppar=ppar_ndd)
        plt.close()