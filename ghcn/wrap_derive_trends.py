"""wrap_derive_trends
Derive trends and relevant information from GHCN and CESM2 data, and
save as a new file.
"""
import sys
sys.path.append(
    "/Users/danielhueholt/Documents/GitHub/" \
    + "exploring-arctic-growing-degree-days/")

import pandas as pd
import polars as pl
from icecream import ic

import classes_gddt as cg
import fun_calc_var as fcv
import fun_ghcn as fg
import fun_process as fproc
import gddt_region_library as g_rlib

include_model_trends = True
data_path = "/Users/danielhueholt/Data/gddt_data/Arctic50N_gt1/"
sc_f = "~/Data/ghcn_data/spancov_lists/spancov_Arctic50N_gt1.csv"
scp = cg.SpanCovParams(f=sc_f, por=10, cov_thr=97, cov_type=">")
ip = cg.IntervalParams(
    span=scp.por, strt_yr=1873, end_yr=2022, type="noverlap")
dp_model = cg.DataParams(
    path='/Users/danielhueholt/Data/gddt_data/gdd/lens2/', tok='*arc*.nc',
    var='gdd5_sum', flag_raw_ds=True, flag_raw_da=True, flag_time_slice=True,
    flag_manage_rlz=True, flag_land_mask=False, flag_roi=True)
setp_model = cg.SetParams(
    area_stat='pass', reg_oi=g_rlib.Arctic50N(), rlz='all')
out_path = '/Users/danielhueholt/Data/gddt_data/trends/'
out_fn = "hist_GHCN_abv50N_cov" + scp.cov_thr_str + "_span" + scp.por_str \
    + "_trend" + ip.strt_yr_str + ip.end_yr_str + '.csv'

out_columns = [
        "interval", "st_id", "lat", "lon", "slope", "percent_change",
        "intercept", "in_nearest_lens2_dist", "nan_trend_bool"]
if include_model_trends:
    rlzs = range(100)
    l_rlzs_str = ['LENS2-' + str(a) for a in rlzs]
    out_columns = pd.Series(out_columns + l_rlzs_str)
else:
    out_columns = pd.Series(out_columns)
l_interval_trends = list()
for loop_count, interval in enumerate(ip.intervals):
    yr_str = str(interval[0]) + "-" + str(interval[1])
    setp_model.yrs = interval
    dict_lens2 = fproc.common_opener(dp=dp_model, setp=setp_model)
    if dict_lens2['roi'] is None:
        raise TypeError('No data matching "roi" input')
    da_lens2 = dict_lens2['roi'].compute()
    years = da_lens2.year
    d_trend = fcv.calc_lin_reg_vec(years, da_lens2)
    da_lens2_grad = d_trend['grad']
    df_match_stations = fg.get_spancov(scp)
    t_slc = (pl.datetime(interval[0], 1, 1), pl.datetime(interval[1], 12, 31))
    stid = df_match_stations["st_id"]
    l_station_trends = list()
    for cc, active_st in enumerate(stid):
        active_fn = active_st + "_refined.csv"
        ldf_active = pl.scan_csv(data_path + active_fn)
        dict_lr = fg.trend_ghcn(ldf_active, scp, t_slc)
        grad = dict_lr["grad"]
        intcpt = dict_lr["intcpt"]
        if grad is None:
            #  Could not calculate (i.e., data not available)
            continue
        else:
            if isinstance(dict_lr['data'], type(None)):
                raise TypeError('"data" key in regression is None')
            ll = ldf_active.select(pl.first('LATITUDE', 'LONGITUDE')).collect()
            ll_lat = ll.select('LATITUDE').item()
            ll_lon = ll.select('LONGITUDE').item()
            ll_lon360 = ll_lon % 360
            np_lens2_grad_ll = da_lens2_grad.sel(
                lat=ll_lat, lon=ll_lon360, method='nearest').data
            check_in_lens2_dist = fproc.check_in_dist(np_lens2_grad_ll, grad)
            if dict_lr['data'][0] != 0:
                lr_trnd = dict_lr['data'][-1] - dict_lr['data'][0]
                percent_change = lr_trnd / dict_lr['data'][0] * 100
                nan_trend = False
            else:
                nan_trend = True
                grad = 100
                percent_change = 100
            active_st_data = [
                yr_str, active_st, ll_lat, ll_lon, grad, percent_change,
                intcpt, check_in_lens2_dist, nan_trend]
            if include_model_trends:
                #  This is horribly inefficient, but I can't immediately
                #  think of a faster approach since embedding an array
                #  as an element isn't available in pandas. It doesn't
                #  run often, so xkcd 1445 applies as well.
                for r in range(len(np_lens2_grad_ll)):
                    active_st_data.append(np_lens2_grad_ll[r])
            active_st_df = pd.DataFrame(columns=out_columns)
            active_st_df.loc[0] = active_st_data
            l_station_trends.append(active_st_df)
    df_interval_trend = pd.concat(l_station_trends)
    l_interval_trends.append(df_interval_trend)
df_all_trend = pd.concat(l_interval_trends)
if include_model_trends:
    out_fn = out_fn.replace('.csv', '_withLENS2.csv')
df_all_trend.to_csv(out_path + out_fn)
ic("Completed!")
