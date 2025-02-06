'''refine_ghcn
Refine Global Historical Climatology Network station data by adding 
average daily temperature, datetimes, growing degree day, and
annual number of NaN days to all files in a region.

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''
import sys
sys.path.append(
    "/Users/dhueholt/Documents/Github/" 
    + "growing-degree-days-treelines-poles/")

from icecream import ic
import numpy as np
import pandas as pd
import xarray as xr

import fun_calc_var as fcv

inv_path = "/Users/dhueholt/Documents/ghcn_data/station_lists/"
inv_name = "ghcn_Arctic50N_TMIN-TMAX_>1.csv"
gdd_d = {
    "base": 5,
    "ndd_flag": True,
}
d_path = "/Volumes/Polycrystal/Data/ghcn/daily-summaries-latest/"
save_flag = True
out_path = "/Volumes/Polycrystal/Data/ghcn/refined-data/Arctic50N_gt1/"

df_inv_roi = pd.read_csv(inv_path + inv_name)
#  Unique IDs ensure no double-saving of data
df_ids_roi = np.unique(df_inv_roi.iloc[:, 0])
notemp_stn_list = list()
#  Loops are not pandas style, but required as it is not reasonable to
#  open all files at once.
for act_st in df_ids_roi:
    ic(act_st)
    act_name = act_st + '.csv'
    df_act = pd.read_csv(d_path + act_name)
    try:
        #  Convert temp data from tenths degree to full deg Celsius
        df_act['TMAX'] = df_act['TMAX'] / 10
        df_act['TMIN'] = df_act['TMIN'] / 10
    except KeyError as ke:
        #  Can't shortcut rest of loop in case TAVG is recorded!
        pass
    try:
        df_act['TAVG'] = df_act['TAVG'] / 10
    except KeyError:
        tavg_msg = 'TAVG not recorded by station ' + act_st + '. Calculating ' \
            + 'from TMAX and TMIN.'
        ic(tavg_msg)
        try:
            df_act['TAVG'] = (df_act['TMAX'] + df_act['TMIN']) / 2
        except KeyError as ke:
            ke_str = str(ke)
            missing_msg = ke_str + ' not recorded by station ' + act_st \
                + ' either. Skipping this station!'
            ic(missing_msg)
            notemp_stn_list.append(act_st)
            continue
    df_act['DATETIME'] = pd.to_datetime(df_act['DATE'])
    df_act.index = df_act['DATETIME']
    last_dt = df_act['DATETIME'].iloc[-1]
    first_dt = df_act['DATETIME'].iloc[0]
    df_act_ri = df_act.reindex(
        pd.date_range(first_dt, last_dt), fill_value=np.nan)
    df_act_ri['DATETIME'] = pd.date_range(first_dt, last_dt)
    #  GDD calc relies on an input DataArray. Profiling indicates minimal
    #  overhead added from this restructuring!
    da_tavg = xr.DataArray(
        data=df_act_ri['TAVG'],
        dims=["time"],
        coords=dict(
            time=df_act_ri['DATETIME']
        )
    )
    da_gdd_sum = fcv.gdd(da_tavg, gdd_d)
    da_ndd_sum = fcv.ndd(da_tavg, gdd_d)
    yrs = df_act_ri['DATETIME'].dt.year
    nonan_yrs = df_act_ri['DATETIME'].dt.year[~np.isnan(yrs)]
    _, topofyr_ind = np.unique(nonan_yrs, return_index=True)
    complete_dims = np.shape(df_act_ri['DATETIME'])
    sum_gdd = np.full(complete_dims, np.nan)
    sum_gdd[topofyr_ind] = da_gdd_sum
    sum_ndd = np.full(complete_dims, np.nan)
    sum_ndd[topofyr_ind] = da_ndd_sum
    df_act_ri['GDD_SUM'] = sum_gdd
    df_act_ri['NDD_SUM'] = sum_ndd
    if save_flag:
        out_name = act_name.replace(".csv", "_refined.csv")
        df_act_ri.to_csv(out_path + out_name)

if save_flag:
    np_notemp = np.array(notemp_stn_list)
    roi_str = inv_name.split('_')[1]
    notemp_stn_f = 'notemp_stn_' + roi_str + '.dat'
    np.savetxt(out_path + notemp_stn_f, np_notemp, delimiter=",", fmt="%s")