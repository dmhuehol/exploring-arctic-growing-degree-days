''' wrap_calc_crossover_grid
Script to calculate crossover year for a given threshold based on
exceedance for gridded data and save as a netCDF file. Output file will
have dimensions lat x lon or rlz x lat x lon.
'''
import random
import sys

from icecream import ic
import numpy as np
import xarray as xr

import fun_calc_var as fcv
import fun_process as fproc
import classes_gddt as cg

random_select = False
dp_exceed = cg.DataParams(
    path='/Users/danielhueholt/data/gddt_data/LENS2/exceedance/',
    tok='allmembers_gexc_1850-2100_base0-2000.nc', var='exceedance',
    flag_raw_ds=True, flag_raw_da=True, flag_time_slice=True,
    flag_manage_rlz=True, flag_land_mask=True, flag_roi=False)
setp_exceed = cg.SetParams(
    area_stat='pass', base_yrs=[0, 2000], mask_flag='none', reg_oi='global',
    rlz='all', window=10, yrs=[1850, 2100])
#  Forced crossover:
#      forced_dict=dict(bool=True, threshold=(80,))
#      member_dict=dict(bool=False,)
#  Member crossover:
#      forced_dict=dict(bool=False, )
#      member_dict=dict(bool=True, threshold=(90,))
#  No-analog state:
#      forced_dict=dict(bool=False, )
#      member_dict=dict(bool=True, threshold=(100,))
ppar_crossover = cg.PlotParams(
    o_bool=True, o_name='',
    o_path='/Users/danielhueholt/data/gddt_data/LENS2/exceedance/crossover/',
    o_prefix='',
    plot_crossover_dict=dict(
        forced_dict=dict(bool=True, threshold=(90,)),
        member_dict=dict(bool=False, threshold=(100,),),),)
assert ppar_crossover.o_path is not None
assert ppar_crossover.o_prefix is not None
assert ppar_crossover.o_name is not None

exceed_dict = fproc.common_opener(dp=dp_exceed, setp=setp_exceed)
da_exceed = exceed_dict["manage_rlz"].compute()
out_attrs = da_exceed.attrs
if random_select is not False:
    random_members = random.sample(range(50), random_select)
    da_exceed = da_exceed.isel(realization=random_members)
    out_attrs["members"] = random_members
    ppar_crossover.o_prefix = 'rs10-smoothed_'
if np.max(da_exceed.data) > 100:
    out_attrs['length_of_base_period'] = exceed_dict[
        'raw_ds'].length_of_base_period
    da_exceed = da_exceed / exceed_dict['raw_ds'].length_of_base_period * 100
else:
    out_attrs['length_of_base_period'] = 'unknown'
if setp_exceed.window is not None:
    da_rolling = da_exceed.rolling(
        year=setp_exceed.window).mean().dropna('year')
    out_attrs['rolling_window'] = setp_exceed.window
else:
    da_rolling = da_exceed
    out_attrs['rolling_window'] = 'none'
lat = da_rolling.lat.data
lon = da_rolling.lon.data
rlzs = da_rolling.realization.data
years = da_rolling.year.data
if ppar_crossover.plot_crossover_dict['member_dict']['bool']:
    np_crossover = np.full((len(rlzs), len(lat), len(lon)), np.nan)
elif ppar_crossover.plot_crossover_dict['forced_dict']['bool']:
    np_crossover = np.full((len(lat), len(lon)), np.nan)
else:
    raise ValueError('Either member or forced crossover must be True')

#  Reversing years for loop ensures the result is the first year beyond
#  the threshold (crossover).
reverse_years = np.flip(years)
list_crossover_rlz = list()
bmb_type = dp_exceed.tok.split('_')[0]
for year_count, year in enumerate(reverse_years):
    if ppar_crossover.plot_crossover_dict['forced_dict']['bool']:
        loop_gexc = da_rolling.sel(year=year).mean(dim='realization')
        if year_count == 0:
            crossover_threshold = ppar_crossover.plot_crossover_dict[
                'forced_dict']['threshold'][0]
            msg_crossover = 'Calculating forced crossover at ' \
                + str(crossover_threshold) + '% threshold'
            ppar_crossover.o_name = bmb_type + '_' + 'forcedcrossover_' \
                + ppar_crossover.o_name + 'threshold' \
                + str(crossover_threshold)
            out_attrs["threshold"] = crossover_threshold
            out_attrs["type_of_crossover"] = 'forced'
            ic(msg_crossover)
        crossover = fcv.calc_crossover(loop_gexc, crossover_threshold)
        np_crossover[crossover] = year
    elif ppar_crossover.plot_crossover_dict['member_dict']['bool']:
        if year_count == 0:
            crossover_threshold = ppar_crossover.plot_crossover_dict[
                'member_dict']['threshold'][0]
            msg_crossover = 'Plotting member crossover at ' \
                + str(crossover_threshold) + '% threshold'
            ppar_crossover.o_name = 'membercrossover_' \
                + ppar_crossover.o_name + 'threshold' \
                + str(crossover_threshold)
            out_attrs["threshold"] = crossover_threshold
            out_attrs["type_of_crossover"] = 'member'
            ic(msg_crossover)
        loop_gexc = da_rolling.sel(year=year)
        crossover = fcv.calc_crossover(loop_gexc, crossover_threshold)
    else:
        ic('Warning: Unknown crossover setting (check output)')
    np_crossover[crossover] = year
out_attrs["units"] = 'year'
#  Assume that data has realization dimension
try:
    ds_crossover = xr.Dataset(
        data_vars = {
            "crossover": (("realization", "lat", "lon"), np_crossover)},
        coords = {
            "realization": rlzs,
            "lat": lat,
            "lon": lon},
        attrs = out_attrs,
    )
#  If no realization dimension
except ValueError:
    ds_crossover = xr.Dataset(
        data_vars = {"crossover": (("lat", "lon"), np_crossover)},
        coords = {
            "lat": lat,
            "lon": lon},
        attrs = out_attrs,
    )

if ppar_crossover.o_bool:
    out_name_full = ppar_crossover.o_path + ppar_crossover.o_prefix \
        + ppar_crossover.o_name
    ic(ds_crossover, out_name_full)
    ds_crossover.to_netcdf(out_name_full + '.nc')
else:
    ic('No data saved')
