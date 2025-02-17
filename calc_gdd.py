''' calc_gdd
Calculate growing degree days from daily temperature data for an input 
base value and save as netCDF.

Configured to work with the following data types:
1) Output from the Community Earth System Model (e.g., LENS2)
2) Station data from the Global Historical Climatology Network (GHCN)
3) Station data from the Danish Meteorological Institute (DMI)
4) Reanalysis from the European Reanalysis 5 (ERA5)

In the work which this codebase accompanies, only 1) and 2) are used.
Additional datasets were part of exploratory analysis and may prove 
useful again in the future.

Sum, number, max, datetimes
'''
########################################################################
import glob
import sys

from icecream import ic
import numpy as np
import xarray as xr

import classes_gddt as cg
import fun_calc_var as fcv

d_d_gen = {
    "p": '/Volumes/Polycrystal/Data/LENS2/daily_TREFHT/refined/',
    "tok": '*.nc',
    "var": 'TREFHT'
}
d_d_hpc_lens2 = {
    "p": '/barnes-engr-scratch1/DATA/CESM2-LE/temp/refined_nomask/',
    "tok": '*.nc',
    "var": 'TREFHT'
}
dp_dmi = cg.DataParams(
    path='/Users/dhueholt/Documents/gddt_data/Greenland_Obs/',
    tok='*.nc', var='mean_temp', flag_raw_ds=True, flag_raw_da=True, 
    flag_time_slice=True, flag_manage_rlz=False, flag_land_mask=False,
    flag_roi=False
)
dp_era5 = cg.DataParams(
    path='/Users/dhueholt/Documents/gddt_data/ERA5/daily_temp/merged/',
    tok='*.nc', var='t2m', flag_raw_ds=True, flag_raw_da=True, 
    flag_time_slice=True, flag_manage_rlz=False, flag_land_mask=False,
    flag_roi=False)
#  ndd_flag False for model output or reanalysis, True for station obs
gddp = cg.GddParams(
    base=5, ndd_flag=True, out_flag=True, out_filename=None,
    out_path='/Users/dhueholt/Documents/gddt_data/gdd/dmi/',
    # out_path='/barnes-engr-scratch1/DATA/CESM2-LE/processed_data/annual/gdd/reproc_20241015/',
    out_var_name='auto')

# out_d = {
#     "flag": True,
#     "p": '/Users/dhueholt/Documents/gddt_data/gdd/era5/',
#     "o_name": None,
#     "var": 'gdd' + str(gdd_d["base"]) + '_sum'
# }
# out_d_hpc = {
#     "flag": True,
#     # "p": '/barnes-engr-scratch1/DATA/CMIP6/processed/CMIP/piControl/day_1/gdd_base5/',
#     "p": '/barnes-engr-scratch1/DATA/CESM2-LE/processed_data/annual/gdd/reproc_20241015/',
#     "o_name": None,
#     "var": 'gdd' + str(gdd_d["base"]) + '_sum'
# }
dp = dp_dmi
# out_d = out_d_hpc

globs = sorted(glob.glob(dp.path + dp.tok))
for f in globs:
    gddp.out_var_name = 'gdd' + str(gddp.base) + '_sum'
    ic(f)
    da_in = xr.open_dataset(f)[dp.var]
    da_gdd_sum = fcv.gdd(da_in, gddp)
    da_ndd_sum = fcv.ndd(da_in)
    try:
        #  Adds station data for GHCN data
        station = da_in.attrs['station']
    except (AttributeError, KeyError):
        station = 'n/a'
    
    if gddp.ndd_flag:
        #  Gridded data: CESM output, ERA5 reanalysis
        try:
            ds_out = xr.Dataset(
                {gddp.out_var_name: (
                    ("year", "lat", "lon"), da_gdd_sum.data),
                    "ndd": (("year", "lat", "lon"), da_ndd_sum.data)},
                coords=dict(
                    lon=da_gdd_sum.lon.data,
                    lat=da_gdd_sum.lat.data,
                    year=da_gdd_sum.year.data),
                attrs=dict(
                    long_name='Growing degree day information',
                )
            )
        #  Station data: GHCN, DMI
        except (IndexError, AttributeError):
            ds_out = xr.Dataset(
                {gddp.out_var_name: (("year",), da_gdd_sum.data),
                    "ndd": (("year",), da_ndd_sum.data)},
                coords=dict(
                    year=da_gdd_sum.year.data),
                attrs=dict(
                    long_name='Growing degree day information',
                    station=station
                )
            )
    else:
        #  Gridded data: CESM output, ERA5 reanalysis
        try:
            ds_out = xr.Dataset(
                {gddp.out_var_name: (
                    ("year", "lat", "lon"), da_gdd_sum.data),},
                coords=dict(
                    lon=da_gdd_sum.lon.data,
                    lat=da_gdd_sum.lat.data,
                    year=da_gdd_sum.year.data),
                attrs=dict(
                    long_name='Growing degree day information',
                )
            )
        #  Station data: GHCN, DMI
        except (IndexError, AttributeError):
            ds_out = xr.Dataset(
                {gddp.out_var_name: (("year",), da_gdd_sum.data),},
                coords=dict(
                    year=da_gdd_sum.year.data),
                attrs=dict(
                    long_name='Growing degree day information',
                    station=station
                )
            )
    
    if gddp.out_flag:
        if gddp.out_filename is None:
            if dp.var != gddp.out_var_name:
                o_name = f.split('/')[-1].replace(dp.var, gddp.out_var_name)
            else:
                o_name = f.split('/')[-1]
        else:
            o_name = gddp.out_filename
        ds_out.to_netcdf(gddp.out_path + o_name)