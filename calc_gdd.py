''' calc_gdd
Calculate growing degree days from daily temperature data for an input 
base value and save as netCDF.

Configured for the following data types:
1) Output from the Community Earth System Model (e.g., LENS2)
2) Station data from the Danish Meteorological Institute (DMI)
3) Reanalysis from the European Reanalysis 5 (ERA5)
GDDs for Global Historical Climatology Network (GHCN) stations 
calculated in calc_gdd_info_ghcn.

Only 1) is directly used in the paper this codebase accompanies. 2) and 
3) were part of exploratory analysis and may prove useful again in the 
future.
'''
import glob
import sys

from icecream import ic
import xarray as xr

import classes_gddt as cg
import fun_calc_var as fcv

dp_cesm = cg.DataParams(
    path='/Volumes/Polycrystal/Data/LENS2/daily_TREFHT/refined/',
    tok='*arc*.nc', var='TREFHT', flag_raw_ds=True, flag_raw_da=True, 
    flag_time_slice=True, flag_manage_rlz=False, flag_land_mask=False,
    flag_roi=False
)
dp_cesm_hpc = cg.DataParams(
    path='/barnes-engr-scratch1/DATA/CESM2-LE/temp/refined_nomask_orig/',
    tok='*arc*.nc', var='TREFHT', flag_raw_ds=True, flag_raw_da=True, 
    flag_time_slice=True, flag_manage_rlz=False, flag_land_mask=False,
    flag_roi=False
)
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
    base=5, ndd_flag=False, out_flag=True, out_filename=None,
    # out_path='/Users/dhueholt/Documents/gddt_data/gdd/lens2/',
    out_path='/barnes-engr-scratch1/DATA/CESM2-LE/processed_data/annual/gdd/reproc_20250218/',
    out_var_name='auto')
dp = dp_cesm_hpc

globs = sorted(glob.glob(dp.path + dp.tok))
for f in globs:
    gddp.out_var_name = 'gdd' + str(gddp.base) + '_sum'
    ic(f)
    da_in = xr.open_dataset(f)[dp.var]
    da_gdd_sum = fcv.gdd(da_in, gddp)
    da_ndd_sum = fcv.ndd(da_in)
    try:
        #  Adds station metadata for GHCN observations
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