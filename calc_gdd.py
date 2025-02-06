''' calc_gdd
Calculate growing degree days for an input base value and save as
netCDF.
Sum, number, max, datetimes
Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''
########################################################################
import glob
import sys

from icecream import ic
import numpy as np
import xarray as xr

import fun_calc_var as fcv

d_d = {
    "p": '/Volumes/Polycrystal/Data/LENS2/daily_TREFHT/refined/',
    "tok": '*.nc',
    "var": 'TREFHT'
}
d_d_hpc_lens2 = {
    "p": '/barnes-engr-scratch1/DATA/CESM2-LE/temp/refined_nomask/',
    "tok": '*.nc',
    "var": 'TREFHT'
}
d_d_dmi = {
    "p": '/Users/dhueholt/Documents/gddt_data/Greenland_Obs/',
    "tok": '*34360_4360*.nc',
    "var": 'mean_temp'
}
d_d_era5 = {
    "p": '/Users/dhueholt/Documents/gddt_data/ERA5/daily_temp/merged/',
    "tok": '*.nc',
    "var": 't2m'
}
gdd_d = {
    "base": 5,
    "ndd_flag": False,
}
out_d = {
    "flag": True,
    "p": '/Users/dhueholt/Documents/gddt_data/gdd/era5/',
    "o_name": None,
    "var": 'gdd' + str(gdd_d["base"]) + '_sum'
}
out_d_hpc = {
    "flag": True,
    # "p": '/barnes-engr-scratch1/DATA/CMIP6/processed/CMIP/piControl/day_1/gdd_base5/',
    "p": '/barnes-engr-scratch1/DATA/CESM2-LE/processed_data/annual/gdd/reproc_20241015/',
    "o_name": None,
    "var": 'gdd' + str(gdd_d["base"]) + '_sum'
}
d_d = d_d_hpc_lens2
out_d = out_d_hpc

globs = sorted(glob.glob(d_d['p'] + d_d['tok']))
for f in globs:
    ic(f)
    da_in = xr.open_dataset(f)[d_d['var']]
    da_gdd_sum = fcv.gdd(da_in, gdd_d)
    da_ndd_sum = fcv.ndd(da_in, gdd_d)
    try:
        station = da_in.attrs['station']
    except:
        station = 'n/a'
    
    if gdd_d['ndd_flag']:
        try:
            ds_out = xr.Dataset( # Gridded data
                {out_d['var']: (("year", "lat", "lon"), da_gdd_sum.data),
                    "ndd": (("year", "lat", "lon"), da_ndd_sum.data)},
                coords=dict(
                    lon=da_gdd_sum.lon.data,
                    lat=da_gdd_sum.lat.data,
                    year=da_gdd_sum.year.data),
                attrs=dict(
                    long_name='Growing degree day information',
                )
            )
        except:
            ds_out = xr.Dataset( # Station data
                {out_d['var']: (("year",), da_gdd_sum.data),
                    "ndd": (("year",), da_ndd_sum.data)},
                coords=dict(
                    year=da_gdd_sum.year.data),
                attrs=dict(
                    long_name='Growing degree day information',
                    station=station
                )
            )
    else:
        try:
            ds_out = xr.Dataset( # Gridded data
                {out_d['var']: (("year", "lat", "lon"), da_gdd_sum.data),},
                coords=dict(
                    lon=da_gdd_sum.lon.data,
                    lat=da_gdd_sum.lat.data,
                    year=da_gdd_sum.year.data),
                attrs=dict(
                    long_name='Growing degree day information',
                )
            )
        except:
            ds_out = xr.Dataset( # Station data
                {out_d['var']: (("year",), da_gdd_sum.data),},
                coords=dict(
                    year=da_gdd_sum.year.data),
                attrs=dict(
                    long_name='Growing degree day information',
                    station=station
                )
            )
    
    if out_d['flag']:
        if out_d['o_name'] is None:
            if d_d['var'] != out_d['var']:
                o_name = f.split('/')[-1].replace(d_d['var'], out_d['var'])
            else:
                o_name = f.split('/')[-1]
        else:
            o_name = out_d['o_name']
        ds_out.to_netcdf(out_d['p'] + o_name)