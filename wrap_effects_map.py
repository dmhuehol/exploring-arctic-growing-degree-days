'''wrap_effects_map
Script to plot a spatial map of some effect size statistic.

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University 
'''
import sys

import cmasher as cmr
from icecream import ic
import matplotlib.colors
from matplotlib import colormaps as cm
import matplotlib.pyplot as plt
import numpy as np
import rioxarray
from shapely import concave_hull
from shapely.geometry import mapping
import xarray as xr

import classes_gddt as cg
import fun_plots as fpl
import fun_process as fproc
import gddt_region_library as g_rlib

# TODO: update input handling so everything happens inside Params objects
aggregate_bool = False
aggregate = g_rlib.atlas_ecoregions2017()
paint_shapefile_bool = True

# cmap = cmr.get_sub_cmap('cmr.lavender_r', 0, 0.9, N=12)
cmap = cmr.get_sub_cmap('cmr.bubblegum_r', 0, 1, N=12)

# cmap = matplotlib.colors.Normalize(vmin=-1, vmax=1, clip=False)

cmn_path = '/Users/dhueholt/Documents/gddt_fig/20241208_furtheraes/'
dp_gdd_gexc = cg.DataParams(
    path='/Users/dhueholt/Documents/gddt_data/LENS2/effects/', 
    tok='gexc_1850-2100_base0-2000.nc', var='gexc', flag_raw_ds=False, flag_raw_da=True, 
    flag_time_slice=True, flag_manage_rlz=True, 
    flag_land_mask=True, flag_roi=False)
dp_altmax = cg.DataParams(
    path='/Users/dhueholt/Documents/gddt_data/LENS2/annual_ALTMAX/',
    tok='*.nc', var='ALTMAX', flag_raw_ds=False, flag_raw_da=True, 
    flag_time_slice=False, flag_manage_rlz=False, 
    flag_land_mask=False, flag_roi=False)
setp_gdd_gexc = cg.SetParams(
    area_stat=None, beat=1999, mask_flag=['land',], 
    effect='gexc', reg_oi='global', rho=50, rlz='all', window=None, 
    yrs=[1850, 2100])
setp_altmax = cg.SetParams(mask_flag='land', reg_oi='global', rlz='all')
ppar = cg.PlotParams(
    cb_bool=True, cb_extent='neither', cb_label='crossover year', cb_vals=[2040, 2100],
    cmap=cmap, dpi=800, figsize=(5, 4), o_bool=True,
    o_name='', o_path=cmn_path, o_prefix='', plot_all=False, 
    proj='Arctic', title='', title_size=10)
ppar_altmask = cg.PlotParams(
    alpha=0.6, cb_bool=False, cb_extent='neither', cb_label='auto', cb_vals=[-10, 0],
    cmap=cm['Greys'], dpi=800, figsize=(5, 4), o_bool=True,
    o_name='', o_path=cmn_path, 
    o_prefix='', 
    plot_all=False, 
    proj='Arctic', set_bad=False, title=': median crossover year', title_size=11)
ppar_sf = cg.PlotParams(
    cmap=cmap,
    o_prefix='testsf_', title='Arctic ecoregions by crossover year', title_size=12)

gexc_d = fproc.common_opener(dp=dp_gdd_gexc, setp=setp_gdd_gexc)
da_all = gexc_d['land_mask'].compute()
ll = (da_all.lat.data, da_all.lon.data)
years = da_all.year.data
np_crossover = np.full((len(ll[0]), len(ll[1])), np.nan)
l_ind = list()
#  Reversing years ensures crossover years aren't progressively overwritten
reverse_years = np.flip(years)
if setp_gdd_gexc.effect == 'gexc':
    if setp_gdd_gexc.beat == 1999:
        ppar.o_name = 'LENS2_ensmnnoanalogyear'
        ppar.title = 'LENS2 no-analog emergence in ensemble mean'
    else:
        ppar.o_name = 'LENS2_ensmnrossoveryear'
        ppar.title = 'LENS2 ensemble mean crossover (>80% of Preindustrial samples)'
    da_gexc_mn = da_all.mean(dim='realization')
    for yrc, yr in enumerate(reverse_years):
        act_gexc = da_gexc_mn.sel(year=yr)
        act_crossover = act_gexc > setp_gdd_gexc.beat
        np_crossover[act_crossover] = yr
elif setp_gdd_gexc.effect == 'robustness':
    ppar.o_name = 'LENS2_mediancrossoveryear'
    ppar.title = 'LENS2 median crossover'
    for yrc, yr in enumerate(reverse_years):
        act_gexc = da_all.sel(year=yr)
        count_robust = fproc.threshold_count(act_gexc, beat=setp_gdd_gexc.beat)
        act_crossover = count_robust > setp_gdd_gexc.rho
        np_crossover[act_crossover] = yr
da_crossover = xr.DataArray(
    data=np_crossover,
    dims=["lat", "lon"],
    coords = {"lat": ll[0], "lon": ll[1],},
    attrs = {"units": 'year'}
)
nd = fproc.namer(gexc_d["raw_da"], setp_gdd_gexc)

if aggregate_bool:
    np_fill = np.full(np.shape(da_crossover.data), np.nan)
    for regc, reg in enumerate(aggregate):
        setp_act = cg.SetParams(reg_oi=reg, area_stat=False)
        da_act, loc_str, loc_title = fproc.manage_area(da_crossover, setp_act)
        setp_act_mn = cg.SetParams(reg_oi=reg, area_stat='mean')
        da_act_val, _, _ = fproc.manage_area(da_crossover, setp_act_mn)
        da_act_val_ll = da_act.where(np.isnan(da_act.data), other=da_act_val.data)
        np_fill[~np.isnan(da_act_val_ll)] = da_act_val.data
        # da_act_val_llnonan = da_act_val_ll.fillna(np.inf).argmin(dim=['lat','lon'])
        # ic(da_act_val_llnonan)
        # sys.exit('STOP')
    da_fill = xr.DataArray(
        data=np_fill,
        dims=["lat", "lon"],
        coords = {"lat": ll[0], "lon": ll[1],},
        attrs = {"units": 'year'}
    )
    plot_this = da_fill
    plt.rcParams.update({'font.family': 'Catamaran'})
    #  'font.weight': normal, bold, heavy, light, ultrabold, ultralight
    fig, ax = fpl.plot_globe(plot_this, ppar)
elif paint_shapefile_bool:
    ppar.o_name = ppar.o_prefix + nd['data_id'] + '_' + ppar.o_name
    geo_arctic, geo_other = fproc.get_arctic_biomes()
    geo_arctic_nps = geo_arctic.to_crs(epsg=4326)
    no_plot_proj = geo_other.to_crs(epsg=3995)
    # ic(geo_arctic_nps)
    da_crossover.coords['lon'] = (da_crossover.coords['lon'] + 180) % 360 - 180
    da_crossover = da_crossover.sortby(da_crossover.lon)
    da_crossover.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
    da_crossover.rio.write_crs("epsg:4326", inplace=True)
    l_crossover = list()
    for ri in np.arange(0, len(geo_arctic_nps)):
        ecoreg = geo_arctic_nps.iloc[ri]
        # ecoreg_hull = concave_hull(ecoreg.geometry)
        # ic(len(geo_arctic))
        # ic(ecoreg)
        da_crossover_ecoreg = da_crossover.rio.clip(
            [ecoreg.geometry], geo_arctic_nps.crs, all_touched=True, drop=False)
        lat_wght = np.cos(np.deg2rad(da_crossover_ecoreg.lat))
        da_wght = da_crossover_ecoreg.weighted(lat_wght)
        da_regmn = da_wght.mean(dim=['lat', 'lon'], skipna=True)
        l_crossover.append(da_regmn.data)
    ppar.color = l_crossover
    plot_this = geo_arctic_nps.to_crs(epsg=3995)
    plt.rcParams.update({'font.family': 'Catamaran'})
    #  'font.weight': normal, bold, heavy, light, ultrabold, ultralight
    plt.rcParams.update({'font.weight': 'heavy'})
    plt.rcParams.update({'font.size': 12})
    fpl.plot_globe_shapefile(plot_this, no_plot_proj, ppar)
else:
    plot_this = da_crossover
    ppar.o_name = ppar.o_prefix + ppar.o_name
    ic(ppar.title, ppar.o_name)
    plt.rcParams.update({'font.family': 'Catamaran'})
    #  'font.weight': normal, bold, heavy, light, ultrabold, ultralight
    plt.rcParams.update({'font.weight': 'light'})
    plt.rcParams.update({'font.size': 10})  
    fig, ax = fpl.plot_globe(plot_this, ppar)
if 'permafrost' in setp_gdd_gexc.mask_flag:
    ic()
    altmax_d = fproc.common_opener(dp=dp_altmax, setp=setp_altmax)
    da_altmax = altmax_d['raw_da'].compute()
    min_altmax = da_altmax.min(dim=['realization', 'time'], skipna=True)
    no_trees = 0.06
    altmask = min_altmax.where(min_altmax < no_trees)
    mute_this = altmask
    ppar_altmask.title = ppar.title
    ppar_altmask.o_name = ppar.o_name
    fpl.plot_globe(mute_this, ppar_altmask, ax=ax, fig=fig)
    # fpl.plot_globe(max_altmax, ppar_altmask)