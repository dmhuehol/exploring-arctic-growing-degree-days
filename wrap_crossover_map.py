"""wrap_crossover_map
Plot crossover year on a map.

Supports analysis of:
    1) Forced crossover
        Defined when the ensemble mean exceeds a certain threshold of
        the preindustrial base period; by default 80% of samples,
        although others can be used.
    2) Member crossover
        Defined for each ensemble member when it exceeds a certain
        threshold of the preindustrial base period; by default 90% of
        samples.
    3) No-analog state
        This is a special case of member crossover for a threshold of
        100% exceedance of the preindustrial base period.

This requires a pre-calculated gridded file of crossover information
with dimensions time, lat, lon (and realizations optional). To make this
file, run wrap_calc_exceedance_grid followed by
wrap_calc_crossover_grid.

The calculation itself addresses the question: "How many samples in one
period exceed the samples from a baseline?" As configured here, this
compares samples from the CESM2 Large Ensemble to the CESM2
Preindustrial control. Beyond a certain threshold, this denotes
"crossover" from a climate state of internal climate variability to one
where the warming trend dominantes.

The underlying statistic is a non-parametric expression of effect size.
The choice of threshold in the definition involves subjectivity. See the
accompanying manuscript Hueholt et al. "Exploring the Influence of
Internal Climate Variability and Forced Change on Arctic Greening" for
discussion of the values used by default here.
"""
import sys
import warnings

import cmasher as cmr
from icecream import ic
from matplotlib import colormaps as cm
import matplotlib.pyplot as plt
import numpy as np
import rioxarray
import xarray as xr

import classes_gddt as cg
import fun_calc_var as fcv
import fun_plots as fpl
import fun_process as fproc

def crossover_map(
    bmb_type, cmn_path, paint_shapefile_bool, dp_crossover, dp_altmax,
    setp_crossover, setp_altmax, ppar_crossover, ppar_altmask, ppar_shapefile):
    """ Plot crossover map
    Arguments:
    bmb_type -- string of biomass burning filename to use
        Valid: "forcingcmip6", "forcingsmoothed", "allmembers"
    cmn_path -- output path for maps
    paint_shapefile_bool -- plot ecoregion shapefile or raw output
    dp_crossover -- DataParams instance for crossover data
    dp_altmax -- DataParams instance for altmax data
    setp_crossover -- SetParams instance for crossover data
    setp_altmax -- SetParams instance for altmax data
    ppar_crossover -- PlotParams instance for crossover data
    ppar_altmask -- PlotParams instance for active layer mask
    ppar_shapefile -- PlotParams instance for ecoregion shapefile
    """
    crossover_dict = fproc.common_opener(dp=dp_crossover, setp=setp_crossover)
    da_crossover = crossover_dict['land_mask'].compute()
    crossover_attrs = crossover_dict['raw_ds'].attrs
    ppar_crossover.o_name = bmb_type + '_' + ppar_crossover.o_name
    if ppar_crossover.plot_crossover_dict['forced_dict']['bool']:
        msg_crossover = 'Plotting forced crossover at ' \
            + str(crossover_attrs['threshold']) + '% threshold'
        ppar_crossover.o_name += 'forcedcrossover_LENS2_threshold' \
            + str(crossover_attrs['threshold'])
        ppar_crossover.title = 'LENS2 forced crossover (ensemble mean >' \
            + str(crossover_attrs['threshold']) + '% of Preindustrial samples)'
    elif ppar_crossover.plot_crossover_dict['member_dict']['bool']:
        msg_crossover = 'Plotting member crossover at quantile ' \
            + str(ppar_crossover.quantile) + ' with ' \
            + str(crossover_attrs['threshold']) + '% threshold'
        ppar_crossover.o_name += 'qoi' + str(ppar_crossover.quantile) \
            + 'membercrossover_LENS2_threshold' \
            + str(crossover_attrs['threshold'])
        ppar_crossover.title = 'LENS2 ' + 'qoi' \
            + str(ppar_crossover.quantile) + ' member crossover (member >' \
            + str(crossover_attrs['threshold']) + '% of Preindustrial samples)'
        if ppar_crossover.quantile is not False:
            da_crossover = da_crossover.quantile(
                ppar_crossover.quantile, dim='realization', skipna=True)
        else:
            raise NotImplementedError(
                'Member crossover only implemented for quantiles.')
    else:
        raise ValueError(
            'Check inputs! Ensure forced OR member crossover are not None.')
    ic(msg_crossover)

    plt.rcParams.update({'font.family': 'Catamaran'})
    #  'font.weight': normal, bold, heavy, light, ultrabold, ultralight
    plt.rcParams.update({'font.weight': 'normal'})
    plt.rcParams.update({'font.size': 12})
    if paint_shapefile_bool:
        ppar_crossover.o_name = ppar_crossover.o_prefix + 'sf_' \
            + ppar_crossover.o_name
        ic(ppar_crossover.o_name, ppar_crossover.title)
        geo_arctic, geo_other = fproc.get_arctic_biomes()
        #  EPSG 4326: World Geodetic | EPSG 3995: North Polar Stereographic
        geo_arctic_proj = geo_arctic.to_crs(epsg=4326)
        no_plot_proj = geo_other.to_crs(epsg=3995)
        da_crossover.coords['lon'] = (
            da_crossover.coords['lon'] + 180) % 360 - 180
        da_crossover = da_crossover.sortby(da_crossover.lon)
        da_crossover.rio.set_spatial_dims(
            x_dim="lon", y_dim="lat", inplace=True)
        da_crossover.rio.write_crs("epsg:4326", inplace=True)
        list_crossover_ecoregion = list()
        for loop_ecoregion in np.arange(0, len(geo_arctic_proj)):
            ecoregion = geo_arctic_proj.iloc[loop_ecoregion]
            da_crossover_ecoregion = da_crossover.rio.clip(
                [ecoregion.geometry], geo_arctic_proj.crs, all_touched=True,
                drop=False)
            lat_weights = np.cos(np.deg2rad(da_crossover_ecoregion.lat))
            da_weighted = da_crossover_ecoregion.weighted(lat_weights)
            da_regionmn = da_weighted.mean(dim=['lat', 'lon'], skipna=True)
            list_crossover_ecoregion.append(da_regionmn.data)
        ppar_crossover.color = list_crossover_ecoregion
        plot_this = geo_arctic_proj.to_crs(epsg=3995)
        fpl.plot_globe_shapefile(plot_this, no_plot_proj, ppar_crossover)
    else:
        ppar_crossover.o_name = ppar_crossover.o_prefix + 'nosf_' \
            + ppar_crossover.o_name
        ic(ppar_crossover.o_name, ppar_crossover.title)
        plot_this = da_crossover
        fig, ax = fpl.plot_globe(plot_this, ppar_crossover)

def make_crossover_map(
    bmb_type='', cmn_path='20260129_crossoverbmb/',
    paint_shapefile_bool=False):
    """ Make crossover map
    Keyword arguments
    bmb_type -- string of biomass burning filename to use
        Valid: "forcingcmip6", "forcingsmoothed", "allmembers"
    cmn_path -- output path for maps with user directory added manually
        in first line of function
    paint_shapefile_bool -- plot ecoregion shapefile or raw output
    """
    cmn_path = '/Users/danielhueholt/Documents/Figures/arc-gdd_fig/' + cmn_path
    dp_crossover = cg.DataParams(
        path='/Users/danielhueholt/Data/gddt_data/LENS2/exceedance/crossover/',
        tok=bmb_type + '_forcedcrossover_threshold80.nc', var='crossover',
        flag_raw_ds=True, flag_raw_da=True, flag_time_slice=False,
        flag_manage_rlz=True, flag_land_mask=True, flag_roi=True)
    #  ALTMAX (active layer depth) used to mask bedrock and permafrost.
    #  This has minimal impact, but it's still useful to note!
    dp_altmax = cg.DataParams(
        path='/Users/danielhueholt/Data/gddt_data/LENS2/annual_ALTMAX/',
        tok='*.nc', var='ALTMAX', flag_raw_ds=False, flag_raw_da=True,
        flag_time_slice=False, flag_manage_rlz=False, flag_land_mask=False,
        flag_roi=False)
    setp_crossover = cg.SetParams(
        area_stat='pass', base_yrs=[0, 2000],
        mask='/Users/danielhueholt/Data/gddt_data/mask/cesm_atm_mask.nc',
        mask_flag='land', reg_oi='global', rlz='all', yrs=[1850, 2100])
    setp_altmax = cg.SetParams(
        area_stat='pass', mask_flag='land', reg_oi='global', rlz='all',
        yrs=[1850, 2100])
    ppar_crossover = cg.PlotParams(
        cb_bool=True, cb_extent='neither', cb_label='crossover year',
        cb_vals=[2000, 2050], cmap=fpl.crossover_n(n=10), dpi=800,
        figsize=(5,4), o_bool=True, o_name='', o_path=cmn_path, o_prefix='',
        plot_crossover_dict=dict(
            forced_dict=dict(bool=True, threshold=(80,)),
            member_dict=dict(bool=False)),
        plot_each_member=False, proj='Arctic', quantile=None, title='',
        title_size=10)
    #  Plot parameters for image muting based on active layer
    ppar_altmask = cg.PlotParams(
        alpha=0.6, cb_bool=False, cb_extent='neither', cb_label='auto',
        cb_vals=[-10, 0], cmap=cm['Greys'], dpi=800, figsize=(5,4),
        o_bool=True, o_name='', o_path=cmn_path, o_prefix='',
        plot_each_member=False, proj='Arctic', set_bad=False, title='',
        title_size=11)
    #  Plot parameters for shapefile. Match cmap here to ppar_crossover
    ppar_shapefile = cg.PlotParams(
        cmap=cmr.get_sub_cmap('cmr.torch_r', 0.15, 0.85, N=10), o_prefix='',
        title='Arctic ecoregions by crossover year', title_size=12)
    crossover_map(
        bmb_type, cmn_path, paint_shapefile_bool, dp_crossover, dp_altmax,
        setp_crossover, setp_altmax, ppar_crossover, ppar_altmask, ppar_shapefile)

def make_noanalog_map(
    bmb_type='', cmn_path='20260129_crossoverbmb/',
    paint_shapefile_bool=False):
    ''' Make no-analog map
    Keyword arguments
    bmb_type -- string of biomass burning filename to use
        Valid: "forcingcmip6", "forcingsmoothed", "allmembers"
    cmn_path -- output path for maps with user directory added manually
        in first line of function
    paint_shapefile_bool -- plot ecoregion shapefile or raw output
    '''
    cmn_path = '/Users/danielhueholt/Documents/Figures/arc-gdd_fig/' + cmn_path
    dp_crossover = cg.DataParams(
        path='/Users/danielhueholt/Data/gddt_data/LENS2/exceedance/crossover/',
        tok=bmb_type + '_membercrossover_threshold100.nc', var='crossover',
        flag_raw_ds=True, flag_raw_da=True, flag_time_slice=False,
        flag_manage_rlz=True, flag_land_mask=True, flag_roi=True)
    dp_altmax = cg.DataParams(
        path='/Users/danielhueholt/Data/gddt_data/LENS2/annual_ALTMAX/',
        tok='*.nc', var='ALTMAX', flag_raw_ds=False, flag_raw_da=True,
        flag_time_slice=False, flag_manage_rlz=False,
        flag_land_mask=False, flag_roi=False)
    setp_crossover = cg.SetParams(
        area_stat='pass', base_yrs=[0, 2000],
        mask='/Users/danielhueholt/Data/gddt_data/mask/cesm_atm_mask.nc',
        mask_flag='land', reg_oi='global', rlz='all', yrs=[1850, 2100])
    setp_altmax = cg.SetParams(
        area_stat='pass', mask_flag='land', reg_oi='global', rlz='all',
        yrs=[1850, 2100])
    ppar_crossover = cg.PlotParams(
        cb_bool=True, cb_extent='neither', cb_label='crossover year',
        cb_vals=[2050, 2100],
        cmap=cmr.get_sub_cmap('cmr.torch_r', 0.15, 0.85, N=10), dpi=800,
        figsize=(5,4), o_bool=True, o_name='', o_path=cmn_path, o_prefix='',
        plot_crossover_dict=dict(
            forced_dict=dict(bool=False, ),
            member_dict=dict(bool=True, threshold=(100,))),
        plot_each_member=False, proj='Arctic', quantile=0.5, title='',
        title_size=10)
    ppar_altmask = cg.PlotParams(
        alpha=0.6, cb_bool=False, cb_extent='neither', cb_label='auto',
        cb_vals=[-10, 0], cmap=cm['Greys'], dpi=800, figsize=(5,4), o_bool=True,
        o_name='', o_path=cmn_path, o_prefix='', plot_each_member=False,
        proj='Arctic', set_bad=False, title='', title_size=11)
    ppar_shapefile = cg.PlotParams(
        cmap=cmr.get_sub_cmap('cmr.torch_r', 0.15, 0.85, N=10), o_prefix='',
        title='Arctic ecoregions by crossover year', title_size=12)
    crossover_map(
        bmb_type, cmn_path, paint_shapefile_bool, dp_crossover, dp_altmax,
        setp_crossover, setp_altmax, ppar_crossover, ppar_altmask, ppar_shapefile)
