'''wrap_crossover_sensitivity
Plot sensitivity of crossover to the choice of threshold value.
'''
import sys

import cmasher as cmr
from icecream import ic
from matplotlib import colormaps as cm
import matplotlib.pyplot as plt
import numpy as np

import classes_gddt as cg
import fun_calc_var as fcv
import fun_plots as fpl
import fun_process as fproc
import gddt_region_library as g_rlib

def sensitivity_map(bmb_type, delta, cmn_path, paint_shapefile_bool):
    ''' Plot sensitivity map
    Inputs
    bmb_type -- string of biomass burning filename to use
    delta -- delta value to perturb from 80 for sensitivity test
    cmn_path -- output path for maps
    paint_shapefil_bool -- plot ecoregion shapefile or raw output
    '''
    dp_crossover_max = cg.DataParams(
        path='/Users/danielhueholt/Data/gddt_data/LENS2/exceedance/crossover/',
        tok=bmb_type + '_forcedcrossover_threshold' + str(80 + delta) + '.nc',
        var='crossover',
        flag_raw_ds=True, flag_raw_da=True, flag_time_slice=False,
        flag_manage_rlz=True, flag_land_mask=True, flag_roi=True)
    dp_crossover_min = cg.DataParams(
        path='/Users/danielhueholt/Data/gddt_data/LENS2/exceedance/crossover/',
        tok=bmb_type + '_forcedcrossover_threshold' + str(80 - delta) + '.nc',
        var='crossover',
        flag_raw_ds=True, flag_raw_da=True, flag_time_slice=False,
        flag_manage_rlz=True, flag_land_mask=True, flag_roi=True)
    setp_crossover = cg.SetParams(
        area_stat='pass', base_yrs=[0, 2000],
        mask='/Users/danielhueholt/Data/gddt_data/mask/cesm_atm_mask.nc',
        mask_flag='land', reg_oi='global', rlz='all', yrs=[1850, 2100])
    setp_altmax = cg.SetParams(
        area_stat='pass', mask_flag='land', reg_oi='global', rlz='all',
        yrs=[1850, 2100])
    #  Plot parameters for base layer (crossover year)
    #  Control which kind of crossover is plotted by using forced_dict and
    #  member_dict within ppar_crossover
    #  Forced crossover:
    #      cb_vals=[2000, 2050]
    #      cmap=fpl.crossover_n(n=10)
    #      forced_dict=dict(bool=True, threshold=(80,))
    #      member_dict=dict(bool=False,)
    #      quantile=None
    #  Median no-analog state:
    #      cb_vals=[2050, 2100]
    #      cmap=cmr.get_sub_cmap('cmr.torch_r', 0.15, 0.85, N=10)
    #      forced_dict=dict(bool=False, )
    #      member_dict=dict(bool=True, threshold=(100,))
    #      quantile=0.5
    ppar_crossover = cg.PlotParams(
        cb_bool=True, cb_extent='neither', cb_label='crossover year',
        cb_vals=[0, 40],
        # cmap=fpl.crossover_n(n=10), dpi=800,
        cmap=cmr.get_sub_cmap('cmr.viola', start=0, stop=1, N=21), dpi=800,
        figsize=(5,4), o_bool=True, o_name='', o_path=cmn_path, o_prefix='',
        plot_crossover_dict=dict(
            forced_dict=dict(bool=True, threshold=(80,)),
            member_dict=dict(bool=False)),
        plot_each_member=False, proj='Arctic', quantile=None, title='',
        title_size=10)
    #  Plot parameters for image muting based on active layer
    ppar_altmask = cg.PlotParams(
        alpha=0.6, cb_bool=False, cb_extent='neither', cb_label='auto',
        cb_vals=[-10, 0], cmap=cm['Greys'], dpi=800, figsize=(5,4), o_bool=True,
        o_name='', o_path=cmn_path, o_prefix='', plot_each_member=False,
        proj='Arctic', set_bad=False, title='', title_size=11)
    #  Plot parameters for shapefile if using shapefile
    #  Match cmap here to whichever is in use above in ppar_crossover
    ppar_shapefile = cg.PlotParams(
        cmap=cmr.get_sub_cmap('cmr.torch_r', 0.15, 0.85, N=10), o_prefix='',
        title='Arctic ecoregions by crossover year', title_size=12)

    crossover_max_dict = fproc.common_opener(dp=dp_crossover_max, setp=setp_crossover)
    crossover_min_dict = fproc.common_opener(dp=dp_crossover_min, setp=setp_crossover)
    da_crossover_max = crossover_max_dict['land_mask'].compute()
    crossover_attrs_max = crossover_max_dict['raw_ds'].attrs
    da_crossover_min = crossover_min_dict['land_mask'].compute()
    crossover_attrs_min = crossover_min_dict['raw_ds'].attrs
    da_sensitivity = da_crossover_max - da_crossover_min
    fcv.check_stats(da_sensitivity)
    if ppar_crossover.plot_crossover_dict['forced_dict']['bool']:
        msg_crossover = 'Plotting sensitivity between ' \
            + str(crossover_attrs_max['threshold']) + '% and ' \
            + str(crossover_attrs_min['threshold']) + '% threshold'
        if paint_shapefile_bool:
            ppar_crossover.o_name = 'sf_' + bmb_type + '_' \
                + 'forcedcrossover_LENS2_threshold' \
                + str(crossover_attrs_max['threshold']) + '-' \
                + str(crossover_attrs_min['threshold'])
        else:
            ppar_crossover.o_name = 'nosf_' + bmb_type + '_' \
                + 'forcedcrossover_LENS2_threshold' \
                + str(crossover_attrs_max['threshold']) + '-' \
                + str(crossover_attrs_min['threshold'])
        ppar_crossover.title = 'LENS2 forced crossover sensitivity between ' \
            + str(crossover_attrs_max['threshold']) + '% and ' \
            + str(crossover_attrs_min['threshold']) + '% threshold'
    elif ppar_crossover.plot_crossover_dict['member_dict']['bool']:
        raise NotImplementedError('Sensitivity for forced only')
        # msg_crossover = 'Plotting member crossover at quantile ' \
        #     + str(ppar_crossover.quantile) + ' with ' \
        #     + str(crossover_attrs['threshold']) + '% threshold'
        # ppar_crossover.o_name = 'qoi' + str(ppar_crossover.quantile) \
        #     + 'membercrossover_LENS2_threshold' \
        #     + str(crossover_attrs['threshold'])
        # ppar_crossover.title = 'LENS2 ' + 'qoi' + str(ppar_crossover.quantile) \
        #     + ' member crossover (member >' + str(crossover_attrs['threshold']) \
        #     + '% of Preindustrial samples'
        # if ppar_crossover.quantile is not False:
        #     da_crossover = da_crossover.quantile(
        #         ppar_crossover.quantile, dim='realization', skipna=True)
        # else:
        #     raise NotImplementedError(
        #         'Member crossover only implemented for quantiles.')
    else:
        raise ValueError(
            'Check inputs! Ensure forced OR member crossover are not None.')
    ic(msg_crossover)
    name_dict = fproc.namer(crossover_max_dict["raw_da"], setp_crossover)

    ppar_crossover.o_name = ppar_crossover.o_prefix + ppar_crossover.o_name
    ic(ppar_crossover.o_name, ppar_crossover.title)
    plt.rcParams.update({'font.family': 'Catamaran'})
    #  'font.weight': normal, bold, heavy, light, ultrabold, ultralight
    plt.rcParams.update({'font.weight': 'normal'})
    plt.rcParams.update({'font.size': 12})
    if paint_shapefile_bool:
        geo_arctic, geo_other = fproc.get_arctic_biomes()
        #  EPSG 4326: World Geodetic | EPSG 3995: North Polar Stereographic
        geo_arctic_proj = geo_arctic.to_crs(epsg=4326)
        no_plot_proj = geo_other.to_crs(epsg=3995)
        da_sensitivity.coords['lon'] = (da_sensitivity.coords['lon'] + 180) % 360 - 180
        da_sensitivity = da_sensitivity.sortby(da_sensitivity.lon)
        da_sensitivity.rio.set_spatial_dims(x_dim="lon", y_dim="lat", inplace=True)
        da_sensitivity.rio.write_crs("epsg:4326", inplace=True)
        list_crossover_ecoregion = list()
        for loop_ecoregion in np.arange(0, len(geo_arctic_proj)):
            ecoregion = geo_arctic_proj.iloc[loop_ecoregion]
            da_sensitivity_ecoregion = da_sensitivity.rio.clip(
                [ecoregion.geometry], geo_arctic_proj.crs, all_touched=True,
                drop=False)
            lat_weights = np.cos(np.deg2rad(da_sensitivity_ecoregion.lat))
            da_weighted = da_sensitivity_ecoregion.weighted(lat_weights)
            da_regionmn = da_weighted.mean(dim=['lat', 'lon'], skipna=True)
            list_crossover_ecoregion.append(da_regionmn.data)
        ppar_crossover.color = list_crossover_ecoregion
        plot_this = geo_arctic_proj.to_crs(epsg=3995)
        fpl.plot_globe_shapefile(plot_this, no_plot_proj, ppar_crossover)
    else:
        plot_this = da_sensitivity
        fig, ax = fpl.plot_globe(plot_this, ppar_crossover)


cmn_path = '/Users/danielhueholt/Documents/Figures/arc-gdd_fig/20260128_crossover/'
l_paint_shapefile_bool = [True, False]
#  Crossover information calculated from wrap_calc_exceedance_grid
#  followed by wrap_calc_crossover_grid
#  DataParams for forced crossover at 80% threshold
l_bmb_type = ['forcingcmip6', 'forcingsmoothed', 'allmembers']
l_delta = [1, 5, 10]

for delta in l_delta:
    for bmb_type in l_bmb_type:
        for paint_shapefile_bool in l_paint_shapefile_bool:
            sensitivity_map(bmb_type, delta, cmn_path, paint_shapefile_bool)
