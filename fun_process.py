''' fun_process
Functions to assist in various processing tasks.

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''
import glob
import sys

import cftime
import geopandas as gp
from icecream import ic
import matplotlib.path as mpth
import numpy as np
import time
from numpy.ma.core import indices
import xarray as xr

import classes_gddt as cg
xr.set_options(keep_attrs=True)


def abrv_var(var_str):
    ''' Return an abbreviated variable name from a netCDF long_name '''
    if ('Sea level pressure' in var_str):
        abv_var = var_str.replace('Sea level pressure', 'SLP')
    elif ('Growing degree day information' in var_str) | ('gdd5' in var_str):
        abv_var = 'GDD'
    else:
        msg = "Unknown variable name, returning original string: " + var_str
        abv_var = var_str

    return abv_var

def apply_mask(in_da, setp):
    ''' Apply land/ocean mask and return DataArray with np.nans used to
    hide masked area.

    Arguments:
    in_da -- input DataArray to be masked
    set_d -- SetParams instance with mask_flag and path to mask specified
    '''
    da_mask_data = in_da
    if ('land' in setp.mask_flag) or ('ocean' in setp.mask_flag):
        ds_mask_org = xr.open_dataset(setp.mask)
        #  Compensate for round-off errors introduced by CDO processing
        ds_mask = ds_mask_org.reindex_like(
            in_da, method='nearest', tolerance=10 ** -13)
        try:
            active_mask = ds_mask.landmask
        except AttributeError:
            active_mask = ds_mask.imask
        if 'land' in setp.mask_flag:
            da_mask_data = in_da.where(active_mask > 0)
        elif 'ocean' in setp.mask_flag:
            da_mask_data = in_da.where(active_mask == 0)
        else:
            ic('WARNING: ocean zeros may affect output if CDO mask was used.')
            da_mask_data = in_da

    return da_mask_data

# def area_weight_stat(in_da, set_d):
#     ''' Return area-weighted mean/sum/etc of entire dataset. More
#     complex region handling will be addressed elsewhere. '''
#     lat_wght = np.cos(np.deg2rad(in_da['lat']))
#     da_wght = in_da.weighted(lat_wght)
#     if set_d['area_stat'] == 'mean':
#         da_stat = da_wght.mean(dim=['lat', 'lon'], skipna=True)
#     elif set_d['area_stat'] == 'sum':
#         da_stat = da_wght.sum(dim=['lat', 'lon'], skipna=True)
#     else:
#         ic("No area-weighted stat applied")
#         da_stat = in_da

#     return da_stat

def calc_anomaly(da, da_full, setp):
    ''' Calculate anomalies by removing mean. Optionally divide by
    standard deviation to obtain z-score.

    Arguments:
    da: DataArray of interest
    da_full: DataArray of base period to remove
    setp: SetParams instance

    Returns:
    da_anom: DataArray of anomalies relative to da_full
    '''
    da_mean = da_full.mean(dim=setp.dims)
    da_stdev = da_full.std(dim=setp.dims)
    da_anom = (da - da_mean)
    da_anom.attrs['long_name'] = da_anom.attrs['long_name'] + ' anomaly'
    if setp.z_flag:
        da_anom = da_anom / da_stdev
        da_anom.attrs['units'] = 'z'

    return da_anom

def calc_trend(da):
    ''' Calculate linear trend using polyfit '''
    da_pf = da.polyfit('year', deg=1)
    slope = da_pf.polyfit_coefficients.sel(degree=1)
    da_guide = slope.squeeze()

    return da_guide


def cliffs_delta(np_x1, np_x2):
    ''' Calculate Cliff's Delta statistic for two samples of realizations. '''
    np_x1_col = np_x1[:, None]
    x1gtx2 = (np_x1_col > np_x2) * 1
    x1ltx2 = (np_x1_col < np_x2) * 1
    count_x1gtx2 = np.count_nonzero(x1gtx2)
    count_x1ltx2 = np.count_nonzero(x1ltx2)
    n1 = len(np_x1)
    try:
        n2 = len(np_x2)
    except TypeError:
        #  For cliffs of ensemble mean
        n2 = 1
    cliff = (count_x1ltx2 - count_x1gtx2) / (n1 * n2)

    return cliff

def common_opener(dp=cg.DataParams(), setp=cg.SetParams()):
    ''' All open functions aggregated into one for ease of use,
    with options to run processing and retain data as requested in
    DataParams instance.

    Keyword arguments:
    dp -- DataParams instance
    setp -- SetParams instance

    Returns:
    open_d: dict with Datasets and DataArrays flagged in dp
    '''
    open_d = dict(
        raw_ds=None, raw_da=None, time_slice=None, manage_rlz=None,
        land_mask=None, flag_roi=None)
    try:
        ds_in = xr.open_mfdataset(
            dp.path + dp.tok, concat_dim='realization', combine='nested',
            chunks={'time': 10000}, coords='minimal')
    except TypeError:
        ds_in = xr.open_mfdataset(
            dp.tok, concat_dim='realization', combine='nested',
            chunks={'time': 10000}, coords='minimal')
    da_in = ds_in[dp.var].squeeze()
    if dp.flag_raw_ds:
        open_d['raw_ds'] = ds_in
    if dp.flag_raw_da:
        open_d['raw_da'] = da_in

    da_toi = da_in
    if dp.flag_time_slice:
        try:
            #  Standard CESM output format
            time_slice = slice(
                cftime.DatetimeNoLeap(setp.yrs[0], 1, 1, 0, 0, 0, 0),
                cftime.DatetimeNoLeap(setp.yrs[1], 12, 31, 14, 24, 0, 0)
            )
            da_toi = da_in.sel(time=time_slice)
        except KeyError:
            #  Processed growing degree days format
            time_slice = slice(setp.yrs[0], setp.yrs[1])
            da_toi = da_in.sel(year=time_slice)
        open_d['time_slice'] = da_toi

    da_rlz = da_toi
    if dp.flag_manage_rlz:
        da_rlz = manage_rlz(da_toi, setp)
        open_d['manage_rlz'] = da_rlz

    da_mask = da_rlz
    if dp.flag_land_mask:
        da_mask = apply_mask(da_rlz, setp)
        open_d['land_mask'] = da_mask

    da_roi = da_mask
    if dp.flag_roi:
        da_roi, loc_str, _ = manage_area(da_mask, setp)
        open_d['roi'] = da_roi
        open_d['loc_str'] = loc_str

    return open_d

def common_calc_exceed(
        da_calc_exceed, da_base_period, years_to_calc, ppar, setp):
    ''' Calculate exceedance and apply moving average. This feels a bit
    conceptually overloaded--but it's a use case that several scripts
    need and modularizing it prevents each script from calling a
    slightly different implementation of this code!
    '''
    np_base_samples = np.ravel(da_base_period.data)
    #  Note that "exceedance" is synonymous with Gexc in robustness from
    #  Hueholt et al. 2022. I think the term "exceedance" is clearer here;
    #  in this calculation the subceedance term is not used so there's no
    #  need for nomenclature that includes both.
    list_exceed = list()
    for loop_yr in years_to_calc:
        da_loop_yr = da_calc_exceed.sel(year=loop_yr)
        np_compare_to_base = da_loop_yr.data
        dict_g = exceed_subceed(np_base_samples, np_compare_to_base)
        list_exceed.append(np.array(dict_g['gexc']))
    if ppar.plot_as_percent:
        np_exceed = np.array(list_exceed) / len(np_base_samples) * 100.
    else:
        np_exceed = np.array(list_exceed)

    if setp.window is not None:
        l_exceed_rolling_avg = list()
        try:
            for rlz in np.arange(0, len(da_calc_exceed.realization)):
                rolling_avg = moving_average(
                    np_exceed[:, rlz], setp.window)
                l_exceed_rolling_avg.append(rolling_avg)
        except AttributeError:
            l_exceed_rolling_avg.append(moving_average(np_exceed, setp.window))
    else:
        l_exceed_rolling_avg = None

    return np_exceed, l_exceed_rolling_avg

def count_threshold(in_data, threshold, type='>', stat=''):
    ''' Count data points past a threshold, return count or percent '''
    if type == '>':
        mask_data = in_data > threshold
    elif type == '<':
        mask_data = in_data < threshold
    else:
        raise ValueError('type must be ">" or "<"')
    count_data = np.count_nonzero(mask_data)
    if stat == 'count':
        out_stat = count_data
    elif stat == 'fraction':
        out_stat = count_data / len(in_data)
    elif stat == 'percent':
        out_stat = count_data / len(in_data) * 100
    else:
        raise ValueError('stat must be "count", "prob", or "percent"')

    return out_stat

def check_in_dist(dist_list, val):
    ''' Check if value is in list. Useful for checking a value against
    a distribution of realizations.

    Arguments:
    dist_list -- list distribution, e.g., large ensemble realizations
    val -- value to check against

    Returns:
    check_in_dist -- bool of whether val is within dist_list bounds
    '''
    dist_max = max(dist_list)
    dist_min = min(dist_list)
    check_in_dist = (val < dist_max) & (val > dist_min)

    return check_in_dist


def get_arctic_biomes(
        sf="/Users/danielhueholt/Data/gddt_data/Ecoregions2017/Ecoregions2017.shp"):
    geo_sf = gp.read_file(sf)
    arctic_biomes = ["Boreal Forests/Taiga", "Tundra", "Rock and Ice"]
    geo_arctic = geo_sf[geo_sf['BIOME_NAME'].isin(arctic_biomes)]
    geo_other = geo_sf[~geo_sf['BIOME_NAME'].isin(arctic_biomes)]

    return geo_arctic, geo_other


def get_params(type='', bmb_type='allmembers', cmn_path=''):
    ''' Get DataParams and cmn_path information based on location '''
    error_msg = 'Obsolete but retained (relevant for processing).'
    raise NotImplementedError(error_msg)
    match type:
        case 'local':
            dp_gdd = cg.DataParams(
                path='/Users/danielhueholt/Data/gddt_data/gdd/lens2/',
                tok='*arc*.nc', var='gdd5_sum', flag_raw_ds=True,
                flag_raw_da=True, flag_time_slice=True, flag_manage_rlz=True,
                flag_land_mask=False, flag_roi=True)
            dp_gdd_roi_alltimes = cg.DataParams(
                path='/Users/danielhueholt/Data/gddt_data/gdd/lens2/',
                tok='*arc*.nc', var='gdd5_sum', flag_raw_ds=False,
                flag_raw_da=False, flag_time_slice=False,
                flag_manage_rlz=False, flag_land_mask=False, flag_roi=True)
            dp_psl = cg.DataParams(
                path='/Users/danielhueholt/Data/gddt_data/LENS2/merge_PSL/AMJJAS/',
                tok='*.nc', var='PSL', flag_raw_ds=True, flag_raw_da=True,
                flag_time_slice=True, flag_manage_rlz=True,
                flag_land_mask=False, flag_roi=False)
            dp_sst = cg.DataParams(
                path='/Users/danielhueholt/Data/gddt_data/LENS2/monthly_SST/AMJJAS/',
                tok=bmb_type + '*.nc', var='SST', flag_raw_ds=True, flag_raw_da=True,
                flag_time_slice=True, flag_manage_rlz=True,
                flag_land_mask=False, flag_roi=False)
            dp_icefrac = None
            if cmn_path == '':
                cmn_path = '/Users/danielhueholt/Documents/Figures/arc-gdd_fig/'
                cmn_path += '20260202_composite-crossover/'
        case 'coe_hpc':
            dp_gdd = cg.DataParams(
                path='/barnes-engr-scratch1/DATA/CESM2-LE/processed_data/annual/gdd/reproc_20250218/',
                tok='*arc*.nc', var='gdd5_sum', flag_raw_ds=True,
                flag_raw_da=True, flag_time_slice=True, flag_manage_rlz=True,
                flag_land_mask=False, flag_roi=True)
            dp_gdd_roi_alltimes = cg.DataParams(
                path='/barnes-engr-scratch1/DATA/CESM2-LE/processed_data/annual/gdd/reproc_20250218/',
                tok='*arc*.nc', var='gdd5_sum', flag_raw_ds=False,
                flag_raw_da=False, flag_time_slice=False,
                flag_manage_rlz=False, flag_land_mask=False, flag_roi=True)
            dp_psl = cg.DataParams(
                path='/barnes-engr-scratch1/DATA/CESM2-LE/processed_data/seasonal/AMJJAS/PSL/',
                tok='*.nc', var='PSL', flag_raw_ds=True, flag_raw_da=True,
                flag_time_slice=True, flag_manage_rlz=True,
                flag_land_mask=False, flag_roi=False)
            dp_sst = cg.DataParams(
                path='/barnes-engr-scratch1/DATA/CESM2-LE/processed_data/seasonal/AMJJAS/SST/',
                tok='*.nc', var='SST', flag_raw_ds=True, flag_raw_da=True,
                flag_time_slice=True, flag_manage_rlz=True,
                flag_land_mask=False, flag_roi=False)
            dp_icefrac = cg.DataParams(
                path='/barnes-engr-scratch1/DATA/CESM2-LE/processed_data/seasonal/AMJJAS/ICEFRAC/',
                tok='*.nc', var='ICEFRAC', flag_raw_ds=True, flag_raw_da=True,
                flag_time_slice=True, flag_manage_rlz=True,
                flag_land_mask=False, flag_roi=False)
            if cmn_path == '':
                cmn_path = '/home/dhueholt/gddt_fig/frames/for_manuscript_202502/'

    return dp_gdd, dp_gdd_roi_alltimes, dp_psl, dp_sst, dp_icefrac, cmn_path


def get_season(da_toi, set_d):
    ''' Obtain a season of interest from a DataArray

    Arguments:
    da_toi -- input DataArray for season extraction
    set_d -- set dict containing months of interest
    '''
    raise NotImplementedError("Season management not yet implemented.")

    return None

def guide(da_ens, gp):
    ''' Guide (select indices corresponding to characteristic) an
    ensemble DataArray. prep_guide should be run first to prepare the
    DataArray.

    Arguments:
    da_ens -- ensemble DataArray, characteristic selected by prep_guide
    gp -- GuideParams object
    '''
    ind_dict = dict()
    np_tmn = da_ens.data
    da_max = da_ens.max(dim='realization').data
    da_min = da_ens.min(dim='realization').data
    ind_sort = np.argsort(np_tmn)
    np_sort = da_ens.sel(realization=ind_sort)
    allq = np.arange(0, 1 + gp.qoi, gp.qoi)
    ind_q = dict()
    for qc, qv in enumerate(allq):
        q_ind = (np_sort <= np.quantile(np_sort, qv)) \
            & (np_sort >= np.quantile(np_sort, allq[qc-1]))
        ind_q[str(qv)] = ind_sort[q_ind]

    ind_positive = ind_sort[da_ens.data[ind_sort] > 0]
    ind_negative = ind_sort[da_ens.data[ind_sort] < 0]
    if "max" in gp.composite_key:
        try:
            noi = int(gp.composite_key.replace('max', ''))
        except ValueError:
            msg = "Assuming max 1 entry."
            ic(msg)
            noi = 1
        maxn = ind_sort[-noi:]
        max_key = "max" + str(noi)
        ind_dict[max_key] = maxn
    elif "min" in gp.composite_key:
        try:
            noi = int(gp.composite_key.replace('min', ''))
        except ValueError:
            msg = "Assuming min 1 entry."
            ic(msg)
            noi = 1
        minn = ind_sort[:noi]
        min_key = "min" + str(noi)
        ind_dict[min_key] = minn
    ind_dict["all_sort"] = ind_sort
    ind_dict["quantiles"] = ind_q
    ind_dict["positive"] = ind_positive
    ind_dict["negative"] = ind_negative

    return ind_dict

def manage_area(darr, setp):
    ''' Manage area operations: obtain global, regional, or pointal output.
    This was copied from SAI-ESM and adapted for style but not all
    functionality has been tested. There may be hidden syntax errors.

    Arguments:
    darr -- DataArray for area selection
    setp -- SetParams instance

    Returns:
    darr: DataArray with selected region
    '''
    if setp.reg_oi == 'global':
        loc_str = 'global'
        loc_title = 'global'
        if setp.area_stat == 'mean':
            lat_wght = np.cos(np.deg2rad(darr['lat']))
            da_wght = darr.weighted(lat_wght)
            darr = da_wght.mean(dim=['lat', 'lon'], skipna=True)
        elif setp.area_stat == 'sum':
            lat_wght = np.cos(np.deg2rad(darr['lat']))
            da_wght = darr.weighted(lat_wght)
            darr = da_wght.sum(dim=['lat', 'lon'], skipna=True)

    elif isinstance(setp.reg_oi, dict): #region_library-style objects
        if len(setp.reg_oi['reg_lats']) == 1: # Point region_library object
            darr = darr.sel(
                lat=setp.reg_oi['reg_lats'], lon=setp.reg_oi['reg_lons'],
                method="nearest")
            ic(darr['lat'].data, darr['lon'].data) #Lat/lon to check 'nearest'
            darr = np.squeeze(darr) #Drop length 1 lat/lon dimensions
            setp.area_stat = None
            loc_str = setp.reg_oi['reg_abv']
            loc_title = setp.reg_oi['reg_str']
            return darr, loc_str, loc_title # Shortcut the rest of the function
        loc_str = setp.reg_oi['reg_abv']
        loc_title = setp.reg_oi['reg_abv']
        lats = darr['lat']
        lons = darr['lon']
        if len(setp.reg_oi['reg_lons'])>2: #non-rectangular region that does not cross Prime Meridian
            grid_mask = make_polygon_mask(
                lats, lons, setp.reg_oi['reg_lats'], setp.reg_oi['reg_lons'])
            darr_mask = darr.where(grid_mask)
        elif isinstance(setp.reg_oi['reg_lons'], tuple): #non-rectangular region that crosses Prime Meridian
            s_grid_mask_list = list()
            for sc in np.arange(0, len(setp.reg_oi['reg_lons'])):
                s_grid_mask = make_polygon_mask(
                    lats, lons, setp.reg_oi['reg_lats'][sc], setp.reg_oi['reg_lons'][sc])
                s_grid_mask_list.append(s_grid_mask)
            grid_mask = np.logical_or.reduce(s_grid_mask_list)
            darr_mask = darr.where(grid_mask)
        else: #rectangular region
            lat_mask = (lats > setp.reg_oi['reg_lats'][0]) & (lats < setp.reg_oi['reg_lats'][1])
            if setp.reg_oi['reg_lons'][0] < setp.reg_oi['reg_lons'][1]: #rectangle does not cross Prime Meridian
                lon_mask = (lons > setp.reg_oi['reg_lons'][0]) & (lons < setp.reg_oi['reg_lons'][1])
            else: #rectangle crosses Prime Meridian
                lon_mask = (lons > setp.reg_oi['reg_lons'][0]) | (lons < setp.reg_oi['reg_lons'][1])
            lats_oi = lats[lat_mask]
            lons_oi = lons[lon_mask]
            darr_mask = darr.sel(lat=lats_oi, lon=lons_oi)
        if setp.area_stat == 'mean':
            lat_wght = np.cos(np.deg2rad(darr_mask['lat']))
            darr_wght = darr_mask.weighted(lat_wght)
            darr = darr_wght.mean(dim=['lat','lon'], skipna=True)
            # darr = darr_mask.mean(dim=['lat','lon'], skipna=True) #No latitude weighting
        elif setp.area_stat == 'sum':
            darr = darr_mask.sum(dim=['lat','lon'], skipna=True)
        elif setp.area_stat is False:
            darr = darr_mask

    else:
        sys.exit('Invalid region!')
    if setp.area_stat == 'sum':
        ic('NaNs in summations result in 0. Re-implement previous workarounds.')
        # try: #spaghetti
        #     cursedZeros = (darr == 0)
        #     darr[cursedZeros] = np.nan
        # except: #spread
        #     for rc in np.arange(0,len(darr['realization'])):
        #         cursedZeros = (darr.isel(realization=rc) == 0)
        #         darr.isel(realization=rc)[cursedZeros] = np.nan

    return darr, loc_str, loc_title


def manage_rlz(da_mod, setp):
    ''' Filter realizations based on input 'rlz' property in SetParams
    instance. Valid 'rlz' values are:
        'all' -- all members
        'all_withmean' -- all members + mean attached as extra member
        'mean' -- only the ensemble mean
        'pass' -- does nothing, return original object
        integer -- selected member (WARNING: ZERO INDEXED)
        list of integers -- selected members (WARNING: ZERO INDEXED)
    If given a list of members, the realizations are bound to the
    DataArray attributes for easy access.

    Arguments:
    da_mod -- Ensemble DataArray or Dataset for filtering
    setp -- SetParams instance, where 'rlz' attribute controls behavior

    Returns:
    da_rlzoi -- Filtered DataArray or Dataset
    '''
    if setp.rlz == 'mean':
        da_rlzoi = da_mod.mean(dim='realization')
    elif (isinstance(setp.rlz, int)) | (isinstance(setp.rlz, list)):
        da_rlzoi = da_mod.isel(realization=setp.rlz)
        da_rlzoi.attrs['members'] = setp.rlz
    elif setp.rlz == 'all_withmean':
        da_mean = da_mod.mean(dim='realization')
        da_rlzoi = xr.concat((da_mod, da_mean), dim='realization')
    elif (setp.rlz == 'pass') | (setp.rlz == 'all'):
        da_rlzoi = da_mod
    else:
        ic('Invalid input for rlz!')

    return da_rlzoi

def make_polygon_mask(lats, lons, reg_lats, reg_lons):
    ''' Make mask for a non-rectangular polygonal region '''
    grid_lon, grid_lat = np.meshgrid(lons, lats) #lonxlat grids of lon/lat
    flat_lon = np.ravel(grid_lon)
    flat_lat = np.ravel(grid_lat)
    flat_latlon = np.transpose(np.vstack((flat_lat, flat_lon))) #Nx2
    reg_poly = np.transpose(np.vstack((reg_lats, reg_lons))) #Polyx2

    reg_path = mpth.Path(reg_poly) #defines the path of the region
    flat_mask = reg_path.contains_points(flat_latlon) #the heavyweight--calculates whether points from the grid are exterior to the path or not
    grid_mask = flat_mask.reshape((len(lats),len(lons)))

    return grid_mask

def match_rlz_extremes(data_rlz, extreme=''):
    ''' Match realizations to absolute extremes (e.g., 'max10' or 'min10'
    indices) '''
    if 'max' in extreme:
        max_n = int(extreme.replace('max',''))
        indices_extreme = np.argsort(data_rlz)[-max_n:]
    elif 'min' in extreme:
        min_n = int(extreme.replace('min',''))
        indices_extreme = np.argsort(data_rlz)[:min_n]
    else:
        raise ValueError('Invalid entry for extreme!')
    msg_extreme = 'Indices matching ' + extreme + ': ' + str(indices_extreme)
    ic(msg_extreme)

    return indices_extreme


def match_rlz_quantiles(data_rlz, quantile, type='equal'):
    ''' Match realizations to quantile (e.g., identify storylines) '''
    if np.isnan(quantile):
        members_quantile = np.isnan(data_rlz) * 1
        indices_quantile = np.nonzero(members_quantile)
    else:
        data_quantile = np.nanquantile(data_rlz, quantile)
        if type == 'equal':
            members_quantile = data_rlz == np.round(data_quantile)
        elif type == 'leq':
            members_quantile = data_rlz <= np.round(data_quantile)
        elif type == 'geq':
            members_quantile = data_rlz >= np.round(data_quantile)
        else:
            raise ValueError('Improper value for type input')
        indices_quantile = np.squeeze(np.nonzero(members_quantile))
    msg_quantile = 'Indices matching ' + str(quantile) + ' quantile: ' \
        + str(indices_quantile)
    ic(msg_quantile)

    return indices_quantile

def moving_average(x, w):
    ''' From https://stackoverflow.com/a/54628145 '''
    return np.convolve(x, np.ones(w), 'valid') / w

def namer(da_plot, setp):
    ''' Return info useful for subsequent naming as well as a suggested
    filename and title corresponding to input DataArray and set dict.

    Arguments:
    da_plot: DataArray of interest
    setp: SetParams instance corresponding to da_plot

    Returns:
    name_dict: dict containing bits useful for filenames, titles, etc.
    '''
    try:
        var_w = abrv_var(da_plot.long_name)
        var_nw = abrv_var(da_plot.long_name).replace(' ', '')
    except AttributeError:
        var_w = abrv_var(da_plot.name)
        var_nw = var_w
    #  Obviously, this is a placeholder that won't work if future datasets are
    #  added! But I'm trying out my "YAGNI" era here :)
    data_id = "LENS2"
    yr_str = str_yrs(setp.yrs)
    rlz_str = str(setp.rlz)
    try:
        reg_str = setp.reg_oi['reg_str']
        reg_abv = setp.reg_oi['reg_abv']
    except KeyError:
        reg_str = ''
        reg_abv = ''
    except TypeError:
        reg_str = 'global'
        reg_abv = 'global'
    build_name = data_id + '_' + var_nw + '_' + reg_abv + '_' + yr_str \
        + '_' + 'rlz' + rlz_str
    build_name_norlz = data_id + '_' + var_nw + '_' + reg_abv + '_' + yr_str
    build_title = data_id + ': ' + var_w + ' ' + yr_str + ' ' + 'rlz ' + \
        rlz_str
    build_title_norlz = data_id + ': ' + var_w + ' ' + yr_str
    name_dict = {
        #  "a_name" & "a_title" useful if input is only dataset plotted
        "a_name": build_name,
        "a_name_norlz": build_name_norlz,
        "a_title": build_title,
        "a_title_norlz": build_title_norlz,
        #  Info below allows recombination for more complex plots
        "var_w": var_w,
        "var_nw": var_nw,
        "data_id": data_id,
        "yr_str": yr_str,
        "reg_str": reg_str,
        "reg_abv": reg_abv,
    }

    return name_dict

# def pieces(fname):
#     ''' Expand filename into bits for plotting '''
#     ex = dict()
#     f_pieces = fname.split('.')[:-1][0].split('_')
#     ex["d_id"] = f_pieces[1].split('/')[-1]
#     if f_pieces[-1] == 'arc':
#         ex["loc"] = 'Arctic'
#     elif f_pieces[-1] == 'antar':
#         ex["loc"] = 'Antarctic'
#     else:
#         ex["loc"] = 'global'
#     if (ex["d_id"] == 'LE2') | (ex["d_id"] == 'LENS2'):
#         ex["d_long"] = 'LENS2'
#     if "gdd5_sum" in fname:
#         ex["var_long"] = 'annual growing degree days'

#     return ex

def prep_guide(da_to_guide, gp):
    ''' Prep DataArray to generate an ensemble guide by obtaining guide
    characteristic and computing the Dask collection (necessary for
    operations in the guide function).

    Arguments:
    da_to_guide -- DataArray for guiding with data as Dask collection
    gp -- GuideParams object

    Returns:
    da_guide_compute -- computed DataArray of guide characteristic
    '''
    if (gp.guide_by == 'mean' or gp.guide_by == 'mn'):
        da_guide = da_to_guide.mean(dim='year')
    elif gp.guide_by == 'median':
        da_guide = da_to_guide.median(dim='year')
    elif (gp.guide_by == 'stdev' or gp.guide_by == 'std'):
        da_guide = da_to_guide.std(dim='year')
    elif gp.guide_by == 'max':
        da_guide = da_to_guide.max(dim='year')
    elif gp.guide_by == 'min':
        da_guide = da_to_guide.min(dim='year')
    elif 'trend' in gp.guide_by:
        da_guide = calc_trend(da_to_guide)
    else:
        ic('Unknown guide_by entry!')
    da_guide_compute = da_guide.compute()

    return da_guide_compute

def exceed_subceed(np_x1, np_x2):
    ''' Calculate exceedance and subceedance (Gexc and Gsub) for two
    samples of realizations. '''
    #  If the only dimensions are lat/lon. This logic is far from
    #  bulletproof, but it's the best I've come up with for now
    #  without explicit dimension checking using xarray objects.
    if len(np.shape(np_x2)) == 2:
        x1gtx2 = np_x1 > np_x2
        x1ltx2 = np_x1 < np_x2
        count_x1gtx2 = np.count_nonzero(x1gtx2, axis=0)
        count_x1ltx2 = np.count_nonzero(x1ltx2, axis=0)
        dict_g = {
            "gexc": count_x1ltx2,
            "gsub": count_x1gtx2
        }
    #  If there is a realization dimension in np_x2
    else:
        list_above = list()
        list_below = list()
        #  Compares each np_x2 element to every np_x1 element. Sadly, I
        #  have yet to think of a way to do this without a loop!
        for x2c, x2_rlz in enumerate(np_x2):
            # if x2c % 10 == 0:
                # ic(cg.TrackProg(cli=x2c, iter=np_x2).message(prefix='gexc '))
            x1gtx2 = np_x1 > x2_rlz
            x1ltx2 = np_x1 < x2_rlz
            count_x1gtx2 = np.count_nonzero(x1gtx2, axis=0)
            count_x1ltx2 = np.count_nonzero(x1ltx2, axis=0)
            list_above.append(count_x1ltx2)
            list_below.append(count_x1gtx2)
        dict_g = {
            "gexc": list_above,
            "gsub": list_below
        }

    return dict_g

def roll_window(list_windows):
    ''' Roll list of windows into a DataArray '''
    da_windows = xr.concat(list_windows, 'window')

    return da_windows

def ensemble_tokens(path=''):
    """Returns tokens for the biomass burning sub-ensembles of LENS2.
    Keyword Arguments:
    path -- path to merge with tokens, must end in forward slash /
    """
    dict_tokens = {
        "forcingcmip6": [
            '1001-001', '1021-002', '1041-003', '1061-004', '1081-005',
            '1101-006', '1121-007', '1141-008', '1161-009', '1181-010',
            '1231-001', '1231-002', '1231-003', '1231-004', '1231-005',
            '1231-006', '1231-007', '1231-008', '1231-009', '1231-010',
            '1251-001', '1251-002', '1251-003', '1251-004', '1251-005',
            '1251-006', '1251-007', '1251-008', '1251-009', '1251-010',
            '1281-001', '1281-002', '1281-003', '1281-004', '1281-005',
            '1281-006', '1281-007', '1281-008', '1281-009', '1281-010',
            '1301-001', '1301-002', '1301-003', '1301-004', '1301-005',
            '1301-006', '1301-007', '1301-008', '1301-009', '1301-010',],
        "forcingsmoothed": [
            '1011-001', '1031-002', '1051-003', '1071-004', '1091-005',
            '1111-006', '1131-007', '1151-008', '1171-009', '1191-010',
            '1231-011', '1231-012', '1231-013', '1231-014', '1231-015',
            '1231-016', '1231-017', '1231-018', '1231-019', '1231-020',
            '1251-011', '1251-012', '1251-013', '1251-014', '1251-015',
            '1251-016', '1251-017', '1251-018', '1251-019', '1251-020',
            '1281-011', '1281-012', '1281-013', '1281-014', '1281-015',
            '1281-016', '1281-017', '1281-018', '1281-019', '1281-020',
            '1301-011', '1301-012', '1301-013', '1301-014', '1301-015',
            '1301-016', '1301-017', '1301-018', '1301-019', '1301-020',
        ],
        "allmembers": '*.nc',
    }
    dict_tokens['forcingcmip6_asterisk'] = [
        '*' + c6 + '*' for c6 in dict_tokens['forcingcmip6']]
    dict_tokens['forcingcmip6_path'] = [
        glob.glob(path + c6)[0] for c6 in dict_tokens['forcingcmip6_asterisk']]
    dict_tokens['forcingsmoothed_asterisk'] = [
        '*' + sm + '*' for sm in dict_tokens['forcingsmoothed']]
    dict_tokens['forcingsmoothed_path'] = [
        glob.glob(
            path + sm)[0] for sm in dict_tokens['forcingsmoothed_asterisk']]
    dict_tokens['allmembers_asterisk'] = '*.nc'
    dict_tokens['allmembers_path'] = path + '*.nc'

    return dict_tokens

def sync_lengths(dict_both, sync_key='', ref_key=''):
    ''' Repeat values indexed with one key to sync with the length
    of values indexed by another key.

    Arguments:
    dict_both: dictionary with both sync_key and ref_key

    Keyword arguments:
    sync_key: key indexing values to be repeated to ref_key length
    ref_key: keywith intended length

    Returns:
    dict_both: input dict, sync_key values repeated to ref_key length
    '''
    if (len(dict_both[sync_key]) == 1) & (len(dict_both[ref_key]) != 1):
        dict_both[sync_key] = dict_both[sync_key] * len(dict_both[ref_key])

    return dict_both


def str_yrs(yrs_list):
    ''' Return useful string for titles from input years '''
    try:
        if yrs_list[0] - yrs_list[1] == 0:
            yr_str = str(yrs_list[0])
        else:
            yr_str = str(yrs_list[0]) + '-' + str(yrs_list[1])
    except IndexError:
        yr_str = str(yrs_list[0])

    return yr_str

def threshold_count(np_count, beat=None):
    ''' Count above beat number '''
    abv_thresh = np.count_nonzero(np_count > beat, axis=0)

    return abv_thresh
