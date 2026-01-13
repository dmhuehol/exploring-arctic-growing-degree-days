''' fun_ghcn
Helper functions for Global Historical Climatology Network (GHCN)
analysis.

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''
import collections
import sys

from icecream import ic
import matplotlib.path as mpth
import numpy as np
import pandas as pd
import polars as pl
from sklearn import linear_model

import fun_calc_var as fcv



def add_dtns(df_ghcn, dt_col='DATETIME'):
    ''' Add ns datetime to GHCN dataframe

    Arguments:
    df -- dataframe of GHCN data
    dt_col -- name of column with datetime information

    Returns:
    df_ghcn_dtns -- df with "DATETIME" column as ns datetime
    '''
    df_ghcn_dtns = df_ghcn.with_columns(
       pl.col(dt_col).str.strptime(pl.Date).cast(pl.Datetime).dt \
       .cast_time_unit("ns"))

    return df_ghcn_dtns

def get_spancov(scp):
    ''' Obtain stations matching span and coverage thresholds.

    Arguments:
    sc_d -- set dict containing thresholds and spancov file
    '''
    df_sc = pl.read_csv(scp.f)
    df_span = df_sc.filter(df_sc["span_days"] > scp.span_thr)
    if scp.cov_type == '>':
        df_span_cov = df_span.filter(df_span["coverage"] > scp.cov_thr)
    elif scp.cov_type == '<':
        df_span_cov = df_span.filter(df_span["coverage"] < scp.cov_thr)
    elif scp.cov_type == '=':
        df_span_cov = df_span.filter(df_span["coverage"] == scp.cov_thr)
    else:
        raise ValueError('Check cov_type input.')

    return df_span_cov

def make_polygon_mask_ghcn(lat_ghcn, lon_ghcn, lat_reg, lon_reg):
    ''' Make mask for a non-rectangular polygonal region, adapted for
    GHCN inventory DataFrame.
    '''
    flat_latlon = np.transpose(np.vstack((lat_ghcn, lon_ghcn)))
    reg_poly = np.transpose(np.vstack((lat_reg, lon_reg)))
    #  Define region as a path
    reg_path = mpth.Path(reg_poly)
    #  The heavyweight: "given path, find points contained inside"
    ll_mask = reg_path.contains_points(flat_latlon)

    return ll_mask

def manage_area_ghcn(df_inv, reg_oi, tol='closest'):
    ''' Operations to filter by area. Based on my previous methods for
    DataArrays, but not the same because GHCN station data is too large
    to hold in memory as a meshgrid. Plus, it doesn't make sense to
    treat station data as a grid anyway!
    tol is tolerance in summed latlon distance for a point. If
    'closest', returns closest station; otherwise, all within tol value.
    '''
    lat = df_inv['lat']
    lon = df_inv['lon'] % 360

    if reg_oi == 'global':
        df_inv_filter = df_inv
    elif isinstance(reg_oi, dict):  # region_library-style objects
        if len(reg_oi['reg_lats']) == 1:  # Point region_library object
            abs_lat_diff = np.abs(lat - reg_oi['reg_lats'])
            abs_lon_diff = np.abs(lon - reg_oi['reg_lons'])
            abs_latlon_diff = abs_lat_diff + abs_lon_diff
            if tol == 'closest':
                in_tol = np.argmin(abs_latlon_diff)
                closest_msg = 'Closest distance: ' + str(
                    np.round(abs_latlon_diff[in_tol], 3))
                ic(closest_msg)
            else:
                in_tol = abs_latlon_diff < tol
                # in_tol = (abs_lat_diff < tol) & (abs_lon_diff < tol)  # Alternate
            try:  # Works if one index
                df_inv_filter = df_inv.iloc[in_tol]
            except:  # Works if multiple indices (pandas indexing is a bit odd!)
                df_inv_filter = df_inv[in_tol]
            return df_inv_filter  # Shortcut rest of function
        elif len(reg_oi['reg_lons']) > 2:  # Non-rectangular not crossing meridian
            latlon_mask = make_polygon_mask_ghcn(
                lat, lon, reg_oi['reg_lats'], reg_oi['reg_lons'])
            try:  # Works if multiple indices
                df_inv_filter = df_inv[latlon_mask]
            except:  # Works if one index (in case only one station in region)
                true_ind = latlon_mask[latlon_mask == True]
                df_inv_filter = df_inv.iloc[true_ind]  # Indexing spelled out for clarity
        elif isinstance(reg_oi['reg_lons'], tuple):  # Non-rectangular crossing meridian
            latlon_mask_list = list()
            for subc, subv in enumerate(reg_oi['reg_lons']):
                s_latlon_mask = make_polygon_mask_ghcn(
                    lat, lon, reg_oi['reg_lats'][subc], subv)
                latlon_mask_list.append(s_latlon_mask)
            latlon_overall_mask = np.logical_or.reduce(latlon_mask_list)
            try:  # Works if multiple indices
                df_inv_filter = df_inv[latlon_overall_mask]
            except:  # Works if one index (in case only one station in region)
                true_ind = latlon_overall_mask[latlon_overall_mask == True]
                df_inv_filter = df_inv.iloc[true_ind]
        #  Rectangular region
        else:
            lat_mask = (lat > reg_oi['reg_lats'][0]) \
                & (lat < reg_oi['reg_lats'][1])
            #  Rectangle does not cross meridian
            if reg_oi['reg_lons'][0] < reg_oi['reg_lons'][1]:
                lon_mask = (lon > reg_oi['reg_lons'][0]) \
                    & (lon < reg_oi['reg_lons'][1])
            #  Rectangle crosses meridian
            else:
                lon_mask = (lon > reg_oi['reg_lons'][0]) \
                    | (lon < reg_oi['reg_lons'][1])
            latlon_mask = lat_mask & lon_mask
            df_inv_filter = df_inv[latlon_mask]
    else:
        sys.exit('Invalid region!')

    return df_inv_filter

def manage_span(df_inv, span_type, span_por):
    ''' Filter GHCN inventory dataframe by span of record '''
    df_span = df_inv['end'] - df_inv['start']
    if span_type == '>':
        df_inv_filter_span = df_inv[df_span > span_por]
    elif span_type == '<':
        ic(df_span)
        df_inv_filter_span = df_inv[df_span < span_por]
        ic(df_inv_filter_span)
    else:
        raise ValueError('Check span_type input!')

    return df_inv_filter_span

def manage_var(df_inv, var_oi):
    ''' Filter GHCN inventory dataframe by variable of interest. If
    multiple variables are input, returns all stations recording either
    variable.
    '''
    df_var = df_inv['var']
    if np.size(var_oi) > 1: # np.size avoids str input issues
        var_mask_list = list()
        for vc, vv in enumerate(var_oi):
            ind_var_mask = df_var == vv
            var_mask_list.append(ind_var_mask)
        var_mask = np.logical_or.reduce(var_mask_list)
    else:
        var_mask = df_var == var_oi
    df_inv_filter_var = df_inv[var_mask]

    return df_inv_filter_var

def open_inv(inv_path, inv_name):
    ''' Open inventory file as dataframe '''
    series_block = pd.read_table(inv_path + inv_name).squeeze()
    regexp = '(?P<id>[^"]{11})' + r'\s+' + r'(?P<lat>-*\d+\.\d{4})' + r'\s+' \
        + r'(?P<lon>-*\d+\.\d{4})' + r'\s' + '(?P<var>[^"]{4})' + r'\s' \
        + r'(?P<start>\d{4})' + r'\s' + r'(?P<end>\d{4})'
    df_inv = series_block.str.extract(regexp)
    df_inv['lat'] = np.round(df_inv['lat'].astype("double"), decimals=4)
    df_inv['lon'] = np.round(df_inv['lon'].astype("double"), decimals=4)
    df_inv['start'] = df_inv['start'].astype("int")
    df_inv['end'] = df_inv['end'].astype("int")

    return df_inv


def report_in_lens2(df_indist, df_total):
    ''' Report stations within LENS2. This comes up in several different
   places so is best as a function, not a separate analysis script.
   Arguments:
       df_indist -- DataFrame within distribution
       df_total -- DataFrame with all data
    '''
    num_indist = len(df_indist)
    num_total = len(df_total.index)
    percent_indist = round(num_indist / num_total * 100, 2)
    msg_in_lens2 = str(percent_indist) + '% of obs (' \
        + str(num_indist) + ' out of ' + str(num_total) + ' stations)' \
        + ' within distribution of LENS2 in interval: '

    return msg_in_lens2


def trend_ghcn(ldf_act, scp, t_slc):
    ''' Calculate trend information for a GHCN station.
    Arguments:
        ldf_act -- a lazy DataFrame of GHCN data
        scp -- SpanCovParams instance
        t_slc -- trend dictionary with years of interest
    '''
    ldf_dtns = add_dtns(ldf_act, dt_col='DATETIME')
    ldf_slice = ldf_dtns.filter(
        pl.col("DATETIME").is_between(t_slc[0], t_slc[1]))
    ldf_gdd_only = ldf_slice.filter(pl.col('GDD_SUM').is_not_nan())
    lgdds = ldf_gdd_only.select(pl.len()).collect().item()
    if lgdds >= scp.por:
        dts = ldf_gdd_only.select('DATETIME').collect().to_series().dt.year().to_numpy()
        gdds = ldf_gdd_only.select('GDD_SUM').collect().to_series().to_numpy()
        d_trend = fcv.calc_lin_reg(dts, gdds)
    else:
        d_trend = {
            #  "data" key inherited on successful fcv.calc_lin_reg
            "grad": None,
            "intcpt": None
        }

    return d_trend
