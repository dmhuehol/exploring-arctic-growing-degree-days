'''fun_calc_var
Calculate variables of interest.

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''
import sys

from icecream import ic
import numpy as np
from sklearn import linear_model
import xarray as xr


def calc_crossover(data_exceedance, threshold):
    ''' Calculate crossover at a given threshold given exceedance
    data. '''
    if threshold == 100:
        all_crossover = data_exceedance == threshold
    elif (threshold < 100) & (threshold > 0):
        all_crossover = data_exceedance > threshold
    else:
        ic(threshold)
        raise ValueError('Threshold must be in range [0, 100].')
    all_crossover = np.array(all_crossover)

    return all_crossover

def calc_lin_reg(x_data, y_data):
    ''' Calculate linear regression and retrieve coefficients '''
    est = linear_model.LinearRegression(fit_intercept=True)
    x_2d = x_data.reshape(-1, 1)
    est.fit(x_2d, y_data)
    temporal_grad = {
        "grad": est.coef_[0],
        "intcpt": est.intercept_,
        "data": est.coef_[0] * x_data + est.intercept_,
    }

    return temporal_grad

def calc_lin_reg_vec(np_x, da_y):
    ''' Calculate linear regression and retrieve coefficients in a
    vectorized way. This saves massive amounts of time for
    multidimensional data, such as large ensembles.
    Derived from implementation at:
    hrishichandanpurkar.blogspot.com/2017/09/vectorized-functions-for-correlation.html

    Arguments:
    np_x -- 1-D numpy array of independent variable
    da_y -- DataArray of dependent variable with 'time' or 'year' dimension

    Returns:
    temporal_grad -- dict containing grad and intercept for each rlz
    '''
    time_entries = len(np_x)
    x_time_mn = np_x.mean()
    x_time_stdev = np_x.std()
    try:
        y_time_mn = da_y.mean(dim='year')
        y_time_stdev = da_y.std(dim='year')
    except ValueError:
        y_time_mn = da_y.mean(dim='time')
        y_time_stdev = da_y.std(dim='time')
    cov = np.sum(
        (np_x - x_time_mn) * (da_y - y_time_mn), axis=0) / (time_entries)
    reg_slp = cov / (x_time_stdev ** 2)
    reg_int = y_time_mn - x_time_mn * reg_slp
    temporal_grad = {
        "grad": reg_slp,
        "intcpt": reg_int
    }

    return temporal_grad


def calc_ransac_reg(x_data, y_data):
    ''' Calculate RANSAC regression and retrieve coefficients '''
    est = linear_model.RANSACRegressor(random_state=13)
    x_2d = x_data.reshape(-1, 1)
    y_rav = y_data.ravel()
    #  Returns model because coefficients are difficult to access
    est.fit(x_2d, y_rav)

    return est

# def calc_theilsen_reg(x_data, y_data):
#     ''' Calculate Theil-Sen regression and retrieve coefficients '''
#     est = linear_model.TheilSenRegressor(random_state=13)
#     x_2d = x_data.reshape(-1, 1)
#     y_rav = y_data.ravel()
#     est.fit(x_2d, y_rav)
#     temporal_grad = {
#         "grad": clf.coef_[0],
#         "intcpt": clf.intercept_
#     }

#     return temporal_grad

def check_stats(in_data):
    ''' Print common statistics for flat input variable '''
    stats_dict = {
        "mean": np.nanmean(in_data),
        "max": np.nanmax(in_data),
        "min": np.nanmin(in_data),
        "median": np.nanmedian(in_data, ),
        "90Q": np.nanquantile(in_data, 0.9),
        "10Q": np.nanquantile(in_data, 0.1),
    }
    ic(stats_dict)
    return None

def gdd(da_mntemp, gddp):
    ''' Calculate growing degree days above a base for an input daily
    mean temperature DataArray.

    Arguments:
    da_mntemp -- DataArray of daily mean temperature
    gddp -- GddParams instance

    Returns:
    da_gdd_sum -- annual count of growing degree days
    '''
    da_daysoi = da_mntemp.where(da_mntemp > gddp.base)
    da_gdd = da_daysoi - gddp.base
    gb_gdd = da_gdd.groupby("time.year")
    da_gdd_sum = gb_gdd.sum()

    return da_gdd_sum

def ndd(da_mntemp):
    ''' Count days with NaNs per year for input daily mean temperature
    DataArray.

    Arguments:
    da_mntemp -- DataArray of daily mean temperature

    Returns:
    da_ndd_sum -- annual count of NaN days
    '''
    da_nandays = np.isnan(da_mntemp)
    gb_ndd = da_nandays.groupby("time.year")
    da_ndd_sum = gb_ndd.sum()

    return da_ndd_sum
