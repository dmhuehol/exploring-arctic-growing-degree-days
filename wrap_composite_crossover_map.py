'''wrap_composite_crossover_map
Script to composite and plot anomalies in one field corresponding to 
crossover in a guiding variable at a given location.

This requires a pre-calculated gridded file of exceedance information 
with dimensions time, lat, lon (and realizations optional). To make this
file, see wrap_calc_exceedance_grid.
'''
import sys

import cftime
from icecream import ic
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

import classes_gddt as cg
import fun_plots as fpl
import fun_process as fproc
import gddt_region_library as g_rlib

#  Exceedance information calculated from wrap_calc_exceedance_grid
dp_crossover = cg.DataParams(
    path='/Users/dhueholt/Documents/gddt_data/LENS2/exceedance/crossover/', 
    tok='membercrossover_threshold90.nc', var='crossover', 
    flag_raw_ds=True, flag_raw_da=True, flag_time_slice=False, 
    flag_manage_rlz=True, flag_land_mask=False, flag_roi=True)
setp_crossover = cg.SetParams(
    area_stat='pass', base_yrs=[0, 2000], mask_flag='none',
    reg_oi=g_rlib.BrooksRange_colonist(), rlz='all', window=10, 
    yrs=[1850, 2100])
_, _, dp_psl, dp_sst, dp_icefrac, cmn_path = fproc.get_params(
    type='local', cmn_path='')
dp_composite = dp_sst
setp_psl = cg.SetParams(
    dims=['time', 'realization'], mask_flag=None, rlz='all', yrs=[1850, 2100], 
    yrs_rel_to=[], z_flag=True)
setp_sst = cg.SetParams(
    dims=['time', 'realization'], mask_flag=None, match_type='geq', rlz='all', 
    yrs=[1850, 2100], yrs_rel_to=[], z_flag=True)
setp_composite = setp_sst
ppar_composite_crossover = cg.PlotParams(
    cb_bool=True, cb_extent='neither', cb_label='z', 
    cb_ticks=[-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6], cb_vals=[-0.6, 0.6],
    cmap=fpl.balance_n(n=18), dpi=800, figsize=(10, 4), o_bool=True, 
    o_name='LENS2_crossover', o_path=cmn_path, o_prefix='', 
    plot_crossover_dict=dict(
        forced_dict=dict(bool=False,),
        member_dict=dict(bool=True, threshold=(90,),),),
    plot_each_member=False, proj='EqualEarth180', quantile=0.9, title='', 
    title_size=10)

crossover_dict = fproc.common_opener(dp=dp_crossover, setp=setp_crossover)
composite_dict = fproc.common_opener(dp=dp_composite, setp=setp_composite)
da_composite_raw = composite_dict['manage_rlz']
da_crossover_roi = crossover_dict['roi']
np_crossover = da_crossover_roi.data.compute()
indices_at_quantile = fproc.match_rlz_quantiles(
    np_crossover, ppar_composite_crossover.quantile, 
    type=setp_composite.match_type)
list_crossover_scenes = list()
for loop_ind, rlz_ind in enumerate(indices_at_quantile):
    time_slice_crossover = slice(
            cftime.DatetimeNoLeap(
                np_crossover[rlz_ind] - 10, 1, 1, 0, 0, 0, 0),
            cftime.DatetimeNoLeap(
                np_crossover[rlz_ind], 12, 31, 14, 24, 0, 0)
    )
    loop_raw = da_composite_raw.isel(
        realization=rlz_ind).sel(time=time_slice_crossover)
    loop_base = da_composite_raw.sel(time=time_slice_crossover)
    loop_anom = fproc.calc_anomaly(loop_raw, loop_base, setp_composite)
    loop_anom_time_mn = loop_anom.mean(dim='time')
    list_crossover_scenes.append(loop_anom_time_mn)
composite = np.mean(list_crossover_scenes, axis=0)
plot_this = xr.DataArray(
    data=composite,
    dims=["lat", "lon"],
    coords={"lat": da_composite_raw.lat, "lon": da_composite_raw.lon},
    attrs={"units": 'year'}
)
ppar_composite_crossover.o_name += '-quantile' \
    + str(ppar_composite_crossover.quantile).replace('.','') + '_' \
    + crossover_dict['loc_str'] + '_' \
    + str(
        indices_at_quantile).replace('[','').replace(']','').replace(' ', '-')
ppar_composite_crossover.title = 'LENS2 ' + dp_composite.var \
    + ' anomaly for ' + str(int(ppar_composite_crossover.quantile * 100)) \
    + 'th percentile crossovers at ' + setp_crossover.reg_oi['reg_str']

ic(ppar_composite_crossover.o_name, ppar_composite_crossover.title)
plt.rcParams.update({'font.family': 'Catamaran'})
#  'font.weight': normal, bold, heavy, light, ultrabold, ultralight
plt.rcParams.update({'font.weight': 'heavy'})
plt.rcParams.update({'font.size': 12})
fpl.plot_globe(plot_this, ppar_composite_crossover)