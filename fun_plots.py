''' fun_plots 
Plotting functions designed to take the same basic inputs of data,
ideally in Numpy array (non-gridded) or xarray DataArray (gridded) 
format, and a PlotParams class instance.
########################################################################
Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''
from glob import glob
import os
import sys

import cartopy.crs as ctcrs
import cartopy.feature as cfeat
import cartopy.util as cutil
import cmasher
import cmocean
import cv2
from icecream import ic
from matplotlib import colorbar
from matplotlib import colormaps
import matplotlib.colors as mcolors
import matplotlib.font_manager as fm
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
import seaborn as sn
import xarray as xr

import fun_process as fproc
########################################################################
####  GLOBAL
########################################################################
#  All plotting functions need access to all system fonts
font_path_local = '/Users/dhueholt/Library/Fonts/'
font_path_coe_hpc = '/home/dhueholt/fonts/'
if os.path.exists(font_path_local):
    for font in fm.findSystemFonts(font_path_local):
        fm.fontManager.addfont(font)
elif os.path.exists(font_path_coe_hpc):
    for font in fm.findSystemFonts(font_path_coe_hpc):
        fm.fontManager.addfont(font)

########################################################################
####  PLOTTING FUNCTIONS
########################################################################
def plot_globe(da_plot, ppar, ax=None, fig=None):
    ''' Make a globe plot. Derived from DrawOnGlobe in SAI-ESM, originally
    written by Elizabeth Barnes. Arctic and regional maps adapted from 
    code by Ariel Morrison.
    
    Arguments:
    da_plot -- DataArray with data to plot
    ppar -- PlotParams instance
    
    Returns:
    ax -- the active axis object
    Saves PNG figure to path specified in ppar
    '''
    if ax is None:
        fig = plt.figure(figsize=ppar.figsize)
        ax = set_proj(fig, ppar)
    d_crs = ctcrs.PlateCarree()
    if ppar.coastline_bool:
        ax.coastlines(linewidth = 1.2, color="#555555")
    try:
        lons = da_plot.longitude
        lats = da_plot.latitude
    except:
        lons = da_plot.lon
        lats = da_plot.lat
    # da_plot = np.nan_to_num(da_plot, nan=-9999)
    # data_cyc, lons_cyc = cutil.add_cyclic_point(da_plot, coord=lons, axis=-1)
    image = ax.pcolormesh(
                lons, lats, da_plot, transform=d_crs, cmap=ppar.cmap, 
                alpha=ppar.alpha)
###############################################################################
    # image = ax.pcolormesh(
    #             lons_cyc, lats, data_cyc, transform=d_crd, cmap=ppar.cmap, 
    #             alpha=1)
    if ppar.cb_bool is True:
        cb = plt.colorbar(
                    image, shrink=.75, orientation="vertical", pad=.02,
                    ticks=ppar.cb_ticks, extend=ppar.cb_extent)
        cb.ax.tick_params(labelsize=8)
        if ppar.cb_label == 'auto':
            cb.set_label(da_plot.attrs['units'], size="small")
        else:
            cb.set_label(ppar.cb_label, size="small")
    try:
        image.set_clim(ppar.cb_vals[0], ppar.cb_vals[1])
    except ValueError:
        max_val = np.nanmax(np.abs(da_plot.data))
        use_val = max_val * 0.9
        auto_msg = 'Automatic colorbar: +/-' + str(np.round(use_val, 2))
        ic(auto_msg)
        image.set_clim(-1 * use_val, use_val)
    if ppar.set_bad:
        ppar.cmap.set_bad('#cfcfcf', alpha=0)
    try:
        plt.title(ppar.title, fontsize=ppar.title_size)
    except ValueError:
        plt.title(ppar.title)
    if ppar.o_bool:
        plt.savefig(ppar.o_path + ppar.o_name + '.png', dpi=ppar.dpi)
    
    return fig, ax


def plot_globe_ng(data_plot, lats, lons, ppar, ax=None):
    ''' Make a globe plot for non-gridded data. Derived from DrawOnGlobe 
    in SAI-ESM, originally written by Elizabeth Barnes. Arctic and 
    regional maps adapted from code by Ariel Morrison. '''
    if ax is None:
        fig = plt.figure(figsize=ppar.figsize)
        ax = set_proj(fig, ppar)
    if ppar.coastline_bool:
        ax.coastlines(linewidth = 1.2, color="#555555")
    d_crs = ctcrs.PlateCarree()
###############################################################################
    try:
        image = ax.scatter(
            lons, lats, s=ppar.marker_size, c=data_plot, transform=d_crs, 
            cmap=ppar.cmap, alpha=ppar.alpha, marker=ppar.marker, 
            edgecolors=ppar.edgecolors, zorder=100)
    except ValueError:
        image = ax.scatter(
            lons, lats, s=ppar.marker_size, c=ppar.cmap, transform=d_crs, 
            alpha=ppar.alpha, marker=ppar.marker, edgecolors=ppar.edgecolors, 
            zorder=100)
    if ppar.color is None:
        image.set_facecolor("none")
    if ppar.cb_bool is True:
        cb = plt.colorbar(
                    image, shrink=.75, orientation="vertical", pad=.02,
                    ticks=ppar.cb_ticks, extend=ppar.cb_extent)
        cb.ax.tick_params(labelsize=8)
        if ppar.cb_label == 'auto':
            cb.set_label(data_plot.attrs['units'], size="small")
        else:
            cb.set_label(ppar.cb_label, size="small")
    try:
        image.set_clim(ppar.cb_vals[0], ppar.cb_vals[1])
    except ValueError:
        max_val = np.nanmax(np.abs(data_plot.data))
        use_val = max_val * 0.9
        auto_msg = 'Colorbar set automatically: +/-' + str(
            np.round(use_val, 2))
        ic(auto_msg)
        image.set_clim(-1 * use_val, use_val)
    # ppar.cmap.set_bad('#cfcfcf', 1.)
    try:
        ax.set_facecolor(ppar.ax_facecolor)
    except ValueError:
        pass
    try:
        plt.title(ppar.title, fontsize=ppar.title_size)
    except ValueError:
        plt.title(ppar.title)
    if ppar.o_bool:
        plt.savefig(ppar.o_path + ppar.o_name + '.png', dpi=400)
    
    return ax

def plot_globe_shapefile(geo_plot, geo_noplot, ppar):
    ''' Plot globe from shapefile '''
    crs = ctcrs.NorthPolarStereo()
    fig, ax = plt.subplots(subplot_kw=dict(projection=crs))
    ax.set_extent([-180, 180, 49, 90], crs=ctcrs.PlateCarree())
    norm = plt.Normalize(ppar.cb_vals[0], ppar.cb_vals[1])
    for rc, geom in enumerate(geo_plot['geometry']):
        acolor = tuple(ppar.cmap(norm(ppar.color[rc])))
        ax.add_geometries(
            geom, crs=crs, facecolor=acolor, edgecolor='#d3d3d3', linewidth=0.2)
    ax.add_geometries(geo_noplot['geometry'], crs=crs, facecolor='#d3d3d3')
    ax.add_feature(cfeat.OCEAN, facecolor='#a9bdc3', edgecolor='none')
    if ppar.cb_bool is True:
        ax_cb = fig.add_axes([0.85, 0.2, 0.02, 0.6])
        cb = colorbar.ColorbarBase(ax_cb, cmap=ppar.cmap, norm=norm, spacing='proportional')
        if ppar.cb_label == 'auto':
            cb.set_label(geo_plot.attrs['units'], size="small")
        else:
            cb.set_label(ppar.cb_label, size="medium")
    ax.set_title(ppar.title, fontsize=ppar.title_size)
    plt.savefig(ppar.o_path + ppar.o_name + '.png', dpi=ppar.dpi)
    
    
def plot_hist(lens2_to_plot, comp_to_plot=None, ppar=None):
    ''' Plot histogram of arbitrary LENS2 data and optional comparison data '''
    fig, ax = plt.subplots()
    if ppar.bw is None:
        bw_da = lens2_to_plot[list(lens2_to_plot.keys())[0]]
        # iqr = stats.iqr(bw_da, nan_policy='omit')
        q75 = np.percentile(bw_da, 75)
        q25 = np.percentile(bw_da, 25)
        iqr = q75 - q25
        bw = 2 * iqr / (len(bw_da) ** (1/3)) # Freedman-Diaconis rule
    else:
        bw = ppar.bw # Manual
    ic(bw)
    colors = ('#296111', '#337915', '#3d9119', '#48aa1e', '#52c222', '#5cda26')
    for perc, per in enumerate(lens2_to_plot.keys()):
        act_dist = np.ravel(lens2_to_plot[per].data)
        ax = sn.histplot(
            data=act_dist, label=per, color=colors[perc % 6], 
            edgecolor=None, stat=ppar.stat, kde=ppar.kde_bool, binwidth=bw)
        if ppar.mn_bool:
            mn_lens2_areas = np.mean(act_dist)
            plt.plot([mn_lens2_areas, mn_lens2_areas], [0, 1000], color=colors[perc % 6], 
                linestyle=ppar.linestyle, linewidth=ppar.lw)
    if comp_to_plot is not None:
        if ppar.comp_hist_bool:
            ax = sn.histplot(
                data=comp_to_plot, label=ppar.comp_label, color=ppar.color_comp, 
                edgecolor=None, stat=ppar.stat, kde=ppar.kde_bool, binwidth=ppar.bw)
        else:
            plt.plot([comp_to_plot, comp_to_plot], [0, 1000], color=ppar.color_comp,
                linewidth=ppar.lw, label=ppar.comp_label)
    # h = 10
    # plt.plot([2, 2], [0, h], lw=2, color='#0cc0aa', label='CESM2 shallowest bedrock')
    # plt.plot([4500, 4500], [0, h], lw=2, color='#f849b6')
    # plt.plot([4750, 4750], [0, h], lw=2, color='#ecac5c')
    # plt.plot([5000, 5000], [0, h], lw=2, color='#830c6f')
    
    plt.annotate(
        "bw: " + str(np.round(bw, 2)), [0.1, 0.9], xytext=[0.13, 0.85],
        xycoords='figure fraction', fontsize=8)
    apply_params_hist(ppar)
    # plt.legend(
    #     fontsize=8, bbox_to_anchor=(0.9, -0.12), 
    #     ncol=len(lens2_to_plot.keys()))
    # plt.tight_layout()
    # if ppar.o_bool:
    #     if ppar.dpi == 'pdf':
    #         plt.savefig(ppar.o_path + ppar.o_name + '.pdf')
    #     else:
    #         plt.savefig(ppar.o_path + ppar.o_name + '.png', dpi=ppar.dpi)
    
    return None
    
    
def plot_matrix(arr, set_plot):
    ''' Plot a matrix, i.e. of pattern correlations '''
    fig, ax = plt.subplots()
    image = plt.imshow(arr, cmap=set_plot['cmap'])
    # for (r, c), val in np.ndenumerate(arr):
        # ax.text(c, r, '{:0.2f}'.format(val), ha='center', va='center')
    if set_plot["cb_bool"] is True:
        cb = plt.colorbar(
                    image, shrink=.75, orientation="vertical", pad=.02,
                    ticks=set_plot['cb_ticks'], extend=set_plot['cb_extent'])
        cb.ax.tick_params(labelsize=8)
        cb.set_label(set_plot['cb_label'], size="small")
    image.set_clim(set_plot['cb_vals'][0], set_plot['cb_vals'][1])
    plt.title(set_plot["title"], fontsize=11)
    plt.savefig(set_plot['o_path'] + set_plot['o_name'], dpi=400)
    
    return None

def plot_timeseries(to_plt, x_d=None, ppar=None):
    ''' Make a simple timeseries of one variable.
    
    Keyword arguments:
    to_plt -- object to plot, ideally an np array
    x_d -- abscissa data
    ppar -- PlotParams instance
    
    Returns: None at present.
    '''
    if x_d is None:
        #  Cut the nonsense and just make a no-frills plot! :)
        plt.plot(to_plt, lw=2, color='k')
    else:
        plt.plot(
            x_d, to_plt, lw=ppar.lw, color=ppar.color, label=ppar.label)

    if ppar is not None:
        if ppar.x_lim == 'auto':
            ppar.x_lim = [x_d[0], x_d[-1]]
        apply_params_ts(ppar=ppar)
    else:
        if ppar.o_flag:
            plt.savefig('auto_out.pdf')
            
    return None

def plot_effect_timeseries(to_plt, x_d=None, lera_gm=None, ppar=None):
    ''' Make a timeseries of an effect size statistic, with optional
    thresholds denoting "crossover."
    
    Keyword arguments:
    to_plt -- object to plot, ideally an np array
    x_d -- abscissa data
    ppar -- PlotParams instance
    
    Returns: None at present.
    '''
    if not ppar.storyline:
        plt.plot(x_d, to_plt, lw=ppar.lw, color=ppar.color_r, label=ppar.label)
        ## HARD CODE
        #db9dbe
        # plt.plot(x_d, to_plt[:, 24], lw=0.9, color='#97215C', label=ppar.label)
        # plt.plot(x_d, to_plt[:, 52], lw=0.9, color='#b66363', label=ppar.label)
        plt.plot(x_d, to_plt[:, 11], lw=0.9, color='k', label=ppar.label)
        plt.plot(x_d, to_plt[:, 52], lw=0.9, color='k', label=ppar.label)
    if ppar.mn_bool:
        to_plt_mn = np.mean(to_plt, axis=1)
        ic(to_plt_mn)
        to_plt_lera = lera_gm / 2000 * 100
        plt.plot(x_d, to_plt_mn, lw=3, color=ppar.color, label='ens mean')
        # plt.plot(x_d, to_plt_lera, lw=3, color='k', label='ens mean')
    #  Manually plot custom vertical line
    if ppar.forced_crossover_bool:
        # cross_thresh = (60, 80, 90, 95, 99, 99.9, 100)
        cross_thresh = (20, 40, 60, 80, 100)
        # cross_thresh = (80, )
        for thr in cross_thresh:
            if thr != 80:
                plt.plot(
                    [1850, 2100], [thr, thr], color='#2098ae',
                    linestyle='--', lw=0.5)
            if thr == 80:
                plt.plot(
                    [1850, 2100], [thr, thr], color='#2098ae',
                    linestyle='--', lw=1)
                crossover = to_plt_mn > thr
                x_crossover = x_d[crossover][0]
                plt.plot(
                    [x_crossover, x_crossover], [-100, 100], color='#2098ae',
                    linestyle='--', lw=1)
            # if thr == 100:
            #     crossover = to_plt_mn == thr
            # else:
            #     crossover = to_plt_mn > thr
            # try:
            #     x_crossover = x_d[crossover][0]
            #     plt.plot(
            #         [x_crossover, x_crossover], [-100, 100], color='#2098ae',
            #         linestyle='--')
            #     str_crossover = 'Crossover ' + str(thr) + '%' + \
            #         ' in ens mean: ' + str(x_crossover)
            #     ppar.title = ppar.title #+ ' ' + str_crossover
            #     ic(str_crossover)
            # except IndexError:
            #     str_crossover = 'No crossover ' + str(thr) + '%' + ' in ens mean'
            # ic(str_crossover)
    if ppar.member_crossover_bool:
        cross_thresh = 90
        if cross_thresh == 100:
            all_crossover = to_plt == cross_thresh
        else:
            all_crossover = to_plt > cross_thresh
        member_crossover = list()
        for ac in all_crossover.T:
            try:
                member_crossover.append(x_d[ac][0])
            except IndexError:
                member_crossover.append(np.nan)
        member_crossover = np.array(member_crossover)
        ic(to_plt)
        ic(member_crossover, member_crossover[11], member_crossover[52])

        if not ppar.storyline:
            mc_10p = np.nanquantile(member_crossover, 0.1)
            mc_90p = np.nanquantile(member_crossover, 0.9)
            ic(mc_10p, mc_90p)
            # plt.plot(
            #     [1850, 2100], [cross_thresh, cross_thresh], color='#d285ae',
            #     linestyle='--')
            plt.fill_between([mc_10p, mc_90p], [-100, -100], [100, 100], color='#d285ae', edgecolor=None, alpha=0.6)
        else:
            if ppar.storyline < 1:
                mc_sp = np.nanquantile(member_crossover, ppar.storyline)
                members_p = member_crossover == np.round(mc_sp)
                ic(mc_sp, members_p.nonzero())
                members_to_plt = to_plt[:, members_p]
                plt.plot(x_d, members_to_plt, lw=ppar.lw, color=ppar.color_r, label=ppar.label)
                plt.plot(
                    [1850, 2100], [cross_thresh, cross_thresh], color='#d285ae',
                    linestyle='--')
            else:
                members_to_plt = to_plt[:, ppar.storyline]
                mc_sp = x_d[members_to_plt > cross_thresh][0]
                ic(mc_sp)
                plt.plot(x_d, members_to_plt, lw=ppar.lw, color=ppar.color_r, label=ppar.label)
                plt.plot(
                    [1850, 2100], [cross_thresh, cross_thresh], color='#d285ae',
                    linestyle='--')
                plt.plot(
                    [np.round(mc_sp), np.round(mc_sp)], [0, 100.05], color='#d285ae',
                    linestyle='--')
                ppar.title = ppar.title + ' crossover: ' + str(np.round(mc_sp))

    plt.gcf().axes[0].spines['top'].set_visible(False)
    plt.gcf().axes[0].spines['right'].set_visible(False)
    if ppar is not None:
        if ppar.x_lim == 'auto':
            ppar.x_lim = [x_d[0], x_d[-1]]
        apply_params_ts(ppar=ppar)
    else:
        if ppar.o_flag:
            plt.savefig('auto_out.pdf')
            
    return None

def plot_timeseries_spaghetti(to_plot, x_d, ppar):
    ''' Make a "spaghetti" timeseries of one variable. Each realization
    is displayed as a separate curve.

    Arguments:
    to_plt -- object to plot, ideally an np array
    x_d -- abscissa data
    ppar -- PlotParams instance
    
    Returns: None at present.
    '''
    for rlz in to_plot.realization:
        act_rlz = to_plot.sel(realization=rlz)
        if rlz == to_plot.realization[-1]:
            plt.plot(x_d, act_rlz, lw=1, color=ppar.color_r, label=ppar.label)
        else:
            plt.plot(x_d, act_rlz, lw=1, color=ppar.color_r)
    if ppar.mn_bool:
        rlz_mn = to_plot.mean(dim='realization', skipna=True)
        plt.plot(x_d, rlz_mn, lw=2, color=ppar.color, label='mean')
    #  Manually add a vertical line
    # plt.plot([1949, 1949], [0, 2000], lw=0.5, color='k')
    if ppar.x_lim == 'auto':
        ppar.x_lim = [x_d[0], x_d[-1]]
    #  Applies params AND saves file
    apply_params_ts(ppar=ppar)
            
    return None
            
def plot_timeseries_spread(da_to_plot, x_data, ppar):
    ''' Make a timeseries of one variable with ensemble variability
    shown as spread. '''
    rlz_max = da_to_plot.max(dim='realization', skipna=True)
    rlz_min = da_to_plot.min(dim='realization', skipna=True)
    rlz_mn = da_to_plot.mean(dim='realization', skipna=True)
    plt.fill_between(
        x_data, rlz_max.data, rlz_min.data, color=ppar.color_r, 
        alpha=0.2, linewidth=0)
    if ppar.mn_bool:
        plt.plot(x_data, rlz_mn, linewidth=2, color=ppar.color)
    ppar.o_name = ppar.o_name.replace('ts', 'ts-spread')
    apply_params_ts(ppar=ppar)
    
    return None
    
            
def images_mp4(
    frame_path, frame_token, mp4_path, mp4_name, fps=1, w=None, h=None):
    ''' Make an mp4 animation from a folder of images '''
    png_frames = sorted(glob(frame_path + frame_token))
    if not png_frames:
        #  Provide clear error message in this common case
        raise FileNotFoundError("No images matching: " + frame_path + frame_token)
    pngs_dur = [(pf, fps) for pf in png_frames]
    for png, dur in pngs_dur:
        frame = cv2.imread(png)
        if w is None:
            h, w, _ = frame.shape
            codec = cv2.VideoWriter_fourcc(*'avc1')
            writer = cv2.VideoWriter(mp4_path + mp4_name, codec, fps, (w, h))
        for repeat in range(round(dur * 1/fps)):
            writer.write(frame)
    writer.release()
    
    return None



########################################################################
####  HELPER FUNCTIONS
########################################################################
########################################################################
###############################################################################
def apply_params_hist(ppar):
    ''' Apply PlotParams for a histogram plot. Specifically, sets the
    title, x-label, y-label, y-limits, y-ticks, legend visibility, and
    saves figure.
    
    Arguments:
    ppar -- PlotParams instance
    
    Returns: 
    None value
    Saves PNG or PDF to path specified in PlotParams if o_bool True
    '''
    plt.title(ppar.title, fontsize=14)
    plt.xlabel(ppar.x_label)
    plt.ylabel(ppar.y_label)
    if ppar.x_lim != 'auto':
        plt.xlim(ppar.x_lim)
    try:
        plt.xticks(ppar.yticks)
    except Exception as e:
        if 'Failed to convert' in e.__str__():
            plt.xticks()
        else:
            sys.exit('Unknown error! ' + e.__str__())
    if ppar.y_lim != 'auto':
        plt.ylim(ppar.y_lim)
    try:
        plt.yticks(ppar.yticks)
    except Exception as e:
        if 'Failed to convert' in e.__str__():
            plt.yticks()
        else:
            sys.exit('Unknown error! ' + e.__str__())
    if ppar.leg_bool:
        plt.legend()
    if ppar.o_bool:
        if ppar.dpi == 'pdf':
            plt.savefig(ppar.o_path + ppar.o_name + '.pdf')
        else:
            plt.savefig(ppar.o_path + ppar.o_name + '.png', dpi=ppar.dpi)
    
    return None


def apply_params_ts(ppar):
    ''' Apply PlotParams for a timeseries plot. Specifically, sets the
    title, x-label, y-label, y-limits, y-ticks, legend visibility, and
    saves figure. 
    
    Arguments:
    ppar -- PlotParams instance
    
    Returns: 
    None value
    Saves PNG or PDF to path specified in PlotParams if o_bool True
    '''
    plt.title(ppar.title, fontsize=14)
    plt.xlabel(ppar.x_label)
    plt.ylabel(ppar.y_label)
    plt.xlim(ppar.x_lim)
    if ppar.y_lim != 'auto':
        if ppar.y_lim == 'fix_nonnegative':
            y_check = plt.ylim()
            if y_check[0] < 0:
                plt.gca().set_ylim(bottom=0)
        else:
            plt.ylim(ppar.y_lim)
    try:
        plt.xticks(ppar.xticks)
    except Exception as e:
        if 'Failed to convert' in e.__str__():
            plt.xticks()
        else:
            sys.exit('Unknown error! ' + e.__str__())
    try:
        plt.yticks(ppar.yticks)
    except Exception as e:
        if 'Failed to convert' in e.__str__():
            plt.yticks()
        else:
            sys.exit('Unknown error! ' + e.__str__())
    if ppar.leg_bool:
        plt.legend()
    if ppar.o_bool:
        if ppar.dpi == 'pdf':
            plt.savefig(ppar.o_path + ppar.o_name + '.pdf')
        else:
            plt.savefig(ppar.o_path + ppar.o_name + '.png', dpi=ppar.dpi)
    
    return None
    
def movie_maker(ppar, intvls, scp, yr_str):
    spn_str = str(intvls[0][0]) + str(intvls[-1][1])
    ppar.anim_d["mp4_path"] = ppar.o_path
    ppar.anim_d["mp4_name"] = ppar.o_name.replace(
        yr_str, spn_str + '_yr' + scp.por_str + '.mp4')
    images_mp4(
        ppar.o_path, ppar.anim_d["frame_tok"], ppar.anim_d["mp4_path"], 
        ppar.anim_d["mp4_name"])
        
def mask_region(region):
    ''' Gets mask for region on a lat/lon grid '''
    lats = np.arange(-90, 91, 1)
    lons = np.arange(0, 360, 1)

    #  Non-rectangular region that does not cross Prime Meridian
    if len(region['reg_lons']) > 2:
        grid_mask = fproc.make_polygon_mask(
            lats, lons, region['reg_lats'], region['reg_lons'])
        lats_to_plot = lats
        lons_to_plot = lons
    #  Non-rectangular region that crosses Prime Meridian
    elif isinstance(region['reg_lons'], tuple):
        l_sgrid_mask = list()
        for sc in np.arange(0, len(region['reg_lons'])):
                sgrid_mask = fproc.make_polygon_mask(
                    lats, lons, region['reg_lats'][sc], region['reg_lons'][sc])
                l_sgrid_mask.append(sgrid_mask)
        grid_mask = np.logical_or.reduce(l_sgrid_mask)
        lats_to_plot = lats
        lons_to_plot = lons
    #  Rectangular regions
    else:
        lat_mask = (lats > region['reg_lats'][0]) & (lats < region['reg_lats'][1])
        #  Rectangle does not cross Prime Meridian
        if region['reg_lons'][0] < region['reg_lons'][1]:
            lon_mask = (lons > region['reg_lons'][0]) & (lons < region['reg_lons'][1])
        #  Rectangle crosses the Prime Meridian
        else:
            lon_mask = (lons > region['reg_lons'][0]) | (lons < region['reg_lons'][1])
        #  Rectangular regions can just grab lat/lon of interest directly
        lats_to_plot = lats[lat_mask]
        lons_to_plot = lons[lon_mask]

    plot_ones = np.ones((len(lats_to_plot), len(lons_to_plot)))
    try:
        #  Non-rectangular regions must retain background grid and mask outside of region
        plot_ones[~grid_mask] = np.nan
    except Exception as e:
        ic(e.__str__())
        pass  
    da_plot_ones = xr.DataArray(
        data=plot_ones,
        dims=["lat", "lon"],
        coords=dict(
            lat=lats_to_plot,
            lon=lons_to_plot
        )
    )

    return da_plot_ones
    

def set_proj(fig, ppar):
    ''' Set up projection and gridlines for a map plot '''
    if isinstance(ppar.proj, str):
        if ppar.proj.lower() == 'arctic':
            ax = fig.add_subplot(
                1, 1, 1, projection=ctcrs.NorthPolarStereo())
            ax.set_extent([-180, 180, 49, 90], crs=ctcrs.PlateCarree())
            gl = ax.gridlines(crs=ctcrs.PlateCarree(), draw_labels=False,
                            linewidth=0.8, color='k', alpha=0.8, linestyle=':',
                            x_inline=False, y_inline=True, rotate_labels=False)
            gl.xlocator = mticker.FixedLocator([-180, -120, -60, 0, 60, 120])
            gl.ylocator = mticker.FixedLocator([60, 70, 80])
            plt.draw()  # Enable the use of `gl._labels`
        elif ppar.proj.lower() == 'antarctic':
            ax = fig.add_subplot(
                1, 1, 1, projection=ctcrs.SouthPolarStereo())
            ax.set_extent([-180, 180, -90, -60], ctcrs.PlateCarree())
            # gl = ax.gridlines(crs=ctcrs.PlateCarree(), draw_labels=False,
            #                   linewidth=0.8, color='k', alpha=0.8, linestyle=':',
            #                   x_inline=False, y_inline=True, rotate_labels=False)
            # gl.xlocator = mticker.FixedLocator([-180, -120, -60, 0, 60, 120])
            # gl.ylocator = mticker.FixedLocator([-60, -70, -80])
            # plt.draw()  # Enable the use of `gl._labels`
        elif 'EqualEarth' in ppar.proj:
            cen_lon = float(ppar.proj.replace('EqualEarth', ''))
            ic(cen_lon)
            ax = fig.add_subplot(
                1, 1, 1, projection=ctcrs.EqualEarth(central_longitude=cen_lon))
    else:
        ax = plt.subplot(1, 1, 1, projection=ppar.proj) #nrow ncol index
        ax.set_global() 
    
    return ax
    
def stitch_images(file1, file2):
    '''Stitch two images together and display them side by side, with 
    the first resized to be the same height as the second. As configured
    for GDDT library, file1 is a composite and file2 is the
    corresponding timeseries.
    Derived from: https://stackoverflow.com/a/34301747
    
    Arguments:
    file1 -- path to first image file
    file2 -- path to second image file
    
    Returns:
    stitched -- Pillow image object with images side by side, resized
    '''
    image1 = Image.open(file1)
    image2 = Image.open(file2)
    (width1, height1) = image1.size
    (width2, height2) = image2.size
    ratio_rs = height2 / height1
    width1_rs = int(width1 * ratio_rs)
    image1_rs = image1.resize((width1_rs, height2))
    stitch_width = width1_rs + width2
    stitched = Image.new('RGB', (stitch_width, height2))
    stitched.paste(im=image1_rs, box=(0, 0))
    stitched.paste(im=image2, box=(width1_rs, 0))

    return stitched
    
    
    
########################################################################
####  CUSTOM COLORS
########################################################################
def balance_n(n=9):
    ''' n-color discrete diverging palette from cmocean balance '''
    d_rgb = cmocean.tools.get_dict(cmocean.cm.balance, N=n)
    cmap_balancen = mcolors.LinearSegmentedColormap(
        'balance_' + str(n), d_rgb, N=n)

    return cmap_balancen

def diff_n(n=9):
    ''' n-color discrete diverging palette from cmocean diff '''
    d_rgb = cmocean.tools.get_dict(cmocean.cm.diff, N=n)
    cmap_diffn = mcolors.LinearSegmentedColormap(
        'diff_' + str(n), d_rgb, N=n)

    return cmap_diffn

def gdd_trend():
    ''' Custom seaborn diverging palette for GDD trends '''
    gdd_trend_map = sn.diverging_palette(
        273, 129, s=100, l=50, sep=1, as_cmap=True)
    
    return gdd_trend_map
    
def rep_color(color='#000000', d=None):
    ''' Repeat input color value for size of data '''
    tup_color = (color, )
    seq_color = tup_color * len(d)
    
    return seq_color

def rgb_to_hex(rgb):
    return '#%02x%02x%02x' % rgb
    
def tarn_n(n=9):
    ''' n-color discrete diverging palette from cmocean tarn '''
    d_rgb = cmocean.tools.get_dict(cmocean.cm.tarn, N=n)
    cmap_tarnn = mcolors.LinearSegmentedColormap(
        'tarn_' + str(n), d_rgb, N=n)

    return cmap_tarnn