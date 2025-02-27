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

import fun_calc_var as fcv
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

def plot_exceed_crossover_ts(to_plt, x_d=None, ppar=None):
    ''' Make a timeseries of exceedance with crossover thresholds 
    and storylines annotated.

    Hand-tuned aesthetics:
        x-limits of exceedance threshold (horizontal lines)
        y-limits of vertical crossover lines
        y-limits of quantile range
        Top and right spines are manually disabled
    
    Arguments:
    to_plt -- object to plot with exceedance data, ideally np array
    
    Keyword arguments:
    x_d -- abscissa data
    ppar -- PlotParams instance

    Returns: None at present
    '''
    # years_span = [x_d[0], x_d[-1]]
    years_span = [1850, 2100]
    to_plt_mn = np.mean(to_plt, axis=1)
    #  Plot all members
    if ppar.rlz_dict['bool']:
        plt.plot(
            x_d, to_plt, color=ppar.rlz_dict['color'], 
            label=ppar.rlz_dict['label'], linestyle=ppar.rlz_dict['linestyle'],
            lw=ppar.rlz_dict['lw'])
    #  Ensemble mean plotting
    if ppar.mn_dict['bool']:
        plt.plot(
            x_d, to_plt_mn, color=ppar.mn_dict['color'], 
            label=ppar.mn_dict['label'], linestyle=ppar.mn_dict['linestyle'],
            lw=ppar.mn_dict['lw'])
    #  Forced crossover calculation and plotting
    if ppar.plot_crossover_dict['forced_dict']['bool']:
        ppar.o_name = 'forced_' + ppar.o_name
        forced_dict = ppar.plot_crossover_dict['forced_dict']
        for sync_key in (
            'exceed_color', 'exceed_linestyle', 'exceed_lw', 'years_alpha',
            'years_color', 'years_linestyle', 'years_lw'):
            forced_dict = fproc.sync_lengths(
                forced_dict, sync_key=sync_key, ref_key='threshold')
        for forced_count, forced_threshold in enumerate(
            forced_dict['threshold']):
            #  Exceedance threshold(s) plotted as horizontal lines
            plt.plot(
                years_span, [forced_threshold, forced_threshold], 
                color=forced_dict['exceed_color'][forced_count],
                linestyle=forced_dict['exceed_linestyle'][forced_count], 
                lw=forced_dict['exceed_lw'][forced_count])
            #  Calculate forced crossover, plot as vertical lines, print 
            #  strings of when this occurs for each requested threshold.
            forced_crossover = fcv.calc_crossover(to_plt_mn, forced_threshold)
            try:
                forced_crossover_year = x_d[forced_crossover][0]
                plt.plot(
                    [forced_crossover_year, forced_crossover_year],
                    [-100, 100], 
                    alpha=forced_dict['years_alpha'][forced_count], 
                    color=forced_dict['years_color'][forced_count], 
                    label=str(forced_crossover_year),
                    linestyle=forced_dict['years_linestyle'][forced_count],
                    lw=forced_dict['years_lw'][forced_count])
                msg_crossover = 'Forced crossover ' + str(forced_threshold) \
                    + '%' + ' threshold: ' + str(forced_crossover_year)
            except IndexError:
                msg_crossover = 'No forced crossover at ' \
                    + str(forced_threshold) + '%' + ' threshold'
            ic(msg_crossover)
    #  Member crossover calculation and plotting
    if ppar.plot_crossover_dict['member_dict']['bool']:
        ppar.o_name = 'member_' + ppar.o_name
        member_dict = ppar.plot_crossover_dict['member_dict']
        for sync_key in (
            'exceed_color', 'exceed_linestyle', 'exceed_lw', 'years_color',
            'years_linestyle', 'years_lw'):
            member_dict = fproc.sync_lengths(
                member_dict, sync_key=sync_key, ref_key='threshold')
        for rlz_threshold_count, rlz_threshold in enumerate(
            member_dict['threshold']):
            if rlz_threshold == 100:
                ppar.o_name = 'noan_' + ppar.o_name
            #  Exceedance threshold(s) plotted as horizontal lines
            plt.plot(
                years_span, [rlz_threshold, rlz_threshold], 
                alpha=member_dict['exceed_alpha'][rlz_threshold_count], 
                color=member_dict['exceed_color'][rlz_threshold_count],
                linestyle=member_dict['exceed_linestyle'][rlz_threshold_count], 
                lw=member_dict['exceed_lw'][rlz_threshold_count])
            l_member_crossover = list()
            all_crossover = fcv.calc_crossover(to_plt, rlz_threshold)
            for member_exceedance in all_crossover.T:
                try:
                    member_crossover = x_d[member_exceedance][0]
                    l_member_crossover.append(member_crossover)
                except IndexError:
                    l_member_crossover.append(np.nan)
            np_member_crossover = np.array(l_member_crossover)
            #  Storylines of member crossover
            try:
                storyline_dict = ppar.storyline_dict
                #  If a quantile is input, find matching members first
                if (storyline_dict['select'][0] < 1) \
                    | (np.isnan(storyline_dict['select'][0])):
                    storyline_dict['select'] = fproc.match_rlz_quantiles(
                        np_member_crossover, storyline_dict['select'])
                #  Plot each storyline
                for story_count, story in enumerate(storyline_dict['select']):
                    to_plt_story = to_plt[:, story]
                    story_crossover = np_member_crossover[story]
                    msg_story_crossover = 'Member crossover ' \
                        + str(rlz_threshold) + '%' \
                        + ' threshold for storyline ' + str(story) + ': ' \
                        + str(story_crossover)
                    ic(msg_story_crossover)
                    story_str = 'story' \
                        + str(
                            storyline_dict['select']).replace(
                                '[','').replace(
                                    ']','').replace(' ','').replace(',','-')
                    ppar.o_name = story_str + '_' + ppar.o_name
                    plt.plot(
                        x_d, to_plt_story, 
                        color=storyline_dict['color'][story_count], 
                        label=storyline_dict['label'][story_count], 
                        linestyle=storyline_dict['linestyle'][story_count],
                        lw=storyline_dict['lw'][story_count])
                    plt.plot(
                        [story_crossover, story_crossover], [-100, 100], 
                        alpha=member_dict['years_alpha'][rlz_threshold_count], 
                        color=member_dict['years_color'][rlz_threshold_count], 
                        label=str(story_crossover),
                        linestyle=member_dict[
                            'years_linestyle'][rlz_threshold_count],
                        lw=member_dict['years_lw'][rlz_threshold_count])     
            except TypeError:
                #  Do nothing if storyline_dict['select'] is None, False
                pass
            #  Plot range of member crossovers
            if member_dict['range_dict']['bool']:
                quantile_lower = np.nanquantile(
                    np_member_crossover, member_dict['range_dict']['range'][0])
                quantile_upper = np.nanquantile(
                    np_member_crossover, member_dict['range_dict']['range'][1])
                msg_range = 'Range at ' + str(rlz_threshold) + '% threshold ' \
                    + str(member_dict['range_dict']['range'][0] * 100) \
                    + ' to ' + str(
                        member_dict['range_dict']['range'][1] * 100) \
                    + ' percentile: ' + str(int(quantile_lower)) + ' to ' \
                    + str(int(quantile_upper))
                ic(msg_range)
                plt.fill_between(
                    [quantile_lower, quantile_upper], [-100, -100], [100, 100], 
                    alpha=member_dict['range_dict']['alpha'],
                    color=member_dict['range_dict']['color'], 
                    edgecolor=member_dict['range_dict']['edgecolor'])
    #  Disable spines, apply timeseries parameters
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
            plt.savefig(ppar.o_path + ppar.o_prefix + ppar.o_name + '.pdf')
        else:
            plt.savefig(
                ppar.o_path + ppar.o_prefix + ppar.o_name + '.png', 
                dpi=ppar.dpi)
    
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