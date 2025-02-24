'''puppets
These functions are messier and would be more intuitive as scripts--but 
handling them modularly allows them to be run with multiple 
configurations simultaneously using multiprocessing to speed up the 
generation of large numbers of output figures.
'''
from glob import glob
import sys
import time

import cftime
from icecream import ic
import matplotlib.pyplot as plt

import fun_process as fproc
import fun_plots as fpl

def paired_composites_ts(
    dp_guide=None, setp_guide=None, dp_composite=None, setp_composite=None, 
    dp_guide_alltimes=None, ip=None, gp=None, ppar=None, ppar_region=None, 
    ppar_super=None, ppar_guide_ts_all=None, ppar_guide_ts=None,
    pair_timeseries=True, super_interval=30, avg_all_composites=True):
    ''' Plot composites of one variable guided by properties of another
    variable for a list of input intervals, with paired timeseries of 
    the variable used to guide the composite.
    Keyword arguments:
    dp_guide: DataParams for guiding variable
    setp_guide: SetParams for guiding variable
    dp_composite: DataParams for composite variable
    setp_composite: SetParams for composite variable
    dp_guide_alltimes: DataParams for guide variable in region of 
        interest used for paired timeseries
    ip: IntervalParams to make figures
    gp: GuideParams instance for composite
    ppar: PlotParams instance for base layer of composite
    ppar_region: PlotParams instance to annotate region of interest on 
        top of base layer
    ppar_super: PlotParams instance for time-average of a list of 
        composites ("supercomposite")
    ppar_guide_ts_all: PlotParams instance for base layer of timeseries,
        the spaghetti of all members and all times
    ppar_guide_ts: PlotParams instance for second layer of timeseries
        emphasizing the members and times guiding the composite
    pair_timeseries: True/False plot paired timeseries (default: True)
    super_interval: Interval to make supercomposites (default: 30)
    avg_all_composites: True/False to make supercomposite of all 
        intervals at end of run (default: True)

    Returns:
    Returns no variables. Saves combination based on inputs of 
        1) composite for each interval (always)
        2) paired timeseries for each interval (if pair_timeseries=True)
        3) time-averaged supercomposite for each super_interval (as long 
               as super_interval > number of intervals)
        4) time-averaged supercomposite for all intervals (if 
            avg_all_composites=True)
    '''
    #  CHECK BEHAVIOR:
    #      ppar_region.o_bool=False: right now, it seems like this just
    #          wouldn't produce any output, based on the if with no 
    #          else? Check on this.
    l_composites_for_all = list()
    l_composites_for_super = list()
    l_intervals_for_super = list()
    total_tic = time.time()
    for loop_count, interval in enumerate(ip.intervals):
        #  For each interval, want to calculate the anomaly relative to
        #  that interval and guide by GDDs for the same period
        setp_composite.yrs = interval
        setp_composite.yrs_rel_to = interval
        setp_guide.yrs = interval
        tic = time.time()
        #### Set up and guide datasets
        guide_dict = fproc.common_opener(dp=dp_guide, setp=setp_guide)
        da_guide = guide_dict["roi"]
        guide_all_dict = fproc.common_opener(
            dp=dp_guide_alltimes, setp=setp_guide)
        da_guide_alltimes = guide_all_dict["roi"]
        composite_dict = fproc.common_opener(
            dp=dp_composite, setp=setp_composite)
        da_composite_raw = composite_dict["manage_rlz"]
        da_prep = fproc.prep_guide(da_guide, gp)
        #  Indices in guide_indices_dict are used to guide the composite
        guide_indices_dict = fproc.guide(da_prep, gp)
        #### Anomaly calculation and composite
        #  Baseline period for composite
        time_slice_base = slice(
            cftime.DatetimeNoLeap(
                setp_composite.yrs_rel_to[0], 1, 1, 0, 0, 0, 0),
            cftime.DatetimeNoLeap(
                setp_composite.yrs_rel_to[1], 12, 31, 14, 24, 0, 0)
        )
        #  da_base includes all members for best baseline estimate
        da_base = composite_dict["raw_da"].sel(time=time_slice_base)
        da_anom = fproc.calc_anomaly(
            da_composite_raw, da_base, setp_composite)
        da_anom_time_mn = da_anom.mean(dim='time')
        loop_guide_indices = guide_indices_dict[gp.composite_key]
        da_composite = da_anom_time_mn.sel(realization=loop_guide_indices)
        rlz_in_composite_msg = 'Number of realizations in composite: ' \
            + str(len(loop_guide_indices))
        ic(rlz_in_composite_msg)
        #### Plotting and related settings
        plot_this = da_composite.compute()
        name_dict_composite = fproc.namer(da_composite, setp_composite)
        name_dict_guide = fproc.namer(guide_dict["raw_ds"], setp_guide)
        yr_str_base = fproc.str_yrs(setp_composite.yrs_rel_to)
        ppar.title = name_dict_composite['data_id'] + ' ' \
            + name_dict_composite['var_w'] + ' ' \
            + name_dict_composite['yr_str'] + ' based on members with\n' \
            + gp.composite_key + ' ' + gp.guide_by + ' ' \
            + name_dict_guide['var_w'] + ' in ' + name_dict_guide['reg_str']
        #  Wait for later plt.savefig call to add file extension
        loop_filename = name_dict_composite['data_id'] + '_' \
            + name_dict_composite['var_nw'] + '_' \
            + name_dict_composite['yr_str'] + '_relto' + yr_str_base + '_bon' \
            + name_dict_guide['var_nw'] + '_' + gp.guide_by \
            + gp.composite_key + '_' + name_dict_guide['yr_str'] + '_' \
            + name_dict_guide['reg_abv']
        plt.rcParams.update({'font.family': 'Catamaran'})
        #  'font.weight': normal, bold, heavy, light, ultrabold, ultralight
        plt.rcParams.update({'font.weight': 'light'})
        plt.rcParams.update({'font.size': 12})       
        if ppar.plot_each_member:
            for r in plot_this.realization.data:
                loop_rlz = plot_this.sel(realization=r)
                rlz_prefix = 'rlz' + str(r) + '_'
                ppar.o_name = rlz_prefix + loop_filename
                ppar.title = ppar.title + ' rlz ' + str(r + 1)
                plt.figure()
                ic(ppar.o_name, ppar.title)
                fpl.plot_globe(loop_rlz, ppar)
                ppar.o_name = ppar.o_name.replace(rlz_prefix, '')
                ppar.title = ppar.title.replace(' rlz ' + str(r + 1), '')
                ppar.o_name = ppar.o_prefix + loop_filename
                plt.close()
        else:
            loop_rlz_mn = plot_this.mean(dim='realization')
            if setp_composite.z_flag:
                ppar_guide_ts.o_name = loop_filename + '_zscore'
            l_composites_for_super.append(loop_rlz_mn)
            l_intervals_for_super.append(interval)
            if ppar_region.o_bool:
                fig, ax_composite = fpl.plot_globe(loop_rlz_mn, ppar)
                ppar_region.title = ppar.title
                ppar_region.o_name = ppar.o_prefix + loop_filename
                #  If based on region, mask region
                if len(setp_guide.reg_oi["reg_lats"]) > 1:
                    reg_ones = fpl.mask_region(setp_guide.reg_oi)
                    fpl.plot_globe(reg_ones, ppar_region, ax=ax_composite)
                #  Otherwise it's a point, just plot that point
                else:
                    fpl.plot_globe_ng(
                        ppar_region.color, setp_guide.reg_oi["reg_lats"], 
                        setp_guide.reg_oi["reg_lons"], ppar_region, 
                        ax=ax_composite)
                plt.close()    
        if pair_timeseries:
            plt.figure()
            rlz_str = str(
                loop_guide_indices).replace(
                    '[', '').replace(']', '').replace(' ', '-').replace('--', '-')
            ppar_guide_ts.title = name_dict_guide["data_id"] + ' ' \
                + name_dict_guide["reg_str"] + ' ' + name_dict_guide["var_w"] \
                + ' ' + name_dict_guide["yr_str"] + ' ' + gp.composite_key
            #  For composites based on a small number of realizations,
            #  include the indices in the title. In my experience, this
            #  becomes impractical beyond ~5 indices so it is omitted
            #  past that number, though still included in the filename.
            if len(loop_guide_indices) < 6:
                ppar_guide_ts.title = ppar_guide_ts.title + ' ' + rlz_str
            ts_filename = name_dict_guide["data_id"] + '_' \
                + name_dict_guide["var_nw"] + '_' \
                + name_dict_guide["reg_abv"] + '_' \
                + name_dict_guide["yr_str"] + '_' + gp.composite_key \
                + '_' + 'rlz' + rlz_str
            ppar_guide_ts.o_name = ppar.o_prefix + ts_filename
            all_guide = da_guide_alltimes.squeeze()
            all_guide_yrs = da_guide_alltimes.year.data
            #  Base layer for timeseries, all members and times
            fpl.plot_timeseries_spaghetti(
                all_guide, all_guide_yrs, ppar_guide_ts_all)
            loop_guide = da_guide.sel(
                realization=loop_guide_indices).squeeze()
            loop_guide_yrs = da_guide.year.data
            #  Top layer for timeseries, members and times in composite
            fpl.plot_timeseries_spaghetti(
                loop_guide, loop_guide_yrs, ppar_guide_ts)
        toc = time.time() - tic
        time_per_pair_msg = 'Seconds per figure pair: ' + str(round(toc, 2))
        ic(time_per_pair_msg)
        #  Use loop_count to make sure a time-average composite is 
        #  produced each super_interval number of iterations
        if ((loop_count + 1) % super_interval == 0) & (loop_count != 0):
            super_info_strs = str(l_intervals_for_super[0][0]) + '-' \
                + str(l_intervals_for_super[-1][-1]) + ip.type
            cut_strs = (
                ('_' + name_dict_composite['yr_str'] + '_relto' + yr_str_base), 
                (name_dict_guide['yr_str'] + '_'))
            da_composite = fproc.roll_window(l_composites_for_super)
            mean_composite = da_composite.mean(dim='window')
            super_filename = loop_filename.replace(
                cut_strs[0], '').replace(cut_strs[1], '') + '_' \
                + super_info_strs
            ppar_super.o_name = ppar_super.o_prefix + super_filename
            ppar_super.title = name_dict_composite['data_id'] + ' ' \
                + name_dict_composite['var_w'] + ' based on members with\n' \
                + ' ' + gp.composite_key + ' ' + gp.guide_by + ' ' \
                + name_dict_guide['var_w'] + ' ' + ip.type + ' ' \
                + ip.span_str + 'yr ' + super_info_strs.replace(ip.type, '') \
                + ' for ' + name_dict_guide['reg_str']
            if ppar_region.o_bool:
                super_ppar_reg = ppar_region
                super_ppar_reg.o_name = ppar_super.o_name
                super_ppar_reg.title = ppar_super.title
                fig, ax_super = fpl.plot_globe(mean_composite, ppar)
                if len(setp_guide.reg_oi["reg_lats"]) > 1:
                    reg_ones = fpl.mask_region(setp_guide.reg_oi)
                    fpl.plot_globe(reg_ones, super_ppar_reg, ax=ax_super)
                else:
                    fpl.plot_globe_ng(
                        ppar_region.color, setp_guide.reg_oi["reg_lats"], 
                        setp_guide.reg_oi["reg_lons"], super_ppar_reg, 
                        ax=ax_super)
            else:
                fpl.plot_globe(mean_composite, ppar_super)
            if avg_all_composites:
                #  Need to append elements, not list
                for cfs in l_composites_for_super:
                    l_composites_for_all.append(cfs)
            #  Reset for next supercomposite
            l_composites_for_super = list()
            l_intervals_for_super = list()
    if avg_all_composites:
        da_composite = fproc.roll_window(l_composites_for_all)
        da_all_composites = da_composite.mean(dim='window')
        all_ppar = ppar_super
        all_ppar.o_name = 'all_' + ppar_super.o_prefix + loop_filename
        ppar_super.title = name_dict_composite['data_id'] + ' ' \
            + name_dict_composite['var_w'] + ' based on members with\n' \
            + ' ' + gp.composite_key + ' ' + gp.guide_by + ' ' \
            + name_dict_guide['var_w'] + ' ' + ip.type + ' ' + ip.span_str \
            + 'yr ' + ip.strt_yr_str + '-' + ip.end_yr_str + ' for ' \
            + name_dict_guide['reg_str']
        fpl.plot_globe(da_all_composites, ppar_super)
        if ppar_region.o_bool:
            fig, ax_comp = fpl.plot_globe(da_all_composites, ppar)
            all_ppar_region = ppar_region
            all_ppar_region.o_name = ppar_super.o_name
            all_ppar_region.title = ppar_super.title
            if len(setp_guide.reg_oi["reg_lats"]) > 1:
                reg_ones = fpl.mask_region(setp_guide.reg_oi)
                fpl.plot_globe(reg_ones, all_ppar_region, ax=ax_comp)
            else:
                fpl.plot_globe_ng(
                    ppar_region.color, setp_guide.reg_oi["reg_lats"], 
                    setp_guide.reg_oi["reg_lons"], all_ppar_region, 
                    ax=ax_comp)
    total_toc = time.time() - total_tic
    time_total_msg = 'Seconds for entire run: ' + str(round(total_toc, 2))
    ic(time_total_msg)
    
def mp4s(stitch_bool, png_path, tok1, tok2, o_path, o_name):
    if stitch_bool:
        l_stitched = list()
        glob1 = sorted(glob(png_path + tok1))
        glob2 = sorted(glob(png_path + tok2))
        ic(len(glob1), len(glob2))
        for fc, fv in enumerate(glob2):
            stitched = fpl.stitch_images(glob1[fc], fv)
            stitched_path = o_path + '/stitched/'
            fname1 = fv.split('/')[-1]
            stitched.save(stitched_path + 'stitched_' + fname1)
        fpl.images_mp4(stitched_path, '*stitched*', o_path, o_name, fps=6)
    else:
        fpl.images_mp4(png_path, tok1, o_path, o_name, fps=6)