'''puppets
This contains functions that in a more reasonable world would be their
own script. In reality, handling them modularly allows them to be
run with multiple configurations simultaneously. 

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University 
'''
from glob import glob
import sys
import time

import cftime
import cmocean as cmo
from icecream import ic
import matplotlib.cm as mcm
import matplotlib.pyplot as plt

import classes_gddt as cg
import fun_process as fproc
import fun_plots as fpl
import gddt_region_library as g_rlib

def composites(
    dp_gde, setp_gde, dp_cmpst, setp_cmpst, dp_gde_roi_alltimes, ip, gp, 
    ppar, ppar_reg, ppar_super, ppar_gdd_ts_gray, ppar_gdd_ts,
    pair_timeseries=True, super_int=30, all_comp=True):
    l_windows = list()
    l_panes = list()
    l_panes_intvl = list()
    total_tic = time.time()
    for intvlc, intvl in enumerate(ip.intervals):
        setp_cmpst.yrs_rel_to = intvl
        setp_cmpst.yrs = intvl
        setp_gde.yrs = intvl
        tic = time.time()
        #### Set up and guide datasets
        gde_d = fproc.common_opener(dp=dp_gde, setp=setp_gde)
        da_gde = gde_d["roi"]
        gde_all_d = fproc.common_opener(dp=dp_gde_roi_alltimes, setp=setp_gde)
        da_roi_alltimes = gde_all_d["roi"]
        cmpst_d = fproc.common_opener(dp=dp_cmpst, setp=setp_cmpst)
        da_cmpst_raw = cmpst_d["manage_rlz"]
        da_prep = fproc.prep_guide(da_gde, gp)
        gde_ind_d = fproc.guide(da_prep, gp)
        #### Anomaly calculation and composite
        time_slice_base = slice(
            cftime.DatetimeNoLeap(setp_cmpst.yrs_rel_to[0], 1, 1, 0, 0, 0, 0),
            cftime.DatetimeNoLeap(setp_cmpst.yrs_rel_to[1], 12, 31, 14, 24, 0, 0)
        )
        #  da_base includes all members
        da_base = cmpst_d["raw_da"].sel(time=time_slice_base)
        da_anom = fproc.calc_anomaly(da_cmpst_raw, da_base, setp_cmpst)
        da_anom_tmn = da_anom.mean(dim='time')
        act_ind = gde_ind_d[gp.cmpst_key]
        da_cmpst = da_anom_tmn.sel(realization=act_ind)
        ic(len(act_ind))
        #### Plotting and related settings
        plot_this = da_cmpst.compute()
        nd_cmpst = fproc.namer(da_cmpst, setp_cmpst)
        nd_gde = fproc.namer(gde_d["raw_ds"], setp_gde)
        yr_str_base = fproc.str_yrs(setp_cmpst.yrs_rel_to)
        ppar.title = nd_cmpst['data_id'] + ' ' + nd_cmpst['var_w'] + ' ' \
            + nd_cmpst['yr_str'] + ' based on ' + nd_gde['var_w'] + ' ' \
            + gp.guide_by + ' ' + gp.cmpst_key + ' ' + nd_gde['yr_str'] \
            + ' ' + 'for ' + nd_gde['reg_str']
        #  Wait for plt.savefig to add file extension to allow further edits
        plot_name = nd_cmpst['data_id'] + '_' + nd_cmpst['var_nw'] + '_' \
            + nd_cmpst['yr_str'] + '_relto' + yr_str_base + '_bon' \
            + nd_gde['var_nw'] + '_' + gp.guide_by + gp.cmpst_key + '_' \
            + nd_gde['yr_str'] + '_' + nd_gde['reg_abv']
        plt.rcParams.update({'font.family': 'Catamaran'})
        #  'font.weight': normal, bold, heavy, light, ultrabold, ultralight
        plt.rcParams.update({'font.weight': 'light'})
        plt.rcParams.update({'font.size': 12})       
        if ppar.plot_all:
            for r in plot_this.realization.data:
                act = plot_this.sel(realization=r)
                rlz_pfx = 'rlz' + str(r) + '_'
                ppar.o_name = rlz_pfx + plot_name
                ppar.title = ppar.title + ' rlz ' + str(r + 1)
                plt.figure()
                ic(ppar.o_name, ppar.title)
                fpl.plot_globe(act, ppar)
                ppar.o_name = ppar.o_name.replace(rlz_pfx, '')
                ppar.title = ppar.title.replace(' rlz ' + str(r + 1), '')
                ppar.o_name = ppar.o_prefix + plot_name
                plt.close()
        else:
            act = plot_this.mean(dim='realization')
            if setp_cmpst.z_flag:
                ppar_gdd_ts.o_name = plot_name + '_zscore'
            l_panes.append(act)
            l_panes_intvl.append(intvl)
            if ppar_reg.o_bool:
                fig, ax_comp = fpl.plot_globe(act, ppar)
                ppar_reg.o_name = ppar.o_prefix + plot_name
                if len(setp_gde.reg_oi["reg_lats"]) > 1:
                    reg_ones = fpl.mask_region(setp_gde.reg_oi)
                    fpl.plot_globe(reg_ones, ppar_reg, ax=ax_comp)
                else:
                    fpl.plot_globe_ng(
                        ppar_reg.color, setp_gde.reg_oi["reg_lats"], 
                        setp_gde.reg_oi["reg_lons"], ppar_reg, ax=ax_comp)
                plt.close()    
        if pair_timeseries:
            plt.figure()
            rlz_str = str(act_ind).replace('[', '').replace(']', '').replace(' ', '-')
            ppar_gdd_ts.title = nd_gde["data_id"] + ' ' + nd_gde["reg_str"] + ' ' \
                + nd_gde["var_w"] + ' ' + nd_gde["yr_str"] + ' ' + gp.cmpst_key
            if len(act_ind) < 6:
                ppar_gdd_ts.title = ppar_gdd_ts.title + ' ' + rlz_str
            ts_name = nd_gde["data_id"] + '_' + nd_gde["var_nw"] + '_' \
                + nd_gde["reg_abv"] + '_' + nd_gde["yr_str"] + '_' + gp.cmpst_key \
                + '_' + 'rlz' + rlz_str
            ppar_gdd_ts.o_name = ppar.o_prefix + ts_name
            all_gdd = da_roi_alltimes.squeeze()
            all_gdd_yrs = da_roi_alltimes.year.data
            fpl.plot_timeseries_spaghetti(all_gdd, all_gdd_yrs, ppar_gdd_ts_gray)
            act_gdd = da_gde.sel(realization=act_ind).squeeze()
            act_gdd_yrs = da_gde.year.data
            fpl.plot_timeseries_spaghetti(act_gdd, act_gdd_yrs, ppar_gdd_ts)
        toc = time.time() - tic; ic(toc)
        if ((intvlc + 1) % super_int == 0) & (intvlc != 0):
            sup_yrs = (str(l_panes_intvl[0][0]), str(l_panes_intvl[-1][-1]))
            sups = sup_yrs[0] + '-' + sup_yrs[1] + ip.type
            cs = (
                ('_' + nd_cmpst['yr_str'] + '_relto' + yr_str_base), 
                (nd_gde['yr_str'] + '_'))
            da_composite = fproc.roll_window(l_panes)
            mean_composite = da_composite.mean(dim='window')
            super_name = plot_name.replace(cs[0], '').replace(cs[1], '') + '_' + sups
            ppar_super.o_name = ppar_super.o_prefix + super_name
            ppar_super.title = nd_cmpst['data_id'] + ' ' + nd_cmpst['var_w'] \
                + ' based on ' + nd_gde['var_w'] + ' ' + gp.guide_by + ' ' \
                + gp.cmpst_key + ' ' + ip.type + ' ' + ip.span_str + 'yr ' \
                + sups.replace(ip.type, '') + ' for ' + nd_gde['reg_str']
            if ppar_reg.o_bool:
                super_ppar_reg = ppar_reg
                super_ppar_reg.o_name = ppar_super.o_name
                super_ppar_reg.title = ppar_super.title
                fig, ax_comp = fpl.plot_globe(mean_composite, ppar)
                if len(setp_gde.reg_oi["reg_lats"]) > 1:
                    reg_ones = fpl.mask_region(setp_gde.reg_oi)
                    fpl.plot_globe(reg_ones, super_ppar_reg, ax=ax_comp)
                else:
                    fpl.plot_globe_ng(
                        ppar_reg.color, setp_gde.reg_oi["reg_lats"], 
                        setp_gde.reg_oi["reg_lons"], super_ppar_reg, ax=ax_comp)
            else:
                fpl.plot_globe(mean_composite, ppar_super)
            if all_comp:
                #  Need to append elements, not list
                for lp in l_panes:
                    l_windows.append(lp)
            #  Reset for next supercomposite
            l_panes = list()
            l_panes_intvl = list()
    if all_comp:
        da_composite = fproc.roll_window(l_windows)
        mean_composite = da_composite.mean(dim='window')
        ppar_super.o_name = 'all_' + ppar_super.o_prefix + plot_name
        ppar_super.title = nd_cmpst['data_id'] + ' ' + nd_cmpst['var_w'] \
            + ' based on ' + nd_gde['var_w'] + ' ' + gp.guide_by + ' ' \
            + gp.cmpst_key + ' ' + ip.type + ' ' + ip.span_str + 'yr ' \
            + ip.strt_yr_str + '-' + ip.end_yr_str + ' for ' + nd_gde['reg_str']
        fpl.plot_globe(mean_composite, ppar_super)
        if ppar_reg.o_bool:
            fig, ax_comp = fpl.plot_globe(mean_composite, ppar)
            all_ppar_reg = ppar_reg
            all_ppar_reg.o_name = ppar_super.o_name
            all_ppar_reg.title = ppar_super.title
            if len(setp_gde.reg_oi["reg_lats"]) > 1:
                reg_ones = fpl.mask_region(setp_gde.reg_oi)
                fpl.plot_globe(reg_ones, all_ppar_reg, ax=ax_comp)
            else:
                fpl.plot_globe_ng(
                    ppar_reg.color, setp_gde.reg_oi["reg_lats"], 
                    setp_gde.reg_oi["reg_lons"], all_ppar_reg, ax=ax_comp)
    total_toc = time.time() - total_tic; ic(total_toc)
    
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