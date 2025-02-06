''' classes_gddt
Custom classes for growing degree days and treeline analysis.

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''

            

class DataParams:
    ''' Collects parameters to import data. At present, configured for
    LENS2 or other model output.
    
     __init__: Initiate an instance of the class with attributes
    '''
    def __init__(
        self, path='', tok='', var='', flag_raw_ds=None, flag_raw_da=None,
        flag_time_slice=None, flag_manage_rlz=None, flag_land_mask=None,
        flag_roi=None):
        ''' Constructor for DataParams class containing parameters
        used for locating and opening data. 
        
        Keyword arguments:
        path: path to directory with data (default: '')
        tok: token to match data files (default: '')
        var: data variable within file (default: '')
        flag_ arguments apply to gridded data using common_opener and 
        have valid inputs True/False/None. All default to None.
        flag_raw_ds: flag to output raw dataset
        flag_raw_da: flag to output raw DataArray
        flag_time_slice: flag to output DataArray isolated by time
        flag_manage_rlz: flag to output DataArray after rlz operations
        flag_land_mask: flag to output masked DataArray
        flag_roi: flag to output specific region of interest
        '''
        self.path = path
        self.tok = tok
        self.var = var
        self.flag_raw_ds = flag_raw_ds
        self.flag_raw_da = flag_raw_da
        self.flag_time_slice = flag_time_slice
        self.flag_manage_rlz = flag_manage_rlz
        self.flag_land_mask = flag_land_mask
        self.flag_roi = flag_roi



class GuideParams:
    ''' Collects parameters to guide composites.
    
    __init__: Initiate an instance of the class with attributes
    '''
    
    def __init__(
        self, cmpst_key='max', guide_by='mean', qoi=0.25):
        ''' Constructor for GuideParams class 
        
        Keyword arguments:
        cmpst_key: composite based on keys in dict from guide_by function
        guide_by: property to classify members
        qoi: quantile of interest
        '''
        self.cmpst_key = cmpst_key
        self.guide_by = guide_by
        self.qoi = qoi
        
        
        
class IntervalParams:
    ''' Collects parameters relevant to intervals, and produces those
    intervals.
    
    __init__: Initiate an instance of the class with attributes
    '''
    def __init__(self, span=10, strt_yr=1850, end_yr=2100, type='noverlap'):
        ''' Constructor for IntervalParams class.
                 
        span: Span of interval
        strt_yr: Starting year of interest
        end_yr: Ending year of interest
        type: Type of plot. '''
        self.span = span
        self.span_str = str(span)
        self.strt_yr = strt_yr
        self.strt_yr_str = str(strt_yr)
        self.end_yr = end_yr
        self.end_yr_str = str(end_yr)
        self.type = type
        self.intervals = self.create_intvls(
            strt_yr=strt_yr, end_yr=end_yr, spn=span, type=type)
  
    def create_intvls(
        self, strt_yr=1873, end_yr=2022, spn=30, type='noverlap'):
        ''' Create non-overlapping intervals of length span between 
        given start/end years.
        
        Keyword arguments:
            strt_yr -- first year of first period (default 1873)
            end_yr -- last year of last period (default 2022)
            spn -- span of each period, inclusive (default 30)
            type -- non-overlapping or rolling (default 'noverlap')
        
        Returns:
            l_spns -- list of [strt_yr, end_yr] non-overlapping or rolling
        
        Example:
            >>> print(create_intvls(strt_yr=1903, end_yr=2022, spn=30, type='noverlap'))
            [[1903, 1932], [1933, 1962], [1963, 1992], [1993, 2022]]
        '''
        if type == 'noverlap':
            all_strt_yrs = range(strt_yr, end_yr, spn)
            all_end_yrs = range(strt_yr + spn, end_yr + spn - 1, spn)
        elif type == 'rolling':
            all_strt_yrs = range(strt_yr, end_yr - spn + 1, 1)
            all_end_yrs = range(strt_yr + spn - 1, end_yr, 1)
        elif type == 'fixed_start':
            all_end_yrs = range(strt_yr, end_yr + 1, 1)
            all_strt_yrs = [strt_yr] * len(all_end_yrs)
        else:
            raise ValueError("Invalid type of interval.")
        l_spns = [[syv, all_end_yrs[syc]] for syc, syv in enumerate(all_strt_yrs)]
        
        return l_spns
    



class PlotParams:
    ''' Collects parameters for plotting.
    
    __init__: Initiate an instance of the class with attributes
    '''
    def __init__(
        self, alpha=1, anim_d=dict(frame_tok='', mp4_path=None, mp4_name=None), 
        anim_flag=False, ax_facecolor='#cccccc', bw=None, cb_bool=True, 
        cb_extent='neither', cb_label='', cb_ticks=None, cb_vals='auto', 
        cmap=None, coastline_bool=True, color='', color_r='#cccccc', 
        color_comp='#ff4444', comp_hist_bool=False, comp_label='', 
        dpi=400, edgecolors=None, figsize=(10, 4), forced_crossover_bool=False, 
        frame_flag=False, kde_bool=False, label=None, leg_bool=False, linestyle='--', lw=2,
        marker_size=8, marker='o', member_crossover_bool=False, mn_bool=False, 
        o_bool=True, o_name=None, o_path='', o_prefix='', plot_all=False, 
        proj=None, set_bad=True, stat='count', storyline=False, title=None, title_size=14, 
        ts_type='spaghetti', x_label='', x_lim='auto', xticks='auto', 
        y_label='', y_lim='auto', yticks='auto',
        ):
        ''' Constructor for PlotParams class, containing all possible 
        parameters for all known plots. At present, this includes GHCN 
        animations (wrap_animate_ghcn_map), GHCN timeseries (plot_ghcn_ts),
        various globes (plot_globe), spaghetti (plot_timeseries_spaghetti),
        and histograms (plot_hist).
        
        Keyword arguments:
        alpha: alpha for objects (default: 1)
        anim_d: dict for mp4 animation settings with keys
            frame_tok: tokens for frame images (default: '')
            mp4_path: out mp4 path (default: None, often auto later)
            mp4_name: out mp4 name (default: None, often auto later)
        anim_flag: mp4 animation True/False (default: False)
        ax_facecolor: axis facecolor for figure (default: '#cccccc')
        bw: binwidth for histogram (default: None)
        cb_bool: true/false plot colorbar (default: True)
        cb_label: label for colorbar (default: '')
        cb_vals: colorbar values or 'auto' (default: 'auto')
        cb_ticks: colorbar ticks or None for auto (default: None)
        cb_extent: set pointed colorbar ends (default: 'neither')
        cmap: colormap (default: None)
        coastline_bool: true/false plot coastlines (default: True)
        color: a single color to be used for plotting (default: '')
        color_comp: color for comparison data on hist (default: '#ff4444')
        color_rlz: a single color for plotting realizations (default: '#cccccc')
        comp_hist_bool: plot comparison data on histogram (default: False)
        comp_label: label for comparison data (default: '')
        dpi: dpi of output image or 'pdf' for pdf (default: 400)
        edgecolors: edgecolor for scatterplot (default: None)
        figsize: figure size (default: (10,4))
        forced_crossover_bool: true/false forced crossover (default: False)
        frame_flag: true/false plot individual frames (default: True)
        kde_bool: true/false kernel density estimate on hist (default: False)
        label: label for legend (default: '')
        leg_bool: true/false to plot legend (default: False)
        linestyle: line style (default: '--')
        lw: linewidth value (default: 2)
        marker_size: marker size for scatterplot (default: 8)
        marker: marker style for scatterplot (default: 'o')
        member_crossover_bool: true/false member crossover (default: False)
        mn_bool: true/false plot ensemble mean on timeseries (default: False)
        o_bool: true/false save plot (default: True)
        o_name: output name (default: None, often auto later)
        o_path: output path (default: '')
        o_prefix: add custom prefix to output file (default: '')
        plot_all: true/false plot members separately (default: False)
        proj: map projection or set_proj keyword (default: None)
        set_bad: true/false special color for NaN
        stat: statistic for histogram (default: count)
        storyline: xth percentile or member to plot or False (default: False)
        title: title for figure (default: None, often auto later)
        title_size: size of title font (default: 14)
        ts_type: 'spread' or 'spaghetti' for timeseries
        x_label: x-axis label string (default: '')
        x_lim: x-axis limits [min, max] (default: 'auto')
        xticks: array of x-axis ticks to use (default: 'auto)
        y_label: y-axis label string (default: '')
        y_lim: y-axis limits [min, max] (default: 'auto')
        yticks: array of y-axis ticks to use (default: 'auto')
        '''
        self.alpha = alpha
        self.anim_d = anim_d
        self.anim_flag = anim_flag
        self.ax_facecolor = ax_facecolor
        self.bw = bw
        self.cb_bool = cb_bool
        self.cb_extent = cb_extent
        self.cb_label = cb_label
        self.cb_ticks = cb_ticks
        self.cb_vals = cb_vals
        self.cmap = cmap
        self.coastline_bool = coastline_bool
        self.color = color
        self.color_comp = color_comp
        self.color_r = color_r
        self.comp_hist_bool = comp_hist_bool
        self.comp_label = comp_label
        self.dpi = dpi
        self.edgecolors = edgecolors
        self.figsize = figsize
        self.forced_crossover_bool = forced_crossover_bool
        self.frame_flag = frame_flag
        self.kde_bool = kde_bool
        self.label = label
        self.leg_bool = leg_bool
        self.linestyle = linestyle
        self.lw = lw
        self.marker_size = marker_size
        self.marker = marker
        self.member_crossover_bool = member_crossover_bool
        self.mn_bool = mn_bool
        self.o_bool = o_bool
        self.o_name = o_name
        self.o_path = o_path
        self.o_prefix = o_prefix
        self.plot_all = plot_all
        self.proj = proj
        self.set_bad = set_bad
        self.stat = stat
        self.storyline = storyline
        self.title = title
        self.title_size = title_size
        self.ts_type = ts_type
        self.x_label = x_label
        self.x_lim = x_lim
        self.xticks = xticks
        self.y_label = y_label
        self.y_lim = y_lim
        self.yticks = yticks
        

    
class Sentinel:
    ''' Watches for conditions that require breaking/continuing
    a loop.
    
    __init__: Initiate an instance of the class with attributes
    span: Check first and last index against span
    '''
    def __init__(self, data=None, match=None):
        ''' Constructor for Sentinel class containing parameters
        used to know when to bypass a loop iteration
        
        Keyword arguments:
        data: Data object of interest
        match: Condition to monitor
        '''
        self.data = data
        self.match = match
    
    def not_in_span(self):
        ''' Check first and last index against span '''
        try:
            oside_first = self.data[0] > self.match[0]
            oside_last = self.data[-1] < self.match[1]
            if oside_first | oside_last:
                return True
        except IndexError:
            return True
        else:
            return False



class SetParams:
    ''' Collects parameters used as settings in processing and analysis
    functions.
    
    __init__: Initiate an instance of the class with attributes
    '''
    
    def __init__(
        self, area_stat=None, base_yrs=1850, beat=4999,
        dims=['time', 'realization'], effect='cliffs',
        mask='/Users/dhueholt/Documents/gddt_data/mask/cesm_atm_mask.nc',
        mask_permafrost='/Users/dhueholt/Documents/gddt_data/LENS2/annual_ALTMAX/', 
        mask_flag=None, reg_oi='global', rho=51, rlz='all', window=None,
        yrs=[1850, 1859], yrs_rel_to=[1850, 1859], z_flag=False
        ):
        ''' Constructor for SetParams class.
        Keyword arguments:
        By default, configured for global data over 1850-1859 relative
        to same period as calculated over time and realizations.
        '''
        self.area_stat = area_stat
        self.base_yrs = base_yrs
        self.beat = beat
        self.mask_permafrost = mask_permafrost
        self.dims = dims
        self.effect = effect
        self.mask = mask
        self.mask_flag = mask_flag
        self.reg_oi = reg_oi
        self.rho = rho
        self.rlz = rlz
        self.window = window
        self.yrs = yrs
        self.yrs_rel_to = yrs_rel_to
        self.z_flag = z_flag



class SpanCovParams:
    ''' Collects parameters for selecting stations.
    
    __init__: Initiate an instance of the class with attributes
    '''
    
    def __init__(self, f='', por=None, cov_thr=None, cov_type='>',):
        ''' Constructor for SpanCovParams class.
        Keyword arguments:
        f: Filename for a spancov csv file
        por: Period of record to filter by
        cov_thr: Coverage threshold to filter by
        cov_type: Coverage type to filter by (valid: >, <, or =)
        '''
        self.por = por
        self.por_str = str(por)
        self.f = f
        self.span_thr = por * 365 + (por / 4).__floor__()
        self.cov_thr = cov_thr
        self.cov_thr_str = str(cov_thr)
        self.cov_type = cov_type



class TrackProg:
    ''' Track progress in a loop. TODO: Folding this functionality into the 
    Sentinel class would be ideal, but not a priority at this moment.
    
    __init__: Initiate an instance of the class with attributes
    message: Message shell about current progress
    '''
    
    def __init__(self, cli=None, iter=None, sf=2):
        ''' Constructor for  __init__ class.
        Keyword arguments:
        cli: Current loop index
        iter: Loop object iterator
        sf: Significant figures to report
        '''
        self.frac = cli / len(iter)
        self.percent = round(self.frac * 100, 2)
        
    def message(self, prefix=''):
        ''' Message about current progress '''
        prog_msg = 'Progress: ' + prefix + str(self.percent) + '%'
        return prog_msg