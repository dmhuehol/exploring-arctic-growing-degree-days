''' wrap_animate_frames
Animate mp4 from directory of still images. Characteristic and region 
are generated automatically for filename; the rest of the information is 
manual. This and assumptions about the directory structure can require a
little doctoring to work from run to run.

Assumes constituent images located in the following directory structure:
cmn_path/RegionName/Variable/Nyr/Characteristic/*.png
For animations where multiple images are stitched together to make a 
frame like the paired map-timeseries plots, this directory must include 
a subfolder named "stitched".

If you see the error: "AVF: AVAssetWriter status: Cannot create file",
this means the filename that fun_plots.images_mp4 is attempting to write
is invalid (often because of an errant /). Check your filenames and 
directory structures, then try again!
'''
import sys

from icecream import ic

import puppets

stitch_bool = True
# cmn_path = '/home/dhueholt/gddt_fig/frames/'
# tup_reg = (
#     'BrooksCoast/PSL/3yr/max20',
#     'EurasianTundra/PSL/3yr/max20', 
#     'Tasiilaq/PSL/3yr/max20',
#     'UpperQuebec/PSL/3yr/max20',
#     'BrooksCoast/PSL/3yr/min20',
#     'EurasianTundra/PSL/3yr/min20', 
#     'Tasiilaq/PSL/3yr/min20',
#     'UpperQuebec/PSL/3yr/min20')
cmn_path = '/Users/dhueholt/Documents/gddt_fig/FramesForAnimations/'
#  No trailing slash on tup_reg
tup_reg = (
    'BrooksRangeColonist/SST/10yr/max20',)
tok1 = '*bonGDD*.png'
tok2 = '*_GDD_*.png'
#  Substrings of path which are cut to automatically create mp4 filename
cut_str = (
    '/arctic_proj', '/min20', '/max20', '/PSL', 'ICEFRAC', 'SST', '/10yr', 
    '/3yr')

for reg in tup_reg:
    reg_only = reg
    for cs in cut_str:
        reg_only = reg_only.replace(cs, '')
    if reg_only[-1] == '/':
        reg_only = reg_only[:-1]
    ic(reg, reg_only)
    png_path = cmn_path + reg + '/'
    o_path = png_path
    if 'min20' in reg:
        o_name = 'rolling_' + reg_only \
            + '_10yr_SSTComposite_min20_18502100.mp4'
    else:
        o_name = 'rolling_' + reg_only \
            + '_10yr_SSTComposite_max20_18502100.mp4'
    ic(png_path)
    puppets.mp4s(stitch_bool, png_path, tok1, tok2, o_path, o_name)
ic('Completed!')