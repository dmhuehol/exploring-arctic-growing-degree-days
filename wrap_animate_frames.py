''' wrap_animate_frames
Animate mp4 from directory of still images. Characteristic
and region are generated automatically for filename, rest of
information is manual.

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
'''
import sys

from icecream import ic

import puppets

stitch_bool = True
cmn_path = '/home/dhueholt/gddt_fig/frames/'
tup_reg = (
    'BrooksCoast/PSL/3yr/max20',
    'EurasianTundra/PSL/3yr/max20', 
    'Tasiilaq/PSL/3yr/max20',
    'UpperQuebec/PSL/3yr/max20',
    'BrooksCoast/PSL/3yr/min20',
    'EurasianTundra/PSL/3yr/min20', 
    'Tasiilaq/PSL/3yr/min20',
    'UpperQuebec/PSL/3yr/min20')
tok1 = '*bonGDD*.png'
tok2 = '*_GDD_*.png'
cut_str = (
    '/arctic_proj', '/min20', '/max20', '/PSL', 'ICEFRAC', '/10yr', '/3yr'
    )

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
            + '_10yr_PSLComposite_min20_18502100.mp4'
    else:
        o_name = 'rolling_' + reg_only \
            + '_10yr_PSLComposite_max20_18502100.mp4'
    ic(png_path)
    puppets.mp4s(stitch_bool, png_path, tok1, tok2, o_path, o_name)
ic('Completed!')