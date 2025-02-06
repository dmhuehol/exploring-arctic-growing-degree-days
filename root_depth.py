from icecream import ic
import matplotlib.pyplot as plt
import numpy as np

import classes_gddt as cg
import fun_plots as fp

ppar_hist = cg.PlotParams(
    bw=2, leg_bool=True,
    o_path='/Users/dhueholt/Documents/gddt_fig/20240925_masks/', 
    o_name='minrootdepth', title='Min root depth from Iversen et al. (2014)', 
    stat='count', x_label='cm', y_lim=[0, 8.5])

root_depth = np.array(
    [20, 50, 3.3, 10, 3.7, 15, 30, 2.9, 25, 35, 40, 30, 3, 40, 30, 20, 1.8, 
     0.7, 10, 12.5, 1.3, 10, 10, 0.9, 20, 10, 2.9, 5, 9.5, 1.6, 100, 30, 7, 
     1.1, 30, 20, 1.5, 4, 5, 40, 4.6, 1.5])
ic(root_depth)
plot_d = {
    'root_depth': root_depth,
}
plt.rcParams.update({'font.family': 'Catamaran'})
#  'font.weight': normal, bold, heavy, light, ultrabold, ultralight
plt.rcParams.update({'font.weight': 'light'})
plt.rcParams.update({'font.size': 12})
fp.plot_hist(plot_d, ppar=ppar_hist)