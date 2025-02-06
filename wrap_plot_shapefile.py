''' wrap_plot_shapefile 
Plot shapefile
'''
import sys

from icecream import ic
import matplotlib.pyplot as plt

import classes_gddt as cg
import fun_process as fproc
import fun_plots as fpl

ppar = cg.PlotParams(
    o_name="4_test-mod_arctic_ecoregions.png",
    o_path="/Users/dhueholt/Documents/gddt_fig/20241010_shapes/",
    title='Arctic ecoregions', title_size=12)

geo_arctic, geo_other = fproc.get_arctic_biomes()
#  EPSG:3995 - North Polar Stereographic
geo_arctic_proj = geo_arctic.to_crs(epsg=3995)
no_plot_proj = geo_other.to_crs(epsg=3995)

plot_this = geo_arctic_proj
plt.rcParams.update({'font.family': 'Catamaran'})
#  'font.weight': normal, bold, heavy, light, ultrabold, ultralight
plt.rcParams.update({'font.weight': 'light'})
plt.rcParams.update({'font.size': 10})
fpl.plot_globe_shapefile(plot_this, no_plot_proj, ppar)
