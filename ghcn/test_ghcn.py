'''test_ghcn
Testing GHCN data.

Written by Daniel Hueholt
Graduate Research Assistant at Colorado State University
########################################################################
###############################################################################
'''
import os
import sys

from icecream import ic
import numpy as np
import pandas as pd

ghcn_path = '/Volumes/Polycrystal/Data/ghcn/refined-data/SthrnAfrica_lt3/'
ghcn_stn = 'SF001807220_refined.csv'

df_ghcn = pd.read_csv(ghcn_path + ghcn_stn)
ic(df_ghcn)
sys.exit('STOP')