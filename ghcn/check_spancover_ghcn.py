''' check_spancover_ghcn 
For input GHCN data files, create dataframe with location metadata,
span of record, total NaN degree days, and percent coverage. This 
information must be calculated explicitly from the data as the
station list file does not contain coverage data and the years are 
unreliable (see e.g., SF001807220).
'''
import collections
from glob import glob
import sys

from icecream import ic
import numpy as np
import pandas as pd

d_path = "/Volumes/Polycrystal/Data/ghcn/refined-data/Arctic50N_gt1/"
#  Adjust token to target specific station(s)
d_tok = "*.csv"
out_path = "/Users/dhueholt/Documents/ghcn_data/spancov_lists/"
out_name = "spancov_Arctic50N_gt1.csv"
save_flag = True

ghcn_files = glob(d_path + d_tok)
ghcn_fnames = [ghcn_str.split('/')[-1] for ghcn_str in ghcn_files]
dtypes_d = collections.defaultdict(lambda: str, {
    #  Default str for everything that is not truly necessary
    "LATITUDE": np.float32,
    "LONGITUDE": np.float32,
    "ELEVATION": np.float32,
    "NDD_SUM": np.float32,
    }
)
lf = collections.defaultdict(list)
for gfc, gfil in enumerate(ghcn_files):
    if gfc % 100 == 0:
        prog_msg = 'Progress: ' + str(
            np.round(gfc / len(ghcn_files) * 100, 2)) + '%'
        ic(prog_msg)
    df_ghcn = pd.read_csv(gfil, dtype=dtypes_d)
    lf['st_id'].append(df_ghcn['STATION'][0])
    lf['lat'].append(df_ghcn['LATITUDE'][0])
    lf['lon'].append(df_ghcn['LONGITUDE'][0])
    lf['elev'].append(df_ghcn['ELEVATION'][0])
    df_ghcn['DATETIME'] = pd.to_datetime(df_ghcn['DATETIME'])
    last_dt = df_ghcn['DATETIME'].iloc[-1]
    first_dt = df_ghcn['DATETIME'].iloc[0]
    act_span_days = (last_dt - first_dt).days
    lf['span_days'].append(act_span_days)
    ndd_sum = df_ghcn['NDD_SUM'].sum()
    pct_bad = ndd_sum / lf['span_days'][gfc] * 100
    lf['coverage'].append(100 - pct_bad)
df_spancov = pd.DataFrame(
    {"st_id": lf['st_id'],
    "lat": lf['lat'],
    "lon": lf['lon'],
    "elev": lf['elev'],
    "span_days": lf['span_days'],
    "coverage": lf['coverage']}
)
if save_flag:
    df_spancov.to_csv(out_path + out_name)

