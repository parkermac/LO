"""
Program to get one year of historical ec data for one river and
add it to the full output from make_historical.py.

The reason for doing this is that sometimes we get gaps in the ec
record because the server times out.  Since it takes many tens of minutes
to run make_historical.py for the full run we instead use this to
fill in gaps as needed.  Hacky but it will do the job.

Before I started running this I made a copy of the file, for safekeeping,
calling it ALL_flow_1980_2020_ORIG.py.

2021.04.16 Using this it was easy for me to fill in all the gaps in the
Fraser River record.
"""
from datetime import datetime, timedelta
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse

from lo_tools import Lfun, zfun
from lo_tools import river_functions as rivf
from importlib import reload
reload(rivf)

Ldir = Lfun.Lstart()

parser = argparse.ArgumentParser()
# required arguments
parser.add_argument('-rn', type=str)   # e.g. fraser
parser.add_argument('-year', type=int) # e.g. 2018
# optional arguments
parser.add_argument('-gtag', type=str, default='cas6_v3') # e.g. cas6_v3
parser.add_argument('-year0', type=int, default=1980) # e.g. 1980
parser.add_argument('-year1', type=int, default=2020) # e.g. 2020

# get the args
args = parser.parse_args()

rn = args.rn
year = args.year
gtag = args.gtag
year0 = args.year0
year1 = args.year1

# Load a dataframe with info for rivers to get
ri_fn = Ldir['LOo'] / 'pre' / 'river' / gtag / 'river_info.csv'
ri_df = pd.read_csv(ri_fn, index_col='rname')

# location for output
in_dir = Ldir['LOo'] / 'pre' / 'river' / gtag / 'Data_historical'
out_name = 'ALL_flow_' + str(year0) + '_' + str(year1) + '.p'
all_df = pd.read_pickle(in_dir / out_name)

rs = ri_df.loc[rn].copy() # a series with info for this river
if pd.notnull(rs.ec):
    rs, qt = rivf.get_ec_data_historical(rs, year)
    if rs['got_data']:
        qt = qt.sort_index() # the ordering of the historical ec data is off
        all_df.loc[qt.index,rn] = qt
        all_df.to_pickle(in_dir / out_name)
        
        plt.close('all')
        all_df[rn].plot()
        plt.show()
        
    else:
        print('Problem getting that year...')
        print(rs)
        
    

