"""
Program to get one year of historical ec data for one river and
add it to the full output from make_historical.py.

The reason for doing this is that sometimes we get gaps in the ec
record because the server times out.  Since it takes many tens of minutes
to run make_historical.py for the full run we instead use this to
fill in gaps as needed.  Hacky but it will do the job.

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
parser.add_argument('-rn', type=str, default='fraser')   # e.g. fraser
parser.add_argument('-year', type=int, default=2018) # e.g. 2018
parser.add_argument('-ctag', type=str, default='lo_base') # e.g. lo_base

# get the args
args = parser.parse_args()

rn = args.rn
year = args.year
ctag = args.ctag

riv_dir0 = Ldir['LOo'] / 'pre' / 'river1' / ctag
# Load a dataframe with info for rivers to get
ri_df_fn = riv_dir0 / 'river_info.p'
ri_df = pd.read_pickle(ri_df_fn)

# location for output
in_dir = riv_dir0 / 'Data_historical'
out_name = 'ALL_flow.p'
all_df = pd.read_pickle(in_dir / out_name)

# save a copy
all_df.to_pickle(in_dir / 'ALL_flow_backup.p')

rs = ri_df.loc[rn].copy() # a series with info for this river
if pd.notnull(rs.ec):
    rs, qt = rivf.get_ec_data_historical(rs, year)
    if rs.got_data:
        qt = qt.sort_index() # the ordering of the historical ec data is off
        all_df.loc[qt.index,rn] = qt
        
        # save the updated DataFrame
        all_df.to_pickle(in_dir / out_name)
        
        # plot to check on things
        plt.close('all')
        all_df[rn].plot()
        plt.show()
        
    else:
        print('Problem getting that year...')
        print(rs)
        
    

