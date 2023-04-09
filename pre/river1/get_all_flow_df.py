"""
A convenience program for loading the flow DataFrame.
"""

from lo_tools import Lfun
Ldir = Lfun.Lstart()

import pandas as pd

gtag = 'cas6_v3'
year0 = 1980
year1 = 2021

# location of historical data to plot
riv_dir = Ldir['LOo'] / 'pre' / 'river' / gtag / 'Data_historical'
all_df = pd.read_pickle(riv_dir / ('ALL_flow_' + str(year0) + '_' + str(year1) + '.p'))