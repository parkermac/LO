"""
Code to make the sta_df.p file for use by extract_casts.py.

This makes a DataFrame with integer Index and columns:
['Station', 'Longitude', 'Latitude', 'Datetime', 'Cruise']

"""

from lo_tools import Lfun
Ldir = Lfun.Lstart()

import pandas as pd
from datetime import datetime, timedelta

in_dir = Ldir['data'] / 'obs' / 'newport_line'
out_dir = Ldir['LOo'] / 'obs' / 'newport_line'
Lfun.make_dir(out_dir)

fn_list = ['sta_df_NH2017.csv', 'sta_df_NH2018.csv', 'sta_df_NH2019.csv']
for fn in fn_list:
    in_fn = in_dir / fn
    a = pd.read_csv(in_fn)
    
    dt_list = []
    for t in a['Datetime']:
        dt_list.append(datetime.strptime(t, '%d-%b-%Y'))
    a['Datetime'] = dt_list
    
    if fn == fn_list[0]:
        A = a.copy()
    else:
        A = pd.concat((A,a))
        
# there are 18 stations that have longitude listed as a positive number, so
# we fix these here, assuming they are "longitude west"
A.loc[A['Longitude']>0, 'Longitude'] = -A.loc[A['Longitude']>0, 'Longitude']
        
A.to_pickle(out_dir / 'sta_df.p')