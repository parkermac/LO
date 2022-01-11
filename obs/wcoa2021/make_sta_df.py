"""
Code to make the sta_df.p file for use by extract_casts.py.

This makes a DataFrame with integer Index and columns:
['Station', 'Longitude', 'Latitude', 'Datetime', 'Cruise']

"""

from lo_tools import Lfun
Ldir = Lfun.Lstart()

import pandas as pd
from datetime import datetime, timedelta

in_dir = Ldir['data'] / 'obs' / 'wcoa2021'
out_dir = Ldir['LOo'] / 'obs' / 'wcoa2021'
Lfun.make_dir(out_dir)

fn = 'WCOA2021_StationLocations_wTimeandDate.xlsx'
in_fn = in_dir / fn
a = pd.read_excel(in_fn)

for ii in a.index:
    a.loc[ii,'Datetime'] = (a.loc[ii,'Date (UTC)'].to_pydatetime()
        + timedelta(days=a.loc[ii,'Time (UTC)'].hour/24))
        
a = a.rename(columns={'Station Name':'Station', 'Lat':'Latitude', 'Long':'Longitude'})
a['Cruise'] = 'WCOA2021'
A = a[['Station', 'Longitude', 'Latitude', 'Datetime', 'Cruise']]
A.to_pickle(out_dir / 'sta_df.p')