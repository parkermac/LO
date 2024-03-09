"""
Code to process the NHL (newport hydrography line) ctd_year data to pickle files.
The ctd data range from year 1997-2021.

"""

import pandas as pd
import numpy as np
import gsw
import sys
from time import time as Time
from lo_tools import Lfun, obs_functions
Ldir = Lfun.Lstart()

# CTD
source = 'NHL'
otype = 'ctd'
in_dir0 = Ldir['data'] / 'obs' / source / otype

testing = False

if testing:
    year_list = [2017]
else:
    year_list = range(1997,2022)

# output location
out_dir = Ldir['LOo'] / 'obs' / source / otype
Lfun.make_dir(out_dir)

col_dict = {
    'datetime': 'time',
    'lat': 'lat',
    'lon': 'lon',
    'cruise': 'cruise',
    'station_name': 'name',
    'pressure (dbar)': 'P (dbar)',
    'temperature (degC)': 'IT',
    'salinity': 'SP',
    'DO (ml/L)': 'DO (ml/L)',
}

tt0 = Time()

for year in year_list:
    ys = str(year)
    print('\n'+ys)

    # name output files
    out_fn = out_dir / (ys + '.p')
    info_out_fn = out_dir / ('info_' + ys + '.p')

    if year in year_list:
        in_fn = in_dir0/f'{year}.csv'
        df = pd.read_csv(in_fn, low_memory=False)
        # Rename the columns using col_dict
        df.rename(columns=col_dict, inplace=True)
        df['time'] = pd.to_datetime(df['time'], format='mixed')
        # missing data is -9999
        df[df==-9999] = np.nan
        df = df.dropna(axis=0, how='all')
        df = df.reset_index(drop=True)
        # Next make derived quantities and do unit conversions

    # (1) Create CT, SA, and z
    # - pull out variables
    SP = df.SP.to_numpy()
    IT = df.IT.to_numpy()
    p = df['P (dbar)'].to_numpy()
    lon = df.lon.to_numpy()
    lat = df.lat.to_numpy()
    # - do the conversions
    SA = gsw.SA_from_SP(SP, p, lon, lat)
    CT = gsw.CT_from_t(SA, IT, p)
    z = gsw.z_from_p(p, lat)
    # - add the results to the DataFrame
    df['SA'] = SA
    df['CT'] = CT
    df['z'] = z
    rho = gsw.rho(SA,CT,p)
    # (2) units
    df['DO (uM)'] = df['DO (ml/L)'] * 44.6596

    df['cid'] = np.nan
    cid = 0
    for time in df.time.unique():
        for name in df.name.unique():
            df.loc[(df.name==name) & (df.time==time),'cid'] = cid
            cid += 1

    # add a cruise column
    df['cruise'] = None

    # Keep only selected columns.
    cols = ['cid','time', 'lat', 'lon', 'name', 'z', 'CT', 'SA', 'DO (uM)','cruise']
    this_cols = [item for item in cols if item in df.columns]
    df = df[this_cols]

    if len(df) > 0:
        # Save the data
        df.to_pickle(out_fn)
        info_df = obs_functions.make_info_df(df)
        info_df.to_pickle(info_out_fn)

print('Total time = %d sec' % (int(Time()-tt0)))

