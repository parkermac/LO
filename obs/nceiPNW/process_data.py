"""
Code to process the neciPNW bottle data to pickle files.
"""

import pandas as pd
import numpy as np
import gsw
import sys
from time import time as Time
from lo_tools import Lfun, obs_functions
Ldir = Lfun.Lstart()

# bottle
source = 'nceiPNW'
otype = 'bottle'
in_dir0 = Ldir['data'] / 'obs' / source / otype

testing = False

if testing:
    year_list = [2017]
else:
    year_list = [1985, 1986, 1989] + list(range(1991, 2008+1))+ [2010, 2013, 2014, 2015, 2016, 2017]

# output location
out_dir = Ldir['LOo'] / 'obs' / source / otype
Lfun.make_dir(out_dir)

for year in year_list:
    ys = str(year)
    print('\n'+ys)

    # name output files
    out_fn = out_dir / (ys + '.p')
    info_out_fn = out_dir / ('info_' + ys + '.p')

    if year in year_list:
        in_fn = in_dir0/f'{year}.csv'
        df = pd.read_csv(in_fn, low_memory=False, parse_dates=['time'])
        # missing data is -999
        df[df==-999] = np.nan
    # This dataset doesn't have unique numbers for each cast, if same stations,
    # force lat and lon to be consistent throughout the given year
    for name in df.name.unique():
        df.loc[df.name==name,'lon'] = df[df.name==name].lon.values[0]
        df.loc[df.name==name,'lat'] = df[df.name==name].lat.values[0]
    # limit Strait of Georgia
    df = df[~((df.lon > -126) & (df.lat > 49))]

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
    for vn in ['DO','NO3', 'PO4', 'SIO4','TA','DIC']:
        if (vn+' (umol/kg)') in df.columns:
            df[vn+' (uM)'] = (rho/1000) * df[vn+' (umol/kg)']

    # (3) retain only selected variables
    cols = ['cid', 'cruise', 'time', 'lat', 'lon', 'name', 'z',
            'CT', 'SA', 'DO (uM)',
            'NO3 (uM)', 'PO4 (uM)', 'SiO4 (uM)',
            'TA (uM)', 'DIC (uM)']
    this_cols = [item for item in cols if item in df.columns]
    df = df[this_cols]

    print(' - processed %d casts' % ( len(df.cid.unique()) ))

    # Renumber cid to be increasing from zero in steps of one.
    df = obs_functions.renumber_cid(df)

    if len(df) > 0:
        # Save the data
        df.to_pickle(out_fn)
        info_df = obs_functions.make_info_df(df)
        info_df.to_pickle(info_out_fn)