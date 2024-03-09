"""
Code to process the LineP CTD data to pickle files.

-it takes 89 seconds for year 1993-2023,
-this code uses the pre-processed ctd files

"""

import pandas as pd
import numpy as np
import gsw
from time import time as Time

from lo_tools import Lfun, obs_functions
Ldir = Lfun.Lstart()

# This is a dict of all the columns after the initial reading.
# We add values to a key for any variable we want to save.
col_dict = {

    'Mission':'',
    'Project':'cruise',
    'Platform':'',
    'Scientist':'',
    'Date [UTC]':'',
    'Time [UTC]':'',
    'Event_number':'cid',
    'Station':'name',
    'Datetime_UTC':'time',
    'Longitude':'lon',
    'Latitude':'lat',
    'Pressure [dbar]':'P (dbar)',
    'Temperature':'IT',
    'Salinity':'SP',
    'Sigma-t [kg/m^3]':'',
    'Oxygen [ml/l]':'',
    'Oxygen [umol/kg]':'DO (umol/kg)',
    'Fluorescence [mg/m^3]':''

}

# CTD
source = 'LineP'
otype = 'ctd'
in_dir0 = Ldir['data'] / 'obs' / source / otype

testing = False

if testing:
    year_list = [2007]
else:
    year_list = list(range(1993,2023+1))
# year 2000 has no data falls into the model domain
# output location
out_dir = Ldir['LOo'] / 'obs' / source / otype
Lfun.make_dir(out_dir)

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
        df['Datetime_UTC']= pd.to_datetime(df['Date [UTC]'] + ' ' + df['Time [UTC]'], format='mixed')
        # Rename the columns using col_dict
        df.rename(columns=col_dict, inplace=True)
        # missing data is -99
        df[df==-99] = np.nan
        df = df.dropna(axis=0, how='all')
        df = df.reset_index(drop=True)

    """
    Line P goes way beyond the model domain, now limit the geographical region to the model domain, most of the time only P4 is with the domain
    # # but P12 is at the edge of the model domain, -130.67. So, include this station as we can use it to test the boundary condition
    """
    df = df[(df.lon>-131) & (df.lon<-122) & (df.lat>42) & (df.lat<52)]

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
    df['DO (uM)'] = df['DO (umol/kg)'] *  rho

    df['cid'] = np.nan
    cid = 0
    for time in df.time.unique():
        for name in df.name.unique():
            df.loc[(df.name==name) & (df.time==time),'cid'] = cid
        cid += 1

    if year > 2006:
        # add a cruise column
        df['cruise'] = 'Line P'

    # Keep only selected columns.
    cols = ['cid','cruise','time', 'lat', 'lon', 'name', 'z', 'CT', 'SA', 'DO (uM)']
    this_cols = [item for item in cols if item in df.columns]
    df = df[this_cols]

    if len(df) > 0:
        # Save the data
        df.to_pickle(out_fn)
        info_df = obs_functions.make_info_df(df)
        info_df.to_pickle(info_out_fn)

print('Total time = %d sec' % (int(Time()-tt0)))