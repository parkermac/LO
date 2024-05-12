"""
Code to process the ecology ctd data.

Takes a 2 minutes to run for 2008-2019.

"""

import pandas as pd
import numpy as np
import gsw
import sys
from time import time as Time

from lo_tools import Lfun, obs_functions
Ldir = Lfun.Lstart()

# BOTTLE
source = 'ecology'
otype = 'ctd'
in_dir0 = Ldir['data'] / 'obs' / source

testing = False

if testing:
    year_list = [2017]
else:
    year_list = range(2008,2020)

# output location
out_dir = Ldir['LOo'] / 'obs' / source / otype
Lfun.make_dir(out_dir)

# We need to associate lat and lon with each station. They are not stored in the bottle
# file, but there are station names, with locations here:
sta_fn = in_dir0 / 'sta_df.p'
sta_df = pd.read_pickle(sta_fn)
xx = sta_df['Long_NAD83 (deg / dec_min)'].values
yy = sta_df['Lat_NAD83 (deg / dec_min)'].values
lon = [-(float(x.split()[0]) + float(x.split()[1])/60) for x in xx]
lat = [(float(y.split()[0]) + float(y.split()[1])/60) for y in yy]
sta_df['lon'] = lon
sta_df['lat'] = lat

tt0 = Time()
load_data = True
for year in year_list:
    ys = str(year)
    print('\n'+ys)
    
    # name output files
    out_fn = out_dir / (ys + '.p')
    info_out_fn = out_dir / ('info_' + ys + '.p')
    
    if year == 2017:
        load_data = True
        ctd_fn = in_dir0 / 'ctd_2017_fixed.p'
    elif year == 2018:
        load_data = True
        ctd_fn = in_dir0 / 'ctd_2018_fixed.p'
    elif year == 2019:
        load_data = True
        ctd_fn = in_dir0 / 'ctd_2019_fixed.p'
    else:
        ctd_fn = in_dir0 / 'ctd_1999_2016_fixed.p'
    
    # DEFAULTS
    data_new_names = ['SP', 'IT', 'DO (mg/L)', 'Chl (mg m-3)']
    
    if year in [2017, 2019]:
        data_original_names = ['Salinity', 'Temp', 'DO_raw', 'Chla_adjusted']
    else:
        data_original_names = ['Salinity', 'Temp', 'DO_adjusted', 'Chla_adjusted']
    
    # read in the data (all stations, all casts)
    if load_data:
        df0 = pd.read_pickle(ctd_fn)
        load_data = False
        
    # add z
    df0['z'] = -df0['Depth']
        
    # select and rename variables
    v_dict = dict(zip(data_original_names, data_new_names))
    v_dict['Date'] = 'time'
    v_dict['Station'] = 'name'
    v_dict['z'] = 'z'
    df1 = pd.DataFrame()
    for v in df0.columns:
        if v in v_dict.keys():
            if len(v_dict[v]) > 0:
                df1[v_dict[v]] = df0[v]
                
    # select one year
    t = pd.DatetimeIndex(df1.time)
    df = df1.loc[t.year==year,:].copy()
    
    # units
    df['DO (uM)'] = (1000/32) * df['DO (mg/L)']
    
    # add lon and lat
    df['lon'] = np.nan
    df['lat'] = np.nan
    for sn in sta_df.index:
        df.loc[df.name==sn,'lon'] = sta_df.loc[sn,'lon']
        df.loc[df.name==sn,'lat'] = sta_df.loc[sn,'lat']
    # and drop stations without good lon, lat, or z
    df = df[df.lon.notna()]
    df = df[df.lat.notna()]
    df = df[df.z.notna()]
    
    # calculated quantities
    SP = df.SP.to_numpy()
    IT = df.IT.to_numpy()
    z = df.z.to_numpy()
    lon = df.lon.to_numpy()
    lat = df.lat.to_numpy()
    p = gsw.p_from_z(z, lat)
    # - do the conversions
    SA = gsw.SA_from_SP(SP, p, lon, lat)
    CT = gsw.CT_from_t(SA, IT, p)
    # - add the results to the DataFrame
    df['SA'] = SA
    df['CT'] = CT
    
    # Keep only selected columns.
    cols = ['time', 'lat', 'lon', 'name', 'z', 'CT', 'SA', 'DO (uM)', 'Chl (mg m-3)']
    this_cols = [item for item in cols if item in df.columns]
    df = df[this_cols]
    
    # Generate the cid (cast ID) field.
    df['cid'] = np.nan
    # This assume that a cast can be identified by (i) happening on a unique day,
    # combined with (ii) having a unique name (station name).
    cid = 0
    for time in df.time.unique():
        for name in df.name.unique():
            df.loc[(df.name==name) & (df.time==time),'cid'] = cid
            cid += 1
            
    # Remove one cast that had obviously bad salinity
    if year == 2017:
        for c in df.cid.unique():
            cdf = df[(df.cid==c)]
            if cdf.loc[(cdf.z<-100),'SA'].mean() <23:
                df = df[df.cid != c]
                break
    
    # Renumber cid to be increasing from zero in steps of one.
    df = obs_functions.renumber_cid(df)

    # check that z is monotonic (top to bottom)
    for c in df.cid.unique():
        if not df.loc[(df.cid==c),'z'].is_monotonic_decreasing:
            print('z problem for ' + str(c))
    # result: they are all monotonic, hooray!
            
    print(' - processed %d casts' % ( len(df.cid.unique()) ))
    
    # add a cruise column
    df['cruise'] = None
    
    if len(df) > 0:
        # Save the data
        df.to_pickle(out_fn)
        info_df = obs_functions.make_info_df(df)
        info_df.to_pickle(info_out_fn)
        
print('Total time = %d sec' % (int(Time()-tt0)))
    
    
    
