"""
Code to process DFO ctd data from the NetCDF version.

Performance: Takes about 18 minutes for all the years (1965-2021)

The processing is very similar to what was done for bottles. So see process_bottle.py
for more extensive comments.
"""
from datetime import datetime, timedelta
import numpy as np
import gsw
import sys
import pandas as pd
from time import time
import xarray as xr
from scipy.stats import binned_statistic

from lo_tools import Lfun, obs_functions
Ldir = Lfun.Lstart()
in_dir = Ldir['data'] / 'obs' / 'dfo1'

testing = False
if testing:
    year_list = [2017]
else:
    year_list = list(range(1965, datetime.now().year + 1))
out_dir = Ldir['LOo'] / 'obs' / 'dfo1' / 'ctd'
Lfun.make_dir(out_dir)

col_dict = {
    'DOXYZZ01': 'DO (ml L-1)',
    'depth':'depth','profile':'cid','latitude':'lat','longitude':'lon',
    'sea_water_practical_salinity':'SP','sea_water_temperature':'TI','time':'tsec'
    }

tt00 = time()
dt_ref = datetime(1970,1,1)
for year in year_list:
    
    if year in range(1965,1993+1):
        name = 'IOS_CTD_Profiles_0.nc'
        if year == 1965 or testing:
            open_ds = True
    elif year in range(1994,2021+1):
        name = 'IOS_CTD_Profiles_1.nc'
        if year == 1994 or testing:
            open_ds = True
    else:
        print('Year out of range')
        
    in_fn = in_dir / name
    if open_ds:
        ds = xr.open_dataset(in_fn, decode_times=False)
        open_ds = False
        
    out_fn = out_dir / (str(year) + '.p')
    info_out_fn = out_dir / ('info_' + str(year) + '.p')
    # limit time range
    dt0 = datetime(year,1,1)
    dt1 = datetime(year+1,1,1)
    tsec0 = (dt0 - dt_ref).total_seconds()
    tsec1 = (dt1 - dt_ref).total_seconds()
    mask = (ds.time.values>=tsec0) & (ds.time.values<tsec1)
    
    if mask.sum() > 0:
    
        # read selected variables into a DataFrame
        df = pd.DataFrame()
        for vn in col_dict.keys():
            df[vn] = ds[vn][mask]
        df = df.rename(col_dict, axis=1)
        
        # Check that there are not two different casts associated with the same profile
        # by looking for large time differences.
        # NOTE: This is the slowest part of all the processing, but it is essential. I have
        # tried a number of things to speed it up, and the solution below gives reasonable
        # performance.
        df = obs_functions.renumber_cid(df) # making the cid's into numbers speeds things up
        cidu = df.cid.unique()
        bad_list = []
        cid_vec = df.cid.to_numpy()
        tsec_vec = df.tsec.to_numpy()
        # tt0 = time()
        # These two methods give identical results, but the first is more complicated
        # so we will use the second
        # 1. first method using binned statistic
        # Cidu = np.concatenate((cidu,np.array([cidu[-1]+1]))) # make these into bin edges
        # amin = binned_statistic(cid_vec, tsec_vec, statistic='min', bins=Cidu).statistic
        # amax = binned_statistic(cid_vec, tsec_vec, statistic='max', bins=Cidu).statistic
        # ii = 0
        # for cid in cidu:
        #     time_diff = amax[ii] - amin[ii]
        #     if time_diff != 0:
        #         #print(' - cid %s has time diff of %d sec' % (str(cid), time_diff))
        #         bad_list.append(cid)
        #     ii += 1
        # 2. second method using boolean mask
        for cid in cidu:
            tvec = tsec_vec[cid_vec==cid]
            time_diff = tvec[-1] - tvec[0]
            if time_diff != 0:
                # print(' - cid %s has time diff of %d sec' % (str(cid), time_diff))
                bad_list.append(cid)
        # print(' -- time to find bad profiles %0.2f' % (time()-tt0))
        # sys.stdout.flush()
        
        # Drop the bad instances
        if len(bad_list) > 0:
            print(' - Dropping %d casts with inconsistent time' % (len(bad_list)))
            for item in bad_list:
                df = df.loc[df.cid!=item,:]
        # RESULT: Mostly the time_diffs are all zero, but we drop about 200 casts over 55 years.
        
        # there are a few fields that had values like 1e37
        for vn in ['depth','TI','SP']:
            df.loc[df[vn]>1e5,vn] = np.nan
        # and a few negative depths and salinities
        df[df['depth']<=0] = np.nan
        df[df['SP']<=0] = np.nan

        # add z and drop rows missing any fundamental variables
        df['z'] = -df.depth
        df = df.dropna(subset=['z','lon','lat'])

        # convert to standard units
        p = gsw.p_from_z(df.z.to_numpy(), df.lat) # no clue why .to_numpy() is needed
        df['SA'] = gsw.SA_from_SP(df.SP, p, df.lon, df.lat)
        df['CT'] = gsw.CT_from_t(df.SA, df.TI, p)
        df['DO (uM)'] = df['DO (ml L-1)'].to_numpy() * 1.42903 * 1000 / 32

        # NOTE: convenient ways to look for nans or out of range values are
        # df.isna().sum() and df.max(), df.min()
        tsec = df.tsec.to_numpy()
        time_list = []
        for item in tsec:
            time_list.append(dt_ref + timedelta(seconds=item))
        df['time'] = time_list

        # clean up columns
        cols = ['cid', 'lon', 'lat', 'time', 'z','SA', 'CT', 'DO (uM)']
        this_cols = [item for item in cols if item in df.columns]
        df = df[this_cols]

        df['name'] = None
        df['cruise'] = None

        # Renumber cid to be increasing from zero in steps of one.
        df = obs_functions.renumber_cid(df)
        
        if len(df) > 0:
            # Save the data
            df.to_pickle(out_fn)
            info_df = obs_functions.make_info_df(df)
            info_df.to_pickle(info_out_fn)
            print('%d: %d casts' % (year,int(info_df.index.max())))
        else:
            print('%d: 0 casts' % (year))
        
    else:
        print('%d: 0 casts' % (year))

print('Elapsed time for processing = %0.1f sec' % (time()-tt00))


# Code to explore contents of the NetCDF files.
if False:
    name_list = ['IOS_CTD_Profiles_0.nc','IOS_CTD_Profiles_1.nc']
    for name in name_list:
        fn = in_dir / name
        print('\n'+fn.name)
        ds = xr.open_dataset(fn,decode_times=False)

        for vn in ds.data_vars:
            try:
                units = ds[vn].units
            except:
                units = ''
            try:
                long_name = ds[vn].long_name
            except:
                long_name = ''
            print('%s: %s [%s]' % (vn,long_name,units))
        
        # time range
        dt_ref = datetime(1970,1,1)
        dt0 = dt_ref + timedelta(seconds=float(ds.time[0].values))
        dt1 = dt_ref + timedelta(seconds=float(ds.time[-1].values))
        print(dt0)
        print(dt1)
            
"""
RESULT:
IOS_CTD_Profiles_0.nc
CNDCST01: Sea Water Electrical Conductivity [S/m]
DOXMZZ01: Oxygen concentration [umol/kg]
*DOXYZZ01: Oxygen concentration [mL/L]
PRESPR01: Pressure [decibar]
PSALST01: Sea Water Practical Salinity [PSS-78]
PSALST02: Sea Water Practical Salinity [PSS-78]
SSALST01: Sea Water Salinity [PPT]
TEMPS601: Sea Water Temperature [degC]
TEMPS602: Sea Water Temperature [degC]
TEMPS901: Sea Water Temperature [degC]
TEMPS902: Sea Water Temperature [degC]
TEMPST01: Sea Water Temperature [degC]
agency:  []
country:  []
*depth: Depth [m]
event_number:  []
filename:  []
geographic_area:  []
instrument_model:  []
instrument_serial_number:  []
instrument_type:  []
*latitude: Latitude [degrees_north]
*longitude: Longitude [degrees_east]
mission_id:  []
platform:  []
*profile: Profile ID []
project:  []
scientist:  []
*sea_water_practical_salinity: Sea Water Practical Salinity [PSS-78]
sea_water_pressure: Pressure [dbar]
*sea_water_temperature: Sea Water Temperature [degC]
*time: Time [seconds since 1970-01-01T00:00:00Z]

Time ranges:

_Profiles_0:
1965-07-05 21:50:00
1993-04-08 22:22:11

_Profiles_1:
1994-02-08 18:57:26
2021-02-03 16:37:06

"""