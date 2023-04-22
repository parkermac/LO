"""
Code to process the ecology bottle data.

Takes 3 minutes to run.

"""

import pandas as pd
import numpy as np
import gsw
import sys
from time import time as Time

from lo_tools import Lfun, zfun, obs_functions
Ldir = Lfun.Lstart()

# BOTTLE
source = 'ecology'
otype = 'bottle'
in_dir0 = Ldir['data'] / 'obs' / source

testing = False

if testing:
    year_list = [2017]
else:
    year_list = range(2008,2018)

# output location
out_dir = Ldir['LOo'] / 'obs' / source / otype
Lfun.make_dir(out_dir)

# This is a dict of all the columns after the initial reading.
# We add a value to a key for any variable we want to save
v_dict = {
'ResultID':'',
'Project':'',
'Date':'time',
'LocalTime':'',
'UTCDateTime':'',
'Year':'',
'Month':'',
'Station':'name',
'Niskin':'',
'Niskin Cast Rep':'',
'NomDepth':'',
'Depth_Matching':'',
'Sampling Depth':'d',
'PO4(uM)D':'PO4 (uM)',
'QC PO4_Lab':'',
'SiOH4(uM)D':'SiO4 (uM)',
'QC SiOH4_Lab':'',
'NO3(uM)D':'NO3 (uM)',
'QC NO3_Lab':'',
'NO2(uM)D':'NO2 (uM)',
'QC NO2_Lab':'',
'NH4(uM)D':'NH4 (uM)',
'QC NH4_Lab':'',
'CTD Cast Rep':'',
'Type':'',
'Replicate?':'',
'LabReplicateNumber':'',
'SampleFieldReplicateNumber':'',
'QA':'',
'Filename':'',
'Bottle Number':'',
'Comments':'',
'Analysis Date':'',
}

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
    
    if (year in range(2006,2018)) and load_data:
        in_fn =  in_dir0 / 'bottle_2006_2017_fixed.p'
        df0 = pd.read_pickle(in_fn)
        # This code was used to generate the text used in v_dict above.
        if False:
            for v in df0.columns:
                print("\'%s\':\'\'," % (v))
        # select and rename variables
        df1 = pd.DataFrame()
        for v in df0.columns:
            if v in v_dict.keys():
                if len(v_dict[v]) > 0:
                    df1[v_dict[v]] = df0[v]
        load_data = False # only load the first time
                
    # select one year
    t = pd.DatetimeIndex(df1.time)
    # df1['time'] = t
    df = df1.loc[t.year==year,:].copy()
    
    df['z'] = -df['d']
    
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
                        
    # a little more cleaning up
    df = df.dropna(axis=0, how='all') # drop rows with no good data
    df = df[df.time.notna()] # drop rows with bad time
    df = df.reset_index(drop=True)
        
    # Generate the cid (cast ID) field.
    df['cid'] = np.nan
    # This assume that a cast can be identified by (i) happening on a unique day,
    # combined with (ii) having a unique name (station name).
    cid = 0
    for time in df.time.unique():
        for name in df.name.unique():
            df.loc[(df.name==name) & (df.time==time),'cid'] = cid
            cid += 1
            
    # Renumber cid to be increasing from zero in steps of one.
    df = obs_functions.renumber_cid(df)
    
    # Go to the processed cast file to get SA and CT, and DO. This is by far
    # the most time-consuming part of the processing, but it is required.
    df['SA'] = np.nan
    df['CT'] = np.nan
    df['DO (uM)'] = np.nan
    ctd_dir = Ldir['LOo'] / 'obs' / source / 'ctd'
    ctd_fn = ctd_dir / (ys + '.p')
    cdf = pd.read_pickle(ctd_fn)
    # Assume that having the same time and cast name means these are the same cast.
    for time in df.time.unique():
        for name in df.name.unique():
            cz = cdf.loc[(cdf.name==name) & (cdf.time==time),'z'].to_numpy()
            cSA = cdf.loc[(cdf.name==name) & (cdf.time==time),'SA'].to_numpy()
            cCT = cdf.loc[(cdf.name==name) & (cdf.time==time),'CT'].to_numpy()
            cDO = cdf.loc[(cdf.name==name) & (cdf.time==time),'DO (uM)'].to_numpy()
            bz = df.loc[(df.name==name) & (df.time==time),'z'].to_numpy()
            iz_list = []
            if (len(cz) > 0) and (len(bz) > 0):
                for zz in bz:
                    iz_list.append(zfun.find_nearest_ind(cz,zz))
                df.loc[(df.name==name) & (df.time==time),'SA'] = cSA[iz_list]
                df.loc[(df.name==name) & (df.time==time),'CT'] = cCT[iz_list]
                df.loc[(df.name==name) & (df.time==time),'DO (uM)'] = cDO[iz_list]
            
    # (3) retain only selected variables
    cols = ['cid', 'cruise', 'time', 'lat', 'lon', 'name', 'z',
        'CT', 'SA', 'DO (uM)',
        'NO3 (uM)', 'NO2 (uM)', 'NH4 (uM)', 'PO4 (uM)', 'SiO4 (uM)',
        'TA (uM)', 'DIC (uM)']
    this_cols = [item for item in cols if item in df.columns]
    df = df[this_cols]
        
    print(' - processed %d casts' % ( len(df.cid.unique()) ))

    df['cruise'] = None
    
    if len(df) > 0:
        # Save the data
        df.to_pickle(out_fn)
        info_df = obs_functions.make_info_df(df)
        info_df.to_pickle(info_out_fn)
        
print('Total time = %d sec' % (int(Time()-tt0)))

