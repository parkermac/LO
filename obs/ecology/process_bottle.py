"""
Code to process the ecology bottle data.

"""

import pandas as pd
import numpy as np
import gsw
import sys

from lo_tools import Lfun
Ldir = Lfun.Lstart()

# BOTTLE
source = 'ecology'
otype = 'bottle'
in_dir0 = Ldir['data'] / 'obs' / source

testing = True

if testing:
    year_list = [2017]
else:
    year_list = range(2008,2019)

# output location
out_dir = Ldir['LOo'] / 'obs' / source / otype
Lfun.make_dir(out_dir)

# This is a dict of all the columns after the initial reading.
# We add values to a key for any variable we want to save
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
'NH4(uM)D':'NH4(uM)',
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
# file, but there are station names, with loccations here:
sta_fn = in_dir0 / 'ParkerMacCreadyCoreStationInfoFeb2018.xlsx'
sta_df = pd.read_excel(sta_fn, index_col='Station')
xx = sta_df['Long_NAD83 (deg / dec_min)'].values
yy = sta_df['Lat_NAD83 (deg / dec_min)'].values
lon = [-(float(x.split()[0]) + float(x.split()[1])/60) for x in xx]
lat = [(float(y.split()[0]) + float(y.split()[1])/60) for y in yy]
sta_df['lon'] = lon
sta_df['lat'] = lat

load_data = True
for year in year_list:
    ys = str(year)
    print('\n'+ys)
    
    # name output files
    out_fn = out_dir / (ys + '.p')
    info_out_fn = out_dir / ('info_' + ys + '.p')
    
    if year in range(2006,2018):
        in_fn =  in_dir0 / 'Parker_2006-2017_Nutrients.xlsx'
        df0 = pd.read_excel(in_fn, sheet_name='2006-2017')
        load_data = False # only load the first time
        
    # for v in df0.columns:
    #     print("\'%s\':\'\'," % (v))

    # select and rename variables
    df1 = pd.DataFrame()
    for v in df0.columns:
        if v in v_dict.keys():
            if len(v_dict[v]) > 0:
                df1[v_dict[v]] = df0[v]
                
    # select one year
    t = pd.DatetimeIndex(df1.time)
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
    
    sys.exit()
                
    # missing data is -999
    df[df==-999] = np.nan
    
    # a little more cleaning up
    df = df.dropna(axis=0, how='all') # drop rows with no good data
    df = df[df.time.notna()] # drop rows with bad time
    df = df.reset_index(drop=True)
    
    # Now proceed with the processing to get a single DataFrame for the year.
    
    # add the "cid" (cast ID) column
    #
    # Note that we will save the field "name" for station number, since this dataset has
    # repeat stations which is helpful for plotting sections. Then we will generate our own
    # cid, a unique one for each cast, being careful to keep them unique for the collection
    # of cruises in this year, even though a station may be repeated on all cruises.
    #
    # We will also save the field "cruise" as a convenient way to select a collection of
    # casts.
    df['cid'] = np.nan
    cid = 0
    for cruise in df.cruise.unique():
        for name in df.name.unique():
            df.loc[(df.name==name) & (df.cruise==cruise),'cid'] = cid
            cid += 1
    for cid in df.cid.unique():
        # Check that there are not two different casts associated with the same station
        # by looking for large time differences. Pretty ad hoc, but it works.
        time_diff = df[df.cid==cid].time.values[-1] - df[df.cid==cid].time.values[0]
        time_diff = pd.to_timedelta(time_diff)
        if time_diff.days > 1 or time_diff.days < -1:
            cruise = df[df.cid==cid].cruise.values[0]
            name = df[df.cid==cid].name.values[0]
            print('Cruise: %s, Station %s has time diff of %d days' % (cruise, str(name), time_diff.days))
            # copy in just the first cast at this repeated station
            dff = df[df.cid==cid].copy()
            dfft = dff.time.values
            Dfft = pd.to_timedelta(dfft - dfft[0])
            dff = dff[Dfft.days==0]
            print('  - length of df before removing repeat cast at this station: %d' % (len(df)))
            df = df[df.cid != cid]
            df = pd.concat((df,dff))
            print('  - length of df before removing repeat cast at this station: %d' % (len(df)))
        # Force certain fields to be the same throughout the cast.
        df.loc[df.cid==cid,'lon'] = df[df.cid==cid].lon.values[0]
        df.loc[df.cid==cid,'lat'] = df[df.cid==cid].lat.values[0]
        df.loc[df.cid==cid,'time'] = df[df.cid==cid].time.values[0]
                    
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
        
    # (3) retain only selected variables
    cols = ['cid', 'cruise', 'time', 'lat', 'lon', 'name', 'z',
        'CT', 'SA', 'DO (uM)',
        'NO3 (uM)', 'NO2 (uM)', 'NH4 (uM)', 'PO4 (uM)', 'SiO4 (uM)',
        'TA (uM)', 'DIC (uM)']
    this_cols = [item for item in cols if item in df.columns]
    df = df[this_cols]
        
    print(' - processed %d casts' % ( len(df.cid.unique()) ))
    cid0 = df.cid.max() + 1
        
    # Sort the result by time, and sort each cast to be bottom to top
    df = df.sort_values(['time','z'], ignore_index=True)
    
    # Rework cid to also be increasing in time
    a = df[['time','cid']].copy()
    a['cid_alt'] = np.nan
    ii = 0
    for t in a.time.unique():
        a.loc[a.time==t,'cid_alt'] = ii
        ii += 1
    df['cid'] = a['cid_alt'].copy()
    
    if len(df) > 0:
        # Save the data
        df.to_pickle(out_fn)

        # Also pull out a dateframe with station info to use for model cast extractions.
        ind = df.cid.unique()
        col_list = ['lon','lat','time','name','cruise']
        info_df = pd.DataFrame(index=ind, columns=col_list)
        for cid in df.cid.unique():
            info_df.loc[cid,col_list] = df.loc[df.cid==cid,col_list].iloc[0,:]
        info_df.index.name = 'cid'
        info_df['time'] = pd.to_datetime(info_df['time'])
        info_df.to_pickle(info_out_fn)
