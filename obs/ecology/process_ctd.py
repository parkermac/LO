"""
Code to process the ecology ctd data.

Takes a couple minutes to run for 2008-2019. Perhaps I could speed it up
by making the first, large, multi-year excel file into a csv.

"""

import pandas as pd
import numpy as np
import gsw
import sys

from lo_tools import Lfun
Ldir = Lfun.Lstart()

# BOTTLE
source = 'ecology'
otype = 'ctd'
in_dir0 = Ldir['data'] / 'obs' / source

testing = True

if testing:
    year_list = [2019]
else:
    year_list = range(2008,2020)

# output location
out_dir = Ldir['LOo'] / 'obs' / source / otype
Lfun.make_dir(out_dir)

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
    
    if year == 2017:
        load_data = True
        ctd_fn = in_dir0 / 'ParkerMacCready2017CTDDataFeb2018.xlsx'
        sheet_name = '2017Provisional_CTDResults'
    elif year == 2018:
        load_data = True
        #ctd_fn = dir0 + 'raw/Parker_2018.xlsx'
        ctd_fn = in_dir0 / 'ParkerMacCready2018CTDDOMar2020.xlsx'
        sheet_name = '2018_CTDDOResults'
    elif year == 2019:
        load_data = True
        ctd_fn = in_dir0 / 'ParkerMacCready2019CTDDataFeb2020.xlsx'
        sheet_name = '2019Provisional_CTDResults'
    else:
        ctd_fn = in_dir0 / 'ParkerMacCready1999-2016CTDDataMay2018.xlsx'
        sheet_name = '1999-2016Finalized_CTDResults'
    
    # DEFAULTS
    date_col_name = 'Date'
    station_col_name = 'Station'
    depth_col_name = 'Depth'
    data_new_names = ['SP', 'IT', 'DO (mg/L)', 'z']
    
    if year in [2017, 2019]:
        data_original_names = ['Salinity', 'Temp', 'DO_raw']
    else:
        data_original_names = ['Salinity', 'Temp', 'DO_adjusted']
    
    # read in the data (all stations, all casts)
    if load_data:
        df0 = pd.read_excel(ctd_fn, sheet_name=sheet_name,
            parse_dates = [date_col_name])
        load_data = False
        
    # add z
    df0['z'] = -df0[depth_col_name]
        
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
    cols = ['time', 'lat', 'lon', 'name', 'z', 'CT', 'SA', 'DO (uM)']
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
            
    # Rework cid to also be increasing from zero in steps of one.
    a = df.cid.values
    au = df.cid.unique() # returns uniques in order
    u_dict = dict(zip(au, np.arange(len(au))))
    b = np.nan * np.ones(len(a))
    for ii in u_dict.keys():
        b[a==ii] = u_dict[ii]
    df['cid'] = b
    
    # Note, the result has casts packed top-to-bottom
        
    print(' - processed %d casts' % ( len(df.cid.unique()) ))
    
    # add a cruise column
    df['cruise'] = None
    
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
    
    
    
