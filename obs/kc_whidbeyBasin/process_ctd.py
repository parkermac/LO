"""
Code to process the King County Water Quality CTD data for Whidbey Basin in Puget Sound.

To process publically-accessible data downloaded by Dakota Mascarenas from:  https://data.kingcounty.gov/Environment-Waste-Management/Whidbey-Basin-CTD-Casts/uz4m-4d96

Station information downloaded by Dakota Mascarenas from: https://data.kingcounty.gov/Environment-Waste-Management/WLRD-Sites/wbhs-bbzf

Initial author date: 2024/04/04

Finalized for public use: 2025/09/05

Last updated: 2025/11/24 to correct conversion unit error for nitrogen species

Written by: Dakota Mascarenas

"""

import pandas as pd
import numpy as np
import gsw

from lo_tools import Lfun, obs_functions
Ldir = Lfun.Lstart()

# source location
source = 'kc_whidbeyBasin'
otype = 'ctd'
in_dir0 = Ldir['data'] / 'obs' / source / otype
year_list = range(2022,2025)

# output location
out_dir = Ldir['LOo'] / 'obs' / source / otype
Lfun.make_dir(out_dir)

# Load big data set and stations.
big_df_raw = pd.read_csv(in_dir0/ 'Whidbey_Basin_CTD_Casts_April2024.csv')
sta_df = pd.read_csv(in_dir0 / 'WLRD_Sites_March2024.csv')

# Merge station data.
big_df = big_df_raw.merge(sta_df[['Locator','Latitude', 'Longitude']], on = 'Locator', how='left')

# Use only downcasts.
big_df_use0 = big_df[big_df['Up Down'] == 'Down']

# Create dictionary for important variable and column names.
v_dict = {'Chlorophyll (µg/L)':'Chl (ug -L)',
          'Dissolved Oxygen (mg/L)':'DO (mg -L)',
          'Nitrate + Nitrite (mg N/L)': 'NO3 (mg -L)', # measured together assuming a 0 NO2
          'Salinity (PSU)':'SP', 
          'Temperature (°C)':'IT' 
          } # not dealing with light/PAR right now...
v_dict_use = {}
for v in v_dict.keys():
    if len(v_dict[v]) > 0:
        v_dict_use[v] = v_dict[v]
v_list = np.array(list(v_dict_use.keys())) #redundant but fine
        
# Clean column names.
big_df_use6 = big_df_use0.copy()
big_df_use6['time'] = pd.DatetimeIndex(big_df_use6['Sample Date'])

# Create unique cast IDs (cid).
big_df_use6['cid'] = np.nan
big_df_use7 = big_df_use6.copy()
big_df_use7['unique_date_location'] = big_df_use7['Locator'] + (big_df_use7['time'].dt.year).astype(str) + (big_df_use7['time'].dt.month).astype(str) + (big_df_use7['time'].dt.day).astype(str)
c = 0
for pid in big_df_use7['unique_date_location'].unique(): # profile ID is unique identifier
    big_df_use7.loc[big_df_use7['unique_date_location'] == pid, 'cid'] = c
    c+=1
    
# Rename some columns in variable dictionary.
v_dict['cid'] = 'cid'
v_dict['time'] = 'time'
v_dict['Latitude'] = 'lat'
v_dict['Longitude'] = 'lon'
v_dict['Depth (meters)'] = 'z' # will be converted to negative later in script
v_dict['Locator'] = 'name'

# Loop through to rename variables and columns, clean the dataset, and produce output dataframes.
df0 = big_df_use7.copy()
for year in year_list:
    ys = str(year)
    print('\n'+ys)
    out_fn = out_dir / (ys + '.p')
    info_out_fn = out_dir / ('info_' + ys + '.p')
    t = pd.DatetimeIndex(df0.time)
    df1 = df0.loc[t.year==year,:].copy()   
    # select and rename variables
    df = pd.DataFrame()
    for v in df1.columns:
        if v in v_dict.keys():
            if len(v_dict[v]) > 0:
                df[v_dict[v]] = df1[v]
    # a little more cleaning up
    df = df.dropna(axis=0, how='all') # drop rows with no good data
    df = df[df.time.notna()] # drop rows with bad time
    df = df.reset_index(drop=True)
    df['z'] = df['z']*-1 # IMPORTANT!!!!!! - from above!
    SP = df.SP.to_numpy()
    IT = df.IT.to_numpy()
    z= df.z.to_numpy()
    lon = df.lon.to_numpy()
    lat = df.lat.to_numpy()
    # do the gsw conversions
    p = gsw.p_from_z(z, lat)
    SA = gsw.SA_from_SP(SP, p, lon, lat)
    CT = gsw.CT_from_t(SA, IT, p)
    # add the results to the DataFrame
    df['SA'] = SA
    df['CT'] = CT
    rho = gsw.rho(SA,CT,p)
    # unit conversions
    if 'DO (mg -L)' in df.columns:
        df['DO (uM)'] = (1000/32) * df['DO (mg -L)']
    if 'NO3 (mg -L)' in df.columns:
        df['NO3 (uM)'] = (1000/14) * df['NO3 (mg -L)']
    if 'Chl (ug -L)' in df.columns:
        df['Chl (mg m-3)'] = df['Chl (ug -L)']
    # retain only selected variables             
    df['cruise'] = ''
    cols = ['cid', 'time', 'lat', 'lon', 'z', 'cruise', 'name',
        'CT', 'SA', 'DO (uM)',
        'NO3 (uM)', 'Chl (mg m-3)']
    this_cols = [item for item in cols if item in df.columns]
    df = df[this_cols]
    # save
    print(' - processed %d casts' % ( len(df.cid.unique()) ))
    if len(df) > 0:
        # Save the data
        df.to_pickle(out_fn)
        info_df = obs_functions.make_info_df(df)
        info_df.to_pickle(info_out_fn)