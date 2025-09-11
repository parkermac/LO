"""
Code to process the King County Water Quality CTD data for Puget Sound.

To process King County QCed dataset received by Dakota Mascarenas from Greg Ikeda, King County, via file transfer on 2024/04/05.

Station information downloaded by Dakota Mascarenas from: https://data.kingcounty.gov/Environment-Waste-Management/WLRD-Sites/wbhs-bbzf

Initial author date: 2023/07/13

Finalized for public use: 2025/09/04

Written by: Dakota Mascarenas

"""

import pandas as pd
import numpy as np
import gsw

import glob

from lo_tools import Lfun, obs_functions
Ldir = Lfun.Lstart()


# source location
source = 'kc'
otype = 'ctd'
in_dir0 = Ldir['data'] / 'obs' / source 
year_list = range(1998,2025)

# output location
out_dir = Ldir['LOo'] / 'obs' / source / otype
Lfun.make_dir(out_dir)

# Load individual data sets.
fn = glob.glob(str(in_dir0) + '/' + otype + '/*.csv')

# Initial clean.
big_df_raw = pd.DataFrame()
for f in fn:
    raw = pd.read_csv(f, encoding='cp1252')
    if 'ï»¿Locator' in raw.columns:
        raw = raw.rename(columns={'ï»¿Locator':'Locator'})
    if big_df_raw.empty:
        big_df_raw = raw
    else:
        big_df_raw = pd.concat([big_df_raw, raw])
    
# Load station data.
sta_df = pd.read_csv(in_dir0 / 'WLRD_Sites_March2024.csv')

# Merge station data.
big_df = big_df_raw.merge(sta_df[['Locator','Latitude', 'Longitude']], on = 'Locator', how='left')

# Use only downcasts.
big_df_use0 = big_df[big_df['Updown'] == 'Down']

# Create dictionary for important variable and column names.
v_dict = {'Chlorophyll, Field (mg/m^3)':'Chl (mg m-3)',
          'Dissolved Oxygen, Field (mg/l ws=2)': 'DO (mg -L)',
          'Nitrite + Nitrate Nitrogen, Field (mg/L)': 'NO3 (mg -L)', # measured together assuming 0 NO2
          'Salinity, Field (PSS)':'SP',
          'Sample Temperature, Field (deg C)':'IT'
          } # not dealing with light/PAR right now...can be added later
v_dict_use = {}
for v in v_dict.keys():
    if len(v_dict[v]) > 0:
        v_dict_use[v] = v_dict[v]
v_list = np.array(list(v_dict_use.keys()))

# Clean column names.
big_df_use6 = big_df_use0.copy()
big_df_use6['time'] = pd.DatetimeIndex(big_df_use6['Sampledate'])

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
v_dict['Depth'] = 'z' # will be converted to negative later in script
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
    if 'NH4 (mg -L)' in df.columns:
        df['NH4 (uM)'] = (1000/18) * df['NH4 (mg -L)']
    if 'NO3 (mg -L)' in df.columns:
        df['NO3 (uM)'] = (1000/62) * df['NO3 (mg -L)']
    if 'SiO4 (mg -L)' in df.columns:
        df['SiO4 (uM)'] = (1000/92) * df['SiO4 (mg -L)']
    if 'PO4 (mg -L)' in df.columns:
        df['PO4 (uM)'] = (1000/95) * df['PO4 (mg -L)']
    if 'Chl (ug -L)' in df.columns:
        df['Chl (mg m-3)'] = df['Chl (ug -L)']
    for vn in ['TA','DIC']:
        if (vn+' (umol -kg)') in df.columns:
            df[vn+' (uM)'] = (rho/1000) * df[vn+' (umol -kg)']
    # retain only selected variables
    df['cruise'] = ''
    cols = ['cid', 'time', 'lat', 'lon', 'z', 'cruise', 'name',
        'CT', 'SA', 'DO (uM)',
        'NO3 (uM)', 'NH4 (uM)', 'PO4 (uM)', 'SiO4 (uM)', #removed NO2...
        'TA (uM)', 'DIC (uM)', 'Chl (mg m-3)']
    this_cols = [item for item in cols if item in df.columns]
    df = df[this_cols]
    # save
    print(' - processed %d casts' % ( len(df.cid.unique()) ))
    if len(df) > 0:
        # Save the data
        df.to_pickle(out_fn)
        info_df = obs_functions.make_info_df(df)
        info_df.to_pickle(info_out_fn)