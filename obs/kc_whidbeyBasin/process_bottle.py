"""
Code to process the King County Water Quality BOTTLE data for Whidbey Basin in Puget Sound.

To process publically-accessible data downloaded by Dakota Mascarenas from: https://data.kingcounty.gov/Environment-Waste-Management/Whidbey-Bottle-Data/vuu8-t6kc

Station information downloaded by Dakota Mascarenas from: https://data.kingcounty.gov/Environment-Waste-Management/WLRD-Sites/wbhs-bbzf

NOTE: This data is replicated in the folder kc bottle data, but is included here since this is the format that King County maintains.

Initial author date: 2026/04/21

Finalized for public use: 2026/04/21

Written by: Dakota Mascarenas

"""

import pandas as pd
import numpy as np
import gsw

from lo_tools import Lfun, obs_functions
Ldir = Lfun.Lstart()

# source location
source = 'kc_whidbeyBasin'
otype = 'bottle'
in_dir0 = Ldir['data'] / 'obs' / source / otype
year_list = range(2022,2027)

# output location
out_dir = Ldir['LOo'] / 'obs' / source / otype
Lfun.make_dir(out_dir)

# Load big data set and stations.
big_df_raw = pd.read_csv(in_dir0 / 'Whidbey_Bottle_Data_20260421.csv', low_memory=False)
sta_df = pd.read_csv(Ldir['data'] / 'obs' / source / 'WLRD_Sites_20260421.csv')

# Merge station data.
big_df = big_df_raw.merge(sta_df[['Locator','Latitude', 'Longitude']], on='Locator', how='left')

# Select only marine offshore stations.
big_df_use0 = big_df[big_df['Site Type'] == 'Marine Offshore']

# Create dictionary for important variable and parameter names.
cols_all = big_df_use0['Parameter'].unique()
v_dict = {col:'' for col in cols_all}
v_dict['Temperature'] = 'IT'
v_dict['Salinity'] = 'SP'
v_dict['Dissolved Oxygen'] = 'DO (mg -L)'
v_dict['Nitrite + Nitrate Nitrogen'] = 'NO3 (mg -L)' # measured together assuming a 0 NO2
v_dict['Ammonia Nitrogen'] = 'NH4 (mg -L)'
v_dict['Orthophosphate Phosphorus'] = 'PO4 (mg -L)'
v_dict['Silica'] = 'SiO4 (mg -L)'
v_dict['Chlorophyll a'] = 'Chl (ug -L)'
v_dict_use = {}
for v in v_dict.keys():
    if len(v_dict[v]) > 0:
        v_dict_use[v] = v_dict[v]
v_list = np.array(list(v_dict_use.keys())) #redundant but fine

# Keep only selected parameters and clean repeated records.
big_df_use1 = big_df_use0[big_df_use0['Parameter'].isin(v_list)].copy()
big_df_use1['Value'] = pd.to_numeric(big_df_use1['Value'], errors='coerce')
big_df_use2 = big_df_use1[['Sample ID', 'Profile ID', 'Collect DateTime', 'Depth (m)', 'Parameter',
                           'Value', 'Replicates', 'Latitude', 'Longitude', 'Locator']].copy()
replicates = big_df_use2['Replicates'].dropna().astype(str).unique()
big_df_use3 = big_df_use2[~big_df_use2['Sample ID'].astype(str).isin(replicates)]
big_df_use4 = big_df_use3[['Profile ID', 'Collect DateTime', 'Depth (m)', 'Parameter',
                           'Value', 'Latitude', 'Longitude', 'Locator']].copy()

# Pivot from long to wide format.
big_df_use5 = big_df_use4.pivot_table(index=['Profile ID', 'Collect DateTime', 'Depth (m)',
                                             'Latitude', 'Longitude', 'Locator'],
                                      columns='Parameter', values='Value').reset_index()

# Clean column names.
big_df_use6 = big_df_use5.copy()
big_df_use6['time'] = pd.to_datetime(big_df_use6['Collect DateTime'])
# Convert from PST (UTC-8, no daylight savings adjustment) to UTC.
big_df_use6['time'] = big_df_use6['time'].dt.tz_localize('Etc/GMT+8').dt.tz_convert('UTC')

# Create unique cast IDs (cid).
big_df_use6['cid'] = np.nan
big_df_use7 = big_df_use6.copy()
c = 0
for pid in big_df_use7['Profile ID'].unique(): # profile ID is unique identifier
    big_df_use7.loc[big_df_use7['Profile ID'] == pid, 'cid'] = c
    c += 1

# Rename some columns in variable dictionary.
v_dict['cid'] = 'cid'
v_dict['time'] = 'time'
v_dict['Latitude'] = 'lat'
v_dict['Longitude'] = 'lon'
v_dict['Depth (m)'] = 'z' # will be converted to negative later in script
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
    df = df[df.z.notna()] # drop rows with missing depth
    df = df.reset_index(drop=True)
    df['z'] = df['z'] * -1 # IMPORTANT!!!!!! - from above!

    # Keep rows where required state variables are available for TEOS-10 conversion.
    if ('SP' in df.columns) and ('IT' in df.columns):
        good = df['SP'].notna() & df['IT'].notna() & df['lat'].notna() & df['lon'].notna()
        df = df.loc[good,:].copy()

    if len(df) == 0:
        print(' - processed 0 casts')
        continue

    SP = df.SP.to_numpy()
    IT = df.IT.to_numpy()
    z = df.z.to_numpy()
    lon = df.lon.to_numpy()
    lat = df.lat.to_numpy()
    # do the gsw conversions
    p = gsw.p_from_z(z, lat)
    SA = gsw.SA_from_SP(SP, p, lon, lat)
    CT = gsw.CT_from_t(SA, IT, p)
    # add the results to the DataFrame
    df['SA'] = SA
    df['CT'] = CT
    rho = gsw.rho(SA, CT, p)
    # unit conversions
    if 'DO (mg -L)' in df.columns:
        df['DO (uM)'] = (1000/32) * df['DO (mg -L)']
    if 'NH4 (mg -L)' in df.columns:
        df['NH4 (uM)'] = (1000/14) * df['NH4 (mg -L)']  # N atomic weight
    if 'NO3 (mg -L)' in df.columns:
        df['NO3 (uM)'] = (1000/14) * df['NO3 (mg -L)']  # N atomic weight
    if 'SiO4 (mg -L)' in df.columns:
        df['SiO4 (uM)'] = (1000/28.0855) * df['SiO4 (mg -L)']  # Si atomic weight
    if 'PO4 (mg -L)' in df.columns:
        df['PO4 (uM)'] = (1000/30.973762) * df['PO4 (mg -L)']  # P atomic weight
    if 'Chl (ug -L)' in df.columns:
        df['Chl (mg m-3)'] = df['Chl (ug -L)']

    # retain only selected variables
    df['cruise'] = ''
    cols = ['cid', 'time', 'lat', 'lon', 'z', 'cruise', 'name',
        'CT', 'SA', 'DO (uM)',
        'NO3 (uM)', 'NH4 (uM)', 'PO4 (uM)', 'SiO4 (uM)',
        'Chl (mg m-3)']
    this_cols = [item for item in cols if item in df.columns]
    df = df[this_cols]
    # save
    print(' - processed %d casts' % ( len(df.cid.unique()) ))
    if len(df) > 0:
        # Save the data
        df.to_pickle(out_fn)
        info_df = obs_functions.make_info_df(df)
        info_df.to_pickle(info_out_fn)
