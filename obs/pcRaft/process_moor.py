"""
Code to process Penn Cove mussel raft sonde data as discussed in Roberts & Carrington (2023): https://doi.org/10.1016/j.jembe.2023.151927

Data received via email on 2026/01/02 and 2026/02/24 from Emily Carrington to Dakota Mascarenas, with metadata clarification also on 2026/02/24.

Initial author date: 2026/02/24

Finalized for public use: 2026/04/13

Written by: Dakota Mascarenas

"""
# %%

import pandas as pd
import numpy as np
import gsw
import xarray as xr

from lo_tools import Lfun, obs_functions
Ldir = Lfun.Lstart()

# source location
source = 'pcRaft'
otype = 'moor'
in_dir0 = Ldir['data'] / 'obs' / source
year_list = range(2014,2020)

# output locations for both output types
out_dir = Ldir['LOo'] / 'obs' / source / otype 
Lfun.make_dir(out_dir)

sn_loc_dict = {
    'raft_main': [-122.705667, 48.220861],
    'raft_hiRes': [-122.705667, 48.220861],
}


# Load big data set and stations.
big_df = pd.read_excel(in_dir0/ 'PennCove-mussel-raft-sonde-data-2014-2019.xlsx', sheet_name='cleaned_data')
big_df['time'] = pd.DatetimeIndex(big_df['TIMESTAMP'])
# Convert from PST (UTC-8) to UTC
big_df['time'] = big_df['time'].dt.tz_localize('Etc/GMT+8').dt.tz_convert('UTC')
big_df['station'] = 'raft_main'
big_df = big_df[['station', 'time', 'depth', 'temp', 'sal', 'chl', 'pH', 'DO']]

# Load in July-December 2017 high resolution sampling (10 minute frequency to be selected only on even hour marks).
little_df = pd.DataFrame()
for depth in [0.5, 1, 2, 3, 4, 5, 6, 7]:
    little_df_temp = pd.read_csv(in_dir0/ ('B8_' + str(depth) + 'm.csv'))
    little_df_temp['depth'] = depth
    little_df_temp = little_df_temp.set_axis(['#', 'TIME', 'temp', 'depth'], axis=1)
    little_df = pd.concat([little_df, little_df_temp])
little_df['time'] = pd.DatetimeIndex(little_df['TIME'])
# Convert from PST (UTC-8) to UTC
little_df['time'] = little_df['time'].dt.tz_localize('Etc/GMT+8').dt.tz_convert('UTC')
little_df = little_df[['time', 'depth', 'temp']]
little_df = (
    little_df
    .set_index('time')
    .groupby('depth')['temp']
    .resample('1h')
    .mean(numeric_only=True)
    .reset_index()
)
little_df['station'] = 'raft_hiRes'

# %%

# Combine dataframes.
big_df_use = pd.concat([big_df, little_df])

# Create dictionary for important variable and column names.
v_dict = {'chl':'Chl (ug -L)',
          'DO':'DO (mg -L)',
          'pH':'PH',
          'sal':'SP', 
          'temp':'IT' 
          }
v_dict_use = {}
for v in v_dict.keys():
    if len(v_dict[v]) > 0:
        v_dict_use[v] = v_dict[v]
v_list = np.array(list(v_dict_use.keys())) #redundant but fine
        
# Clean column names and add coordinates.
df0 = big_df_use.copy()
    

# Rename some columns in variable dictionary.
v_dict['time'] = 'time'
v_dict['depth'] = 'z' # will be converted to negative later in script
# %%
# Loop through to rename variables and columns, clean the dataset, and produce output dataframes.

for sn in ['raft_main', 'raft_hiRes']:
    out_fn = out_dir / (sn + '.nc')
    df1 = df0[df0['station'] == sn].copy()
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
    df['lon'] = sn_loc_dict[sn][0]
    df['lat'] = sn_loc_dict[sn][1]
    IT = df.IT.to_numpy()
    z= df.z.to_numpy()
    lon = df.lon.to_numpy()
    lat = df.lat.to_numpy()
    if sn == 'raft_main':
        SP = df.SP.to_numpy()
        # do the gsw conversions
        p = gsw.p_from_z(z, lat)
        SA = gsw.SA_from_SP(SP, p, lon, lat)
        CT = gsw.CT_from_t(SA, IT, p)
        # add the results to the DataFrame
        df['SA'] = SA
        df['CT'] = CT
    elif sn == 'raft_hiRes':
        # For each depth in raft_hiRes, find the closest depth in raft_main and use its salinity
        df_main = df0[df0['station'] == 'raft_main'].copy()
        # Rename depth to z to match the processed dataframe
        df_main['z'] = df_main['depth'] * -1
        # Remove timezone info for comparison (times are already in UTC)
        df_main['time'] = df_main['time'].dt.tz_localize(None)
        df['time'] = df['time'].dt.tz_localize(None)
        SP = np.zeros_like(IT)
        
        for idx, (time_val, depth_val) in enumerate(zip(df['time'].values, z)):
            # Find matching time in raft_main data
            main_time_match = df_main[df_main['time'] == time_val]
            
            if len(main_time_match) > 0:
                # Find closest depth in raft_main for this time
                available_depths = main_time_match['z'].values
                closest_depth_idx = np.argmin(np.abs(available_depths - depth_val))
                SP[idx] = main_time_match.iloc[closest_depth_idx]['sal']
            else:
                # If no exact time match, interpolate from nearest time
                time_diff = np.abs(df_main['time'].values - time_val)
                nearest_time_idx = np.argmin(time_diff)
                nearest_time = df_main.iloc[nearest_time_idx]['time']
                
                # Find closest depth at that time
                main_at_time = df_main[df_main['time'] == nearest_time]
                available_depths = main_at_time['z'].values
                closest_depth_idx = np.argmin(np.abs(available_depths - depth_val))
                SP[idx] = main_at_time.iloc[closest_depth_idx]['sal']
        
        # do the gsw conversions using matched SP data
        p = gsw.p_from_z(z, lat)
        SA = gsw.SA_from_SP(SP, p, lon, lat)
        CT = gsw.CT_from_t(SA, IT, p)
        # add the results to the DataFrame
        df['SA'] = np.nan # set SA to NaN since it's not directly measured in raft_hiRes
        df['CT'] = CT

    # unit conversions
    if 'DO (mg -L)' in df.columns:
        df['DO (uM)'] = (1000/32) * df['DO (mg -L)']
    if 'NO3 (mg -L)' in df.columns:
        df['NO3 (uM)'] = (1000/62) * df['NO3 (mg -L)']
    if 'Chl (ug -L)' in df.columns:
        df['Chl (mg m-3)'] = df['Chl (ug -L)']
    
    # retain only selected variables
    cols = ['time', 'lat', 'lon', 'z',
        'SA', 'CT', 'DO (uM)',
        'PH', 'Chl (mg m-3)']
    
    # Select and filter columns that exist in the dataframe
    this_cols = [item for item in cols if item in df.columns]
    df = df[this_cols]
    
    # Remove duplicate records
    if 'z' in df.columns:
        df = df[~df.duplicated(subset=["time", "z"], keep=False)]
    
    # Save to NetCDF in correct format (hourly)
    if len(df) > 0:
        print(f'Processing {sn}: {len(df)} records')
        
        # Convert to numpy arrays for xarray Dataset creation
        time_arr = df['time'].unique()
        # Remove timezone info for NetCDF compatibility (times are already in UTC)
        time_arr = pd.DatetimeIndex(time_arr).tz_localize(None)
        z_arr = np.sort(df['z'].unique())
        NT = len(time_arr)
        NZ = len(z_arr)
        
        # Create a 2D array for each variable (time, z)
        def create_2d_array(var_name):
            """Create 2D array (time, z) from dataframe"""
            arr = np.full((NT, NZ), np.nan)
            # Remove timezone for matching
            df_time_naive = df.copy()
            df_time_naive['time'] = df_time_naive['time'].dt.tz_localize(None)
            for t_idx, t in enumerate(time_arr):
                time_data = df_time_naive[df_time_naive['time'] == t]
                for z_idx, z in enumerate(z_arr):
                    match = time_data[time_data['z'] == z]
                    if len(match) > 0:
                        arr[t_idx, z_idx] = match[var_name].values[0]
            return arr
        
        # Create coordinate arrays
        coords = {'time': ('time', time_arr), 'z': ('z', z_arr)}
        
        # Create data_vars dict with variables that exist
        data_vars = {}
        data_vars['CT'] = (('time', 'z'), create_2d_array('CT'))
        
        if 'SA' in df.columns:
            data_vars['SA'] = (('time', 'z'), create_2d_array('SA'))
        if 'DO (uM)' in df.columns:
            data_vars['DO (uM)'] = (('time', 'z'), create_2d_array('DO (uM)'))
        if 'PH' in df.columns:
            data_vars['PH'] = (('time', 'z'), create_2d_array('PH'))
        if 'Chl (mg m-3)' in df.columns:
            data_vars['Chl (mg m-3)'] = (('time', 'z'), create_2d_array('Chl (mg m-3)'))
        
        # Create Dataset in correct format
        lon = float(df['lon'].iloc[0])
        lat = float(df['lat'].iloc[0])
        ds = xr.Dataset(coords=coords, 
                       attrs={'Station Name': sn, 'lon': lon, 'lat': lat})
        
        # Add variables with attributes
        for var_name, (dims, data) in data_vars.items():
            ds[var_name] = xr.DataArray(data, dims=dims)
        
        # Set variable attributes
        ds['CT'].attrs = {'units': 'degC', 'long_name': 'Conservative Temperature'}
        if 'SA' in ds.data_vars:
            ds['SA'].attrs = {'units': 'g kg-1', 'long_name': 'Absolute Salinity'}
        if 'DO (uM)' in ds.data_vars:
            ds['DO (uM)'].attrs = {'units': 'uM', 'long_name': 'Dissolved Oxygen'}
        if 'PH' in ds.data_vars:
            ds['PH'].attrs = {'units': 'NBS scale', 'long_name': 'pH'}
        if 'Chl (mg m-3)' in ds.data_vars:
            ds['Chl (mg m-3)'].attrs = {'units': 'mg m-3', 'long_name': 'Chlorophyll'}
        
        ds.to_netcdf(out_fn)
        ds.close()
        print(f'  Saved to: {out_fn}')
    