"""
Code to process WA Dept. of Ecology historical monitoring BOTTLE data for Puget Sound.

To process data received via records request by Dakota Mascarenas on 2024/05/31.

Initial author date: 2024/07/10

Finalized for group use: 2025/09/03

Last updated: 2025/11/25 to correct conversion unit errors for N, P, and Si

Written by: Dakota Mascarenas

NOTE: Despite the labeling on the Excel files received from WA Dept. of Ecology records request, we consider all data in this set bottle (discrete) data. Email requests for information were sent on 2025/05/06 with follow up on 2025/09/03.

"""

import pandas as pd
import numpy as np
import gsw

from lo_tools import Lfun, obs_functions
Ldir = Lfun.Lstart()


# source location
source = 'ecology_his'
otype = 'bottle'
in_dir0 = Ldir['data'] / 'obs' / source
year_list = range(1973,1999)

# output location
out_dir = Ldir['LOo'] / 'obs' / source / otype
Lfun.make_dir(out_dir)

# Load big data sets and stations.
big_df_raw0 = pd.read_excel(in_dir0/ 'Aug1973toOct1989CTDandDiscrete.xlsx', parse_dates=['Date'], skiprows=[0])
big_df_raw1 = pd.read_excel(in_dir0/ 'Nov1989toDec1998Discrete.xlsx', parse_dates=['Date'], skiprows=[0])
sta_df = pd.read_excel(in_dir0 / 'ParkerMacCreadyCoreStationInfoFeb2018.xlsx') #station list copied from apogee: dat1/parker/LO_data/ecology/

# Parse station lat/lon.
xx = sta_df['Long_NAD83 (deg / dec_min)'].values
yy = sta_df['Lat_NAD83 (deg / dec_min)'].values
lon = [-(float(x.split()[0]) + float(x.split()[1])/60) for x in xx]
lat = [(float(y.split()[0]) + float(y.split()[1])/60) for y in yy]
sta_df['lon'] = lon
sta_df['lat'] = lat

# Make column naming consistent for variables across datasets.
big_df_raw0 = big_df_raw0.rename(columns = {'NO2+NO3':'NO2+NO3 (total size fraction)',
                                            'NO2':'NO2 (total size fraction)',
                                            'NH4':'NH4 (total size fraction)',
                                            'TP':'TP (total size fraction)',
                                            'PO4':'PO4 (total size fraction)',
                                            'TOC':'TOC (total size fraction)'})
big_df_raw1 = big_df_raw1.rename(columns = {'NO2+NO3':'NO2+NO3 (dissolved size fraction)',
                                            'NO2':'NO2 (dissolved size fraction)',
                                            'NH4':'NH4 (dissolved size fraction)',
                                            'PO4':'PO4 (dissolved size fraction)',
                                            'NO2+NO3.1':'NO2+NO3 (total size fraction)',
                                            'NO2.1':'NO2 (total size fraction)', 
                                            'NH4.1':'NH4 (total size fraction)', 
                                            'PO4.1':'PO4 (total size fraction)', 
                                            'TP':'TP (total size fraction)'})

# Concatenate two datasets.
big_df_raw = pd.concat([big_df_raw0, big_df_raw1])

# Merge station data.
big_df = big_df_raw.merge(sta_df[['Station','lat', 'lon']], on = 'Station', how='left')

# Create unique cast IDs (cid).
big_df_use = big_df.copy()
big_df_use['cid'] = np.nan
big_df_use['unique_date_location'] = big_df_use['Station'] + big_df_use['Date'].dt.strftime('%Y-%m-%d')
c = 0
for pid in big_df_use['unique_date_location'].unique(): # profile ID is unique identifier
    big_df_use.loc[big_df_use['unique_date_location'] == pid, 'cid'] = c
    c+=1

# Create dictionary for important variable and column names.
data_original_names = ['Temp', 'Salinity', 'DO', 'Chla', 'NO2 (total size fraction)', 'NO2+NO3 (total size fraction)', 'NH4 (total size fraction)', 'PO4 (total size fraction)']
data_new_names = ['IT', 'SP', 'DO (mg -L)', 'Chl (ug -L)', 'NO2 (mg -L)', 'NO2+NO3 (mg -L)', 'NH4 (mg -L)', 'PO4 (mg -L)']
v_dict = dict(zip(data_original_names, data_new_names))
v_dict['cid'] = 'cid'
v_dict['Date'] = 'time'
v_dict['lat'] = 'lat'
v_dict['lon'] = 'lon'
v_dict['DepthInterval'] = 'z' # will be converted to negative later in script
v_dict['Station'] = 'name'

# Loop through to rename variables and columns, clean the dataset, and produce output dataframes.
df0 = big_df_use.copy()
for year in year_list:
    ys = str(year)
    print('\n'+ys)
    out_fn = out_dir / (ys + '.p')
    info_out_fn = out_dir / ('info_' + ys + '.p')
    t = pd.DatetimeIndex(df0['Date'])
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
    # add the results to the dataframe
    df['SA'] = SA
    df['CT'] = CT
    rho = gsw.rho(SA,CT,p)
    # unit coversions
    if 'DO (mg -L)' in df.columns:
        df['DO (uM)'] = (1000/32) * df['DO (mg -L)']
    if 'NH4 (mg -L)' in df.columns:
        df['NH4 (uM)'] = (1000/14) * df['NH4 (mg -L)']
    if 'NO3 (mg -L)' in df.columns:
        df['NO3 (uM)'] = (1000/14) * df['NO3 (mg -L)']
    if 'SiO4 (mg -L)' in df.columns:
        df['SiO4 (uM)'] = (1000/28.0855) * df['SiO4 (mg -L)']
    if 'PO4 (mg -L)' in df.columns:
        df['PO4 (uM)'] = (1000/30.973762) * df['PO4 (mg -L)']
    if 'Chl (ug -L)' in df.columns:
        df['Chl (mg m-3)'] = df['Chl (ug -L)']
    if 'NO2+NO3 (mg -L)' in df.columns:
        df['NO3 (uM)'] = (1000/62) * (df['NO2+NO3 (mg -L)'] - df['NO2 (mg -L)'])
    for vn in ['TA','DIC']:
        if (vn+' (umol -kg)') in df.columns:
            df[vn+' (uM)'] = (rho/1000) * df[vn+' (umol -kg)']
    # retain only selected variables
    df['cruise'] = ''
    cols = ['cid', 'time', 'lat', 'lon', 'z', 'cruise', 'name',
        'CT', 'SA', 'DO (uM)',
        'NO3 (uM)', 'NH4 (uM)', 'PO4 (uM)', 'SiO4 (uM)',
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