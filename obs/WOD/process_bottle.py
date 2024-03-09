"""
Code to process WOD (World Ocean Database) bottle data from the NetCDF version, for years (1993-2018)
year 2011-2013 doesn't have bottle data within the model domain
year 2015-2018 doesn't have bottle data within the model domain either
"""
from datetime import datetime, timedelta
import numpy as np
import gsw
import sys
import pandas as pd
from time import time
import glob
import xarray as xr

from lo_tools import Lfun, obs_functions
Ldir = Lfun.Lstart()

# This is a dict of most of the variables after the initial reading.
# We add values to a key for any variable we want to save.
# be sure about the units.
v_dict = {
    'country': '',
    'WOD_cruise_identifier': 'cruise',
    'originators_station_identifier': '',
    'originators_cruise_identifier': '',
    'wod_unique_cast': 'cid',
    'date': '',
    'GMT_time': '',
    'Access_no': '',
    'Platform': '',
    'Cast_Tow_number': '',
    'Orig_Stat_Num': '',
    'dataset': 'name',
    'dbase_orig': '',
    'origflagset': '',
    'z_origflag': '',
    'z_WODflag': '',
    'z_sigfig': '',
    'Temperature': 'IT',
    'Temperature_sigfigs': '',
    'Temperature_WODflag': '',
    'Temperature_origflag': '',
    'Temperature_WODprofileflag': '',
    'Temperature_Scale': '',
    'Temperature_Instrument': '',
    'Salinity': 'SP',
    'Salinity_sigfigs': '',
    'Salinity_WODflag': '',
    'Salinity_origflag': '',
    'Salinity_WODprofileflag': '',
    'Oxygen': 'DO (umol/kg)',
    'Oxygen_sigfigs': '',
    'Oxygen_WODflag': '',
    'Oxygen_origflag': '',
    'Oxygen_WODprofileflag': '',
    'Phosphate': 'PO4 (umol/kg)',
    'Phosphate_sigfigs': '',
    'Phosphate_WODflag': '',
    'Phosphate_origflag': '',
    'Phosphate_WODprofileflag': '',
    'Silicate': 'SiO4 (umol/kg)',
    'Silicate_sigfigs': '',
    'Silicate_WODflag': '',
    'Silicate_origflag': '',
    'Silicate_WODprofileflag': '',
    'Alkalinity': 'TA (mmol/l)',
    'Alkalinity_sigfigs': '',
    'Alkalinity_WODflag': '',
    'Alkalinity_origflag': '',
    'Alkalinity_WODprofileflag': '',
    'Alkalinity_Original_units': '',
    'tCO2': 'DIC (mmol/l)',
    'tCO2_sigfigs': '',
    'tCO2_WODflag': '',
    'tCO2_origflag': '',
    'tCO2_Original_units': '',
    'Nitrate': 'NO3 (umol/kg)',
    'Nitrate_sigfigs': '',
    'Nitrate_WODflag': '',
    'Nitrate_origflag': '',
    'Nitrate_WODprofileflag': '',
    'NO2NO3': 'NO3 (umol/kg)',
    'NO2NO3_sigfigs': '',
    'NO2NO3_WODflag': '',
    'NO2NO3_origflag': '',
    'NO2NO3_Original_units': '',
    'Pressure': 'P (dbar)',
    'Pressure_sigfigs': '',
    'Chlorophyll': 'Chl (ug/l)',
    'Chlorophyll_sigfigs': '',
    'Chlorophyll_WODflag': '',
    'Chlorophyll_WODprofileflag': '',
    'crs': '',
    'WODf': '',
    'WODfp': '',
    'WODfd': '',
    'Oflag': '',
    'Odflag': '',

}

# bottle
source = 'WOD'
otype = 'bottle'
# output location
out_dir = Ldir['LOo'] / 'obs' / source / otype
Lfun.make_dir(out_dir)

testing = False

if testing:
    year_list = [1999]
else:
    year_list = list(range(1993, 2010+1)) + [2014]

# input location
for year in year_list:
    ys = str(year)
    print('\n'+ys)
    # name output files
    out_fn = out_dir / (ys + '.p')
    info_out_fn = out_dir / ('info_' + ys + '.p')
    df = pd.DataFrame()
    # Read files
    files = f'{Ldir["data"]}/obs/{source}/{year}/OSD/wod*.nc'
    fn_list = glob.glob(files)
    for fn in fn_list:
        # print(fn)
        ds = xr.open_dataset(fn)
        # select and rename variables
        df1 = pd.DataFrame()
        for v in list(ds.data_vars):
            # print(f"Variable: {v}")
            if v in v_dict.keys():
                if len(v_dict[v]) > 0:
                    df1[v_dict[v]] = ds[v].values
        df1['z'] = ds['z'].values
        df1['lat'] = ds['lat'].values
        df1['lon'] = ds['lon'].values
        df1['time'] = ds['time'].values
        df1['cid'] = ds['wod_unique_cast'].values
        # Extract the element and decode it
        df1['name'] = ds['dataset'].values.item().decode('utf-8').rstrip('\x00')
        df1['cruise'] = ds['WOD_cruise_identifier'].values.item().decode('utf-8').rstrip('\x00')

        if len(df1['z']) != 1:
            df= pd.concat([df, df1], ignore_index=True)
        else:
            pass
    for col in df.columns:
        if col != 'lon' and col != 'cruise' and col != 'time' and col != 'name':
            df.loc[df[col] <= 0, col] = np.nan
    # drop rows with no good data
    df = df.dropna(axis=0, how='all')
    df = df.reset_index(drop=True)

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
    # - add the results to the DataFrame
    df['SA'] = SA
    df['CT'] = CT
    df['z'] = -1 * df['z']
    rho = gsw.rho(SA,CT,p)

    # (2) units
    for vn in ['DO','NO3', 'PO4', 'SIO4','TA','DIC']:
        if (vn+' (umol/kg)') in df.columns:
            df[vn+' (uM)'] = (rho/1000) * df[vn+' (umol/kg)']
        elif (vn+' (mmol/l)') in df.columns:
            df[vn+' (uM)'] = (1000) * df[vn+' (mmol/l)']

    # (3) retain only selected variables
    cols = ['cid', 'cruise', 'time', 'lat', 'lon', 'name', 'z',
            'CT', 'SA', 'DO (uM)',
            'NO3 (uM)', 'PO4 (uM)', 'SiO4 (uM)',
            'TA (uM)', 'DIC (uM)']
    this_cols = [item for item in cols if item in df.columns]
    df = df[this_cols]

    print(' - processed %d casts' % ( len(df.cid.unique()) ))
    # Renumber cid to be increasing from zero in steps of one.
    df = obs_functions.renumber_cid(df)

    if len(df) > 0:
        # Save the data
        df.to_pickle(out_fn)
        info_df = obs_functions.make_info_df(df)
        info_df.to_pickle(info_out_fn)
