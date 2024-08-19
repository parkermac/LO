"""
Code to process the NCEI Coastal bottle data.

This takes just [] to process 2008-2018.

"""

import pandas as pd
import numpy as np
import gsw
import sys

from lo_tools import Lfun, obs_functions
Ldir = Lfun.Lstart()

# BOTTLE
source = 'nceiCoastal'
otype = 'bottle'
in_dir0 = Ldir['data'] / 'obs' / source
year_list = range(2008,2020)

# output location
out_dir = Ldir['LOo'] / 'obs' / source / otype
Lfun.make_dir(out_dir)

# This is a dict of all the columns after the initial reading.
# We add values to a key for any variable we want to save. I looked
# at units_dict (created below) to be sure about the units.
v_dict = {
    'time':'time',
    'Accession':'',
    'EXPOCODE':'',
    'Cruise_flag':'',
    'Cruise_ID':'cruise',
    'Observation_type':'',
    'Profile_number':'cid',
    'Station_ID':'name',
    'Cast_number':'',
    'Niskin_ID':'',
    'Niskin_flag':'',
    'Sample_ID':'',
    'Year_UTC':'',
    'Month_UTC':'',
    'Day_UTC':'',
    'Time_UTC':'',
    'Latitude':'lat',
    'Longitude':'lon',
    'Depth_bottom':'',
    'Max_sample_depth':'',
    'CTDPRES':'P (dbar)',
    'Depth':'',
    'CTDTEMP_ITS90':'IT',
    'CTDTEMP_flag':'',
    'CTDSAL_PSS78':'',
    'CTDSAL_flag':'',
    'Salinity_PSS78':'',
    'Salinity_flag':'',
    'recommended_Salinity_PSS78':'SP',
    'recommended_Salinity_flag':'',
    'CTDOXY':'',
    'CTDOXY_flag':'',
    'Oxygen':'',
    'Oxygen_flag':'',
    'recommended_Oxygen':'DO (umol/kg)',
    'recommended_Oxygen_flag':'',
    'AOU':'',
    'AOU_flag':'',
    'DIC':'DIC (umol/kg)',
    'DIC_flag':'',
    'TALK':'TA (umol/kg)',
    'TALK_flag':'',
    'pH_TS_measured':'',
    'TEMP_pH':'',
    'pH_flag':'',
    'pH_TS_insitu_measured':'',
    'pH_TS_insitu_calculated':'',
    'fCO2_measured':'',
    'TEMP_fCO2':'',
    'fCO2_flag':'',
    'fCO2_insitu_measured':'',
    'fCO2_insitu_calculated':'',
    'Carbonate_measured':'',
    'TEMP_Carbonate':'',
    'Carbonate_flag':'',
    'Carbonate_insitu_measured':'',
    'Carbonate_insitu_calculated':'',
    'Aragonite':'',
    'Calcite':'',
    'Revelle_Factor':'',
    'Silicate':'SiO4 (umol/kg)',
    'Silicate_flag':'',
    'Phosphate':'PO4 (umol/kg)',
    'Phosphate_flag':'',
    'Nitrate':'NO3 (umol/kg)',
    'Nitrate_flag':'',
    'Nitrite':'NO2 (umol/kg)',
    'Nitrite_flag':'',
    'Nitrate_and_Nitrite':'',
    'Nitrate_and_Nitrite_flag':'',
    'recommended_Nitrate_and_Nitrite':'',
    'recommended_Nitrate_and_Nitrite_flag':'',
    'Ammonium':'NH4 (umol/kg)',
    'Ammonium_flag':'',
}

load_data = True
for year in year_list:
    ys = str(year)
    print('\n'+ys)
    
    # name output files
    out_fn = out_dir / (ys + '.p')
    info_out_fn = out_dir / ('info_' + ys + '.p')
    
    """
    NOTE: There was a garbled Time_UTC column in the csv version, so I used the excel instead.
    This was slow so I converted it to a csv. Also, since the second row was the units, I pulled
    this into a separate csv (and dropped that row from the data csv). Reading in the csv is
    MUCH faster.
    Commands to do this:
    in_fn =  in_dir0 / 'CODAP_NA_v2021.xlsx'
    df00 = pd.read_excel(in_fn, skiprows=[1])
    df00.to_csv(in_dir0 / 'converted_from_excel.csv')
    dfu = pd.read_excel(in_fn, nrows=1)
    dfu.to_csv(in_dir0 / 'units_converted_from_excel.csv')
    """
    
    in_fn = in_dir0 / 'converted_from_excel.csv'

    if load_data:
        # only load this once, because it has all the years through 2017
        df0 = pd.read_csv(in_fn, low_memory=False)

    # a little massaging to get the time column
    df00 = df0[['Year_UTC','Month_UTC','Day_UTC']].copy()
    df000 = df00.rename(columns={'Year_UTC':'year','Month_UTC':'month','Day_UTC':'day'})
    hms = df0['Time_UTC'].copy().to_list()
    hms_hour = [item[:2] for item in hms]
    hms_minute = [item[3:5] for item in hms]
    hms_second = [item[6:] for item in hms]
    df000['hour'] = hms_hour
    df000['minute'] = hms_minute
    df000['second'] = hms_second
    df0['time'] = pd.to_datetime(df000)
    
    units_fn = in_dir0 / 'units_converted_from_excel.csv'
    units_df = pd.read_csv(units_fn)
    keys = list(units_df.columns)
    values = list(units_df.loc[0,])
    units_dict = { k:v for (k,v) in zip(keys, values)}
    
    load_data = False # only load the first time
            
    # select one year
    t = pd.DatetimeIndex(df0.time)
    df1 = df0.loc[t.year==year,:].copy()
    
    # select and rename variables
    df = pd.DataFrame()
    for v in df1.columns:
        if v in v_dict.keys():
            if len(v_dict[v]) > 0:
                df[v_dict[v]] = df1[v]
                
    # missing data is -999
    df[df==-999] = np.nan
    
    # a little more cleaning up
    df = df.dropna(axis=0, how='all') # drop rows with no good data
    df = df[df.time.notna()] # drop rows with bad time
    df = df.reset_index(drop=True)
        
    # Force certain fields to be the same throughout the cast. This dataset already
    # has unique numbers for each cast (Profile_number), which is convenient.
    for cid in df.cid.unique():
        df.loc[df.cid==cid,'lon'] = df[df.cid==cid].lon.values[0]
        df.loc[df.cid==cid,'lat'] = df[df.cid==cid].lat.values[0]
        df.loc[df.cid==cid,'time'] = df[df.cid==cid].time.values[0]
    
    # limit the geographical region to the model domain
    df = df[(df.lon>-130) & (df.lon<-122) & (df.lat>42) & (df.lat<52)]
    
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

    # (2) units
    for vn in ['DO','NO3', 'NO2', 'NH4', 'PO4', 'SIO4','TA','DIC']:
        if (vn+' (umol/kg)') in df.columns:
            df[vn+' (uM)'] = (rho/1000) * df[vn+' (umol/kg)']
        
    # (3) retain only selected variables
    cols = ['cid', 'cruise', 'time', 'lat', 'lon', 'name', 'z',
        'CT', 'SA', 'DO (uM)',
        'NO3 (uM)', 'NO2 (uM)', 'NH4 (uM)', 'PO4 (uM)', 'SiO4 (uM)',
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
