"""
Code to process the NCEI Coastal bottle data.

This takes just [] to process 2008-2018.

"""

import pandas as pd
import numpy as np
import gsw
import sys

from lo_tools import Lfun
Ldir = Lfun.Lstart()

# BOTTLE
source = 'nceiCoastal'
otype = 'bottle'
in_dir0 = Ldir['data'] / 'obs' / source
year_list = [2017] #range(2008,2019)

# output location
out_dir = Ldir['LOo'] / 'obs' / source / otype
Lfun.make_dir(out_dir)

# This is a dict of all the columns after the initial reading.
# We add values to a key for any variable we want to save
# v_dict = {
#     'time':'time',
#     'record':'',
#     'EXPOCODE':'',
#     'CRUISE_ID':'cruise',
#     'DATE_LOCAL':'',
#     'TIME_LOCAL':'',
#     'LONGITUDE_DEC':'lon',
#     'LATITUDE_DEC':'lat',
#     'STATION_NO':'name',
#     'NISKIN_NO':'',
#     'CTDPRS_DBAR':'P (dbar)',
#     'CTDTMP_DEG_C_ITS90':'IT', # in situ temperature, deg C
#     'CTDTMP_FLAG_W':'',
#     'CTDSAL_PSS78':'SP',
#     'CTDSAL_FLAG_W':'',
#     'SIGMATHETA_KG_M3':'',
#     'CTDOXY_UMOL_KG_ADJ':'DO (umol/kg)',
#     'CTDOXY_UMOL_KG':'',
#     'CTDOXY_MG_L_1':'',
#     'CTDOXY_MG_L_2':'',
#     'CTDOXY_FLAG_W':'',
#     'OXYGEN_UMOL_KG':'',
#     'OXYGEN_MG_L_1':'',
#     'OXYGEN_MG_L_2':'',
#     'OXYGEN_MG_L_3':'',
#     'OXYGEN_FLAG_W':'',
#     'TA_UMOL_KG':'TA (umol/kg)',
#     'DIC_UMOL_KG':'DIC (umol/kg)',
#     'TA_FLAG_W':'',
#     'DIC_FLAG_W':'',
#     'NITRATE_UMOL_KG':'',
#     'NITRATE_UMOL_L':'NO3 (uM)',
#     'NITRITE_UMOL_KG':'',
#     'NITRITE_UMOL_L':'NO2 (uM)',
#     'AMMONIA_UMOL_KG':'',
#     'AMMONIUM_UMOL_L':'NH4 (uM)',
#     'PHOSPHATE_UMOL_KG':'',
#     'PHOSPHATE_UMOL_L':'PO4 (uM)',
#     'SILICATE_UMOL_KG':'',
#     'SILICATE_UMOL_L':'Si4 (uM)',
#     'NUTRIENTS_FLAG_W':'',
# }

v_dict = {
    'Accession':'',
    'EXPOCODE':'',
    'Cruise_flag':'',
    'Cruise_ID':'',
    'Observation_type':'',
    'Profile_number':'',
    'Station_ID':'',
    'Cast_number':'',
    'Niskin_ID':'',
    'Niskin_flag':'',
    'Sample_ID':'',
    'Year_UTC':'',
    'Month_UTC':'',
    'Day_UTC':'',
    'Time_UTC':'',
    'Latitude':'',
    'Longitude':'',
    'Depth_bottom':'',
    'Max_sample_depth':'',
    'CTDPRES':'',
    'Depth':'',
    'CTDTEMP_ITS90':'',
    'CTDTEMP_flag':'',
    'CTDSAL_PSS78':'',
    'CTDSAL_flag':'',
    'Salinity_PSS78':'',
    'Salinity_flag':'',
    'recommended_Salinity_PSS78':'',
    'recommended_Salinity_flag':'',
    'CTDOXY':'',
    'CTDOXY_flag':'',
    'Oxygen':'',
    'Oxygen_flag':'',
    'recommended_Oxygen':'',
    'recommended_Oxygen_flag':'',
    'AOU':'',
    'AOU_flag':'',
    'DIC':'',
    'DIC_flag':'',
    'TALK':'',
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
    'Silicate':'',
    'Silicate_flag':'',
    'Phosphate':'',
    'Phosphate_flag':'',
    'Nitrate':'',
    'Nitrate_flag':'',
    'Nitrite':'',
    'Nitrite_flag':'',
    'Nitrate_and_Nitrite':'',
    'Nitrate_and_Nitrite_flag':'',
    'recommended_Nitrate_and_Nitrite':'',
    'recommended_Nitrate_and_Nitrite_flag':'',
    'Ammonium':'',
    'Ammonium_flag':'',
}

load_data = True
for year in year_list:
    ys = str(year)
    print('\n'+ys)
    
    # name output files
    out_fn = out_dir / (ys + '.p')
    info_out_fn = out_dir / ('info_' + ys + '.p')
    
    if year in range(2008,2019):
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
        # in_fn =  in_dir0 / 'CODAP_NA_v2021.xlsx'
        # units_df = pd.read_excel(in_fn, nrows=1
        # keys = list(units_df.columns)
        # values = list(units_df.loc[0,])
        # units_dict = { k:v for (k,v) in zip(keys, values)}
        # df0 = pd.read_excel(in_fn, skiprows=[1], nrows=10, parse_dates={'time':['Year_UTC','Month_UTC','Day_UTC','Time_UTC']})
        
        in_fn = in_dir0 / 'converted_from_excel.csv'
        df0 = pd.read_csv(in_fn, low_memory=False, parse_dates={'time':['Year_UTC','Month_UTC','Day_UTC','Time_UTC']})
        
        units_fn = in_dir0 / 'units_converted_from_excel.csv'
        units_df = pd.read_csv(units_fn)
        keys = list(units_df.columns)
        values = list(units_df.loc[0,])
        units_dict = { k:v for (k,v) in zip(keys, values)}
        
        load_data = False # only load the first time
        
    sys.exit()
    
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

    # (2) units
    if 'DO (mg/L)' in df.columns:
        df['DO (uM)'] = (1000/32) * df['DO (mg/L)']
    for vn in ['DO','TA','DIC']:
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
