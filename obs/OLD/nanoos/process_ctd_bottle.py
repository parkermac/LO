"""
Code for initial processing of ctd and bottle data from nanoos cruises.

Hack notes:

2021
- has three cruises (April: RC0051, July: RC0058, September: RC0063)
- none have DIC or TA, but samples were taken. They just were slow to be
  processed by PMEL. I may get draft values from Anna Boyar.
- the original September bottle file had a column alignment type to it
  appeared to have near-zero nitrate. Anna fixed this and the copy I have
  in LO_data is corrected. Not sure about the website.
- decimal longitudes are positive (fixed with hack)
- April has bad rows at the end of the xlsx with only lat, lon (fixed)
- April bottle file has a typo for the latitude of station 5 (repeats that of 28)
  (fixed with a hack)
"""

from datetime import datetime, timedelta
import numpy as np
import gsw
import sys
import pandas as pd
from time import time

from lo_tools import Lfun
Ldir = Lfun.Lstart()

testing = False

# CTD
# Not written yet

# BOTTLE
source = 'nanoos'
otype = 'bottle'

in_dir0 = Ldir['data'] / 'obs' / source

# process bottles
out_dir = Ldir['LOo'] / 'obs' / source / otype
if testing:
    Lfun.make_dir(out_dir, clean=True)
else:
    Lfun.make_dir(out_dir)

year_list = [2017]

for year in year_list:
    ys = str(year)
    
    # name output files
    out_fn = out_dir / (str(year) + '.p')
    info_out_fn = out_dir / ('info_' + str(year) + '.p')
    
    in_dirs = list(in_dir0.glob('*'+ys+'*'))
    
    fn_list = []
    for in_dir in in_dirs:
        fn_list.append(list(in_dir.glob('*labupcast*'))[0])
        
    if testing:
        fn_list = [fn_list[1]]
    
    cid0 = 0
    DF = pd.DataFrame()
    for fn in fn_list:
        print(fn.name)
        # selected columns to ingest
        if year in [2021]:
            """
            Available columns:
            ['record no', 'CRUISE_ID', 'DATE_UTC', 'TIME_UTC', 'DATE_LOCAL',
            'TIME_LOCAL', 'LATITUDE_DEG', 'LONGITUDE_DEG', 'LATITUDE_DEC',
            'LONGITUDE_DEC', 'STATION_NO', 'NISKIN_NO', 'NISKIN_NO_FLAG_W',
            'CTDPRS_DBAR', 'DEPTH (M)', 'CTDTMP_DEG_C_ITS90', 'CTDTMP2_DEC_C_ITS90',
            'CTDTMP_FLAG_W', 'CTD/TEMP_COMMENTS', 'CTDSAL_PSS78', 'CTDSAL2_PSS78',
            'CTDSAL_FLAG', 'CTD/SAL_COMMENTS', 'SIGMATHETA_KG_M3',
            'SIGMATHETA2_KG_M3', 'CTDOXY_MG_L_1', 'CTDOXY_MG_L_2',
            'CTDOXY_MG_L_AVG', 'CTDOXY_FLAG', 'CTD/O2_COMMENTS', 'OXYGEN_MG_L_1',
            'OXYGEN_MG_L_2', 'OXYGEN_MG_L_3', 'OXYGEN_avg_mg_L', 'OXYGEN_UMOL_KG',
            'OXYGEN_FLAG_W', 'OXYGEN COMMENTS', 'CTD_PH', 'SALINITY_PSS78',
            'SALINITY_FLAG', 'Nutrient lab temperature', 'NITRATE_UMOL_L',
            'NITRITE_UMOL_L', 'AMMONIUM_UMOL_L', 'PHOSPHATE_UMOL_L',
            'SILICATE_UMOL_L', 'NUTRIENTS_FLAG_W', 'CTD FLU (mg/m3)', 'CHLA (ug/l)',
            'CHLA 2 (ug/l)', 'CHLA avg (ug/l)', 'CHLA_FLAG', 'PHAEOPIGMENT (ug/l)',
            'PHAEOPIGMENT 2 (ug/l)', 'PHAEOPIGMENT avg (ug/l)', 'PHAEOPIGMENT_FLAG',
            'TA_UMOL_KG', 'TA_FLAG_W', 'DIC_UMOL_KG', 'DIC_FLAG_W',
            'SECCHI DEPTH (m)']
            """
            cols = ['CRUISE_ID', 'DATE_UTC', 'TIME_UTC',
            'LATITUDE_DEC', 'LONGITUDE_DEC', 'STATION_NO',
            'CTDPRS_DBAR', 'CTDTMP_DEG_C_ITS90', 'CTDSAL_PSS78',
            'OXYGEN_avg_mg_L',
            'NITRATE_UMOL_L', 'NITRITE_UMOL_L', 'AMMONIUM_UMOL_L', 'PHOSPHATE_UMOL_L',
            'SILICATE_UMOL_L', 'CTD FLU (mg/m3)',
            'CHLA avg (ug/l)','PHAEOPIGMENT avg (ug/l)',
            'TA_UMOL_KG', 'DIC_UMOL_KG', 'SECCHI DEPTH (m)']
            # load the data to a DataFrame
            df = pd.read_excel(fn, usecols=cols, dtype={'DATE_UTC':str, 'TIME_UTC':str})
            # NOTE: this would have been a lot easier using
            # parse_dates={'time':['DATE_UTC', 'TIME_UTC']} in the pd.read_excel() call,
            # but there are some bad rows with missing data, so we slog through doing it
            # from scratch.
            dstr = df.DATE_UTC
            tstr = df.TIME_UTC
            tlist = []
            for ii in range(len(dstr)):
                try:
                    tlist.append(pd.to_datetime(dstr[ii][:10] + ' ' + tstr[ii][:8]))
                except:
                    tlist.append(np.nan)
            df['time'] = tlist
            df = df.dropna(axis=0, how='all') # drop rows with no good data
            df = df[df.time.notna()] # drop rows with bad time
            
        elif year in [2017]:
            # this only works for the first two cruises!
            cols = ['Station number',
               'Date_UTC', 'Time_UTC',
               ' CruiseID',
               ' Potemp090C', ' Sal00', ' PrDM',
               ' Latitude', ' Longitude',
               'CTDOXY_UMOL_KG_ADJ',
               'NITRATE_UMOL_L', 'NITRITE_UMOL_L',
               'AMMONIUM_UMOL_L', 'PHOSPHATE_UMOL_L', 'SILICATE_UMOL_L',
               'CHLA avg (ug/l)',
               'PHAEOPIGMENT avg (ug/l)']
            df = pd.read_excel(fn, usecols=cols, parse_dates={'time':['Date_UTC', 'Time_UTC']})
            df = df.dropna(axis=0, how='all') # drop rows with no good data
        
        # if we got any good rows, proceed with processing
        if len(df) > 0:
            # rename variables for convenience (some will be final names)
            if year in [2021]:
                df = df.rename({'STATION_NO':'name', 'CRUISE_ID':'cruise',
                    'LONGITUDE_DEC':'lon', 'LATITUDE_DEC':'lat',
                    'CTDPRS_DBAR':'P (dbar)', 'CTDTMP_DEG_C_ITS90':'PT', 'CTDSAL_PSS78':'SP',
                    'OXYGEN_avg_mg_L':'DO (mg/L)',
                    'NITRATE_UMOL_L':'NO3 (uM)', 'NITRITE_UMOL_L':'NO2 (uM)',
                    'AMMONIUM_UMOL_L':'NH4 (uM)',
                    'PHOSPHATE_UMOL_L':'PO4 (uM)', 'SILICATE_UMOL_L':'SiO4 (uM)',
                    'CTD FLU (mg/m3)':'Fluor (ug/L)','CHLA avg (ug/l)':'ChlA (ug/L)',
                    'PHAEOPIGMENT avg (ug/l)':'Phaeo (ug/L)',
                    'TA_UMOL_KG':'TA (umol/kg)', 'DIC_UMOL_KG':'DIC (umol/kg)',
                    'SECCHI DEPTH (m)':'Secchi (m)',
                }, axis=1)
            elif year in [2017]:
                df = df.rename({'Station number':'name',
                   ' CruiseID':'cruise',
                   ' Potemp090C':'PT', ' Sal00':'SP', ' PrDM':'P (dbar)',
                   ' Latitude':'lat', ' Longitude':'lon',
                   'CTDOXY_UMOL_KG_ADJ':'DO (umol/kg)',
                   'NITRATE_UMOL_L':'NO3 (uM)', 'NITRITE_UMOL_L':'NO2 (uM)',
                   'AMMONIUM_UMOL_L':'NH4 (uM)', 'PHOSPHATE_UMOL_L':'PO4 (uM)', 'SILICATE_UMOL_L':'SiO4 (uM)',
                   'CHLA avg (ug/l)':'ChlA (ug/L)',
                   'PHAEOPIGMENT avg (ug/l)':'Phaeo (ug/L)'
                }, axis=1)
            
            # add the "cid" (cast ID) column
            #
            # Note that we will save the field "name" for station number, since this dataset has
            # repeat stations which is helpful for plotting sections. Then we will generate our own
            # cid, a unique one for each cast, being careful to keep them unique for the collection
            # of cruises in this year, even though a station may be repeated on all cruises.
            #
            # We will also save the field "cruise" as a convenient way to select a collection of
            # casts from a given month.
            df['cid'] = np.nan
            cid = cid0
            for name in df.name.unique():
                df.loc[df.name==name,'cid'] = cid
                cid += 1
            for cid in df.cid.unique():
                # Check that there are not two different casts associated with the same station
                # by looking for large time differences. Pretty ad hoc.
                time_diff = df[df.cid==cid].time.values[-1] - df[df.cid==cid].time.values[0]
                time_diff = pd.to_timedelta(time_diff)
                if time_diff.days > 1 or time_diff.days < -1:
                    print('Station %d has time diff of %d days' % (cid, time_diff.days))
                # RESULT: the time_diffs are all zero, so it appears that in this database
                # the Station field is a unique cast identifier.

                # Force certain fields to be the same throughout the cast.
                df.loc[df.cid==cid,'lon'] = - df[df.cid==cid].lon.values[0] # also fix sign
                df.loc[df.cid==cid,'lat'] = df[df.cid==cid].lat.values[0]
                df.loc[df.cid==cid,'time'] = df[df.cid==cid].time.values[0]
                
            # Next make derived quantities and do unit conversions
            
            # (1) Create CT, SA, and z
            # - pull out variables
            SP = df.SP.to_numpy()
            pt = df.PT.to_numpy()
            p = df['P (dbar)'].to_numpy()
            lon = -(np.abs(df.lon.to_numpy())) # force sign to be negative
            df['lon'] = lon
            lat = df.lat.to_numpy()
            # - do the conversions
            SA = gsw.SA_from_SP(SP, p, lon, lat)
            CT = gsw.CT_from_pt(SA, pt)
            z = gsw.z_from_p(p, lat)
            # - add the results to the DataFrame
            df['SA'] = SA
            df['CT'] = CT
            df['z'] = z
            
            # (2) units
            if 'DO (gg/L)' in df.columns:
                df['DO (uM)'] = (1000/32) * df['DO (mg/L)']
            elif 'DO (umol/kg)' in df.columns:
                rho = gsw.rho(SA,CT,p)
                df['DO (uM)'] = (rho/1000) * df['DO (umol/kg)']
            
            # (3) retain only selected variables
            cols = ['cid', 'cruise', 'time', 'lat', 'lon', 'name', 'z',
                'CT', 'SA', 'DO (uM)',
                'NO3 (uM)', 'NO2 (uM)', 'NH4 (uM)', 'PO4 (uM)', 'SiO4 (uM)',
                'Fluor (ug/L)', 'ChlA (ug/L)', 'Phaeo (ug/L)', 'TA (umol/kg)',
                'DIC (umol/kg)', 'Secchi (m)']
            this_cols = [item for item in cols if item in df.columns]
            df = df[this_cols]
            DF = pd.concat([DF,df])
            
            print(' - processed %d casts' % ( len(df.cid.unique()) ))
            cid0 = df.cid.max() + 1
            
    # Sort the result by time, and sort each cast to be bottom to top
    DF = DF.sort_values(['time','z'], ignore_index=True)
        
    # Rework cid to also be increasing in time
    a = DF[['time','cid']].copy()
    a['cid_alt'] = np.nan
    ii = 0
    for t in a.time.unique():
        a.loc[a.time==t,'cid_alt'] = ii
        ii += 1
    DF['cid'] = a['cid_alt'].copy()

    # clean up spaces in cruise names
    a = DF['cruise'].to_list()
    aa = [item.strip() for item in a]
    DF['cruise'] = aa
        
    # Hack: fix a location typo
    DF.loc[(DF.cruise=='RC0051') & (DF.name==5), 'lat'] = 47.883

    # Save the data
    DF.to_pickle(out_fn)
    
    # Also pull out a dateframe with station info to use for model cast extractions.
    ind = DF.cid.unique()
    col_list = ['lon','lat','time','name','cruise']
    info_df = pd.DataFrame(index=ind, columns=col_list)
    for cid in DF.cid.unique():
        info_df.loc[cid,col_list] = DF.loc[DF.cid==cid,col_list].iloc[0,:]
    info_df.index.name = 'cid'
    info_df['time'] = pd.to_datetime(info_df['time'])
    info_df.to_pickle(info_out_fn)

        

