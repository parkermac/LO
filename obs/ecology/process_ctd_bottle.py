"""
Code for initial processing of WA Dept. of Ecology ctd and bottle files.

These have somewhat different formats over the years.

Notes 2020-2021:
- DO is only a field in the ctd file, but it is always missing.
"""

from datetime import datetime, timedelta
import numpy as np
import gsw
import sys
import pandas as pd
from time import time

from lo_tools import Lfun
Ldir = Lfun.Lstart()

source = 'ecology'
in_dir0 = Ldir['data'] / 'obs' / source


testing = False

do_ctd = True
do_bottle = False

# CTD
if do_ctd:
    otype = 'ctd'
    
    # process ctd casts
    out_dir = Ldir['LOo'] / 'obs' / source / otype
    if testing:
        Lfun.make_dir(out_dir, clean=True)
    else:
        Lfun.make_dir(out_dir)
        
    if testing:
        year_list = [2021]
    else:
        year_list = [2020,2021]
    
    for year in year_list:
        tt0 = time()
        ys = str(year)
    
        if year in [2020, 2021]:
            in_dir = in_dir0 / 'ecy_2020_2021'
    
            # name output files
            out_fn = out_dir / (str(year) + '.p')
            info_out_fn = out_dir / ('info_' + ys + '.p')
    
            # info about stations
            sta_fn = in_dir / 'Sites.xlsx'
            sta_df = pd.read_excel(sta_fn)
            # columns are:
            # ['SiteID', 'SiteCode', 'SiteName', 'SiteDescription', 'SiteTypeCode',
            # 'LatitudeDecimal', 'LongitudeDecimal']
            sta_df = sta_df.set_index('SiteCode')
        
            # data
            # tt00 = time()
            ctd_fn = in_dir / '2020_to_2021_CTD_data.csv'
            ctd_df = pd.read_csv(ctd_fn, parse_dates={'time':['Date']})
            ctd_df = ctd_df[ctd_df.Year==year].copy()
            # columns are:
            # ['Date', 'Year', 'Month', 'Station', 'Rep', 'Depth', 'Temperature',
            #        'Salinity', 'Density', 'DOAdj', 'FluorAdj', 'XTrans', 'Turbidity', 'pH',
            #        'PAR', 'SUNA']
            
            ctd_df['z'] = - ctd_df.Depth
            # NOTE: We are making some guesses about units here!
            ctd_df = ctd_df.rename({'Station':'name', 'Temperature':'PT', 'Salinity':'SP',
                'DOAdj':'DO (uM)'}, axis=1)
            
            # NOTE: The times in this file are just dates. Here we add a bit like this,
            # assuming noon local time is a better estimate of the sampling time.
            ctd_df.time += timedelta(hours=20) # hour 20 UTC is noon PST
            # print('time for initial read = %0.2f sec' % ((time()-tt00)))
            # sys.stdout.flush()
    
            DF = pd.DataFrame()
    
            # identify a cast by (i) a day, and (ii) a station name
            cid = 0
            for t in ctd_df.time.unique():
                for name in sta_df.index:
                    df = ctd_df[(ctd_df.name==name) & (ctd_df.time==t)].copy()
            
                    if len(df) > 0:
                        # tt00 = time()
                        df['cid'] = cid
                        df['lon'] = sta_df.loc[name, 'LongitudeDecimal']
                        df['lat'] = sta_df.loc[name, 'LatitudeDecimal']
                
                        # retain only selected variables (copied from obs/nanoos)
                        cols = ['cid', 'cruise', 'time', 'lat', 'lon', 'name', 'z',
                            'PT', 'SP', 'DO (uM)']
                        this_cols = [item for item in cols if item in df.columns]
                        df = df[this_cols]
                        # print(' - time for one cast df = %0.2f sec' % ((time()-tt00)))
                        # sys.stdout.flush()
                        
                        # tt00 = time()
                        DF = pd.concat([DF,df])
                        # print(' - time for concat = %0.2f sec' % ((time()-tt00)))
                        # sys.stdout.flush()
                                        
                        cid += 1
                        
                if testing and (cid >= 3):
                    break
            print('CTD: %s - processed %d casts in %d sec' % (ys, len(DF.cid.unique()),int(time()-tt0) ))
    
        # Sort the result by time, and sort each cast to be bottom to top
        DF = DF.sort_values(['time','cid','z'], ignore_index=True)
        
        # calculate CT and SA
        z = DF.z.to_numpy()
        PT = DF.PT.to_numpy()
        SP = DF.SP.to_numpy()
        lon = DF.lon.to_numpy()
        lat = DF.lat.to_numpy()
        P = gsw.p_from_z(z, lat)
        SA = gsw.SA_from_SP(SP, P, lon, lat)
        CT = gsw.CT_from_pt(SA, PT)
        # - add the results to the DataFrame
        DF['SA'] = SA
        DF['CT'] = CT
        DF.drop(['PT','SP'], axis=1)
        
        # Save the data
        DF.to_pickle(out_fn)
    
        # Also pull out a dateframe with station info to use for model cast extractions.
        ind = DF.cid.unique()
        col_list = ['lon','lat','time','name']
        info_df = pd.DataFrame(index=ind, columns=col_list)
        for cid in DF.cid.unique():
            info_df.loc[cid,col_list] = \
                DF.loc[DF.cid==cid,col_list].iloc[0,:]
        info_df.index.name = 'cid'
        info_df['time'] = pd.to_datetime(info_df['time'])
        info_df.to_pickle(info_out_fn)

# BOTTLE
if do_bottle:
    otype = 'bottle'
    # process bottles
    out_dir = Ldir['LOo'] / 'obs' / source / otype
    if testing:
        Lfun.make_dir(out_dir, clean=True)
    else:
        Lfun.make_dir(out_dir)

    year_list = [2020, 2021]

    for year in year_list:
        ys = str(year)
    
        if year in [2020, 2021]:
            in_dir = in_dir0 / 'ecy_2020_2021'
    
            # name output files
            out_fn = out_dir / (str(year) + '.p')
            info_out_fn = out_dir / ('info_' + ys + '.p')
    
            # info about stations
            sta_fn = in_dir / 'Sites.xlsx'
            sta_df = pd.read_excel(sta_fn)
            # columns are:
            # ['SiteID', 'SiteCode', 'SiteName', 'SiteDescription', 'SiteTypeCode',
            # 'LatitudeDecimal', 'LongitudeDecimal']
            sta_df = sta_df.set_index('SiteCode')
        
            # data
            nuts_fn = in_dir / '2020_to_2021_Nutrients_and_Chl_data.csv'
            nuts_df = pd.read_csv(nuts_fn, parse_dates={'time':['Date']})
            nuts_df = nuts_df[nuts_df.Year==year].copy()
            # columns are:
            # ['time', 'Year', 'Month', 'Station', 'Rep', 'Depth', 'NO3', 'NO2', 'NH4',
            # 'PO4', 'SiOH4', 'Chla']
            nuts_df['z'] = - nuts_df.Depth
            # NOTE: We are making some guesses about units here!
            nuts_df = nuts_df.rename({'Station':'name', 'NO3':'NO3 (uM)', 'NO2':'NO2 (uM)', 'NH4':'NH4 (uM)',
                    'PO4':'PO4 (uM)', 'SiOH4':'SiO4 (uM)', 'Chla':'ChlA (ug/L)'}, axis=1)
                
            # NOTE: The times in this file are just dates. Here we add a bit like this,
            # assuming noon local time is a better estimate of the sampling time.
            nuts_df.time += timedelta(hours=20) # hour 20 UTC is noon PST
        
            DF = pd.DataFrame()
        
            # identify a cast by (i) a day, and (ii) a station name
            cid = 0
            for t in nuts_df.time.unique():
                for name in sta_df.index:
                    df = nuts_df[(nuts_df.name==name) & (nuts_df.time==t)].copy()
                
                    # NOTE: until I know what dept = 999 means I will drop this column.
                    df = df[df.z != -999]
                    if len(df) > 0:
                        df['cid'] = cid
                        df['lon'] = sta_df.loc[name, 'LongitudeDecimal']
                        df['lat'] = sta_df.loc[name, 'LatitudeDecimal']
                    
                        # retain only selected variables (copied from obs/nanoos)
                        cols = ['cid', 'cruise', 'time', 'lat', 'lon', 'name', 'z',
                            'CT', 'SA', 'DO (uM)',
                            'NO3 (uM)', 'NO2 (uM)', 'NH4 (uM)', 'PO4 (uM)', 'SiO4 (uM)',
                            'Fluor (ug/L)', 'ChlA (ug/L)', 'Phaeo (ug/L)', 'TA (umol/kg)',
                            'DIC (umol/kg)', 'Secchi (m)']
                        this_cols = [item for item in cols if item in df.columns]
                        df = df[this_cols]
                        DF = pd.concat([DF,df])
                    
                        cid += 1
            print('Bottles: %s - processed %d casts' % (ys, len(DF.cid.unique()) ))
        
        # Sort the result by time, and sort each cast to be bottom to top
        DF = DF.sort_values(['time','cid','z'], ignore_index=True)
        
        # Save the data
        DF.to_pickle(out_fn)
    
        # Also pull out a dateframe with station info to use for model cast extractions.
        ind = DF.cid.unique()
        col_list = ['lon','lat','time','name']
        info_df = pd.DataFrame(index=ind, columns=col_list)
        for cid in DF.cid.unique():
            info_df.loc[cid,col_list] = \
                DF.loc[DF.cid==cid,col_list].iloc[0,:]
        info_df.index.name = 'cid'
        info_df['time'] = pd.to_datetime(info_df['time'])
        info_df.to_pickle(info_out_fn)