"""
Code for initial processing of ctd and bottle data from nanoos cruises.
"""

from datetime import datetime, timedelta
import numpy as np
import gsw
import sys
import pandas as pd
from time import time

from lo_tools import Lfun
Ldir = Lfun.Lstart()

testing = True

# BOTTLES
source = 'nanoos'
otype = 'bottle'

in_dir0 = Ldir['data'] / 'obs' / source

# process bottles
out_dir = Ldir['LOo'] / 'obs' / source / otype
if testing:
    Lfun.make_dir(out_dir, clean=True)
else:
    Lfun.make_dir(out_dir)

year_list = [2021]

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
        fn_list = [fn_list[0]]
        
    for fn in fn_list:
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
        # selected columns
        cols = ['DATE_UTC', 'TIME_UTC',
        'LATITUDE_DEC', 'LONGITUDE_DEC', 'STATION_NO',
        'CTDPRS_DBAR', 'CTDTMP_DEG_C_ITS90', 'CTDSAL_PSS78',
        'OXYGEN_avg_mg_L', 'SALINITY_PSS78',
        'NITRATE_UMOL_L', 'NITRITE_UMOL_L', 'AMMONIUM_UMOL_L', 'PHOSPHATE_UMOL_L',
        'SILICATE_UMOL_L', 'CTD FLU (mg/m3)','CHLA avg (ug/l)','PHAEOPIGMENT avg (ug/l)',
        'TA_UMOL_KG', 'DIC_UMOL_KG', 'SECCHI DEPTH (m)']
        
        df = pd.read_excel(fn, usecols=cols, parse_dates={'time':['DATE_UTC', 'TIME_UTC']})
        
        cid0 = 0
        if len(df) > 0:
            df = df.rename({'STATION_NO':'name',
                'LONGITUDE_DEC':'lon', 'LATITUDE_DEC':'lat',
                'CTDPRS_DBAR':'P (dbar)', 'CTDTMP_DEG_C_ITS90':'PT', 'CTDSAL_PSS78':'SP',
                'OXYGEN_avg_mg_L':'DO (mg/L)'
                'NITRATE_UMOL_L':'NO3 (uM)', 'NITRITE_UMOL_L':'NO2 (uM)', 'AMMONIUM_UMOL_L':'NH4 (uM)',
                'PHOSPHATE_UMOL_L':'PO4 (uM)', 'SILICATE_UMOL_L':'SiO4 (uM)',
                'CTD FLU (mg/m3)':'Fluor (ug/L)','CHLA avg (ug/l)':'ChlA (ug/L)',
                'PHAEOPIGMENT avg (ug/l)':'Phaeo (ug/L)',
                'TA_UMOL_KG':'TA (umol/kg)', 'DIC_UMOL_KG':'DIC (umol/kg)', 'SECCHI DEPTH (m)':'Secchi (m)'
            }, axis=1)
            # Note that we will save the field "name" for station number, since this dataset has
            # repeat stations which is helpful for plotting sections. Then we will generate our own
            # cid, a unique one for each cast, being careful not to keep them unique for the collection
            # of cruises in this year.
            df['cid'] = np.nan
            cid = cid0
            for name in df.name.unique():
                df.loc[df.name==name,'cid'] = cid
                cid += 1
            for cid in df.cid.unique():
                # Check that there are not two different casts associated with the same Station
                # by looking for large time differences.
                time_diff = df[df.cid==cid].time.values[-1] - df[df.cid==cid].time.values[0]
                time_diff = pd.to_timedelta(time_diff)
                if time_diff.days > 1 or time_diff.days < -1:
                    print('Station %d has time diff of %d days' % (cid, time_diff.days))
                # RESULT: the time_diffs are all zero, so it appears that in this database
                # the Station field is a unique cast identifier.

                # Force certain fields to be the same throughout the cast.
                df.loc[df.cid==cid,'lon'] = df[df.cid==cid].lon.values[0]
                df.loc[df.cid==cid,'lat'] = df[df.cid==cid].lat.values[0]
                df.loc[df.cid==cid,'time'] = df[df.cid==cid].time.values[0]
                
            # derived quantities and unit conversions
        
        

