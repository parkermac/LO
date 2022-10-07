"""
Code to TEST the processing of DFO bottle data.

This is the first use of the new (2022.07.14) cast data organization
scheme.

Runs in about a minute.
"""
from datetime import datetime, timedelta
import numpy as np
import gsw
import sys
import pandas as pd

import loadDFO
from importlib import reload
reload(loadDFO)

from lo_tools import Lfun
Ldir = Lfun.Lstart()

testing = True

in_dir = Ldir['data'] / 'obs' / 'dfo'
# out_dir = Ldir['LOo'] / 'obs' / 'dfo'
# Lfun.make_dir(out_dir)

if testing:
    year_list = [2019]
else:
    year_list = list(range(1930, datetime.now().year + 1))

for year in year_list:
    print(year)
    datelims = (datetime(year,1,1), datetime(year,12,31))
    # out_fn = out_dir / ('bottles_' + str(year) + '.p')
    # info_out_fn = out_dir / ('bottles_info_' + str(year) + '.p')
    
    df = loadDFO.loadDFO_bottle(basedir=str(in_dir), datelims=datelims)
    """
    A typical Dataframe will have columns:
    
    ['Station', 'Lat', 'Lon', 'Chl', 'Chl_units', 'N', 'N_units', 'Si',
           'Si_units', 'SA', 'CT', 'DO', 'DO_units', 'Z', 'dtUTC']
    
    and a mish-mash of units:
    
    In [5]: df.DO_units.unique()
    Out[5]: array(['umol/kg', None], dtype=object) [also mL/L]

    In [6]: df.Chl_units.unique()
    Out[6]: array(['mg/m^3', None], dtype=object)

    In [7]: df.N_units.unique()
    Out[7]: array(['umol/L', None], dtype=object)
    
    In [71]: df.Si_units.unique()
    Out[71]: array(['umol/L', None], dtype=object)
    
    """
    
    # process to LO standard form
    if df is None:
        print(' - no data')
    else:
        df = df.dropna(axis=1, how='all') # drop empty rows
        if len(df) > 0:
            df = df.rename({'Station':'cid', 'Lon':'lon', 'Lat':'lat', 'dtUTC':'time',
                'Z':'z'}, axis=1)
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
            
            # Check for outliers in the units. Use {sets}
            if 'DO_units' in df.columns:
                DO_units = {item for item in df.DO_units.unique() if item != None}
                if not DO_units.issubset({'umol/kg', 'mL/L'}):
                    print(' DO units problem: %s' % (DO_units))
            if 'Chl_units' in df.columns:
                Chl_units = {item for item in df.Chl_units.unique() if item != None}
                if not Chl_units.issubset({'mg/m^3'}):
                    print(' Chl units problem: %s' % (Chl_units))
            if 'N_units' in df.columns:
                N_units = {item for item in df.N_units.unique() if item != None}
                if not N_units.issubset({'umol/L'}):
                    print(' N units problem: %s' % (N_units))
            if 'Si_units' in df.columns:
                Si_units = {item for item in df.Si_units.unique() if item != None}
                if not Si_units.issubset({'umol/L','mmol/m**3'}):
                    print(' Si units problem: %s' % (Si_units))
            # RESULT: there are no outliers beyond out lists, except for 1975.
            
            # fix NO3 units
            df['NO3 (uM)'] = np.nan
            if 'N' in df.columns:
                df['NO3 (uM)'] = df['N']
        
            # fix DO units
            df['DO (uM)'] = np.nan
            if 'DO' in df.columns:
                df['DO (uM)'] = df['DO']
                mask1 = df.DO_units=='umol/kg'
                mask2 = df.DO_units=='mL/L'
                # print(' mask1:%d mask2:%d' % (mask1.sum(),mask2.sum()))
                if mask1.sum() > 0:
                    # calculate in situ density
                    p = gsw.p_from_z(df.z.values, df.lat.values)
                    rho = gsw.rho(df['SA'].values, df['CT'].values, p)
                    df['rho'] = rho
                    df.loc[mask1,'DO (uM)'] *= df.loc[mask1,'rho']/1000
                if mask2.sum() > 0:
                    df.loc[mask2,'DO (uM)'] *= 1.42903
                
            # Fix Si name and handle a few units outliers
            df['Si (uM)'] = np.nan
            if 'Si' in df.columns:
                mask = (df['Si_units']=='umol/L') | (df['Si_units']=='mmol/m**3')
                df.loc[mask,'Si (uM)'] = df.loc[mask,'Si']
            # Note that in 1961-1975 Si_units also showed up as "microg-at/l" but
            # we just won't use those.
            
            # Fix Chl name
            df['Chl (mg m-3)'] = np.nan
            if 'Chl' in df.columns:
                df['Chl (mg m-3)'] = df['Chl']
            
            # clean up columns
            df = df[['cid', 'lon', 'lat', 'time', 'z',
                'SA', 'CT',
                'DO (uM)', 'NO3 (uM)', 'Si (uM)', 'Chl (mg m-3)']]
            df['name'] = None

            # save
            # df.to_pickle(out_fn)
            
            # # Also pull out a dateframe with station info to use for model cast extractions.
            ind = df.cid.unique()
            info_df = pd.DataFrame(index=ind, columns=['lon','lat','time','name'])
            for cid in df.cid.unique():
                info_df.loc[cid,['lon','lat','time']] = df.loc[df.cid==cid,['lon','lat','time']].iloc[0,:]
            info_df.name = None
            info_df.index.name = 'cid'
            info_df['time'] = pd.to_datetime(info_df['time'])
            # info_df.to_pickle(info_out_fn)
    



