"""
Code to process DFO bottle data.

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

testing = False

in_dir = Ldir['data'] / 'obs' / 'dfo'
out_dir = Ldir['LOo'] / 'obs' / 'dfo'
Lfun.make_dir(out_dir)

if testing:
    year_list = [1930, 1975]
else:
    #year_list = list(range(1930, datetime.now().year + 1))
    year_list = list(range(1930, datetime.now().year + 1))

for year in year_list:
    print(year)
    datelims = (datetime(year,1,1), datetime(year,12,31))
    out_fn = out_dir / ('bottles_' + str(year) + '.p')
    
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
            df = df.rename({'Station':'sta', 'Lon':'lon', 'Lat':'lat', 'dtUTC':'time',
                'SA':'salt (SA g kg-1)', 'CT':'temp (CT degC)', 'Z':'z'}, axis=1)
            for sta in df.sta.unique():
                
                # Check that there are not two different casts associated with the same Station
                # by looking for large time differences.
                time_diff = df[df.sta==sta].time.values[-1] - df[df.sta==sta].time.values[0]
                time_diff = pd.to_timedelta(time_diff)
                if time_diff.days > 1 or time_diff.days < -1:
                    print('Station %d has time diff of %d days' % (sta, time_diff.days))
                # RESULT: the time_diffs are all zero, so it appears that in this database
                # the Station field is a unique cast identifier.

                # Force certain fields to be the same throughout the cast.
                df.loc[df.sta==sta,'lon'] = df[df.sta==sta].lon.values[0]
                df.loc[df.sta==sta,'lat'] = df[df.sta==sta].lat.values[0]
                df.loc[df.sta==sta,'time'] = df[df.sta==sta].time.values[0]
            
            # Check for outliers in the units.
            if 'DO_units' in df.columns:
                DO_units = [item for item in df.DO_units.unique() if item != None]
                if set(DO_units) > set(['umol/kg', 'mL/L']):
                    print('DO units problem: %s' % (DO_units))
            if 'Chl_units' in df.columns:
                Chl_units = [item for item in df.Chl_units.unique() if item != None]
                if set(Chl_units) > set(['mg/m^3']):
                    print('Chl units problem: %s' % (Chl_units))
            if 'N_units' in df.columns:
                N_units = [item for item in df.N_units.unique() if item != None]
                if set(N_units) > set(['umol/L']):
                    print('N units problem: %s' % (N_units))
            if 'Si_units' in df.columns:
                Si_units = [item for item in df.Si_units.unique() if item != None]
                if set(Si_units) > set(['umol/L','mmol/m**3']):
                    print('Si units problem: %s' % (Si_units))
            # RESULT: there are no outliers beyond out lists.
            
            # fix NO3 units
            if 'N' in df.columns:
                df['NO3 (uM)'] = df['N']
            else:
                df['NO3 (uM)'] = np.nan
        
            # fix DO units
            if 'DO' in df.columns:
                # calculate in situ density
                p = gsw.p_from_z(df.z.values, df.lat.values)
                rho = gsw.rho(df['salt (SA g kg-1)'].values, df['temp (CT degC)'].values, p)
                df['rho'] = rho
                df['DO (uM)'] = df['DO']
                df.loc[df.DO_units=='umol/kg','DO (uM)'] *= df.loc[df.DO_units=='umol/kg','rho']/1000
                df.loc[df.DO_units=='mL/L','DO (uM)'] *= 1.42903
            else:
                df['DO (uM)'] = np.nan
                
            # Fix Si name and handle a few units outliers
            df['Si (uM)'] = np.nan
            if 'Si' in df.columns:
                mask = (df['Si_units']=='umol/L') | (df['Si_units']=='mmol/m**3')
                df.loc[mask,'Si (uM)'] = df.loc[mask,'Si']
            # Note that in 1975 Si_units also showed up as "microg-at/l" but
            # we just won't use those.
            
            # Fix Chl name
            if 'Chl' in df.columns:
                df['Chl (mg m-3)'] = df['Chl']
            else:
                df['Chl (mg m-3)'] = np.nan
            
            # clean up columns
            df = df[['sta', 'lon', 'lat', 'time', 'z',
                'salt (SA g kg-1)', 'temp (CT degC)',
                'DO (uM)', 'NO3 (uM)', 'Si (uM)', 'Chl (mg m-3)']]

            # save
            df.to_pickle(out_fn)
    



