"""
Code to process DFO cast data.

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
    year_list = [2017]
else:
    year_list = list(range(2013, datetime.now().year + 1))

for year in year_list:
    print(year)
    datelims = (datetime(year,1,1), datetime(year,12,31))
    out_fn = out_dir / ('casts_' + str(year) + '.p')
    
    df = loadDFO.loadDFO_CTD(basedir=str(in_dir), datelims=datelims)
    """
    Expected Fields:
    ['Station', 'Lat', 'Lon', 'Z', 'SA', 'CT', 'Fluor', 'DO_mLL',
           'DO_umolkg', 'dtUTC']
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
                # the Station field is a unique cast identifier. ??

                # Force certain fields to be the same throughout the cast.
                df.loc[df.sta==sta,'lon'] = df[df.sta==sta].lon.values[0]
                df.loc[df.sta==sta,'lat'] = df[df.sta==sta].lat.values[0]
                df.loc[df.sta==sta,'time'] = df[df.sta==sta].time.values[0]
            
            # Check for outliers in the units.
            if 'DO_units' in df.columns:
                DO_units = [item for item in df.DO_units.unique() if item != None]
                if set(DO_units) > set(['DO_mLL','DO_umolkg']):
                    print('DO units problem: %s' % (DO_units))
            # RESULT: No outliers.
            
            # fix DO units
            df['DO (uM)'] = np.nan
            df['DO (uM) alt'] = np.nan
            if 'DO_mLL' in df.columns:
                df['DO (uM)'] = (1000/32) * 1.42903 * df['DO_mLL']
            # We will skip DO_umolkg this because it appears that (i) it is basically identical
            # to what we get from DO_mLL, and (ii) there are always more missing values for DO_umolkg.
            # if 'DO_umolkg' in df.columns:
            #     # calculate in situ density
            #     p = gsw.p_from_z(df.z.values, df.lat.values)
            #     rho = gsw.rho(df['salt (SA g kg-1)'].values, df['temp (CT degC)'].values, p)
            #     df['rho'] = rho
            #     df['DO (uM) alt'] = df['DO_umolkg'] * df['rho']/1000
            #     # check on how much there is
            #     print(' - isnull DO_mLL:DO_umolkg = %d:%d' %
            #         (df['DO (uM)'].isnull().sum(), df['DO (uM) alt'].isnull().sum()))
        
            # clean up columns
            df = df[['sta', 'lon', 'lat', 'time', 'z',
                'salt (SA g kg-1)', 'temp (CT degC)',
                'DO (uM)', 'Fluor']]
            
            # save
            df.to_pickle(out_fn)