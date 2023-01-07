"""
Code to process DFO ctd data.

Runs in about [] minutes.

NEED to update to the same standards as process_bottle.py.

"""
from datetime import datetime, timedelta
import numpy as np
import gsw
import sys
import pandas as pd
from time import time

import loadDFO
from importlib import reload
reload(loadDFO)

from lo_tools import Lfun
Ldir = Lfun.Lstart()
in_dir = Ldir['data'] / 'obs' / 'dfo'

testing = True

if testing:
    year_list = [2017]
else:
    year_list = list(range(2013, datetime.now().year + 1))
out_dir = Ldir['LOo'] / 'obs' / 'dfo' / 'ctd'
Lfun.make_dir(out_dir)

tt0 = time()
for year in year_list:
    print(year)
    datelims = (datetime(year,1,1), datetime(year,12,31))
    out_fn = out_dir / (str(year) + '.p')
    info_out_fn = out_dir / ('info_' + str(year) + '.p')

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
                # the Station field is a unique cast identifier. ??

                # Force certain fields to be the same throughout the cast.
                df.loc[df.cid==cid,'lon'] = df[df.cid==cid].lon.values[0]
                df.loc[df.cid==cid,'lat'] = df[df.cid==cid].lat.values[0]
                df.loc[df.cid==cid,'time'] = df[df.cid==cid].time.values[0]
        
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
            #     rho = gsw.rho(df['SA'].values, df['CT'].values, p)
            #     df['rho'] = rho
            #     df['DO (uM) alt'] = df['DO_umolkg'] * df['rho']/1000
            #     # check on how much there is
            #     print(' - isnull DO_mLL:DO_umolkg = %d:%d' %
            #         (df['DO (uM)'].isnull().sum(), df['DO (uM) alt'].isnull().sum()))
    
            # clean up columns
            df = df[['cid', 'lon', 'lat', 'time', 'z',
                'SA', 'CT',
                'DO (uM)', 'Fluor']]
            df['name'] = None
        
            # save
            df.to_pickle(out_fn)
        
            # # Also pull out a dateframe with station info to use for model cast extractions.
            ind = df.cid.unique()
            info_df = pd.DataFrame(index=ind, columns=['lon','lat','time','name'])
            for cid in df.cid.unique():
                info_df.loc[cid,['lon','lat','time']] = df.loc[df.cid==cid,['lon','lat','time']].iloc[0,:]
            info_df.name = None
            info_df.index.name = 'cid'
            info_df['time'] = pd.to_datetime(info_df['time'])
            info_df.to_pickle(info_out_fn)
        
print('Elapsed time for ctd = %0.1f sec' % (time()-tt0))



