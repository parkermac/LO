"""
Code to process the CTD data collected from Olympic coast national marine
sanctuary, WA.

This takes just seconds to process 2004-2023.

"""

import pandas as pd
import numpy as np
import gsw
from time import time as Time
from pathlib import Path
import xarray as xr

from lo_tools import Lfun, obs_functions
Ldir = Lfun.Lstart()

# CTD
source = 'ocnms_ctd'
folder_list = ['ctd01','ctd02'] # from two separate sources
otype = 'ctd'

testing = False
if testing:
    year_list = [2005]
    sn_list = ['CA015','CA042','CA065','Cape_Alava']

else:

    year_list = range(2004,2024)
    sn_list = ['CA015', 'CA042', 'CA065', 'CE015', 'CE042', 'CE065',
               'KL015', 'KL027', 'KL050', 'MB015', 'MB042', 'TH015',
               'TH042', 'TH065', 'Cape_Alava', 'Cape_Elizabeth',
               'Hoh_Head', 'Moclips', 'Raft_River', 'Teahwhit_Head']

# output location
out_dir = Ldir['LOo'] / 'obs' / source / otype

Lfun.make_dir(out_dir)

# process data
tt0 = Time()
for year in year_list:
    ys = str(year)
    print('\n'+ys)
    # name output files
    out_fn = out_dir / (ys + '.p')
    info_out_fn = out_dir / ('info_' + ys + '.p')
    fn_list = []
    for sn in sn_list:
        # print(sn)
        for folder in folder_list:
            # print(folder)
            in_dir0 = Ldir['data'] / 'obs' / source / folder/ 'netcdf_files'/ 'binned_profiles'
            in_dir = Path(in_dir0 / sn)
            # print(in_dir)
            in_file = list(in_dir.glob(f'*{year}*.nc'))
            if not in_file:
                pass
                # print("The list is empty.")
            else:
                fn_list.extend(in_file)
    # print(fn_list)
    df0 = []
    n = 0
    for file in fn_list:
        # print(file)
        # get station name
        sn = Path(file).parent.stem
        # print(sn)
        ds = xr.open_dataset(file)
        # check if year is matching
        time = ds['time'].values
        yr = time.astype('datetime64[Y]').astype(int) + 1970
        if yr == year:
            p = ds.pressure.values
            cid = np.repeat(n, len(p))
            cruise = np.repeat(source, len(p))
            time = np.repeat(ds['time'].values, len(p))
            lat = np.repeat(ds.latitude.values, len(p))
            lon = np.repeat(ds.longitude.values, len(p))
            name = np.repeat(sn, len(p))
            z = ds.depth.values * -1
            IT = ds['temperature'].values
            SP = ds['salinity'].values
            # - do the conversions
            DO = ds['dissolved_oxygen'].values * 44.661 # from ml/L to umol/L
            SA = gsw.SA_from_SP(SP, p, lon, lat)
            CT = gsw.CT_from_t(SA, IT, p)
            lon = lon - 360
            data = {'cid': cid,
                    'cruise': cruise,
                    'time': time,
                    'lat': lat,
                    'lon': lon,
                    'name': name,
                    'z': z,
                    'DO (uM)': DO,
                    'SA': SA,
                    'CT': CT
                    }
            df1 = pd.DataFrame(data)
            # missing data is 9999.0
            df1[df1 == 9999.0] = np.nan
            # a little more cleaning up
            df1  = df1.dropna(axis=0, how='all') # drop rows with no good data
            df1  = df1[df1.time.notna()] # drop rows with bad time
            df1 = df1.reset_index(drop=True)
            df0.append(df1)
            n += 1
    print(f'processed {n+1} casts')
    df = pd.concat(df0, ignore_index=True)
    if len(df) > 0:
        # Save the data
        df.to_pickle(out_fn)
        info_df = obs_functions.make_info_df(df)
        info_df.to_pickle(info_out_fn)

print('Total time = %d sec' % (int(Time()-tt0)))


