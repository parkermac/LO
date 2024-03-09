"""
NOTE: There are multiple files for a given year, we want to combine those files to yearly for easier process
one file in 2008 is malformed, so I manually copied it to a new .csv file, but did not change any data, it works.
I manually changed the following files to have consistent data variables at that year while using Pandas:

18DD20150210_export_whp.csv
18DD20150607_hy1.csv
18DD20160817_hy1.csv
18DD20030527_hy1.csv
18DD20030831_hy1.csv
18DD20010618_hy1.csv
18DD20010823_hy1.csv

original files are in the folder "about"

"""

import pandas as pd
import glob
from lo_tools import Lfun
Ldir = Lfun.Lstart()

# This is a dict of all the columns after the initial reading.
# We add values to a key for any variable we want to save. I looked
# at units_dict (created below) to be sure about the units.
v_dict = {
    'time':'time',
    'EXPOCODE':'',
    'SECT_ID':'cruise',
    'SECT':'cruise',
    'DATE':'',
    'TIME':'',
    'LATITUDE':'lat',
    'LONGITUDE':'lon',
    'STNNBR':'name',
    'STNNB':'name',
    'CAST':'cid',
    'CASTNO':'cid',
    'BTLNBR':'',
    'BTLNBR_FLAG_W':'',
    'SAMPNO':'',
    'DEPTH':'',
    'CTDPRS':'P (dbar)',
    'CTDTMP':'IT',
    'CTDSAL':'SP',
    'CTDSAL_FLAG_W':'',
    'CTDOXY':'',
    'SALNTY':'',
    'SALNTY_FLAG':'',
    'OXYGEN':'DO (umol/kg)',
    'OXYGEN_FLAG':'',
    'OXYGEN_FLAG_W':'',
    'NO2+NO3':'NO3 (umol/kg)',
    'NO3+NO2':'NO3 (umol/kg)',
    'NO2+NO3_2':'NO3 (umol/kg)',
    'NO2+NO3_FLAG':'',
    'NO2+NO3_FLAG_W':'',
    'NITRAT':'NO3 (umol/kg)',
    'NITRAT_FLAG_W':'',
    'SILCAT':'SiO4 (umol/kg)',
    'SILCAT_2':'SiO4 (umol/kg)',
    'SILCAT_FLAG_W':'',
    'PHSPHT':'PO4 (umol/kg)',
    'PHSPHT_2':'PO4 (umol/kg)',
    'PHSPHT_FLAG_W':'',
    'TCARBN':'DIC (umol/kg)',
    'TCARBN_FLAG':'',
    'TCARBN_FLAG_W':'',
    'ALKALI':'TA (umol/kg)',
    'ALKALI_FLAG':'',
    'ALKALI_FLAG_W':'',
    'PH_TOT':'PH',
    'PH_TOT_FLAG_W':'',
}

# bottle
source = 'nceiPNW'
otype = 'bottle'
in_dir = Ldir['data'] / 'obs' / source

testing = False

if testing:
    year_list = [2001]
else:
    # 1985-2000 have same patterns
    year_list = [1985, 1986, 1989] + list(range(1991, 2008+1))+ [2010, 2013, 2014, 2015, 2016, 2017]
# output location
out_dir = Ldir['data'] / 'obs' / source / otype
Lfun.make_dir(out_dir)

for year in year_list:
    ys = str(year)
    print('\n'+ys)

    in_fn = in_dir.glob(f'*18DD{year}*.csv')
    df0 = pd.DataFrame()

    for files in in_fn:
        skiprows = []
        drop_last_row = True
        cruise = False
        time = False
        date = False
        sal = False
        if year ==2001:
            skiprows = [0]
            date = True
        elif year == 2002:
            skiprows = range(0, 8)
        elif year ==2003:
            date = True
            time = True
            drop_last_row = False
        elif year == 2010:
            skiprows = range(0, 21)
            sal = True
        elif year == 2014:
            skiprows = range(0, 7)
            sal = True
        elif year == 2013 or year ==2015 or year == 2016 or year ==2017:
            drop_last_row = False
            time = True
        elif year in range(2004, 2009):
            skiprows = range(0, 8)
            cruise = True
        else:
            skiprows = [0]

        df = pd.read_csv(files, skiprows=skiprows, dtype={'DATE': 'object', 'TIME': 'object'})
        df = df.drop([0])  # drop the first row for unit
        if drop_last_row:
            df = df[:-1]  # drop the last row {END_DATA}
        if date:
            df['DATE'] = pd.to_datetime(df['DATE'], format='mixed')
            df['DATE'] = df['DATE'].dt.strftime('%Y%m%d')
        if cruise:
            df['SECT'] = 'Line P'
        if sal:
            df['CTDSAL'] = df['SALNTY']
        datetime_str = df['DATE'] + ' ' + df['TIME'].str.zfill(4)

        if time:
            df['time'] = pd.to_datetime(datetime_str, format='%Y%m%d %H:%M')
        else:
            df['time'] = pd.to_datetime(datetime_str, format='%Y%m%d %H%M')
            df['time'] = df['time'].dt.strftime('%Y-%m-%dT%H:%M')
        df0 = pd.concat([df0, df], ignore_index=True)

    # select and rename variables
    df1 = pd.DataFrame()
    for v in df0.columns:
        if v in v_dict.keys():
            if len(v_dict[v]) > 0:
                df1[v_dict[v]] = df0[v]
    # re-order variables
    columns = list(df1.columns)
    df1 = df1[['time', 'cruise'] + [col for col in df1.columns if col not in ['time', 'cruise']]]
    df1['lon'] = df1.lon.astype('float64')
    df1['lat'] = df1.lat.astype('float64')
    df1 = df1[(df1.lon > -130) & (df1.lon < -124) & (df1.lat > 42) & (df1.lat < 52)]
    df1.to_csv(out_dir/(f'{year}.csv'), index=False)
