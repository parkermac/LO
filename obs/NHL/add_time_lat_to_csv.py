"""
NOTE: There are multiple files for a given year, and the datatime is in the "name" of each file, not within the file as a variable,
so you have to extract the datatime info while reading the filename.
'lat' variable is missing (which should be 44.652, because NHL is a fixed line)
'cruise' should be 'NH'
Also, you want to combine those monthly files to yearly
"""


import pandas as pd
import glob
from datetime import datetime
from lo_tools import Lfun
Ldir = Lfun.Lstart()

# CTD
source = 'NHL'
otype = 'ctd'
in_dir0 = Ldir['data'] / 'obs' / source /'newport_line_ctd_raw'
# output location
out_dir0 = Ldir['data'] / 'obs' / source / 'ctd_biweekly'
Lfun.make_dir(out_dir0)

v_dict = {
    'NHL station number':'station_name',
    ' longitude (degW)':'lon',
    ' pressure (dbar)':'pressure (dbar)',
    'temperature (degC)':'temperature (degC)',
    ' practical salinity':'salinity',
    ' dissolved oxygen (ml/L)':'DO (ml/L)',
}

in_fn= list(in_dir0.glob('*.csv'))
for file in in_fn:
    df0 = pd.read_csv(file, skiprows=1, usecols=list(v_dict.keys()))
    df0.rename(columns=v_dict, inplace=True)
    datetime_str = file.stem
    df0['datetime'] = datetime.strptime(datetime_str, '%Y%m%dT%H%M%S')
    df0['lat'] = 44.652
    df0['cruise'] ='NH'
    columns = list(df0.columns)
    this_cols = ['datetime', 'lat', 'lon'] + [col for col in columns if col not in ['datetime', 'lat', 'lon']]
    df0 = df0[this_cols]
    df0.to_csv(out_dir0/(file.stem + '.csv'), index=False)

testing = False

if testing:
    year_list = [2017]
else:
    year_list = range(1997,2022)
in_dir = out_dir0
# output location
out_dir = Ldir['data'] / 'obs' / source / otype
Lfun.make_dir(out_dir)

for year in year_list:
    ys = str(year)
    print('\n'+ys)

    if year in year_list:
        #NOTE: There are multiple files for a given year
        in_fn = in_dir.glob(f'*{year}*.csv')
        df0 = pd.DataFrame()
        for files in in_fn:
            df = pd.read_csv(files)
            df0 = pd.concat([df0, df], ignore_index=True)
            df0.to_csv(out_dir/(f'{year}' + '.csv'), index=False)