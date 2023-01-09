"""
Code to process the ecology ctd data.

"""

import pandas as pd
import numpy as np
import gsw
import sys

from lo_tools import Lfun
Ldir = Lfun.Lstart()

# BOTTLE
source = 'ecology'
otype = 'ctd'
in_dir0 = Ldir['data'] / 'obs' / source

testing = True

if testing:
    year_list = [2017]
else:
    year_list = range(2008,2019)

# output location
out_dir = Ldir['LOo'] / 'obs' / source / otype
Lfun.make_dir(out_dir)

# This is a dict of all the columns after the initial reading.
# We add values to a key for any variable we want to save
v_dict = {

}

# We need to associate lat and lon with each station. They are not stored in the bottle
# file, but there are station names, with loccations here:
sta_fn = in_dir0 / 'ParkerMacCreadyCoreStationInfoFeb2018.xlsx'
sta_df = pd.read_excel(sta_fn, index_col='Station')
xx = sta_df['Long_NAD83 (deg / dec_min)'].values
yy = sta_df['Lat_NAD83 (deg / dec_min)'].values
lon = [-(float(x.split()[0]) + float(x.split()[1])/60) for x in xx]
lat = [(float(y.split()[0]) + float(y.split()[1])/60) for y in yy]
sta_df['lon'] = lon
sta_df['lat'] = lat

load_data = True
for year in year_list:
    ys = str(year)
    print('\n'+ys)
    
    # name output files
    out_fn = out_dir / (ys + '.p')
    info_out_fn = out_dir / ('info_' + ys + '.p')
    
    if year == 2017:
        load_data = True
        ctd_fn = dir0 + 'raw/ParkerMacCready2017CTDDataFeb2018.xlsx'
        sheet_name = '2017Provisional_CTDResults'
    elif year == 2018:
        load_data = True
        #ctd_fn = dir0 + 'raw/Parker_2018.xlsx'
        ctd_fn = dir0 + 'raw/ParkerMacCready2018CTDDOMar2020.xlsx'
        sheet_name = '2018_CTDDOResults'
    elif year == 2019:
        load_data = True
        ctd_fn = dir0 + 'raw/ParkerMacCready2019CTDDataFeb2020.xlsx'
        sheet_name = '2019Provisional_CTDResults'
    else:
        ctd_fn = dir0 + 'raw/ParkerMacCready1999-2016CTDDataMay2018.xlsx'
        sheet_name = '1999-2016Finalized_CTDResults'

    # data long names; we retain only these fields
    
    # DEFAULTS
    date_col_name = 'Date'
    station_col_name = 'Station'
    depth_col_name = 'Depth'
    data_new_names = ['Salinity', 'Temperature', 'DO', 'Z']
    
    if year == 2017:
        data_original_names = ['Salinity', 'Temp', 'DO_raw', 'Z']
    elif year == 2018: # missing Chl
        #date_col_name = 'UTCDate'
        #station_col_name = 'SiteCode'
        #depth_col_name = 'ActualDepthDecimal'
        data_original_names = ['Salinity', 'Temp', 'DO_adjusted', 'Z']
        
    elif year == 2019:
        data_original_names = ['Salinity', 'Temp', 'DO_raw', 'Z']
    else:
        data_original_names = ['Salinity', 'Temp', 'DO_adjusted', 'Z']
    
    # read in the data (all stations, all casts)
    if load_data:
        all_casts = pd.read_excel(ctd_fn, sheet_name=sheet_name,
            parse_dates = [date_col_name])
        load_data = False
    
