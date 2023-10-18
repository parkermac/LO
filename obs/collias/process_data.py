"""
Code to process the Collias data.

The first time you run it use this flag:
run process_data -preprocess True

After that you can just do:
run process_data
(much faster)

Here are the unique values of 'Result_Parameter_Name' and 'Result_Value_Units':
                      Salinity [ppt]
    Alkalinity, Total as CaCO3 [mg/L]
               Ortho-Phosphate [mg/L]
                       Nitrate [mg/L]
            Temperature, water [deg C]
              Dissolved Oxygen [mL/L]
         Silicic acid (H4SiO4) [mg/L]
              Nitrogen dioxide [mg/L]
                       Density [mg/L]

Years with data are 1932 to 1975, with some gaps

"""

import pandas as pd
import numpy as np
import gsw
import sys

from lo_tools import Lfun, obs_functions
Ldir = Lfun.Lstart()

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-preprocess', default=False, type=Lfun.boolean_string)
parser.add_argument('-testing', default=False, type=Lfun.boolean_string)
args = parser.parse_args()
testing = args.testing

# BOTTLE
source = 'collias'
otype = 'bottle'
in_dir = Ldir['data'] / 'obs' / source
# year_list = range(2008,2019)

# output location
out_dir = Ldir['LOo'] / 'obs' / source / otype
Lfun.make_dir(out_dir)

cols = [
    'Location_ID',
    'Field_Collection_Start_Date_Time',
    'Field_Collection_Upper_Depth',
    # 'Field_Collection_Depth_Units', # this is [m]
    # 'Sample_ID',
    'Result_Parameter_Name',
    'Result_Value',
    'Result_Value_Units',
    'Calculated_Latitude_Decimal_Degrees_NAD83HARN',
    'Calculated_Longitude_Decimal_Degrees_NAD83HARN'
]
 
# dict for renaming columns
col_dict = {
    'Location_ID':'name',
    'Field_Collection_Start_Date_Time':'time',
    'Field_Collection_Upper_Depth':'depth',
    'Result_Parameter_Name':'vn',
    'Result_Value':'val',
    'Result_Value_Units':'units',
    'Calculated_Latitude_Decimal_Degrees_NAD83HARN':'lat',
    'Calculated_Longitude_Decimal_Degrees_NAD83HARN':'lon'
}

if args.preprocess:
    fn = in_dir / 'EIMDiscreteResults_2023May17_456416.csv'
    df = pd.read_csv(fn, low_memory=False, usecols=cols,
        parse_dates = ['Field_Collection_Start_Date_Time'])
    df = df.rename(columns=col_dict)
    df['z'] = -df.depth
    df = df.drop(columns=['depth'])
    df.to_pickle(in_dir / 'preprocessed.p')
else:
    df = pd.read_pickle(in_dir / 'preprocessed.p')
    
# ORGANIZE INTO CASTS

# set index so we can loop over years
df = df.set_index('time')

sta_list = list(df.name.unique())

# loop over years
if testing:
    year_list = [1966]
else:
    year_list = range(1932,1976)
    
for year in year_list:
    
    # initialize a DataFrame to hold all data for a year
    ydf = df[df.index.year==year].copy()
    Ydf = pd.DataFrame()
    
    if len(ydf) == 0:
        print('No data in ' + str(year))
        continue
    else:
        print('Working on ' + str(year))

    # loop over stations
    for sn in sta_list:
        
        sdf = ydf[ydf.name==sn].copy()
        if len(sdf) == 0:
            continue
        
        Sdf = pd.DataFrame()
        
        # We add time column because we drop the index (time) below.
        sdf['time'] = sdf.index.copy()
        
        # Initialize DataFrame that holds all cleaned casts at this station.
        # Cleaned means no repeat depths
        cdf = pd.DataFrame()
        
        print_info = True
        ncast = 0
        # loop over dates (meaning individual casts)
        for dd in sdf.index.unique(): # Note: the index is still time at this point
            sdft = sdf[sdf.index==dd].copy() # just data from this time
            if len(sdft) == 0:
                continue
                
            sdfz = sdft.set_index('z')
            # drop repeated values (duplicate depths)
            sdfzu = sdfz[~sdfz.index.duplicated()]
            sdfzu = sdfzu.sort_index() # sort deepest to shallowest
            if (len(sdfz) != len(sdfzu)) and print_info:
                print('%s %s dropped %d repeat bottles' % (sn, str(dd), len(sdfz)-len(sdfzu)))
                
            if len(sdfzu) == 0:
                continue
            Sdf = pd.concat((Sdf,sdfzu)) # this has "z" as its index
            ncast += 1
            
        if print_info and (len(sdfzu)>0):
            print('*** %s There were %d casts at this station ***\n' % (sn, ncast))

        Sdf['z'] = Sdf.index.copy() # save because we drop it in the concat below
        Ydf = pd.concat((Ydf, Sdf), ignore_index=True, sort=False)
        # The final DataFrame just has numbers for its index.
        
    # save result
    # Bottles.to_pickle(dir0 + 'Bottles_' + str(year) + '.p')


if False:
    
    # Check on data and units
    rpn_set = set(df.vn.to_list())
    for rpn in rpn_set:
        unit_set = set(df.loc[df.vn==rpn,'units'].to_list())
        for unit in unit_set:
            print('%30s [%s]' % (rpn,unit))
    
    # Just look at one variable, trying to reconstruct things
    vn = 'Temperature, water'
    this_df = df.loc[df.vn==vn,:]

    lon = this_df.lon.to_numpy()
    lat = this_df.lat.to_numpy()
    z = this_df.z.to_numpy()
    t = this_df.index
    # tu = t.sort_values().unique() # could this be used to define "casts"?
    
    # Plotting
    import matplotlib.pyplot as plt
    from lo_tools import plotting_functions as pfun
    plt.close('all')
    pfun.start_plot(figsize=(12,8))

    fig = plt.figure()

    ax = fig.add_subplot(211)
    ax.plot(lon, lat,'.b')
    pfun.add_coast(ax)
    pfun.dar(ax)

    ax.axis([-125, -122, 47, 49])
    ax = fig.add_subplot(212)
    ax.plot(t,z,'.b')

    plt.show()
    pfun.end_plot()
