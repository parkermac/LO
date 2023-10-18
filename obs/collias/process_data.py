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
parser.add_argument('-test','--testing', default=False, type=Lfun.boolean_string)
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

# dict for renaming data variables (only save this list)
vn_dict = {
    'Salinity': 'salt (ppt)',
    'Nitrate': 'NO3 (mg/L)',
    'Temperature, water': 'temp (degC)',
    'Dissolved Oxygen': 'DO (mL/L)',
    'Silicic acid (H4SiO4)': 'SiO4 (mg/L)',
    'Nitrogen dioxide': 'NO2 (mg/L)'
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

# default lists
year_list = range(1932,1976)
sta_list = list(df.name.unique())

if testing:
    print_info = True
    year_list = [1966]
    # sta_list = ['LCH551']
    pd.set_option('display.max_rows', 500)
    pd.set_option('display.max_columns', 500)
    pd.options.display.width = 0 # auto-detect full display width
    
# loop over years
for year in year_list:
    
    # initialize cast ID
    cid = 0
    
    # Initialize a DataFrame to hold all data for a year
    Ydf = pd.DataFrame()
    
    # pull out data for this year
    ydf = df[df.index.year==year].copy()
    
    if len(ydf) == 0:
        print('No data in ' + str(year))
        continue # this handy command just skips to the next item in the for loop
    else:
        print('Working on ' + str(year))

    # loop over stations
    for sn in sta_list:
        
        # Initialize a DataFrame to hold all data for this station
        Sdf = pd.DataFrame()
        
        # pull out data for this station
        sdf = ydf[ydf.name==sn].copy()
        if len(sdf) == 0:
            continue
        
        # We add a duplicate time column because we swap to z for the index below.
        sdf['time'] = sdf.index.copy()
                
        ncast = 0
        # Loop over dates
        # NOTE: we are ASSUMING an individual date corresponds to an individual cast!
        # How can we tell if this is true?
        for dd in sdf.index.unique(): # Note: the index is still time at this point
            sdft = sdf[sdf.index==dd].copy() # just data from this time
            if len(sdft) == 0:
                continue
                
            sdfz = sdft.set_index('z')
                
            A = pd.DataFrame(index=sdfz.index.unique().sort_values())
            for vn in vn_dict.keys():
                # pull in each variable as a Series
                aa = sdfz.loc[sdfz.vn==vn,'val'].copy()
                if len(aa) > 0:
                    aa = aa[~aa.index.duplicated()]
                    aa = aa.sort_index()
                    A[vn_dict[vn]] = aa
                else:
                    A[vn_dict[vn]] = np.nan
            
            # Then add remaining columns
            bb = sdfz.iloc[0,:]
            for vn in ['name','time','lon','lat']:
                A.loc[:,vn] = bb[vn]
                
            if len(A) == 0:
                continue
            else:
                A['cid'] = cid
                ncast += 1
                Sdf = pd.concat((Sdf,A)) # this has "z" as its index
                cid += 1
            
        if print_info and (len(A)>0):
            print('*** %s There were %d casts at this station ***\n' % (sn, ncast))

        Sdf['z'] = Sdf.index.copy() # save because we drop it in the concat below
        
        # Append this station's casts to the year DataFrame
        Ydf = pd.concat((Ydf, Sdf), ignore_index=True, sort=False)
        # The final DataFrame just has numbers for its index.
        
    Ydf['cruise'] = None
    
    # still neet to fix units and generate CT and SA
    
    # if len(Ydf) > 0:
    #     # Save the data
    #     df.to_pickle(out_fn)
    #     info_df = obs_functions.make_info_df(df)
    #     info_df.to_pickle(info_out_fn)
    #
    #
    # # save result
    # Ydf.to_pickle(dir0 / ('Bottles_' + str(year) + '.p'))

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
