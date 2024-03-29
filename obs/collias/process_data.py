"""
Code to process the Collias data.

The first time you run it use this flag:
run process_data -preprocess True

After that you can just do:
run process_data
(much faster)

Performance: 5 minutes for all years.

"""

import pandas as pd
import numpy as np
import gsw
from time import time

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

# output location
out_dir = Ldir['LOo'] / 'obs' / source / otype
Lfun.make_dir(out_dir)

# Columns to keep
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

# Read in data
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
print_info = False

if testing:
    print_info = True
    year_list = [1932]
    sta_list = ['DAB523']
    pd.set_option('display.max_rows', 500)
    pd.set_option('display.max_columns', 500)
    pd.options.display.width = 0 # auto-detect full display width
    
# Loop over years
tt00 = time()
for year in year_list:
    tt0 = time()
    ys = str(year)
    print('\n'+ys)
    
    # name output files
    out_fn = out_dir / (ys + '.p')
    info_out_fn = out_dir / ('info_' + ys + '.p')
    
    # initialize cast ID
    cid = 0
    
    # Initialize a DataFrame to hold all data for a year
    Ydf = pd.DataFrame()
    
    # pull out data for this year
    ydf = df[df.index.year==year].copy()
    
    if len(ydf) == 0:
        print('- No data in ' + str(year))
        continue # this handy command just skips to the next item in the for loop

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
            
            # The initial file has all variables in one column. We use the DataFrame "A"
            # to reorganize it so that each data variable has its own column, and
            # all data at a given depth is on one row.
            A = pd.DataFrame(index=sdfz.index.unique().sort_values())
            for vn in vn_dict.keys():
                # Pull in each variable as a Series
                aa = sdfz.loc[sdfz.vn==vn,'val'].copy()
                if len(aa) > 0:
                    # Remove duplicate depths, and sort by z (deep to shallow)
                    aa = aa[~aa.index.duplicated()]
                    aa = aa.sort_index()
                    A[vn_dict[vn]] = aa
                else:
                    A[vn_dict[vn]] = np.nan
                    
            # remove rows with no data
            A = A.dropna(axis=1,how='all')
            if len(A) == 0:
                continue
            else:
                # ==============================================================
                # NOTE: Deal with an issue where in many years, at some stations
                # (especially in Hood Canal)
                # the DO column was mistakenly a reapeat of the salinity column
                if ('salt (ppt)' in A.columns) and ('DO (mL/L)' in A.columns):
                    mask = A['salt (ppt)'] == A['DO (mL/L)']
                    A.loc[mask,'DO (mL/L)'] = np.nan
                # ==============================================================
                # NOTE: Deal with an issue where in some years, at some stations
                # the SiO4 and NO2 columns are the same
                if ('SiO4 (mg/L)' in A.columns) and ('NO2 (mg/L)' in A.columns):
                    mask = A['SiO4 (mg/L)'] == A['NO2 (mg/L)']
                    A.loc[mask,'NO2 (mg/L)'] = np.nan
                    A.loc[mask,'SiO4 (mg/L)'] = np.nan
                # ==============================================================
                
                # Then add remaining columns. This ensures that in the final product
                # each cast has these values tha same for all depths.
                bb = sdfz.iloc[0,:]
                for vn in ['name','time','lon','lat']:
                    A.loc[:,vn] = bb[vn]
                # Add a unique cast ID (cid) for each cast. This will only be unique
                # within the year being processed.
                A['cid'] = cid
                Sdf = pd.concat((Sdf,A))
                # Note: this has "z" as its index, with many repeats
                ncast += 1
                cid += 1
            
        if print_info and (len(Sdf)>0):
            print('*** %s There were %d casts at this station ***' % (sn, ncast))

        Sdf['z'] = Sdf.index.copy()
        # Duplicate z as a column because we drop it in the concat below
        
        # Append this station's casts to the year DataFrame
        Ydf = pd.concat((Ydf, Sdf), ignore_index=True, sort=False)
        # The final DataFrame just has numbers for its index.
        
    Ydf['cruise'] = None # no cruise info for this data
    
    # Fix units and generate CT and SA
    # SA and CT
    SP = Ydf['salt (ppt)'].to_numpy() # assume this is close enough to psu
    IT = Ydf['temp (degC)'].to_numpy() # I assume this is in-situ
    z = Ydf.z.to_numpy()
    lon = Ydf.lon.to_numpy()
    lat = Ydf.lat.to_numpy()
    p = gsw.p_from_z(z, lat)
    # - do the conversions
    SA = gsw.SA_from_SP(SP, p, lon, lat)
    CT = gsw.CT_from_t(SA, IT, p)
    # - add the results to the DataFrame
    Ydf['SA'] = SA
    Ydf['CT'] = CT
    # Other variables
    fac_dict = {
        'DO (mL/L)': 1.42903 * 1000 / 32,
        'NO3 (mg/L)': 1000 / 14,
        'NO2 (mg/L)': 1000 / 14,
        'SiO4 (mg/L)': 1000 / 28.0855
    }
    name_dict = {
        'DO (mL/L)': 'DO (uM)',
        'NO3 (mg/L)': 'NO3 (uM)',
        'NO2 (mg/L)': 'NO2 (uM)',
        'SiO4 (mg/L)': 'SiO4 (uM)'
    }
    for vnn in fac_dict.keys():
        if vnn in Ydf.columns:
            Ydf[name_dict[vnn]] = Ydf[vnn].to_numpy() * fac_dict[vnn]
        else:
            pass
    
    # Retain only selected variables
    cols = ['cid', 'cruise', 'time', 'lat', 'lon', 'name', 'z',
        'CT', 'SA', 'DO (uM)',
        'NO3 (uM)', 'NO2 (uM)', 'NH4 (uM)', 'PO4 (uM)', 'SiO4 (uM)',
        'TA (uM)', 'DIC (uM)']
    this_cols = [item for item in cols if item in Ydf.columns]
    Ydf = Ydf[this_cols]
    
    # Save results
    if (len(Ydf) > 0) and (testing==False):
        Ydf.to_pickle(out_fn)
        info_df = obs_functions.make_info_df(Ydf)
        info_df.to_pickle(info_out_fn)
        
    print('- took %0.1f sec to process' % (time()-tt0))

total_minutes = (time()-tt00)/60
print('Total processing time %0.1f minutes' % (total_minutes))
 