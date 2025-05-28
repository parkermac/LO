"""
Note that this script is only compatible with trapsD01 (listed in traps_data_ver.csv)

This script compiles all of the excel files from Mohamedali et al. (2020)
and the csv files from Wasielewski et al. (2024) into three netCDF files
which are saved in LO_data/trapsD01/processed_data:
    river_data_mohamedali_etal_2020.nc
    wwtp_data_mohamedali_etal_2020.nc
    wwtp_data_wasielewski_etal_2024.nc

In theory, this script only needs to be run once.
Then, the netCDF files can be referenced to generate climatologies

Takes about 15 minutes to run on my local machine.

To run from ipython:
    run rawdata_2netCDF.py
"""

#################################################################################
#                              Import packages                                  #
#################################################################################
from lo_tools import Lfun
Ldir = Lfun.Lstart()

import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import cmocean
import os
import datetime
import xarray as xr
from pathlib import Path
from dateutil.relativedelta import relativedelta

from lo_tools import plotting_functions as pfun

plt.close('all')

#################################################################################
#                              Helper functions                                 #
#################################################################################

# Helper functions are designed to specifically process data from two datasets.
# Those designed to process data from Mohamedali et al. (2020) are prefaced by "moh20"
# Those designed to process data from Wasielewski et al. (2024) are prefaced by "was24"

def start_ds(date,source,numdates,Nsources,riv_or_wwtp,data_source):
    '''
    Initialize dataset to store loading data
    '''
    if data_source == 'moh20':
        dataset = 'Mohamedali et al. (2020)'
    elif data_source == 'was24':
        dataset = 'Wasielewski et al. (2024)'
    else:
        raise Exception('data_source must be either `moh20` or `was24`')

    ds = xr.Dataset(data_vars=dict(
        ID   =(['source'], ['placeholder placeholder placeholder']*Nsources),
        lon  =(['source'], np.ones((Nsources,))),
        lat  =(['source'], np.ones((Nsources,))),
        name =(['source'], ['placeholder placeholder placeholder']*Nsources),
        flow =(['source', 'date'], np.zeros((Nsources, numdates))),
        temp =(['source', 'date'], np.zeros((Nsources, numdates))),
        NO3  =(['source', 'date'], np.zeros((Nsources, numdates))),
        NH4  =(['source', 'date'], np.zeros((Nsources, numdates))),
        TIC  =(['source', 'date'], np.zeros((Nsources, numdates))),
        Talk =(['source', 'date'], np.zeros((Nsources, numdates))),
        DO   =(['source', 'date'], np.zeros((Nsources, numdates))),),
    coords=dict(source=source, date=date,),
    attrs=dict(description=dataset+' data for '+riv_or_wwtp),)
    
    return ds

def add_metadata(ds):
    '''
    Create metadata for loading data
    '''
    ds['lon'].attrs['long_name'] = 'source longitude'
    ds['lat'].attrs['long_name'] = 'source latitude'

    ds['flow'].attrs['long_name'] = 'discharge rate'
    ds['flow'].attrs['units'] = 'm3/s'

    ds['temp'].attrs['long_name'] = 'discharge temperature'
    ds['temp'].attrs['units'] = 'C'

    ds['NO3'].attrs['long_name'] = 'nitrate+nitrite concentration'
    ds['NO3'].attrs['units'] = 'mmol/m3'

    ds['NH4'].attrs['long_name'] = 'ammonium concentration'
    ds['NH4'].attrs['units'] = 'mmol/m3'

    ds['TIC'].attrs['long_name'] = 'total inorganic carbon'
    ds['TIC'].attrs['units'] = 'mmol/m3'

    ds['Talk'].attrs['long_name'] = 'total alkalinity'
    ds['Talk'].attrs['units'] = 'meq/m3'

    ds['DO'].attrs['long_name'] = 'dissolved oxygen concentration'
    ds['DO'].attrs['units'] = 'mmol/m3'
    return ds

def moh20_monthly2daily(df, start_date, end_date):
    '''
    turn a monthly dataframe into daily data (for the length of the raw data timeseries)
    '''
    # duplicate last row
    double_lr_df = pd.concat([df, df.iloc[-1:]], ignore_index=True)
    dates = pd.date_range(start_date, end_date, freq='MS')
    # Replace month column with things that look like dates
    double_lr_df['Month'] = dates
    double_lr_df = double_lr_df.set_index('Month')
    # Change monthly to daily, and fill daily values with repeats of monthly values
    double_lr_daily_df = double_lr_df.resample('D').ffill()
    # delete last row (Jan 1 on the next year)
    daily_df = double_lr_daily_df[:-1]
    # reset index
    daily_df.reset_index(inplace=True)
    return daily_df

def was24_monthly2daily(df, start_date, end_date):
    '''
    turn a monthly dataframe into daily data (for the length of the raw data timeseries)
    '''
    # duplicate last row
    double_lr_df = pd.concat([df, df.iloc[-1:]], ignore_index=True)
    # get dates (need to include one month beyond actual date range)
    dates = pd.date_range(start_date, end_date, freq='MS')
    # Replace month column with things that look like dates
    double_lr_df['date'] = dates
    double_lr_df = double_lr_df.set_index('date')
    # Change monthly to daily, and fill daily values with repeats of monthly values
    daily_df = double_lr_df.resample('D').ffill()
    # delete last row (Jan 1 on the next year)
    daily_df = daily_df[:-1]
    # reset index
    daily_df.reset_index(inplace=True)
    return daily_df

def moh20_wwtp_add_data(WWTP_index, ds, source_ID, source_name, moh20_latlon_df, moh20_wwtp_daily_df):
    '''
    Add Mohamedali et al. (2020) WWTP data to datasets, and convert to units that LO uses
    '''

    i = WWTP_index

    # Add source ID, name, and lat/lon
    ds['ID'][i] = source_ID
    ds['name'][i] = source_name
    ds['lat'][i] = moh20_latlon_df.loc[moh20_latlon_df['ID'] == source_ID, 'Lat'].values[0]
    ds['lon'][i] = moh20_latlon_df.loc[moh20_latlon_df['ID'] == source_ID, 'Lon'].values[0]
    
    # Add physics and biology data to dataset
    ds.flow[i,:] = moh20_wwtp_daily_df['Flow(m3/s)']
    ds.temp[i,:] = moh20_wwtp_daily_df['Temp(C)']
    ds.NO3[i,:]  = moh20_wwtp_daily_df['NO3+NO2(mg/L)'] * 71.4 # convert to mmol/m3
    ds.NH4[i,:]  = moh20_wwtp_daily_df['NH4(mg/L)']     * 71.4 # convert to mmol/m3
    ds.TIC[i,:]  = moh20_wwtp_daily_df['DIC(mmol/m3)']
    ds.Talk[i,:] = moh20_wwtp_daily_df['Alk(mmol/m3)']
    ds.DO[i,:]   = moh20_wwtp_daily_df['DO(mg/L)']      * 31.26 # convert to mmol/m3

    return ds

def get_moh20_wwtp_clim(moh20_wwtp_ds):
    '''
    Calculate climatology for temp, TIC, Talk, and DO using data from
    Mohamedali et al. (2020) for WWTPs

    output:
    list of climatologies in order: [temp, TIC, Talk, DO]
    '''

    # initialize empty list
    clims = []

    # create climatology for all variables, just arbitratily picking the first WWTP listed
    for var in ['temp','TIC','Talk','DO']:
        # generate climatology for leap years (366 days)
        leapyear_clim = moh20_wwtp_ds[var][0,:].groupby('date.dayofyear').mean(dim='date').values
        # remove leap day to get non-leap year climatology (365 days)
        feb29_index = 59
        nonleapyear_clim = np.delete(leapyear_clim, feb29_index)
        # correctly alternate leap and non-leap years to create artificial 2005-2020 time series
        leap = leapyear_clim   # renamed for convenience
        non = nonleapyear_clim # renamed for convenience
        clim_2005to2020 = np.concatenate((
        #   2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020
            non , non , non , leap, non , non , non , leap, non , non , non , leap, non , non , non , leap))
        # add 2005 - 2020 climatology to dataset
        clims.append(clim_2005to2020)

    return clims

def was24_wwtp_add_data(WWTP_index, ds, source_ID, source_name, lat, lon, was24_daily_df, moh20_clim_list):
    '''
    Add Wasielewski et al. (2024) WWTP data to datasets, and convert to units that LO uses
    '''

    i = WWTP_index

    # Add source ID, name, and lat/lon
    ds['ID'][i] = source_ID
    ds['name'][i] = source_name
    ds['lat'][i] = lat
    ds['lon'][i] = lon
    
    # Add flow and nutrients to dataset (and replace '.' with zeros)
    ds.flow[i,:] = was24_daily_df['FLOW_MGD'].replace('.', 0).astype(float)     * 0.0438126364 # convert from millions gallons per day to m3/s
    ds.NO3[i,:]  = was24_daily_df['NO2NO3N_MG_L'].replace('.', 0).astype(float) * 71.4 # convert to mmol/m3
    ds.NH4[i,:]  = was24_daily_df['NH4N_MG_L'].replace('.', 0).astype(float)    * 71.4 # convert to mmol/m3

    # Add other variables using climatology data from WWTPs in Mohamedali et al. (2020)
    for var_index,var in enumerate(['temp','TIC','Talk','DO']):
        ds[var][i,:] = moh20_clim_list[var_index] # note that variables are already in the correct units

    return ds

def moh20_riv_add_data(river_index, ds, source_ID, source_name, moh20_latlon_df, moh20_riv_df):
    '''
    Add Mohamedali et al. (2020) river data to datasets, and convert to units that LO uses
    '''

    i = river_index

    # Add source ID and name
    ds['ID'][i] = source_ID
    ds['name'][i] = source_name

    # Add source lat/lon
    # NOTE: we take the mean here because in SSM, large rivers are spread across two grid cell
    #       meaning that they have two lat/lon coordinates. We average to consolidate into one
    #       lat/lon coordinate. For rivers that are already in a single grid cell, the average of
    #       itself is itself.
    ds['lat'][i] = np.mean(moh20_latlon_df.loc[moh20_latlon_df['ID'] == source_ID, 'Lat'].values)
    ds['lon'][i] = np.mean(moh20_latlon_df.loc[moh20_latlon_df['ID'] == source_ID, 'Lon'].values)
    
    # Add physics and biology data to dataset
    ds.flow[i,:] = moh20_riv_df['Flow(m3/s)']
    ds.temp[i,:] = moh20_riv_df['Temp(C)']
    ds.NO3[i,:]  = moh20_riv_df['NO3+NO2(mg/L)'] * 71.4 # convert to mmol/m3
    ds.NH4[i,:]  = moh20_riv_df['NH4(mg/L)']     * 71.4 # convert to mmol/m3
    ds.TIC[i,:]  = moh20_riv_df['DIC(mmol/m3)']
    ds.Talk[i,:] = moh20_riv_df['Alk(mmol/m3)']
    ds.DO[i,:]   = moh20_riv_df['DO(mg/L)']      * 31.26 # convert to mmol/m3

    return ds

#################################################################################
#                               Get data location                               #
#################################################################################

# read dataset version (i.e. trapsP## listed in traps_data_ver.csv)
this_dir = Path(__file__).absolute().parent
with open(this_dir / 'traps_data_ver.csv','r') as f:
    for ver in f:
        trapsD = ver

#################################################################################
#                     Old WWTP data (Mohamedali et al., 2020)                   #
#################################################################################

# date range of data
moh20_start_date = datetime.date(1999, 1, 1)
moh20_end_date = datetime.date(2017, 8, 1) # one day after end of data time span
# get number of days in dataset time span
moh20_numdates = len(pd.date_range(moh20_start_date, moh20_end_date, freq='D')) - 1
# create date values
moh20_dates = pd.Series(pd.date_range(moh20_start_date, moh20_end_date, freq='D').values[0:-1])

# location of Mohamedali et al. (2020) data to process
moh20_pointsource_dir = Ldir['data'] / trapsD / 'mohamedali_etal2020' / 'point_sources'
moh20_pointsource_fns = os.listdir(moh20_pointsource_dir)

# Mohamedali et al. (2020) metadata with lat/lon coordinates
moh20_meta_fn = Ldir['data'] / trapsD / 'SSM_source_info.xlsx'
moh20_latlon_df = pd.read_excel(moh20_meta_fn,usecols='D,E,F,G,N,O')

# names of industrial facilities (which get omitted from data processing)
moh20_indu = ['BP Cherry Point',
            'Conoco Phillips',
            'Intalco',
            'Kimberly_Clark',
            'Nippon Paper',
            'Port Townsend Paper',
            'Shell Oil',
            'Tesoro',
            'US Oil & Refining',
            'West Rock']

# names of WWTPs that are in both moh20 and was24 (and get omitted from moh20 processing)
# open dataset of WWTP names for both datasets
same_wwtp_names_fn = Ldir['data'] / trapsD / 'wwtp_names.xlsx'
same_WWTPs_df = pd.read_excel(same_wwtp_names_fn)
# get list of WWTPs that are in both datasets (using naming convention in Mohamedali et al., 2020)
wwtp_in_both_moh20names = same_WWTPs_df['Mohamedali et al., 2020'][same_WWTPs_df['Wasielewski et al., 2024'].notna()].tolist()

# get number of WWTPs
# (total point sources in dataset minus number of industrial facilities
# and minus number of WWTPs that are also in Wasielewski et al., 2024)
NWWTP = np.shape(moh20_pointsource_fns)[0] - len(moh20_indu) - len(wwtp_in_both_moh20names)

# Start Dataset (with empty data)
source = np.arange(1,NWWTP+1)
moh20_wwtp_ds = start_ds(moh20_dates,source,moh20_numdates,NWWTP,'WWTPs','moh20')

# Add dataset metadata
moh20_wwtp_ds = add_metadata(moh20_wwtp_ds)

print('\nLooping through Mohamedali et al. (2020) WWTP files...')
# initialize WWTP_index for keeping track of WWTPs count
WWTP_index = 0

# Loop through all WWTPs and add data to dataset
for i,fn in enumerate(moh20_pointsource_fns):

    # get ID and source name
    source_ID = int(fn.split('_', 1)[0])
    source_name_xlsx = fn.split('_', 1)[1]
    source_name = source_name_xlsx.split('.', 1)[0]

    # do not add industrial facilities
    if source_name in moh20_indu:
        continue
    # do not add WWTPs that are also in Wasielewski et al. (2024)
    if source_name in wwtp_in_both_moh20names:
        continue

    # print counter
    print('    {}/{}: {}'.format(WWTP_index+1,NWWTP,source_name))

    # load data as a dataframe
    moh20_wwtp_fp = str(moh20_pointsource_dir)  + '/' + fn
    moh20_wwtp_monthly_df = pd.read_excel(moh20_wwtp_fp, skiprows=[0]) 
    
    # rename columns so that they are standardized
    # I have previously verified that Ecology's .xlsx files all have the same parameters
    moh20_wwtp_monthly_df = moh20_wwtp_monthly_df.set_axis(['Date', 'Year', 'Month', 'Day',
                            'Hour', 'Minute', 'Bin1', 'Flow(m3/s)',
                            'Temp(C)','Salt(ppt)','NH4(mg/L)',
                            'NO3+NO2(mg/L)', 'PO4(mg/L)', 'DO(mg/L)',
                            'pH', 'DON(mg/L)', 'PON(mg/L)', 'DOP(mg/L)',
                            'POP(mg/L)', 'POCS(mg/L)', 'POCF(mg/L)',
                            'POCR(mg/L)', 'DOCS(mg/L)', 'DOCF(mg/L)',
                            'Diatoms', 'Dinoflag', 'Chl', 'DIC(mmol/m3)',
                            'Alk(mmol/m3)'], axis=1)#, inplace=False)
    
    # point source data is monthly. Convert to daily
    moh20_wwtp_daily_df = moh20_monthly2daily(moh20_wwtp_monthly_df, moh20_start_date, moh20_end_date)

    # Add Mohamedali et al. (2020) WWTP data to dataset
    moh20_wwtp_ds = moh20_wwtp_add_data(WWTP_index, moh20_wwtp_ds, source_ID, source_name, moh20_latlon_df, moh20_wwtp_daily_df)

    # advance WWTP counter
    WWTP_index += 1

# save dataset as .nc file in LO_data
out_fn = '../../../LO_data/' + trapsD + '/processed_data/wwtp_data_mohamedali_etal_2020.nc'
moh20_wwtp_ds.to_netcdf(out_fn)
moh20_wwtp_ds.close()
print('Mohamedali et al. (2020) WWTPs complete --------------------------------------------\n')

#################################################################################
#                   New WWTP data (Wasielewski et al., 2020)                    #
#################################################################################

# date range of data
was24_start_date = datetime.date(2005, 1, 1)
was24_end_date = datetime.date(2021, 1, 1) # one day after end of data time span
# get number of days in dataset time span
was24_numdates = len(pd.date_range(was24_start_date, was24_end_date, freq='D')) - 1

# create date values
was24_dates = pd.Series(pd.date_range(was24_start_date, was24_end_date, freq='D').values[0:-1])

# location of point source data to process
was24_pointsource_dir = Ldir['data'] / trapsD / 'wasielewski_etal2024' / 'point_sources'
# Wasielewski et al. (2024)
was24_pointsource_meta = was24_pointsource_dir / 'fac_attributes.csv'
was24_pointsource_loads = was24_pointsource_dir / 'nutrient_loads.csv'
# get point source data
was24_meta_df = pd.read_csv(was24_pointsource_meta)
was24_loads_df = pd.read_csv(was24_pointsource_loads)

# get a list of all of the source IDs
all_was24_sources = was24_meta_df['FAC_ID']

# get list of wwtp names which are in both datasets (using naming convention from Wasielewski et al., 2024) 
wwtp_in_both_was24names = same_WWTPs_df['Wasielewski et al., 2024'][same_WWTPs_df['Wasielewski et al., 2024'].notna()].tolist()
# get associated IDs
wwtp_in_both_was24IDs = was24_meta_df[was24_meta_df['FAC_NAME'].isin(wwtp_in_both_was24names)]['FAC_ID'].values

# get climatology for temp, TIC, Talk, DO from Mohamedali et al. (2020)
moh20_clim_list = get_moh20_wwtp_clim(moh20_wwtp_ds)

# initialize a list of source IDs to process
was24_wwtp_IDs = []

# loop through all of the sources and keep only WWTPs that are also listed in Mohamedali et al. (2020)
for ID in all_was24_sources:
        # remove anything that is not a WWTP (ie omit fish hatcheries and industrial facilities)
        # get facility type
        fac_type = was24_meta_df.loc[was24_meta_df['FAC_ID'] == ID, 'FAC_TYPE'].values[0]
        if 'sic_4952' in fac_type:
            # only include WWTPs if they are also listed in Mohamedali et al. (2020)
            if ID in wwtp_in_both_was24IDs:
                was24_wwtp_IDs = was24_wwtp_IDs + [ID]  

# remove WWTPs that stopped operation prior to 2012 (earliest LiveOcean can run)
was24_wwtp_IDs.remove('WA0020567-001') # OAK HARBOR STP which ends mid-way through 2010
was24_wwtp_IDs.remove('WA0020893-thru2012') # LAKE STEVENS SEWER DISTRICT which was moved at the start of 2013

# Start Dataset (with empty data)
NWWTP = len(was24_wwtp_IDs)
source = np.arange(1,NWWTP+1)
was24_wwtp_ds = start_ds(was24_dates,source,was24_numdates,NWWTP,'WWTPs','was24')

# Add dataset metadata
was24_wwtp_ds = add_metadata(was24_wwtp_ds)

print('\nLooping through Wasielewski et al. (2024) WWTPs...')
# initialize WWTP_index for keeping track of WWTPs count
WWTP_index = 0

# Loop through all WWTPs and add data to dataset
for i,FAC_ID in enumerate(was24_wwtp_IDs):

    # source name
    source_name = was24_meta_df.loc[was24_meta_df['FAC_ID'] == FAC_ID, 'FAC_NAME'].values[0]

    # print counter
    print('    {}/{}: {}'.format(i+1,NWWTP,source_name))

    # use lat/lon positions from Mohamedali et al. (2020) datasset
    # get name in Mohamedali et al. (2020) dataset
    if source_name == 'Everett Water Pollution Control Facility':
            if FAC_ID == 'WA0024490_Gardner':
                    moh20name = 'OF100'
                    source_name = 'Gardner - Everett Water Pollution Control Facility'
            elif FAC_ID == 'WA0024490_Snohomish':
                    moh20name = 'Everett Snohomish'
                    source_name = 'Snohomish - Everett Water Pollution Control Facility'
    elif source_name == 'OAK HARBOR STP':
            if FAC_ID == 'WA0020567-002':
                    moh20name = 'Oak Harbor Lagoon'
    else:
        moh20name = same_WWTPs_df.loc[same_WWTPs_df['Wasielewski et al., 2024']==source_name]['Mohamedali et al., 2020'].values[0]
    # get corresponding lat/lon coordinate
    lat = moh20_latlon_df.loc[moh20_latlon_df['Name'] == moh20name, 'Lat'].values[0]
    lon = moh20_latlon_df.loc[moh20_latlon_df['Name'] == moh20name, 'Lon'].values[0]

    # crop loading data to just the current point source
    was24_monthly_df = was24_loads_df[was24_loads_df['FAC_ID']==FAC_ID]

    # special handling case to combine both Lake Stevens Sewer District WWTPs into one time series
    if FAC_ID == 'WA0020893':
         # get old Lake Stevens WWTP data
         old_LakeStevens_df = was24_loads_df[was24_loads_df['FAC_ID']=='WA0020893-thru2012']
         # concatenate new and old WWTP datasets
         was24_monthly_df = pd.concat([old_LakeStevens_df,was24_monthly_df])

    # test if any WWTPs have missing data,
    # and if there are any, plot what data points are missing
    if len(was24_monthly_df['FAC_ID']) < 192:
        # create figure to show what dates are missing
        fig, ax = plt.subplots(1,1,figsize = (5,5))
        ax.grid(True,color='silver',linewidth=1,linestyle='--',axis='both')
        # add data
        ax.scatter(was24_monthly_df['YEAR'],was24_monthly_df['MONTH'],zorder=3)
        ax.set_ylabel('Month')
        ax.set_xlabel('Year')
        ax.set_title('{}: {}'.format(FAC_ID,source_name))
        ax.set_ylim([0,13])
        ax.set_xlim([2004,2021])
        plt.show()

    # convert monthly to daily data
    was24_daily_df = was24_monthly2daily(was24_monthly_df, was24_start_date, was24_end_date)

    # add data to dataset
    was24_wwtp_ds = was24_wwtp_add_data(WWTP_index, was24_wwtp_ds, FAC_ID, source_name,
                                        lat, lon, was24_daily_df, moh20_clim_list)
        
    # increment WWTP counter
    WWTP_index += 1


# save dataset as .nc file in LO_data
out_fn = '../../../LO_data/' + trapsD + '/processed_data/wwtp_data_wasielewski_etal_2024.nc'
was24_wwtp_ds.to_netcdf(out_fn)
was24_wwtp_ds.close()
print('Wasielewski et al. (2024) WWTPs complete --------------------------------------------\n')

#################################################################################
#                     River data (Mohamedali et al., 2020)                      #
#################################################################################

# location of Mohamedali et al. (2020) data to process
moh20_river_dir = Ldir['data'] / trapsD / 'mohamedali_etal2020' / 'nonpoint_sources'
moh20_river_fns = os.listdir(moh20_river_dir)
NTRIV = np.shape(moh20_river_fns)[0]

# Start Dataset (with empty data)
source = np.arange(1,NTRIV+1)
moh20_riv_ds = start_ds(moh20_dates,source,moh20_numdates,NTRIV,'Rivers','moh20')

# Add dataset metadata
moh20_riv_ds = add_metadata(moh20_riv_ds)

print('\nLooping through Mohamedali et al. (2020) river files...')
# initialize river_index for keeping track of river count
river_index = 0

# Loop through all rivers and add data to dataset
for i,fn in enumerate(moh20_river_fns):

    # get ID and source name
    source_ID = int(fn.split('_', 1)[0])
    source_name_xlsx = fn.split('_', 1)[1]
    source_name = source_name_xlsx.split('.', 1)[0]
    print('    {}/{}: {}'.format(i+1,NTRIV,source_name))

    # load data as a dataframe
    moh20_riv_fp = str(moh20_river_dir)  + '/' + fn
    moh20_riv_df = pd.read_excel(moh20_riv_fp, skiprows=[0]) 
    
    # rename columns so that they are standardized
    # I have previously verified that Ecology's .xlsx files all have the same parameters
    # note that river data are already daily values (not monthly)
    moh20_riv_df = moh20_riv_df.set_axis(['Date', 'Year', 'Month', 'Day',
                            'Hour', 'Minute', 'Bin1', 'Flow(m3/s)',
                            'Temp(C)','Salt(ppt)','NH4(mg/L)',
                            'NO3+NO2(mg/L)', 'PO4(mg/L)', 'DO(mg/L)',
                            'pH', 'DON(mg/L)', 'PON(mg/L)', 'DOP(mg/L)',
                            'POP(mg/L)', 'POCS(mg/L)', 'POCF(mg/L)',
                            'POCR(mg/L)', 'DOCS(mg/L)', 'DOCF(mg/L)',
                            'Diatoms', 'Dinoflag', 'Chl', 'DIC(mmol/m3)',
                            'Alk(mmol/m3)'], axis=1)#, inplace=False)

    # Add Ecology data to dataset
    moh20_riv_ds = moh20_riv_add_data(river_index, moh20_riv_ds, source_ID, source_name, moh20_latlon_df, moh20_riv_df)

    # increment index
    river_index += 1

# save dataset as .nc file in LO_data
out_fn = '../../../LO_data/' + trapsD + '/processed_data/river_data_mohamedali_etal_2020.nc'
moh20_riv_ds.to_netcdf(out_fn)
moh20_riv_ds.close()
print('Mohamedali et al. (2020) rivers complete --------------------------------------------\n')

# print(moh20_wwtp_ds)
# print(was24_wwtp_ds)
# print(moh20_riv_ds)