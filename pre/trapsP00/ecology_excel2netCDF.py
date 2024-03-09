"""
This script compiles all of Ecology's excel
loading data into two netCDF files:
one for point sources
one for rivers

Looks at data stored in LO_data/trapsD##
(To change the Ecology data version, modify traps_data_ver.csv)

In theory, this script only needs to be run once.
Then, the netCDF files can be referenced to generate climatologies

Takes about 5-8 minutes to run on my local machine.

To run from ipython:
run ecology_excel2netCDF.py
"""

#################################################################################
#                              Import packages                                  #
#################################################################################
from lo_tools import Lfun
Ldir = Lfun.Lstart()

import pandas as pd
import numpy as np
import os
import datetime
import xarray as xr
from pathlib import Path

#################################################################################
#                              Helper functions                                 #
#################################################################################

def monthly2daily(df):
    '''
    turn a monthly dataframe into daily data (for the lenght of Ecology's timeseries)
    '''
    # duplicate last row
    double_lr_df = pd.concat([df, df.iloc[-1:]], ignore_index=True)
    # picking arbitrary year to fill with daily data (but includes a leap year)
    start_date = datetime.date(1999, 1, 1)
    end_date = datetime.date(2017, 8, 1)
    dates = pd.date_range(start_date, end_date, freq='MS')
    # Replace month column with things that look like dates
    double_lr_df['Month'] = dates
    double_lr_df = double_lr_df.set_index('Month')
    # Change monthly to daily
    double_lr_daily_df = double_lr_df.resample('D').ffill()
    # delete last row (1/1 on the next year)
    daily_df = double_lr_daily_df[:-1]
    # make index start from 1 and go to 366
    daily_df.reset_index(inplace=True)
    return daily_df

def start_ds(date,source,numdates,Nsources,source_type):
    '''
    Initialize dataset to story Ecology's loading data
    '''
    ds = xr.Dataset(data_vars=dict(ID=(['source'], np.ones((Nsources,), dtype=int)),
        lon=(['source'], np.ones((Nsources,))),
        lat=(['source'], np.ones((Nsources,))),
        name=(['source'], ['placeholder placeholder placeholder']*Nsources),
        flow=(['source', 'date'], np.zeros((Nsources, numdates))),
        temp=(['source', 'date'], np.zeros((Nsources, numdates))),
        NO3=(['source', 'date'], np.zeros((Nsources, numdates))),
        NH4=(['source', 'date'], np.zeros((Nsources, numdates))),
        TIC=(['source', 'date'], np.zeros((Nsources, numdates))),
        Talk=(['source', 'date'], np.zeros((Nsources, numdates))),
        DO=(['source', 'date'], np.zeros((Nsources, numdates))),),
    coords=dict(source=source, date=date,),
    attrs=dict(description='Ecology data for '+source_type),)
    
    return ds

def add_metadata(ds):
    '''
    Create metadata for dataset of Ecology loading data
    '''
    ds['ID'].attrs['long_name'] = 'source ID used in Salish Sea Model'
    ds['lon'].attrs['long_name'] = 'point source longitude'
    ds['lat'].attrs['long_name'] = 'point source latitude'
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

def add_data(ds, source_ID, source_name, latlon_df, ecologydata_df):
    '''
    Add Ecology's data to datasets, and convert to units that LO uses
    '''
    # Add source ID and name
    ds['ID'][i] = source_ID
    ds['name'][i] = source_name

    # Add source lat/lon
    # NOTE: we take the mean here because in SSM, large rivers are spread across two grid cell
    #       meaning that they have two lat/lon coordinates. We average to consolidate into one
    #       lat/lon coordinate. For rivers that are already in a single grid cell, the average of
    #       itself is itself.
    ds['lat'][i] = np.mean(latlon_df.loc[latlon_df['ID'] == source_ID, 'Lat'].values)
    ds['lon'][i] = np.mean(latlon_df.loc[latlon_df['ID'] == source_ID, 'Lon'].values)
    
    # Add physics and biology data to dataset
    ds.flow[i,:] = ecologydata_df['Flow(m3/s)']
    ds.temp[i,:] = ecologydata_df['Temp(C)']
    ds.NO3[i,:]  = ecologydata_df['NO3+NO2(mg/L)'] * 71.4 # convert to mmol/m3
    ds.NH4[i,:]  = ecologydata_df['NH4(mg/L)']     * 71.4 # convert to mmol/m3
    ds.TIC[i,:]  = ecologydata_df['DIC(mmol/m3)']
    ds.Talk[i,:] = ecologydata_df['Alk(mmol/m3)']
    ds.DO[i,:]   = ecologydata_df['DO(mg/L)']      * 31.26 # convert to mmol/m3

    return ds

#################################################################################
#                              Get path to data                                 #
#################################################################################

# read Ecology data version (i.e. trapsP## listed in traps_data_ver.csv)
this_dir = Path(__file__).absolute().parent
with open(this_dir / 'traps_data_ver.csv','r') as f:
    for ver in f:
        trapsD = ver

# location of historical data to process
wwtp_dir = Ldir['data'] / trapsD / 'point_sources'
wwtp_fns = os.listdir(wwtp_dir)
NWWTP = np.shape(wwtp_fns)[0]

# location of historical data to process
riv_dir = Ldir['data'] / trapsD / 'nonpoint_sources'
riv_fns = os.listdir(riv_dir)
NTRIV = np.shape(riv_fns)[0]

# SSM metadata with lat/lon coordinates
trapsll_fn = Ldir['data'] / trapsD / 'SSM_source_info.xlsx'
latlon_df = pd.read_excel(trapsll_fn,usecols='D,E,F,G,N,O')

# Start with one point source to get date information
# note that need to use river data because river data is daily, wwtp is only monthly
riv_fp = str(riv_dir)  + '/' + riv_fns[0]
df_example = pd.read_excel(riv_fp, skiprows=[0])
numdates = len(df_example['Date'])

#################################################################################
#                         Create point source dataset                           #
#################################################################################

# Start Dataset (with empty data)
date = df_example['Date']
source = np.arange(1,NWWTP+1)
pointsource_ds = start_ds(date,source,numdates,NWWTP,'point sources')

# Add dataset metadata
pointsource_ds = add_metadata(pointsource_ds)

print('Looping through point source files...')

# Loop through all WWTPs and add data to dataset
for i,fn in enumerate(wwtp_fns):

    # get ID and source name
    source_ID = int(fn.split('_', 1)[0])
    source_name_xlsx = fn.split('_', 1)[1]
    source_name = source_name_xlsx.split('.', 1)[0]
    print('{}/{}: {}'.format(i+1,NWWTP,source_name))

    # load data as a dataframe
    wwtp_fp = str(wwtp_dir)  + '/' + fn
    wwtp_monthly_df = pd.read_excel(wwtp_fp, skiprows=[0]) 
    
    # rename columns so that they are standardized
    # I have previously verified that Ecology's .xlsx files all have the same parameters
    wwtp_monthly_df = wwtp_monthly_df.set_axis(['Date', 'Year', 'Month', 'Day',
                            'Hour', 'Minute', 'Bin1', 'Flow(m3/s)',
                            'Temp(C)','Salt(ppt)','NH4(mg/L)',
                            'NO3+NO2(mg/L)', 'PO4(mg/L)', 'DO(mg/L)',
                            'pH', 'DON(mg/L)', 'PON(mg/L)', 'DOP(mg/L)',
                            'POP(mg/L)', 'POCS(mg/L)', 'POCF(mg/L)',
                            'POCR(mg/L)', 'DOCS(mg/L)', 'DOCF(mg/L)',
                            'Diatoms', 'Dinoflag', 'Chl', 'DIC(mmol/m3)',
                            'Alk(mmol/m3)'], axis=1, inplace=False)
    
    # point source data is monthly. Convert to daily
    wwtp_df = monthly2daily(wwtp_monthly_df)

    # Add Ecology data to dataset
    pointsource_ds = add_data(pointsource_ds, source_ID, source_name, latlon_df, wwtp_df)

# save dataset as .nc file in LO_data
out_fn = '../../../LO_data/' + trapsD + '/all_point_source_data.nc'
pointsource_ds.to_netcdf(out_fn)
pointsource_ds.close()
print('Point sources complete --------------------------------------------\n')

#################################################################################
#                       Create nonpoint source dataset                          #
#################################################################################

# Start Dataset (with empty data)
date = df_example['Date']
source = np.arange(1,NTRIV+1)
nonpointsource_ds = start_ds(date,source,numdates,NTRIV,'nonpoint sources')

# Add dataset metadata
nonpointsource_ds = add_metadata(nonpointsource_ds)

print('Looping through nonpoint source files...')

# Loop through all rivers and add data to dataset
for i,fn in enumerate(riv_fns):

    # get ID and source name
    source_ID = int(fn.split('_', 1)[0])
    source_name_xlsx = fn.split('_', 1)[1]
    source_name = source_name_xlsx.split('.', 1)[0]
    print('{}/{}: {}'.format(i+1,NTRIV,source_name))

    # load data as a dataframe
    riv_fp = str(riv_dir)  + '/' + fn
    riv_df = pd.read_excel(riv_fp, skiprows=[0]) 
    
    # rename columns so that they are standardized
    # I have previously verified that Ecology's .xlsx files all have the same parameters
    riv_df = riv_df.set_axis(['Date', 'Year', 'Month', 'Day',
                            'Hour', 'Minute', 'Bin1', 'Flow(m3/s)',
                            'Temp(C)','Salt(ppt)','NH4(mg/L)',
                            'NO3+NO2(mg/L)', 'PO4(mg/L)', 'DO(mg/L)',
                            'pH', 'DON(mg/L)', 'PON(mg/L)', 'DOP(mg/L)',
                            'POP(mg/L)', 'POCS(mg/L)', 'POCF(mg/L)',
                            'POCR(mg/L)', 'DOCS(mg/L)', 'DOCF(mg/L)',
                            'Diatoms', 'Dinoflag', 'Chl', 'DIC(mmol/m3)',
                            'Alk(mmol/m3)'], axis=1, inplace=False)

    # Add Ecology data to dataset
    nonpointsource_ds = add_data(nonpointsource_ds, source_ID, source_name, latlon_df, riv_df)

# save dataset as .nc file in LO_data
out_fn = '../../../LO_data/' + trapsD + '/all_nonpoint_source_data.nc'
nonpointsource_ds.to_netcdf(out_fn)
nonpointsource_ds.close()
print('Nonpoint sources complete ---------------------------------------')

# print(pointsource_ds)
# print(nonpointsource_ds)