"""
Code to process the ocnms mooring data from 2011 - 2023 

This imports pre-processed hourly matlab data and places to LO format 
in xarray Datasets, and then saves as a netcdf.

Takes a ~10 seconds to run on Kate's mac

todos: (1) when get new 2024 data, will need to update & add so can place
place straight to LO/obs/... the naming of files from OCNMS changes 
+ sometimes how the files are saved too. This is a work in progress, 
but someone can run this code after minimal matlab processing
(2) might want to add flags so that can tell when removed data that was flagged suspect. 
eg the late deployment O2 data from 2021 at MB042 

"""

from scipy.io import loadmat
import pandas as pd
import xarray as xr
import numpy as np
import gsw
import datetime as dt
from time import time

from lo_tools import Lfun
Ldir = Lfun.Lstart()

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-test','--testing', default=False, type=Lfun.boolean_string)
args = parser.parse_args()
testing = args.testing

# input location
source = 'ocnms'
otype = 'moor' 
in_dir = Ldir['data'] / 'obs' / 'ocnms' / 'ocnms_mooring' / 'datasets'

# output location
out_dir = Ldir['LOo'] / 'obs' / source / otype

sn_name_dict = {
    'MB015':'Makah Bay 15m',
    'MB042':'Makah Bay 42m',
    'CA015':'Cape Alava 15m',
    'CA042':'Cape Alava 42m',
    'TH015':'Teahwhit Head 15m',
    'TH042':'Teahwhit Head 42m',
    'KL015':'Kalaloch 15m',
    'KL027':'Kalaloch 27m',
    'CE015':'Cape Elizabeth 15m',
    'CE042':'Cape Elizabeth 42m'
}
    
if testing:
    sn_list = ['CE042']
else:    
    sn_list = list(sn_name_dict.keys())
    
for sn in sn_list:
    print(sn)
    in_fn = in_dir / (sn + '_2011_2023_hourly.mat')
    out_fn = out_dir / (sn + '_2011_2023_hourly.nc')
    
    mat_dict = loadmat(in_fn,squeeze_me=True)
    
    # convert matlab struct data to LO standard units 
    
    # import times import matlab datenums: we round to the nearest second and 
    # 719529 is the datenum value of the Unix epoch start (1970-01-01)
    # i.e. in matlab: datenum(datetime(1970,01,01)) = 719529
    matlab_datetimes = mat_dict['timestamp_UTC']                               # an array of datenums from matlab 
    
    # if save with the timezone when we set up the dataset it saves 
    # the time coordinate as an object and not a datetime index. Changed to remove 
    # utc = True, and then everything runs fine when saving xr.Dataset
    # tt = pd.to_datetime(matlab_datetimes-719529,unit='d',utc=True).round('s')  # datetime index
      
    tt = pd.to_datetime(matlab_datetimes-719529,unit='d').round('s')  # datetime index
    NT = len(tt)
    
    z = mat_dict['Z_est']  # set up shallow to deep; depths are estimated 
    NZ = len(z)
    Z = z.reshape((1,NZ)) * np.ones((NT,1))
    
    # salt and temp
    SP = mat_dict['SAL'] 
    IT = mat_dict['IT']
    lon = mat_dict['longitude']
    lat = mat_dict['latitude']
    P = gsw.p_from_z(Z,lat)
    SA = gsw.SA_from_SP(SP, P, lon, lat)
    CT = gsw.CT_from_t(SA, IT, P)
    
    # oxygen converted from mg/L to uM 
    DO = (1000/32) * mat_dict['OXY']
    
    # compute sigma0
    # potential density relative to 0 dbar, minus 1000 kg m-3
    SIG0 = gsw.sigma0(SA,CT)
    
    #initialize new dataset and fill
    coords = {'time':('time',tt),'z':('z',z)}
    ds = xr.Dataset(coords=coords, attrs={'Station Name':sn_name_dict[sn],'lon':lon,'lat':lat})
    
    ds['SA'] = xr.DataArray(SA, dims=('time','z'),
        attrs={'units':'g kg-1', 'long_name':'Absolute Salinity'})
    ds['SP'] = xr.DataArray(SA, dims=('time','z'),
        attrs={'units':' ', 'long_name':'Practical Salinity'})
    ds['IT'] = xr.DataArray(IT, dims=('time','z'),
        attrs={'units':'degC', 'long_name':'Insitu Temperature'})
    ds['CT'] = xr.DataArray(CT, dims=('time','z'),
        attrs={'units':'degC', 'long_name':'Conservative Temperature'})
    ds['DO (uM)'] = xr.DataArray(DO, dims=('time','z'),
        attrs={'units':'uM', 'long_name':'Dissolved Oxygen'})
    ds['SIG0'] = xr.DataArray(SIG0, dims=('time','z'),
        attrs={'units':'kg m-3', 'long_name':'Sigma0'})
        
    if not testing:
        ds.to_netcdf(out_fn, unlimited_dims='time')
    

    
    
    
    
    
    
