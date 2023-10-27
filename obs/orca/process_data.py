"""
Code to process the ORCA mooring data.

This is just a light reformatting of the extensive pre-processing Erin Broatch did
to get the ORCA mooring data into clean daily values in xarray Datasets.

This also is the first instance of formatting MOORING data in the LO standard
for LO_output/obs.

Performance: Takes only seconds to run

"""

import pandas as pd
import xarray as xr
import numpy as np
import gsw
from time import time

from lo_tools import Lfun
Ldir = Lfun.Lstart()

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-test','--testing', default=False, type=Lfun.boolean_string)
args = parser.parse_args()
testing = args.testing

# input location
source = 'orca'
otype = 'moor' # introducing a new "otype" beyond ctd and bottle
in_dir = Ldir['data'] / 'obs' / 'ORCA' / 'orca_profiles' / 'datasets'

# output location
out_dir = Ldir['LOo'] / 'obs' / source / otype
if not testing:
    Lfun.make_dir(out_dir, clean=True)

sn_name_dict = {
    'CI':'Carr Inlet',
    'PW':'Point Wells',
    'NB':'Hansville',
    'DB':'Dabob Bay',
    'HP':'Hoodsport',
    'TW':'Twanoh'
}

if testing:
    sn_list = ['CI']
else:
    sn_list = list(sn_name_dict.keys())
    
sn_loc_dict = {
    'CI': [-122.7300, 47.2800],
    'PW': [-122.3972, 47.7612],
    'NB': [-122.6270, 47.9073],
    'DB': [-122.8029, 47.8034],
    'HP': [-123.1126, 47.4218],
    'TW': [-123.0083, 47.3750]
}

for sn in sn_list:
    print(sn)
    in_fn = in_dir / (sn + '_ds_daily.nc')
    out_fn = out_dir / (sn + '_daily.nc')
    ds0 = xr.open_dataset(in_fn)
    """
    Original form:
    Dimensions:    (time: 6181, z: 29)
    Coordinates:
      * time       (time) datetime64[ns] 2005-01-20 2005-01-21 ... 2021-12-22
        pressure   (z) float64 ...
        depth      (z) float64 ...
    Dimensions without coordinates: z
    Data variables:
        sal        (z, time) float64 ...
        temp       (z, time) float64 ...
        oxy        (z, time) float64 ...
        nitrate    (z, time) float64 ...
        fluor      (z, time) float64 ...
        par        (z, time) float64 ...
        depth_obs  (z, time) float64 ...
    
    And the units of the data_vars, from doing:
    for vn in ds0.data_vars:
        print('%s (%s)' % (vn,ds0[vn].units))
    sal (psu)
    temp (degC)
    oxy (mg/L)
    nitrate (umol)
    fluor (mg/m^3)
    par (uEinstein/m^2s)
    depth_obs (m)
    """
    # convert to LO standard units
    z = -ds0.depth.values # packed shallow to deep
    zz = z[::-1] # packed deep to shallow
    NZ = len(zz)
    t = ds0.time.values
    tt = pd.DatetimeIndex(t)
    NT = len(tt)
    Z = zz.reshape((1,NZ)) * np.ones((NT,1))
    # salt and temp
    SP = np.transpose(ds0.sal.values[::-1,:])
    IT = np.transpose(ds0.temp.values[::-1,:])
    lon = sn_loc_dict[sn][0]
    lat = sn_loc_dict[sn][1]
    P = gsw.p_from_z(Z, lat)
    # - do the conversions
    SA = gsw.SA_from_SP(SP, P, lon, lat)
    CT = gsw.CT_from_t(SA, IT, P)
    # oxygen converted mg/L to uM
    DO = (1000/32) * np.transpose(ds0.oxy.values[::-1,:])
    # nitrate is already uM
    NO3 = np.transpose(ds0.nitrate.values[::-1,:])
    # no conversion of fluor or par (is the data any good?)
    FLUOR = np.transpose(ds0.fluor.values[::-1,:])
    PAR = np.transpose(ds0.par.values[::-1,:])
    # compute sigma0
    # potential density relative to 0 dbar, minus 1000 kg m-3
    SIG0 = gsw.sigma0(SA,CT)
    
    # initialize new Dataset and fill
    coords = {'time':(('time'),tt),'z':(('z'),zz)}
    ds = xr.Dataset(coords=coords, attrs={'Station Name':sn_name_dict[sn],'lon':lon,'lat':lat})
    ds['SA'] = xr.DataArray(SA, dims=('time','z'),
        attrs={'units':'g kg-1', 'long_name':'Absolute Salinity'})
    ds['CT'] = xr.DataArray(CT, dims=('time','z'),
        attrs={'units':'degC', 'long_name':'Conservative Temperature'})
    ds['DO (uM)'] = xr.DataArray(DO, dims=('time','z'),
        attrs={'units':'uM', 'long_name':'Dissolved Oxygen'})
    ds['NO3 (uM)'] = xr.DataArray(NO3, dims=('time','z'),
        attrs={'units':'uM', 'long_name':'Nitrate'})
    ds['FLUOR'] = xr.DataArray(FLUOR, dims=('time','z'),
        attrs={'units':'mg m-3', 'long_name':'Fluoresence'})
    ds['PAR'] = xr.DataArray(PAR, dims=('time','z'),
        attrs={'units':'uEinstein m-2 s-1',
            'long_name':'Photosynthetically Active Radiation'})
    # add sigma0
    ds['SIG0'] = xr.DataArray(SIG0, dims=('time','z'),
        attrs={'units':'kg m-3', 'long_name':'Sigma0'})
    
    if not testing:
        ds.to_netcdf(out_fn)
