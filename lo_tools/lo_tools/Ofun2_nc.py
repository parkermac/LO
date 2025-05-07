"""
Functions for writing saved fields to ROMS NetCDF.

Rewritten to use xarray instead of netCDF4

RESULT: this works great with ROMS as long as we:
(i) transform the arrays that have nan's
(ii) write them to NetCDF with a _FillValue
(iii) and the compression makes them MUCH smaller!
"""

import os
import xarray as xr
import numpy as np
from lo_tools import Lfun, zrfun
import pickle

def make_ini_file(in_fn, out_fn):
    """
    Create the ini file from the first time of the clm file.
    """
    ds0 = xr.open_dataset(in_fn, decode_times=False)
    ot_vec = ds0.zeta_time.values
    ds = xr.Dataset()
    ds['ocean_time'] = (('ocean_time',), [ot_vec[0]])
    ds['ocean_time'].attrs['units'] = Lfun.roms_time_units
    for vn in ds0.data_vars:
        vinfo = zrfun.get_varinfo(vn)
        ndims = len(ds0[vn].dims)
        if ndims == 3:
            ds[vn] = (('ocean_time',) + ds0[vn].dims[1:], ds0[vn].values[[0],:,:])
            # Note: we use [0] instead of 0 to retain the singleton dimension
        elif ndims == 4:
            ds[vn] = (('ocean_time',) + ds0[vn].dims[1:], ds0[vn].values[[0],:,:,:])
        ds[vn].attrs['units'] = vinfo['units']
        ds[vn].attrs['long_name'] = vinfo['long_name']
    ds0.close()
    Enc_dict = {vn:zrfun.enc_dict for vn in ds.data_vars}
    # and save to NetCDF
    ds.to_netcdf(out_fn, encoding=Enc_dict)
    ds.close()
    
def make_bry_file(in_fn, out_fn):
    """
    Create the bry file from the edges of the clm file.
    """
    ds0 = xr.open_dataset(in_fn, decode_times=False)
    ot_vec = ds0.salt_time.values
    ds = xr.Dataset()
    for vn in ds0.data_vars:
        dm = ds0[vn].dims
        ndims = len(dm)
        for D in ['north', 'south', 'east', 'west']:
            # rename variable
            # We have to deal with some inconsistency in the naming of bio boundary variables.
            bvn_dict = {'phytoplankton': 'phyt', 'zooplankton': 'zoop', 'chlorophyll': 'chlo',
                    'LdetritusN': 'LDeN', 'SdetritusN': 'SDeN',
                    'LdetritusC': 'LdeC', 'SdetritusC': 'SdeC'}
            if vn in bvn_dict.keys():
                vn_new = bvn_dict[vn]
            else:
                vn_new = vn
            Vn = vn_new + '_' + D
            try:
                vinfo = zrfun.get_varinfo(Vn, vartype='climatology')
            except Exception as e:
                print(Vn)
                print(e)
            tname = vinfo['time']
            
            # create time coordinate
            ds[tname] = ((tname,), ot_vec)
            ds[tname].attrs['units'] = Lfun.roms_time_units
            
            # trim dimensions
            if D in ['east','west']:
                Dm = tuple(item for item in dm if 'xi_' not in item)
            elif D in ['north','south']:
                Dm = tuple(item for item in dm if (('eta_' not in item) or ('zeta' in item)))
            
            # replace time dimension
            Dm = tuple(tname if '_time' in item else item for item in Dm) 

            # write boundary arrays
            if ndims == 3:
                if D == 'north':
                    ds[Vn] = (Dm, ds0[vn].values[:,-1,:])
                elif D == 'south':
                    ds[Vn] = (Dm, ds0[vn].values[:,0,:])
                elif D == 'east':
                    ds[Vn] = (Dm, ds0[vn].values[:,:,-1])
                elif D == 'west':
                    ds[Vn] = (Dm, ds0[vn].values[:,:,0])
            elif ndims == 4:
                if D == 'north':
                    ds[Vn] = (Dm, ds0[vn].values[:,:,-1,:])
                elif D == 'south':
                    ds[Vn] = (Dm, ds0[vn].values[:,:,0,:])
                elif D == 'east':
                    ds[Vn] = (Dm, ds0[vn].values[:,:,:,-1])
                elif D == 'west':
                    ds[Vn] = (Dm, ds0[vn].values[:,:,:,0])
                    
            # add attributes
            ds[Vn].attrs['units'] = vinfo['units']
            ds[Vn].attrs['long_name'] = vinfo['long_name']
    
    ds0.close()
    Enc_dict = {vn:zrfun.enc_dict for vn in ds.data_vars}
    # and save to NetCDF
    ds.to_netcdf(out_fn, encoding=Enc_dict)
    ds.close()

