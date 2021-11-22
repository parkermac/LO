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
from lo_tools import Lfun
import pickle

enc_dict = {'zlib':True, 'complevel':1, '_FillValue':1e20}
# Using compression (zlib=True, complevel=1) results in files that are just 2% of the
# uncompressed files (for hc0, which has a lot of nan's).
# Using complevel=9 makes the files half as big as complevel=1, but takes much about 10x longer.

def make_clm_file(temp_out_fn_list, NT, out_fn):
    # associate variables with dimenstions
    vn_dict = {
            'zeta': ('zeta_time', 'eta_rho', 'xi_rho'),
            'ubar': ('v2d_time', 'eta_u', 'xi_u'),
            'vbar': ('v2d_time', 'eta_v', 'xi_v'),
            'salt': ('salt_time', 's_rho', 'eta_rho', 'xi_rho'),
            'temp': ('temp_time', 's_rho', 'eta_rho', 'xi_rho'),
            'u': ('v3d_time', 's_rho', 'eta_u', 'xi_u'),
            'v': ('v3d_time', 's_rho', 'eta_v', 'xi_v')
            }
    # assign attributes to variables
    attrs_dict = {
            'zeta': {'long_name': 'sea surface height climatology', 'units':'meter'},
            'ubar': {'long_name': 'vertically averaged u-momentum climatology', 'units':'meter second-1'},
            'vbar': {'long_name': 'vertically averaged v-momentum climatology', 'units':'meter second-1'},
            'salt': {'long_name': 'salinity climatology', 'units':'g kg-1'},
            'temp': {'long_name': 'potential temperature climatology', 'units':'Celsius'},
            'u': {'long_name': 'u-momentum component climatology', 'units':'meter second-1'},
            'v': {'long_name': 'v-momentum component climatology', 'units':'meter second-1'}
            }
    # write fields to the Dataset
    ds = xr.Dataset()
    ot_vec = np.nan * np.ones(NT)
    ii = 0
    for temp_out_fn in temp_out_fn_list:
        dd = pickle.load(open(temp_out_fn, 'rb'))
        if ii == 0:
            # initialize the variables in the Dataset
            for vn in dd.keys():
                if vn == 'ocean_time':
                    pass
                else:
                    ds[vn] = (vn_dict[vn] , np.nan * np.ones((NT,) + dd[vn].shape))
                    ds[vn].attrs = attrs_dict[vn]
        for vn in dd.keys():
            ddv = dd[vn]
            if vn == 'ocean_time':
                ot_vec[ii] = ddv
            else:
                if ddv.ndim == 2:
                    ds[vn][ii, :, :] = ddv
                elif ddv.ndim == 3:
                    ds[vn][ii, :, :, :] = ddv
                else:
                    print('problem with ndim?')
        ii += 1
    # time coordinates
    ds['ocean_time'] = (('ocean_time',), ot_vec)
    ds['ocean_time'].attrs['units'] = Lfun.roms_time_units
    ds['ocean_time'].attrs['long_name'] = 'ocean time'
    for vn in vn_dict.keys():
        vnt = vn_dict[vn][0]
        ds[vnt] = ((vnt,), ot_vec)
        ds[vnt].attrs['units'] = Lfun.roms_time_units
    # and save to NetCDF
    Enc_dict = {vn:enc_dict for vn in vn_dict.keys()}
    ds.to_netcdf(out_fn, encoding=Enc_dict)
    ds.close()

def make_ini_file(in_fn, out_fn):
    """
    Create the ini file from the first time of the clm file.
    """
    ds0 = xr.open_dataset(in_fn)
    ds = xr.Dataset()
    for vn in ds0.data_vars:
        ndims = len(ds0[vn].dims)
        if ndims == 3:
            ds[vn] = (ds0[vn].dims, ds0[vn].values[[0],:,:])
            # Note: we use [0] instead of 0 to retain the singleton dimension
        elif ndims == 4:
            ds[vn] = (ds0[vn].dims, ds0[vn].values[[0],:,:,:])
        Attrs = ds0[vn].attrs
        Attrs['long_name'] = Attrs['long_name'].replace('climatology','').strip()
        ds[vn].attrs = Attrs
    for cn in ds0.coords:
        ds.coords[cn] = ds0.coords[cn][[0]]
    ds0.close()
    Enc_dict = {vn:enc_dict for vn in ds.data_vars}
    # and save to NetCDF
    ds.to_netcdf(out_fn, encoding=Enc_dict)
    ds.close()
    
def make_bry_file(in_fn, out_fn):
    """
    Create the bry file from the edges of the clm file.
    """
    ds0 = xr.open_dataset(in_fn)
    ds = xr.Dataset()
    for vn in ds0.data_vars:
        dm = ds0[vn].dims
        ndims = len(dm)
        for D in ['north', 'south', 'east', 'west']:
            # rename variable
            Vn = vn + '_' + D
            # trim dimensions
            if D in ['east','west']:
                Dm = tuple(item for item in dm if 'xi_' not in item)
            elif D in ['north','south']:
                Dm = tuple(item for item in dm if (('eta_' not in item) or ('zeta' in item)))

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
            Attrs = ds0[vn].attrs
            Attrs['long_name'] = Attrs['long_name'].replace('climatology','').strip()
            ds[Vn].attrs = Attrs
    for cn in ds0.coords:
        ds.coords[cn] = ds0.coords[cn]
    ds0.close()
    Enc_dict = {vn:enc_dict for vn in ds.data_vars}
    # and save to NetCDF
    ds.to_netcdf(out_fn, encoding=Enc_dict)
    ds.close()

