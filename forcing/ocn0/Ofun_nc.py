#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 13:51:06 2017

@author: PM5

Functions for writing saved fields to ROMS NetCDF.
"""

import os
import netCDF4 as nc
import numpy as np
import Lfun

ncformat = 'NETCDF3_64BIT_OFFSET'

def make_clm_file(Ldir, nc_dir, fh_dir, c_dict, dt_list, S, G):
    # name output file
    clm_fn = nc_dir + 'ocean_clm.nc'
    # get rid of the old version, if it exists
    try:
        os.remove(clm_fn)
    except OSError:
        pass # assume error was because the file did not exist
    foo = nc.Dataset(clm_fn, 'w', format=ncformat)
    # create dimensions
    for vn in ['salt', 'temp', 'v3d', 'v2d', 'zeta', 'ocean']:
        foo.createDimension(vn+'_time', len(dt_list)) # could use None
    foo.createDimension('s_rho', S['N'])
    for tag in ['rho', 'u', 'v']:
        foo.createDimension('eta_'+tag, G['lat_'+tag].shape[0])
        foo.createDimension('xi_'+tag, G['lon_'+tag].shape[1])
    # add time data
    dtm_list = []
    for dt in dt_list:
        dtm_list.append(Lfun.datetime_to_modtime(dt))        
    for vn in ['salt', 'temp', 'v3d', 'v2d', 'zeta', 'ocean']:
        vv = foo.createVariable(vn+'_time', float, (vn+'_time',))
        vv.units = 'seconds since 1970.01.01 UTC'
        vv[:] = np.array(dtm_list)
    # add 2d field data
    vv = foo.createVariable('zeta', float, ('zeta_time', 'eta_rho', 'xi_rho'))
    vv.long_name = 'sea surface height climatology'
    vv.units = 'meter'
    for t in c_dict.keys():
        c = c_dict[t]
        vv[t,:,:] = c['ssh']
    vv = foo.createVariable('ubar', float, ('v2d_time', 'eta_u', 'xi_u'))
    vv.long_name = 'vertically averaged u-momentum climatology'
    vv.units = 'meter second-1'
    for t in c_dict.keys():
        c = c_dict[t]
        vv[t,:,:] = c['ubar']
    vv = foo.createVariable('vbar', float, ('v2d_time', 'eta_v', 'xi_v'))
    vv.long_name = 'vertically averaged v-momentum climatology'
    vv.units = 'meter second-1'
    for t in c_dict.keys():
        c = c_dict[t]
        vv[t,:,:] = c['vbar']
    # add 3d field data
    vv = foo.createVariable('u', float, ('v3d_time', 's_rho', 'eta_u', 'xi_u'))
    vv.long_name = 'u-momentum component climatology'
    vv.units = 'meter second-1'
    for t in c_dict.keys():
        c = c_dict[t]
        vv[t,:,:,:] = c['u3d']
    vv = foo.createVariable('v', float, ('v3d_time', 's_rho', 'eta_v', 'xi_v'))
    vv.long_name = 'v-momentum component climatology'
    vv.units = 'meter second-1'
    for t in c_dict.keys():
        c = c_dict[t]
        vv[t,:,:,:] = c['v3d']
    vv = foo.createVariable('salt', float, ('salt_time', 's_rho', 'eta_rho', 'xi_rho'))
    vv.long_name = 'salinity climatology'
    vv.units = 'PSU'
    for t in c_dict.keys():
        c = c_dict[t]
        vv[t,:,:,:] = c['s3d']
    vv = foo.createVariable('temp', float, ('temp_time', 's_rho', 'eta_rho', 'xi_rho'))
    vv.long_name = 'potential temperature climatology'
    vv.units = 'Celsius'
    for t in c_dict.keys():
        c = c_dict[t]
        vv[t,:,:,:] = c['theta']
    foo.close()
    print('-Writing ocean_clm.nc')
        
def make_ini_file(nc_dir):
    # Initial condition, copied from first time of ocean_clm.nc
    ds1 = nc.Dataset(nc_dir + 'ocean_clm.nc', mode='r')
    # name output file
    ini_fn = nc_dir + 'ocean_ini.nc'
    # get rid of the old version, if it exists
    try:
        os.remove(ini_fn)
    except OSError:
        pass # assume error was because the file did not exist
    ds2 = nc.Dataset(ini_fn, 'w', format=ncformat)
    # Copy dimensions
    for dname, the_dim in ds1.dimensions.items():
        if 'time' in dname:
            ds2.createDimension(dname, 1)
        else:
            ds2.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)
    # Copy variables
    for v_name, varin in ds1.variables.items():
        outVar = ds2.createVariable(v_name, varin.datatype, varin.dimensions)
        # Copy variable attributes, {} is a dict comprehension, cool!
        outVar.setncatts({k: varin.getncattr(k).replace('climatology','').strip() for k in varin.ncattrs()})
        if varin.ndim > 1:
            outVar[:] = varin[0,:]
        else:
            outVar[:] = varin[0]
    ds1.close()
    ds2.close()
    print('-Writing ocean_ini.nc')

def make_bry_file(nc_dir):
    # Boundary conditions, copied from edges of ocean_clm.nc
    ds1 = nc.Dataset(nc_dir + 'ocean_clm.nc', mode='r')
    # name output file
    bry_fn = nc_dir + 'ocean_bry.nc'
    # get rid of the old version, if it exists
    try:
        os.remove(bry_fn)
    except OSError:
        pass # assume error was because the file did not exist
    ds2 = nc.Dataset(bry_fn, 'w', format=ncformat)
    # Copy dimensions
    for dname, the_dim in ds1.dimensions.items():
        ds2.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)
    # Copy parts of variables
    for v_name, varin in ds1.variables.items():
        if varin.ndim in [3,4]: 
            for bname in ['north', 'south', 'east', 'west']:
                # hacks to get the right boundary names for a few BGD variables
                # (really this should be regularized in the ROMS NZD code)
                if v_name == 'phytoplankton':
                    v_name = 'phyt'
                if v_name == 'zooplankton':
                    v_name = 'zoop'
                if v_name == 'alkalinity':
                    v_name = 'Talk'
                outname = v_name + '_' + bname
                if bname in ['north', 'south']:
                    outdims = tuple([item for item in varin.dimensions if item[:3] != 'eta'])
                elif bname in ['west', 'east']:
                    outdims = tuple([item for item in varin.dimensions if item[:2] != 'xi'])
                outVar = ds2.createVariable(outname, varin.datatype, outdims)    
                outVar.setncatts({k: varin.getncattr(k).replace('climatology','').strip() for k in varin.ncattrs()})
                if varin.ndim == 4:
                    if bname == 'north':
                        outVar[:] = varin[:,:,-1,:]
                    elif bname == 'south':
                        outVar[:] = varin[:,:,0,:]
                    elif bname == 'east':
                        outVar[:] = varin[:,:,:,-1]
                    elif bname == 'west':
                        outVar[:] = varin[:,:,:,0]
                elif varin.ndim == 3:
                    if bname == 'north':
                        outVar[:] = varin[:,-1,:]
                    elif bname == 'south':
                        outVar[:] = varin[:,0,:]
                    elif bname == 'east':
                        outVar[:] = varin[:,:,-1]
                    elif bname == 'west':
                        outVar[:] = varin[:,:,0]
        else:
            outname = v_name
            outdims = tuple([item for item in varin.dimensions])
            outVar = ds2.createVariable(outname, varin.datatype, outdims)    
            outVar.setncatts({k: varin.getncattr(k).replace('climatology','').strip() for k in varin.ncattrs()})
            outVar[:] = varin[:]    
    ds1.close()
    ds2.close()
    print('-Writing ocean_bry.nc')
