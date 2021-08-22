#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code to convert bathymetry files from matlab to netcdf.

Takes only a few seconds to run.
"""

import netCDF4 as nc
import numpy as np

import os
dir0 = os.environ.get('HOME') + '/Documents/'


from importlib import reload
import gfun_utility as gfu
reload(gfu)

tdir0 = dir0 + 'ptools_data/topo/'

tdict = {'srtm15': ['topo15.grd'],
         'cascadia': ['cascadia_gridded.mat'],
         'psdem': ['PS_183m.mat', 'PS_91m.mat', 'PS_27m.mat'],
         'ttp_patch': ['TTP_Regional_27m_patch.mat']}

for k in tdict.keys():
    
    for in_fn in tdict[k]:
        print(50*'=')
        print(in_fn)
        t_fn = tdir0 + k + '/' + in_fn
        if k == 'srtm15':
            lon_vec = np.array([-135, -122])
            lat_vec = np.array([35, 53])
            lon, lat, z = gfu.load_bathy2(t_fn, lon_vec, lat_vec)
        else:
            lon, lat, z = gfu.load_bathy(t_fn)
        z = np.ma.masked_where(np.isnan(z), z)
        # write the data to NetCDF
        out_fn = t_fn[:-4] + '.nc'
        foo = nc.Dataset(out_fn, 'w')
        L = len(lon)
        M = len(lat)
        foo.createDimension('x', L)
        foo.createDimension('y', M)
        lon_var = foo.createVariable('lon', float, ('x'))
        lat_var = foo.createVariable('lat', float, ('y'))
        z_var = foo.createVariable('z', float, ('y', 'x'))
        lon_var[:] = lon
        lat_var[:] = lat
        z_var[:] = z
        foo.close()

