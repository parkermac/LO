# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 15:34:17 2016

@author: PM5

Code to initialize the creation of a ROMS grid file.

NOTE: the gridname is set in gfun.gstart().

Throughout this code I try to use ROMS naming conventions, except that
when manipulating or plotting I refer to [lon,lat]_rho as [lon,lat],
and [lon,lat]_psi_ex as [plon,plat].
"""

from importlib import reload
import gfun
reload(gfun)
Gr =gfun.gstart()
import gfun_utility as gfu
reload(gfu)

import numpy as np
import pickle

import Lfun
import zfun

Lfun.make_dir(Gr['gdir'], clean=True)

fn = 'grid_m00_r00_s00_x00.nc'
out_fn = Gr['gdir'] + fn
print(50*'*')
print(out_fn)

dch = gfun.default_choices(Gr)

# GRID DEFINITIONS

# vectors to define the plaid grid
# start with cell corners (like an extended psi grid)
    
if Gr['gridname'] == 'sal0':
    # start of a salish nest grid
    aa = [-124, -122, 47, 49]
    res = 300 # target resolution (m)
    plon_vec, plat_vec = gfu.simple_grid(aa, res)
    dch['t_list'] = ['cascadia/cascadia_gridded.nc',
             'psdem/PS_183m.nc',
             'ttp_patch/TTP_Regional_27m_patch.nc']
    dch['nudging_edges'] = ['north', 'west']

elif Gr['gridname'] == 'hc1':
    # mid Hood Canal nest
    aa = [-123, -122.55, 47.5, 47.9]
    res = 100 # target resolution (m)
    plon_vec, plat_vec = gfu.simple_grid(aa, res)
    dch['t_list'] = ['psdem/PS_27m.nc']
    dch['nudging_edges'] = ['north', 'west']
    dch['nudging_days'] = (0.1, 1.0)

elif Gr['gridname'] == 'cas3': # a stretched MoSSea-like grid, became cas4/5
    maxres = 1500
    minres = 500
    lon_list = [-127.4, -124, -122]
    x_res_list = [maxres, minres, minres]
    lat_list = [42, 47, 49, 50.3]
    y_res_list = [maxres, minres, minres, maxres]
    plon_vec, plat_vec = gfu.stretched_grid(lon_list, x_res_list,
                                        lat_list, y_res_list)
    dch['nudging_edges'] = ['south', 'west']
    # new: a list of good masks to work from
    dch['maskfiles'] = ['cas2/grid_m05_r01_s01_x01.nc', 'sal0/grid_m06_r03_s05_x02.nc']

elif Gr['gridname'] == 'cas6': # an extended version of the excellent cas4/5
    maxres = 1500
    minres = 500
    extres = 3000
    lon_list = [-130, -127.4, -124, -122]
    x_res_list = [extres, maxres, minres, minres]
    lat_list = [42, 47, 49, 50.3, 52]
    y_res_list = [maxres, minres, minres, maxres, extres]
    plon_vec, plat_vec = gfu.stretched_grid(lon_list, x_res_list,
                                        lat_list, y_res_list)
    dch['nudging_edges'] = ['north', 'south', 'west']
    # new: a list of good masks to work from
    dch['maskfiles'] = ['cas5/grid_m05_r01_s02_x02.nc']
    
elif Gr['gridname'] == 'sj0':
    # San Juan Islands nest, first version
    aa = [-123.3, -122.65, 48.3, 48.8]
    res = 100 # target resolution (m)
    plon_vec, plat_vec = gfu.simple_grid(aa, res)
    dch['t_list'] = ['cascadia/cascadia_gridded.nc']
    dch['nudging_edges'] = ['north', 'south', 'east', 'west']
    dch['nudging_days'] = (0.1, 1.0)

elif Gr['gridname'] == 'aestus1': # idealized model
    lon_list = [-1, 0, 1, 2, 3]
    x_res_list = [5000, 1000, 1000, 5000, 5000]
    lat_list = [44, 44.9, 45.1, 46]
    y_res_list = [5000, 1000, 1000, 5000]
    plon_vec, plat_vec = gfu.stretched_grid(lon_list, x_res_list,
                                        lat_list, y_res_list)
    dch['analytical'] = True
    dch['nudging_edges'] = ['north', 'south', 'west']
    
elif Gr['gridname'] == 'aestus2':
    # idealized model, higher resolution and larger domain than aestus1
    lon_list = [-2, 0, 1, 2, 3]
    x_res_list = [2500, 500, 500, 2500, 2500]
    lat_list = [43, 44.9, 45.1, 47]
    y_res_list = [2500, 500, 500, 2500]
    plon_vec, plat_vec = gfu.stretched_grid(lon_list, x_res_list,
                                        lat_list, y_res_list)
    dch['analytical'] = True
    dch['nudging_edges'] = ['north', 'south', 'west']
    
elif Gr['gridname'] == 'aestus3':
    # idealized model, like aestus2 but cut off just inside the estuary mouth,
    # with the goal of more manipulative forcing of the exchange flow
    lon_list = [-2, 0, .1]
    x_res_list = [2500, 500, 500]
    lat_list = [43, 44.9, 45.1, 47]
    y_res_list = [2500, 500, 500, 2500]
    plon_vec, plat_vec = gfu.stretched_grid(lon_list, x_res_list,
                                        lat_list, y_res_list)
    dch['analytical'] = True
    dch['nudging_edges'] = ['north', 'south', 'west']

# save the default choices for use by other code
pickle.dump(dch, open(Gr['gdir'] + 'choices.p', 'wb'))

plon, plat = np.meshgrid(plon_vec, plat_vec)
#ax_lims = (plon_vec[0], plon_vec[-1], plat_vec[0], plat_vec[-1])

# make box centers
lon_vec = plon_vec[:-1] + np.diff(plon_vec)/2
lat_vec = plat_vec[:-1] + np.diff(plat_vec)/2
lon, lat = np.meshgrid(lon_vec, lat_vec)
NR, NC = lon.shape

# initialize the final bathymetry array
z = np.nan * lon
if dch['analytical']==True:
    if Gr['gridname'] in ['aestus1', 'aestus2', 'aestus3']:
        # make grid and bathymetry by hand
        z = np.zeros(lon.shape)
        x, y = zfun.ll2xy(lon, lat, 0, 45)
        zshelf = x * 1e-3
        zestuary = -20 + 20*x/1e5 + 20/(1e4)*np.abs(y)
        z = zshelf
        mask = zestuary < z
        z[mask] = zestuary[mask]
        m = np.ones_like(lon)
else:
    # add bathymetry automatically from files
    # m is the start of a mask: 1=water, 0=land
    m = np.ones_like(lon)
    for t_file in dch['t_list']:
        t_fn = dch['t_dir'] + t_file
        print('\nOPENING BATHY FILE: ' + t_file)
        tlon_vec, tlat_vec, tz = gfu.load_bathy_nc(t_fn)
        if isinstance(tz, np.ma.masked_array):
            tz1 = tz.data
            tz1[tz.mask==True] = np.nan
            tz = tz1
        tlon, tlat = np.meshgrid(tlon_vec, tlat_vec)
        z_part = zfun.interp2(lon, lat, tlon, tlat, tz)
        # put good values of z_part in z
        z[~np.isnan(z_part)] = z_part[~np.isnan(z_part)]
    if dch['use_z_offset']:
        z = z + dch['z_offset']

#%% save the output to NetCDF
gfu.make_nc(out_fn, plon, plat, lon, lat, z, m, dch)



