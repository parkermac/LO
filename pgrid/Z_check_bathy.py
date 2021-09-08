#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 11:31:05 2017

@author: PM5

Code to check bathymetry parameters such as average depth and volume.

"""

from importlib import reload
import gfun
reload(gfun)
Gr =gfun.gstart()
import pfun
reload(pfun)
import gfun_plotting as gfp

import gfun_utility as gfu
reload(gfu)

import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np

import zfun

using_old_grid = False
# Set this to True to look at grids we have already created,
# e.g. ones currently in use for LiveOcean.
# Set it to False when interacting with grids from pgrid_output.

if using_old_grid==True:
    fn = gfun.select_file(Gr, using_old_grid=True)
    in_fn = fn
    gridname = in_fn.split('/')[-2]
elif using_old_grid==False:
    fn = gfun.select_file(Gr)
    in_fn = Gr['gdir'] + fn
    gridname = Gr['gridname']
    
dch =  gfun.default_choices(Gr)


#%% load a model grid
ds = nc.Dataset(in_fn)
z = -ds.variables['h'][:]
mask_rho = ds.variables['mask_rho'][:]

plon, plat = gfp.get_plon_plat(using_old_grid, ds)

zm = np.ma.masked_where(mask_rho == 0, z)

dch['t_list'] = ['cascadia/cascadia_gridded.nc',
         'psdem/PS_183m.nc']

lon = ds['lon_rho'][:]
lat = ds['lat_rho'][:]

#%% re-create the same bathymetry without smoothing
z0 = np.nan * lon
# m is the start of a mask: 1=water, 0=land
m = np.ones_like(lon)
for t_file in dch['t_list']:
    t_fn = dch['t_dir'] + t_file
    #print('\nOPENING BATHY FILE: ' + t_file)
    tlon_vec, tlat_vec, tz = gfu.load_bathy_nc(t_fn)
    if isinstance(tz, np.ma.masked_array):
        tz1 = tz.data
        tz1[tz.mask==True] = np.nan
        tz = tz1
    tlon, tlat = np.meshgrid(tlon_vec, tlat_vec)
    z_part = zfun.interp2(lon, lat, tlon, tlat, tz)
    # put good values of z_part in z
    z0[~np.isnan(z_part)] = z_part[~np.isnan(z_part)]
if dch['use_z_offset']:
    z0 = z0 + dch['z_offset']
    
z0[z0>0] = 0

z0m = np.ma.masked_where(mask_rho==0, z0)

dz = z0m - zm
#%% plotting

#plt.close('all')

fig = plt.figure(figsize=(12,12))
ax = fig.add_subplot(111)
vv = 50
cs = plt.pcolormesh(plon, plat, dz, vmin=-vv, vmax=vv, cmap='rainbow')
aa = [-124, -122, 47, 49]
ax.axis(aa)
pfun.dar(ax)
fig.colorbar(cs)
ax.set_title(gridname)

plt.show()

#%% calculations
i0 = zfun.find_nearest_ind(lon[0,:], aa[0])
i1 = zfun.find_nearest_ind(lon[0,:], aa[1])
j0 = zfun.find_nearest_ind(lat[:,0], aa[2])
j1 = zfun.find_nearest_ind(lat[:,0], aa[3])

dx = 1/ds['pm'][:]
dy = 1/ds['pn'][:]
da = dx*dy

Z = zm[j0:j1+1, i0:i1+1]
Z0 = z0m[j0:j1+1, i0:i1+1]
DA = da[j0:j1+1, i0:i1+1]
M = mask_rho[j0:j1+1, i0:i1+1]
DAM = np.ma.masked_where(M==0, DA)

V = np.sum(-Z*DA)
V0 = np.sum(-Z0*DA)
A = np.sum(DAM)
H = V/A
H0 = V0/A

print('\n*********** ' + gridname + ' ***********')
print('Original volume = %d km' % (int(V0/1e9)))
print('   Final volume = %d km' % (int(V/1e9)))
print('Original mean depth = %d m' % (int(H0)))
print('   Final mean depth = %d m' % (int(H)))
