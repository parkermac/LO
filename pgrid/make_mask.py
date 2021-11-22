"""
Code to create an initial mask for a grid.
"""

import numpy as np
import pickle
import xarray as xr

from lo_tools import zfun
from lo_tools import plotting_functions as pfun

import gfun

Gr =gfun.gstart()
# select and increment grid file
in_fn = gfun.select_file(Gr)
out_fn = gfun.increment_filename(in_fn, '_m')

# get the grid from NetCDF
ds = xr.open_dataset(in_fn)
lon = ds.lon_rho.values
lat = ds.lat_rho.values
plon, plat = pfun.get_plon_plat(lon, lat)
z = -ds.h.values
mask_rho_orig = ds.mask_rho.values
plon_vec = plon[0,:]
plat_vec = plat[:,0]
lon_vec = lon[0,:]
lat_vec = lat[:,0]

# load the default choices
dch = pickle.load(open(Gr['gdir'] / 'choices.p', 'rb'))

def mask_from_existing(in_fn, maskfiles, pgdir):
    ds = xr.open_dataset(in_fn)
    lonr = ds['lon_rho'].values
    latr = ds['lat_rho'].values
    ds.close()
    
    m10 = np.ones(lonr.shape)
        
    for mf in maskfiles:
        print(' - using ' + pgdir + mf)
        ds = xr.open_dataset(pgdir + mf)
        xx = ds['lon_rho'].values
        yy = ds['lat_rho'].values
        mm = ds['mask_rho'].values
        # mask_rho = 1. = water, and 0. = land
        m_part = zfun.interp2(lonr, latr, xx, yy, mm)
        m10[m_part < .5] = 0
        ds.close()
    m = m10 == 0 # boolean array (False = water, True = land)
    return m

# PROCESSING

# Create a boolean mask array (True where masked = land)
# following the numpy masked array convention.
# Note that this is the opposite of the ROMS convention
# where mask_rho = 1. over water, and 0. over land.
if mask_rho_orig.all() == 1:    
    print('Original mask all ones')
    if len(dch['maskfiles']) == 0:
        # set z position of initial dividing line (positive up)
        m = z >= dch['z_land']
    else:
        print('using maskfiles')
        m = mask_from_existing(in_fn, dch['maskfiles'], Gr['pgdir'])
        m[z >= dch['z_land']] = True
        
        print(m.shape)
        print(np.sum(m))
        print(np.sum(~m))
        
# unmask the coast
if dch['unmask_coast']:
    # This unmasks it in the places where the
    # coastline crosses a tile, to facilitate wetting-drying
    cx, cy = pfun.get_coast()
    cmask = np.isnan(cx)
    cx = cx[~cmask]
    cy = cy[~cmask]
    ii0, ii1, ifr = zfun.get_interpolant(cx, lon_vec)
    jj0, jj1, jfr = zfun.get_interpolant(cy, lat_vec)
    # Don't unmask extrapolated points.
    ii0 = ii0[~np.isnan(ifr) & ~np.isnan(jfr)]
    jj0 = jj0[~np.isnan(ifr) & ~np.isnan(jfr)]
    m[jj0, ii0] = False
      
# remove islands and lakes
if dch['remove_islands']:
    # What this does is mask any water point that has land on 3 sides
    # or any land point that has water on three sides. By doing this repeatedly
    # you get rid of stray channels or peninsulas.
    # The number in range() determines how long of a feature is removed.
    # What the algorithm will not do, for example, is get rid of
    # a square lake of 4 cells.    
    for ii in range(7): # was range(5)
        NR, NC = m.shape
        mm = m[1:-1, 1:-1]
        mn = m[2:, 1:-1]
        ms = m[:-2, 1:-1]
        me = m[1:-1, 2:]
        mw = m[1:-1, :-2]
        # remove islands of ocean
        MMo = ~mm & mn & ms & me
        mm[MMo] = True
        MMo = ~mm & mn & ms & mw
        mm[MMo] = True
        MMo = ~mm & mn & me & mw
        mm[MMo] = True
        MMo = ~mm & ms & me & mw
        mm[MMo] = True
        # remove islands of land
        MMl = mm & ~mn & ~ms & ~me
        mm[MMl] = False
        MMl = mm & ~mn & ~ms & ~mw
        mm[MMl] = False
        MMl = mm & ~mn & ~me & ~mw
        mm[MMl] = False
        MMl = mm & ~ms & ~me & ~mw
        mm[MMl] = False
        m[1:-1, 1:-1] = mm
        
# Save the output file

# create the new mask_rho
# 1 = water
# 0 = land
mask_rho = np.ones(mask_rho_orig.shape)
mask_rho[m == True] = 0

if not np.all(mask_rho == mask_rho_orig):
    print('Creating ' + str(out_fn))
    ds.update({'mask_rho':(('eta_rho', 'xi_rho'), mask_rho)})
    ds.to_netcdf(out_fn)
    ds.close()
else:
    print('No change to mask')
