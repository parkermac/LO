"""
Code to initialize the creation of a ROMS grid file.

NOTE: the gridname is set in gfun.gstart().

Throughout this code I try to use ROMS naming conventions, except that
when manipulating or plotting I refer to: [lon_rho,lat_rho] as [lon,lat].

Also [plon,plat] is just like [lon_psi, lat_psi] but extended by one on all
directions so that it is box corners around all rho-grid points.

"""
import numpy as np
import pickle

import gfun
import gfun_utility as gfu
from importlib import reload
reload(gfun)
reload(gfu)
from lo_tools import Lfun, zfun

Gr =gfun.gstart()

Lfun.make_dir(Gr['gdir'], clean=True)

fn = 'grid_m00_r00_s00_x00.nc'
out_fn = Gr['gdir'] / fn
print(60*'*')
print(str(out_fn).center(60,'-'))

dch = gfun.default_choices(Gr)

# GRID DEFINITIONS

# vectors to define the plaid grid
# start with cell corners (like an extended psi grid)
    
if Gr['gridname'] == 'sal0':
    # A Salish Sea grid, used as an example.
    aa = [-124, -122, 47, 49]
    res = 600 # target resolution (m)
    Lon_vec, Lat_vec = gfu.simple_grid(aa, res)
    dch['nudging_edges'] = ['north', 'west']

# save the default choices for use by other code
pickle.dump(dch, open(Gr['gdir'] / 'choices.p', 'wb'))

# the rho grid
lon, lat = np.meshgrid(Lon_vec, Lat_vec)
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
    for t_fn in dch['t_list']:
        print('\nOPENING BATHY FILE: ' + t_fn.name)
        tlon_vec, tlat_vec, tz = gfu.load_bathy_nc(t_fn)
        tlon, tlat = np.meshgrid(tlon_vec, tlat_vec)
        z_part = zfun.interp2(lon, lat, tlon, tlat, tz)
        # put good values of z_part in z
        z[~np.isnan(z_part)] = z_part[~np.isnan(z_part)]
    if dch['use_z_offset']:
        z = z + dch['z_offset']

# save the output to NetCDF
gfu.make_nc(out_fn, lon, lat, z, m, dch)



