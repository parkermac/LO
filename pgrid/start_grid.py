"""
Code to initialize the creation of a ROMS grid file.

NOTE: the gridname is set in gfun.gstart().

Throughout this code I try to use ROMS naming conventions, except that
when manipulating or plotting I refer to [lon,lat]_rho as [lon,lat],
and [lon,lat]_psi_ex as [plon,plat].
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
    # start of a salish nest grid
    aa = [-124, -122, 47, 49]
    res = 600 # target resolution (m)
    plon_vec, plat_vec = gfu.simple_grid(aa, res)
    # t_dir = Ldir['data'] / 'topo'
    # dch['t_list'] = [
    #           t_dir / 'cascadia' / 'cascadia_gridded.nc',
    #          t_dir / 'psdem' / 'PS_183m.nc',
    #          t_dir / 'ttp_patch' / 'TTP_Regional_27m_patch.nc']
    dch['nudging_edges'] = ['north', 'west']

# save the default choices for use by other code
pickle.dump(dch, open(Gr['gdir'] / 'choices.p', 'wb'))

plon, plat = np.meshgrid(plon_vec, plat_vec)

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
    for t_fn in dch['t_list']:
        print('\nOPENING BATHY FILE: ' + t_fn.name)
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



