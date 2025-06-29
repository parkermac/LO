"""
Smooth the grid.

"""
import numpy as np
import pickle
import xarray as xr
from time import time

from lo_tools import zfun
from lo_tools import plotting_functions as pfun

import gfun
import gfun_utility as gfu

Gr =gfun.gstart()
# select and increment grid file
in_fn = gfun.select_file(Gr)
out_fn = gfun.increment_filename(in_fn, '_s')

# get the grid from NetCDF
ds = xr.open_dataset(in_fn)
plon, plat = pfun.get_plon_plat(ds.lon_rho.values, ds.lat_rho.values)
h = ds.h.values
mask_rho = ds.mask_rho.values
plon_vec = plon[0,:]
plat_vec = plat[:,0]
dx = 1/ds.pm.values
dy = 1/ds.pn.values

# load the default choices
dch = pickle.load(open(Gr['gdir'] / 'choices.p', 'rb'))

# Smoothing set up
MSK = mask_rho.copy()
Hobs = h.copy()
rx0max = dch['rx0max']

# Make sure the entire grid is not shallower than min_depth.
# Note that ROMS, and the smoothing algorithm, both require that the
# depth of unmasked cells be > 0.
if dch['use_min_depth']:
    if dch['min_depth'] > 0:
        pass
    elif dch['min_depth'] <= 0:
        print('min_depth must be > 0')
        sys.exit()
    Hobs[Hobs < dch['min_depth']] = dch['min_depth']

# create the area matrix
AreaMatrix = dx * dy

# create smoothed bathymetry
tt0 = time()
# New scheme to make sure the intertidal is smooth even at low tide.
# We will do the smoothing several times, starting by shifting it up
# and gradually moving to zero shift, using a temporary mask each time.
Htemp = Hobs.copy()
for shift in [6,5,4,3,2,1,0]:
    Htemp = Htemp - shift
    MSKtemp=MSK.copy()
    # Mask convention:
    # 1 = water
    # 0 = land
    MSKtemp[Htemp < dch['min_depth']] = 0
    Hnew = gfu.GRID_PlusMinusScheme_rx0(MSKtemp, Htemp, rx0max, AreaMatrix,
            fjord_cliff_edges=True)
    Htemp = Hnew + shift
print('Smoothing took %0.1f seconds' % (time() - tt0))

# Make sure the entire grid is not shallower than min_depth, again.
if dch['use_min_depth']:
    if dch['min_depth'] > 0:
        pass
    elif dch['min_depth'] <= 0:
        print('min_depth must be > 0')
        sys.exit()
    Hnew[Hnew < dch['min_depth']] = dch['min_depth']

# save the updated mask and h
ds.update({'mask_rho': (('eta_rho', 'xi_rho'), MSK)})
ds.update({'h': (('eta_rho', 'xi_rho'), Hnew)})
ds.to_netcdf(out_fn)
ds.close()

