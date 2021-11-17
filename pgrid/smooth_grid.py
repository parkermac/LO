"""
Smooth the grid.

This is slow-ish. I wish I could vectorize the
gfu.GRID_PlusMinusScheme_rx0().  But still it only takes a
minute for a grid with 3 million points, so not bad.
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

# Make sure that anything NOT MASKED is not shallower than min_depth.
if dch['use_min_depth']:
    # Hobs[(MSK==1) & (Hobs < dch['min_depth'])] = dch['min_depth']
    Hobs[Hobs < dch['min_depth']] = dch['min_depth']

# create the area matrix
AreaMatrix = dx * dy

# create smoothed bathymetry
tt0 = time()
# Shift the whole depth grid for the case where we have active bathymetry
# on land (e.g. when using wet_dry), because otherwise the smoothing
# code will fail.  We shift it back at the end.  All the shifting is done in
# gfu.GRID_PlusMinusScheme_rx0().
if dch['min_depth'] > 0:
    shift = 0
elif dch['min_depth'] <= 0:
    print('min_depth must be > 0')
    sys.exit()
    
# Do the smoothing.
Hnew = gfu.GRID_PlusMinusScheme_rx0(MSK, Hobs, rx0max, AreaMatrix,
            fjord_cliff_edges=True, shift=shift)

print('Smoothing took %0.1f seconds' % (time() - tt0))

# Again, make sure that anything NOT MASKED is not shallower than min_depth.
if dch['use_min_depth']:
    Hnew[(MSK==1) & (Hnew < dch['min_depth'])] = dch['min_depth']

# save the updated mask and h
ds.update({'mask_rho': (('eta_rho', 'xi_rho'), MSK)})
ds.update({'h': (('eta_rho', 'xi_rho'), Hnew)})
ds.to_netcdf(out_fn)
ds.close()

