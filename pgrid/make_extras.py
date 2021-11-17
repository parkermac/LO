"""
Add the other masks.

"""
import numpy as np
import pickle
import xarray as xr

import gfun

Gr =gfun.gstart()
# select and increment grid file
in_fn = gfun.select_file(Gr)
out_fn = gfun.increment_filename(in_fn, '_x')

# get the grid from NetCDF
ds = xr.open_dataset(in_fn)
h = ds.h.values
mask_rho = ds.mask_rho.values

# load the default choices
dch = pickle.load(open(Gr['gdir'] / 'choices.p', 'rb'))

# enforce min depth
if dch['use_min_depth']:
    # set min depth everywhere
    h[ h <= dch['min_depth'] ] = dch['min_depth']
    # also set it to min_depth where masked
    h[mask_rho == 0] = dch['min_depth']

# Make the masks.

mask_u_bool = (mask_rho[:, 1:] == 0) | (mask_rho[:, :-1] == 0)
mask_u = np.ones(mask_u_bool.shape)
mask_u[mask_u_bool] = 0

mask_v_bool = (mask_rho[1:, :] == 0) | (mask_rho[:-1, :] == 0)
mask_v = np.ones(mask_v_bool.shape)
mask_v[mask_v_bool] = 0

mask_psi_bool = ( (mask_rho[1:, 1:] == 0) | (mask_rho[:-1, :-1] == 0) |
                (mask_rho[1:, :-1] == 0) | (mask_rho[:-1, 1:] == 0) )
mask_psi = np.ones(mask_psi_bool.shape)
mask_psi[mask_psi_bool] = 0

# save the updated mask and h
ds.update({'h': (('eta_rho', 'xi_rho'), h)})
ds.update({'mask_u': (('eta_u', 'xi_u'), mask_u)})
ds.update({'mask_v': (('eta_v', 'xi_v'), mask_v)})
ds.update({'mask_psi': (('eta_psi', 'xi_psi'), mask_psi)})

enc_dict = dict()
for vn in ds.data_vars:
    if vn not in ['spherical', 'sph']:
        enc_dict[vn] = {'_FillValue': 1e20}
    else:
        enc_dict[vn] = {'dtype': 'S1'}
ds.to_netcdf(out_fn, encoding=enc_dict)
ds.close()



