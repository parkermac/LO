"""
Code to work with ocnN.  It gets the data_dict entries for one time and
Saves the result.  Meant to work as a parallel subprocess
"""

from lo_tools import zfun
import xarray as xr
import numpy as np
import sys
from scipy.spatial import cKDTree
import argparse
import pickle

parser = argparse.ArgumentParser()
parser.add_argument('-grid_fn', type=str) # full path to nest grid.nc file
parser.add_argument('-his_fn', type=str) # full path to history file to interpolate from
parser.add_argument('-start_type', type=str) # continuation or new
parser.add_argument('-out_fn', type=str) # full path to output pickle file
args = parser.parse_args()

def get_bounds(x_big, y_big, x_small, y_small, pad=3):
    """
    This function takes two pairs of plaid, 2-D, lon, lat arrays:
    - one pair bigger (that we hope to nest inside) and
    - one pair smaller (the grid of the nest)
    and returns the indices to use for making a trimmed version of
    the bigger grid that the smaller grid still fits inside.
    We add "pad" around the edges to make sure things fit comfortably.
    """
    # First: error checking
    if (x_small[0,0] < x_big[0,pad]) or (x_small[0,-1] > x_big[0,-pad]):
        print('ERROR: lon out of bounds ')
        sys.exit()
    if (y_small[0,0] < y_big[pad,0]) or (y_small[-1,0] > y_big[-pad,0]):
        print('ERROR: lat out of bounds ')
        sys.exit()
    # Second: get indices
    ix0 = zfun.find_nearest_ind(x_big[0,:], x_small[0,0]) - pad
    ix1 = zfun.find_nearest_ind(x_big[0,:], x_small[0,-1]) + pad
    iy0 = zfun.find_nearest_ind(y_big[:,0], y_small[0,0]) - pad
    iy1 = zfun.find_nearest_ind(y_big[:,0], y_small[-1,0]) + pad
    return ix0, ix1, iy0, iy1

tag_list = ['rho', 'u', 'v']

# the new grid
ds = xr.open_dataset(args.grid_fn)

xx = {}; yy = {}; mm = {}; xynew = {}
for tag in tag_list:
    xx[tag] = ds['lon_' + tag].values
    yy[tag] = ds['lat_' + tag].values
    mm[tag] = ds['mask_' + tag].values
    
ds.close()

if args.start_type == 'continuation':
    pad = 20
elif args.start_type == 'new':
    pad = 0
else:
    print('Error: Unrecognized start_type')
    sys.exit()
if pad > 0:
    # mask out the inside of the nest fields, since we only use
    # the edges (unless Ldir['start_type']=='new')
    for tag in tag_list:
        mm[tag][pad:-pad, pad:-pad] = 0 # this speeds things up
    
for tag in tag_list:
    xynew[tag] = np.array((xx[tag][mm[tag]==1],yy[tag][mm[tag]==1])).T

# Create 2-D search trees for the old grid
ds = xr.open_dataset(args.his_fn)
N = len(ds.s_rho.values)
xtrim = {}; ytrim = {}; mtrim = {}; xyT = {}
ix0 = {}; ix1 = {}; iy0 = {}; iy1 = {}
for tag in tag_list:
    x = ds['lon_' + tag].values
    y = ds['lat_' + tag].values
    m = ds['mask_' + tag].values # 1=water
    # trim the old grid before making the search tree
    ix0[tag], ix1[tag], iy0[tag], iy1[tag] = get_bounds(x, y, xx[tag], yy[tag])
    xtrim[tag] = x[iy0[tag]:iy1[tag], ix0[tag]:ix1[tag]]
    ytrim[tag] = y[iy0[tag]:iy1[tag], ix0[tag]:ix1[tag]]
    mtrim[tag] = m[iy0[tag]:iy1[tag], ix0[tag]:ix1[tag]]
    xyorig = np.array((xtrim[tag][mtrim[tag]==1],ytrim[tag][mtrim[tag]==1])).T
    xyT[tag] = cKDTree(xyorig)
ds.close()

# associate variables to process with grids
vn_dict = {'salt':('rho',3), 'temp':('rho',3), 'zeta':('rho',2),
        'u':('u',3), 'v':('v',3), 'ubar':('u',2), 'vbar':('v',2)}

# create blank arrays for results
data_dict = dict()
for vn in vn_dict.keys():
    tag = vn_dict[vn][0]
    dm = vn_dict[vn][1]
    NR, NC = xx[tag].shape
    if dm == 2:
        data_dict[vn] = np.nan* np.ones((NR,NC))
    elif dm == 3:
        data_dict[vn] = np.nan* np.ones((N,NR,NC))

# Interpolate to fill all data arrays for new grid.
ds = xr.open_dataset(args.his_fn, decode_times=False)
data_dict['ocean_time'] = ds.ocean_time.values[0]
for vn in vn_dict.keys():
    tag = vn_dict[vn][0]
    dm = vn_dict[vn][1]
    if dm == 2:
        vtrim = ds[vn][0,iy0[tag]:iy1[tag], ix0[tag]:ix1[tag]].values
        vv = np.nan * np.ones(xx[tag].shape) 
        vv[mm[tag]==1] = vtrim[mtrim[tag]==1][xyT[tag].query(xynew[tag], workers=-1)[1]]
        # note that "workers" has replaced "n_jobs"
        data_dict[vn][:, :] = vv
    elif dm == 3:
        for nn in range(N):
            vtrim = ds[vn][0,nn,iy0[tag]:iy1[tag], ix0[tag]:ix1[tag]].values
            vv = np.nan * np.ones(xx[tag].shape) 
            vv[mm[tag]==1] = vtrim[mtrim[tag]==1][xyT[tag].query(xynew[tag], workers=-1)[1]]
            data_dict[vn][nn, :, :] = vv
ds.close()
    
pickle.dump(data_dict, open(args.out_fn,'wb'))