"""
Code to test the nearest neighbor search algorithm.

RESULT: the algorithm appears to work perfectly, finding the exact
points requested (with nudge=0).  With nudge != 0 the error is equal
to the nudge I added to the search points.

"""

import os, sys
sys.path.append(os.path.abspath('../../LiveOcean/alpha'))
import Lfun
import zrfun

from scipy.spatial import cKDTree
from time import time
import pickle
import numpy as np
import argparse

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# optional command line arguments, can be input in any order
parser = argparse.ArgumentParser()
parser.add_argument('-gridname', default='cas6', type=str)
parser.add_argument('-gtagex', default='cas6_v3_lo8b', type=str)
parser.add_argument('-ds', default='2019.07.04', type=str)
# set which roms output directory to look in (refers to Ldir['roms'] or Ldir['roms2'])
parser.add_argument('-rd', '--roms_dir', default='roms', type=str)
# valid arguments to pass are: roms, roms2 (see alpha/get_lo_info.sh)
args = parser.parse_args()
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Ldir = Lfun.Lstart()
Ldir['roms'] = Ldir[args.roms_dir]
fn = Ldir['roms'] + 'output/' + args.gtagex + '/f' + args.ds + '/ocean_his_0001.nc'
outdir0 = Ldir['LOo'] + 'tracker_trees/'
Lfun.make_dir(outdir0)
outdir = outdir0 + args.gridname + '/'

xyzT_rho = pickle.load(open(outdir + 'xyzT_rho.p', 'rb'))

indir = outdir
in_fn = indir + 'xyzT_rho.p'

G, S, T = zrfun.get_basic_info(fn)
h = G['h']

z = zrfun.get_z(h, 0*h, S, only_rho=True)

tag = 'rho'
x = G['lon_' + tag]
y = G['lat_' + tag]
mask = G['mask_' + tag]

N,M,L = z.shape
X = np.tile(x.reshape(1,M,L),[N,1,1])
Y = np.tile(y.reshape(1,M,L),[N,1,1])
H = np.tile(h.reshape(1,M,L),[N,1,1])
Z = z/H # fractional depth (-1 to 0)
# we use these to make the arrays to query

Mask = np.tile(mask.reshape(1,M,L),[N,1,1])

# we subsample this to find points to look for
xyz = np.array((X[Mask],Y[Mask],Z[Mask])).T
# xyz.shape gives (12834660, 3), or about 13 million points for each field

# find grid spacing
dx = np.diff(x[0,:]).max()
dy = np.diff(y[:,0]).max()
dz = np.diff(Z[:,0,0]).max()
# I don't use these but could use them to generate a little variation
# to the query field xyz_q.  One interesting thing is that they are all similar
# in size dx = 0.04, dy = 0.03, dz = 0.06, meaning that the search space is
# not terribly flat in any direction, I think

# these are the points at which we will query our arrays
nudge = 0.001 # intended to be less than a grid spacing
xyz_q = xyz[::1000,:] + nudge # subsample to find points to query, and add a little nudge

# these are the arrays to query to see if the search is working
xm = X[Mask].data
ym = Y[Mask].data
zm = Z[Mask].data

# these are the answers to the query, for all three axes
x_a = xm[xyzT_rho.query(xyz_q, n_jobs=-1)[1]]
y_a = ym[xyzT_rho.query(xyz_q, n_jobs=-1)[1]]
z_a = zm[xyzT_rho.query(xyz_q, n_jobs=-1)[1]]
# and they represent the nearest positions to the query

# checking on the results
x_err = (xyz_q[:,0] - x_a).max()
y_err = (xyz_q[:,1] - y_a).max()
z_err = (xyz_q[:,2] - z_a).max()
# if it is working these should all be identical to nudge

print('Errors: x=%0.6f, y=%0.6f, z = %0.6f' % (x_err, y_err, z_err))



