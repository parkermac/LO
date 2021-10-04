"""
Creates and saves KDTrees for a given model run.

NOTE: for some reason I have to make the 2D trees individually (meaning
not in a loop like I do for the 3D trees).  No idea why.

Usage:

python make_KDTrees.py (uses defaults)

python make_KDTrees.py -gridname cas6 -gtagex cas6_v3_lo8b -ds 2019.07.04

(** you have to be pointing it to a history file that exists **)

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
Lfun.make_dir(outdir, clean=True)

G, S, T = zrfun.get_basic_info(fn)
h = G['h']

# 2D trees
X = G['lon_rho']; Y = G['lat_rho']
Maskr = G['mask_rho'] # True over water
xy = np.array((X[Maskr],Y[Maskr])).T
xyT_rho = cKDTree(xy)
xy = np.array((X.flatten(),Y.flatten())).T
xyT_rho_un = cKDTree(xy) # unmasked version
X = G['lon_u']; Y = G['lat_u']
Masku = G['mask_u'] # True over water
xy = np.array((X[Masku],Y[Masku])).T
xyT_u = cKDTree(xy)
X = G['lon_v']; Y = G['lat_v']
Maskv = G['mask_v'] # True over water
xy = np.array((X[Maskv],Y[Maskv])).T
xyT_v = cKDTree(xy)

pickle.dump(xyT_rho, open(outdir + 'xyT_rho.p', 'wb'))
pickle.dump(xyT_rho_un, open(outdir + 'xyT_rho_un.p', 'wb'))
pickle.dump(xyT_u, open(outdir + 'xyT_u.p', 'wb'))
pickle.dump(xyT_v, open(outdir + 'xyT_v.p', 'wb'))

if True:
    # 3D trees
    for tag in ['w', 'rho', 'u', 'v']:
        # prepare fields to make the tree
        tt0 = time()

        if tag == 'u':
            hh = (h[:,:-1] + h[:,1:])/2
        elif tag == 'v':
            hh = (h[:-1,:] + h[1:,:])/2
        elif tag in ['rho', 'w']:
            hh = h.copy()

        if tag in ['rho', 'u', 'v']:
            z = zrfun.get_z(hh, 0*hh, S, only_rho=True)
            x = G['lon_' + tag]
            y = G['lat_' + tag]
            mask = G['mask_' + tag]
        elif tag == 'w':
            z = zrfun.get_z(hh, 0*hh, S, only_w=True)
            x = G['lon_rho']
            y = G['lat_rho']
            mask = G['mask_rho']
        
        N,M,L = z.shape
        X = np.tile(x.reshape(1,M,L),[N,1,1])
        Y = np.tile(y.reshape(1,M,L),[N,1,1])
        H = np.tile(hh.reshape(1,M,L),[N,1,1])
        Z = z/H # fractional depth (-1 to 0)
    
        Mask = np.tile(mask.reshape(1,M,L),[N,1,1])
    
        xyz = np.array((X[Mask],Y[Mask],Z[Mask])).T
        
        print('Prepare fields to make tree %0.2f sec' % (time()-tt0))
        # create the nearest neighbor Tree objects
        tt0 = time()
    
        if tag == 'rho':
            xyzT_rho = cKDTree(xyz)
            pickle.dump(xyzT_rho, open(outdir + 'xyzT_rho.p', 'wb'))
        elif tag == 'u':
            xyzT_u = cKDTree(xyz)
            pickle.dump(xyzT_u, open(outdir + 'xyzT_u.p', 'wb'))
        elif tag == 'v':
            xyzT_v = cKDTree(xyz)
            pickle.dump(xyzT_v, open(outdir + 'xyzT_v.p', 'wb'))
        elif tag == 'w':
            xyzT_w = cKDTree(xyz)
            pickle.dump(xyzT_w, open(outdir + 'xyzT_w.p', 'wb'))

        print('Create 3D tree for %s: %0.2f sec' % (tag, time()-tt0))
        sys.stdout.flush()
    
