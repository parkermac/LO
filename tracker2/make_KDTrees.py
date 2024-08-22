"""
Creates and saves KDTrees for a given model run.

NOTE: for some reason I have to make the 2D trees individually (meaning
not in a loop like I do for the 3D trees).  No idea why.

Usage:

python make_KDTrees.py (uses defaults)

python make_KDTrees.py -gtx cas7_t0_x4b -d 2017.07.04 -ro 0

(** you have to be pointing it to a history file that exists **)

PERFORMANCE: takes about 30 sec on my mac for the cas7 grid.

"""

from lo_tools import Lfun, zrfun, zfun

import sys
from scipy.spatial import cKDTree
from time import time
import pickle
import numpy as np
import argparse

# optional command line arguments, can be input in any order
parser = argparse.ArgumentParser()
parser.add_argument('-gtx','--gtagex', default='cas7_t0_x4b', type=str)
parser.add_argument('-d', '--date_string', default='2017.07.04', type=str)
parser.add_argument('-ro', '--roms_out_num', type=int, default=0) # 1 = Ldir['roms_out1'], etc.

# get the args and put into Ldir
args = parser.parse_args()
argsd = args.__dict__
gridname, tag, ex_name = args.gtagex.split('_')
# get the dict Ldir
Ldir = Lfun.Lstart(gridname=gridname, tag=tag, ex_name=ex_name)
# add more entries to Ldir
for a in argsd.keys():
    if a not in Ldir.keys():
        Ldir[a] = argsd[a]
# set where to look for model output
if Ldir['roms_out_num'] == 0:
    pass
elif Ldir['roms_out_num'] > 0:
    Ldir['roms_out'] = Ldir['roms_out' + str(Ldir['roms_out_num'])]

fn = Ldir['roms_out'] / Ldir['gtagex'] / ('f' + Ldir['date_string']) / 'ocean_his_0002.nc'
outdir = Ldir['LOo'] / 'tracker2_trees' / Ldir['gridname']
Lfun.make_dir(outdir, clean=True)

G, S, T = zrfun.get_basic_info(fn)
h = G['h']

# 2D trees
X = G['lon_rho']; Y = G['lat_rho']
Maskr = G['mask_rho']==1 # True over water
xy = np.array((X[Maskr],Y[Maskr])).T
xyT_rho = cKDTree(xy)
xy = np.array((X.flatten(),Y.flatten())).T
xyT_rho_un = cKDTree(xy) # unmasked version
X = G['lon_u']; Y = G['lat_u']
Masku = G['mask_u']==1 # True over water
xy = np.array((X[Masku],Y[Masku])).T
xyT_u = cKDTree(xy)
X = G['lon_v']; Y = G['lat_v']
Maskv = G['mask_v']==1 # True over water
xy = np.array((X[Maskv],Y[Maskv])).T
xyT_v = cKDTree(xy)

pickle.dump(xyT_rho, open(outdir / 'xyT_rho.p', 'wb'))
pickle.dump(xyT_rho_un, open(outdir / 'xyT_rho_un.p', 'wb'))
pickle.dump(xyT_u, open(outdir / 'xyT_u.p', 'wb'))
pickle.dump(xyT_v, open(outdir / 'xyT_v.p', 'wb'))

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
        mask = G['mask_' + tag]==1
    elif tag == 'w':
        z = zrfun.get_z(hh, 0*hh, S, only_w=True)
        x = G['lon_rho']
        y = G['lat_rho']
        mask = G['mask_rho']==1
    
    N,M,L = z.shape
    X = np.tile(x.reshape(1,M,L),[N,1,1])
    Y = np.tile(y.reshape(1,M,L),[N,1,1])
    H = np.tile(hh.reshape(1,M,L),[N,1,1])
    Z = z/H # fractional depth (-1 to 0)

    Mask = np.tile(mask.reshape(1,M,L),[N,1,1])

    xyz = np.array((X[Mask],Y[Mask],Z[Mask])).T
    
    #------- jx
    # make a tree in the unit of meter
    x_c = np.mean(x)  # domain center
    y_c = np.mean(y)
    x_m, y_m = zfun.ll2xy(x, y, x_c, y_c)
    X_m = np.tile(x_m.reshape(1,M,L),[N,1,1])
    Y_m = np.tile(y_m.reshape(1,M,L),[N,1,1])
    Z_m = z  # surface = 0m, bottom = - local depth, surface elevation is not included
    xyz_m = np.array((X_m[Mask],Y_m[Mask],Z_m[Mask])).T
    #------- jx
    
    print('Prepare fields to make tree %0.2f sec' % (time()-tt0))
    # create the nearest neighbor Tree objects
    tt0 = time()

    if tag == 'rho':
        xyzT_rho = cKDTree(xyz)
        pickle.dump(xyzT_rho, open(outdir / 'xyzT_rho.p', 'wb'))
        xyzT_rho_m = cKDTree(xyz_m)  #jx
        pickle.dump(xyzT_rho_m, open(outdir / 'xyzT_rho_m.p', 'wb')) #jx
    elif tag == 'u':
        xyzT_u = cKDTree(xyz)
        pickle.dump(xyzT_u, open(outdir / 'xyzT_u.p', 'wb'))
        xyzT_u_m = cKDTree(xyz_m)  #jx
        pickle.dump(xyzT_u_m, open(outdir / 'xyzT_u_m.p', 'wb'))  #jx
    elif tag == 'v':
        xyzT_v = cKDTree(xyz)
        pickle.dump(xyzT_v, open(outdir / 'xyzT_v.p', 'wb'))
        xyzT_v_m = cKDTree(xyz_m) #jx
        pickle.dump(xyzT_v_m, open(outdir / 'xyzT_v_m.p', 'wb'))  #jx
    elif tag == 'w':
        xyzT_w = cKDTree(xyz)
        pickle.dump(xyzT_w, open(outdir / 'xyzT_w.p', 'wb'))
        xyzT_w_m = cKDTree(xyz_m) #jx
        pickle.dump(xyzT_w_m, open(outdir / 'xyzT_w_m.p', 'wb'))  #jx

    print('Create 3D tree for %s: %0.2f sec' % (tag, time()-tt0))
    sys.stdout.flush()
    
