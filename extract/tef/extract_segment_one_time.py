"""
This code extracts all needed segment data for one history file,
looping over all variables and all segments.

Performance: []
"""
from pathlib import Path
import sys
from datetime import datetime, timedelta
import argparse
import numpy as np
import netCDF4 as nc
from time import time
import pickle
import pandas as pd

pth = Path(__file__).absolute().parent.parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import Lfun
import zfun, zrfun

import tef_fun

# command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-pth', type=str) # path to LO_roms
parser.add_argument('-out_dir', type=str) # path to temporary directory for output
parser.add_argument('-gtagex', type=str) # gtagex
parser.add_argument('-d', type=str) # date string like 2019.07.04
parser.add_argument('-nhis', type=int, default=1) # history file number, 1-25
parser.add_argument('-get_bio', type=Lfun.boolean_string, default=True)
parser.add_argument('-testing', type=Lfun.boolean_string, default=False)
args = parser.parse_args()

Ldir = Lfun.Lstart(gridname=args.gtagex.split('_')[0])

nhiss = ('0000' + str(args.nhis))[-4:]
fn = Path(args.pth) / args.gtagex / ('f' + args.d) / ('ocean_his_' + nhiss + '.nc')
out_dir = Path(args.out_dir)
Lfun.make_dir(out_dir)
out_fn  = out_dir / ('A_' + args.d + '_' + nhiss + '.p')
out_fn.unlink(missing_ok=True)

# ---
# get grid info
G = zrfun.get_basic_info(fn, only_G=True)
S = zrfun.get_basic_info(fn, only_S=True)
h = G['h']
DA = G['DX'] * G['DY']
DA3 = DA.reshape((1,G['M'],G['L']))

# get segment info
vol_dir = Ldir['LOo'] / 'extract' / 'tef' / ('volumes_' + Ldir['gridname'])
v_df = pd.read_pickle(vol_dir / 'volumes.p')
j_dict = pickle.load(open(vol_dir / 'j_dict.p', 'rb'))
i_dict = pickle.load(open(vol_dir / 'i_dict.p', 'rb'))
seg_list = list(v_df.index)

# set list of variables to extract
if args.get_bio:
    vn_list = tef_fun.vn_list
else:
    vn_list = ['salt']

tt0 = time()
print(fn)
    
ds = nc.Dataset(fn)

vn_dict = {}
for vn in vn_list:
    vn_dict[vn] = ds[vn][0,:,:,:]
zeta = ds['zeta'][0,:,:]
ds.close()

# find the volume and other variables for each segment, at this time
A = pd.DataFrame(index=seg_list)
AA = dict()
for seg_name in seg_list:
    
    jjj = j_dict[seg_name]
    iii = i_dict[seg_name]
    z_r, z_w = zrfun.get_z(h[jjj,iii], zeta[jjj,iii], S)
    dz = np.diff(z_w, axis=0)
    dzr = np.diff(z_r, axis=0)
    DV = dz * DA3[0,jjj,iii]
    DVR = dzr * DA3[0,jjj,iii]
    volume = DV.sum()
    AA['volume'] = volume
    
    for vn in vn_list:
        AA[vn] = (vn_dict[vn][:,jjj,iii] * DV).sum()/volume
    # store results
    for vn in vn_list + ['volume']:
        A.loc[seg_name, vn] = AA[vn]
            
print('  ** took %0.1f sec' % (time()-tt0))
sys.stdout.flush()

pickle.dump(A, open(out_fn, 'wb'))
print('Time to extract all segment data = %0.2f sec' % (time()-tt0))
sys.stdout.flush()

