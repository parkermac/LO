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
DXu = (G['DX'][:,1:]+G['DX'][:,:-1])/2
DX3u = DXu.reshape((1,G['M'],G['L']-1))
DYv = (G['DY'][1:,:]+G['DY'][:-1,:])/2
DY3v = DYv.reshape((1,G['M']-1,G['L']))

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
salt = ds['salt'][0,:,:,:]
AKs = ds['AKs'][0,:,:,:]
KH = float(ds['nl_tnu2'][0].data)
zeta = ds['zeta'][0,:,:]
ot = ds['ocean_time'][:]
ds.close()

# calculate horizontal salinity gradient for hmix
# Erin's results say that centered differencing is fine, and this is used here.
sx2 = np.square(np.diff(salt,axis=2)/DX3u)
SX2 = 0*salt
SX2[:,:,1:-1] = (sx2[:,:,1:]+sx2[:,:,:-1])/2
sy2 = np.square(np.diff(salt,axis=1)/DY3v)
SY2 = 0*salt
SY2[:,1:-1,:] = (sy2[:,1:,:]+sy2[:,:-1,:])/2
    
dt = Lfun.modtime_to_datetime(ot.data[0])

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
    net_salt = (salt[:,jjj,iii] * DV).sum()
    AA['mean_salt'] = net_salt/volume
    net_salt2 = (salt[:,jjj,iii] * salt[:,jjj,iii] * DV).sum()
    AA['mean_salt2'] = net_salt2/volume
    
    dsdz = (salt[1:,jjj,iii] - salt[:-1,jjj,iii])/dzr
    AA['mix'] = -2*(AKs[1:-1,jjj,iii] * dsdz * dsdz * DVR).sum()
    
    AA['hmix'] = -2 * KH * ((SX2[:,jjj,iii] + SY2[:,jjj,iii]) * DV).sum()
    
    # store results
    for vn in ['mean_salt', 'mean_salt2', 'mix', 'hmix', 'volume']:
        A.loc[seg_name, vn] = AA[vn]
    
    if args.testing:
        print('%3s: Mean Salinity = %0.4f, Volume  = %0.4f km3' %
            (seg_name, AA['mean_salt'], AA['volume']/1e9))
        print('%3s: Mean Salinity Squared = %0.4f, Volume  = %0.4f km3' %
            (seg_name, AA['mean_salt2'], AA['volume']/1e9))
            
print('  ** took %0.1f sec' % (time()-tt0))
sys.stdout.flush()

pickle.dump(A, open(out_fn, 'wb'))
print('Time to extract all segment data = %0.2f sec' % (time()-tt0))
sys.stdout.flush()

