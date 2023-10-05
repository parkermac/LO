"""
This code extracts all needed segment data for one history file,
looping over all variables and all segments.

To test on mac:
...

"""
from pathlib import Path
import sys
import argparse
import numpy as np
import xarray as xr
from time import time
import pandas as pd

from lo_tools import Lfun, zfun, zrfun
#import tef_fun

# command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('-in_fn', type=str) # path to a history file
parser.add_argument('-out_dir', type=str) # path to temporary directory for output
parser.add_argument('-file_num', type=str) # for numbering the output file
parser.add_argument('-seg_fn', type=str) # path to pickled seg_info_dict
parser.add_argument('-get_bio', type=Lfun.boolean_string, default=True)
parser.add_argument('-testing', type=Lfun.boolean_string, default=False)
args = parser.parse_args()

fn = Path(args.in_fn)
out_dir = Path(args.out_dir)
out_str = ('000000' + args.file_num)[-6:]

Ldir = Lfun.Lstart()

seg_info_dict = pd.read_pickle(args.seg_fn)

out_dir = Path(args.out_dir)
out_fn  = out_dir / ('A_' + out_str + '.p')
print(out_fn)
out_fn.unlink(missing_ok=True)

# ---
# get grid info
G, S, T = zrfun.get_basic_info(fn)
h = G['h']
DA = G['DX'] * G['DY']
DA3 = DA.reshape((1,G['M'],G['L']))

j_dict = dict(); i_dict = dict()
seg_list = seg_info_dict.keys()
for seg in seg_list:
    ji_list = seg_info_dict[seg]['ji_list']
    # make index vectors for fancy indexing
    jj = []; ii = []
    for ji in ji_list:
        jj.append(ji[0])
        ii.append(ji[1])
    JJ = np.array(jj,dtype=int)
    II = np.array(ii,dtype=int)
    j_dict[seg] = JJ
    i_dict[seg] = II

# set list of variables to extract
if args.get_bio:
    # NEED to deal with this!
    vn_list = 'junk'#tef_fun.vn_list
else:
    vn_list = ['salt']

tt0 = time()
print(fn)
    
ds = xr.open_dataset(fn)

# This pre-loading seems like a bad idea when I am doing a lot of bio vars!
# Also soon we will have to deal with other things like air-sea fluxes.
vn_dict = {}
for vn in vn_list:
    vn_dict[vn] = ds[vn][0,:,:,:].values
zeta = ds['zeta'][0,:,:].values
ds.close()

# find the volume and other variables for each segment, at this time
A = pd.DataFrame(index=seg_list)
AA = dict()
for seg in seg_list:
    
    jjj = j_dict[seg]
    iii = i_dict[seg]
    z_w = zrfun.get_z(h[jjj,iii], zeta[jjj,iii], S, only_w=True)
    dz = np.diff(z_w, axis=0)
    DV = dz * DA3[0,jjj,iii]
    volume = DV.sum()
    AA['volume'] = volume
    
    for vn in vn_list:
        AA[vn] = (vn_dict[vn][:,jjj,iii] * DV).sum()/volume
    # store results
    for vn in vn_list + ['volume']:
        A.loc[seg, vn] = AA[vn]
            
print('  ** took %0.1f sec' % (time()-tt0))
sys.stdout.flush()

A.to_pickle(out_fn)
print('Time to extract all segment data = %0.2f sec' % (time()-tt0))
sys.stdout.flush()

