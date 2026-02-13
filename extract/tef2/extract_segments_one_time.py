"""
This code is is the main worker for extract_segments.py.

"""
from pathlib import Path
import sys
import argparse
import numpy as np
import xarray as xr
from time import time
import pandas as pd

from lo_tools import Lfun, zfun, zrfun
import tef_fun

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

tt0 = time()
print(fn)
    
ds = xr.open_dataset(fn)
# we get zeta for the DV calculation below
zeta = ds['zeta'][0,:,:].values

# Other 2-D quantities we will want for budgets:
two_d = dict()

# EminusP
#
# standard_name:   surface_upward_water_flux
# long_name:       modeled surface net freshwater flux, (E-P)/rhow
# units:           meter second-1
# negative_value:  upward flux, freshening (net precipitation)
# positive_value:  downward flux, salting (net evaporation)
if 'EminusP' in ds.data_vars:
    two_d['EminusP'] = ds.EminusP[0,:,:].values

# Surface salinity, to use with EminusP
two_d['salt_surf'] = ds.salt[0,-1,:,:].values

# shflux
#
# standard_name:   surface_downward_heat_flux_in_sea_water
# long_name:       surface net heat flux
# units:           watt meter-2
# negative_value:  upward flux, cooling
# positive_value:  downward flux, heating
if 'shflux' in ds.data_vars:
    two_d['shflux'] = ds.shflux[0,:,:].values

# set list of variables to extract
if args.get_bio:
    if 'NH4' in ds.data_vars:
        vn_list = tef_fun.vn_list
    else:
        # old roms version
        vn_list = ['salt', 'temp', 'oxygen',
            'NO3', 'phytoplankton', 'zooplankton', 'detritus', 'Ldetritus',
            'TIC', 'alkalinity']
else:
    vn_list = ['salt']
    
# Trim vn_list to only have variable in ds
for vn in vn_list:
    if vn not in ds.data_vars:
        vn_list.remove(vn)
    
# add custom 3-D variables, like salt-squared
vn_list = vn_list + ['salt2']

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
    this_DA = DA[jjj,iii]
    area = this_DA.sum()
    # store results
    A.loc[seg, 'volume'] = volume
    A.loc[seg, 'area'] = area
# 3-D tracers
for vn in vn_list:
    if vn == 'salt2':
        fld = ds.salt[0,:,:,:].values * ds.salt[0,:,:,:].values
    else:
        fld = ds[vn][0,:,:,:].values
    for seg in seg_list:
        jjj = j_dict[seg]
        iii = i_dict[seg]
        z_w = zrfun.get_z(h[jjj,iii], zeta[jjj,iii], S, only_w=True)
        dz = np.diff(z_w, axis=0)
        DV = dz * DA3[0,jjj,iii]
        volume = DV.sum()
        AA[vn] = (fld[:,jjj,iii] * DV).sum()/volume
        # store results
        A.loc[seg, vn] = AA[vn]
# 2-D properties, e.g. for surface fluxes
for vn in two_d.keys():
    fld = two_d[vn]
    for seg in seg_list:
        jjj = j_dict[seg]
        iii = i_dict[seg]
        this_DA = DA[jjj,iii]
        area = this_DA.sum()
        AA[vn] = (fld[jjj,iii] * this_DA).sum()/area
        # store results
        A.loc[seg, vn] = AA[vn]
ds.close()
            
print('  ** took %0.1f sec' % (time()-tt0))
sys.stdout.flush()

A.to_pickle(out_fn)
print('Time to extract all segment data = %0.2f sec' % (time()-tt0))
sys.stdout.flush()

