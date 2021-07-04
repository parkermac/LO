"""
A tool to extract hourly time series of volume and selected tracers and derived quantities
for budgets in the segments.

To test on mac:
run extract_segment_tracers -g cas6 -t v3 -x lo8b -ro 2 -0 2019.07.04 -1 2019.07.06 -test True

"""

from pathlib import Path
import sys
from datetime import datetime, timedelta

pth = Path(__file__).absolute().parent.parent.parent / 'alpha'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import extract_argfun as exfun
Ldir = exfun.intro() # this handles the argument passing

from time import time
import Lfun
import numpy as np
import netCDF4 as nc
import pickle
import pandas as pd
from subprocess import Popen as Po
from subprocess import PIPE as Pi
import zrfun
import tef_fun
if Ldir['testing']:
    from importlib import reload
    reload(tef_fun)
    
# set list of variables to extract
if Ldir['get_bio']:
    vn_list = tef_fun.vn_list
else:
    vn_list = ['salt']

ds0 = Ldir['ds0']
ds1 = Ldir['ds1']

tt00 = time()
print(' Doing segment extraction for '.center(60,'='))
print(' gtagex = ' + Ldir['gtagex'])
outname = 'segment_' + ds0 + '_' + ds1
print(' outname = ' + outname)

# make sure the output directory exists
out_dir = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef' / outname
Lfun.make_dir(out_dir, clean=True)

# make the scratch directory for holding temporary files
temp_dir = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / ('segment_temp_' + ds0 + '_' + ds1)
Lfun.make_dir(temp_dir, clean=True)

# get the DataFrame of all sections
gridname=Ldir['gtagex'].split('_')[0]
sect_df = tef_fun.get_sect_df(gridname)

# get list of history files to process
fn_list = Lfun.get_fn_list('hourly', Ldir, ds0, ds1)
NT = len(fn_list)

# -------------------- OLD ------------------

# get grid info
fn = fn_list[0]
G = zrfun.get_basic_info(fn, only_G=True)
S = zrfun.get_basic_info(fn, only_S=True)
h = G['h']
DA = G['DX'] * G['DY']
DA3 = DA.reshape((1,G['M'],G['L']))
DXu = (G['DX'][:,1:]+G['DX'][:,:-1])/2
DX3u = DXu.reshape((1,G['M'],G['L']-1))
DYv = (G['DY'][1:,:]+G['DY'][:-1,:])/2
DY3v = DYv.reshape((1,G['M']-1,G['L']))

# get volume info
vol_dir = Ldir['LOo'] / 'extract' / 'tef' / ('volumes_' + Ldir['gridname'])
v_df = pd.read_pickle(vol_dir / 'volumes.p')
bathy_dict = pickle.load(open(vol_dir / 'bathy_dict.p', 'rb'))
ji_dict = pickle.load(open(vol_dir / 'ji_dict.p', 'rb'))

seg_list = list(v_df.index)

if Ldir['testing']:
    verbose = True
    #seg_list = seg_list[:2]
else:
    verbose = False

j_dict = {}; i_dict = {}
for seg_name in seg_list:
    jj = []; ii = []
    ji_list_full = ji_dict[seg_name]
    for ji in ji_list_full:
        jj.append(ji[0])
        ii.append(ji[1])
    jjj = np.array(jj)
    iii = np.array(ii)
    j_dict[seg_name] = jjj
    i_dict[seg_name] = iii

s_df = pd.DataFrame(columns=seg_list)
s2_df = pd.DataFrame(columns=seg_list)
mix_df = pd.DataFrame(columns=seg_list)
hmix_df = pd.DataFrame(columns=seg_list)
v_df = pd.DataFrame(columns=seg_list)

for fn in fn_list:
    
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
    
    # find the volume and volume-mean salinity
    for seg_name in seg_list:
        
        jjj = j_dict[seg_name]
        iii = i_dict[seg_name]
        z_r, z_w = zrfun.get_z(h[jjj,iii], zeta[jjj,iii], S)
        dz = np.diff(z_w, axis=0)
        dzr = np.diff(z_r, axis=0)
        DV = dz * DA3[0,jjj,iii]
        DVR = dzr * DA3[0,jjj,iii]
        volume = DV.sum()
        net_salt = (salt[:,jjj,iii] * DV).sum()
        mean_salt = net_salt/volume
        net_salt2 = (salt[:,jjj,iii] * salt[:,jjj,iii] * DV).sum()
        mean_salt2 = net_salt2/volume
        
        dsdz = (salt[1:,jjj,iii] - salt[:-1,jjj,iii])/dzr
        mix = -2*(AKs[1:-1,jjj,iii] * dsdz * dsdz * DVR).sum()
        
        hmix = -2 * KH * ((SX2[:,jjj,iii] + SY2[:,jjj,iii]) * DV).sum()
        
        # store results
        s_df.loc[dt, seg_name] = mean_salt
        s2_df.loc[dt, seg_name] = mean_salt2
        mix_df.loc[dt, seg_name] = mix
        hmix_df.loc[dt, seg_name] = hmix
        v_df.loc[dt, seg_name] = volume
        
        if verbose:
            print('%3s: Mean Salinity = %0.4f, Volume  = %0.4f km3' %
                (seg_name, mean_salt, volume/1e9))
            print('%3s: Mean Salinity Squared = %0.4f, Volume  = %0.4f km3' %
                (seg_name, mean_salt2, volume/1e9))
                
    print('  ** took %0.1f sec' % (time()-tt0))
    sys.stdout.flush()

s_out_fn = out_dir / 'hourly_segment_salinity.p'
s2_out_fn = out_dir / 'hourly_segment_salinity2.p'
mix_out_fn = out_dir / 'hourly_segment_mix.p'
hmix_out_fn = out_dir / 'hourly_segment_hmix.p'
v_out_fn = out_dir / 'hourly_segment_volume.p'

s_df.to_pickle(s_out_fn)
s2_df.to_pickle(s2_out_fn)
mix_df.to_pickle(mix_out_fn)
hmix_df.to_pickle(hmix_out_fn)
v_df.to_pickle(v_out_fn)
    
        
    
